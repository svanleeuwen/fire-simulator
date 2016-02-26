#define DEBUG
#define VIENNACL_WITH_EIGEN

#define G -9.81f
#define AMBIENT_TEMP 273.0f
#define AIR_DENSITY 1.3f

/* Guessed this value based on Figure 6 from 
 * "Effective Density Characterization of Soot
 * Agglomerates from Various Sources and
 * Comparison to Aggregation Theory"
 */
#define SOOT_DENSITY 400.0f

#define ALPHA 17.0f
#define BETA 1.0f

#include "mac_grid.hpp"
#include "viennacl/linalg/cg.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <stdlib.h>

using std::cout;
using std::cerr;
using std::endl;

using Eigen::Vector4f;

/* *********** First Fire Paper ************
 *
 * TODO: Read Bridson and paper
 *
 */

/* *********** Major Goals ****************
 *
 * TODO: First Fire Paper
 * TODO: Second Fire Paper
 *
 * TODO: Convert to 3D 
 *           - Vorticy Confinement (can't do it in 2D)
 *
 * TODO: Set up GitHub repository
 * TODO: Put a .pdf of references in the code directory for github
 *
 */

/* *********** Extra Stuff **************
 *
 * TODO: Consider rk3 (90% sure you should be doing this since rk3 helps preserve vortices)
 * TODO: May be worth implementing your own solver for speedup
 * TODO: Handle solids
 * TODO: Maybe switch linear system in project() 
 *       to double precision
 * TODO: Precondition the matrix
 * TODO: Maybe use particle methods for rendering
 * TODO: Set up UI to tweak alpha/beta for buoyency and phi for vorticy confinment
 *
 */

/********** Optional smoke stuff ************
 *
 * TODO: Update cubic interpolator using Appendix B of Fedkiw or Bridson's course notes
 * TODO: Decay (maybe look to Fedkiw to see if this is used in theirs)
 * TODO: Diffusion (p100 Bridson) (Look to viscous fluids if
 *you want to do this properly)
 * TODO: Variable Density (not necessary for the most part)
 * TODO: Divergence solves (only needed if variable density is implemented)
 **
 */

//**************** Testing *******************
#ifdef DEBUG
void MacGrid::setTestValues()
{
    for(uint i = 2; i < m_press.size() - 2; ++i)
    {
        for(uint j = 2; j < m_press[0].size() - 2; ++j)
        {
            m_press[i][j] = 1;
            m_mat[i][j] = Material::Fluid;
        }
    }

    for(uint i = 2; i < m_vel.u_.size() - 2; ++i) {
        for(uint j = 2; j < m_vel.u_[0].size() - 2; ++j) {
            m_vel.u_[i][j] = 1.0f;
        }
    }

    for(uint i = 2; i < m_vel.v_.size() - 2; ++i) {
        for(uint j = 2; j < m_vel.v_[0].size() - 2; ++j) {
            m_vel.v_[i][j] = 1.0f;
        }
    }

    for(int i = 2; i < m_nx - 2;  ++i) {
        for(int j = 2; j < m_ny - 2; ++j) {
    //        m_temp[i].push_back((i + j) % 3);
        }
    }

    m_src[(int)(m_nx / 2.0f)- 1][2] = true;
    m_src[(int)(m_nx / 2.0f)][2] = true;

    /*cout << "Initialized grids" << endl << endl;
    cout << "P :" << endl << m_press;
    cout << endl << m_vel;
    cout << endl << "T :" << endl << m_temp;
    cout << endl << "S :" << endl << m_soot;
    cout << endl;*/
}

void MacGrid::testSolver()
{
    setTestValues();

    int size = (m_nx-PADDING_WIDTH*2)*(m_ny-PADDING_WIDTH*2);
    
    SparseMatrix<float, Eigen::RowMajor> sm(size, size);
    sm.reserve(5 * size);
    
    VectorXf rhs(size);
    buildSystem(sm, rhs);
    
//    cout << sm << endl;
//    cout << rhs << endl;

    VectorXf result(size);
    solveSystem(sm, rhs, result);

//    cout << result << endl;
}
#endif

//**************** Public ******************
MacGrid::MacGrid(int res, float dt) :
    m_vel(res, res), m_press(res, res, 0.0f), m_temp(res, res, AMBIENT_TEMP), 
    m_soot(res, res, 0.0f)
{
    m_nx = res + PADDING_WIDTH * 2;
    m_ny = res + PADDING_WIDTH * 2;
    m_dt = dt;
    m_dx = 1.0 * 1.0f / res;

    m_mat.resize(m_nx);
    m_src.resize(m_nx);

    for(int i = 0; i < m_nx; ++i) 
    {
        m_mat[i].resize(m_ny, Material::Empty);
        m_src[i].resize(m_ny, false);
    }

    setTestValues();
}

void MacGrid::step() 
{
    addSmoke();
    advect();

    /*cout << "Advected grids" << endl << endl;
    cout << "P :" << endl << m_press;
    cout << endl << m_vel;
    cout << endl << "T :" << endl << m_temp;
    cout << endl << "S :" << endl << m_soot;
    cout << endl;*/

    applyForces();
    project();

/*    cout << "Projected grids" << endl << endl;
    cout << "P :" << endl << m_press;
    cout << endl << m_vel;
    cout << endl << "T :" << endl << m_temp;
    cout << endl << "S :" << endl << m_soot;
    cout << endl;*/
}

float MacGrid::getDensity(int i, int j)
{
    int a = i + PADDING_WIDTH;
    int b = j + PADDING_WIDTH;

    return m_soot[a][b];
}

void MacGrid::addSmoke()
{
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            if(m_src[i][j])
            {
                 m_temp[i][j] += 10.0f;
                 m_soot[i][j] += 0.50f;
            }
        }
    }
}

//******************* Advection **********************

void MacGrid::advect() 
{
    m_temp.advect(m_vel, m_dt);
    m_soot.advect(m_vel, m_dt);
    m_vel.advect(m_dt);
}

//********************* Forces ***************************

void MacGrid::applyForces() 
{
    // alpha and beta have values from p.103 of Bridson
    float alpha = ALPHA;// (SOOT_DENSITY - AIR_DENSITY) / AIR_DENSITY;
    float beta = BETA; // 1.0f / AMBIENT_TEMP;

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; 
                j < m_ny + 1 - PADDING_WIDTH; ++j)
        {
            float soot = (m_soot[i][j-1] 
                    + m_soot[i][j]) / 2.0f;
            float temp = (m_temp[i][j-1] 
                    + m_temp[i][j]) / 2.0f;
            float buoy = (alpha * soot - beta 
                    * (temp - AMBIENT_TEMP)) * G;

            m_vel.v_[i][j] = m_vel.v_[i][j] + m_dt * buoy;
        }
    }
}

//******************* Projection **************************

void MacGrid::project() 
{
    int size = (m_nx-PADDING_WIDTH*2)*(m_ny-PADDING_WIDTH*2);
    
    SparseMatrix<float, Eigen::RowMajor> sm(size, size);
    sm.reserve(5 * size);
    
    VectorXf rhs(size);
    buildSystem(sm, rhs);

    VectorXf result(size);
    solveSystem(sm, rhs, result);

    updatePressure(result);
    applyPressure();
}

void MacGrid::addTriplets(const Mat &A_diag, const Mat &A_x, 
        const Mat &A_y, vector<Trip> &trip)
{
    int width = A_diag.size();
    int height = A_diag[0].size();

    for(int i = 0; i < width; ++i)
    {
        for(int j = 0; j < height; ++j)
        {
            int index = i*height + j;
            
            if(A_diag[i][j] > EPSILON)
            {
                trip.push_back(Trip(index, index, A_diag[i][j]));
            }

            if(A_x[i][j] < EPSILON)
            {
                if(index < (width-1)*height)
                {
                    trip.push_back(Trip(index, index + height,
                                A_x[i][j]));
                    trip.push_back(Trip(index + height, index,
                                A_x[i][j]));
                }
            }

            if(A_y[i][j] < EPSILON)
            {
                if(index < width*height - 1)
                {
                    trip.push_back(Trip(index, index + 1,
                                A_y[i][j]));
                    trip.push_back(Trip(index + 1, index,
                                A_y[i][j]));
                }
            }
        }
    }
}

// Algorithm from p.78 of Bridson, ignoring solids
void MacGrid::addEntries(Mat &A_diag, Mat &A_x, Mat &A_y, 
            const Vector2i &g_point)
{
    int i = g_point[0];
    int j = g_point[1];

    int a = i + PADDING_WIDTH;
    int b = j + PADDING_WIDTH;

    float scale = m_dt / (AIR_DENSITY * m_dx * m_dx);
    A_diag[i][j] = 4.0f * scale;

    if(m_mat[a+1][b] == Material::Fluid)
    {
        A_x[i][j] = -scale;
    }
    
    if(m_mat[a][b+1] == Material::Fluid)
    {
        A_y[i][j] = -scale;
    }
}

// Divergence calculation from p.72 of Bridson, ignoring scale,
//      to build rhs
void MacGrid::buildSystem(
        SparseMatrix<float, Eigen::RowMajor> &sm, 
        VectorXf &rhs)
{
    int width = m_nx - PADDING_WIDTH*2;
    int height = m_ny - PADDING_WIDTH*2;
    float scale = 1.0f / m_dx;
   
    Mat A_diag(width, 
            vector<float>(height, 0.0f));
    Mat A_x = A_diag;
    Mat A_y = A_diag;
 
    vector<Trip> trip;
    trip.reserve(5 * width * height);

    for(int i = 0; i < width; ++i)
    {
        int a = i + PADDING_WIDTH;

        for(int j = 0; j < height; ++j)
        {
            int index = i * (height) + j;
            int b = j + PADDING_WIDTH;

            if(m_mat[a][b] == Material::Fluid)
            {
                rhs[index] = -scale * 
                    ((m_vel.u_[a+1][b] - m_vel.u_[a][b])
                        + (m_vel.v_[a][b+1] - m_vel.v_[a][b]));

                addEntries(A_diag, A_x, A_y, Vector2i(i, j));
            }
            else
            {
                rhs[index] = 0.0f;
            }
        }
    }

    addTriplets(A_diag, A_x, A_y, trip);
    sm.setFromTriplets(trip.begin(), trip.end());
}

void MacGrid::solveSystem(
        const SparseMatrix<float, Eigen::RowMajor> &sm, 
        const VectorXf &rhs, VectorXf &result)
{
    result = viennacl::linalg::solve(sm, rhs, viennacl::linalg::cg_tag());
}

void MacGrid::updatePressure(const VectorXf &result)
{
    int size = (m_nx-PADDING_WIDTH*2)*(m_ny-PADDING_WIDTH*2);
    int height = m_ny - PADDING_WIDTH * 2;

    for(int k = 0; k < size; ++k)
    {
        int i = (k / height) + PADDING_WIDTH;
        int j = (k % height) + PADDING_WIDTH;

        m_press[i][j] = result[k];
    }
}

// Algorithm from p.71 of Bridson, ignoring possibility of solid
// cells
void MacGrid::applyPressure()
{
    float scale = m_dt / (AIR_DENSITY * m_dx);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            m_vel.u_[i][j] -= scale * (m_press[i][j] - m_press[i-1][j]);
            m_vel.v_[i][j] -= scale * (m_press[i][j] - m_press[i][j-1]);
        }
    }
}
