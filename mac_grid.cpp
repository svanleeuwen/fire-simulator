#define CIRCLE 
#define VIENNACL_WITH_EIGEN

#define G -9.81f
#define AIR_DENSITY 1.3f // kg/m^3

/* Guessed this value based on Figure 6 from 
 * "Effective Density Characterization of Soot
 * Agglomerates from Various Sources and
 * Comparison to Aggregation Theory"
 *
 * #define SOOT_DENSITY 400.0f // kg/m^3
 * Need to really reduce density for this to work
 *
 */

#define SOOT_DENSITY FLAME_DENSITY

// alpha and beta have values from p.103 of Bridson
#define ALPHA (SOOT_DENSITY - AIR_DENSITY) / AIR_DENSITY
#define BETA 1.0f / AMBIENT_TEMP

#define VORTICITY_EPS 1.0f

#include "mac_grid.hpp"
#include "viennacl/linalg/cg.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <stdlib.h>

using std::cout;
using std::cerr;
using std::endl;

using Eigen::Vector4f;


/* *********** Major Goals ****************
 *
 * TODO: Write goals for second fire paper
 *
 * TODO: Convert to 3D 
 * TODO: Vorticy Confinement
 *
 * TODO: Set up GitHub repository
 * TODO: Put a .pdf of references in the code directory for github
 *
 */

/************ Fire Goals ********************
 *
 * TODO: Variable density for the matrix: smoke + p136
 *
 */


/********** Smoke goals ************
 *
 * TODO: Diffusion (p100 Bridson) (Look to viscous fluids if
 *       you want to do this properly)
 * TODO: Variable Density (not necessary for the most part) <<<< Probably necessary for fire
 * TODO: Divergence solves (only needed if variable density is implemented)
 **
 */

/* *********** Extra Stuff **************
 *
 * TODO: Handle solids
 * TODO: Precondition the matrix
 * TODO: Set up UI to tweak alpha/beta for buoyency and phi for vorticy confinment
 *
 */


//**************** Testing *******************
void MacGrid::setTestValues()
{
    for(uint i = 2; i < m_press.size() - 2; ++i)
    {
        for(uint j = 2; j < m_press[0].size() - 2; ++j)
        {
            m_press[i][j] = 1;
            m_mat[i][j] = Material::Fire;
        }
    }
/*
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
*/
    for(int i = 2; i < m_nx - 2;  ++i) {
        for(int j = 2; j < m_ny - 2; ++j) {
    //        m_temp[i].push_back((i + j) % 3);
        }
    }

#ifdef FLAMETHROWERS
    for(int i = PADDING_WIDTH + 3; 
            i < PADDING_WIDTH + 8; ++i)
    {
        for(int j = PADDING_WIDTH + 13; 
                j < PADDING_WIDTH + 18; ++j)
        {
            m_src[i][j] = true;
       //     m_src[m_nx - i][j] = true;
        }
    }
#endif

#ifdef CIRCLE
    for(int i = PADDING_WIDTH + 8; i < m_nx - PADDING_WIDTH - 8; i += 8)
    {
        m_interface->addCircle(2.0f, m_vel, Vector2f(i, 15), 
                &m_src);
    }               
#endif

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
    
    SparseMatrix<double, Eigen::RowMajor> sm(size, size);
    sm.reserve(5 * size);
    
    VectorXd rhs(size);
    buildSystem(sm, rhs);
    
//    cout << sm << endl;
//    cout << rhs << endl;

    VectorXd result(size);
    solveSystem(sm, rhs, result);

//    cout << result << endl;
}

//**************** Public ******************
MacGrid::MacGrid(int res, float dt, bool fire) :
    m_vel(res, res), m_press(res, res, 0.0f), 
    m_temp(res, res, AMBIENT_TEMP), 
    m_soot(res, res, 0.0f)
{
    m_nx = res + PADDING_WIDTH * 2;
    m_ny = res + PADDING_WIDTH * 2;
    m_dt = dt;
    m_dx = 1.0f / res;

    m_mat.resize(m_nx);
    m_src.resize(m_nx);

    for(int i = 0; i < m_nx; ++i) 
    {
        m_mat[i].resize(m_ny, Material::Empty);
        m_src[i].resize(m_ny, false);
    }

    if(fire)
    {
        m_interface = new LevelSet(res, res);
    }
    else
    {
        m_interface = NULL;
    }

    setTestValues();
}

void MacGrid::step() 
{
    addFuel();
    advect();

/*    cout << "Advected grids" << endl << endl;
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

bool MacGrid::isFuel(int i, int j)
{
    int a = i + PADDING_WIDTH;
    int b = j + PADDING_WIDTH;
    
    return m_interface->at(a)[b] < 0;
}

float MacGrid::getTemp(int i, int j)
{
    int a = i + PADDING_WIDTH;
    int b = j + PADDING_WIDTH;

    return m_temp[a][b];
}

void MacGrid::addSmoke()
{
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            if(m_src[i][j])
            {
                 m_temp[i][j] += 100.0f;
                 m_soot[i][j] += 0.10f;
            }
        }
    }
}

void MacGrid::addFuel()
{
#ifdef CIRCLE
    for(int i = PADDING_WIDTH + 8; i < m_nx - PADDING_WIDTH - 8; i += 8)
    {
        m_interface->addCircle(2.0f, m_vel, Vector2f(i, 15), 
                &m_src);
    }               
/*
    m_interface->addCircle(4.0f, m_vel, 
            Vector2f(m_nx/2.0f, 15));
    
    m_interface->addCircle(4.0f, m_vel, 
            Vector2f(m_nx/2.0f + 20, 15));
    m_interface->addCircle(4.0f, m_vel, 
            Vector2f(m_nx/2.0f - 20, 15));*/
#endif
    
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            if(m_src[i][j])
            {
#ifdef FLAMETHROWERS
                m_interface->at(i)[j] = -1;

                if(i < (int)(m_nx / 2.0))
                {
                    m_vel.u_[i][j] = 100;
                }
                else
                {
                    m_vel.u_[i][j] = -100;
                }
                m_vel.v_[i][j] = 0;
#endif

#ifdef CIRCLE
                m_vel.v_[i][j] = 30.0f;
#endif
            }
        }
    }

    m_interface->redistance();
}

//******************* Advection **********************

static void clampSoot(QGrid &soot, int nx, int ny)
{
    for(int i = 0; i < nx; ++i)
    {
        for(int j = 0; j < ny; ++j)
        {
            if(soot[i][j] < 0.0f)
            {
                soot[i][j] = 0.0f;
            }
            else if(soot[i][j] > 1.0f)
            {
                soot[i][j] = 1.0f;
            }
        }
    }
}

void MacGrid::advect() 
{
    if(m_interface != NULL)
    {
        m_interface->advect(m_dt, m_vel);
        m_interface->redistance();
    }
    m_temp.advect(m_vel, m_dt, m_interface);

    m_soot.advect(m_vel, m_dt, m_interface);
    clampSoot(m_soot, m_nx, m_ny);

    m_vel.advect(m_dt, m_interface);
}

//********************* Forces ***************************

void MacGrid::applyForces()
{
    applyBuoyancy();
    applyVorticity();
}

void MacGrid::applyBuoyancy() 
{
    float alpha = ALPHA;
    float beta = BETA;

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

            m_vel.v_[i][j] += m_dt * buoy;
        }
    }
}

// 2D version of vorticity from p167 of Bridson
void MacGrid::applyVorticity()
{
    Grid omega(m_nx - PADDING_WIDTH, m_ny - PADDING_WIDTH, 0.0f);
    
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; 
                j < m_ny - PADDING_WIDTH; ++j)
        {
            omega[i][j] = getOmega(i, j);
        }
    }

    QGrid vortX(m_nx - PADDING_WIDTH, 
            m_ny - PADDING_WIDTH, 0.0f);
    QGrid vortY(m_nx - PADDING_WIDTH, 
            m_ny - PADDING_WIDTH, 0.0f);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; 
                j < m_ny - PADDING_WIDTH; ++j)
        {
            Vector2f gradOmega = getGradOmega(omega, i, j);

            float gradOmegaNorm = sqrt(gradOmega[0] * 
                    gradOmega[0] + gradOmega[1] * gradOmega[1]);
            float denom = gradOmegaNorm + 10e-20 * (1.0f / 
                    (m_dx * m_dt));

            Vector2f N = {gradOmega[0] / denom,
                gradOmega[1] / denom};

            // Here, we take the cross product of
            //      (0, 0, omega) and (N[0], N[1], 0)
            Vector2f f = VORTICITY_EPS * m_dx * 
                Vector2f(-omega[i][j] * N[1], omega[i][j] * N[0]);

            vortX[i][j] = f[0];
            vortY[i][j] = f[1];
        }
    }

    for(int i = PADDING_WIDTH; i < m_nx + 1 - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; 
                j < m_ny + 1 - PADDING_WIDTH; ++j)
        {
            if(i < m_nx - PADDING_WIDTH)
            {
                m_vel.u_[i][j] += m_dt * vortX.getQuantity(
                        Vector2f(i-0.5f, j));
            }
            if(j < m_ny - PADDING_WIDTH)
            {
                m_vel.v_[i][j] += m_dt * vortY.getQuantity(
                        Vector2f(i, j-0.5f));
            }
        }
    }
}

float MacGrid::getOmega(int i, int j)
{
    float t2;
    if(i == m_nx - PADDING_WIDTH - 1)
    {
        t2 =  (m_vel.getVelocity(Vector2f(i, j))[0] - 
            m_vel.getVelocity(Vector2f(i-1, j))[0]) 
            / (2.0f*m_dx);
    }
    else if(i == PADDING_WIDTH)
    {
        t2 =  (m_vel.getVelocity(Vector2f(i+1, j))[0] - 
            m_vel.getVelocity(Vector2f(i, j))[0]) 
            / (2.0f*m_dx);
    }
    else
    {
        t2 =  (m_vel.getVelocity(Vector2f(i+1, j))[0] - 
            m_vel.getVelocity(Vector2f(i-1, j))[0]) 
            / (2.0f*m_dx);
    }

    float t1;
    if(j == m_ny - PADDING_WIDTH - 1)
    {
        t1 = (m_vel.getVelocity(Vector2f(i, j))[1] - 
            m_vel.getVelocity(Vector2f(i, j-1))[1]) 
            / (2.0f*m_dx);
    }
    else if(j == 0)
    {
        t1 = (m_vel.getVelocity(Vector2f(i, j+1))[1] - 
            m_vel.getVelocity(Vector2f(i, j))[1]) 
            / (2.0f*m_dx);
    }
    else
    {
        t1 = (m_vel.getVelocity(Vector2f(i, j+1))[1] - 
            m_vel.getVelocity(Vector2f(i, j-1))[1]) 
            / (2.0f*m_dx);
    }

    return t1 - t2;
}

Vector2f MacGrid::getGradOmega(Grid &omega, int i, int j)
{
    float gradOmegaX;
    if(i == m_nx - PADDING_WIDTH - 1)
    {
        gradOmegaX = (omega[i][j] - omega[i-1][j])
            / (2.0f * m_dx);
    }
    else if(i == PADDING_WIDTH)
    {
        gradOmegaX = (omega[i+1][j] - omega[i][j])
            / (2.0f * m_dx);
    }
    else
    {
        gradOmegaX = (omega[i+1][j] - omega[i-1][j])
            / (2.0f * m_dx);
    }

    float gradOmegaY;
    if(j == m_ny - PADDING_WIDTH - 1)
    {
        gradOmegaY = (omega[i][j] - omega[i][j-1])
            / (2.0f * m_dx);
    }
    else if(j == 0)
    {
        gradOmegaY = (omega[i][j+1] - omega[i][j])
            / (2.0f * m_dx);
    }
    else
    {
        gradOmegaY = (omega[i][j+1] - omega[i][j-1])
            / (2.0f * m_dx);
    }

    return Vector2f(gradOmegaX, gradOmegaY);
}

//******************* Projection **************************

void MacGrid::project() 
{
    int size = (m_nx-PADDING_WIDTH*2)*(m_ny-PADDING_WIDTH*2);
    
    SparseMatrix<double, Eigen::RowMajor> sm(size, size);
    sm.reserve(5 * size);
    
    VectorXd rhs(size);
    buildSystem(sm, rhs);

    VectorXd result(size);
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

    if(m_mat[a+1][b] == Material::Fire)
    {
        A_x[i][j] = -scale;
    }
    
    if(m_mat[a][b+1] == Material::Fire)
    {
        A_y[i][j] = -scale;
    }
}

// Divergence calculation from p.72 of Bridson, ignoring scale,
//      to build rhs
void MacGrid::buildSystem(
        SparseMatrix<double, Eigen::RowMajor> &sm, 
        VectorXd &rhs)
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

            if(m_mat[a][b] == Material::Fire)
            {

                rhs[index] = -scale 
                    * getCorrectedRhs(Vector2i(a, b));

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

float MacGrid::getCorrectedRhs(Vector2i pos)
{
    int a = pos[0];
    int b = pos[1];

    float l_jump = 0.0f;
    float r_jump = 0.0f;
    float t_jump = 0.0f;
    float b_jump = 0.0f;

    if(m_interface != NULL)
    {
        float dv = (FUEL_DENSITY 
                / FLAME_DENSITY - 1.0f) * BURN_RATE;

        bool isFuel = m_interface->at(a)[b] < 0;

        bool isLeftFuel = 
            m_interface->lerp(Vector2f(a-0.5f, b)) < 0; 
        
        if(isFuel && !isLeftFuel)
        {
            l_jump -= dv 
                * m_interface->getGradientX(Vector2f(a, b));
        }
        else if(!isFuel && isLeftFuel)
        {
            l_jump += dv 
                * m_interface->getGradientX(Vector2f(a, b));
        }

        bool isRightFuel = 
            m_interface->lerp(Vector2f(a+0.5f, b)) < 0; 

        if(isFuel && !isRightFuel)
        {
            r_jump -= dv 
                * m_interface->getGradientX(Vector2f(a+1, b));
        }
        else if(!isFuel && isRightFuel)
        {
            r_jump += dv 
                * m_interface->getGradientX(Vector2f(a+1, b));
        }

        bool isBottomFuel = 
            m_interface->lerp(Vector2f(a, b-0.5f)) < 0; 

        if(isFuel && !isBottomFuel)
        {
            b_jump -= dv
                * m_interface->getGradientY(Vector2f(a, b));
        }
        else if(!isFuel && isBottomFuel)
        {
            b_jump += dv
                * m_interface->getGradientY(Vector2f(a, b));
        }

        bool isTopFuel = 
            m_interface->lerp(Vector2f(a, b+0.5f)) < 0; 
        
        if(isFuel && !isTopFuel)
        {
            t_jump -= dv 
                * m_interface->getGradientY(Vector2f(a, b+1));
        }
        else if(!isFuel && isTopFuel)
        {
            t_jump += dv 
                * m_interface->getGradientY(Vector2f(a, b+1));
        }
    }

    float left = m_vel.u_[a][b] + l_jump;
    float right = m_vel.u_[a+1][b] + r_jump;
    float top = m_vel.v_[a][b+1] + t_jump;
    float bottom = m_vel.v_[a][b] + b_jump;

    return (right - left + top - bottom);
}

void MacGrid::solveSystem(
        const SparseMatrix<double, Eigen::RowMajor> &sm, 
        const VectorXd &rhs, VectorXd &result)
{
    result = viennacl::linalg::solve(sm, rhs, viennacl::linalg::cg_tag());
}

void MacGrid::updatePressure(const VectorXd &result)
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
            m_vel.u_[i][j] -= scale * (m_press[i][j] 
                    - m_press[i-1][j]);
            m_vel.v_[i][j] -= scale * (m_press[i][j] 
                    - m_press[i][j-1]);
        }
    }
}
