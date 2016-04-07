#define CIRCLE
#define VIENNACL_WITH_EIGEN

// Constants defined in Hong's paper
/*
#define DCJ 0.2
#define C_1 10
#define C_2 0.3
#define C_3 300
#define C_4 0
#define MU_THETA 2
#define C_5_THETA 2.5
*/

#define DCJ 0.01
#define C_1 0.1
#define C_2 0.01
#define C_3 0.1
#define C_4 0
#define MU_THETA 0.01
#define C_5_THETA 0.5

#define GRAVITY -9.81f
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
#include <Eigen/Dense>
#include <stdlib.h>

using std::cout;
using std::cerr;
using std::endl;

using Eigen::Vector4f;
using Eigen::Vector3d;


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

//**************** Public ******************
MacGrid::MacGrid(int res, float dt, float scale) :
    m_dt(dt),
    m_scale(scale),
    m_dx(scale/res),
    
    m_interface(res),
    m_curvature(res, 0.0f, m_dx),
   
    m_burn(res, BURN_RATE, m_dx),
    m_dburn(res, 0.0f, m_dx),
    
    m_vel(res, res, res, m_dx, &m_interface, &m_burn), 
    m_press(res, 1.0f, m_dx), 
    m_temp(res, AMBIENT_TEMP, m_dx), 
    m_soot(res, 0.0f, m_dx)
{
    m_nx = res + PADDING_WIDTH * 2;
    m_ny = res + PADDING_WIDTH * 2;
    m_nz = res + PADDING_WIDTH * 2;

    m_src.resize(m_nx);

    for(int i = 0; i < m_nx; ++i) 
    {
        m_src[i].resize(m_ny);

        for(int j = 0; j < m_ny; ++j)
        {
            m_src[i][j].resize(m_nz, false);
        }
    }

    setTestValues();
}

//**************** Testing *******************
void MacGrid::setTestValues()
{
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
   //             m_press[i][j][k] = 1.0f;
            }
        }
    }

#ifdef FLAMETHROWERS
    for(int i = PADDING_WIDTH + 3; 
            i < PADDING_WIDTH + 8; ++i)
    {
        for(int j = PADDING_WIDTH + 40; 
                j < PADDING_WIDTH + 50; ++j)
        {
            for(int k = m_nx / 2 - 10; k < m_nx / 2 + 10; ++k)
            {
                m_src[i][j] = true;
            }
        }
    }
#endif

#ifdef CIRCLE
    m_interface.addCircle((8.0f/100.0f)/(m_dx * m_scale), 
            Vector3f(m_nx/2.0f, (15/100.0f)/(m_dx * m_scale),
                m_nz/2.0f), &m_src);

/*    for(int i = PADDING_WIDTH + 24; i < m_nx - PADDING_WIDTH - 24; i += 16)
    {
        m_interface.addCircle(4.0f, Vector2f(i, 15), 
                &m_src);
    }*/
#endif

    addFuel();
/*    for(int i = 0; i < m_nx; ++i)
    {
        for(int j = 0; j < m_ny; ++j)
        {
            for(int k = 0; k < m_nz; ++k)
            {
                m_interface[i][j][k] = 1;
            }
        }
    }*/
 
#ifdef DSD
    m_interface.updateMeanCurvature(m_curvature);
    m_curvature.extrapolate(m_interface);
#endif

    /*cout << "Initialized grids" << endl << endl;
    cout << "P :" << endl << m_press;
    cout << endl << m_vel;
    cout << endl << "T :" << endl << m_temp;
    cout << endl << "S :" << endl << m_soot;
    cout << endl;*/
}

void MacGrid::step() 
{
    cout << "Adding Fuel" << endl;
    addFuel();
    cout << "Advecting Level Set" << endl;
    updateBurn();

    cout << "Advecting" << endl;
    advect();
    cout << "Applying Forces" << endl;
    applyForces();
    cout << "Projecting" << endl;
    project();
    cout << "Done" << endl;
}

float MacGrid::getDensity(int i, int j, int k)
{
    int a = i + PADDING_WIDTH;
    int b = j + PADDING_WIDTH;
    int c = k + PADDING_WIDTH;

    return m_soot[a][b][c];
}

bool MacGrid::isFuel(int i, int j, int k)
{
    int a = i + PADDING_WIDTH;
    int b = j + PADDING_WIDTH;
    int c = k + PADDING_WIDTH;
    
    return m_interface.at(a)[b][c] < 0;
}

float MacGrid::getTemp(int i, int j, int k)
{
    int a = i + PADDING_WIDTH;
    int b = j + PADDING_WIDTH;
    int c = k + PADDING_WIDTH;

    return m_temp[a][b][c];
}

void MacGrid::addSmoke()
{
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH; 
                    ++k)
            {
                if(m_src[i][j][k])
                {
                    m_temp[i][j][k] = 100.0f + AMBIENT_TEMP;
                    m_soot[i][j][k] = 0.10f;
                }
            }
        }
    }
}

void MacGrid::addFuel()
{
#ifdef CIRCLE
/*    for(int i = PADDING_WIDTH + 24; i < m_nx - PADDING_WIDTH - 24; i += 16)
    {
        m_interface.addCircle(4.0f, Vector2f(i, 15));
    }*/
    m_interface.addCircle((8.0f/100.0f)/(m_dx * m_scale), 
        Vector3f(m_nx/2.0f, (15/100.0f)/(m_dx * m_scale),
            m_nz/2.0f));
/*    
    m_interface.addCircle(4.0f, 
            Vector2f(m_nx/2.0f + 20, 15));
    m_interface.addCircle(4.0f, 
            Vector2f(m_nx/2.0f - 20, 15));*/
#endif
    
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                if(m_src[i][j][k])
                {
                    m_burn[i][j][k] = BURN_RATE;
#ifdef FLAMETHROWERS
                    m_interface.at(i)[j][k] = -1;

                    if(i < (int)(m_nx / 2.0))
                    {
                        m_vel.u_[i][j][k] = 100;
                    }
                    else
                    {
                        m_vel.u_[i][j][k] = -100;
                    }
                    m_vel.v_[i][j][k] = 0;
                    m_vel.w_[i][j][k] = 0;
#endif

#ifdef CIRCLE
                    m_vel.v_[i][j][k] = (30.0f/100.0f) 
                        / (m_scale * m_dx);
#endif
                }
            }
        }
    }

    m_interface.redistance();
}

//******************* Advection **********************

static void clampSoot(QGrid &soot, int nx, int ny, int nz)
{
    for(int i = 0; i < nx; ++i)
    {
        for(int j = 0; j < ny; ++j)
        {
            for(int k = 0; k < nz; ++k)
            {
                if(soot[i][j][k] < 0.0f)
                {
                    soot[i][j][k] = 0.0f;
                }
                else if(soot[i][j][k] > 1.0f)
                {
                    soot[i][j][k] = 1.0f;
                }
            }
        }
    }
}

void MacGrid::advect() 
{
    // Note: some advection happens in updateBurn()
    m_temp.advect(m_vel, m_dt, &m_interface);

    m_soot.advect(m_vel, m_dt, &m_interface);
    clampSoot(m_soot, m_nx, m_ny, m_nz);

    m_vel.advect(m_dt);
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
            for(int k = PADDING_WIDTH;
                    k < m_nz - PADDING_WIDTH; ++k)
            {
                float soot = (m_soot[i][j-1][k]
                        + m_soot[i][j][k]) / 2.0f;
                float temp = (m_temp[i][j-1][k]
                        + m_temp[i][j][k]) / 2.0f;
                float buoy = (alpha * soot - beta 
                        * (temp - AMBIENT_TEMP)) * GRAVITY;

                m_vel.v_[i][j][k] += m_dt * buoy;
            }
        }
    }
}

// Algorithm from p.167 of Bridson
void MacGrid::applyVorticity()
{
    Grid omegaX(m_nx - PADDING_WIDTH, 
            0.0f, m_dx);
    Grid omegaY(m_nx - PADDING_WIDTH, 
            0.0f, m_dx);
    Grid omegaZ(m_nx - PADDING_WIDTH, 
             0.0f, m_dx);
    
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = 0; k < PADDING_WIDTH; ++k)
            {
                Vector3f omega = getOmega(i, j, k);
                omegaX[i][j][k] = omega[0];
                omegaY[i][j][k] = omega[1];
                omegaZ[i][j][k] = omega[2];
            }
        }
    }

    QGrid vortX(m_nx - PADDING_WIDTH, 
            0.0f, m_dx);
    QGrid vortY(m_nx - PADDING_WIDTH, 
            0.0f, m_dx);
    QGrid vortZ(m_nx - PADDING_WIDTH, 
            0.0f, m_dx);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                Vector3f fIndex(i, j, k);
                Vector3i iIndex(i, j, k);
                Vector3f gradOmega;
                
                gradOmega[0] = omegaX.getGradient(fIndex,
                        Grid::Axis::X);
                gradOmega[1] = omegaY.getGradient(fIndex,
                        Grid::Axis::Y);
                gradOmega[2] = omegaZ.getGradient(fIndex,
                        Grid::Axis::Z);

                float gradOmegaNorm = gradOmega.norm();
                float denom = gradOmegaNorm + 10e-20 * (1.0f / 
                        (m_dx * m_dt));

                Vector3d N(gradOmega[0] / denom,
                    gradOmega[1] / denom,
                    gradOmega[2] / denom);
                
                Vector3d omega(omegaX.get(iIndex),
                    omegaY.get(iIndex),
                    omegaZ.get(iIndex));

                Vector3d f = VORTICITY_EPS * m_dx * 
                    N.cross(omega);

                vortX[i][j][k] = f[0];
                vortY[i][j][k] = f[1];
                vortZ[i][j][k] = f[2];
            }
        }
    }

    for(int i = PADDING_WIDTH; i < m_nx + 1 - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny + 1 - PADDING_WIDTH;
                ++j)
        {
            for(int k = PADDING_WIDTH; 
                    k < m_nz + 1 - PADDING_WIDTH; ++k)
            {
                if(j < m_ny - PADDING_WIDTH && 
                        k < m_nz - PADDING_WIDTH)
                {
                    m_vel.u_[i][j][k] += m_dt * 
                        vortX.getQuantity(
                            Vector3f(i-0.5f, j, k));
                }

                if(i < m_nx - PADDING_WIDTH && 
                        k < m_nz - PADDING_WIDTH)
                {
                    m_vel.v_[i][j][k] += m_dt * 
                        vortY.getQuantity(
                            Vector3f(i, j-0.5f, k));
                }

                if(i < m_nx - PADDING_WIDTH &&
                        j < m_ny - PADDING_WIDTH)
                {
                    m_vel.w_[i][j][k] += m_dt * 
                        vortZ.getQuantity(
                            Vector3f(i, j, k-0.5f));
                }
            }
        }
    }
}

Vector3f MacGrid::getOmega(int i, int j, int k)
{
    Vector3f index(i, j, k);
    
    float omegaX = m_vel.getGradient(index, Grid::Axis::Z) -
        m_vel.getGradient(index, Grid::Axis::Y);

    float omegaY = m_vel.getGradient(index, Grid::Axis::X) -
        m_vel.getGradient(index, Grid::Axis::Z);

    float omegaZ = m_vel.getGradient(index, Grid::Axis::Y) -
        m_vel.getGradient(index, Grid::Axis::X);

    return Vector3f(omegaX, omegaY, omegaZ);
}

//******************* Projection **************************

static void solveSystem(
        const SparseMatrix<double, Eigen::RowMajor> &sm, 
        const VectorXd &rhs, VectorXd &result)
{
    result = viennacl::linalg::solve(sm, rhs, viennacl::linalg::cg_tag());
}

void MacGrid::project() 
{
    int size = (m_nx-PADDING_WIDTH*2)*(m_ny-PADDING_WIDTH*2)*
        (m_nz-PADDING_WIDTH*2);
    
    SparseMatrix<double, Eigen::RowMajor> sm(size, size);
    sm.reserve(7 * size);
    
    VectorXd rhs(size);
    buildSystem(sm, rhs);

    VectorXd result(size);
    solveSystem(sm, rhs, result);

    updatePressure(result);
    applyPressure();
}

// Algorithm from p.78 of Bridson, everywhere is fluid
void MacGrid::addTriplets(vector<Trip> &trip)
{
    int width = m_nx - PADDING_WIDTH*2;
    int height = m_ny - PADDING_WIDTH*2;
    int depth = m_nz - PADDING_WIDTH*2;
    
    float scale = m_dt / (AIR_DENSITY * m_dx * m_dx);

    for(int i = 0; i < width; ++i)
    {
        for(int j = 0; j < height; ++j)
        {
            for(int k = 0; k < depth; ++k)
            {
                int index = i*height*depth + j*depth + k;
                trip.push_back(Trip(index, index, 6.0*scale));

                // X-axis values
                if(i + 1 < width)
                {
                    trip.push_back(Trip(index, 
                                index + height*depth, -scale));
                    trip.push_back(Trip(index + height*depth, 
                                index, -scale));
                }

                // Y-axis values
                if(j + 1 < height)
                {
                    trip.push_back(Trip(index, index + depth,
                                -scale));
                    trip.push_back(Trip(index + depth, index,
                                -scale));
                }

                // Z-axis values
                if(k + 1 < depth)
                {
                    trip.push_back(Trip(index, 
                                index + 1, -scale));
                    trip.push_back(Trip(index + 1 ,index, 
                                -scale));
                }
            }
        }
    }
}

void MacGrid::buildSystem(
        SparseMatrix<double, Eigen::RowMajor> &sm, 
        VectorXd &rhs)
{
    int width = m_nx - PADDING_WIDTH*2;
    int height = m_ny - PADDING_WIDTH*2;
    int depth = m_nz - PADDING_WIDTH*2;

    float scale = 1.0f / m_dx;
   
    vector<Trip> trip;
    trip.reserve(7 * width * height * depth);

    for(int i = 0; i < width; ++i)
    {
        int a = i + PADDING_WIDTH;

        for(int j = 0; j < height; ++j)
        {
            int b = j + PADDING_WIDTH;

            for(int k = 0; k < depth; ++k)
            {
                int c = k + PADDING_WIDTH;
                int index = i * depth * height + j * depth + k; 

                rhs[index] = -scale 
                    * getCorrectedRhs(Vector3i(a, b, c));
            }
        }
    }

    addTriplets(trip);
    sm.setFromTriplets(trip.begin(), trip.end());
}

// Uses divergence calculation from p.72 of Bridson and jump
// conditions from Nguyen 2002 to get RHS
float MacGrid::getCorrectedRhs(Vector3i pos)
{
    int a = pos[0];
    int b = pos[1];
    int c = pos[2];

    float dv = (FUEL_DENSITY 
            / FLAME_DENSITY - 1.0f) * m_burn[a][b][c];

    bool isFuel = m_interface.at(a)[b][c] < 0;

    float rhsVal = 0.0f;

    for(int ax = 0; ax < DIM; ++ax)
    {
        for(int side = 0; side <= 1; ++side)
        {
            Vector3i adjacentIndex = pos;
            adjacentIndex[ax] += side;
            
            float termVal = 0.0f;

            switch(ax)
            {
                case Grid::Axis::X:
                    termVal = m_vel.u_.get(adjacentIndex);
                    break;

                case Grid::Axis::Y:
                    termVal = m_vel.v_.get(adjacentIndex);
                    break;

                case Grid::Axis::Z:
                    termVal = m_vel.w_.get(adjacentIndex);
                    break;
            }

            Vector3f adjacentPos(a, b, c);
            adjacentPos[ax] += (float)side - 0.5f;
            
            bool isAdjacentFuel = 
                m_interface.lerp(adjacentPos) < 0;

            if(isFuel && !isAdjacentFuel)
            {
                termVal -= dv * 
                    m_interface.getGradient(adjacentPos, 
                            (Grid::Axis)ax);
            }
            else if(!isFuel && isAdjacentFuel)
            {
                termVal += dv *
                    m_interface.getGradient(adjacentPos, 
                            (Grid::Axis)ax);
            }

            float sign = side ? 1.0f : -1.0f;
            rhsVal += sign * termVal;
        }
    }

    return rhsVal;
}

void MacGrid::updatePressure(const VectorXd &result)
{
    int width = m_nx - PADDING_WIDTH*2;
    int height = m_ny - PADDING_WIDTH*2;
    int depth = m_nz - PADDING_WIDTH*2;

    for(int i = 0; i < width; ++i)
    {
        int a = i + PADDING_WIDTH;

        for(int j = 0; j < height; ++j)
        {
            int b = j + PADDING_WIDTH;

            for(int k = 0; k < depth; ++k)
            {
                int c = k + PADDING_WIDTH;
                int index = i * depth * height + j * depth + k; 

                m_press[a][b][c] = result[index];
            }
        }
    }
}


// Algorithm from p.71 of Bridson, ignoring possibility of solid
// cells
void MacGrid::applyPressure()
{
    float scale = m_dt / (AIR_DENSITY * m_dx);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH + 1; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH + 1; 
                ++j)
        {
            for(int k = PADDING_WIDTH; 
                    k < m_nz - PADDING_WIDTH + 1; ++k)
            {
                m_vel.u_[i][j][k] -= scale * (m_press[i][j][k]
                        - m_press[i-1][j][k]);
                m_vel.v_[i][j][k] -= scale * (m_press[i][j][k]
                        - m_press[i][j-1][k]);
                m_vel.w_[i][j][k] -= scale * (m_press[i][j][k]
                        - m_press[i][j][k-1]);
            }
        }
    }
}

// Update_D_With_DSD from Hong 2007
void MacGrid::updateBurn()
{
    m_interface.advect(m_dt, m_vel, m_burn);
    m_interface.redistance();

#ifdef DSD
    m_burn.advect(m_vel, m_dt, &m_interface);
    m_dburn.advect(m_vel, m_dt, &m_interface);
    m_curvature.advect(m_vel, m_dt, &m_interface);

    DGrid oldCurvature(m_nx - PADDING_WIDTH*2, 
            0.0f, m_dx);

    std::swap(m_curvature, oldCurvature);
    m_interface.updateMeanCurvature(m_curvature);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                if(m_interface.inFuelRegion(i, j, k))
                {
                    float dcurv = (m_curvature[i][j][k]
                            - oldCurvature[i][j][k])/m_dt;

                    float delta = m_burn[i][j][k] - DCJ;
                    float alpha = exp(MU_THETA * delta);
                    float lcj = log(fabs(1 + C_5_THETA * 
                                m_curvature[i][j][k] / alpha));

                    float ddburn = -C_1 * alpha * alpha * delta 
                        - C_2 * alpha * m_dburn[i][j][k]
                        - C_3 * alpha * alpha * lcj 
                        - C_4 * dcurv;

                    m_dburn[i][j][k] += ddburn * m_dt;
                    m_burn[i][j][k] += m_dburn[i][j][k] * m_dt;
                }
            }
        }
    }

    m_burn.extrapolate(m_interface);
    m_dburn.extrapolate(m_interface);
    m_curvature.extrapolate(m_interface);
#endif
}
