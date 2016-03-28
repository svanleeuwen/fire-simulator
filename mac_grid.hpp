#ifndef MAC_GRID_HPP
#define MAC_GRID_HPP

// From Nguyen 2002
#define FLAME_DENSITY 0.05f // kg/m^3
#define FUEL_DENSITY 1.0f // kg/m^3
#define BURN_RATE 0.1f // m/s

#include "s_grid.hpp"
#include "t_grid.hpp"
#include "v_grid.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

using Eigen::Vector2i;
using Eigen::Vector2f;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;

class MacGrid {

public:
    MacGrid(int res, float dt, bool fire = false);

    void step();
    float getDensity(int i, int j);
    bool isFuel(int i, int j);
    float getTemp(int i, int j);

    void testSolver();

private:
    typedef vector< vector<float> > Mat;
    typedef Triplet<float> Trip;

    enum Material { Fire, Solid, Empty };

//************** Testing ******************
    void setTestValues();

//**************** Stuff ****************
    void addSmoke();
    void addFuel();
    void advect();

    void applyForces();
    void applyBuoyancy();
    void applyVorticity();
   
    float getOmega(int i, int j);
    Vector2f getGradOmega(Grid &omega, int i, int j);
    Vector2f getVorticity(int i, int j);

//*************** Projection *****************
    void project();
    
    void addTriplets(const Mat &A_diag, const Mat &A_x, 
            const Mat &A_y, vector<Trip> &trip);
    void addEntries(Mat &A_diag, Mat &A_x, Mat &A_y, 
            const Vector2i &g_point);
    void buildSystem(SparseMatrix<double, Eigen::RowMajor> &sm,
            VectorXd &rhs);

    float getCorrectedRhs(Vector2i pos);
    void solveSystem(
            const SparseMatrix<double, Eigen::RowMajor> &sm,
            const VectorXd &rhs, VectorXd &result);
    void updatePressure(const VectorXd &result);
    void applyPressure();

//**************** Variables ******************
    int m_nx;
    int m_ny;
    float m_dt;
    float m_dx;
    
    VGrid m_vel;
    Grid m_press;
    
    TGrid m_temp;
    SGrid m_soot;

    LevelSet *m_interface;

    vector< vector<Material> > m_mat;
    vector< vector<bool> > m_src;
};

#endif
