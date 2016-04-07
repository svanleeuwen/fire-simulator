#ifndef MAC_GRID_HPP
#define MAC_GRID_HPP

// From Nguyen 2002
#define FLAME_DENSITY 0.01f // kg/m^3
#define FUEL_DENSITY 1.0f // kg/m^3
#define BURN_RATE 0.2f // m/s

#include "s_grid.hpp"
#include "t_grid.hpp"
#include "v_grid.hpp"
#include "d_grid.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

using Eigen::Vector2i;
using Eigen::Vector2f;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;

class MacGrid {

public:
    MacGrid(int res, float dt, float scale);

    void step();
    float getDensity(int i, int j, int k);
    bool isFuel(int i, int j, int k);
    float getTemp(int i, int j, int k);

private:
    typedef Triplet<float> Trip;

    enum Material { Fire, Solid, Empty };

//************** Testing ******************
    void setTestValues();
    void addSmoke();
    void addFuel();
    
//************** Advection *****************
    void advect();

//************ Update Forces ****************

    void applyForces();
    void applyBuoyancy();
    void applyVorticity();
   
    Vector3f getOmega(int i, int j, int k);

//*************** Projection *****************
    void project();
    
    void addTriplets(vector<Trip> &trip);
    void buildSystem(SparseMatrix<double, Eigen::RowMajor> &sm,
            VectorXd &rhs);

    float getCorrectedRhs(Vector3i pos);
    void updatePressure(const VectorXd &result);
    void applyPressure();

//****************** DSD ********************
    void updateBurn();

//**************** Variables ******************
    int m_nx;
    int m_ny;
    int m_nz;

    float m_dt;
    float m_scale;
    float m_dx;
    
    LevelSet m_interface;
    DGrid m_curvature;
    
    DGrid m_burn;
    DGrid m_dburn;
    
    VGrid m_vel;
    Grid m_press;
    
    TGrid m_temp;
    SGrid m_soot;

    vector< vector< vector<bool> > > m_src;
};

#endif
