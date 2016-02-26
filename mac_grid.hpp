#ifndef MAC_GRID_HPP
#define MAC_GRID_HPP

#include "q_grid.hpp"
#include "v_grid.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

using Eigen::Vector2i;
using Eigen::Vector2f;
using Eigen::VectorXf;
using Eigen::SparseMatrix;
using Eigen::Triplet;

class MacGrid {

public:
    MacGrid(int res, float dt);

    void step();
    float getDensity(int i, int j);

    void testSolver();

private:
    typedef vector< vector<float> > Mat;
    typedef Triplet<float> Trip;

    enum Material { Fluid, Solid, Empty };
    enum Quantity { Temp, Soot };

//************** Testing ******************
    void setTestValues();

//**************** Stuff ****************
    void addSmoke();
    void advect();
    void applyForces();

//*************** Projection *****************
    void project();
    
    void addTriplets(const Mat &A_diag, const Mat &A_x, 
            const Mat &A_y, vector<Trip> &trip);
    void addEntries(Mat &A_diag, Mat &A_x, Mat &A_y, 
            const Vector2i &g_point);
    void buildSystem(SparseMatrix<float, Eigen::RowMajor> &sm,
            VectorXf &rhs);
    void solveSystem(
            const SparseMatrix<float, Eigen::RowMajor> &sm,
            const VectorXf &rhs, VectorXf &result);
    void updatePressure(const VectorXf &result);
    void applyPressure();

//**************** Variables ******************
    int m_nx;
    int m_ny;
    float m_dt;
    float m_dx;
    
    VGrid m_vel;
    
    Grid m_press;
    
    QGrid m_temp;
    QGrid m_soot;
    
    vector< vector<Material> > m_mat;
    vector< vector<bool> > m_src;
};

#endif
