#ifndef V_GRID_HPP
#define V_GRID_HPP

#include "grid.hpp"
#include "level_set.hpp"
#include "d_grid.hpp"

using Eigen::Vector2i;

class VGrid
{
public:
    VGrid(int nx, int ny, int nz, float dx,
            LevelSet *interface, DGrid *burn);

    // Must advect quantities before velocities
    void advect(float dt);
   
    Vector3f getVelocity(const Vector3f &x, 
            bool advecting = false);

    Vector3f rk3(const Vector3f &x_g, const Vector3f &vel, 
            float dt, bool advecting = false);

    float getGradient(Vector3f x, Grid::Axis ax);

    Grid u_;
    Grid v_;
    Grid w_;

private:
    void advect(float dt, Grid::Axis ax, Grid &ret);
    void advectHandler(Vector2i bounds, float dt, Grid::Axis ax,
            Grid &ret);
    float applyJumpConditions(Vector3i index, Grid::Axis ax,
            float val);

    float lerp(Grid &g, Vector3f x, Grid::Axis ax);
   
    float getValue(Grid &g, const Vector3i &index, 
            const Vector3f &s, Grid::Axis ax, bool advecting);
    float getValue2D(Grid &g, int i, const Vector2i &index, 
        const Vector2f &s, Grid::Axis ax);
    
    float getVelocity(const Vector3f &x, Grid::Axis ax, 
            bool advecting = false);

    Vector3f clamp_pos(const Vector3f &pos);

    friend std::ostream& operator<<(std::ostream &out, VGrid &grid);

    int m_nx;
    int m_ny;
    int m_nz;

    float m_dx;
    bool m_advectingFlame;

    LevelSet *m_interface;
    DGrid *m_burn;
};
#endif
