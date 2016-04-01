#ifndef V_GRID_HPP
#define V_GRID_HPP

#include "grid.hpp"
#include "level_set.hpp"
#include "d_grid.hpp"

using Eigen::Vector2i;

class VGrid
{
public:
    VGrid(int nx, int ny, LevelSet *interface, DGrid *burn);

    // Must advect quantities before velocities
    void advect(float dt);
   
    Vector2f getVelocity(const Vector2f &x, 
            bool advecting = false);
    Vector2f rk2(const Vector2f &x_g, const Vector2f &vel, 
            float dt, bool advecting = false);

    Grid u_;
    Grid v_;

protected:
    Vector2f rk3(const Vector2f &x_g, const Vector2f &vel, 
            float dt, bool advecting);

private:
    void advect(float dt, Grid::Axis ax, Grid &ret);
    float applyJumpConditions(Vector2i x, Grid::Axis ax,
            float val);

    float getValue(Grid &g, const Vector2i &x, const Vector2f &s,
            Grid::Axis ax, bool advecting);
    
    float getVelocityU(const Vector2f &x, bool advecting = false);
    float getVelocityV(const Vector2f &x, bool advecting = false);

    Vector2f clamp_pos(const Vector2f &pos);

    friend std::ostream& operator<<(std::ostream &out, VGrid &grid);

    int m_nx;
    int m_ny;

    bool advectingFlame;

    LevelSet *m_interface;
    DGrid *m_burn;
};
#endif
