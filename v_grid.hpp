#ifndef V_GRID_HPP
#define V_GRID_HPP

#include "grid.hpp"
#include "level_set.hpp"

using Eigen::Vector2i;

class VGrid
{
public:
    VGrid(int nx, int ny);

    // Must advect quantities before velocities
    void advect(float dt, LevelSet *ls = NULL);
   
    Vector2f getVelocity(const Vector2f &x);
    Vector2f rk2(const Vector2f &x_g, const Vector2f &vel, 
            float dt);

    Grid u_;
    Grid v_;

protected:
    Vector2f getVelocity(const Vector2f &x, LevelSet *ls);
    Vector2f rk2(const Vector2f &x_g, const Vector2f &vel, 
            float dt, LevelSet *ls);

    Vector2f rk3(const Vector2f &x_g, const Vector2f &vel, 
            float dt, LevelSet *ls);

private:
    void advect(float dt, Grid::Axis ax, Grid &ret, 
            LevelSet *ls = NULL);
    float applyJumpConditions(Vector2i x, Grid::Axis ax,
            float val, LevelSet *ls);

    float getValue(Grid &g, const Vector2i &x, const Vector2f &s,
            Grid::Axis ax, LevelSet *ls);
    
    float getVelocityU(const Vector2f &x, LevelSet *ls = NULL);
    float getVelocityV(const Vector2f &x, LevelSet *ls = NULL);

    Vector2f clamp_pos(const Vector2f &pos);

    friend std::ostream& operator<<(std::ostream &out, VGrid &grid);

    int m_nx;
    int m_ny;

    bool advectingFlame;
};
#endif
