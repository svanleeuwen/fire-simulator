#ifndef V_GRID_HPP
#define V_GRID_HPP

#include "grid.hpp"

using Eigen::Vector2i;

class VGrid
{
public:
    VGrid(int nx, int ny);

    // Must advect quantities before velocities
    void advect(float dt);
    
    // Gets the velocity at an arbitrary point x
    Vector2f getVelocity(const Vector2f &x);
    Vector2f rk2(const Vector2f &x_g, const Vector2f &vel, 
            float dt);

    Grid u_;
    Grid v_;

private:
    enum Axis { U, V };
    
    void advect(float dt, Axis ax, Grid &ret);
    
    float getVelocityU(const Vector2f &x);
    float getVelocityV(const Vector2f &x);

    Vector2f clamp_pos(const Vector2f &pos);

    friend std::ostream& operator<<(std::ostream &out, VGrid &grid);

    int m_nx;
    int m_ny;
};
#endif
