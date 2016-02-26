#define OUTER_VELOCITY 0.0f

#include "v_grid.hpp"

#include <stdlib.h>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;
using Eigen::Vector4f;

VGrid::VGrid(int nx, int ny) :
    u_(nx+1, ny, OUTER_VELOCITY), 
    v_(nx, ny+1, OUTER_VELOCITY) 
{
    m_nx = nx + PADDING_WIDTH * 2;
    m_ny = ny + PADDING_WIDTH * 2;
}

void VGrid::advect(float dt)
{
    int nx = m_nx - PADDING_WIDTH * 2;
    int ny = m_ny - PADDING_WIDTH * 2;

    Grid u_copy(nx+1, ny, OUTER_VELOCITY);
    Grid v_copy(nx, ny+1, OUTER_VELOCITY);

    advect(dt, Axis::U, u_copy);
    advect(dt, Axis::V, v_copy);

    u_ = std::move(u_copy);
    v_ = std::move(v_copy);
}

Vector2f VGrid::getVelocity(const Vector2f &x)
{
    return Vector2f(getVelocityU(x), getVelocityV(x));
}

// Runge-Kutta 2 from Bridson p.32
Vector2f VGrid::rk2(const Vector2f &x_g, const Vector2f &vel,
        float dt) 
{
    Vector2f x_mid = clamp_pos(x_g - 0.5f * dt * vel);
    return clamp_pos(x_g - dt * getVelocity(x_mid));
}

void VGrid::advect(float dt, Axis ax, Grid &ret)
{
    int width = m_nx + (ax == Axis::U);
    int height = m_ny + (ax == Axis::V);

    for(int i = PADDING_WIDTH; i < width - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < height - PADDING_WIDTH; ++j)
        {
            Vector2f x_g(i, j);
            x_g[ax] -= 0.5f;
           
            Vector2f vel = getVelocity(x_g);
            Vector2f x_p = rk2(x_g, vel, dt);

            switch(ax)
            {
                case U:
                    ret[i][j] = getVelocityU(x_p);
                    break;

                case V:
                    ret[i][j] = getVelocityV(x_p);
                    break;
               
                default:
                    cerr << "Unimplemented axis" << endl;
                    exit(1);
            }
        }
    }
}

float VGrid::getVelocityU(const Vector2f &x)
{
    int i = (int)(x[0] + 0.5f) - 1;
    int j = (int)(x[1]) - 1;
    float s_x = x[0] - i - 0.5f;
    float s_y = x[1] - j - 1.0f;

    return u_.getValue(Vector2i(i, j), Vector2f(s_x, s_y));
}

float VGrid::getVelocityV(const Vector2f &x)
{
    int i = (int)(x[0]) - 1;
    int j = (int)(x[1] + 0.5f) - 1;
    float s_x = x[0] - i - 1.0f;
    float s_y = x[1] - j - 0.5f;

    return v_.getValue(Vector2i(i, j), Vector2f(s_x, s_y));
}

Vector2f VGrid::clamp_pos(const Vector2f &pos)
{
    Vector2f x = pos;

    float low = PADDING_WIDTH - 0.5f;
    float x_high = m_nx - PADDING_WIDTH - 0.5f;
    float y_high = m_ny - PADDING_WIDTH - 0.5f;

    if(x[0] < low) 
    {
        x[0] = low; 
    }
    else if(x[0] > x_high)
    {
        x[0] = x_high;
    }

    if(x[1] < low) 
    {
        x[1] = low;
    }
    else if(x[1] > y_high)
    {
        x[1] = y_high;
    }

    return x;
}

std::ostream& operator<<(std::ostream &out, VGrid &grid)
{
    out << "U:" << endl << grid.u_ << endl;
    out << "V:" << endl << grid.v_ << endl;
    return out;
}
