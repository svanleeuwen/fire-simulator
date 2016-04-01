#define OUTER_VELOCITY 0.0f

#include "v_grid.hpp"
#include "mac_grid.hpp"

#include <stdlib.h>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;
using Eigen::Vector4f;

VGrid::VGrid(int nx, int ny, LevelSet *interface, DGrid *burn) :
    u_(nx+1, ny, OUTER_VELOCITY), 
    v_(nx, ny+1, OUTER_VELOCITY),
    m_interface(interface),
    m_burn(burn)
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

    advect(dt, Grid::Axis::X, u_copy);
    advect(dt, Grid::Axis::Y, v_copy);

    u_ = std::move(u_copy);
    v_ = std::move(v_copy);
}

Vector2f VGrid::getVelocity(const Vector2f &x, bool advecting)
{
    return Vector2f(getVelocityU(x, advecting), 
            getVelocityV(x, advecting));
}

// Runge-Kutta 2 from Bridson p.32
Vector2f VGrid::rk2(const Vector2f &x_g, const Vector2f &vel,
        float dt, bool advecting) 
{
/*    Vector2f x_mid = clamp_pos(x_g - 0.5f * dt * vel);
    return clamp_pos(x_g - dt * getVelocity(x_mid, ls));*/

    return rk3(x_g, vel, dt, advecting);
}

Vector2f VGrid::rk3(const Vector2f &x_g, const Vector2f &vel,
        float dt, bool advecting) 
{
    Vector2f k1 = vel;
    Vector2f k2 = getVelocity(clamp_pos(
                x_g - 0.5f * dt * k1), advecting);
    Vector2f k3 = getVelocity(clamp_pos(
                x_g - 0.75f * dt * k2), advecting);

    Vector2f ret = x_g - dt * (1.0f/9.0f) * 
        (2.0f * k1 + 3.0f * k2 + 4.0f * k3);

    return clamp_pos(ret);
}


void VGrid::advect(float dt, Grid::Axis ax, Grid &ret)
{
    int width = m_nx + (ax == Grid::Axis::X);
    int height = m_ny + (ax == Grid::Axis::Y);

    for(int i = PADDING_WIDTH; i < width - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < height - PADDING_WIDTH; ++j)
        {
            Vector2f x_g(i, j);
            x_g[ax] -= 0.5f;

            advectingFlame = m_interface->lerp(x_g) > 0;
                       
            Vector2f vel = getVelocity(x_g, true);
            Vector2f x_p = rk2(x_g, vel, dt, true);

            switch(ax)
            {
                case Grid::Axis::X:
                    ret[i][j] = getVelocityU(x_p, true);
                    
                    break;

                case Grid::Axis::Y:
                    ret[i][j] = getVelocityV(x_p, true);
                    break;
               
                default:
                    cerr << "Unimplemented axis" << endl;
                    exit(1);
            }
        }
    }
}

// Jump condition p135 of Bridson
float VGrid::applyJumpConditions(Vector2i x, Grid::Axis ax,
            float val)
{
    int i = x[0];
    int j = x[1];

    float dv = (FUEL_DENSITY / FLAME_DENSITY - 1.0f)
        * m_burn->at(i)[j];
    float n;

    // Flame advecting from fuel
    if(m_interface->at(i)[j] > 0 && !advectingFlame)
    {
        switch(ax)
        {
            case Grid::Axis::X:
                n = m_interface->getGradientX(Vector2f(i, j));
                break;
            
            case Grid::Axis::Y:
                n = m_interface->getGradientY(Vector2f(i, j));
                break;
            
            default:
                cerr << "Unimplemented axis" << endl;
                exit(1);
        }

        return val + dv * n;
    }

    // Fuel advecting from flame
    else if(m_interface->at(i)[j] < 0 && advectingFlame)
    {
        switch(ax)
        {
            case Grid::Axis::X:
                n = m_interface->getGradientX(Vector2f(i, j));
                break;
            
            case Grid::Axis::Y:
                n = m_interface->getGradientY(Vector2f(i, j));
                break;
            
            default:
                cerr << "Unimplemented axis" << endl;
                exit(1);
        }
        
        return val - dv * n;
    }

    return val;
}

// Monotonic cubic interpolation from Fedkiw 2001 for smoke
float VGrid::getValue(Grid &g, const Vector2i &x, 
        const Vector2f &s, Grid::Axis ax, bool advecting) 
{
    if(!advecting)
    {
        return g.getValue(x, s);
    }

    int i = x[0];
    int j = x[1];

    Vector4f q = {0.0f, 0.0f, 0.0f, 0.0f};

    for(int a = 0; a < 4; ++a) 
    {
        if(s[1] < EPSILON) 
        {
            if(s[0] < EPSILON)
            {
                return applyJumpConditions(Vector2i(i+1, j+1),
                        ax, g.at(i+1)[j+1]);
            }

            q[a] = applyJumpConditions(Vector2i(i+a,j+1),
                    ax, g.at(i+a)[j+1]);
        }
        else
        {
            if(s[0] < EPSILON)
            {
                a = 1;
            }

            Vector4f p;
            for(int c = 0; c < 4; ++c)
            {
                p[c] = applyJumpConditions(Vector2i(i+a,j+c),
                        ax, g.at(i+a)[j+c]);
            }

            Vector4f w_y = g.getW(p);
            
            for(int b = 0; b < 4; ++b) 
            {
                q[a] += w_y[b] * pow(s[1], b);
            }

            if(s[0] < EPSILON)
            {
                break;
            }
        }
    }

    float val = 0.0f;
    Vector4f w_x = g.getW(q);

    for(int a = 0; a < 4; ++a) 
    {
        val += w_x[a] * pow(s[0], a);
    }

    return val;
}


float VGrid::getVelocityU(const Vector2f &x, bool advecting)
{
    int i = (int)(x[0] + 0.5f) - 1;
    int j = (int)(x[1]) - 1;
    float s_x = x[0] - i - 0.5f;
    float s_y = x[1] - j - 1.0f;

    return getValue(u_, Vector2i(i, j), Vector2f(s_x, s_y),
            Grid::Axis::X, advecting);
}

float VGrid::getVelocityV(const Vector2f &x, bool advecting)
{
    int i = (int)(x[0]) - 1;
    int j = (int)(x[1] + 0.5f) - 1;
    float s_x = x[0] - i - 1.0f;
    float s_y = x[1] - j - 0.5f;

    return getValue(v_, Vector2i(i, j), Vector2f(s_x, s_y),
            Grid::Axis::Y, advecting);
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
