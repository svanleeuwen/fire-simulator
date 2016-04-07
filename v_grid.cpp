#define OUTER_VELOCITY 0.0f

#include "v_grid.hpp"
#include "mac_grid.hpp"

#include <stdlib.h>
#include <iostream>
#include <thread>

using std::cout;
using std::cerr;
using std::endl;

using std::thread;
using Eigen::Vector4f;

VGrid::VGrid(int nx, int ny, int nz, float dx,
        LevelSet *interface, DGrid *burn) :
    u_(nx+1, ny, nz, OUTER_VELOCITY, dx), 
    v_(nx, ny+1, nz, OUTER_VELOCITY, dx),
    w_(nx, ny, nz+1, OUTER_VELOCITY, dx),
    m_interface(interface),
    m_burn(burn)
{
    m_nx = nx + PADDING_WIDTH * 2;
    m_ny = ny + PADDING_WIDTH * 2;
    m_nz = nz + PADDING_WIDTH * 2;

    m_dx = dx;
}

void VGrid::advect(float dt)
{
    int nx = m_nx - PADDING_WIDTH * 2;
    int ny = m_ny - PADDING_WIDTH * 2;
    int nz = m_nz - PADDING_WIDTH * 2;

    Grid u_copy(nx+1, ny, nz, OUTER_VELOCITY, m_dx);
    Grid v_copy(nx, ny+1, nz, OUTER_VELOCITY, m_dx);
    Grid w_copy(nx, ny, nz+1, OUTER_VELOCITY, m_dx);

    advect(dt, Grid::Axis::X, u_copy);
    advect(dt, Grid::Axis::Y, v_copy);
    advect(dt, Grid::Axis::Z, w_copy);

    u_ = std::move(u_copy);
    v_ = std::move(v_copy);
    w_ = std::move(w_copy);

    u_.resetGradients();
    v_.resetGradients();
    w_.resetGradients();
}

Vector3f VGrid::getVelocity(const Vector3f &x, bool advecting)
{
    return Vector3f(getVelocity(x, Grid::Axis::X, advecting), 
            getVelocity(x, Grid::Axis::Y, advecting),
            getVelocity(x, Grid::Axis::Z, advecting));
}

Vector3f VGrid::rk3(const Vector3f &x_g, const Vector3f &vel,
        float dt, bool advecting) 
{
    Vector3f k1 = vel;
    Vector3f k2 = getVelocity(clamp_pos(
                x_g - 0.5f * dt * k1), advecting);
    Vector3f k3 = getVelocity(clamp_pos(
                x_g - 0.75f * dt * k2), advecting);

    Vector3f ret = x_g - dt * (1.0f/9.0f) * 
        (2.0f * k1 + 3.0f * k2 + 4.0f * k3);

    return clamp_pos(ret);
}

float VGrid::getGradient(Vector3f x, Grid::Axis ax)
{
    x[ax] += 0.5f;

    switch(ax)
    {
        case Grid::Axis::X :
            return u_.getGradient(x, ax);
            
        case Grid::Axis::Y :
            return v_.getGradient(x, ax);
        
        case Grid::Axis::Z :
            return w_.getGradient(x, ax);

        default:
            cerr << "Unimplemented axis" << endl;
            exit(1);
    }
}

void VGrid::advect(float dt, Grid::Axis ax, Grid &ret)
{
    int nx = ax == Grid::Axis::X ? m_nx + 1 : m_nx;

    vector<Vector2i> bounds;
    int bWidth = (nx - PADDING_WIDTH*2) / THREAD_COUNT;

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        bounds.push_back(Vector2i(i*bWidth + PADDING_WIDTH, 
                    (i+1)*bWidth + PADDING_WIDTH));
    }
    
    bounds[THREAD_COUNT-1][1] = nx - PADDING_WIDTH;

    thread tt[THREAD_COUNT];

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i] = thread(&VGrid::advectHandler, this, 
                std::ref(bounds[i]),
                dt, ax, std::ref(ret));
    }

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i].join();
    }
}

void VGrid::advectHandler(Vector2i bounds, float dt, 
        Grid::Axis ax, Grid &ret)
{
    Vector3i dim(m_nx, m_ny, m_nz);
    dim[ax] += 1;

    for(int i = bounds[0]; i < bounds[1]; ++i)
    {
        for(int j = PADDING_WIDTH; j < dim[1] - PADDING_WIDTH; 
                ++j)
        {
            for(int k = PADDING_WIDTH; k < dim[2] 
                    - PADDING_WIDTH; ++k)
            {
                Vector3f x_g(i, j, k);
                x_g[ax] -= 0.5f;

                m_advectingFlame = m_interface->lerp(x_g) > 0;

                Vector3f vel = getVelocity(x_g, true);
                Vector3f x_p = rk3(x_g, vel, dt, true);

                ret[i][j][k] = getVelocity(x_p, ax, true);
            }
        }
    }
}

// Jump condition p135 of Bridson
float VGrid::applyJumpConditions(Vector3i index, Grid::Axis ax,
            float val)
{
    float dv = (FUEL_DENSITY / FLAME_DENSITY - 1.0f)
        * m_burn->get(index);
    float n;

    Vector3f findex(index[0], index[1], index[2]);

    // Flame advecting from fuel
    if(m_interface->get(index) > 0 && !m_advectingFlame)
    {
        n = m_interface->getGradient(findex, ax);
        return val + dv * n;
    }

    // Fuel advecting from flame
    else if(m_interface->get(index) < 0 && m_advectingFlame)
    {
        n = m_interface->getGradient(findex, ax);
        return val - dv * n;
    }

    return val;
}

float VGrid::lerp(Grid &g, Vector3f x, Grid::Axis ax)
{
    int i = (int)x[0];
    int j = (int)x[1];
    int k = (int)x[2];

    float s_x = x[0] - i;
    float s_y = x[1] - j;
    float s_z = x[2] - k;

    Vector2f p_x = { 1.0f - s_x, s_x };
    Vector2f p_y = { 1.0f - s_y, s_y };
    Vector2f p_z = { 1.0f - s_z, s_z };

    float val = 0.0f;

    for(int a = 0; a < 2; ++a)
    {
        for(int b = 0; b < 2; ++b)
        {
            for(int c = 0; c < 2; ++c)
            {
                val += p_x[a] * p_y[b] * p_z[c] 
                    * applyJumpConditions(Vector3i(i+a, j+b, k+c)
                            , ax, g.at(i+a)[j+b][k+c]);
            }
        }
    }

    return val;
}

// Monotonic cubic interpolation from Fedkiw 2001 for smoke
float VGrid::getValue(Grid &g, const Vector3i &index, 
            const Vector3f &s, Grid::Axis ax, bool advecting)
{
    if(!advecting)
    {
        return g.getValue(index, s);
    }
    
    int i = index[0];
    int j = index[1];
    int k = index[2];

    Vector4f q(0.0f, 0.0f, 0.0f, 0.0f);

    for(int a = 0; a < 4; ++a)
    {
        q[a] = getValue2D(g, i+a, Vector2i(j, k), 
                Vector2f(s[1], s[2]), ax);
    }

    float val = 0.0f;
    Vector4f w_x = g.getW(q);

    for(int a = 0; a < 4; ++a) 
    {
        val += w_x[a] * pow(s[0], a);
    }

    return val;
}

float VGrid::getValue2D(Grid &g, int i, const Vector2i &index, 
        const Vector2f &s, Grid::Axis ax) 
{
    int j = index[0];
    int k = index[1];

    Vector4f q = {0.0f, 0.0f, 0.0f, 0.0f};

    for(int b = 0; b < 4; ++b) 
    {
        Vector4f p;
        for(int c = 0; c < 4; ++c)
        {
            p[c] = applyJumpConditions(Vector3i(i,j+b,k+c),
                    ax, g.at(i)[j+b][k+c]);
        }

        Vector4f w_z = g.getW(p);

        for(int c = 0; c < 4; ++c) 
        {
            q[b] += w_z[c] * pow(s[1], c);
        }
    }

    float val = 0.0f;
    Vector4f w_y = g.getW(q);

    for(int a = 0; a < 4; ++a) 
    {
        val += w_y[a] * pow(s[0], a);
    }

    return val;
}

float VGrid::getVelocity(const Vector3f &x, Grid::Axis ax,
        bool advecting)
{
    Vector3f shift(0.0f, 0.0f, 0.0f);
    shift[ax] += 0.5f;

    int i = (int)(x[0] + shift[0]) - 1;
    int j = (int)(x[1] + shift[1]) - 1;
    int k = (int)(x[2] + shift[2]) - 1;

    float s_x = x[0] - i - 1.0f + shift[0];
    float s_y = x[1] - j - 1.0f + shift[1];
    float s_z = x[2] - k - 1.0f + shift[2];

    Grid *vel[3] = {&u_, &v_, &w_};

    return getValue(*vel[ax], Vector3i(i, j, k), 
            Vector3f(s_x, s_y, s_z),
            ax, advecting);
}

Vector3f VGrid::clamp_pos(const Vector3f &pos)
{
    Vector3f x = pos;

    float low = PADDING_WIDTH - 0.5f;
    Vector3f high = { m_nx - PADDING_WIDTH - 0.5f,
        m_ny - PADDING_WIDTH - 0.5f,
        m_nz - PADDING_WIDTH - 0.5f };

    for(int ax = 0; ax < DIM; ++ax)
    {
        if(x[ax] < low) 
        {
            x[ax] = low; 
        }
        else if(x[ax] > high[ax])
        {
            x[ax] = high[ax];
        }
    }

    return x;
}

std::ostream& operator<<(std::ostream &out, VGrid &grid)
{
    out << "U:" << endl << grid.u_ << endl;
    out << "V:" << endl << grid.v_ << endl;
    out << "W:" << endl << grid.w_ << endl;
    return out;
}
