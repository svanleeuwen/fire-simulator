#define SWEEP_COUNT 2

#include "level_set.hpp"
#include "mac_grid.hpp"
#include "v_grid.hpp"
#include "d_grid.hpp"

#include <cfloat>
#include <iostream>
#include <stdlib.h>

using std::cerr;
using std::endl;
using std::cout;

LevelSet::LevelSet(int nx, int ny) :
    Grid(nx, ny, FLT_MAX)
{
    m_dx = 1.0f / m_nx;
}

void LevelSet::redistance()
{
    redistanceAdjacent();

    for(int k = 0; k < SWEEP_COUNT; ++k)
    {
        for(int i = 0; i < DIM; ++i)
        {
            for(int j = 0; j < 2; ++j)
            {
                sweep((bool)j, (Axis)i);
            }
        }
    }
}

void LevelSet::addCircle(float radius,
        Vector2f center,
        vector< vector<bool> > *src)
{
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {

            float dist = sqrt(pow((i-center[0]),2) 
                        + pow((j-center[1]),2)) - radius;

            if(dist < 0)
            {
                at(i)[j] = dist;

/*                Vector2f dir = {i-center[0], j-center[1]};
                if(dir[0] != 0 || dir[1] != 0)
                    dir.normalize();

                dir *= 30;

                vel.u_[i][j] = dir[0];
                vel.v_[i][j] = dir[1];
*/
                if(src != NULL)
                {
                    src->at(i)[j] = true;
                }
            }
        }
    }
}

void LevelSet::advect(float dt, VGrid &velocities, DGrid &burn)
{
    int nx = m_nx - PADDING_WIDTH * 2;
    int ny = m_ny - PADDING_WIDTH * 2;

    LevelSet update(nx, ny);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            Vector2f x_g(i, j);
            Vector2f vel = velocities.getVelocity(x_g);

            Vector2f x_p = velocities.rk2(x_g, vel, dt);
            update[i][j] = lerp(x_p) + dt * burn.lerp(x_p);
            
//            Vector2f n = getGradient(Vector2f(i, j));
//            Vector2f uf = velocities.getVelocity(Vector2f(i, j));

//            Vector2f w = uf + BURN_RATE * n;
//            Vector2f phi = getUpwind(w, Vector2i(i, j));
            
//            update[i][j] = at(i)[j] + dt * w.dot(phi);
        }
    }

    *this = std::move(update);
}

Vector2f LevelSet::getGradient(Vector2f x)
{
    return Vector2f(getGradientX(x), getGradientY(x));
}

float LevelSet::getGradientX(Vector2f x)
{
    int i = (int)(x[0] - 0.5f);
    int j = (int)x[1];

    float bl = (at(i+1)[j] - at(i)[j])/m_dx;
    float br = (at(i+2)[j] - at(i+1)[j])/m_dx;
    float tl = (at(i+1)[j+1] - at(i)[j+1])/m_dx;
    float tr = (at(i+2)[j+1] - at(i+1)[j+1])/m_dx;

    float s_x = x[0] - (i + 0.5f);
    float s_y = x[1] - j;

    return (1.0f - s_x) * (1.0f - s_y) * bl 
        + s_x * (1.0f - s_y) * br 
        + (1.0f - s_x) * s_y * tl 
        + s_x * s_y * tr;
}

float LevelSet::getGradientY(Vector2f x)
{
    int i = (int)x[0];
    int j = (int)(x[1] - 0.5f);

    float bl = (at(i)[j+1] - at(i)[j])/m_dx;
    float br = (at(i+1)[j+1] - at(i+1)[j])/m_dx;
    float tl = (at(i)[j+2] - at(i)[j+1])/m_dx;
    float tr = (at(i+1)[j+2] - at(i+1)[j+1])/m_dx;

    float s_x = x[0] - i;
    float s_y = x[1] - (j + 0.5f);
    
    return (1.0f - s_x) * (1.0f - s_y) * bl 
        + s_x * (1.0f - s_y) * br 
        + (1.0f - s_x) * s_y * tl 
        + s_x * s_y * tr;

}

bool LevelSet::inFuelRegion(int i, int j)
{
    return (at(i)[j] <= 0 
            || at(i-1)[j] < 0 
            || at(i+1)[j] < 0 
            || at(i)[j-1] < 0 
            || at(i)[j+1] < 0);
}

// This is suboptimal
void LevelSet::redistanceAdjacent()
{
    LevelSet closest(m_nx - PADDING_WIDTH*2, 
            m_ny - PADDING_WIDTH*2);

    for(int i = 0; i < m_nx; ++i)
    {
        for(int j = 0; j < m_ny; ++j)
        {
            float dist = FLT_MAX;

            float val = at(i)[j];
            float sign = (val > 0.0f) ? 1.0f : -1.0f;
            
            if(i > 0)
            {
                if(sign * at(i-1)[j] < 0)
                {
                    float theta = val / (val - at(i-1)[j]);
                    dist = fmin(dist, theta * m_dx);
                }
            }

            if(j > 0)
            {
                if(sign * at(i)[j-1] < 0)
                {
                    float theta = val / (val - at(i)[j-1]);
                    dist = fmin(dist, theta * m_dx);
                }
            }

            if(i < m_nx - 1)
            {
                if(sign * at(i+1)[j] < 0)
                {
                    float theta = val / (val - at(i+1)[j]);
                    dist = fmin(dist, theta * m_dx);
                }
            }

            if(j < m_ny - 1)
            {
                if(sign * at(i)[j+1] < 0)
                {
                    float theta = val / (val - at(i)[j+1]);
                    dist = fmin(dist, theta * m_dx);
                }
            }

            closest[i][j] = sign * dist;
        }
    }

    std::swap(*this, closest);
}

static Vector2f sortDistances(Vector2f phi)
{
    if(phi[1] < phi[0])
    {
        return Vector2f(phi[1], phi[0]);
    }
    else
    {
        return phi;
    }
}

void LevelSet::sweep(bool increasing, Axis ax)
{
    for(int a = 0; a < m_nx; ++a)
    {
        for(int b = 0; b < m_ny; ++b)
        {
            int i, j;
            switch(ax)
            {
                case X:
                    i = (ax == X && !increasing) ? 
                        m_nx - 1 - a : a;
                    j = (ax == Y && !increasing) ? 
                        m_ny - 1 - b : b;
                    break;

                case Y:
                    j = (ax == X && !increasing) ? 
                        m_nx - 1 - a : a;
                    i = (ax == Y && !increasing) ? 
                        m_ny - 1 - b : b;
                    break;

                default:
                    cerr << "Unimplemented axis" << endl;
                    exit(1);
            }

            Vector2f phi;
            float sign = (at(i)[j] > 0.0f) ? 1.0f : -1.0f;

            if(i == 0)
            {
                phi[0] = fabs(at(i+1)[j]);
            }
            else if(i == m_nx - 1)
            {
                phi[0] = fabs(at(i-1)[j]);
            }
            else
            {
                phi[0] = fmin(fabs(at(i-1)[j]), 
                        fabs(at(i+1)[j]));
            }

            if(j == 0)
            {
                phi[1] = fabs(at(i)[j+1]);
            }
            else if(j == m_nx - 1)
            {
                phi[1] = fabs(at(i)[j-1]);
            }
            else
            {
                phi[1] = fmin(fabs(at(i)[j-1]), 
                        fabs(at(i)[j+1]));
            }

            if(phi[0] == FLT_MAX && phi[1] == FLT_MAX)
            {
                continue;
            }

            phi = sortDistances(phi);
            at(i)[j] = sign 
                * fmin(fabs(at(i)[j]),updateDistance(phi));
        }
    }
}

// Figure 4.4 from p56 of Bridson
float LevelSet::updateDistance(Vector2f phi)
{
    float d = phi[0] + m_dx;

    if(d > phi[1])
    {
        d = 0.5f * (phi[0] + phi[1] + 
                sqrt(2.0f * m_dx * m_dx 
                    - pow(phi[1] - phi[0], 2.0f)));
    }

    return d;
}

Vector2f LevelSet::getUpwind(Vector2f w, Vector2i pos)
{
    Vector2f phi;

    int i = pos[0];
    int j = pos[1];

    if(w[0] > 0)
    {
        phi[0] = (at(i)[j] - at(i-1)[j])/m_dx;
    }
    else
    {
        phi[0] = (at(i+1)[j] - at(i)[j])/m_dx;
    }

    if(w[1] > 0)
    {
        phi[0] = (at(i)[j] - at(i)[j-1])/m_dx;
    }
    else
    {
        phi[0] = (at(i)[j+1] - at(i)[j])/m_dx;
    }

    return phi;
}
