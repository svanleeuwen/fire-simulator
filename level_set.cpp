#define SWEEP_COUNT 2

#include "level_set.hpp"
#include "mac_grid.hpp"
#include "v_grid.hpp"
#include "d_grid.hpp"

#include <cfloat>
#include <iostream>
#include <stdlib.h>
#include <thread>

using std::cerr;
using std::endl;
using std::cout;
using std::thread;

LevelSet::LevelSet(int res) :
    Grid(res, FLT_MAX, 1.0f / res)
{}

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

    resetGradients();
}

void LevelSet::addCircle(float radius,
        Vector3f center,
        vector< vector< vector<bool> > > *src)
{
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                float dist = sqrt(pow((i-center[0]),2) 
                        + pow((j-center[1]),2) 
                        + pow((k-center[2]),2)) - radius;

                if(dist < 0)
                {
                    at(i)[j][k] = dist;

                    if(src != NULL)
                    {
                        src->at(i)[j][k] = true;
                    }
                }
            }
        }
    }
}

void LevelSet::advect(float dt, VGrid &v, DGrid &burn)
{
    int res = m_nx - PADDING_WIDTH * 2;

    vector<Vector2i> bounds;
    int bWidth = res / THREAD_COUNT;

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        bounds.push_back(Vector2i(i*bWidth + PADDING_WIDTH, 
                    (i+1)*bWidth + PADDING_WIDTH));
    }
    
    bounds[THREAD_COUNT-1][1] = res + PADDING_WIDTH;

    LevelSet update(res);
    thread tt[THREAD_COUNT];

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i] = thread(&LevelSet::advectHandler, this, 
                std::ref(bounds[i]), std::ref(v), dt, 
                std::ref(burn), std::ref(update));
    }

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i].join();
    }
    
    *this = std::move(update);
    resetGradients();
}

void LevelSet::advectHandler(Vector2i bounds, VGrid &v, 
            float dt, DGrid &burn, LevelSet &update)
{
    for(int i = bounds[0]; i < bounds[1]; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                Vector3f x_g(i, j, k);
                Vector3f vel = v.getVelocity(x_g);

                Vector3f x_p = v.rk3(x_g, vel, dt);
                update[i][j][k] = 
                    lerp(x_p) + dt * burn.lerp(x_p);
            }
        }
    }
}

bool LevelSet::inFuelRegion(int i, int j, int k)
{
    Vector2i val = {-1, 1};
    bool inRegion = at(i)[j][k] <= 0;

    for(int ax = 0; ax < DIM; ++ax)
    {
        for(int b = 0; b < 2; ++b)
        {
            Vector3i index = {i, j, k};
            index[ax] += val[b];

            inRegion = inRegion || get(index) < 0;
        }
    }

    return inRegion;
}

void LevelSet::updateMeanCurvature(DGrid &curvature)
{
    for(int ax = 0; ax < DIM; ++ax)
    {
        if(m_grad[ax] == NULL)
        {
            initGradient(Grid::Axis(ax));
        }
    }

    // Adapted divergence calculation from Bridson p72
    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                Vector3f index(i, j, k);
                curvature[i][j][k] = 
                    -(m_grad[0]->getGradient(index, 
                                Grid::Axis::X)
                            + m_grad[1]->getGradient(index,
                                Grid::Axis::Y)
                            + m_grad[2]->getGradient(index,
                                Grid::Axis::Z));
            }
        }
    }

    curvature.resetGradients();
}

// This is suboptimal
void LevelSet::redistanceAdjacent()
{
    int res = m_nx - PADDING_WIDTH*2;
    LevelSet closest(res); 

    for(int i = 0; i < m_nx; ++i)
    {
        for(int j = 0; j < m_ny; ++j)
        {
            for(int k = 0; k < m_nz; ++k)
            {
                float dist = FLT_MAX;

                float val = at(i)[j][k];
                float sign = (val > 0.0f) ? 1.0f : -1.0f;
                
                if(i > 0)
                {
                    if(sign * at(i-1)[j][k] < 0)
                    {
                        float theta = val / 
                            (val - at(i-1)[j][k]);
                        dist = fmin(dist, theta * m_dx);
                    }
                }

                if(j > 0)
                {
                    if(sign * at(i)[j-1][k] < 0)
                    {
                        float theta = val / 
                            (val - at(i)[j-1][k]);
                        dist = fmin(dist, theta * m_dx);
                    }
                }

                if(k > 0)
                {
                    if(sign * at(i)[j][k-1] < 0)
                    {
                        float theta = val / 
                            (val - at(i)[j][k-1]);
                        dist = fmin(dist, theta * m_dx);
                    }
                }

                if(i < m_nx - 1)
                {
                    if(sign * at(i+1)[j][k] < 0)
                    {
                        float theta = val /
                           (val - at(i+1)[j][k]);
                        dist = fmin(dist, theta * m_dx);
                    }
                }

                if(j < m_ny - 1)
                {
                    if(sign * at(i)[j+1][k] < 0)
                    {
                        float theta = val / 
                            (val - at(i)[j+1][k]);
                        dist = fmin(dist, theta * m_dx);
                    }
                }

                if(k < m_nz - 1)
                {
                    if(sign * at(i)[j][k+1] < 0)
                    {
                        float theta = val / 
                            (val - at(i)[j][k+1]);
                        dist = fmin(dist, theta * m_dx);
                    }
                }

                closest[i][j][k] = sign * dist;
            }
        }
    }

    *this = std::move(closest);
}

static Vector3f sortDistances(Vector3f phi)
{
    bool flag = true;

    while(flag)
    {
        flag = false;

        for(int i = 0; i <=1; ++i)
        {
            if(phi[i+1] < phi[i])
            {
                flag = true;

                float temp = phi[i];
                phi[i] = phi[i+1];
                phi[i+1] = temp;
            }
        }
    }

    return phi;
}

void LevelSet::sweep(bool increasing, Axis ax)
{
    int res = m_nx;

    for(int a = 0; a < m_nx; ++a)
    {
        for(int b = 0; b < m_ny; ++b)
        {
            for(int c = 0; c < m_nz; ++c)
            {
                int i, j, k;
                switch(ax)
                {
                    case X:
                        i = increasing ? a : res - 1 - a;
                        j = b;
                        k = c;
                        break;

                    case Y:
                        i = b;
                        j = increasing ? a : res - 1 - a;
                        k = c;
                        break;

                    case Z:
                        i = b;
                        j = c;
                        k = increasing ? a : res - 1 - a;
                        break;

                    default:
                        cerr << "Unimplemented axis" << endl;
                        exit(1);
                }

                Vector3f phi;
                float sign = (at(i)[j][k] > 0.0f) ? 1.0f : -1.0f;

                Vector3i index(i, j, k);
                for(int ax_it = 0; ax_it < DIM; ++ax_it)
                {
                    Vector3i less = index;
                    less[ax_it] -= 1;

                    Vector3i more = index;
                    more[ax_it] += 1;

                    if(index[ax_it] == 0)
                    {
                        phi[ax_it] = fabs(get(more));
                    }
                    else if(index[ax_it] == res - 1)
                    {
                        phi[ax_it] = fabs(get(less));
                    }
                    else
                    {
                        phi[ax_it] = fmin(fabs(get(less)), 
                                fabs(get(more)));
                    }
                }
                
                if(phi[0] == FLT_MAX && phi[1] == FLT_MAX
                        && phi[2] == FLT_MAX)
                {
                    continue;
                }

                phi = sortDistances(phi);
                at(i)[j][k] = sign *
                    fmin(fabs(get(index)),updateDistance(phi));
            }
        }
    }
}

// Figure 4.4 from p56 of Bridson
float LevelSet::updateDistance(Vector3f phi)
{
    float d = phi[0] + m_dx;

    if(d > phi[1])
    {
        d = 0.5f * (phi[0] + phi[1] + 
                sqrt(2.0f * m_dx * m_dx 
                    - pow(phi[1] - phi[0], 2.0f)));

        if(d > phi[2])
        {
            d = (1.0f/3.0f) * (phi[0] + phi[1] + phi[2] +
                    sqrt(fmax(0.0f,
                            pow(phi[0] + phi[1] + phi[2], 2.0)
                            - 3.0f * (phi[0] * phi[0]
                                + phi[1] * phi[1]
                                + phi[2] * phi[2]
                                - m_dx * m_dx))));
        }
    }

    return d;
}
