#include "d_grid.hpp"
#include "v_grid.hpp"

#include <iostream>
#include <queue>

#include <climits>
#include <thread>

using std::cout;
using std::cerr;
using std::endl;
using std::queue;

using Eigen::Vector2i;
using std::thread;

DGrid::DGrid(int res, float ambient, float dx)
    : QGrid(res, ambient, dx)
{}

void DGrid::advect(VGrid &v, float dt, LevelSet *ls)
{
    if(ls == NULL)
    {
        cerr << "Level set shouldn't be NULL" << endl;
        exit(1);
    }

    int res = m_nx - PADDING_WIDTH * 2;

    vector<Vector2i> bounds;
    int bWidth = res / THREAD_COUNT;

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        bounds.push_back(Vector2i(i*bWidth + PADDING_WIDTH, 
                    (i+1)*bWidth + PADDING_WIDTH));
    }
    
    bounds[THREAD_COUNT-1][1] = res + PADDING_WIDTH;

    DGrid update(res, m_ambient, m_dx);
    thread tt[THREAD_COUNT];

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i] = thread(&DGrid::advectHandler, this, 
                std::ref(bounds[i]), std::ref(v), dt, 
                ls, std::ref(update));
    }

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i].join();
    }
    
    *this = std::move(update);
    resetGradients();
}

void DGrid::advectHandler(Vector2i bounds, VGrid &v, 
            float dt, LevelSet *ls, DGrid &update)
{
    for(int i = bounds[0]; i < bounds[1]; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                if(ls->inFuelRegion(i, j, k))
                {
                    Vector3f x_g(i, j, k);
                    Vector3f vel = v.getVelocity(x_g);

                    Vector3f x_p = v.rk3(x_g, vel, dt);
                    update[i][j][k] = getQuantity(x_p);
                }
            }
        }
    }
}


bool DGrid::adjacentZero(int i, int j, int k, 
        IGrid &dist)
{
    Vector3i index(i, j, k);
    Vector2i shift(-1, 1);

    bool adjacent = false;
    int res = m_nx;

    for(int ax = 0; ax < DIM; ++ax)
    {
        Vector3i less = index;
        less[ax] -= 1;
        Vector3i more = index;
        more[ax] += 1;

        adjacent = adjacent ||
            (index[ax] > 0 && 
             dist.get(less) == 0) ||
            (index[ax] < res - 1 && 
             dist.get(more) == 0);
    }

    return adjacent;
}

// Extrapolation algorithm from p65 of Bridson
void DGrid::extrapolate(LevelSet &ls)
{
    int res = m_nx;
    IGrid dist(res, ls);

    queue<Vector3i> wavefront;

    // Initialize first wavefront
    for(int i = 0; i < m_nx; ++i) 
    {
        for(int j = 0; j < m_ny; ++j) 
        {
            for(int k = 0; k < m_nz; ++k)
            {
                if(dist[i][j][k] != 0 && 
                        adjacentZero(i, j, k, dist))
                {
                    dist[i][j][k] = 1;
                    wavefront.push(Vector3i(i, j, k));
                }
            }
        }
    }
    
    // Main loop
    while(wavefront.size() > 0)
    {
        Vector3i index = wavefront.front();
        wavefront.pop();

        float sum = 0.0f;
        int total = 0;

        for(int ax = 0; ax < DIM; ++ax)
        {
            Vector3i less = index;
            less[ax] -= 1;
            Vector3i more = index;
            more[ax] += 1;

            if(index[ax] > 0 && dist.get(less) < dist.get(index))
            {
                sum += get(less);
                ++total;
            }
            else if(index[ax] > 0 && dist.get(less) == INT_MAX)
            {
                dist.set(less, dist.get(index) + 1);
                wavefront.push(less);
            }

            if(index[ax] < res - 1 && 
                    dist.get(more) < dist.get(index))
            {
                sum += get(more);
                ++total;
            }
            else if(index[ax] < res - 1 && 
                    dist.get(more) == INT_MAX)
            {
                dist.set(more, dist.get(index) + 1);
                wavefront.push(more);
            }

        }
       
        int i = index[0];
        int j = index[1];
        int k = index[2];

        at(i)[j][k] = sum / (float)total;
    }

    resetGradients();
}

DGrid::IGrid::IGrid(int res, LevelSet &ls) :
    m_nx(res),
    m_ny(res),
    m_nz(res)
{
    resize(m_nx);

    for(int i = 0; i < m_nx; ++i) 
    {
        at(i).resize(m_ny);

        for(int j = 0; j < m_ny; ++j)
        {
            at(i)[j].resize(m_nz, INT_MAX);
        }
    }


    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH;
                ++j)
        {
            for(int k = PADDING_WIDTH;
                    k < m_nz - PADDING_WIDTH; ++k)
            {
                if(ls.inFuelRegion(i, j, k))
                {
                    at(i)[j][k] = 0;
                }
            }
        }
    }
}

inline int DGrid::IGrid::get(Vector3i index)
{
    return at(index[0])[index[1]][index[2]];
}

inline void DGrid::IGrid::set(Vector3i index, int val)
{
    at(index[0])[index[1]][index[2]] = val;
}
