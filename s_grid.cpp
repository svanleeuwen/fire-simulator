#include "s_grid.hpp"
#include "v_grid.hpp"

#include <iostream>
#include <thread>

#define MAX_SOOT 0.5f

using std::cerr;
using std::endl;
using std::thread;

SGrid::SGrid(int res, float ambient, float dx) :
    QGrid(res, ambient, dx)
{}

void SGrid::advect(VGrid &v, float dt, LevelSet *ls)
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

    SGrid update(res, m_ambient, m_dx);
    thread tt[THREAD_COUNT];

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i] = thread(&SGrid::advectHandler, this, 
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
void SGrid::advectHandler(Vector2i bounds, VGrid &v, 
            float dt, LevelSet *ls, SGrid &update)
{
    for(int i = bounds[0]; i < bounds[1]; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                if(ls->at(i)[j][k] < 0)
                {
                    update[i][j][k] = 0.0f;
                }
                else
                {
                    Vector3f x_g(i, j, k);
                    Vector3f vel = v.getVelocity(x_g);

                    Vector3f x_p = v.rk3(x_g, vel, dt);

                    if(ls->lerp(x_p) < 0)
                    {
                        update[i][j][k] = MAX_SOOT;
                    }
                    else
                    {
                        update[i][j][k] = getQuantity(x_p);
                    }
                }
            }
        }
    }
}
