#include "t_grid.hpp"
#include "v_grid.hpp"

#include <iostream>
#include <thread>

using std::cerr;
using std::endl;
using std::thread;

TGrid::TGrid(int res, float ambient, float dx)
    : QGrid(res, ambient, dx)
{}

void TGrid::advect(VGrid &v, float dt, LevelSet *ls)
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

    TGrid update(res, m_ambient, m_dx);
    thread tt[THREAD_COUNT];

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i] = thread(&TGrid::advectHandler, this, 
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

void TGrid::advectHandler(Vector2i bounds, VGrid &v, 
            float dt, LevelSet *ls, TGrid &update)
{
    for(int i = bounds[0]; i < bounds[1]; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            for(int k = PADDING_WIDTH; k < m_nz - PADDING_WIDTH;
                    ++k)
            {
                if(ls != NULL && ls->at(i)[j][k] < 0)
                {
                    update[i][j][k] = IGNITION_TEMP;
                    continue;
                }

                Vector3f x_g(i, j, k);
                Vector3f vel = v.getVelocity(x_g);

                Vector3f x_p = v.rk3(x_g, vel, dt);

                if(ls != NULL && ls->lerp(x_p) < 0 && 
                        ls->at(i)[j][k] > 0)
                {
                    update[i][j][k] = MAX_TEMP;
                }
                else
                {
                    float temp = getQuantity(x_p);
                    update[i][j][k] = decay(temp, dt);
                }
            }
        }
    }
}

// Uses analytical equation from p.137 of Bridson
float TGrid::decay(float temp, float dt)
{
    if(temp <= AMBIENT_TEMP + 1.0f)
    {
        return AMBIENT_TEMP;
    }

    float a = 1.0f / pow(temp - AMBIENT_TEMP, 3);
    float b = (3.0f * COOLING_CONSTANT * dt) /
        pow(MAX_TEMP - AMBIENT_TEMP, 4);

    return AMBIENT_TEMP + pow(a + b, -1.0f/3.0f);
}
