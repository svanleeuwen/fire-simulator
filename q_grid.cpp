#include "q_grid.hpp"
#include "v_grid.hpp"

#include <iostream>
#include <thread>

using std::cout;
using std::endl;
using std::thread;

QGrid::QGrid(int res, float ambient, float dx) :
    Grid(res, ambient, dx) 
{}

void QGrid::advect(VGrid &v, float dt, LevelSet *ls)
{
    (void)ls;
    int res = m_nx - PADDING_WIDTH * 2;

    vector<Vector2i> bounds;
    int bWidth = res / THREAD_COUNT;

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        bounds.push_back(Vector2i(i*bWidth + PADDING_WIDTH, 
                    (i+1)*bWidth + PADDING_WIDTH));
    }
    
    bounds[THREAD_COUNT-1][1] = res + PADDING_WIDTH;

    QGrid update(res, m_ambient, m_dx);
    thread tt[THREAD_COUNT];

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i] = thread(&QGrid::advectHandler, this, 
                std::ref(bounds[i]), std::ref(v), dt, 
                std::ref(update));
    }

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i].join();
    }
    
    *this = std::move(update);
    resetGradients();
}

float QGrid::getQuantity(const Vector3f &x)
{
    int i = (int)(x[0]) - 1;
    int j = (int)(x[1]) - 1;
    int k = (int)(x[2]) - 1;
    float s_x = x[0] - i - 1.0f;
    float s_y = x[1] - j - 1.0f;
    float s_z = x[2] - k - 1.0f;

    return getValue(Vector3i(i,j,k), Vector3f(s_x, s_y, s_z));
}

void QGrid::advectHandler(Vector2i bounds, VGrid &v, 
            float dt, QGrid &update)
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
                update[i][j][k] = getQuantity(x_p);
            }
        }
    }
}
