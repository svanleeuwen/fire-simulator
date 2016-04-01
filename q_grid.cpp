#include "q_grid.hpp"
#include "v_grid.hpp"

#include <iostream>

using std::cout;
using std::endl;

QGrid::QGrid(int nx, int ny, float ambient) :
    Grid(nx, ny, ambient) 
{}

void QGrid::advect(VGrid &v, float dt, LevelSet *ls)
{
    (void) ls;

    int nx = m_nx - PADDING_WIDTH * 2;
    int ny = m_ny - PADDING_WIDTH * 2;

    QGrid update(nx, ny, m_ambient);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j) 
        {
            Vector2f x_g(i, j);
            Vector2f vel = v.getVelocity(x_g);

            Vector2f x_p = v.rk2(x_g, vel, dt);
            update[i][j] = getQuantity(x_p);
        }
    }

    *this = std::move(update);
}

float QGrid::getQuantity(const Vector2f &x)
{
    int i = (int)(x[0]) - 1;
    int j = (int)(x[1]) - 1;
    float s_x = x[0] - i - 1.0f;
    float s_y = x[1] - j - 1.0f;

    return getValue(Vector2i(i,j), Vector2f(s_x, s_y));
}
