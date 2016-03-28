#include "s_grid.hpp"

#define MAX_SOOT 0.5f

SGrid::SGrid(int nx, int ny, float ambient) :
    QGrid(nx, ny, ambient)
{}

void SGrid::advect(VGrid &v, float dt, LevelSet *ls)
{
    if(ls == NULL)
    {
        QGrid::advect(v, dt, ls);
    }

    int nx = m_nx - PADDING_WIDTH * 2;
    int ny = m_ny - PADDING_WIDTH * 2;

    SGrid update(nx, ny, m_ambient);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j) 
        {
            if(ls->at(i)[j] < 0)
            {
                update[i][j] = 0.0f;
            }
            else
            {
                Vector2f x_g(i, j);
                Vector2f vel = v.getVelocity(x_g);

                Vector2f x_p = v.rk2(x_g, vel, dt);

                if(ls->lerp(x_p) < 0)
                {
                    update[i][j] = MAX_SOOT;
                }
                else
                {
                    update[i][j] = getQuantity(x_p);
                }
            }
        }
    }

    *this = std::move(update);
}
