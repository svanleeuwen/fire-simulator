#include "t_grid.hpp"

TGrid::TGrid(int nx, int ny, float ambient)
    : QGrid(nx, ny, ambient)
{}

void TGrid::advect(VGrid &v, float dt, LevelSet *ls)
{
    int nx = m_nx - PADDING_WIDTH * 2;
    int ny = m_ny - PADDING_WIDTH * 2;

    TGrid update(nx, ny, m_ambient);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j) 
        {
            if(ls != NULL && ls->at(i)[j] < 0)
            {
                update[i][j] = IGNITION_TEMP;
                continue;
            }

            Vector2f x_g(i, j);
            Vector2f vel = v.getVelocity(x_g);

            Vector2f x_p = v.rk2(x_g, vel, dt);
            
            if(ls != NULL && ls->lerp(x_p) < 0 && ls->at(i)[j] > 0)
            {
                update[i][j] = MAX_TEMP;
            }
            else
            {
                float temp = getQuantity(x_p);
                update[i][j] = decay(temp, dt);
            }
        }
    }

    *this = std::move(update);
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
