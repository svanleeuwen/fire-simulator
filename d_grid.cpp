#include "d_grid.hpp"
#include "v_grid.hpp"

#include <iostream>
#include <queue>

#include <climits>

using std::cout;
using std::cerr;
using std::endl;
using std::queue;

DGrid::DGrid(int nx, int ny, float ambient)
    : QGrid(nx, ny, ambient)
{}

void DGrid::advect(VGrid &v, float dt, LevelSet *ls)
{
    if(ls == NULL)
    {
        cerr << "Level set shouldn't be NULL" << endl;
        exit(1);
    }

    int nx = m_nx - PADDING_WIDTH * 2;
    int ny = m_ny - PADDING_WIDTH * 2;

    DGrid update(nx, ny, m_ambient);

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i) 
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j) 
        {
            if(ls->inFuelRegion(i, j))
            {
                Vector2f x_g(i, j);
                Vector2f vel = v.getVelocity(x_g);

                Vector2f x_p = v.rk2(x_g, vel, dt);
                update[i][j] = getQuantity(x_p);
            }
        }
    }

    *this = std::move(update);
}

bool DGrid::adjacentZero(int i, int j, vector<vector<int>> &dist)
{
    return (i > 0 && dist[i-1][j] == 0)
        || (i < m_nx - 1 && dist[i+1][j] == 0)
        || (j > 0 && dist[i][j-1] == 0)
        || (j < m_ny - 1 && dist[i][j+1] == 0);
}

static void printDist(vector<vector<int>> &grid)
{
    for(int j = grid[0].size() - 1; j >= 0; --j) 
    {
        for(uint i = 0; i < grid.size(); ++i) 
        {
            cout << grid[i][j] << ", ";
        }

        cout << endl;
    }
}

// Extrapolation algorithm from p65 of Bridson
void DGrid::extrapolate(LevelSet &ls)
{
    vector<vector<int>> dist;
    dist.resize(m_nx);

    for(int i = 0; i < m_nx; ++i) 
    {
        dist[i].resize(m_ny, INT_MAX);
    }

    for(int i = PADDING_WIDTH; i < m_nx - PADDING_WIDTH; ++i)
    {
        for(int j = PADDING_WIDTH; j < m_ny - PADDING_WIDTH; ++j)
        {
            if(ls.inFuelRegion(i, j))
            {
                dist[i][j] = 0;
            }
        }
    }

    queue<Vector2i> wavefront;

    // Initialize first wavefront
    for(int i = 0; i < m_nx; ++i) 
    {
        for(int j = 0; j < m_ny; ++j) 
        {
            if(dist[i][j] != 0 && adjacentZero(i, j, dist))
            {
                dist[i][j] = 1;
                wavefront.push(Vector2i(i, j));
            }
        }
    }
    
    // Main loop
    while(wavefront.size() > 0)
    {
        int i = wavefront.front()[0];
        int j = wavefront.front()[1];
        wavefront.pop();

        float sum = 0.0f;
        int total = 0;

        if(i > 0 && dist[i-1][j] < dist[i][j])
        {
            sum += at(i-1)[j];
            ++total;
        }
        else if(i > 0 && dist[i-1][j] == INT_MAX)
        {
            dist[i-1][j] = dist[i][j] + 1;
            wavefront.push(Vector2i(i-1, j));
        }

        if(i < m_nx - 1 && dist[i+1][j] < dist[i][j])
        {
            sum += at(i+1)[j];
            ++total;
        }
        else if(i < m_nx - 1 && dist[i+1][j] == INT_MAX)
        {
            dist[i+1][j] = dist[i][j] + 1;
            wavefront.push(Vector2i(i+1, j));
        }

        if(j > 0 && dist[i][j-1] < dist[i][j])
        {
            sum += at(i)[j-1];
            ++total;
        }
        else if(j > 0 && dist[i][j-1] == INT_MAX)
        {
            dist[i][j-1] = dist[i][j] + 1;
            wavefront.push(Vector2i(i, j-1));
        }

        if(j < m_ny - 1 && dist[i][j+1] < dist[i][j])
        {
            sum += at(i)[j+1];
            ++total;
        }
        else if(j < m_ny - 1 && dist[i][j+1] == INT_MAX)
        {
            dist[i][j+1] = dist[i][j] + 1;
            wavefront.push(Vector2i(i, j+1));
        }

        at(i)[j] = sum / (float)total;
    }
}
