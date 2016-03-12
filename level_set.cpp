#include "level_set.hpp"

#include <cfloat>
#include <iostream>
#include <stdlib.h>

using std::cerr;
using std::endl;

LevelSet::LevelSet(int nx, int ny) :
    Grid(nx, ny, FLT_MAX)
{
    m_dx = 1.0f / m_nx;
}

void LevelSet::redistance()
{
    redistanceAdjacent();

    for(int i = 0; i < DIM; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            sweep((bool)j, (Axis)i);
        }
    }
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
