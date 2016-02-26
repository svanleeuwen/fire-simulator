#define THIRD 1.0f/3.0f
#define SIXTH 1.0f/6.0f

#include "grid.hpp"

#include <iostream>

using std::cout;
using std::endl;

using Eigen::Vector4f;

Grid::Grid(int nx, int ny, float ambient) :
    m_ambient(ambient)
{
    m_nx = nx + PADDING_WIDTH * 2;
    m_ny = ny + PADDING_WIDTH * 2;

    resize(m_nx);

    for(int i = 0; i < m_nx; ++i) 
    {
        at(i).resize(m_ny, m_ambient);
    }
}

// w0 = w_-1 from Bridson p.40
inline static float w0(float s)
{
    return -THIRD * s * (1.0f + 0.5f * s * ( - 3.0 + s));
}

inline static float w1(float s) 
{
    return 1.0f - s * ( 0.5f + s * (1.0f - 0.5f * s));
}

inline static float w2(float s) 
{
    return s * (1.0f + 0.5f * (s - s*s));
}

inline static float w3(float s) 
{
    return SIXTH * s * (s*s - 1.0f);
}

float Grid::getValue(const Vector2i &x, const Vector2f &s) {
    int i = x[0];
    int j = x[1];

    Vector4f w_x = {w0(s[0]), w1(s[0]), w2(s[0]), w3(s[0])};
    Vector4f w_y = {w0(s[1]), w1(s[1]), w2(s[1]), w3(s[1])};

    vector<float> q;
    q.resize(4, 0.0f);

    for(int a = 0; a < 4; ++a) 
    {
        if(s[1] < EPSILON) 
        {
            q[a] = at(i + a)[j + 1];
        }
        else
        {
            for(int b = 0; b < 4; ++b) 
            {
                /*cout << "m_nx: " << m_nx << endl
                    << "m_ny: " << m_ny << endl
                    << "i_x: " << i + a << endl
                    << "i_y: " << j + b << endl; */
                q[a] += w_y[b] * at(i + a)[j + b];
            }
        }
    }

    float val = 0.0f;

    for(int a = 0; a < 4; ++a) 
    {
        val += w_x[a] * q[a];
    }

    return val;
}

std::ostream& operator<<(std::ostream &out, Grid &grid)
{
    for(int j = grid[0].size() - 1; j >= 0; --j) 
    {
        for(uint i = 0; i < grid.size(); ++i) 
        {
            out << grid[i][j] << ", ";
        }

        out << endl;
    }

    return out;
}
