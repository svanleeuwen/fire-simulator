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
/*
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
*/

Vector4f Grid::getW(const Vector4f &q)
{
    Vector4f w;

    w[0] = q[1];

    Vector2f d;
    d[0] = (q[2] - q[0]) / 2.0f;
    d[1] = (q[3] - q[1]) / 2.0f;

    float delta = q[2] - q[1];

    if(delta*d[0] <= 0)
    {
        d[0] = 0.0f;
    }

    if(delta*d[1] <= 0)
    {
        d[1] = 0.0f;
    }

    w[1] = d[0];
    w[2] = 2.0f*delta - 2*d[0] - d[1];
    w[3] = d[0] + d[1] - delta;

    return w;
}

// Monotonic cubic interpolation from Fedkiw 2001 for smoke
float Grid::getValue(const Vector2i &x, const Vector2f &s) {
    int i = x[0];
    int j = x[1];

    Vector4f q = {0.0f, 0.0f, 0.0f, 0.0f};

    for(int a = 0; a < 4; ++a) 
    {
        if(s[1] < EPSILON) 
        {
            q[a] = at(i + a)[j + 1];
        }
        else
        {
            Vector4f p;
            for(int c = 0; c < 4; ++c)
            {
                p[c] = at(i+a)[j+c];
            }

            Vector4f w_y = getW(p);
            
            for(int b = 0; b < 4; ++b) 
            {
                q[a] += w_y[b] * pow(s[1], b);
            }
        }
    }

    float val = 0.0f;
    Vector4f w_x = getW(q);

    for(int a = 0; a < 4; ++a) 
    {
        val += w_x[a] * pow(s[0], a);
    }

    return val;
}

float Grid::lerp(const Vector2f &x)
{
    int i = (int)x[0];
    int j = (int)x[1];

    float s_x = x[0] - i;
    float s_y = x[1] - j;

    return (1.0f - s_x) * (1.0f - s_y) * at(i)[j] 
        + s_x * (1.0f - s_y) * at(i+1)[j]
        + (1.0f - s_x) * s_y * at(i)[j+1]
        + s_x * s_y * at(i+1)[j+1];
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
