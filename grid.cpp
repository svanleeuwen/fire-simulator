#define THIRD 1.0f/3.0f
#define SIXTH 1.0f/6.0f

#include "grid.hpp"

#include <iostream>

using std::cout;
using std::endl;

using Eigen::Vector2f;
using Eigen::Vector4f;

Grid::Grid(int res, float ambient, float dx) :
    Grid(res, res, res, ambient, dx)
{}

Grid::Grid(int nx, int ny, int nz, float ambient, float dx) :
    m_ambient(ambient),
    m_dx(dx)
{
    m_nx = nx + PADDING_WIDTH * 2;
    m_ny = ny + PADDING_WIDTH * 2;
    m_nz = nz + PADDING_WIDTH * 2;

    resize(m_nx);

    for(int i = 0; i < m_nx; ++i) 
    {
        at(i).resize(m_ny);

        for(int j = 0; j < m_ny; ++j)
        {
            at(i)[j].resize(m_nz, m_ambient);
        }
    }

    m_grad[Axis::X] = NULL;
    m_grad[Axis::Y] = NULL;
    m_grad[Axis::Z] = NULL;
}

Grid::~Grid()
{
    resetGradients();
}

float Grid::get(Vector3i index)
{
    return at(index[0])[index[1]][index[2]];
}

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

float Grid::getValue(Vector3f x)
{
    int i = (int)(x[0]) - 1;
    int j = (int)(x[1]) - 1;
    int k = (int)(x[2]) - 1;
    float s_x = x[0] - i - 1.0f;
    float s_y = x[1] - j - 1.0f;
    float s_z = x[2] - k - 1.0f;

//    Vector3f  temp= x + Vector3f(1, 1,1);
//    return lerp(temp);
    return getValue(Vector3i(i,j,k), Vector3f(s_x, s_y, s_z));
}

// Monotonic cubic interpolation from Fedkiw 2001 for smoke
float Grid::getValue(const Vector3i &index, const Vector3f &s) 
{
//    Vector3f temp(index[0] +1+ s[0],index[1] +1+ s[1],index[2] +1+ s[2]);
//    return lerp(temp);

    int i = index[0];
    int j = index[1];
    int k = index[2];

    Vector4f r(0.0f, 0.0f, 0.0f, 0.0f);
    
    for(int a = 0; a < 4; ++a) 
    {
        Vector4f q = Vector4f(0.0f, 0.0f, 0.0f, 0.0f);

        for(int b = 0; b < 4; ++b)
        {
            Vector4f p;
            for(int c = 0; c < 4; ++c)
            {
                p[c] = at(i+a)[j+b][k+c];
            }

            Vector4f w_z = getW(p);

            for(int d = 0; d < 4; ++d)
            {
                q[b] += w_z[d] * pow(s[2], d);
            }
        }

        Vector4f w_y = getW(q);

        for(int b = 0; b < 4; ++b)
        {
            r[a] += w_y[b] * pow(s[1], b);
        }
    }

    float val = 0.0f;
    Vector4f w_x = getW(r);

    for(int a = 0; a < 4; ++a) 
    {
        val += w_x[a] * pow(s[0], a);
    }

    return val;
}
/*
float Grid::getValue(const Vector3i &index, const Vector3f &s)
{
    int i = index[0];
    int j = index[1];
    int k = index[2];

    Vector4f q(0.0f, 0.0f, 0.0f, 0.0f);

    for(int a = 0; a < 4; ++a)
    {
        q[a] = getValue2D(i+a, Vector2i(j, k), 
                Vector2f(s[1], s[2]));
    }

    float val = 0.0f;
    Vector4f w_x = getW(q);

    for(int a = 0; a < 4; ++a) 
    {
        val += w_x[a] * pow(s[0], a);
    }

    return val;
}*/

float Grid::getValue2D(int i, const Vector2i &index, 
        const Vector2f &s) 
{
    int j = index[0];
    int k = index[1];

    Vector4f q = {0.0f, 0.0f, 0.0f, 0.0f};

    for(int a = 0; a < 4; ++a) 
    {
        Vector4f p;
        for(int c = 0; c < 4; ++c)
        {
            p[c] = at(i)[j+a][k+c];
        }

        Vector4f w_y = getW(p);
        
        for(int b = 0; b < 4; ++b) 
        {
            q[a] += w_y[b] * pow(s[1], b);
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

float Grid::lerp(Vector3f x)
{
    int i = (int)x[0];
    int j = (int)x[1];
    int k = (int)x[2];

    float s_x = x[0] - i;
    float s_y = x[1] - j;
    float s_z = x[2] - k;

    Vector2f p_x = { 1.0f - s_x, s_x };
    Vector2f p_y = { 1.0f - s_y, s_y };
    Vector2f p_z = { 1.0f - s_z, s_z };

    float val = 0.0f;

    for(int a = 0; a < 2; ++a)
    {
        for(int b = 0; b < 2; ++b)
        {
            for(int c = 0; c < 2; ++c)
            {
                val += p_x[a] * p_y[b] * p_z[c] 
                    * at(i+a)[j+b][k+c];
            }
        }
    }

    return val;
}

Vector3f Grid::getGradient(Vector3f x)
{
    return Vector3f(getGradient(x, Axis::X),
            getGradient(x, Axis::Y),
            getGradient(x, Axis::Z));
}

float Grid::getGradient(Vector3f x, Axis ax)
{
    if(m_grad[ax] == NULL)
    {
        initGradient(ax);
    }

    return m_grad[ax]->lerp(x);
}

void Grid::resetGradients()
{
    for(int i = 0; i < DIM; ++i)
    {
        if(m_grad[i] != NULL)
        {
            delete m_grad[i];
            m_grad[i] = NULL;
        }
    }
}

void Grid::initGradient(Axis ax)
{
    m_grad[ax] = new Grid(m_nx - PADDING_WIDTH*2, 
            m_ny - PADDING_WIDTH*2, 
            m_nz - PADDING_WIDTH*2, 0.0f, m_dx);
    Vector3i dim(m_nx, m_ny, m_nz);

    for(int i = 0; i < m_nx; ++i)
    {
        for(int j = 0; j < m_ny; ++j)
        {
            for(int k = 0; k < m_nz; ++k)
            {
                Vector3i index(i, j, k);

                Vector3i less = index;
                less[ax] -= 1;
                
                Vector3i more = index;
                more[ax] += 1;

                if(index[ax] == dim[ax] - 1)
                {
                    (m_grad[ax])->at(i)[j][k] = 
                        (get(index) - get(less))
                        / m_dx;
                }
                else if(index[ax] == 0)
                {
                    (m_grad[ax])->at(i)[j][k] = 
                        (get(more) - get(index))
                        / m_dx;
                }
                else
                {
                    (m_grad[ax])->at(i)[j][k] = 
                        (get(more) - get(less))
                        / (2.0f * m_dx);
                }
            }
        }
    }
}

std::ostream& operator<<(std::ostream &out, Grid &grid)
{
    for(uint k = 0; k < grid[0][0].size(); ++k)
    {
        out << "z = " << k << endl;
        for(int j = grid[0].size() - 1; j >= 0; --j) 
        {
            for(uint i = 0; i < grid.size(); ++i) 
            {
                out << grid[i][j][k] << ", ";
            }

            out << endl;
        }

        out << endl;
    }

    return out;
}
