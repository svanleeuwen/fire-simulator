#ifndef GRID_HPP
#define GRID_HPP

#define EPSILON 1.0e-10f
#define PADDING_WIDTH 2
#define DIM 2

#include <vector>
#include <Eigen/Core>

using Eigen::Vector2f;
using Eigen::Vector2i;

using std::vector;

class Grid : public vector< vector<float> >
{
public:
    enum Axis { X, Y };
    
    Grid(int nx, int ny, float ambient);
    float getValue(const Vector2i &x, const Vector2f &s);

protected:
    int m_nx;
    int m_ny;
    float m_ambient;

private:
    friend std::ostream& operator<<(std::ostream &out, Grid &g);
};
#endif
