#ifndef GRID_HPP
#define GRID_HPP

#define EPSILON 1.0e-10f
#define PADDING_WIDTH 4
#define DIM 3
#define THREAD_COUNT 8

#include <vector>
#include <Eigen/Core>

using Eigen::Vector3f;
using Eigen::Vector3i;
using Eigen::Vector4f;
using Eigen::Vector2i;
using Eigen::Vector2f;

using std::vector;

class Grid : public vector< vector< vector<float> > >
{
public:
    enum Axis { X, Y, Z };
   
    Grid(int res, float ambient, float dx);
    Grid(int nx, int ny, int nz, float ambient, float dx);
    virtual ~Grid();
   
    float get(Vector3i index);

    float getValue(Vector3f x);
    float getValue(const Vector3i &index, const Vector3f &s);
    float getValue2D(int i, const Vector2i &index, 
            const Vector2f &s);
    float lerp(Vector3f x);

    Vector4f getW(const Vector4f &q);

    Vector3f getGradient(Vector3f x);
    float getGradient(Vector3f x, Axis ax);
    void resetGradients();

protected:
    void initGradient(Axis ax);
    
    int m_nx;
    int m_ny;
    int m_nz;

    float m_ambient;
    float m_dx;

    Grid *m_grad [3];

private:
    friend std::ostream& operator<<(std::ostream &out, Grid &g);
};
#endif
