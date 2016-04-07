#ifndef Q_GRID_HPP
#define Q_GRID_HPP

#include "grid.hpp"
#include "level_set.hpp"

#include <vector>

using std::vector;
using Eigen::Vector2i;

class VGrid;

class QGrid : public Grid 
{
public:
    QGrid(int res, float ambient, float dx);

    // Must advect quantities before velocities
    virtual void advect(VGrid &v, float dt, LevelSet *ls = NULL);
    virtual float getQuantity(const Vector3f &x);

private:
    void advectHandler(Vector2i bounds, VGrid &v, 
            float dt, QGrid &update);
};
#endif
