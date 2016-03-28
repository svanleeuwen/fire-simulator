#ifndef Q_GRID_HPP
#define Q_GRID_HPP

#include "grid.hpp"
#include "v_grid.hpp"

#include <vector>

using std::vector;

class QGrid : public Grid 
{
public:
    QGrid(int nx, int ny, float ambient);

    // Must advect quantities before velocities
    virtual void advect(VGrid &v, float dt, LevelSet *ls = NULL);
    virtual float getQuantity(const Vector2f &x);
};
#endif
