#ifndef S_GRID_HPP
#define S_GRID_HPP

#include "q_grid.hpp"

class SGrid : public QGrid
{
public:
    SGrid(int nx, int ny, float ambient);
    virtual void advect(VGrid &v, float dt, LevelSet *ls = NULL);
};

#endif
