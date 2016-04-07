#ifndef S_GRID_HPP
#define S_GRID_HPP

#include "q_grid.hpp"

class SGrid : public QGrid
{
public:
    SGrid(int res, float ambient, float dx);
    virtual void advect(VGrid &v, float dt, LevelSet *ls = NULL);

private:
    void advectHandler(Vector2i bounds, VGrid &v, 
            float dt, LevelSet *ls, SGrid &update);
};

#endif
