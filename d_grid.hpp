#ifndef D_GRID_HPP
#define D_GRID_HPP

#include "q_grid.hpp"

class DGrid : public QGrid
{
public:
    DGrid(int nx, int ny, float ambient);
    virtual void advect(VGrid &v, float dt, LevelSet *ls);
    
    void extrapolate(LevelSet &ls);

private:
    bool adjacentZero(int i, int j, vector<vector<int>> &dist);
};

#endif
