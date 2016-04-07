#ifndef T_GRID_HPP
#define T_GRID_HPP

#define AMBIENT_TEMP 273.0f // K

#define MAX_TEMP 2000.0f // K
#define IGNITION_TEMP 1000.0f // K

#define COOLING_CONSTANT MAX_TEMP * MAX_TEMP * 0.0001f

#include "q_grid.hpp"

class TGrid : public QGrid
{
public:
    TGrid(int res, float ambient, float dx);
    virtual void advect(VGrid &v, float dt, LevelSet *ls = NULL);

private:
    void advectHandler(Vector2i bounds, VGrid &v, 
            float dt, LevelSet *ls, TGrid &update);

    float decay(float temp, float dt);
};

#endif
