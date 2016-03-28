#ifndef T_GRID_HPP
#define T_GRID_HPP

#define AMBIENT_TEMP 273.0f // degK

#define MAX_TEMP 2000.0f
#define IGNITION_TEMP 1000.0f

// Made up...might have to change
#define COOLING_CONSTANT MAX_TEMP * MAX_TEMP * 0.0001f

#include "q_grid.hpp"

class TGrid : public QGrid
{
public:
    TGrid(int nx, int ny, float ambient);
    virtual void advect(VGrid &v, float dt, LevelSet *ls = NULL);

protected:
    using QGrid::getQuantity;
    float getQuantity(const Vector2f &x, LevelSet *ls);
    float getValue(const Vector2i &x, const Vector2f &s,
        LevelSet *ls);

private:
    float decay(float temp, float dt);
};

#endif
