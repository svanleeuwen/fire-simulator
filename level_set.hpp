#ifndef LEVEL_SET_HPP 
#define LEVEL_SET_HPP

#include "grid.hpp"

class VGrid;
class DGrid;

class LevelSet : public Grid
{
public:
    LevelSet(int nx, int ny);

    void redistance();
    void addCircle(float radius,
            Vector2f center,
            vector< vector<bool> > *src = NULL);
    void advect(float dt, VGrid &velocities, DGrid &burn);

    Vector2f getGradient(Vector2f x);
    float getGradientX(Vector2f x);
    float getGradientY(Vector2f x);

    bool inFuelRegion(int i, int j);

private:
    void redistanceAdjacent();
    void sweep(bool increasing, Axis ax);
    float updateDistance(Vector2f phi);

    Vector2f getUpwind(Vector2f w, Vector2i pos);

    float m_dx;
};

#endif
