#ifndef LEVEL_SET_HPP 
#define LEVEL_SET_HPP

#include "grid.hpp"

class LevelSet : public Grid
{
public:
    LevelSet(int nx, int ny);

    void redistance();

private:
    void redistanceAdjacent();
    void sweep(bool increasing, Axis ax);
    float updateDistance(Vector2f phi);

    float m_dx;
};

#endif
