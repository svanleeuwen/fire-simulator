#ifndef LEVEL_SET_HPP 
#define LEVEL_SET_HPP

#include "grid.hpp"

class VGrid;
class DGrid;

using Eigen::Vector2i;

class LevelSet : public Grid
{
public:
    LevelSet(int res);

    void redistance();
    void addCircle(float radius,
            Vector3f center,
            vector< vector< vector<bool> > > *src = NULL);
    void advect(float dt, VGrid &velocities, DGrid &burn);

    bool inFuelRegion(int i, int j, int k);
    void updateMeanCurvature(DGrid &curvature);

private:
    void redistanceAdjacent();
    void sweep(bool increasing, Axis ax);
    float updateDistance(Vector3f phi);

    void advectHandler(Vector2i bounds, VGrid &v, 
            float dt, DGrid &burn, LevelSet &update);
};

#endif
