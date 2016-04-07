#ifndef D_GRID_HPP
#define D_GRID_HPP

#include "q_grid.hpp"

class DGrid : public QGrid
{
public:
    DGrid(int res, float ambient, float dx);
    virtual void advect(VGrid &v, float dt, LevelSet *ls);
    
    void extrapolate(LevelSet &ls);

private:
    class IGrid : public vector< vector< vector<int> > >
    {
        public:
        IGrid(int res, LevelSet &ls);
        inline int get(Vector3i index);
        inline void set(Vector3i index, int val);
        
        private:
        int m_nx;
        int m_ny;
        int m_nz;
    };

    bool adjacentZero(int i, int j, int k,
        IGrid &dist);
    void initDistance(vector< vector< vector<int> > > &dist,
        LevelSet &ls);

    void advectHandler(Vector2i bounds, VGrid &v, 
            float dt, LevelSet *ls, DGrid &update);
};

#endif
