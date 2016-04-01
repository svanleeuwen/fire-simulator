#include "d_grid.hpp"
#include "level_set.hpp"

#include <iostream>

#define RES 10

using std::cout;
using std::endl;

int main()
{
    LevelSet ls(RES, RES);
    DGrid g(RES, RES, 0.0f);

    for(int i = 3; i < 6; ++i)
    {
        for(int j = 3; j < 6; ++j)
        {
            ls[i][j] = -1.0f;
        }
    }

    for(int i = 2; i < RES + PADDING_WIDTH; ++i)
    {
        for(int j = 2; j < RES + PADDING_WIDTH; ++j)
        {
            if(ls.inFuelRegion(i, j))
                g[i][j] = 3.0f;
        }
    }

    ls.redistance();
    g.extrapolate(ls);

    cout << g << endl;     
}
