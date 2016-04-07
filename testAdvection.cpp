#define RES 10

#include "v_grid.hpp"
#include "level_set.hpp"
#include "d_grid.hpp"
#include "mac_grid.hpp"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>

using std::cout;
using std::endl;

using Eigen::Triplet;

int main()
{
    LevelSet ls(RES);
    DGrid burn(RES, 0.08, 1.0f/RES);
    VGrid vel(RES, RES, RES, 1.0f/RES, &ls, &burn);

    MacGrid mac(RES, 0.25f, 1.0f);

    for(int i = PADDING_WIDTH; i < RES + PADDING_WIDTH + 1; ++i)
    {
        for(int j = PADDING_WIDTH; j < RES + PADDING_WIDTH + 1; ++j)
        {
            for(int k = PADDING_WIDTH; k < RES + PADDING_WIDTH + 1; ++k)
            {
                if(j < RES + PADDING_WIDTH && k < RES + PADDING_WIDTH)
                    vel.u_[i][j][k] = 1.0f;
                if(i < RES + PADDING_WIDTH && k < RES + PADDING_WIDTH)
                    vel.v_[i][j][k] = 1.0f;
                if(j < RES + PADDING_WIDTH && i < RES + PADDING_WIDTH)
                    vel.w_[i][j][k] = 1.0f;
            }
        }
    }

    
    for(int i = 0; i < 5; ++i)
    {
        mac.step();
    }
}
