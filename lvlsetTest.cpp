#define SIZE 12

#include "level_set.hpp"
#include "v_grid.hpp"

#include <iostream>

using std::cout;
using std::endl;

int main()
{
    LevelSet ls(SIZE, SIZE);

    for(int i = 0; i < SIZE + 2*PADDING_WIDTH; ++i)
    {
        for(int j = 0; j < SIZE + 2*PADDING_WIDTH; ++j)
        {
            ls[i][j] = 1;
        }
    }

    for(int i = 5; i < 9; ++i)
    {
        for(int j = 5; j < 9; ++j)
        {
            ls[i][j] = -1;
        }
    }

    cout << ls << endl;

    ls.redistance();

    cout << ls << endl;

    VGrid m_vel(SIZE, SIZE);

    for(uint i = 2; i < m_vel.u_.size() - 2; ++i) {
        for(uint j = 2; j < m_vel.u_[0].size() - 2; ++j) {
            m_vel.u_[i][j] = 1.0f;
        }
    }

    for(uint i = 2; i < m_vel.v_.size() - 2; ++i) {
        for(uint j = 2; j < m_vel.v_[0].size() - 2; ++j) {
            m_vel.v_[i][j] = 1.0f;
        }
    }

    for(int i = 0; i < 2; ++i)
    {
        ls.advect(0.25f, m_vel);
        ls.redistance();

        m_vel.advect(0.25, &ls);
    }
    
    cout << "V:" << endl << m_vel << endl;
    cout << "LS:" << endl << ls << endl;
}
