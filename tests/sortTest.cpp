#include <Eigen/Core>
#include <iostream>

using Eigen::Vector3f;
using std::cout;
using std::endl;

static Vector3f sortDistances(Vector3f phi)
{
    bool flag = true;

    while(flag)
    {
        flag = false;

        for(int i = 0; i <=1; ++i)
        {
            if(phi[i+1] < phi[i])
            {
                flag = true;

                float temp = phi[i];
                phi[i] = phi[i+1];
                phi[i+1] = temp;
            }
        }
    }

    return phi;
}

int main()
{
    Vector3f a(0, 1, -2);
    Vector3f b(-3, 1, 2);
    Vector3f c(0, -3, 2);

    cout << "a: " << a << endl << "sorted:" << sortDistances(a)
        << endl;
    cout << "b: " << b << endl << "sorted:" << sortDistances(b)
        << endl;
    cout << "c: " << c << endl << "sorted:" << sortDistances(c)
        << endl;
}
