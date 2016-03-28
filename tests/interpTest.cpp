#include <Eigen/Core>
#include <iostream>

using Eigen::Vector4f;
using Eigen::Vector2f;

using std::cout;
using std::endl;

Vector4f getW_x(const Vector4f &q)
{
    Vector4f w;

    w[0] = q[1];

    Vector2f d;
    d[0] = (q[2] - q[0]) / 2.0f;
    d[1] = (q[3] - q[1]) / 2.0f;

    float delta = q[2] - q[1];

    if(delta*d[0] <= 0)
    {
        d[0] = 0.0f;
    }

    if(delta*d[1] <= 0)
    {
        d[1] = 0.0f;
    }

    w[1] = d[0];
    w[2] = 2.0f*delta - 2*d[0] - d[1];
    w[3] = d[0] + d[1] - delta;

    return w;
}

int main()
{
    bool c = true;

    while(c)
    {
        Vector4f q(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX);
        cout << "q: " << q << endl;

        float s = rand() / (float)RAND_MAX;
        cout << "s: " << s << endl;

        Vector4f w = getW_x(q);
        cout << "w: " << w << endl;

        float val = w[0] + s * w[1] + s * s * w[2] + s * s * s * w[3];
        cout << "interp : " << val << endl;;

        if(val < 0 || val > 1)
        {
            c = false;
        }
    }
}
