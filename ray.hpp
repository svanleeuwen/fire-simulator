#ifndef RAY_HPP
#define RAY_HPP

#include <Eigen/Core>

using Eigen::Vector3f;

class Ray
{
public:
    Ray(Vector3f origin, Vector3f dir);
    Ray(const Ray &other);
    Ray& operator=(const Ray &other);

    inline Vector3f operator() (float t) const
    {
        return m_origin + t*m_dir;    
    }

    Vector3f m_origin;
    Vector3f m_dir;
};
#endif
