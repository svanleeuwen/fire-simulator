#include "ray.hpp"

Ray::Ray(Vector3f origin, Vector3f dir) :
    m_origin(origin), m_dir(dir)
{}

Ray::Ray(const Ray &other)
{
    this->m_origin = other.m_origin;
    this->m_dir = other.m_dir;
}

Ray& Ray::operator=(const Ray &other)
{
    this->m_origin = other.m_origin;
    this->m_dir = other.m_dir;

    return *this;
}
