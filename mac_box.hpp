#ifndef MAC_BOX_HPP
#define MAC_BOX_HPP

#include "mac_grid.hpp"
#include "ray.hpp"

#include <QImage>

using Eigen::Vector4i;

class MacBox
{
public:
    MacBox(int res, float dt, float scale);
    ~MacBox();

    void step();
    void computeImage(QImage &img);

private:
    void computeHandler(Vector4i bounds, QImage &img);
    uint traceSample(const Ray &r);
    
    Vector3f next(Ray r);
    bool inside(Vector3f point);

    MacGrid *m_mac;

    int m_res;
    float m_dx;

    int m_imgWidth;
    int m_imgHeight;

    Vector3f m_bbox_min;
    Vector3f m_bbox_max;
};

#endif
