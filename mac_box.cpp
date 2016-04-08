#include "mac_box.hpp"
#include "grid.hpp"
#include "colourUtil.hpp"

#include <thread>
#include <cfloat>

#include <iostream>

#define SAMPLE_WIDTH 1

using std::thread;

using std::cout;
using std::endl;

MacBox::MacBox(int res, float dt, float scale)
{
    m_mac = new MacGrid(res, dt, scale);
    m_res = res;
    m_dx = 2.0f / res;

    m_bbox_min = Vector3f(-1.0f, -1.0f, 1.0f);
    m_bbox_max = Vector3f(1.0f, 1.0f, 3.0f);
}

MacBox::~MacBox()
{
    delete m_mac;
}

void MacBox::step()
{
    m_mac->step();
}

void MacBox::computeImage(QImage &img)
{
    m_imgWidth = img.width();
    m_imgHeight = img.height();

    vector<Vector4i> bounds;
    int dx = m_imgWidth / THREAD_COUNT; 

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        bounds.push_back(Vector4i(i*dx, (i+1)*dx, 0, 
                    m_imgHeight));
    }

    bounds[THREAD_COUNT-1][1] = m_imgWidth;

    thread tt[THREAD_COUNT];

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i] = thread(&MacBox::computeHandler, this, 
                std::ref(bounds[i]),
                std::ref(img));
    }

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        tt[i].join();
    }
}

void MacBox::computeHandler(Vector4i bounds, QImage &img)
{
    Vector3f origin(0.0f, 0.0f, 0.0f);

    for(int i = bounds[0]; i < bounds[1]; ++i)
    {
        for(int j = bounds[2]; j < bounds[3]; ++j)
        {
            uint red = 0;
            uint green = 0;
            uint blue = 0;

            for(int a = 0; a < SAMPLE_WIDTH; ++a)
            {
                for(int b = 0; b < SAMPLE_WIDTH; ++b)
                {
                    float x = ((float)i + 
                            (float)a / SAMPLE_WIDTH) / 
                            m_imgWidth * 2.0f - 1.0f;

                    float y = ((float)j + 
                                (float)b / SAMPLE_WIDTH) / 
                            m_imgHeight * 2.0f - 1.0f;

                    Vector3f dir(x, y, 2.0f);
                    Ray ray(origin, dir);

                    uint colour = traceSample(ray);

                    red += (colour >> 16) % 256;
                    green += (colour >> 8) % 256;
                    blue += colour % 256;
                }
            }

            red /= (SAMPLE_WIDTH * SAMPLE_WIDTH);
            green /= (SAMPLE_WIDTH * SAMPLE_WIDTH);
            blue /= (SAMPLE_WIDTH * SAMPLE_WIDTH);

            uint colour = (red << 16) + (green << 8) + blue;

            img.setPixel(i, m_imgHeight - j - 1, colour);
        }
    }
}

uint MacBox::traceSample(const Ray &r)
{
    int z = Grid::Axis::Z;
    float first_t = m_bbox_min[z] / r.m_dir[z];

    Vector3f intersect = r(first_t);

    while(inside(intersect))
    {
        Vector3f trans = intersect + Vector3f(1.0f, 1.0f, -1.0f);

        int i = trans[0] / m_dx;
        int j = trans[1] / m_dx;
        int k = trans[2] / m_dx;
    
        if(trans[2] > 1.66f)
            return 0;

        if(!m_mac->isFuel(i, j, k) && 
                m_mac->getTemp(i, j, k) > IGNITION_TEMP &&
                m_mac->getDensity(i, j, k) > EPSILON)
        {
            return getBlackbodyRGB(m_mac->getTemp(i, j, k) * 
                    2.0); 
        }
        else if(m_mac->isFuel(i, j, k))
        {
            return 255;
        }

        intersect = next(Ray(intersect, r.m_dir));
    }

    return 0;
}

Vector3f MacBox::next(Ray r)
{
    Vector3f offset(1.0f, 1.0f, -1.0f);
    Vector3f trans = r.m_origin + offset;
        
    int i = trans[0] / m_dx;
    int j = trans[1] / m_dx;
    int k = trans[2] / m_dx;

    Vector3i index(i, j, k);
    float minDist = FLT_MAX;

    for(int ax = 0; ax < DIM; ++ax)
    {
        float sign = r.m_dir[ax] < 0 ? -1.0f : 1.0f;
        float next = ((index[ax] + sign) * m_dx - offset[ax]);

        float dist = (next - r.m_origin[ax] + m_dx / 2.0f) / r.m_dir[ax];

        minDist = fmin(dist, minDist);
    }

    return r(minDist);
}

bool MacBox::inside(Vector3f point)
{
    for(int ax = 0; ax < DIM; ++ax)
    {
        if(point[ax] < m_bbox_min[ax] || 
                point[ax] > m_bbox_max[ax])
        {
            return false;
        }
    }

    return true;
}

