#define FRAMERATE 10
#define STEPS_PER_FRAME 3

#include <QtGui>
#include "paintcanvas.hpp"
#include "colourUtil.hpp"

#define FIRE_SIMULATION true 
#define MACSIZE 100

PaintCanvas::PaintCanvas(QWidget *parent) :
    QWidget(parent)
{
    m_mac = new MacGrid(MACSIZE, 1.0f / (FRAMERATE * STEPS_PER_FRAME), FIRE_SIMULATION);
    m_img = QImage(width(), height(), QImage::Format_RGB32);

    m_updateTimer = new QTimer(this);
    connect(m_updateTimer, SIGNAL(timeout()), this, 
            SLOT(refresh()));
    m_updateTimer->start(1000/FRAMERATE);
}

PaintCanvas::~PaintCanvas() 
{
}

QSize PaintCanvas::minimumSizeHint() const 
{
    return QSize(50, 50);
}   

QSize PaintCanvas::sizeHint() const 
{
    return QSize(600, 600);
}

void PaintCanvas::paintEvent(QPaintEvent* event) 
{
    (void) event;

    QPainter painter(this);
    painter.drawImage(QRect(0, 0, width(), height()), m_img);
}

void PaintCanvas::refresh()
{
    for(int i = 0; i < STEPS_PER_FRAME; ++i)
    {
        m_mac->step();
    }
    m_img = QImage(width(), height(), QImage::Format_RGB32);

    for(int i = 0; i < width(); ++i)
    {
        for(int j = 0; j < height(); ++j)
        {
            int grid_width = ceil(width() / (float)MACSIZE);
            int grid_height = ceil(height() / (float)MACSIZE);

            uint colour;

            if(!FIRE_SIMULATION)
            {
                colour = 255 * std::max(0.0f, m_mac->getDensity(
                    i / grid_width,
                    j / grid_height));
                colour = colour + (colour << 8) + (colour << 16);
            }
            else
            {
                int a = i / grid_width;
                int b = j / grid_height;

                if(!m_mac->isFuel(a, b) && m_mac->getTemp(a, b) > IGNITION_TEMP)
                {
                    colour = getBlackbodyRGB( m_mac->getTemp(
                        i / grid_width,
                        j / grid_height) * 2.0); 
                   /* *
                        std::max(0.0f, m_mac->getDensity(
                        i / grid_width,
                        j / grid_height));*/
                }
                else if(!m_mac->isFuel(a, b))
                {
                    colour = 255 * std::max(0.0f, m_mac->getDensity(
                        i / grid_width,
                        j / grid_height));
                    colour = colour + (colour << 8) + (colour << 16);
                }
                else
                {
                    colour = 255;
                }
            }

            m_img.setPixel(i, height() - j - 1, colour);
        }
    }

    update();
}
