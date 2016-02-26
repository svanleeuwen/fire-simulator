#define FRAMERATE 8

#include <QtGui>
#include "paintcanvas.hpp"

// TODO: Draw when the timer runs out, don't start the computation
//          when it runs out

#define MACSIZE 100

PaintCanvas::PaintCanvas(QWidget *parent) :
    QWidget(parent)
{
    m_mac = new MacGrid(MACSIZE, 1.0f / FRAMERATE);
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
    m_mac->step();
    m_img = QImage(width(), height(), QImage::Format_RGB32);

    for(int i = 0; i < width(); ++i)
    {
        for(int j = 0; j < height(); ++j)
        {
            int grid_width = ceil(width() / (float)MACSIZE);
            int grid_height = ceil(height() / (float)MACSIZE);

            uint colour = 255 * std::max(0.0f, m_mac->getDensity(
                    i / grid_width,
                    j / grid_height));
            colour = colour + (colour << 8) + (colour << 16);

            m_img.setPixel(i, height() - j - 1, colour);
        }
    }

    update();
}
