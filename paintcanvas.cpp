#define FRAMERATE 20
#define STEPS_PER_FRAME 1

#include <QtGui>
#include <iostream>

#include "paintcanvas.hpp"
#include "colourUtil.hpp"

#define MACSIZE 100
#define SCALE 1.0f

using std::cout;
using std::endl;

PaintCanvas::PaintCanvas(QWidget *parent) :
    QWidget(parent)
{
    m_mac = new MacBox(MACSIZE, 1.0f / (FRAMERATE * STEPS_PER_FRAME), SCALE);
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
    m_mac->computeImage(m_img);

    update();
}
