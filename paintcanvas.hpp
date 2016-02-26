#ifndef PAINTCANVAS_HPP
#define PAINTCANVAS_HPP

#include <list>
#include <QWidget>
#include <QPainter>

#include "mac_grid.hpp"

class PaintCanvas : public QWidget {

    Q_OBJECT

public:
    PaintCanvas(QWidget *parent = 0);
    virtual ~PaintCanvas();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

protected:
    virtual void paintEvent(QPaintEvent* event);

private:
    MacGrid *m_mac;
    
    QTimer *m_updateTimer;
    QImage m_img;

private slots:
    void refresh();
};
#endif
