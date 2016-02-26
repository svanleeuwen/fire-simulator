#ifndef PAINTWINDOW_HPP
#define PAINTWINDOW_HPP

#include <QMainWindow>
#include <QMenu>
#include <QAction>
#include <QtWidgets>
#include "paintcanvas.hpp"

class PaintWindow : public QMainWindow
{
    Q_OBJECT

public:
    PaintWindow();
    virtual ~PaintWindow() {}

private:
    virtual void keyPressEvent(QKeyEvent* event);

    PaintCanvas* m_canvas;
};

#endif
