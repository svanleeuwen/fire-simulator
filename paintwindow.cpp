#include "paintwindow.hpp"

PaintWindow::PaintWindow() 
{
    setWindowTitle("Fire Simulator");

    QVBoxLayout *layout = new QVBoxLayout;

    m_canvas = new PaintCanvas(this);
    layout->addWidget(m_canvas);

    setCentralWidget(new QWidget);
    centralWidget()->setLayout(layout);
}

void PaintWindow::keyPressEvent(QKeyEvent *event) {
    switch(event->key()) {
        case Qt::Key_Q: 
            QCoreApplication::instance()->quit();
            return;
        default: 
            QWidget::keyPressEvent(event);
            return;
    }
}
