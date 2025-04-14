#include "mainwindow.h"
#include <QHBoxLayout>

MainWindow::MainWindow()
{
    glWidget = new GLWidget();

    QHBoxLayout *container = new QHBoxLayout;
    container->addWidget(glWidget);
    this->setLayout(container);

    // Set a larger window size for better visualization
    this->resize(1024, 768);
}

MainWindow::~MainWindow()
{
    delete glWidget;
}
