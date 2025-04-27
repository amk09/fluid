#pragma once

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include "simulation.h"
#include "graphics/camera.h"
#include "graphics/shader.h"

#include <QOpenGLWidget>
#include <QElapsedTimer>
#include <QTimer>
#include <memory>

class GLWidget : public QOpenGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = nullptr);
    ~GLWidget();

private:
    static const int FRAMES_TO_AVERAGE = 30;

private:
    // Basic OpenGL Overrides
    void initializeGL()         override;
    void paintGL()              override;
    void resizeGL(int w, int h) override;

    // Event Listeners
    void mousePressEvent  (QMouseEvent *event) override;
    void mouseMoveEvent   (QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent       (QWheelEvent *event) override;
    void keyPressEvent    (QKeyEvent   *event) override;
    void keyReleaseEvent  (QKeyEvent   *event) override;

private:
    QElapsedTimer m_deltaTimeProvider; // For measuring elapsed time
    QTimer        m_intervalTimer;     // For triggering timed events

    Simulation m_sim;
    Camera     m_camera;
    Shader    *m_shader;

    int m_forward;
    int m_sideways;
    int m_vertical;
    bool m_pause = false;
    int m_lastX;
    int m_lastY;

    bool m_capture;

    // Variables to track interaction modes
    bool m_cameraControl = false;
    bool m_fluidInteract = false;

private slots:

    // Physics Tick
    void tick();
};
