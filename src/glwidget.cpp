#include "glwidget.h"

#include <QApplication>
#include <QKeyEvent>
#include <iostream>

#define SPEED 1.5
#define ROTATE_SPEED 0.0025

using namespace std;

GLWidget::GLWidget(QWidget *parent) :
    QOpenGLWidget(parent),
    m_deltaTimeProvider(),
    m_intervalTimer(),
    m_sim(),
    m_camera(),
    m_shader(),
    m_forward(),
    m_sideways(),
    m_vertical(),
    m_lastX(),
    m_lastY(),
    m_capture(false)
{
    // GLWidget needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    QApplication::setOverrideCursor(Qt::ArrowCursor);

    // GLWidget needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // Function tick() will be called once per interva
    connect(&m_intervalTimer, SIGNAL(timeout()), this, SLOT(tick()));
}

GLWidget::~GLWidget()
{
    if (m_shader != nullptr) delete m_shader;
}

// ================== Basic OpenGL Overrides


void GLWidget::initializeGL()
{
    // Initialize GL extension wrangler
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK) fprintf(stderr, "Error while initializing GLEW: %s\n", glewGetErrorString(err));
    fprintf(stdout, "Successfully initialized GLEW %s\n", glewGetString(GLEW_VERSION));

    // Set clear color to black
    glClearColor(0, 0, 0, 1);

    // Enable depth-testing and backface culling
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Initialize the shader and simulation
    m_shader = new Shader(":/resources/shaders/shader.vert", ":/resources/shaders/shader.frag");
    m_sim.init();

    // Initialize camera with a closer position
    Eigen::Vector3f eye    = {0, 2, -8};
    Eigen::Vector3f target = {0, 1,  0};
    m_camera.lookAt(eye, target);
    m_camera.setOrbitPoint(target);

    m_camera.setPerspective(120, width() / static_cast<float>(height()), 0.1, 50);

    m_deltaTimeProvider.start();
    m_intervalTimer.start(1000 / 60);
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_BLEND); // Enable transparency blending
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    m_shader->bind();
    // m_shader->setUniform("proj", m_camera.getProjection());
    // m_shader->setUniform("view", m_camera.getView());
    Eigen::Matrix4f view = m_camera.getView();
    Eigen::Matrix4f proj = m_camera.getProjection();
    Eigen::Matrix4f invViewProj = (proj * view).inverse();
    m_shader->setUniform("invViewProj", invViewProj);
    m_sim.draw(m_shader);
    m_shader->unbind();

    glDisable(GL_BLEND);
}

void GLWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    m_camera.setAspect(static_cast<float>(w) / h);
}

// ================== Event Listeners

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    m_lastX = event->position().x();
    m_lastY = event->position().y();

    // Right button controls camera rotation
    if (event->button() == Qt::RightButton) {
        m_cameraControl = true;
    }
    // Left button controls fluid interaction
    else if (event->button() == Qt::LeftButton) {
        m_fluidInteract = true;
        // Immediately trigger fluid interaction
        m_sim.handleMousePress(m_lastX, m_lastY, width(), height());
    }
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    int currX = event->position().x();
    int currY = event->position().y();

    // Right button drag - camera control
    if (m_cameraControl) {
        int deltaX = currX - m_lastX;
        int deltaY = currY - m_lastY;

        if (deltaX != 0 || deltaY != 0) {
            m_camera.rotate(deltaY * ROTATE_SPEED, -deltaX * ROTATE_SPEED);
        }
    }
    // Left button drag - fluid interaction
    else if (m_fluidInteract) {
        m_sim.handleMouseMove(currX, currY, width(), height());
    }

    m_lastX = currX;
    m_lastY = currY;
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::RightButton) {
        m_cameraControl = false;
    }
    else if (event->button() == Qt::LeftButton) {
        m_fluidInteract = false;
    }
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
    float zoom = 1 - event->pixelDelta().y() * 0.1f / 120.f;
    m_camera.zoom(zoom);
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  += SPEED; break;
    case Qt::Key_S: m_forward  -= SPEED; break;
    case Qt::Key_A: m_sideways -= SPEED; break;
    case Qt::Key_D: m_sideways += SPEED; break;
    case Qt::Key_F: m_vertical -= SPEED; break;
    case Qt::Key_R: m_vertical += SPEED; break;
    case Qt::Key_C: m_camera.toggleIsOrbiting(); break;
    case Qt::Key_T: m_sim.toggleWire(); break;

    // Added color and rendering mode switch keys
    case Qt::Key_1:
        m_sim.fluidCube.setColorMap(1);
        break;
    case Qt::Key_2:
        m_sim.fluidCube.setColorMap(2);
        break;
    case Qt::Key_3:
        m_sim.fluidCube.setColorMap(3);
        break;
    case Qt::Key_4:
        m_sim.fluidCube.setColorMap(4);
        break;
    case Qt::Key_0:
        m_sim.fluidCube.setColorMap(0);
        break;

    case Qt::Key_M:
        // Toggle shell rendering
        m_sim.fluidCube.setRenderMode(
            m_sim.fluidCube.getRenderMode() == 0 ? 1 : 0);
        std::cout << "Render Mode: " << (m_sim.fluidCube.getRenderMode() == 0 ? "Volume" : "Shell") << std::endl;
        break;

    // Add vorticity strength control
    case Qt::Key_V: // Increase vorticity
        m_sim.fluidCube.setVorticityStrength(
            m_sim.fluidCube.getVorticityStrength() + 0.2f);
        std::cout << "Vorticity Strength: " << m_sim.fluidCube.getVorticityStrength() << std::endl;
        break;

    case Qt::Key_B: // Decrease vorticity
        m_sim.fluidCube.setVorticityStrength(
            std::max(0.0f, m_sim.fluidCube.getVorticityStrength() - 0.2f));
        std::cout << "Vorticity Strength: " << m_sim.fluidCube.getVorticityStrength() << std::endl;
        break;

    case Qt::Key_Escape: QApplication::quit();
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  -= SPEED; break;
    case Qt::Key_S: m_forward  += SPEED; break;
    case Qt::Key_A: m_sideways += SPEED; break;
    case Qt::Key_D: m_sideways -= SPEED; break;
    case Qt::Key_F: m_vertical += SPEED; break;
    case Qt::Key_R: m_vertical -= SPEED; break;
    case Qt::Key_P: m_pause = !m_pause; break;
    }
}

// ================== Physics Tick

void GLWidget::tick()
{
    // Run at 60 FPS
    float accumulatedTime = 0.0f;
    const float fixedDelta = 1.0f / 120.0f;

    float deltaSeconds = m_deltaTimeProvider.restart() / 1000.f;
    accumulatedTime += deltaSeconds;

    while (accumulatedTime >= fixedDelta && !m_pause) {
        m_sim.update(fixedDelta);
        accumulatedTime -= fixedDelta;
    }

    // Move camera
    auto look = m_camera.getLook();
    look.y() = 0;
    look.normalize();
    Eigen::Vector3f perp(-look.z(), 0, look.x());
    Eigen::Vector3f moveVec = m_forward * look.normalized() + m_sideways * perp.normalized() + m_vertical * Eigen::Vector3f::UnitY();
    moveVec *= deltaSeconds;
    m_camera.move(moveVec);

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
