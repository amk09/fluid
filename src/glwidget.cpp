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

    Eigen::Vector3f eye    = {0, 1.5, -10};
    Eigen::Vector3f target = {0.5, 0.5,  0.5};
    m_camera.lookAt(eye, target);

    // Auto Camera OrbitPoint
    m_camera.setOrbitPoint(target);
    m_camera.setIsOrbiting(true);

    m_camera.setPerspective(120, width() / static_cast<float>(height()), 0.1, 50);

    m_deltaTimeProvider.start();
    m_intervalTimer.start(1000 / 60);
}

//void GLWidget::dumpFrame()
//{
//    static int frame = 0;                     // or a member variable
//    QImage img = grabFramebuffer();           // <- RGBA8888 by default
//    // Make sure the directory exists; Qt will create parent dirs for you
//    const QString path = QStringLiteral("captures/frame_%1.png")
//                             .arg(frame++, 6, 10, QLatin1Char('0'));
//    img.save(path);                           // writes PNG
//}

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
//    if (m_dumpFrames) {
//        dumpFrame();
//    }

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
        std::cout << "Cyan-blue" << std::endl;
        break;
    case Qt::Key_2:
        m_sim.fluidCube.setColorMap(2);
        std::cout << "Purple Electric " << std::endl;
        break;
    case Qt::Key_3:
        m_sim.fluidCube.setColorMap(3);
        std::cout << "Cyan-yellow" << std::endl;
        break;
    case Qt::Key_4:
        m_sim.fluidCube.setColorMap(4);
        std::cout << "Orange-grey" << std::endl;
        break;
    case Qt::Key_5:
        m_sim.fluidCube.setColorMap(5); // New fire gradient
        std::cout << "Fire gradient" << std::endl;
        break;
    case Qt::Key_6:
        m_sim.fluidCube.setColorMap(6); // New ocean gradient
        std::cout << "Ocean gradient" << std::endl;
        break;
    case Qt::Key_7:
        m_sim.fluidCube.setColorMap(7); // New plasma gradient
        std::cout << "Plasma gradient" << std::endl;
        break;
    case Qt::Key_8:
        m_sim.fluidCube.setColorMap(8); // New Rainbow gradient
        std::cout << "Rainbow" << std::endl;
        break;
    case Qt::Key_9:
        m_sim.fluidCube.setColorMap(9); // New Aurora gradient
        std::cout << "Aurora" << std::endl;
        break;
    case Qt::Key_L:
        m_sim.fluidCube.setColorMap(10); // New Lava gradient
        std::cout << "Lava" << std::endl;
        break;
    case Qt::Key_0:
        m_sim.fluidCube.setColorMap(0);
        std::cout << "Default Water - Color type set to: " <<
            m_sim.fluidCube.getCurrentColorType() << std::endl;
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
    case Qt::Key_O:
        m_sim.addObstacle();
        std::cout << "Add Obstacle" << endl;
            break;
    case Qt::Key_I:
        m_sim.clearObstacle();
        std::cout << "Clear Obstacle" << endl;
        break;
    case Qt::Key_Up:
        m_sim.moveUp();
        std::cout << "Move Up" << std::endl;
        break;

    case Qt::Key_Down:
        m_sim.moveDown();
        std::cout << "Move Down" << std::endl;
        break;

    case Qt::Key_Left:
        m_sim.moveLeft();
        std::cout << "Move Left" << std::endl;
        break;

    case Qt::Key_Right:
        m_sim.moveRight();
        std::cout << "Move Right" << std::endl;
        break;

    case Qt::Key_Equal:
        m_sim.moveForward();
        std::cout << "Move Forward (Z++)" << std::endl;
        break;

    case Qt::Key_Minus:
        m_sim.moveBackward();
        std::cout << "Move Backward (Z--)" << std::endl;
        break;

    // Add clear fluid shortcuts
    case Qt::Key_Space:
    case Qt::Key_Delete:
    case Qt::Key_Return:
        m_sim.clearAllFluids();
        cout << "Cleared all fluids" << std::endl;
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
    static float accumulatedTime = 0.0f;
    const float fixedDelta = 1.0f / m_sim.fluidCube.frameRate;

    float deltaSeconds = m_deltaTimeProvider.restart() / 1000.f;
    accumulatedTime += deltaSeconds;

    while (accumulatedTime >= fixedDelta && !m_pause) {
        m_sim.update(fixedDelta);
        accumulatedTime -= fixedDelta;
    }

    //Auto Camera
    if (m_camera.getIsOrbiting()) {
        Eigen::Vector3f target = {0, 0,  0};
        Eigen::AngleAxisf rotation_y(M_PI/ 360, Eigen::Vector3f::UnitY());
        m_camera.setPosition( rotation_y.toRotationMatrix() * m_camera.getPosition());
        m_camera.lookAt(m_camera.getPosition(), target);
    }

    //Move camera
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
