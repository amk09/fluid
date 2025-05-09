#pragma once

#include "Eigen/Dense"

class Camera
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Camera();

    void setPosition(const Eigen::Vector3f      &position);
    void move       (const Eigen::Vector3f &deltaPosition);

    void setRotation(float      pitch, float      yaw);
    void rotate     (float deltaPitch, float deltaYaw);

    void lookAt(const Eigen::Vector3f &eye, const Eigen::Vector3f &target);

    void setOrbitPoint(const Eigen::Vector3f &target);
    bool getIsOrbiting();
    void setIsOrbiting(bool orbit);
    void toggleIsOrbiting();
    void zoom(float zoomMultiplier);

    const Eigen::Matrix4f &getView();
    const Eigen::Matrix4f &getProjection();
    const Eigen::Vector3f &getLook();

    void setPerspective(float fovY, float aspect, float near, float far);
    void setAspect(float aspect);

    const Eigen::Vector3f &getPosition() const { return m_position; }
    const Eigen::Vector3f &getOrbitPoint() const { return m_orbitPoint; }
    float getZoom() const { return m_zoom; }


private:
    void updateLook();
    void updatePitchAndYaw();

    // Do not mess with the order of these variables. Some Eigen voodoo will cause an inexplicable crash.

    Eigen::Vector3f m_position;

    float m_pitch, m_yaw;
    Eigen::Vector3f m_look;

    Eigen::Vector3f m_orbitPoint;
    bool m_isOrbiting;

    Eigen::Matrix4f m_view, m_proj;
    bool m_viewDirty, m_projDirty;

    float m_fovY, m_aspect, m_near, m_far;
    float m_zoom;
};
