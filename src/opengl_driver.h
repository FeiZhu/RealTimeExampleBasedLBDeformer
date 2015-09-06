/*
 * @file: opengl_driver.h
 * @brief: UI interface for simulation
 * @author: Fei Zhu
 *
 */

#ifndef OPENGL_DRIVER_H_
#define OPENGL_DRIVER_H_

#include <string>

class GLUI;
class SphericalCamera;
class Lighting;

namespace RTLB{

class RealTimeExampleBasedDeformer;

class OpenGLDriver
{
public:
    explicit OpenGLDriver(const std::string &config_file_name);
    ~OpenGLDriver();

private:
    static OpenGLDriver* activeInstance();
    //init functions
    void initGLUT();
    void initGLUI();
    void initConfigurations(const std::string &config_file_name);
    void initSimulation();
    void initCamera();
    void initGraphics();
    //callback functions
    static void displayFunction(void);
    static void idleFunction(void);
    static void reshapeFunction(int width, int height);
    static void keyboardFunction(unsigned char key, int x, int y);
    static void specialFunction(int key, int x, int y);
    static void motionFunction(int x, int y);
    static void mouseFunction(int button, int state, int x, int y);
    //GLUI callback methods
    static void updateRenderMesh(int code);
    static void updateCurrentExample(int code);
    static void changeSimulationMode(int code);
    static void loadObjectEigenfunctions(int code);
    static void saveObjectEigenfunctions(int code);
    static void loadExampleEigenfunctions(int code);
    static void saveExampleEigenfunctions(int code);
    static void loadReducedBasis(int code);
    static void loadObjectCubicaData(int code);
    static void exitApplication(int code);
    //misc
    void drawAxis(double axis_length) const;
private:
    static OpenGLDriver *active_instance_;
    //simulation
    RealTimeExampleBasedDeformer *simulator_ = NULL;
    bool pause_simulation_ = true;
    std::string simulation_mesh_file_name_ = std::string("None");
    std::string example_file_name_prefix_ = std::string("None");
    unsigned int example_num_ = 0;
    std::string visual_mesh_file_name_ = std::string("None");
    std::string reduced_basis_file_name_ = std::string("None");
    std::string object_eigen_file_name_ = std::string("None");
    std::string example_eigen_file_name_ = std::string("None");
    std::string plane_file_name_ = std::string("None");
    std::string correspondence_file_name_ = std::string("None");
    double gravity_ = -9.8;
    double time_step_ = 1.0/30;
    //window
    std::string window_name_ = std::string("Example-based Simulator");
    int window_id_ = -1;
    unsigned int window_width_ = 800;
    unsigned int window_height_ = 600;
    //camera
    double znear_ = 0.01;
    double zfar_ = 10.0;
    double camera_radius_ = 17.5;
    double focus_position_[3];
    double camera_longitude_ = -60.0;
    double camera_lattitude_ = 20.0;
    SphericalCamera *camera_ = NULL;
    //light
    Lighting *lighting_ = NULL;
    std::string lighting_config_file_name_ = std::string("None");
    //mouse
    int mouse_pos_[2];
    bool left_button_down_ = false;
    bool middle_button_down_ = false;
    bool right_button_down_ = false;
    //render switches
    bool render_axis_ = true;
    bool render_vertices_ = false;
    bool render_wireframe_ = true;
    bool render_fixed_vertices_ = true;
    bool render_eigenfunction_ = false;
    enum RenderMeshType{
        VISUAL_MESH = 0,
        OBJECT_EIGEN_MESH,
        EXAMPLE_MESH
    };
    RenderMeshType render_mesh_type_ = VISUAL_MESH;
    //glui controls
    GLUI *glui_ = NULL;
};

}  //namespace RTLB

#endif //OPENGL_DRIVER_H_
