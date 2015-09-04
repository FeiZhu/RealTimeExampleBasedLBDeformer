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
    void exitApplication(int code);
    //misc
    void drawAxis(double axis_length) const;
private:
    static OpenGLDriver *active_instance_;
    //simulation
    RealTimeExampleBasedDeformer *simulator_;
    bool pause_simulation_;
    std::string simulation_mesh_file_name_;
    std::string example_file_name_prefix_;
    unsigned int example_num_;
    std::string visual_mesh_file_name_;
    std::string reduced_basis_file_name_;
    std::string object_eigen_file_name_;
    std::string example_eigen_file_name_;
    std::string plane_file_name_;
    std::string correspondence_file_name_;
    double gravity_;
    double time_step_;
    //window
    std::string window_name_;
    int window_id_;
    GLUI *glui_;
    unsigned int window_width_;
    unsigned int window_height_;
    //camera
    double znear_;
    double zfar_;
    double camera_radius_;
    double focus_position_[3];
    double camera_longitude_;
    double camera_lattitude_;
    SphericalCamera *camera_;
    //light
    Lighting *lighting_;
    //mouse
    int mouse_pos_[2];
    bool left_button_down_;
    bool middle_button_down_;
    bool right_button_down_;
    //render
    bool render_axis_;
    bool render_vertices_;
    bool render_wireframe_;
    bool render_fixed_vertices_;
    bool render_eigenfunction_;
    //glui controls
};

}  //namespace RTLB

#endif //OPENGL_DRIVER_H_
