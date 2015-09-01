/*
 * @file: opengl_driver.cpp
 * @brief: UI interface for simulation
 * @author: Fei Zhu
 *
 */

#include "GL/glui.h"
#include "camera.h"
#include "lighting.h"
#include "real_time_example_based_deformer.h"
#include "opengl_driver.h"

namespace RTLB{

OpenGLDriver* OpenGLDriver::active_instance_ = NULL;

OpenGLDriver::OpenGLDriver(const std::string &config_file_name)
    :simulator_(NULL), pause_simulation_(true), example_num_(0),
     gravity_(-9.8), time_step_(1.0/30), window_id_(-1), glui_(NULL),
     window_width_(800),window_height_(600),
     znear_(0.01),zfar_(10.0),camera_radius_(17.5),
     camera_longitude_(-60.0),camera_lattitude_(20.0),
     camera_(NULL), lighting_(NULL), left_button_down_(false),
     middle_button_down_(false), right_button_down_(false),
     render_axis_(true), render_vertices_(false), render_wireframe_(true),
     render_fixed_vertices_(true), render_eigenfunction_(false),
{
    if(active_instance_)
        delete active_instance_;
    active_instance_ = this;
    //TO DO: init everything and enter mainloop
}

OpenGLDriver::~OpenGLDriver()
{
    exitApplication(0);
}

void OpenGLDriver::initConfigurations(const std::string &config_file_name)
{

}

OpenGLDriver* OpenGLDriver::activeInstance()
{
    return active_instance_;
}

void OpenGLDriver::initGLUT()
{

}

void OpenGLDriver::initGLUI()
{

}

void OpenGLDriver::initSimulation()
{

}

void OpenGLDriver::displayFunction()
{

}

void OpenGLDriver::idleFunction()
{

}

void OpenGLDriver::reshapeFunction(int width, int height)
{

}

void OpenGLDriver::keyboardFunction(unsigned char key, int x, int y)
{

}

void OpenGLDriver::specialFunction(int key, int x, int y)
{

}

void OpenGLDriver::motionFunction(int x, int y)
{

}

void OpenGLDriver::mouseFunction(int button, int state, int x, int y)
{

}

void OpenGLDriver::exitApplication(int code)
{
    //TO DO: release memories
}

}  //namespace RTLB
