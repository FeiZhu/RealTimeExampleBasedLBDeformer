/*
 * @file: opengl_driver.cpp
 * @brief: UI interface for simulation
 * @author: Fei Zhu
 *
 */

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include "GL/freeglut.h"
#include "GL/glui.h"
#include "camera.h"
#include "lighting.h"
#include "real_time_example_based_deformer.h"
#include "opengl_driver.h"

namespace RTLB{

OpenGLDriver* OpenGLDriver::active_instance_ = NULL;

OpenGLDriver::OpenGLDriver(const std::string &config_file_name)
{
    if(active_instance_)
        delete active_instance_;
    active_instance_ = this;
    //TO DO: init everything and enter mainloop
    initGLUT();
    initGLUI();
    initGraphics();
    initConfigurations(config_file_name);
    initSimulation();
    glutMainLoop();
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
    int argc = 0;
    char *argv = "dummy argv\n";
    glutInit(&argc, &argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH|GLUT_STENCIL);
    glutInitWindowSize(window_width_,window_height_);
    window_id_ = glutCreateWindow(window_name_.c_str());
    glutDisplayFunc(displayFunction);
    GLUI_Master.set_glutIdleFunc(idleFunction);
    GLUI_Master.set_glutKeyboardFunc(keyboardFunction);
    GLUI_Master.set_glutReshapeFunc(reshapeFunction);
    GLUI_Master.set_glutMouseFunc(mouseFunction);
    glutMotionFunc(motionFunction);
}

void OpenGLDriver::initGLUI()
{
    glui_ = GLUI_Master.create_glui("Controls",0,window_width_+52,0);
    //init GLUI control panel
    GLUI_Panel *display_mode_panel = glui_->add_panel("Display Mode",GLUI_PANEL_EMBOSSED);
    display_mode_panel->set_alignment(GLUI_ALIGN_LEFT);
    //switch render mesh
    //GLUI_RadioGroup *display_mode_radio_group = glui_->add_radiogroup_to_panel(display_mode_panel,)

    //exit button
    glui_->add_button("Exit",0,exitApplication);
    glui_->sync_live();
    glui_->set_main_gfx_window(window_id_);
}

void OpenGLDriver::initSimulation()
{
    //TO DO: setup simulation
}

void OpenGLDriver::initCamera()
{
    double MY_PI = 3.141592653589793238462643;
    double up_dir[3] = {0,1,0};
    znear_ = camera_radius_*0.01;
    zfar_ = camera_radius_*100;
    camera_ = new SphericalCamera(camera_radius_,camera_longitude_/360*2*MY_PI,camera_lattitude_/360*2*MY_PI,
                                  focus_position_,up_dir);
    camera_->SetOrigin(focus_position_);
}

void OpenGLDriver::initGraphics()
{
    // clear to white
    glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);
    // clear to light blue
    //glClearColor(233.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);
    // clear to gray
    //glClearColor(196.0 / 256, 196.0 / 256, 196.0 / 256, 0.0);
    // clear to brown
    //glClearColor(255.0 / 256, 156.0 / 256, 17.0 / 256, 0.0);
    // clear to medical blue
    //glClearColor(148.0 / 256, 199.0 / 256, 211.0 / 256, 0.0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_LINE_SMOOTH);

    reshapeFunction(window_width_,window_height_);
    //init camera
    initCamera();

    //printf ("Graphics initialization complete.\n")
}

void OpenGLDriver::displayFunction()
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    active_instance->camera_->Look();
}

void OpenGLDriver::idleFunction()
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    glutSetWindow(active_instance->window_id_);
}

void OpenGLDriver::reshapeFunction(int width, int height)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    glViewport(0,0,width,height);
    active_instance->window_width_ = width;
    active_instance->window_height_ = height;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f,1.0f*width/height,active_instance->znear_,active_instance->zfar_);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void OpenGLDriver::keyboardFunction(unsigned char key, int x, int y)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    switch(key)
    {
    case 27: //ESC, exit
        exitApplication(0);
        break;
    case 'a': //render axis
        active_instance->render_axis_ = !(active_instance->render_axis_);
        break;
    case 'w': //render wireframe
        active_instance->render_wireframe_ = !(active_instance->render_wireframe_);
        break;
    case 'b': //render fixed vertices
        active_instance->render_fixed_vertices_ = !(active_instance->render_fixed_vertices_);
        break;
    case 'r': //reset camera
        active_instance->camera_->Reset();
        break;
    case 32: //space button, pause simulation
        active_instance->pause_simulation_ = !(active_instance->pause_simulation_);
        break;
    case 'v': //render vertices
        active_instance->render_vertices_ = !(active_instance->render_vertices_);
        break;
    case 'e': //render eigenfunctions
        active_instance->render_eigenfunction_ = !(active_instance->render_eigenfunction_);
        break;
    }
}

void OpenGLDriver::specialFunction(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_LEFT:
        break;
    case GLUT_KEY_RIGHT:
        break;
    case GLUT_KEY_DOWN:
        break;
    case GLUT_KEY_UP:
        break;
    case GLUT_KEY_PAGE_UP:
        break;
    case GLUT_KEY_PAGE_DOWN:
        break;
    case GLUT_KEY_HOME:
        break;
    case GLUT_KEY_END:
        break;
    case GLUT_KEY_INSERT:
        break;
    default:
        break;
    }
}

void OpenGLDriver::motionFunction(int x, int y)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    int mouse_delta_x = x - active_instance->mouse_pos_[0];
    int mouse_delta_y = y - active_instance->mouse_pos_[1];
    active_instance->mouse_pos_[0] = x;
    active_instance->mouse_pos_[1] = y;
    SphericalCamera *camera = active_instance->camera_;
    if(active_instance->left_button_down_)
    {}
    if(active_instance->right_button_down_) //camera rotation
    {
        const double factor = 0.2;
        camera->MoveRight(factor*mouse_delta_x);
        camera->MoveUp(factor*mouse_delta_y);
    }
    if(active_instance->middle_button_down_) //zoom in/out
    {
        const double factor = 0.1;
        camera->ZoomIn(active_instance->camera_radius_*factor*mouse_delta_y);
    }
}

void OpenGLDriver::mouseFunction(int button, int state, int x, int y)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);

}

void OpenGLDriver::updateRenderMesh(int code)
{
    //TO DO
}

void OpenGLDriver::updateCurrentExample(int code)
{
    //TO DO
}

void OpenGLDriver::changeSimulationMode(int code)
{
    //TO DO
}

void OpenGLDriver::loadObjectEigenfunctions(int code)
{
    //TO DO
}

void OpenGLDriver::saveObjectEigenfunctions(int code)
{
    //TO DO
}

void OpenGLDriver::loadExampleEigenfunctions(int code)
{
    //TO DO
}

void OpenGLDriver::saveExampleEigenfunctions(int code)
{
    //TO DO
}

void OpenGLDriver::loadReducedBasis(int code)
{
    //TO DO
}

void OpenGLDriver::loadObjectCubicaData(int code)
{
    //TO DO
}

void OpenGLDriver::exitApplication(int code)
{
    //TO DO: release memories
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(active_instance->simulator_)
        delete active_instance->simulator_;
    if(active_instance->camera_)
        delete active_instance->camera_;
    if(active_instance->lighting_)
        delete active_instance->lighting_;
    exit(0);
}

void OpenGLDriver::drawAxis(double axis_length) const
{
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glBegin(GL_LINES);
    for(unsigned int i = 0; i < 3; ++i)
    {
        float color[3] = {0,0,0};
        color[i] = 1.0f;
        glColor3fv(color);
        float vertex[3] = {0,0,0};
        vertex[i] = axis_length;
        glVertex3fv(vertex);
        glVertex3f(0,0,0);
    }
    glEnd();
}

}  //namespace RTLB
