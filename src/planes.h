/*
  Code author: Fei Zhu
  zhuf@graphics.pku.edu.cn

  A class that constructs arbitrary number of planes from a text configuration file.
  Usage: read the configuration file during initialization (using the constructor)

  The configuration file syntax supports the following parameters (listed with their types):

  Plane-specific parameters:
  --------------------------

  enable (int)
  bounce (double)
  center_X (double)
  center_Y (double)
  center_Z (double)
  normal_X (double)
  normal_Y (double)
  normal_Z (double)
  size (double)
  color_R (double)
  color_G (double)
  color_B (double)
  ambientIntensity (double)
  diffuseIntensity (double)
  specularIntensity (double)
  shininessIntensity (double)

*/

#ifndef _PLANES_H_
#define _PLANES_H_

#include <vector>
#include "vec3d.h"

#ifdef WIN32
#include <windows.h>
#endif

#include "openGL-headers.h"
using std::vector;

class ObjMesh;

class Planes
{
public:
    //read plane configurations from a configuration file
    Planes(char *config_file_name,int plane_number);
    ~Planes();

    //call this inside you OpenGL display routine, rendering the enabled planes
    void render();

    //resolve contact between the enabled planes and given mesh
    //return penalty forces, the memory for forces has be allocated outside the function
    void resolveContact(const ObjMesh *mesh/*,double *forces*/,const double *vel,double *u_new,double *vel_new);

    //set the enable status of the plane with given index (start with 1)
    //if the index is out of range, nothing happens
    void setEnableStatus(bool status,int plane_idx);

protected:
    void renderPlane(int plane_index);

protected:
    int plane_number;
    vector<int> plane_enabled;
    vector<double> plane_bounce;
    vector<Vec3d> plane_center;
    vector<Vec3d> plane_normal;
    vector<double> plane_size;
    vector<Vec3d> plane_color;
    vector<double> plane_ambient;
    vector<double> plane_diffuse;
    vector<double> plane_specular;
    vector<double> plane_shininess;
    GLuint displaylist_index;
};


#endif
