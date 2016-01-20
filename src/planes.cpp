#include <cstdio>
#include <cstring>
#include "configFile.h"
#include "objMesh.h"
#include "planes.h"
#include "volumetricMesh.h"

Planes::Planes(const char* config_file_name,unsigned int plane_number):plane_number(plane_number)
{
    ConfigFile config_file;

    plane_enabled.resize(plane_number);
    plane_bounce.resize(plane_number);
    plane_center.resize(plane_number);
    plane_normal.resize(plane_number);
    plane_size.resize(plane_number);
    plane_color.resize(plane_number);
    plane_ambient.resize(plane_number);
    plane_diffuse.resize(plane_number);
    plane_specular.resize(plane_number);
    plane_shininess.resize(plane_number);

    for(int plane_index=0;plane_index<plane_number;++plane_index)
    {
	char option_name[128];
	char plane_index_ch='0'+plane_index;

	sprintf(option_name,"planeEnabled_%c",plane_index_ch);
	config_file.addOption(option_name,&plane_enabled[plane_index]);
	sprintf(option_name,"planeBounce_%c",plane_index_ch);
	config_file.addOption(option_name,&plane_bounce[plane_index]);
	sprintf(option_name,"planeCenter_%c_X",plane_index_ch);
	config_file.addOption(option_name,&plane_center[plane_index][0]);
	sprintf(option_name,"planeCenter_%c_Y",plane_index_ch);
	config_file.addOption(option_name,&plane_center[plane_index][1]);
	sprintf(option_name,"planeCenter_%c_Z",plane_index_ch);
	config_file.addOption(option_name,&plane_center[plane_index][2]);
	sprintf(option_name,"planeNormal_%c_X",plane_index_ch);
	config_file.addOption(option_name,&plane_normal[plane_index][0]);
	sprintf(option_name,"planeNormal_%c_Y",plane_index_ch);
	config_file.addOption(option_name,&plane_normal[plane_index][1]);
	sprintf(option_name,"planeNormal_%c_Z",plane_index_ch);
	config_file.addOption(option_name,&plane_normal[plane_index][2]);
	sprintf(option_name,"planeSize_%c",plane_index_ch);
	config_file.addOption(option_name,&plane_size[plane_index]);
	sprintf(option_name,"planeColor_%c_X",plane_index_ch);
	config_file.addOption(option_name,&plane_color[plane_index][0]);
	sprintf(option_name,"planeColor_%c_Y",plane_index_ch);
	config_file.addOption(option_name,&plane_color[plane_index][1]);
	sprintf(option_name,"planeColor_%c_Z",plane_index_ch);
	config_file.addOption(option_name,&plane_color[plane_index][2]);
	sprintf(option_name,"planeAmbient_%c",plane_index_ch);
	config_file.addOption(option_name,&plane_ambient[plane_index]);
	sprintf(option_name,"planeDiffuse_%c",plane_index_ch);
	config_file.addOption(option_name,&plane_diffuse[plane_index]);
	sprintf(option_name,"planeSpecular_%c",plane_index_ch);
	config_file.addOption(option_name,&plane_specular[plane_index]);
	sprintf(option_name,"planeShininess_%c",plane_index_ch);
	config_file.addOption(option_name,&plane_shininess[plane_index]);
    }
    printf("Parsing planes from %s...\n",config_file_name);
    if(config_file.parseOptions(config_file_name)!=0)
	   throw 1;
    config_file.printOptions();
    //build display list for the enabled planes
    displaylist_index=glGenLists(1);
    if(displaylist_index!=0)
    {
    	glNewList(displaylist_index,GL_COMPILE);
    	for(int plane_index=0;plane_index<plane_number;++plane_index)
    	{
    	    if(plane_enabled[plane_index])
    		renderPlane(plane_index);
    	}
    	glEndList();
    }
    else
	   printf("Warning: Failed to create display list for planes.\n");
}

Planes::~Planes()
{
    glDeleteLists(displaylist_index,1);
}

void Planes::render()
{
    if(displaylist_index)
	   glCallList(displaylist_index);
    else
    {
    	for(int plane_index=0;plane_index<plane_number;++plane_index)
    	{
    	    if(plane_enabled[plane_index])
    		renderPlane(plane_index);
    	}
    }
}
void Planes::resolveContact(ObjMesh *mesh,double *forces)
{
    // std::cout<<".............BEGIN.\n";
    int vert_num=mesh->getNumVertices();
    double threshold=1.0/*mesh->getDiameter()*/;//the threshold to start resolve contact
    memset(forces,0.0,sizeof(double)*3*vert_num);
    for(int plane_index=0;plane_index<plane_number;++plane_index)
    {
    	if(!plane_enabled[plane_index])
    	    continue;
    	Vec3d unit_plane_normal=norm(plane_normal[plane_index]);//normalize the plane normal
        for(int vert_index=0;vert_index<vert_num;++vert_index)
    	{
    	    Vec3d rel_vec=mesh->getPosition(vert_index)-plane_center[plane_index];
        //     if(vert_index==0)
        // std::cout<<"mesh->getPosition(vert_index):"<<mesh->getPosition(vert_index)<<","<<plane_center[plane_index]<<","<<rel_vec[1]<<"\n";
        // getchar();
    	    double dist_vec=dot(rel_vec,unit_plane_normal);
            // std::cout<<dist_vec<<",";
        // std::cout<<dist_vec<<"\n";
    	    if(dist_vec<threshold)//close than a threshold or penetrated
    	    {
                // std::cout<<"kkkkkkkkkkk\n";
        		if(dist_vec<0)
        		    dist_vec=-dist_vec;
                // std::cout<<dist_vec<<"\n";
        		forces[3*vert_index+0]+=plane_bounce[plane_index]*dist_vec*unit_plane_normal[0];
        		forces[3*vert_index+1]+=plane_bounce[plane_index]*dist_vec*unit_plane_normal[1];
        		forces[3*vert_index+2]+=plane_bounce[plane_index]*dist_vec*unit_plane_normal[2];
            // std::cout<<forces[3*vert_index+1]<<"\n";
    	    }
    	}
    }
    // std::cout<<".............END.\n";
}
void Planes::resolveContact(const ObjMesh *mesh/*,double *forces*/,const double *vel,double *u_new,double *vel_new)
{
    int vert_num=mesh->getNumVertices();
    double threshold=0.03*mesh->getDiameter();//the threshold to start resolve contact

	//memset(forces,0.0,sizeof(double)*3*vert_num);
	memset(u_new,0.0,sizeof(double)*3*vert_num);
	memset(vel_new,0.0,sizeof(double)*3*vert_num);
    for(int plane_index=0;plane_index<plane_number;++plane_index)
    {
		if(!plane_enabled[plane_index])
			continue;
		Vec3d unit_plane_normal=norm(plane_normal[plane_index]);//normalize the plane normal
		for(int vert_index=0;vert_index<vert_num;++vert_index)
		{
			Vec3d rel_vec=mesh->getPosition(vert_index)-plane_center[plane_index];
			double dist_vec=dot(rel_vec,unit_plane_normal);
		//	Vec3d vert_pos=mesh->getPosition(vert_index);
			Vec3d vert_vel(0.0);
			for(unsigned int i=0;i<3;++i)
			{
				vert_vel[i]=vel[3*vert_index+i];
			}
			double vel_dot = dot(vert_vel,unit_plane_normal); 
			if((vel_dot<0)&&(dist_vec<threshold))
			{
				//v_new = v - dot(v,n)n
                                for(int i=0;i<3;++i)
				    vel_new[3*vert_index+i]=vert_vel[i]- vel_dot*unit_plane_normal[i];
                                std::cout<<"....\n";
			}
			//change x to handle penetrated
			if(dist_vec<0)
			{
				for(unsigned int i = 0; i < 3; ++i)
				    u_new[3*vert_index+i]=-dist_vec*unit_plane_normal[i];
			}

		}
    }
}

void Planes::setEnableStatus(bool status,int plane_idx)
{
    if(plane_idx<=0||plane_idx>plane_number)
	return;
    plane_enabled[plane_idx-1]=status?1:0;
    //update display list
    glDeleteLists(displaylist_index,1);
    //build display list for the enabled planes
    displaylist_index=glGenLists(1);
    if(displaylist_index!=0)
    {
	glNewList(displaylist_index,GL_COMPILE);
	for(int plane_index=0;plane_index<plane_number;++plane_index)
	{
	    if(plane_enabled[plane_index])
		renderPlane(plane_index);
	}
	glEndList();
    }
}

void Planes::renderPlane(int plane_index)
{
    Vec3d current_center=plane_center[plane_index];
    Vec3d current_normal=plane_normal[plane_index];
    Vec3d current_color=plane_color[plane_index];
    double current_size=plane_size[plane_index];
    double current_ambient=plane_ambient[plane_index];
    double current_diffuse=plane_diffuse[plane_index];
    double current_specular=plane_specular[plane_index];
    double current_shininess=plane_shininess[plane_index];

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0,1.0);

    GLfloat ambient[4]={GLfloat(current_ambient*current_color[0]),GLfloat(current_ambient*current_color[1]),GLfloat(current_ambient*current_color[2]),1.0f};
    GLfloat diffuse[4]={GLfloat(current_diffuse*current_color[0]),GLfloat(current_diffuse*current_color[1]),GLfloat(current_diffuse*current_color[2]),1.0f};
    GLfloat specular[4]={GLfloat(current_specular*current_color[0]),GLfloat(current_specular*current_color[1]),GLfloat(current_specular*current_color[2]),1.0f};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, current_shininess);

    //the transformation between the plane and the x-z plane
    //transform x-z plane to the plane
    float rotate_around_y[16],rotate_around_z[16],translate[16];
    memset(rotate_around_y,0.0f,sizeof(float)*16);
    memset(rotate_around_z,0.0f,sizeof(float)*16);
    memset(translate,0.0f,sizeof(float)*16);

    translate[0]=translate[5]=translate[10]=translate[15]=1.0f;
    translate[12]=current_center[0];
    translate[13]=current_center[1];
    translate[14]=current_center[2];

    rotate_around_y[5]=rotate_around_y[15]=1.0f;
    Vec3d normal_xz=Vec3d(current_normal[0],0.0,current_normal[2]);
    Vec3d x_axis=Vec3d(1.0,0.0,0.0);
    double normal_xz_len=len(normal_xz);
    double cos_y,sin_y;
    if(normal_xz_len>0)
    {
        cos_y=dot(normal_xz,x_axis)/normal_xz_len;
        sin_y=-normal_xz[2]/normal_xz_len;
    }
    else
    {
        cos_y=1.0;
        sin_y=0.0;
    }
    rotate_around_y[0]=rotate_around_y[10]=cos_y;
    rotate_around_y[2]=-sin_y;
    rotate_around_y[8]=sin_y;

    rotate_around_z[10]=rotate_around_z[15]=1.0f;
    Vec3d normal_xy=Vec3d(normal_xz_len,current_normal[1],0.0);
    double normal_xy_len=len(normal_xy);
    double cos_z,sin_z;
    if(normal_xy_len>0)
    {
        cos_z=normal_xy[1]/normal_xy_len;
        sin_z=-normal_xy[0]/normal_xy_len;
    }
    else
    {
        cos_z=1.0;
        sin_z=0.0;
    }
    rotate_around_z[0]=rotate_around_z[5]=cos_z;
    rotate_around_z[1]=sin_z;
    rotate_around_z[4]=-sin_z;

    const int plane_resolution=1;
    double plane_increment=current_size/plane_resolution;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixf(translate);
    glMultMatrixf(rotate_around_z);
    glMultMatrixf(rotate_around_y);
    glNormal3f(current_normal[0],current_normal[1],current_normal[2]);
    for(int i=0;i<plane_resolution;++i)
        for(int j=0;j<plane_resolution;++j)
        {
            glBegin(GL_QUADS);
            glVertex3f(-current_size/2.0+i*plane_increment,0.0,-current_size/2.0+j*plane_increment);
            glVertex3f(-current_size/2.0+i*plane_increment,0.0,-current_size/2.0+(j+1)*plane_increment);
            glVertex3f(-current_size/2.0+(i+1)*plane_increment,0.0,-current_size/2.0+(j+1)*plane_increment);
            glVertex3f(-current_size/2.0+(i+1)*plane_increment,0.0,-current_size/2.0+j*plane_increment);
            glEnd();
        }
    glPopMatrix();
    glDisable(GL_POLYGON_OFFSET_FILL);
}
