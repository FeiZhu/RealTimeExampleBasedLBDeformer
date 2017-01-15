/*
 * @file: main.cpp
 * @brief: entry of the program
 * @author: Fei Zhu,Mirror
 *
 */

#include <unistd.h>
#include <string>
#include "opengl_driver.h"

int main(int argc,char *argv[])
{
   //std::cout<<"aaaaaaaaaaaaaaa\n";
    int args_num=3;
    if(argc<args_num)
    {
        std::cout<<"Real-time Example-based deformable object simulator.\n";
        std::cout<<"Usage: "<<argv[0]<<" [working directory] [config file]\n";
        return 1;
    }
    if(chdir(argv[1])!=0)
    {
        std::cout<<"Error changeing working directory to: "<<argv[1]<<std::endl;
        return 1;
    }
    //parse command line option
    char *config_file_name_c=argv[2];
    opt_t opttable[]=
    {
        {NULL,0,NULL}
    };
    argv+=args_num-1;
    argc-=args_num-1;
    int optup=getopts(argc,argv,opttable);
    if(optup!=argc)
        std::cout<<"Error parsing optins. Error at option "<<argv[optup]<<".\n";
    std::cout<<"Starting application.\n";

    std::string config_file_name=std::string(config_file_name_c);
    //std::string config_file_name=std::string("examples/bar/bar.config");
    std::cout<<"Loading scene configuration from "<<config_file_name.c_str()<<".\n";
    RTLB::OpenGLDriver driver(config_file_name);
    //std::cout<<"c";
    return 0;
}
