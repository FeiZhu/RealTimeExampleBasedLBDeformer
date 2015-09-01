/*
 * @file: main.cpp
 * @brief: entry of the program
 * @author: Fei Zhu
 *
 */

#include <string>
#include "opengl_driver.h"

int main()
{
    std::string config_file_name("test.config");
    RTLB::OpenGLDriver(config_file_name);
    return 0;
}
