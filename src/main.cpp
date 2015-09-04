/*
 * @file: main.cpp
 * @brief: entry of the program
 * @author: Fei Zhu,Mirror
 *
 */

#include <string>
#include "opengl_driver.h"

int main()
{
    std::string config_file_name("test.config");
    RTLB::OpenGLDriver driver(config_file_name);
    return 0;
}
