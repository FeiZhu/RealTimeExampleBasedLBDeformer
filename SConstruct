#################################################################
# Scons script for real-time example-based materials in laplace-
# beltrami shape space project
# @author: Fei Zhu, 08/23/2015
# Usage: enter root directory of the project in terminal and
#        enter "scons"
#################################################################

######
import fnmatch
import os
from os.path import basename
import platform

#OS TYPE
os_name=platform.system()
os_architecture=platform.architecture()[0]

#BUILD TYPE
build_type='release'
#build_type='debug'

#SRC ROOT DIRECTORY
src_dir = './src/'

#GL LIBRARY
gl_include_path = src_dir+'GL/'
gl_library_path = src_dir+'GL/lib/'

#VegaFEM LIBRARY
vega_include_path = src_dir+'VegaFEM-lib/include/'
vega_library_path = src_dir+'VegaFEM-lib/lib/'

#NLopt LIBRARY
nlopt_include_path = src_dir+'nlopt/include/'
nlopt_library_path = src_dir+'nlopt/lib/'

#Opt++ LIBRARY
optpp_include_path = src_dir+'optpp/include/'
optpp_library_path = src_dir+'optpp/lib/'

#ALGLIB LIBRARY
alglib_include_path = src_dir+'alglib/'

#MKL LIBRARY
mkl_library_path = '/opt/intel/mkl/lib/intel64/'

#PATHS MADE PLATFORM SPECIFIC
if os_name=='Linux':
    vega_library_path=vega_library_path+'Linux/'
    nlopt_library_path=nlopt_library_path+'Linux/'
    optpp_library_path=optpp_library_path+'Linux/'
    gl_library_path=gl_library_path+'Linux/'
elif os_name=='Darwin':
    vega_library_path=vega_library_path+'Apple/'
    nlopt_library_path=nlopt_library_path+'Apple/'
    optpp_library_path=optpp_library_path+'Apple/'
    gl_library_path=gl_library_path+'Apple/'

if os_architecture=='32bit':
    vega_library_path=vega_library_path+'X86/'
    nlopt_library_path=nlopt_library_path+'X86/'
    optpp_library_path=optpp_library_path+'X86/'
else:
    vega_library_path=vega_library_path+'X64/'
    nlopt_library_path=nlopt_library_path+'X64/'
    optpp_library_path=optpp_library_path+'X64/'

#exampleBasedDeformableSimulator
source_filename=Glob(src_dir+'*.cpp',True,False,True)
source_filename.append(Glob(alglib_include_path+'*.cpp',True,False,True))

target_filename='RealTimeExampleBasedLBDeformer'+build_type

#LIB FILES
lib_files=[]
#VegaFEM LIBS
vega_libs = 'sceneObject integrator elasticForceModel forceModel loadList insertRows lighting performanceCounter \
             configFile renderVolumetricMesh volumetricMesh getopts camera graph isotropicHyperelasticFEM stvk corotationalLinearFEM \
             polarDecomposition massSpringSystem objMesh sparseSolver matrix sparseMatrix minivector'
lib_files.append(Split(vega_libs))
#NLOPT,OPT++ LIBS
lib_files.append('nlopt')
lib_files.append('opt')
lib_files.append('newmat')
#GL LIBS
if os_architecture=='32bit':
    lib_files.append('glui32')
else:
    lib_files.append('glui64')
lib_files.append('glut')
lib_files.append('GLU')
lib_files.append('GL')
#MKL
if os_name=='Linux':
   lib_files.append('mkl_core')

#COMPILER OPTIONS
CC='g++'
CXX='g++'
tools=['gcc', 'g++', 'gnulink']
CPPPATH=[gl_include_path,vega_include_path,alglib_include_path,nlopt_include_path,optpp_include_path]
LIBPATH=[gl_library_path,vega_library_path,nlopt_library_path,optpp_library_path]
RPATH=[gl_library_path]
LIBS=lib_files
ENV={'PATH':os.environ['PATH']}
if build_type=='release':
    CCFLAGS=['-O3','-fno-strict-aliasing','-std=gnu++11','-DNDEBUG','-DHAVE_STD','-DHAVE_NAMESPACES']
else:
    CCFLAGS=['-std=gnu++11','-fno-strict-aliasing','-g','-DHAVE_STD','-DHAVE_NAMESPACES']

env=Environment(CC=CC,CXX=CXX,tools=tools,CCFLAGS=CCFLAGS,CPPPATH=CPPPATH,LIBPATH=LIBPATH,RPATH=RPATH,LIBS=LIBS,ENV=ENV)

#BUILD
env.Program(target_filename,source_filename)
