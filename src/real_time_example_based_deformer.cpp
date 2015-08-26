/*
 * @file: real_time_example_based_deformer.cpp
 * @brief: simulation driver for real-time simulation of example-based materials in
 *         Laplace-Beltrami shape space
 * @author: Fei Zhu
 *
 */

#include "volumetricMesh.h"
#include "sceneObjectDeformable.h"
#include "planes.h"
#include "real_time_example_based_deformer.h"

RealTimeExampleBasedDeformer::RealTimeExampleBasedDeformer()
{
    //TO DO
}

RealTimeExampleBasedDeformer::~RealTimeExampleBasedDeformer()
{
    //TO DO
}

void RealTimeExampleBasedDeformer::setupSimulation()
{
    //TO DO
}

void RealTimeExampleBasedDeformer::advanceStep()
{
    //TO DO
}

void RealTimeExampleBasedDeformer::loadSimulationMesh(const std::string &file_name)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::loadExamples(const std::string &file_name_prefix, unsigned int example_num)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::loadVisualMesh(const std::string &file_name)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::loadReducedBasis(const std::string &file_name)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::loadObjectEigenfunctions(const std::string &file_name)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::loadExampleEigenFunctions(const std::string &file_name_prefix)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::loadPlanesInScene(const std::string &file_name, unsigned int plane_num)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::saveSimulationMesh(const std::string &file_name) const
{
    //TO DO
}

void RealTimeExampleBasedDeformer::saveExamples(const std::string &file_name_prefix) const
{
    //TO DO
}

void RealTimeExampleBasedDeformer::saveVisualMesh(const std::string &file_name) const
{
    //TO DO
}

void RealTimeExampleBasedDeformer::saveObjectEigenfunctions(const std::string &file_name) const
{
    //TO DO
}

void RealTimeExampleBasedDeformer::saveExampleEigenfunctions(const std::string &file_name_prefix) const
{
    //TO DO
}

void RealTimeExampleBasedDeformer::loadCorrespondenceData(const std::string &file_name)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::registerEigenfunctions()
{
    //TO DO
}

void RealTimeExampleBasedDeformer::projectOnEigenFunctions(const VolumetricMesh *mesh, const double *displacement, const double *vertex_volume,
                                                           const double **eigenfunctions, const double *eigenvalues, unsigned int eigenfunction_num,
                                                           Vec3d *eigencoefs)
{
    //TO DO
}

void RealTimeExampleBasedDeformer::reconstructFromEigenCoefs(const double **eigenfunctions, const double *eigenvalues,const Vec3d *eigencoefs,
                                                             int eigenfunction_num, int vert_num, Vec3d *vert_pos)
{
    //TO DO
}
