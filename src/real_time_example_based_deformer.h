/*
 * @file: real_time_example_based_deformer.h
 * @brief: simulation driver for real-time simulation of example-based materials in
 *         Laplace-Beltrami shape space
 * @author: Fei Zhu
 *
 */

#ifndef REAL_TIME_EXAMPLE_BASED_DEFORMER_H_
#define REAL_TIME_EXAMPLE_BASED_DEFORMER_H_

#include <string>

class Vec3d;
class VolumetricMesh;
class SceneObjectDeformable;
class Planes;
class CoupledQuasiHarmonics;

class RealTimeExampleBasedDeformer
{
public:
    RealTimeExampleBasedDeformer();
    ~RealTimeExampleBasedDeformer();

    void setupSimulation();  //preprocess
    void advanceStep();  //one time step

    //basic load && save
    bool loadSimulationMesh(const std::string &file_name);//.veg
    bool loadExamples(const std::string &file_name_prefix, unsigned int example_num);//.veg
    bool loadVisualMesh(const std::string &file_name);//.obj
    bool loadReducedBasis(const std::string &file_name);//.basis
    bool loadObjectEigenfunctions(const std::string &file_name);//.eigen
    bool loadExampleEigenFunctions(const std::string &file_name_prefix);//.eigen
    bool loadPlanesInScene(const std::string &file_name, unsigned int plane_num);
    bool saveSimulationMesh(const std::string &file_name) const;
    bool saveExamples(const std::string &file_name_prefix) const;
    bool saveVisualMesh(const std::string &file_name) const;
    bool saveObjectEigenfunctions(const std::string &file_name) const;
    bool saveExampleEigenfunctions(const std::string &file_name_prefix) const;

    //get && set
    unsigned int exampleNum() const{return example_num_;}
    double gravity() const {return gravity_;}
    void setGravity(double gravity){ gravity_ = gravity;}
    double timeStep() const {return time_step_;}
    void setTimeStep(double dt){time_step_ = dt;}
    unsigned int reducedBasisNum() const{return reduced_basis_num_;}
    unsigned int objectEigenfunctionNum() const{return object_eigenfunction_num_;}
    unsigned int exampleEigenfunctionNum() const{return example_eigenfunction_num_;}
    const Planes* planesInScnene() const {return planes_;}

    //registration of eigenfunctions
    bool loadCorrespondenceData(const std::string &file_name);
    bool registerEigenfunctions();

private:
    //project && unproject with eigenfunctions
    void projectOnEigenFunctions(const VolumetricMesh *mesh, const double *displacement, const double *vertex_volume,
                                 const double **eigenfunctions, const double *eigenvalues, unsigned int eigenfunction_num,
                                 Vec3d *eigencoefs);
    void reconstructFromEigenCoefs(const double **eigenfunctions, const double *eigenvalues,const Vec3d *eigencoefs,
                                   int eigenfunction_num, int vert_num, Vec3d *vert_pos);

private:
    //volumetric meshes
    VolumetricMesh *simulation_mesh_ = NULL;
    VolumetricMesh **examples_ = NULL;
    unsigned int example_num_ = 0;
    //visual mesh for rendering
    SceneObjectDeformable *visual_mesh_ = NULL;
    //simulation data
    double *displacement_ = NULL;
    double *velocity_ = NULL;
    double *external_force_ = NULL;
    double gravity_ = -9.8;
    double time_step_ = 1.0/30;
    //reduced simulation data
    unsigned int reduced_basis_num_ = 0;
    double **reduced_basis_ = NULL;
    double *reduced_displacement_ = NULL;
    double *reduced_velocity_ = NULL;
    //eigenfunction data
    double **object_eigenfunctions_ = NULL;
    double *object_eigenvalues_ = NULL;
    unsigned int object_eigenfunction_num_ = 0;
    double ***example_eigenfunctions_ = NULL;
    double **example_eigenvalues_ = NULL;
    unsigned int example_eigenfunction_num_ = 0;
    //shapes in LB shape space
    Vec3d *object_eigencoefs_ = NULL;
    Vec3d **example_eigencoefs_ = NULL;
    Vec3d *target_eigencoefs_ = NULL;
    //total volume and per-vertex volume
    double object_volume_ = 0;
    double *object_vertex_volume_ = NULL;
    double *example_volume_ = NULL;
    double **example_vertex_volume_ = NULL;
    //for eigenfunction registration
    double **object_corresponding_functions_ = NULL;
    double ***example_corresponding_functions_ = NULL;
    unsigned int corresponding_function_num_ = 0;
    bool is_region_based_correspondence_ = false;
    //planes in scene, for contact
    Planes *planes_ = NULL;
    unsigned int plane_num_ = 0;
};

#endif //REAL_TIME_EXAMPLE_BASED_DEFORMER_H_
