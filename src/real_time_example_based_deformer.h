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
#include "optimization.h"
#include "mat3d.h"
using alglib::real_1d_array;

class Vec3d;
class VolumetricMesh;
class SceneObjectDeformable;
class Planes;
class CoupledQuasiHarmonics;
class ConfigFile;

namespace RTLB{

class RealTimeExampleBasedDeformer
{
public:
    RealTimeExampleBasedDeformer();
    ~RealTimeExampleBasedDeformer();
private:
    static RealTimeExampleBasedDeformer* activeInstance();
public:
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
    bool loadFixedVertices(const std::string &file_name);
    bool loadObjectCubicaData(const std::string &file_name);//tetID : 0-indexed
    bool loadExampleCubicaData(const std::string &file_name_prefix);
    bool saveSimulationMesh(const std::string &file_name) const;
    bool saveExamples(const std::string &file_name_prefix) const;
    bool saveVisualMesh(const std::string &file_name) const;
    bool saveObjectEigenfunctions(const std::string &file_name) const;
    bool saveExampleEigenfunctions(const std::string &file_name_prefix) const;

    //get && set
    unsigned int exampleNum() const{return example_num_;}
    void setExampleNum(int num) {example_num_=num;}
    double gravity() const {return gravity_;}
    void setGravity(double gravity){ gravity_ = gravity;}
    double timeStep() const {return time_step_;}
    void setTimeStep(double dt){time_step_ = dt;}
    unsigned int reducedBasisNum() const{return reduced_basis_num_;}
    unsigned int objectEigenfunctionNum() const{return object_eigenfunction_num_;}
    unsigned int exampleEigenfunctionNum() const{return example_eigenfunction_num_;}
    unsigned int correspondingFunctionNum() const{return corresponding_function_num_;}
    const Planes* planesInScnene() const {return planes_;}
    unsigned int fixedVertexNum() const{return fixed_vertex_num_;}
    unsigned int* fixedVertexPtr() const{return fixed_vertices_;}
    VolumetricMesh* simulationMesh() const{return simulation_mesh_;}
    VolumetricMesh* exampleMesh(unsigned int example_idx) const;
    const SceneObjectDeformable* visualMesh() const{return visual_mesh_;}
    double* objectVertexVolume() const{return object_vertex_volume_;}
    double** exampleVertexVolume() const{return example_vertex_volume_;}
    double** objectEigenFunctions() const{return object_eigenfunctions_;}
    double* objectEigenValues() const{return object_eigenvalues_;}
    double*** exampleEigenFunctions() const{return example_eigenfunctions_;}
    double** exampleEigenValues() const{return example_eigenvalues_;}
    Vec3d* objectEigencoefs() const{return object_eigencoefs_;}
    Vec3d** exampleEigencoefs() const{return example_eigencoefs_;}
    unsigned int objectCubicaEleNum() const{return object_cubica_ele_num_;}
    unsigned int* objectCubicaElements() const{return object_cubica_elements_;}
    double* objectCubicaWeights() const{return object_cubica_weights_;}
    unsigned int* exampleCubicaEleNum() const{return example_cubica_ele_num_;}
    unsigned int** exampleCubicaElements() const{return example_cubica_elements_;}
    double** exampleCubicaWeights() const{return example_cubica_weights_;}
    void setEnableEigenWeightControl(bool enable_value){enable_eigen_weight_control_=enable_value;}

    //registration of eigenfunctions
    bool loadCorrespondenceData(const std::string &file_name);//the vertex is 1-indexed;
    bool registerEigenfunctions();

    void initDisplacementMatrixOnElement(VolumetricMesh *mesh);
    //project && unproject with eigenfunctions
    void projectOnEigenFunctions(VolumetricMesh *mesh, double *displacement, double *vertex_volume,
                                 double **eigenfunctions, double *eigenvalues, unsigned int eigenfunction_num,
                                 Vec3d *eigencoefs);
    void reconstructFromEigenCoefs(double **eigenfunctions, double *eigenvalues,const Vec3d *eigencoefs,
                                   int eigenfunction_num, int vert_num, double *vert_pos);
    void projectOnExampleManifold(Vec3d *object_eigencoefs, Vec3d *target_eigencoefs);
    static void evaluateObjectiveAndGradient1(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr);
    static void evaluateObjectiveAndGradient2(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr);
    Mat3d computeDeformGradient(const Mat3d &init_matrix,const Mat3d &deformed_matrix/*Vec3d *init_pos,Vec3d *deform_pos*/);
    //compute the deformed energy from source pos to deformed pos
    //we use optimizing cubature method to compute energy on reduced space
    //solving course need the initial displacement matrix Dm, and the deformed displacement matrix Ds
    // to consider solving time we compute the object and examples initial displacement Dm when we load volumetric meshes
    //the last two parameters is to identify the energy we solved is for example mesh deformation or object deformation
    void computeReducedEnergyAndGradient(VolumetricMesh *mesh,const double *init_pos,const double *displacement, const unsigned int cubica_num,
    									const unsigned int *cubica_elements, const double *cubica_weights,const unsigned int example_flag,
    									const unsigned int dis_ex_idx,double &energy,double *grad);
    void computeForceOnReducedSubSpace(VolumetricMesh *mesh,const double *init_pos,const double *displacement,const unsigned int example_flag,
    	 								const unsigned int dis_ex_idx,double *g);
    void computeForceForSingleElementOnReducedSubSpace(const VolumetricMesh *mesh,const unsigned int ele,const double *ele_dis,unsigned int example_flag,
						                const unsigned int dis_ex_idx,const Mat3d &H,double *g);
private:
    int ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const;    //modified SVD for inversion handling
    // given a vector, find a unit vector that is orthogonal to it
    void FindOrthonormalVector(Vec3d & v, Vec3d & result) const;
private:
    static RealTimeExampleBasedDeformer *active_instance_;
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
    double epsilon_=1.0e-12;
    double integrator_epsilon_=1.0e-6;
    double last_initial_weight_=1.5;//useful when explicit weight control enable
    unsigned int fixed_vertex_num_ = 0;
    unsigned int *fixed_vertices_ = NULL;
    bool enable_eigen_weight_control_=false;
    bool pure_example_linear_interpolation_=false;
    Mat3d *object_init_element_dis_matrix_ = NULL;
    Mat3d **example_init_element_dis_matrix_ = NULL;

    //reduced simulation data
    unsigned int reduced_basis_num_ = 0;
    double **reduced_basis_ = NULL;
    double *reduced_displacement_ = NULL;
    double *reduced_velocity_ = NULL;
    //cubica data
    unsigned int object_cubica_ele_num_ = 0;
    unsigned int *object_cubica_elements_ = NULL;
    double *object_cubica_weights_ = NULL;
    unsigned int *example_cubica_ele_num_ = NULL;
    unsigned int **example_cubica_elements_ = NULL;
    double **example_cubica_weights_ = NULL;
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

} //namespace RTLB

#endif //REAL_TIME_EXAMPLE_BASED_DEFORMER_H_
