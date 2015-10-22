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
//#include "NLF.h"
#include "matrix.h"
using alglib::real_1d_array;
//the NEWMAT::Matrix<double> index from 1
//using NEWMAT::Matrix<double>;

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
    //bool loadExampleCubicaData(const std::string &file_name_prefix);
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
    double** reducedBasis() const{return reduced_basis_;}
    void setInterpolateEigenfunctionNum(int num) {interpolate_eigenfunction_num_=num;}
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
    void setEnableEigenWeightControl(bool enable_value){enable_eigen_weight_control_=enable_value;}

    //registration of eigenfunctions
    bool loadCorrespondenceData(const std::string &file_name);//the vertex is 1-indexed;
    bool registerEigenfunctions();

    void initDisplacementMatrixOnElement(VolumetricMesh *mesh);
    //project && unproject with eigenfunctions
    void projectOnEigenFunctions(VolumetricMesh *mesh, double *displacement, double *vertex_volume,
                                 double **eigenfunctions, double *eigenvalues, unsigned int eigenfunction_num,
                                 Vec3d *eigencoefs);
    void reconstructFromEigenCoefs(double **eigenfunctions,double *eigenvalues,Vec3d *eigencoefs,Vec3d *target_eigencoefs,
                                 const int &eigenfunction_num,const int &input_reconstruct_eigenfunction_num,const int &vert_num, double *vert_pos);
    void projectOnExampleManifold(Vec3d *object_eigencoefs, Vec3d *target_eigencoefs);
    static void evaluateObjectiveAndGradient1(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr);
    static void evaluateObjectiveAndGradient2(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr);
    Mat3d computeDeformationGradient(const Mat3d &init_matrix,const Mat3d &deformed_matrix/*Vec3d *init_pos,Vec3d *deform_pos*/);
    //compute the deformed energy from source pos to deformed pos
    //we use optimizing cubature method to compute energy on reduced space
    //solving course need the initial displacement matrix Dm, and the deformed displacement matrix Ds
    // to consider solving time we compute the object and examples initial displacement Dm when we load volumetric meshes
    //the last two parameters is to identify the energy we solved is for example mesh deformation or object deformation
    void computeReducedEnergy(const Vec3d *reduced_dis,double &energy) const;
    void computeReducedInternalForce(const Vec3d *reduced_dis,double *forces) const;
    void computeReducedStiffnessMatrix(const Vec3d *reduced_dis,Matrix<double> &reduced_K) const;
    void testObjectiveGradients();
private:
    void preComputeForReducedSimulation();
    Matrix<double> vertexSubBasis(const int &vert_idx) const;//3*r
    Matrix<double> tetSubBasis(const int &ele) const;//12*r
    Mat3d computeDs(const double *reduced_dis) const;
    Mat3d computeDmInv(const int &ele) const;
    void computeF(const Vec3d *reduced_dis) const;
    Mat3d firstPiolaKirchhoff(Mat3d &F) const;
    Mat3d computeF_gradient(const int &ele,const int &vert_idx,const int &vert_idx_dim) const;
    Mat3d computeP_gradient(const int &ele,const Mat3d F,const int &vert_idx,const int &vert_idx_dim) const;

    int ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const;    //modified SVD for inversion handling
    // given a vector, find a unit vector that is orthogonal to it
    void FindOrthonormalVector(Vec3d & v, Vec3d & result) const;
    void projectOnSubBasis(VolumetricMesh *mesh,
                           double **eigenfunctions, unsigned int eigenfunction_num,
                           Vec3d *eigencoefs);
private:
    static RealTimeExampleBasedDeformer *active_instance_;
    //volumetric meshes
    VolumetricMesh *simulation_mesh_ = NULL;
    VolumetricMesh **examples_ = NULL;
    unsigned int example_num_ = 0;
    //visual mesh for rendering
    SceneObjectDeformable *visual_mesh_ = NULL;
    //simulation data
    double **ex_dis_ = NULL;
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
    double *reduced_basis_values_ = NULL;
    double *reduced_displacement_ = NULL;
    double *reduced_velocity_ = NULL;
    //cubica data
    unsigned int object_cubica_ele_num_ = 0;
    unsigned int *object_cubica_elements_ = NULL;
    double *object_cubica_weights_ = NULL;
    //eigenfunction data
    double **object_eigenfunctions_ = NULL;
    double *object_eigenvalues_ = NULL;
    unsigned int interpolate_eigenfunction_num_ = 0;
    double ***example_eigenfunctions_ = NULL;
    double **example_eigenvalues_ = NULL;
    unsigned int reconstruct_eigenfunction_num_ = 0;
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
    bool isPreComputeReducedData_=true;
    //used for reduced cubica element Computation
    Vec3d *q_;
    Mat3d *F_;
    double **restpos_;//compute rest position for cubica elements
    //material
    double mu_=0.0;
    double lamda_=0.0;
};

} //namespace RTLB

#endif //REAL_TIME_EXAMPLE_BASED_DEFORMER_H_
