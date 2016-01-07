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
#include <cfloat>
#include <map>
#include "VegaHeaders.h"
#include "optimization.h"
#include "reducedNeoHookeanForceModel.h"
#include "reducedStVKCubatureForceModel.h"
#include "rigidBody.h"
using alglib::real_1d_array;

class Planes;
class CoupledQuasiHarmonics;

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
    bool loadObjectCubicaData(const std::string &file_name);//tetID : 1-indexed
    bool loadObjectLBCubicaData(const std::string &file_name);//tetID : 1-indexed
    // bool loadObjectMass(const std::string &file_name);
    bool loadMassmatrix(const std::string &file_name);
    //bool loadExampleCubicaData(const std::string &file_name_prefix);
    bool saveSimulationMesh(const std::string &file_name) const;
    bool saveExamples(const std::string &file_name_prefix) const;
    bool saveVisualMesh(const std::string &file_name) const;
    bool saveObjectEigenfunctions(const std::string &file_name) const;
    bool saveExampleEigenfunctions(const std::string &file_name_prefix) const;

    //get && set
    // double** objectMass() const{return mass_;}
    unsigned int exampleNum() const{return example_num_;}
    void setExampleNum(int num) {example_num_=num;}
    void enableGravity(bool value){add_gravity_=value;}
    void setEnableExampleBasedSimulation(bool value){enable_example_simulation_ = value;}
    double gravity() const {return gravity_;}
    void setGravity(double gravity){ gravity_ = gravity;}
    double timeStep() const {return time_step_;}
    void setTimeStep(double dt){time_step_ = dt;}
    void setFrameRate(double frame_rate){frame_rate_ = frame_rate;}
    void setTotalFrames(int num){total_frames_ = num;}
    void setDampingMassCoef(double value){damping_mass_coef_ = value;}
    void setDampingStiffnessCoef(double value){damping_stiffness_coef_ = value;}
    void setExampleStiffnessScale(double value){example_stiffness_scale_ = value;}
    void setPrincipalStretchThreshold(double value){principal_stretch_threshold_ = value;}

    unsigned int reducedBasisNum() const{return r_;}
    double** reducedBasis(){return reduced_basis_;}
    void setInterpolateEigenfunctionNum(int num) {interpolate_eigenfunction_num_=num;}
    unsigned int correspondingFunctionNum() const{return corresponding_function_num_;}
    const Planes* planesInScnene() const {return planes_;}
    int fixedVertexNum() const{return fixed_vertex_num_;}
    int* fixedVertexPtr() const{return fixed_vertices_;}
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
    ModalMatrix* getModalmatrix() const{return modal_matrix_;}
    double* getq(){return q_;}
    double* getu(){return u_;}
    void setu(double *u);
    enum SimulationMode{
        FULLSPACE,
        REDUCEDSPACE
    };
    SimulationMode simulation_mode_ = FULLSPACE;
    enum SolverType{
        IMPLICITNEWMARK,
        IMPLICITBACKWARDEULER,
        EULER,
        CENTRALDIFFERENCES,
        REDUCEDCENTRALDIFFERENCES,
        REDUCEDIMPLICITNEWMARK,
        REDUCEDIMPLICITBACKWARDEULER,
        UNKNOWN
    };
    SolverType solver_type_ = UNKNOWN;
    enum MaterialMode{
        REDUCED_NEOHOOKEAN,
        REDUCED_STVK
    };
    MaterialMode material_mode_ = REDUCED_STVK;
    void setSimulationType(const std::string &simulation_type);
    void setMaterialType(const std::string &material_type);
    void setSolverType(std::string solver){solver_method_=solver;}
    void setEnableEigenWeightControl(bool enable_value){enable_eigen_weight_control_=enable_value;}
    void setExternalForces(double *ext_forces);
    void setReducedExternalForces(double *ext_forces);
    void setGravity(bool add_gravity,double gravity);
    // void setReducedSimulationMesh(SceneObjectReduced *mesh){reduced_simulation_mesh_=mesh};

    //registration of eigenfunctions
    bool loadCorrespondenceData(const std::string &file_name);//the vertex is 1-indexed;
    bool registerEigenfunctions();
    void testEnergyGradients();
    void testObjectiveGradients();
    void projectOnEigenFunctions(VolumetricMesh *mesh, double *displacement, double *vertex_volume,
                                double **eigenfunctions, double *eigenvalues, unsigned int eigenfunction_num,
                                Vec3d *eigencoefs);
    void reconstructFromEigenCoefs(Vec3d *target_eigencoefs,double *vert_pos,int flag=0);
    void saveReconstructMesh(double *vert_pos);
private:
    void preAllocateLocalFrameCorrespondingVertices();
    void rigidBodyPreComputation();
    void generateLocalDetailVector(const double *vert_pos);
    void generateNewDetailVector(const double *vert_pos);
    void initDisplacementMatrixOnElement(VolumetricMesh *mesh);
    //project && unproject with eigenfunctions

    void projectReducedSubspaceOnLBSubspace(double *reduced_dis, Vec3d *eigencoefs);
    void projectLBSubspaceOnReducedShpae(Vec3d *eigencoefs, double *reduced_dis);
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


    void fullspaceSimulation();
    void reducedspaceSimulation();
    Mat3d computeDs(const double *reduced_dis) const;
    Mat3d computeDmInv(const int &ele) const;
    Mat3d computeF(const int &cubica_idx,const Vec3d *reduced_dis) const;
    Mat3d firstPiolaKirchhoff(Mat3d &F) const;
    Mat3d computeF_gradient(const int &ele,const int &vert_idx,const int &vert_idx_dim) const;
    Mat3d computeP_gradient(const int &ele,const Mat3d &F,const int &vert_idx,const int &vert_idx_dim) const;

    int ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const;    //modified SVD for inversion handling
    // given a vector, find a unit vector that is orthogonal to it
    void FindOrthonormalVector(Vec3d & v, Vec3d & result) const;
    void projectOnSubBasis(VolumetricMesh *mesh,
                           double **eigenfunctions, unsigned int eigenfunction_num,
                           Vec3d *eigencoefs);

private:
    static RealTimeExampleBasedDeformer *active_instance_;
    // double **mass_ = NULL;
    //volumetric meshes
    VolumetricMesh *simulation_mesh_ = NULL;
    TetMesh *tet_mesh_=NULL;
    int simulation_vertices_num_=0;
    VolumetricMesh **examples_ = NULL;
    unsigned int example_num_ = 0;
    Graph *mesh_graph_ = NULL;
    SparseMatrix *laplacian_matrix_ = NULL;
    SparseMatrix *laplacian_damping_matrix_ = NULL;
    //visual mesh for rendering
    SceneObjectDeformable *visual_mesh_ = NULL;
    //simulation data
    double **ex_dis_ = NULL;
    double *displacement_ = NULL;
    double *velocity_ = NULL;
    double *external_force_ = NULL;
    double time_step_ = 1.0/30;
    double frame_rate_ = 30.0;
    int total_frames_ = 0;
    int time_step_counter_=0;
    int total_steps_=0;
    double epsilon_=1.0e-12;
    double integrator_epsilon_=1.0e-6;
    double last_initial_weight_=1.5;//useful when explicit weight control enable
    int fixed_vertex_num_ = 0;
    int *fixed_vertices_ = NULL;
    bool enable_eigen_weight_control_=false;
    bool pure_example_linear_interpolation_=false;
    double principal_stretch_threshold_=-DBL_MAX;
    float newmark_beta_=0.25;
    float newmark_gamma_=0.5;
    float damping_mass_coef_=0.0;
    float damping_stiffness_coef_=0.0;
    float damping_laplacian_coef_=0.0;
    double example_stiffness_scale_=1.0;//the stiffness used to compute example force is scaled
    bool add_gravity_ = false;
    double gravity_=9.8;
    int max_iterations_=1;
    int solver_threads_num_=4;//number of threads used for integration solver
    int positive_definite_=0;
    int central_difference_tangential_damping_update_mode_=1;
    // Mat3d *object_init_element_dis_matrix_ = NULL;
    // Mat3d **example_init_element_dis_matrix_ = NULL;
    int fixed_dofs_num_=0;
    int *fixed_dofs_=NULL;
    bool enable_example_simulation_=false;

    //reduced simulation data
    unsigned int r_ = 0;
    double **reduced_basis_ = NULL;
    double *reduced_basis_values_ = NULL;
    //cubica data for reduced simulation
    unsigned int object_cubica_ele_num_ = 0;
    unsigned int *object_cubica_elements_ = NULL;
    double *object_cubica_weights_ = NULL;
    //cubica data for LB space:eigenfunctions
    unsigned int object_LB_cubica_ele_num_ = 0;
    unsigned int *object_LB_cubica_elements_ = NULL;
    double *object_LB_cubica_weights_ = NULL;
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
    Vec3d *object_current_eigencoefs_ = NULL;
    Vec3d *target_eigencoefs_ = NULL;
    Vec3d *target_eigencoefs_diff_ = NULL;
    double *target_deformation_ = NULL;
    double *example_guided_deformation_ = NULL;
    double *example_based_LB_forces_ = NULL;
    double *example_based_forces_ = NULL;
    double *example_based_q_ = NULL;
    double *example_based_fq_ = NULL;
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
    // bool isPreComputeReducedData_=true;
    bool isload_cubica_ = false;
    bool isload_LB_cubica_ = false;
    bool isload_reduced_basis_ = false;
    // double **LB_object_eigenfunctions_ = NULL;
    double ***cubica_LB_tetsubBasis_ = NULL;
    //material
    double mu_=0.0;
    double lamda_=0.0;
    std::string simulation_type_;
    std::string material_type_;
    std::string solver_method_;
    SparseMatrix *mass_matrix_= NULL;
    SparseMatrixOutline *mass_matrix_outline;
    IsotropicMaterial *isotropic_material_ = NULL;
    IsotropicHyperelasticFEM *isotropic_hyperelastic_fem_ = NULL;
    ForceModel *force_model_ = NULL;
    IntegratorBase *integrator_base_ = NULL;
    ImplicitNewmarkSparse *implicit_newmark_sparse_ = NULL;
    IntegratorBaseSparse *integrator_base_sparse_ = NULL;

    ReducedForceModel *reduced_force_model_=NULL;
    ReducedNeoHookeanForceModel *reduced_neoHookean_force_model_=NULL;
    ReducedStVKCubatureForceModel *reduced_stvk_cubature_force_model_= NULL;
    IntegratorBaseDense *integrator_base_dense_ = NULL;
    ImplicitNewmarkDense *implicit_newmark_dense_ = NULL;
    ImplicitBackwardEulerDense *implicit_backward_euler_dense_ = NULL;
    CentralDifferencesDense *central_differences_dense_ = NULL;

    //fullspace Simulation
    double *u_=NULL;
    double *vel_=NULL;
    double *u_initial_=NULL;
    double *vel_initial_=NULL;
    double *f_ext_=NULL;
    double *f_col_=NULL;
    // double *full_drag_force_=NULL;
    //used for reduced cubica element Computation
    double **restpos_ = NULL;//compute rest position for cubica elements
    double *deformed_ = NULL;
    double **LB_restpos_ = NULL;//compute rest position for LB_cubica elements
    double *q_=NULL;
    double *qvel_=NULL;
    double *fq_=NULL;
    double *fqBase_=NULL;
    double *fq_ext_=NULL;
    double *reduced_mass_matrix_ = NULL;
    ModalMatrix *modal_matrix_ = NULL;
    double *U_ = NULL;
    VolumetricMesh *temp_mesh_ = NULL;
    double *examples_deformation_=NULL;//temp
    // double *reduced_drag_force_=NULL;
    // SceneObjectReduced *reduced_simulation_mesh_=NULL;
    // SceneObjectReducedCPU *reduced_simulation_mesh_cpu_=NULL;

    //temp mesh_eigen_skeleton
    double *initial_eigen_skeleton_=NULL;
	double *initial_detail_vector_=NULL;
    double *local_detail_vector_=NULL;
    double *deformed_detail_vector_=NULL;
    std::map<int,int> vert_vertex1,vert_vertex2;
    double mass_=0.0;
    Vec3d rigid_center_;
};

} //namespace RTLB

#endif //REAL_TIME_EXAMPLE_BASED_DEFORMER_H_
