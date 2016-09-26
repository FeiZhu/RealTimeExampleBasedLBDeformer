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
#include "rigidBody_generalTensor.h"

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
    bool loadInertiaTensor(const std::string &file_name);
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
    bool getExampleSimulationState(){return enable_example_simulation_;}
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
    void setDampingExampleStiffness(double value){damping_example_stiffness_ = value;}
    void setExampleBiasFactor(double value){example_bias_ = value;}
    void setPrincipalStretchThreshold(double value){principal_stretch_threshold_ = value;}
    void setCollisionNumLimited(int value){col_limited_num_=value;}
    void setu(double *u);
    void setVelAfterCollision(double *vel);
    void setCollisionNum(int num){collide_vert_num_=num;}
    void setRigidInitialVel(double x,double y,double z){initial_rigid_vel_[0]=x;initial_rigid_vel_[1]=y;initial_rigid_vel_[2]=z;}
    void setInitialForceFilename(const std::string &file_name){force_loads_file_name_=file_name;}
    void setInitialVelFilename(const std::string &file_name){initial_velocity_file_name_=file_name;}
    void setInitialPosFilename(const std::string &file_name){initial_position_file_name_=file_name;}
    void setInitialTetMeshFilename(const std::string &file_name){initial_tetmesh_file_name_=file_name;}
    void setObjectAffectedVerticesFilename(const std::string &file_name){object_affected_vertices_file_name_=file_name;}
    void setExampleAffectedVerticesFilebase(const std::string &file_name){example_affected_vertices_file_base_=file_name;}
    void setTorqueCoef(double value){torque_coef_=value;}
    // void setInitialVel(double *vel);
    // void setInitialDis(double *pos);
    // void setInitialForceLoad(double *force);
    // void setInitialForceLoadNum(int num){num_force_loads_=num;}

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
    double* getvel(){return vel_;}
    double* getExternalForce(){return f_ext_;}
    RigidBody_GeneralTensor* getRigid() const{return rigid_;}
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
    SolverType solver_type_=IMPLICITNEWMARK;
    enum MaterialMode{
        REDUCED_NEOHOOKEAN,
        REDUCED_STVK,
        INV_STVK,
        INV_NEOHOOKEAN
    };
    MaterialMode material_mode_ = INV_STVK;
    void setSimulationType(const std::string &simulation_type);
    void setMaterialType(const std::string &material_type);
    void setSolverType(const std::string &solver);
    void setEnableEigenWeightControl(bool enable_value){enable_eigen_weight_control_=enable_value;}
    void setExternalForces(double *ext_forces);
    void setReducedExternalForces(double *ext_forces);
    void setGravity(bool add_gravity,double gravity);
    void setConstrains(bool with_constrains){with_constrains_=with_constrains;}
    // void setReducedSimulationMesh(SceneObjectReduced *mesh){reduced_simulation_mesh_=mesh};

    //registration of eigenfunctions
    bool loadCorrespondenceData(const std::string &file_name);//the vertex is 1-indexed;
    bool registerEigenfunctions();
    void testEnergyGradients();
    void testObjectiveGradients();
    void projectOnEigenFunctions(VolumetricMesh *mesh, double *displacement, double *vertex_volume,
                                double **eigenfunctions, double *eigenvalues, unsigned int eigenfunction_num,
                                Vec3d *eigencoefs,int affected_vertices_num=0,int *affected_vertices=NULL);
    void projectOnEigenFunctions1(double *displacement, double *vertex_volume,
                                double **eigenfunctions, double *eigenvalues, unsigned int eigenfunction_num,Vec3d *eigencoefs);
    void reconstructFromEigenCoefs(Vec3d *target_eigencoefs,double *vert_pos);//full space
    void reconstructFromEigenCoefs(Vec3d *target_eigencoefs,int flag=0);//reduced_space
    void saveReconstructMesh(double *vert_pos);
    void saveReconstructMesh1(double *vert_pos);
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
    static void evaluateObjectiveAndGradient3(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr);
    Mat3d computeDeformationGradient(const Mat3d &init_matrix,const Mat3d &deformed_matrix/*Vec3d *init_pos,Vec3d *deform_pos*/);
    //compute the deformed energy from source pos to deformed pos
    //we use optimizing cubature method to compute energy on reduced space
    //solving course need the initial displacement matrix Dm, and the deformed displacement matrix Ds
    // to consider solving time we compute the object and examples initial displacement Dm when we load volumetric meshes
    //the last two parameters is to identify the energy we solved is for example mesh deformation or object deformation

    void computeEnergy(const Vec3d *reduced_dis,double &energy) const;
    void computeEnergyGrad(const Vec3d *reduced_dis,double *forces) const;
    void computeReducedEnergy(const Vec3d *reduced_dis,double &energy) const;
    void computeReducedInternalForce(const Vec3d *reduced_dis,double *forces) const;
    void computeReducedStiffnessMatrix(const Vec3d *reduced_dis,Matrix<double> &reduced_K) const;


    void fullspaceSimulation();
    void reducedspaceSimulationWithConstraints();
    void reducedspaceSimulationWithoutConstraints();
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
    double damping_example_stiffness_=0.0;//damping stiffness force coefficients
    double example_bias_ = 0.8;
    bool add_gravity_ = true;
    double gravity_=9.8;
    int max_iterations_=1;
    int solver_threads_num_=4;//number of threads used for integration solver
    int positive_definite_=0;
    int central_difference_tangential_damping_update_mode_=1;
    double torque_coef_=1.0e-5;
    int col_limited_num_=20;
    // Mat3d *object_init_element_dis_matrix_ = NULL;
    // Mat3d **example_init_element_dis_matrix_ = NULL;
    int fixed_dofs_num_=0;
    int *fixed_dofs_=NULL;
    bool enable_example_simulation_=false;
    std::string force_loads_file_name_="none";
    std::string initial_velocity_file_name_="none";
    std::string initial_position_file_name_="none";
    std::string initial_tetmesh_file_name_="none";
    int force_loads_num_=0;
    double *force_loads_=NULL;

    std::string object_affected_vertices_file_name_="none";
    std::string example_affected_vertices_file_base_="none";
    int object_affected_vertices_num_=0;
    int *example_affected_vertices_num_=NULL;
    int *object_affected_vertices_=NULL;
    int **example_affected_vertices_=NULL;

    //reduced simulation data
    unsigned int r_ = 0;
    double **reduced_basis_ = NULL;
    double *reduced_basis_values_ = NULL;
    double **LB_to_reducespace_ = NULL;//dimension:rxm(r:reduced basis num,m:interpolate eigenfunction num)
    double **reducespace_to_LB_ = NULL;//dimension:rxm(r:reduced basis num,m:interpolate eigenfunction num)
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
    Vec3d example_length_;
    Vec3d **object_example_eigencoefs_=NULL;
    Vec3d *object_current_eigencoefs_ = NULL;
    Vec3d *target_eigencoefs_ = NULL;
    Vec3d *target_eigencoefs_diff_ = NULL;
    Vec3d *energy_grad_ = NULL;
    double *target_deformation_ = NULL;
    double *example_guided_deformation_ = NULL;
    double *temp_example_guided_deformation_=NULL;
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
    bool with_constrains_ = true;
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
    double *collide_vel_=NULL;
    double *u_initial_=NULL;
    double *vel_initial_=NULL;
    double *f_ext_=NULL;
    double *example_f_=NULL;
    // double *full_drag_force_=NULL;
    //used for reduced cubica element Computation
    double **restpos_ = NULL;//compute rest position for cubica elements
    double *deformed_ = NULL;
    double **LB_restpos_ = NULL;//compute rest position for LB_cubica elements
    double *q_=NULL;
    double *temp_q_=NULL;
    double *temp_grad_=NULL;
    double *qvel_=NULL;
    double *qaccel_=NULL;
    double *fq_=NULL;
    double *fqBase_=NULL;
    double *fq_ext_=NULL;
    double *reduced_mass_matrix_ = NULL;
    ModalMatrix *modal_matrix_ = NULL;
    double *U_ = NULL;
    VolumetricMesh *temp_mesh_ = NULL;

    //temp mesh_eigen_skeleton
    double *examples_deformation_=NULL;//temp
    // double *initial_eigen_skeleton_=NULL;
	// double *initial_detail_vector_=NULL;
    double *local_detail_vector_=NULL;
    double *target_reconstruction_=NULL;
    Vec3d *temp_eigencoefs_=NULL;
    double *examples_deformation0_=NULL;//temp
    Vec3d *temp_eigencoefs0_=NULL;
    std::map<int,int> vert_vertex1,vert_vertex2;

    //rigid body simulation
    RigidBody_GeneralTensor *rigid_=NULL;
    // RigidBody *rigid_=NULL;
    double total_mass_=0.0;
    double *vert_mass_=NULL;
    Vec3d rigid_center_;
    double initial_inertia_tensor_[9];
    // double inertia_tensor_[9];
    Vec3d linear_velocity_;
    Vec3d angular_velocity_;
    Vec3d new_linear_velocity_;
    Vec3d new_angular_velocity_;
    Vec3d rigid_external_force_;
    Vec3d initial_rigid_vel_;
    Vec3d am_;
    // Vec3d linear_acc_;
    // Vec3d angular_acc_;
    // double *f_cor_=NULL;
    // double *f_ine_=NULL;
    // double *f_enl_=NULL;
    // double *f_cen_=NULL;
    // double *vert_f_ext_=NULL;
    double R_[9];
    double new_R_[9];
    double t_[3];
    double *local_reference_=NULL;
    double *local_u_=NULL;
    double *global_u_=NULL;//not include the displacement generated by inertial force
    double *Uqvel_=NULL;
    double *Uqaccel_=NULL;
    int collide_vert_num_=1;
    double total_projection_time_=0.0;
    double total_target_time_=0.0;
    double total_reconstruction_time_=0.0;
    double total_time_=0.0;
    int timer_sample_interval_=100;
    // Matrix<double> elastic_tensor_(6,6,0.0);

};

} //namespace RTLB

#endif //REAL_TIME_EXAMPLE_BASED_DEFORMER_H_
