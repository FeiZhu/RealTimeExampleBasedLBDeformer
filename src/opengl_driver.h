/*
 * @file: opengl_driver.h
 * @brief: UI interface for simulation
 * @author: Fei Zhu
 *
 */

#ifndef OPENGL_DRIVER_H_
#define OPENGL_DRIVER_H_

#include <string>
#include <cfloat>
#include "planes.h"
#include "VegaHeaders.h"
#include "reducedNeoHookeanForceModel.h"

class GLUI;
class GLUI_StaticText;
class GLUI_Spinner;
class GLUI_Button;
class SphericalCamera;
class Lighting;
class SceneObjectDeformable;
class ConfigFile;
class planes;
class Graph;
class ReducedNeoHookeanForceModel;
class ModalMatrix;

namespace RTLB{

class RealTimeExampleBasedDeformer;


class OpenGLDriver
{
public:
    explicit OpenGLDriver(const std::string &config_file_name);
    ~OpenGLDriver();
private:
    static OpenGLDriver* activeInstance();
    //init functions
    void initGLUT();
    void initGLUI();
    void initConfigurations(const std::string &config_file_name);
    void initSimulation();
    void initCamera();
    void initGraphics();
    //callback functions
    static void displayFunction(void);
    static void idleFunction(void);
    static void reshapeFunction(int width, int height);
    static void keyboardFunction(unsigned char key, int x, int y);
    static void specialFunction(int key, int x, int y);
    static void motionFunction(int x, int y);
    static void mouseFunction(int button, int state, int x, int y);
    static void addGravitySwitch(bool add_gravity);
    //GLUI callback methods
    static void updateRenderMesh(int code);
    static void updateCurrentExample(int code);
    static void changeCurrentEigenIndex(int code);
    static void changeSimulationMode(int code);
    static void enableReducedSimulation(int code);
    static void loadObjectEigenfunctions(int code);
    static void saveObjectEigenfunctions(int code);
    static void loadExampleEigenfunctions(int code);
    static void saveExampleEigenfunctions(int code);
    static void loadReducedBasis(int code);
    static void loadObjectCubicaData(int code);
    static void loadExampleCubicaData(int code);
    static void exitApplication(int code);
    static void loadCorrespondenceData(int code);
    static void registerEigenfunctions(int code);
    static void resetDeformation(int code);

    static void saveTetMesh(int code);
    //misc
    void drawAxis(double axis_length) const;
    void drawIndexColorTable() const;
private:
    static const unsigned int string_length=1024;
    static OpenGLDriver *active_instance_;
    //simulation
    RealTimeExampleBasedDeformer *simulator_ = NULL;
    bool pause_simulation_ = true;
    int example_num_ = 0;
    int interpolate_eigenfunction_num_ = 0;
    int reconstruct_eigenfunction_num_ = 0;
    //int coupled_eigenfunction_num_ = 0;
    char simulation_mesh_file_name_[string_length];
    char example_file_name_prefix_[string_length];
    char visual_mesh_file_name_[string_length];
    char volumetric_surface_mesh_file_name_[string_length];
    char reduced_basis_file_name_[string_length];
    char object_eigen_file_name_[string_length];
    char example_eigen_file_name_prefix_[string_length];
    char output_object_eigen_file_name_[string_length];
    char output_eigen_file_name_prefix_[string_length];
    char object_interpolation_file_name_[string_length];
    char example_cubica_file_name_prefix_[string_length];
    char object_cubica_file_name_[string_length];
    char fixed_vertex_file_name_[string_length];
    char force_loads_file_name_[string_length];
    char initial_position_file_name_[string_length];
    char initial_velocity_file_name_[string_length];
    Planes *planes_=NULL;
    char plane_file_name_[string_length];
    int plane_num_ = 0;
    char mass_matrix_file_name_[string_length];
    char corresponding_file_name_[string_length];
    char deformable_model_[string_length];
    char invertible_material_[string_length];
    //char objectEigen
    char object_affected_vertices_file_name_[string_length];
    char example_affected_vertices_file_name_[string_length];
    char extra_objects_file_name_base_[string_length];

    bool add_gravity_=false;
    double gravity_ = 9.8;
    double time_step_ = 1.0/30;
    int time_step_counter_=0;
    float newmark_beta_=0.25;
    float newmark_gamma_=0.5;
    float damping_mass_coef_=0.0;
    float damping_stiffness_coef_=0.0;
    float damping_laplacian_coef_=0.0;
    double epsilon_=1.0e-12;//numerical accuracy in this file
    double integrator_epsilon_=1.0E-6;//numerical accuracy used in integrator
    double deformable_object_compliance_=1.0;
    double example_stiffness_scale_=1.0;//the stiffness used to compute example force is scaled
    int max_iterations_=1;
    int force_neighbor_size_=3;//the influence range of the integration force
    int solver_threads_num_=1;//number of threads used for integration solver
    double principal_stretch_threshold_=-DBL_MAX;
    bool enable_eigen_weight_control_=false;    //enable the explicit weight control
    int corotational_linearFEM_warp_=1;
    int central_difference_tangential_damping_update_mode_=1;
    int positive_definite_=0;
    double last_initial_weight_=1.5;//useful when explicit weight control enabled

    VolumetricMesh *simulation_mesh_=NULL;
    VolumetricMesh **example_mesh_=NULL;
    VolumetricMesh *current_example_mesh_=NULL;
    TetMesh *tet_mesh_=NULL;
    unsigned int simulation_vertices_num_ = 0;
    RenderVolumetricMesh *render_volumetric_mesh_=NULL;
    SceneObjectDeformable *render_surface_mesh_=NULL;
    SceneObjectDeformable *visual_mesh_=NULL;
    double *u_render_surface_=NULL;
    Vec3d **example_eigencoefs_ = NULL; //the coefficients of the example geometry projected onto the eigenfunctions
    Vec3d *initial_object_eigencoefs_ = NULL;//the coefficients of the object projected onto all eigenfunctions
    Vec3d *deformed_object_eigencoefs_ = NULL;//the coefficients of the object projected onto the interpolation eigenfunctions
    Vec3d *target_eigencoefs_ = NULL;
    double *initial_object_configurations_ = NULL;//initial simulation mesh position,size:3n*1
    double *deformed_object_configurations_ = NULL;//deformed simulation mesh position for each step ,size:3n*1
    double *temp_deformed_object_dis_ = NULL;//displacement(current_configuration-target_configuration)
    //deformation
    //double *object_elastic_subspace_force_ = NULL;
    //double *object_elastic_fullspace_force_ = NULL;
    double *example_guided_subspace_force_ = NULL;
    double *example_guided_fullspace_force_ = NULL;
    double *example_guided_deformation_ = NULL;
    double *target_initial_deformation_ = NULL;

    Graph *mesh_graph_ = NULL;
    SparseMatrix *mass_matrix_= NULL;
    SparseMatrix *laplacian_matrix_ = NULL;
    SparseMatrix *laplacian_damping_matrix_ = NULL;
    StVKElementABCD *precomputed_integrals_ = NULL;
//    StVKInternalForces *stvk_internal_force_ = NULL,*example_stvk_internal_force_ = NULL;
//    StVKStiffnessMatrix *stvk_stiffness_matrix_ = NULL,*example_stvk_stiffness_matrix_=NULL;
//    CorotationalLinearFEM *corotational_linear_fem_ = NULL,*example_corotational_linear_fem_ = NULL;
    IsotropicMaterial *isotropic_material_ = NULL,*example_isotropic_material_ = NULL;
    IsotropicHyperelasticFEM *isotropic_hyperelastic_fem_ = NULL,*example_isotropic_hyperelastic_fem_ = NULL;
    ForceModel *force_model_ = NULL;
    IntegratorBase *integrator_base_ = NULL;
    ImplicitNewmarkSparse *implicit_newmark_sparse_ = NULL;
    IntegratorBaseSparse *integrator_base_sparse_ = NULL;
    IntegratorBaseDense *integrator_base_dense_ = NULL;
    ImplicitNewmarkDense *implicit_newmark_dense_ = NULL;
    ImplicitBackwardEulerDense *implicit_backward_euler_dense_ = NULL;
    CentralDifferencesDense *central_differences_dense_ = NULL;


    int pulled_vertex_=-1; //the index of vertex pulled by user
    double *u_=NULL;
    double *vel_=NULL;
    double *u_initial_=NULL;
    double *vel_initial_=NULL;
    double *f_ext_=NULL;
    double *f_col_=NULL;
    int fixed_vertices_num_=0;
    int *fixed_vertices_=NULL;
    int fixed_dofs_num_=0;
    int *fixed_dofs_=NULL;
    int force_loads_num_=0;
    double *force_loads_=NULL;
    int object_interpolation_element_vertices_num_=0;
    int *object_interpolation_vertices_=NULL;
    double *object_interpolation_weights_=NULL;
    int extra_objects_num_=0;
    SceneObjectDeformable **extra_objects_=NULL;
    //window
    std::string window_name_ = std::string("Real-time example-based Simulator");
    int window_id_ = -1;
    int window_width_ = 800;
    int window_height_ = 600;
    //camera
    double znear_ = 0.01;
    double zfar_ = 10.0;
    double camera_radius_ = 17.5;
    double focus_position_[3];
    double camera_longitude_ = -60.0;
    double camera_lattitude_ = 20.0;
    SphericalCamera *camera_ = NULL;
    //light
    Lighting *lighting_ = NULL;
    char lighting_config_file_name_[string_length];
    //mouse
    int mouse_pos_[2];
    bool left_button_down_ = false;
    bool middle_button_down_ = false;
    bool right_button_down_ = false;
    int shift_pressed_=0;
    int alt_pressed_=0;
    int ctrl_pressed_=0;
    int drag_start_x_,drag_start_y_;
    double render_velocity_scale_ = 1.0;
    //render switches
    unsigned int enable_textures_=0;
    bool render_axis_ = true;
    bool render_vertices_ = false;
    bool render_wireframe_ = true;
    bool render_fixed_vertices_ = true;
    bool render_velocity_ = false;
    bool render_dis_ = false;
    bool render_eigenfunction_ = false;
    bool isrender_surface_mesh_ = true;
    bool isrender_volumetric_mesh_ = false;
    bool isload_example_eigen_ = false;
    bool isload_object_eigen_ = false;
    bool isload_correspondence_data_ = false;
    bool isload_cubica_ = false;
    bool isload_example_cubica_ = false;
    bool enable_example_simulation_ = false;
    bool save_tet_mesh_ = false;
    char solver_method_[string_length];
    enum RenderMeshType{
        VISUAL_MESH = 0,
        OBJECT_EIGEN_MESH,
        EXAMPLE_MESH
    };
    RenderMeshType render_mesh_type_ = VISUAL_MESH;
    enum SolverType{
        IMPLICITNEWMARK,
        IMPLICITBACKWARDEULER,
        EULER,
        SYMPLECTICEULER,
        CENTRALDIFFERENCES,
        REDUCEDCENTRALDIFFERENCES,
        REDUCEDIMPLICITNEWMARK,
        REDUCEDIMPLICITBACKWARDEULER,
        UNKNOWN
    };
    SolverType solver_type_ = UNKNOWN;
    enum InvertibleMaterialType{
        INV_STVK,
        INV_NEOHOOKEAN,
        REDUCED_STVK,
        REDUCED_INV_NEOHOOKEAN,
        INV_NONE
    };
    InvertibleMaterialType invertible_material_type_ = INV_NONE;
    enum DeformableObjectType{
        INVERTIBLEFEM,
        REDUCEDFEM,
        UNSPECIFIED
    };
    DeformableObjectType deformable_object_type_ = UNSPECIFIED;
    enum SimulationMode{
        FULLSPACE,
        REDUCEDSPACE
    };
    SimulationMode simulation_mode_ = FULLSPACE;
    //glui controls
    GLUI *glui_ = NULL;
    GLUI_StaticText *glui_object_surface_eigenfunctions_loaded_;
    GLUI_StaticText *glui_example_with_eigenfunctions_loaded_;
    GLUI_StaticText *glui_current_example_eigenfunctions_loaded_;
    GLUI_StaticText *glui_rendering_eigenfunctions_enabled_;
    GLUI_Spinner *render_eigen_index_spinner_;
    GLUI_Button *change_simulation_mode_button_;
    GLUI_Button *reduced_simulation_button_;
    unsigned int current_example_index_ = 0;
    unsigned int current_render_eigen_idx_=0;
    unsigned int save_eigenfunction_num_ = 0;
    ConfigFile config_file_;
    //reuced simulation
    bool enable_reduced_simulation_ = false;
    SceneObjectReducedCPU *render_reduced_surface_mesh_=NULL;
    ReducedForceModel *reduced_force_model_=NULL;
    ReducedNeoHookeanForceModel *reduced_neoHookean_force_model_=NULL;
    int r_ = 0;
    double *reduced_mass_matrix_ = NULL;
    double *reduced_stiffness_matrix_ = NULL;
    double *reduced_f_ext_ = NULL;
    Vec3d *vec_q_ = NULL;
    double *q_ = NULL;
    bool reduced_simulation_ = false;
    ModalMatrix *modal_matrix_ = NULL;
    double *fq_ = NULL;
    double *fqBase_ = NULL;
    double *U_ = NULL;
};

}  //namespace RTLB

#endif //OPENGL_DRIVER_H_
