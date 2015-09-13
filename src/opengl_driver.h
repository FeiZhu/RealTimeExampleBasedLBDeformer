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
#include "configFile.h"
#include "planes.h"
#include "VegaHeaders.h"

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
class SparseMatrix;
class StVKElementABCD;
class IsotropicMaterial;
class IsotropicHyperelasticFEM;
class CorotationalLinearFEM;
class RenderVolumetricMesh;

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
    //GLUI callback methods
    static void updateRenderMesh(int code);
    static void updateCurrentExample(int code);
    static void changeSimulationMode(int code);
    static void loadObjectEigenfunctions(int code);
    static void saveObjectEigenfunctions(int code);
    static void loadExampleEigenfunctions(int code);
    static void saveExampleEigenfunctions(int code);
    static void loadReducedBasis(int code);
    static void loadObjectCubicaData(int code);
    static void exitApplication(int code);
    static void loadCorrespondenceData(int code);
    static void registerEigenfunctions(int code);
    //misc
    void drawAxis(double axis_length) const;
private:
    static const unsigned int string_length=1024;
    static OpenGLDriver *active_instance_;
    //simulation
    RealTimeExampleBasedDeformer *simulator_ = NULL;
    bool pause_simulation_ = true;
    int example_num_ = 0;
    int example_eigenfunction_num_ = 0;
    int object_eigenfunction_num_ = 0;
    char simulation_mesh_file_name_[string_length];
    char example_file_name_prefix_[string_length];
    char visual_mesh_file_name_[string_length];
    char reduced_basis_file_name_[string_length];
    char object_eigen_file_name_[string_length];
    char example_eigen_file_name_[string_length];
    char object_interpolation_file_name_[string_length];
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

    bool add_gravity_=false;
    double gravity_ = -9.8;
    double time_step_ = 1.0/30;
    float newmark_beta_=0.25;
    float newmark_gamma_=0.5;
    float damping_mass_coef_=0.0;
    float damping_stiffness_coef_=0.0;
    float damping_laplacian_coef_=0.0;
    double integrator_epsilon_=1.0e-6;
    double deformable_object_compliance_=1.0;
    double example_stiffness_scale_=1.0;//the stiffness used to compute example force is scaled
    int max_iterations_=1;
    int force_neighbor_size_=3;
    int solver_threads_num_=1;
    double principal_stretch_threshold_=-DBL_MAX;
    bool enable_eigen_weight_control_=false;    //enable the explicit weight control
    int corotational_linearFEM_warp_=1;
    int central_difference_tangential_damping_update_mode_=1;
    int positive_definite_=0;

    VolumetricMesh *simulation_mesh_=NULL;
    VolumetricMesh **example_mesh_=NULL;
    VolumetricMesh *current_simulation_mesh_=NULL;
    TetMesh *tet_mesh_=NULL;
    unsigned int simulation_vertices_num_ = 0;
    RenderVolumetricMesh *render_volumetric_mesh_=NULL;
    SceneObjectDeformable *rendering_mesh_=NULL;
    SceneObjectDeformable *visual_mesh_=NULL;
    Graph *mesh_graph_ = NULL;
    SparseMatrix *mass_matrix_=NULL;
    SparseMatrix *laplacian_matrix_=NULL;
    SparseMatrix *laplacian_damping_matrix_=NULL;
    StVKElementABCD *precomputed_integrals_=NULL,*example_precomputed_integrals_=NULL;
    StVKInternalForces *stvk_internal_force_=NULL,*example_stvk_internal_force_=NULL;
    StVKStiffnessMatrix *stvk_stiffness_matrix_=NULL,*example_stvk_stiffness_matrix_=NULL;
    CorotationalLinearFEM *corotational_linear_fem_=NULL,*example_corotational_linear_fem_=NULL;
    IsotropicMaterial *isotropic_material_=NULL,*example_isotropic_material_=NULL;
    IsotropicHyperelasticFEM *isotropic_hyperelastic_fem_=NULL,*example_isotropic_hyperelastic_fem_=NULL;
    ForceModel *force_model_;
    IntegratorBase *integrator_base_=NULL;
    IntegratorBaseSparse *integrator_base_sparse_=NULL;

    int pull_vertex_=-1; //the index of vertex pulled by user
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
    //window
    std::string window_name_ = std::string("Example-based Simulator");
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
    //render switches
    unsigned int enable_textures_=0;
    bool render_axis_ = true;
    bool render_vertices_ = false;
    bool render_wireframe_ = true;
    bool render_fixed_vertices_ = true;
    bool render_eigenfunction_ = false;
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
        UNKNOWN
    };
    SolverType solver_type_ = UNKNOWN;
    enum InvertibleMaterialType{
        INV_STVK,
        INV_NEOHOOKEAN,
        INV_MOONEYRIVLIN,
        INV_NONE
    };
    InvertibleMaterialType invertible_material_type_ = INV_NONE;
    enum DeformableObjectType{
        STVK,
        COROTLINFEM,
        LINFEM,
        MASSSPRING,
        INVERTIBLEFEM,
        UNSPECIFIED
    };
    DeformableObjectType deformable_object_type_ = UNSPECIFIED;
    //glui controls
    GLUI *glui_ = NULL;
    GLUI_StaticText *glui_object_surface_eigenfunctions_loaded;
    GLUI_StaticText *glui_example_with_eigenfunctions_loaded;
    GLUI_StaticText *glui_current_example_eigenfunctions_loaded;
    GLUI_StaticText *glui_rendering_eigenfunctions_enabled;
    GLUI_Spinner *render_eigen_index_spinner;
    GLUI_Button *change_simulation_mode_button;
    unsigned int current_example_index_ = 1;
    unsigned int example_with_eigen_num_ = 0;
    unsigned int save_eigenfunction_num_ = 0;
    ConfigFile config_file_;
};

}  //namespace RTLB

#endif //OPENGL_DRIVER_H_
