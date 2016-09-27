/*
 * @file: opengl_driver.cpp
 * @brief: UI interface for simulation
 * @author: Fei Zhu
 *
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <set>
#include "GL/freeglut.h"
#include "GL/glui.h"
#include "camera.h"
#include "lighting.h"
#include "real_time_example_based_deformer.h"
#include "opengl_driver.h"
#include "colorTable.h"
using namespace std;
using std::set;

namespace RTLB{

OpenGLDriver* OpenGLDriver::active_instance_ = NULL;

OpenGLDriver::OpenGLDriver(const std::string &config_file_name)
{
    if(active_instance_)
        delete active_instance_;
    active_instance_ = this;
    simulator_=new RTLB::RealTimeExampleBasedDeformer();
    //TO DO: init everything and enter mainloop
    initConfigurations(config_file_name);
    initGLUT();
    initGLUI();
    initGraphics();
    initSimulation();
    glutMainLoop();
}

OpenGLDriver::~OpenGLDriver()
{
    exitApplication(0);
}

void OpenGLDriver::initConfigurations(const std::string &config_file_name)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    std::cout<<"Parsing configuration file "<<config_file_name<<std::endl;
    config_file_.addOptionOptional("windowWidth",&window_width_,window_width_);
    config_file_.addOptionOptional("windowHeight",&window_height_,window_height_);
    config_file_.addOptionOptional("cameraRadius",&camera_radius_,camera_radius_);
    config_file_.addOptionOptional("cameraLongitude",&camera_longitude_,camera_longitude_);
    config_file_.addOptionOptional("cameraLattitude",&camera_lattitude_,camera_lattitude_);
    config_file_.addOptionOptional("focusPositionX",&focus_position_[0],focus_position_[0]);
    config_file_.addOptionOptional("focusPositionY",&focus_position_[1],focus_position_[1]);
    config_file_.addOptionOptional("focusPositionZ",&focus_position_[2],focus_position_[2]);
    config_file_.addOption("lightingConfigFilename",lighting_config_file_name_);

    //planes config
    config_file_.addOptionOptional("planeNumber",&plane_num_,plane_num_);
    config_file_.addOptionOptional("planesFilename",plane_file_name_,plane_file_name_);
    config_file_.addOptionOptional("extraObjectsFilenameBase",extra_objects_file_name_base_,"extraObject_");
    //cubica FILE
    config_file_.addOptionOptional("isLoadCubica",&isload_cubica_,isload_cubica_);
    config_file_.addOptionOptional("isLoadLBCubica",&isload_LB_cubica_,isload_LB_cubica_);
    config_file_.addOptionOptional("withConstrains",&with_constrains_,true);
    config_file_.addOptionOptional("objectCubicaFilename",object_cubica_file_name_,"none");
    config_file_.addOptionOptional("objectLBCubicaFilename",object_LB_cubica_file_name_,"none");
    config_file_.addOptionOptional("exampleCubicaFilenameBase",example_cubica_file_name_prefix_,"none");
    //solver and materials config
    config_file_.addOptionOptional("simulationType",simulation_type_,"UNKNOWN");
    config_file_.addOptionOptional("invertibleMaterial",invertible_material_type_,"reducedStVK");
    config_file_.addOptionOptional("solver",solver_method_,"UNKNOWN");
    config_file_.addOptionOptional("deformableModel",deformable_model_,"neoHookean");
    config_file_.addOptionOptional("invertibleMaterial",invertible_material_,"none");
    config_file_.addOptionOptional("principalStretchThreshold",&principal_stretch_threshold_,principal_stretch_threshold_);

    //simulation object
    config_file_.addOptionOptional("simulationMeshFilename",simulation_mesh_file_name_,"none");
    config_file_.addOptionOptional("objectRenderSurfaceMeshFilename",visual_mesh_file_name_,"none");
    config_file_.addOptionOptional("objectVolumetricSurfaceMeshFilename",volumetric_surface_mesh_file_name_,"none");
    config_file_.addOptionOptional("objectInterpolationFilename",object_interpolation_file_name_,"none");
    config_file_.addOptionOptional("objectInitialTetMesh",initial_tetmesh_file_name_,"none");
    config_file_.addOptionOptional("initialRigidVelX",&initial_rigidvel_x_,initial_rigidvel_x_);
    config_file_.addOptionOptional("initialRigidVelY",&initial_rigidvel_y_,initial_rigidvel_y_);
    config_file_.addOptionOptional("initialRigidVelZ",&initial_rigidvel_z_,initial_rigidvel_z_);
    //examples
    config_file_.addOptionOptional("exampleFilenameBase",example_file_name_prefix_,"none");
    //config_file_.addOptionOptional("exampleSurfaceMeshFilenameBase",example_surface_mesh_file_name_,"none");
    //config_file_.addOptionOptional("exampleInterpolationFilenameBase",example_interpolation_file_name_,"none");

    //eigen configuration
    config_file_.addOption("exampleNum",&example_num_);
    config_file_.addOption("interpolateEigenfunctionNum",&interpolate_eigenfunction_num_);
    config_file_.addOption("reconstructEigenfunctionNum",&reconstruct_eigenfunction_num_);
    config_file_.addOptionOptional("exampleEigenFunctionFilenameBase",example_eigen_file_name_prefix_,"none");
    config_file_.addOptionOptional("objectEigenFunctionFilename",object_eigen_file_name_,"none");
    config_file_.addOptionOptional("savedExampleEigenFunctionFilenameBase",output_eigen_file_name_prefix_,"none");
    config_file_.addOptionOptional("savedObjectEigenFunctionFilename",output_object_eigen_file_name_,"none");
    config_file_.addOptionOptional("saveObjMeshFilenameBase",output_objmesh_file_name_base_,"none");
    config_file_.addOptionOptional("enableSaveSurfaceMtl",&enable_save_surfacemtl_,enable_save_surfacemtl_);
    config_file_.addOptionOptional("correspondingFunctionFilename",corresponding_file_name_,"none");
    config_file_.addOptionOptional("reducedBasisFilename",reduced_basis_file_name_,"none");
    //output
    //config_file_.addOptionOptional("outputFilenameBase",output_file_name_,"output");

    //Simulation
    config_file_.addOptionOptional("forceNeighborSize",&force_neighbor_size_,force_neighbor_size_);
    config_file_.addOptionOptional("enableEigenWeightControl",&enable_eigen_weight_control_,enable_eigen_weight_control_);
    config_file_.addOptionOptional("timestep",&time_step_,time_step_);
    config_file_.addOptionOptional("framerate",&frame_rate_,frame_rate_);
    config_file_.addOptionOptional("totalFrame",&total_frames_,total_frames_);
    config_file_.addOptionOptional("massMatrixFilename",mass_matrix_file_name_,"none");
    config_file_.addOptionOptional("inertiaTensorFilename",inertia_tensor_file_name_,"none");
    config_file_.addOptionOptional("dampingMassCoef",&damping_mass_coef_,damping_mass_coef_);
    config_file_.addOptionOptional("dampingStiffnessCoef",&damping_stiffness_coef_,damping_stiffness_coef_);
    config_file_.addOptionOptional("dampingLaplacianCoef",&damping_laplacian_coef_,damping_laplacian_coef_);
    config_file_.addOptionOptional("integratorEpsilon",&integrator_epsilon_,integrator_epsilon_);
    config_file_.addOptionOptional("deformableObjectCompliance",&deformable_object_compliance_,deformable_object_compliance_);
    config_file_.addOptionOptional("exampleBias",&example_bias_,example_bias_);
    config_file_.addOptionOptional("exampleStiffnessScale",&example_stiffness_scale_,example_stiffness_scale_);
    config_file_.addOptionOptional("dampingExampleStiffness",&damping_example_stiffness_,damping_example_stiffness_);
    config_file_.addOptionOptional("maxIterations",&max_iterations_,max_iterations_);
    config_file_.addOptionOptional("torqueCoef",&torque_coef_,torque_coef_);
    config_file_.addOptionOptional("collisionNumLimited",&col_limited_num_,col_limited_num_);

    //initial conditions
    config_file_.addOptionOptional("fixedVerticesFilename",fixed_vertex_file_name_,"none");
    config_file_.addOptionOptional("forceLoadsFilename",force_loads_file_name_,"none");
    config_file_.addOptionOptional("initialPositionFilename",initial_position_file_name_,"none");
    config_file_.addOptionOptional("initialVelocityFilename",initial_velocity_file_name_,"none");
    //gravity
    config_file_.addOptionOptional("addGravity",&add_gravity_,add_gravity_);
    config_file_.addOptionOptional("gravity",&gravity_,gravity_);
    //local Example-based
    config_file_.addOptionOptional("objectAffectedVerticesFilename",object_affected_vertices_file_name_,"none");
    config_file_.addOptionOptional("exampleAffectedVerticesFilenameBase",example_affected_vertices_file_base_,"none");

    //parse the configuration file
    if(config_file_.parseOptions((char*)config_file_name.c_str())!=0)
    {
        std::cout<<"Error: parsing options failed.\n";
        exit(0);
    }
    //print the variables that were just parsed
    config_file_.printOptions();
    simulator_->setTimeStep(time_step_);
    simulator_->setFrameRate(frame_rate_);
    simulator_->setTotalFrames(total_frames_);
    simulator_->setDampingMassCoef(damping_mass_coef_);
    simulator_->setDampingStiffnessCoef(damping_stiffness_coef_);
    simulator_->setExampleStiffnessScale(example_stiffness_scale_);
    simulator_->setDampingExampleStiffness(damping_example_stiffness_);
    simulator_->setExampleBiasFactor(example_bias_);
    std::cout<<"example_bias_:"<<example_bias_<<".........................\n";
    simulator_->setTorqueCoef(torque_coef_);
    simulator_->setCollisionNumLimited(col_limited_num_);
    simulator_->setPrincipalStretchThreshold(principal_stretch_threshold_);
    //enable eigen weight control
    simulator_->setEnableEigenWeightControl(enable_eigen_weight_control_);
    simulator_->enableGravity(add_gravity_);
    simulator_->setGravity(gravity_);
    simulator_->setSolverType(solver_method_);
    simulator_->setConstrains(with_constrains_);
    //set simulation mode
    if(strcmp(simulation_type_,"fullspace")==0)
    {
        simulation_mode_=FULLSPACE;
    }
    if(strcmp(simulation_type_,"reducedspace")==0)
    {
        simulation_mode_=REDUCEDSPACE;
    }
    if(strcmp(simulation_type_,"UNKNOWN")==0)
    {
        std::cout<<"Error:unknown simulation mode specified."<<std::endl;
        exit(0);
    }
    // if(strcmp(force_loads_file_name_,"none")!=0)
    //     simulator_->setInitialForceFilename(force_loads_file_name_);
    // if(strcmp(initial_velocity_file_name_,"none")!=0)
    //     simulator_->setInitialVelFilename(initial_velocity_file_name_);
    // if(strcmp(initial_position_file_name_,"none")!=0)
    //     simulator_->setInitialPosFilename(initial_position_file_name_);
    simulator_->setSimulationType(simulation_type_);
    simulator_->setMaterialType(invertible_material_type_);
    simulator_->setRigidInitialVel(initial_rigidvel_x_,initial_rigidvel_y_,initial_rigidvel_z_);
    std::cout<<plane_num_<<"\n";
}

OpenGLDriver* OpenGLDriver::activeInstance()
{
    return active_instance_;
}

void OpenGLDriver::initGLUT()
{
    int argc = 0;
    char *argv = "dummy argv\n";
    glutInit(&argc, &argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH|GLUT_STENCIL);
    glutInitWindowSize(window_width_,window_height_);
    window_id_ = glutCreateWindow(window_name_.c_str());
    glutDisplayFunc(displayFunction);
    GLUI_Master.set_glutIdleFunc(idleFunction);
    GLUI_Master.set_glutKeyboardFunc(keyboardFunction);
    GLUI_Master.set_glutReshapeFunc(reshapeFunction);
    GLUI_Master.set_glutMouseFunc(mouseFunction);
    glutMotionFunc(motionFunction);
}

void OpenGLDriver::initGLUI()
{
    glui_ = GLUI_Master.create_glui("Controls",0,window_width_+52,0);
    //init GLUI control panel
    GLUI_Panel *display_mode_panel = glui_->add_panel("Display Mode",GLUI_PANEL_EMBOSSED);
    display_mode_panel->set_alignment(GLUI_ALIGN_LEFT);
    //switch render mesh
    GLUI_RadioGroup *display_mode_radio_group = glui_->add_radiogroup_to_panel(display_mode_panel,
                                                        (int*)&render_mesh_type_,0,updateRenderMesh);
    glui_->add_radiobutton_to_group(display_mode_radio_group,"Visual Mesh");
    glui_->add_radiobutton_to_group(display_mode_radio_group,"Object Eigen Mesh");
    //example_num_=simulator_->exampleNum();
    if(example_num_ > 0)
        glui_->add_radiobutton_to_group(display_mode_radio_group,"Example Mesh");
    //example pose information
    GLUI_Panel *examples_panel = glui_->add_panel("Examples",GLUI_PANEL_EMBOSSED);
    examples_panel->set_alignment(GLUI_ALIGN_LEFT);
    std::stringstream adaptor;
    adaptor<<example_num_;
    std::string example_pose_num_str, static_text_content("Number of examples: ");
    adaptor>>example_pose_num_str;
    static_text_content += example_pose_num_str;
    glui_->add_statictext_to_panel(examples_panel,static_text_content.c_str());
    if(example_num_ > 0)
    {
        GLUI_Spinner *example_index_spinner = glui_->add_spinner_to_panel(examples_panel,"Current example index: ",
                                                  GLUI_SPINNER_INT,&current_example_index_,0,updateCurrentExample);
        example_index_spinner->set_int_limits(1,example_num_,GLUI_LIMIT_CLAMP);
    }
    //eigenfunction information
    GLUI_Panel *eigen_panel = glui_->add_panel("Laplace-Beltrami Eigenfunctions",GLUI_PANEL_EMBOSSED);
    eigen_panel->set_alignment(GLUI_ALIGN_LEFT);
    std::string eigenfunction_num_str;
    std::string reconstruct_eigenfunction_num_str;
    if(example_num_ > 0)
    {
        adaptor.clear();
        adaptor<<interpolate_eigenfunction_num_;
        adaptor>>eigenfunction_num_str;
        static_text_content.clear();
        static_text_content = std::string("Number of eigenfunctions for examples: ");
        static_text_content += eigenfunction_num_str;
        adaptor.clear();
        adaptor<<reconstruct_eigenfunction_num_;
        adaptor>>reconstruct_eigenfunction_num_str;
        static_text_content += "/" + reconstruct_eigenfunction_num_str;
        glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
        adaptor.clear();
        glui_current_example_eigenfunctions_loaded_=glui_->add_statictext_to_panel(eigen_panel,"Eigenfunctions for current example loaded:No");
        glui_->add_button_to_panel(eigen_panel,"Load Eigenfunctions for All Examples",0,loadExampleEigenfunctions);
        glui_->add_separator_to_panel(eigen_panel);
    }
    adaptor.clear();

    std::string interpolate_eigenfunction_num_str;
    adaptor<<interpolate_eigenfunction_num_;
    adaptor>>interpolate_eigenfunction_num_str;
    adaptor.clear();
    adaptor<<reconstruct_eigenfunction_num_;
    adaptor>>reconstruct_eigenfunction_num_str;
    static_text_content.clear();
    static_text_content="Number of eigenfunctions for simulation object:";
    static_text_content+=interpolate_eigenfunction_num_str;
    static_text_content+="/"+reconstruct_eigenfunction_num_str;
    glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    static_text_content.clear();
    static_text_content="Simulation object surface eigenfunctions loaded:No";
    glui_object_surface_eigenfunctions_loaded_=glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    glui_->add_button_to_panel(eigen_panel,"Load Eigenfunctions for Simulation Object Surface",1,loadObjectEigenfunctions);
    glui_->add_separator_to_panel(eigen_panel);
    //glui_->add_button_to_panel(eigen_panel,"Load reduced Basis",0,loadReducedBasis);
    static_text_content.clear();
    static_text_content="Rendering eigenfunctions enable:No";
    glui_rendering_eigenfunctions_enabled_=glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    render_eigen_index_spinner_=glui_->add_spinner_to_panel(eigen_panel,"Current eigenfunction index: ",GLUI_SPINNER_INT,&current_render_eigen_idx_,0,changeCurrentEigenIndex);
    render_eigen_index_spinner_->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);

    //coupled quasi-harmonic bases
    if(example_num_>0)
    {
        GLUI_Panel *coupled_panel=glui_->add_panel("Coupled quasi-harmonics",GLUI_PANEL_EMBOSSED);
        coupled_panel->set_alignment(GLUI_ALIGN_LEFT);
        //GLUI_Spinner *coupled_base_num_spinner=glui_->add_spinner_to_panel(coupled_panel,"Coupled Eigenfunction Number:",GLUI_SPINNER_INT,&coupled_eigenfunction_num_);
        //coupled_base_num_spinner->set_int_limits(1,(interpolate_eigenfunction_num_<interpolate_eigenfunction_num_?interpolate_eigenfunction_num_:interpolate_eigenfunction_num_));
        glui_->add_button_to_panel(coupled_panel,"Load All Corresponding Functions",0,loadCorrespondenceData);
        glui_->add_button_to_panel(coupled_panel,"register Eigenfunctions",0,registerEigenfunctions);
    }
    //save eigenfunctions
    GLUI_Panel *save_panel=glui_->add_panel("Save Eigenfunctions",GLUI_PANEL_EMBOSSED);

    if(example_num_>0)
    {
        GLUI_Spinner *saved_base_num_spinner=glui_->add_spinner_to_panel(save_panel,"Saved Eigenfunction Number: ",GLUI_SPINNER_INT,&save_eigenfunction_num_);
        saved_base_num_spinner->set_int_limits(1,(interpolate_eigenfunction_num_<interpolate_eigenfunction_num_?interpolate_eigenfunction_num_:interpolate_eigenfunction_num_));
        //glui_->add_button_to_panel(save_panel,"Save Eigenfunctions for Current Example",0,/*saveCurrentEigenfunctions*/saveExampleEigenfunctions);
        glui_->add_button_to_panel(save_panel,"Save Eigenfunctions for All Examples",0,saveExampleEigenfunctions);
        glui_->add_button_to_panel(save_panel,"Save color for all example eigenfunctions",0,saveExampleEigenfunctionsColor);
    }
    glui_->add_button_to_panel(save_panel,"Save Eigenfunctions for Simulation Object Surface",1,saveObjectEigenfunctions);
    glui_->add_button_to_panel(save_panel,"Save object color for all eigenfunctions",0,saveObjectEigenfunctionsColor);

    //simulation
    GLUI_Panel *sim_panel=glui_->add_panel("Simulation",GLUI_PANEL_EMBOSSED);
    //add_rollout("Simulation",true);
    sim_panel->set_alignment(GLUI_ALIGN_LEFT);
    glui_->add_button_to_panel(sim_panel,"Reset",0,/*resetDeformation*/changeSimulationMode);
    //glui_->add_button_to_panel(sim_panel,"Save Current object Surface",0,/*saveCurrentObjectSurface*/changeSimulationMode);
    if(example_num_>0)
        change_simulation_mode_button_=glui_->add_button_to_panel(sim_panel,"Enable Example-based Simulation",0,changeSimulationMode);
//
    //reduced_simulation_button_=glui_->add_button_to_panel(sim_panel,"Enable Reduced Simulation",0,enableReducedSimulation);

    //exit button
    glui_->add_button("Exit",0,exitApplication);
    glui_->sync_live();
    glui_->set_main_gfx_window(window_id_);
}

void OpenGLDriver::initSimulation()
{
    std::cout<<"initSimulation:\n";
    //
    // getchar();
    total_steps_=(int)((1.0/time_step_)/frame_rate_)*total_frames_;
    try
    {
        lighting_ = new Lighting(lighting_config_file_name_);
    }
    catch(int exception_code)
    {
        std::cout<<"Error: "<<exception_code<<" reading lighting information from "<<lighting_config_file_name_<<".\n";
        exit(0);
    }
    //init planes in scene
    if(plane_num_>0)
    {
        if(strcmp(plane_file_name_,"none")!=0)
        {
            planes_=new Planes(plane_file_name_,plane_num_);
            // simulator_->loadPlanesInScene(plane_file_name_,plane_num_);
        }
        else
        {
            std::cout<<"Error: configuration file for planes in scene not specified.\n";
            exit(1);
        }
    }
    //load simulation mesh
    std::cout<<"Loading object volumetric mesh from file "<<simulation_mesh_file_name_<<".\n";
    if(simulator_->loadSimulationMesh(simulation_mesh_file_name_))
    {
        simulation_mesh_=simulator_->simulationMesh();
    }
    else
    {
        std::cout<<"Error: unable to load the object simulation mesh from "<<simulation_mesh_file_name_<<".\n";
        exit(1);
    }
    simulation_vertices_num_=simulation_mesh_->getNumVertices();
    mesh_graph_=GenerateMeshGraph::Generate(simulation_mesh_);

    //init surface mesh of the volumetric Mesh
    if(strcmp(visual_mesh_file_name_,"none")==0)
    {
        std::cout<"Error: the surface mesh of the volumetric mesh was not specified.\n";
        exit(1);
    }
    visual_mesh_=new SceneObjectDeformable(visual_mesh_file_name_);

    //load cubica file
    if(isload_cubica_)
    {
        if(strcmp(object_cubica_file_name_,"none")==0)
        {
            std::cout<<"Error: object cubica file unloaded.\n";
            exit(0);
        }
        std::cout<<"load cubica:\n";
        loadObjectCubicaData(0);
        std::cout<<"load cubica done\n";
    }
    loadReducedBasis(0);
    //init object render surface Mesh
    isrender_volumetric_mesh_=false;
    isrender_surface_mesh_=true;
    visual_mesh_=visual_mesh_;

    vel_initial_=new double[3*simulation_vertices_num_];
    memset(vel_initial_,0.0,sizeof(double)*3*simulation_vertices_num_);
    u_=new double[3*simulation_vertices_num_];
    memset(u_,0.0,sizeof(double)*3*simulation_vertices_num_);
    collide_u_=new double[3*simulation_vertices_num_];
    memset(collide_u_,0.0,sizeof(double)*3*simulation_vertices_num_);
    collide_vel_=new double[3*simulation_vertices_num_];
    memset(collide_vel_,0.0,sizeof(double)*3*simulation_vertices_num_);
    f_ext_=new double[3*simulation_vertices_num_];
    memset(f_ext_,0.0,sizeof(double)*3*simulation_vertices_num_);
    f_col_=new double[3*simulation_vertices_num_];
    memset(f_col_,0.0,sizeof(double)*3*simulation_vertices_num_);
    r_=simulator_->reducedBasisNum();
    q_= new double[r_];
    memset(q_,0.0,sizeof(double)*r_);
    fq_= new double[r_];
    memset(fq_,0.0,sizeof(double)*r_);
    fqBase_= new double[r_];
    memset(fqBase_,0.0,sizeof(double)*r_);
    fq_plane_=new double[r_];
    memset(fq_plane_,0.0,sizeof(double)*r_);
    collide_vert_num_=new int[1];
    //load mass matrix
    if(strcmp(mass_matrix_file_name_,"none")==0)
    {
        std::cout<<"Error: mass matrix for the deformable model not specified.\n";
        exit(1);
    }
    std::cout<<"Loading the mass matrix from file "<<mass_matrix_file_name_<<".\n";
    loadMassmatrix(0);
    if((simulation_mode_==REDUCEDSPACE)&&(!with_constrains_))
        loadInertiaTensor(0);
    // if((simulation_mode_==REDUCEDSPACE)&&(with_constrains_))
    // {
    //     // if(with_constrains_)
    //     // {
    //         render_reduced_surface_mesh_cpu_ = new SceneObjectReducedCPU(volumetric_surface_mesh_file_name_,simulator_->getModalmatrix());
    //         render_reduced_surface_mesh_=render_reduced_surface_mesh_cpu_;
    //     // }
    //     // else
    //     // {
    //     //     render_surface_mesh_=new SceneObjectDeformable(volumetric_surface_mesh_file_name_);
    //     // }
    // }
    // else
    // {
        render_surface_mesh_= new SceneObjectDeformable(volumetric_surface_mesh_file_name_);
    // }
    if(enable_textures_)
        visual_mesh_->SetUpTextures(SceneObject::MODULATE,SceneObject::NOMIPMAP);
    visual_mesh_->ResetDeformationToRest();
    visual_mesh_->BuildNeighboringStructure();
    visual_mesh_->BuildNormals();

    u_render_surface_=new double[3*visual_mesh_->Getn()];//---------------to check
    memset(u_render_surface_,0.0,sizeof(double)*3*visual_mesh_->Getn());
    render_volumetric_mesh_=new RenderVolumetricMesh();
    //load interpolation structure for object
    if(strcmp(object_interpolation_file_name_,"none")==0)
    {
        std::cout<<"Error: no object interpolation filename specified.\n";
        exit(1);
    }
    object_interpolation_element_vertices_num_=VolumetricMesh::getNumInterpolationElementVertices(object_interpolation_file_name_);
    if(object_interpolation_element_vertices_num_<0)
    {
        std::cout<<"Error: unable to open file "<<object_interpolation_file_name_<<".\n";
        exit(1);
    }
    std::cout<<"Num interpolation element vertices: "<<object_interpolation_element_vertices_num_<<".\n";
    VolumetricMesh::loadInterpolationWeights(object_interpolation_file_name_,visual_mesh_->Getn(),object_interpolation_element_vertices_num_,
                                            &object_interpolation_vertices_,&object_interpolation_weights_);
    //load fixed vertices only for fullspace
    if(simulation_mode_==FULLSPACE)
    {
        if(strcmp(fixed_vertex_file_name_,"none")==0)
        {
            std::cout<<"Error: no fixed vertices file name specified.\n";
        //    exit(1);
        }
        else
        {
            std::cout<<"load fixed vertex file:\n";
            if(simulator_->loadFixedVertices(fixed_vertex_file_name_))
            {
                std::cout<<"load fixed vertices succeed.\n";
            }
        }
		//load initial vel
		if(strcmp(initial_velocity_file_name_,"none")!=0)
		{
            simulator_->setInitialVelFilename(initial_velocity_file_name_);
		}
		//load initial velocity
		if(strcmp(initial_position_file_name_,"none")!=0)
		{
            simulator_->setInitialPosFilename(initial_position_file_name_);
		}

    }
    //load initial force
    if(strcmp(force_loads_file_name_,"none")!=0)
    {
        simulator_->setInitialForceFilename(force_loads_file_name_);
    }
    //load local example affected vertices file for object and examples
    if(strcmp(object_affected_vertices_file_name_,"none")!=0)
    {
        simulator_->setObjectAffectedVerticesFilename(object_affected_vertices_file_name_);
    }
    //load inital mesh
    if(strcmp(initial_tetmesh_file_name_,"none")!=0)
    {
        simulator_->setInitialTetMeshFilename(initial_tetmesh_file_name_);
    }
    if(example_num_>0)
    {
        //load local example affected vertices file for examples
        if(strcmp(example_affected_vertices_file_base_,"none")!=0)
        {
            simulator_->setExampleAffectedVerticesFilebase(example_affected_vertices_file_base_);
        }
    }

    //load example volumetric meshes
    if(example_num_>0)
    {
        simulator_->setExampleNum(example_num_);
        example_mesh_=new VolumetricMesh *[example_num_];
        if(simulator_->loadExamples(example_file_name_prefix_,example_num_))
        {
            for(unsigned int i=0;i<example_num_;++i)
                example_mesh_[i]=simulator_->exampleMesh(i);
        }
        else
        {
            std::cout<<"Error: unable to load the example simulation mesh from "<<simulation_mesh_file_name_<<".\n";
            exit(1);
        }
        current_example_index_=1;
        current_example_mesh_=example_mesh_[0];
        //set interpolate_eigenfunction num for simulator;
        simulator_->setInterpolateEigenfunctionNum(interpolate_eigenfunction_num_);
        //set LB cubica file name
        if(isload_LB_cubica_)
        {
            if(strcmp(object_LB_cubica_file_name_,"none")==0)
            {
                std::cout<<"Error: object cubica file for LB subspace basis unloaded.\n";
                exit(0);
            }
            loadLBObjectCubicaData(0);
        }
    }

    simulator_->setupSimulation();
    std::cout<<"init Simulation finish.\n";
}

void OpenGLDriver::initCamera()
{
    // std::cout<<"initCamera function:\n";
    double MY_PI = 3.141592653589793238462643;
    double up_dir[3] = {0,1,0};
    znear_ = camera_radius_*0.01;
    zfar_ = camera_radius_*100;
    camera_ = new SphericalCamera(camera_radius_,camera_longitude_/360*2*MY_PI,camera_lattitude_/360*2*MY_PI,
                                  focus_position_,up_dir);
    camera_->SetOrigin(focus_position_);
}

void OpenGLDriver::initGraphics()
{
    std::cout<<"initGraphics function:\n";
    // clear to white
    //glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.5);
    // clear to light blue
    //    glClearColor(233.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);
    // clear to gray
    //glClearColor(196.0 / 256, 196.0 / 256, 196.0 / 256, 0.0);
    // clear to brown
    //glClearColor(255.0 / 256, 156.0 / 256, 17.0 / 256, 0.0);
    // clear to medical blue
    glClearColor(148.0 / 256, 199.0 / 256, 211.0 / 256, 0.0);
    //clear to black
    //glClearColor(0.0,0.0,0.0,0.0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    reshapeFunction(window_width_,window_height_);

    //init camera
    initCamera();
    std::cout<<"Graphics initialization complete.\n";
}

void OpenGLDriver::displayFunction()
{
    PerformanceCounter display_counter;
    display_counter.StartCounter();
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    active_instance->camera_->Look();
    //set OpenGL Lighting
    if(active_instance->isrender_surface_mesh_)
    {
        if(active_instance->render_mesh_type_==VISUAL_MESH)
        {
            // if(active_instance->simulation_mode_==FULLSPACE)
                active_instance->visual_mesh_->SetLighting(active_instance->lighting_);
            // else
            //     active_instance->render_reduced_surface_mesh_->SetLighting(active_instance->lighting_);
        }
        glEnable(GL_LIGHTING);
        glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE);
        glStencilFunc(GL_ALWAYS,0,~(0u));
    }
    //show frame_rate on the left top of the window
    if(active_instance->render_fps_)
    {
        std::stringstream adaptor;
        adaptor.precision(2);
        std::string str;
        adaptor<<active_instance->fps_;
        str=std::string("FPS:") + adaptor.str();
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glDisable(GL_LIGHTING);
        glDisable(GL_TEXTURE_2D);
        gluOrtho2D(0,active_instance->window_width_,0,active_instance->window_height_);
        glColor3f(1.0,0.0,0.0);
        glRasterPos2i(5,active_instance->window_height_-19);
        for(int i=0;i<str.length();++i)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,str[i]);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        glEnable(GL_LIGHTING);
    }
    //render eigenfunction
    //  if(active_instance->render_eigenfunction_)
    //  {
    //     active_instance->drawIndexColorTable();//draw color table at left bottom corner of the window
    //  }

    //render extra objects
    if(active_instance->extra_objects_num_>0)//render the extra objects in sceneObjectDeformable
    {
        for(int i=0;i<active_instance->extra_objects_num_;++i)
            active_instance->extra_objects_[i]->Render();
    }

    glStencilFunc(GL_ALWAYS,1,~(0u));
    //render different mode
    if(active_instance->isrender_surface_mesh_&&active_instance->render_mesh_type_==VISUAL_MESH)
    {//object surface mesh
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0,1.0);
        //  active_instance->visual_mesh_->Render();
        // if(active_instance->simulation_mode_==FULLSPACE)
                active_instance->visual_mesh_->Render();
        // else
        // {
            // if(active_instance->with_constrains_)
            //     active_instance->render_reduced_surface_mesh_->Render();
            // else
                // active_instance->render_surface_mesh_->Render();
        // }

        if(active_instance->render_vertices_)
        {
            glDisable(GL_LIGHTING);
            glColor3f(0.5,0.0,0.0);
            glPointSize(8.0);
        //    if(active_instance->simulation_mode_==FULLSPACE)
                active_instance->visual_mesh_->RenderVertices();
            // else
            //     active_instance->render_reduced_surface_mesh_->RenderVertices();
            glEnable(GL_LIGHTING);
        }
        if(active_instance->render_wireframe_)
        {
            glDisable(GL_LIGHTING);
            glColor3f(0.0,0.0,0.0);
            // if(active_instance->simulation_mode_==FULLSPACE)
                // active_instance->visual_mesh_->RenderEdges();
            // else
            //     active_instance->render_reduced_surface_mesh_->RenderEdges();
                active_instance->render_surface_mesh_->RenderEdges();
            glEnable(GL_LIGHTING);
        }
        glDisable(GL_BLEND);
        active_instance->isrender_volumetric_mesh_=false;
    }
    else if(active_instance->isrender_volumetric_mesh_&&active_instance->render_mesh_type_==OBJECT_EIGEN_MESH)
    {//object eigen volumetric mesh
        if(active_instance->simulation_mesh_==NULL)
        {
            std::cout<<"Error: object simulation mesh is null.\n";
            exit(1);
        }
        else
        {
            // memcpy(active_instance->u_,active_instance->integrator_base_->Getq(),sizeof(double)*3*active_instance->simulation_vertices_num_);

            if((active_instance->render_eigenfunction_)&&(active_instance->isload_object_eigen_))
            {
                active_instance->render_volumetric_mesh_->RenderVertexColorMap(active_instance->simulation_mesh_,
                                active_instance->simulator_->objectEigenFunctions()[active_instance->current_render_eigen_idx_-1],active_instance->u_);
            }
            else
            {
                glDisable(GL_LIGHTING);
                glColor3f(0.0,0.5,0.0);
                active_instance->render_volumetric_mesh_->RenderDeformation(active_instance->simulation_mesh_,active_instance->u_);
                // active_instance->render_volumetric_mesh_->Render(active_instance->simulation_mesh_);
                if(active_instance->render_vertices_)
                {
                    glDisable(GL_LIGHTING);
                    glColor3f(0.5,0.0,0.0);
                    glPointSize(8.0);
                    for(int i=0;i<active_instance->simulation_vertices_num_;++i)
                        active_instance->render_volumetric_mesh_->RenderVertexDeformed(active_instance->simulation_mesh_,i,active_instance->u_);
                    // active_instance->render_volumetric_mesh_->RenderVertices(active_instance->simulation_mesh_);
                    glEnable(GL_LIGHTING);
                }
                if(active_instance->render_wireframe_)
                {
                    glDisable(GL_LIGHTING);
                    glColor3f(0.0,0.0,0.0);
                    active_instance->render_volumetric_mesh_->RenderSolidAndWireframeDeformation(active_instance->simulation_mesh_,active_instance->u_);
                    // active_instance->render_volumetric_mesh_->RenderWireframe(active_instance->simulation_mesh_);
                    glEnable(GL_LIGHTING);
                }
                glDisable(GL_BLEND);
            }
        }
        active_instance->isrender_surface_mesh_=false;
    }
    else if(active_instance->isrender_volumetric_mesh_&&active_instance->render_mesh_type_==EXAMPLE_MESH)
    {//EXAMPLE_MESH
        if(active_instance->example_mesh_[active_instance->current_example_index_-1]==NULL)
        {
            std::cout<<"Error: object simulation mesh is null.\n";
            exit(1);
        }
        else
        {
            if((active_instance->render_eigenfunction_)&&(active_instance->isload_example_eigen_))
            {
                active_instance->render_volumetric_mesh_->RenderVertexColorMap(active_instance->example_mesh_[active_instance->current_example_index_-1],
                                active_instance->simulator_->exampleEigenFunctions()[active_instance->current_example_index_-1][active_instance->current_render_eigen_idx_-1]);
            }
            else
            {
                glDisable(GL_LIGHTING);
                glColor3f(0.0,0.5,0.0);
                active_instance->render_volumetric_mesh_->Render(active_instance->example_mesh_[active_instance->current_example_index_-1]);
                //active_instance->render_volumetric_mesh_->RenderDeformation(active_instance->example_mesh_[active_instance->current_example_index_-1],active_instance->simulator_->exampleDis()[active_instance->current_example_index_-1]);
                if(active_instance->render_vertices_)
                {
                    glDisable(GL_LIGHTING);
                    glColor3f(0.5,0.0,0.0);
                    glPointSize(8.0);
                    active_instance->render_volumetric_mesh_->RenderVertices(active_instance->example_mesh_[active_instance->current_example_index_-1]);
                    // for(int i=0;i<active_instance->example_mesh_[active_instance->current_example_index_-1]->getNumVertices();++i)
                    //     active_instance->render_volumetric_mesh_->RenderVertexDeformed(active_instance->example_mesh_[active_instance->current_example_index_-1],i,active_instance->simulator_->exampleDis()[active_instance->current_example_index_-1]);
                    glEnable(GL_LIGHTING);
                }
                if(active_instance->render_wireframe_)
                {
                    glDisable(GL_LIGHTING);
                    glColor3f(0.0,0.0,0.5);
                    active_instance->render_volumetric_mesh_->RenderWireframe(active_instance->example_mesh_[active_instance->current_example_index_-1]);
                    //active_instance->render_volumetric_mesh_->RenderSolidAndWireframeDeformation(active_instance->example_mesh_[active_instance->current_example_index_-1],active_instance->simulator_->exampleDis()[active_instance->current_example_index_-1]);
                    glEnable(GL_LIGHTING);
                }
                glDisable(GL_BLEND);
            }
        }
        active_instance->isrender_surface_mesh_=false;
    }
    glStencilFunc(GL_ALWAYS,0,~(0u));
    glDisable(GL_TEXTURE_2D);
    //render planes
    if(active_instance->plane_num_>0)
    {
        active_instance->planes_->render();
    }
    glDisable(GL_LIGHTING);
    glStencilFunc(GL_ALWAYS,0,~(0u));
    glColor3f(1.0,0.0,0.0);
    glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP);
    //render axis
    if(active_instance->render_axis_)
    {
        glLineWidth(1.0);
        active_instance->drawAxis(10.0);
    }

    //render the currently pulled vertex
    if(active_instance->pulled_vertex_!=-1)
    {
        glColor3f(1.0,0.3,0.0);
        Vec3d pulled_vertex_pos;
        // if(active_instance->simulation_mode_==FULLSPACE)
            // active_instance->visual_mesh_->GetSingleVertexPositionFromBuffer(active_instance->pulled_vertex_,
            //                                 &pulled_vertex_pos[0],&pulled_vertex_pos[1],&pulled_vertex_pos[2]);
            active_instance->render_surface_mesh_->GetSingleVertexPositionFromBuffer(active_instance->pulled_vertex_,
                                            &pulled_vertex_pos[0],&pulled_vertex_pos[1],&pulled_vertex_pos[2]);
        // else
        //     active_instance->render_reduced_surface_mesh_->GetSingleVertexPositionFromBuffer(active_instance->pulled_vertex_,
        //                                     &pulled_vertex_pos[0],&pulled_vertex_pos[1],&pulled_vertex_pos[2]);
        // for(int i=0;i<3;++i)
        //     pulled_vertex_pos[i]=(*active_instance->render_reduced_surface_mesh_->getVertex(active_instance->pulled_vertex_))[i];

        glEnable(GL_POLYGON_OFFSET_POINT);
        glPolygonOffset(-1.0,-1.0);
        glPointSize(7.0);
        glBegin(GL_POINTS);
        glVertex3f(pulled_vertex_pos[0],pulled_vertex_pos[1],pulled_vertex_pos[2]);
        glEnd();
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
    //render fixed vertices--for volumetric mesh
    if(active_instance->render_fixed_vertices_)
    {
        glColor3f(1,0,0);
        for(int i=0;i<active_instance->simulator_->fixedVertexNum();++i)
        {
            Vec3d * vertex = active_instance->simulation_mesh_->getVertex(active_instance->simulator_->fixedVertexPtr()[i]);
            glEnable(GL_POLYGON_OFFSET_POINT);
            glPolygonOffset(-1.0,-1.0);
            glPointSize(10.0);
            glBegin(GL_POINTS);
            glVertex3f((*vertex)[0],(*vertex)[1],(*vertex)[2]);
            glEnd();
            glDisable(GL_POLYGON_OFFSET_FILL);
        }
    }

    //glStencilFunc(GL_ALWAYS,1,~(0u));
    //render vertex velocity
    if(active_instance->render_velocity_)
    {
        glDisable(GL_LIGHTING);
        for(int i=0;i<active_instance->simulation_vertices_num_;++i)
        {
            Vec3d vert_pos,vert_new_pos;
            for(int j=0;j<3;++j)
            {
                vert_pos[j]=(*active_instance->simulation_mesh_->getVertex(i))[j]+active_instance->u_[3*i+j];
                vert_new_pos[j]=vert_pos[j]+active_instance->simulator_->getvel()[3*i+j]*active_instance->render_velocity_scale_;
            }
            glColor3f(1.0,0.3,0.0);
            glLineWidth(1.0);
            glBegin(GL_LINES);
            glVertex3f(vert_pos[0],vert_pos[1],vert_pos[2]);
            glVertex3f(vert_new_pos[0],vert_new_pos[1],vert_new_pos[2]);
            glEnd();
        }
        glEnable(GL_LIGHTING);
    }
    //render plane collision force
    if(active_instance->render_ext_force_)
    {
        glDisable(GL_LIGHTING);
        for(int i=0;i<active_instance->simulation_vertices_num_;++i)
        {
            Vec3d vert_pos,vert_new_pos;
            for(int j=0;j<3;++j)
            {
                vert_pos[j]=(*active_instance->simulation_mesh_->getVertex(i))[j]+active_instance->u_[3*i+j];
                vert_new_pos[j]=vert_pos[j]+active_instance->simulator_->getExternalForce()[3*i+j]*0.001;
            }
            glColor3f(1.0,0.3,1.0);
            glLineWidth(2.0);
            glBegin(GL_LINES);
            glVertex3f(vert_pos[0],vert_pos[1],vert_pos[2]);
            glVertex3f(vert_new_pos[0],vert_new_pos[1],vert_new_pos[2]);
            glEnd();
        }
        glEnable(GL_LIGHTING);
    }
    //render vertex displacement
    // if(active_instance->render_dis_)
    // {
    //     glDisable(GL_LIGHTING);
    //     for(int i=0;i<active_instance->example_mesh_[0]->getNumVertices();++i)
    //     {
    //         Vec3d vert_pos,vert_new_pos;
    //         for(int j=0;j<3;++j)
    //         {
    //             vert_pos[j]=(*active_instance->example_mesh_[0]->getVertex(i))[j];
    //             vert_new_pos[j]=vert_pos[j]+active_instance->simulator_->exampleDis()[active_instance->current_example_index_-1][3*i+j]*2.0;
    //         }
    //         glColor3f(1.0,0.3,1.0);
    //         glLineWidth(1.0);
    //         glBegin(GL_LINES);
    //         glVertex3f(vert_pos[0],vert_pos[1],vert_pos[2]);
    //         glVertex3f(vert_new_pos[0],vert_new_pos[1],vert_new_pos[2]);
    //         glEnd();
    //     }
    //     glEnable(GL_LIGHTING);
    // }
    //render example guided surface deformation---to do
    display_counter.StopCounter();
    active_instance->total_display_time_+=display_counter.GetElapsedTime();
    if(active_instance->time_step_counter_==(active_instance->total_steps_-1))
    {
        std::cout<<"Total display time is: "<<active_instance->total_display_time_<<" s.";
        std::cout<<"Average is: "<<(active_instance->total_display_time_)/(active_instance->total_frames_)<<" s/frame. \n";
    }

    glutSwapBuffers();
}
void OpenGLDriver::idleFunction()
{
    // std::cout<<"idleFunction\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    // memcpy(active_instance->u_,active_instance->simulator_->getu(),sizeof(double)*3*active_instance->simulation_vertices_num_);
    // getchar();
	// std::cout<<"6\n";
    // getchar();
    glutSetWindow(active_instance->window_id_);
	// std::cout<<"7\n";
    // getchar();
    if(active_instance->time_step_counter_<active_instance->total_steps_)
    {
        //time statistics by frame
        PerformanceCounter each_frame_performance_counter;
        double each_frame_time=0.0;
        each_frame_performance_counter.StartCounter();
    	// std::cout<<"8\n";
        // getchar();
        memcpy(active_instance->u_,active_instance->simulator_->getu(),sizeof(double)*3*active_instance->simulation_vertices_num_);
    	// std::cout<<"9\n";
        // getchar();
        if(!active_instance->pause_simulation_)
        {
            // std::cout<<"reduced:testEnergyGradients:\n";
            // active_instance->reduced_neoHookean_force_model_->testEnergyGradients();
            // getchar();
            // std::cout<<"testInternalForceGradients:\n";
            // active_instance->reduced_neoHookean_force_model_->testObjectiveGradients();
            // std::cout<<"testInternalForceGradients-end:\n";
            // getchar();
            // active_instance->simulator_->testEnergyGradients();
            // getchar();
            // active_instance->simulator_->testObjectiveGradients();
            // getchar();
            // std::cout<<"-------------"<<active_instance->simulation_mode_<<"----------\n";

            memset(active_instance->f_ext_,0.0,sizeof(double)*3*active_instance->simulation_vertices_num_);

        	// std::cout<<"10\n";
            // getchar();
            //time statistics by step, time_step_intervals is 100
        	PerformanceCounter step_counter;
            step_counter.StartCounter();
            // if(active_instance->simulation_mode_==REDUCEDSPACE)
            // {
            //     if(active_instance->left_button_down_)
            //     {
            //         std::cout<<"pulled_vertex_:"<<active_instance->pulled_vertex_<<"\n";
            //         if(active_instance->pulled_vertex_!=-1)
            //         {
            //             double force_x=active_instance->mouse_pos_[0]-active_instance->drag_start_x_;
            //             double force_y=(-1.0)*(active_instance->mouse_pos_[1]-active_instance->drag_start_y_);
            //             double external_force[3];
            //             active_instance->camera_->CameraVector2WorldVector_OrientationOnly3D(force_x,force_y,0,external_force);
            //             std::cout<<active_instance->pulled_vertex_<<" fx: "<<force_x<<",fy: "<<force_y<<" | "<<external_force[0]<<",";
            //             std::cout<<"external_force:"<<external_force[0]<<","<<external_force[1]<<","<<external_force[2]<<std::endl;
            //             for(int i=0;i<3;++i)
            //             {
            //                 external_force[i]*=active_instance->deformable_object_compliance_;
            //             }
            //             //register force on the pulled vertex
            //             active_instance->f_ext_[3*(active_instance->pulled_vertex_)+0]=external_force[0];
            //             active_instance->f_ext_[3*(active_instance->pulled_vertex_)+1]=external_force[1];
            //             active_instance->f_ext_[3*(active_instance->pulled_vertex_)+2]=external_force[2];
            //             //distributing force over the neighboring vertices
            //             set<int> affected_vertices;
            //             set<int> last_layer_vertices;
            //             affected_vertices.insert(active_instance->pulled_vertex_);
            //             last_layer_vertices.insert(active_instance->pulled_vertex_);
            //             for(int j=1;j<active_instance->force_neighbor_size_;++j)
            //             {
            //                 //linear kernel
            //                 double force_magnitude=1.0*(active_instance->force_neighbor_size_-j)/active_instance->force_neighbor_size_;
            //                 set<int> new_affected_vertices;
            //                 for(set<int>::iterator iter=last_layer_vertices.begin();iter!=last_layer_vertices.end();++iter)
            //                 {
            //                     //traverse all neighbors and check if they were already previously inserted
            //                     int vtx=*iter;
            //                     int deg=active_instance->mesh_graph_->GetNumNeighbors(vtx);
            //                     for(int k=0;k<deg;++k)
            //                     {
            //                         int vtx_neighbor=active_instance->mesh_graph_->GetNeighbor(vtx,k);
            //                         if(affected_vertices.find(vtx_neighbor)==affected_vertices.end())
            //                             new_affected_vertices.insert(vtx_neighbor);
            //                     }
            //                 }
            //                 last_layer_vertices.clear();
            //                 for(set<int>::iterator iter=new_affected_vertices.begin();iter!=new_affected_vertices.end();iter++)
            //                 {
            //                     //apply forces
            //                     active_instance->f_ext_[3*(*iter)+0]+=force_magnitude*external_force[0];
            //                     active_instance->f_ext_[3*(*iter)+1]+=force_magnitude*external_force[1];
            //                     active_instance->f_ext_[3*(*iter)+2]+=force_magnitude*external_force[2];
            //                     //generate new layers
            //                     last_layer_vertices.insert(*iter);
            //                     affected_vertices.insert(*iter);
            //                 }
            //             }
            //         }
            //     }
            //
            //     //plane--
            //
            // 	// std::cout<<"11\n";
                // if(!active_instance->with_constrains_)
                // {
                    if(active_instance->plane_num_>0)
                	{
                        //plane:collision-method 1
                		active_instance->planes_->resolveContact(active_instance->render_surface_mesh_->GetMesh(),active_instance->f_col_,active_instance->collide_vert_num_);
                        for(int i=0;i<3*active_instance->simulation_vertices_num_;++i)
                        {
                            active_instance->f_ext_[i]+=active_instance->f_col_[i];
                            // if(active_instance->f_col_[i]>1.0e-6)
                            //     std::cout<<active_instance->f_col_[i]<<",";
                        }
                        //plane:return new dis and new vel-method 2
                        // active_instance->planes_->resolveContact(active_instance->render_surface_mesh_->GetMesh(),//active_instance->f_col_,
                        //                                 active_instance->simulator_->getvel(),active_instance->collide_u_,
                        //                                 active_instance->collide_vel_,active_instance->collide_vert_num_);
                        // for(int i=0;i<3*active_instance->simulation_mesh_->getNumVertices();++i)
                        // {
                        //     // active_instance->u_[i]+=active_instance->collide_u_[i];
                        //     active_instance->f_ext_[i]=(active_instance->collide_vel_[i]-active_instance->simulator_->getvel()[i])/active_instance->time_step_;
                        // }
                        // active_instance->simulator_->setu(active_instance->u_);
                        // active_instance->simulator_->setVelAfterCollision(active_instance->collide_vel_);
                        active_instance->simulator_->setCollisionNum(active_instance->collide_vert_num_[0]);
                    }
                // }
            //
            // 	// std::cout<<"12\n";
                // active_instance->simulator_->setExternalForces(active_instance->f_ext_);
            // 	// std::cout<<"13\n";
                // active_instance->simulator_->advanceStep();
            //
            // 	// std::cout<<"14\n";
            //     // getchar();
            // }
            // else
            // {
                //drag force
                if(active_instance->left_button_down_)
                {
                    std::cout<<"pulled_vertex_:"<<active_instance->pulled_vertex_<<"\n";
                    if(active_instance->pulled_vertex_!=-1)
                    {
                        double force_x=active_instance->mouse_pos_[0]-active_instance->drag_start_x_;
                        double force_y=(-1.0)*(active_instance->mouse_pos_[1]-active_instance->drag_start_y_);
                        double external_force[3];
                        active_instance->camera_->CameraVector2WorldVector_OrientationOnly3D(force_x,force_y,0,external_force);
                        std::cout<<"external_force:"<<external_force[0]<<","<<external_force[1]<<","<<external_force[2]<<std::endl;

                        for(int i=0;i<3;++i)
                        {
                            external_force[i]*=active_instance->deformable_object_compliance_;
                        }
                        std::cout<<active_instance->pulled_vertex_<<" fx: "<<force_x<<",fy: "<<force_y<<" | "<<external_force[0]<<",";

                        std::cout<<"external_force:"<<external_force[0]<<","<<external_force[1]<<","<<external_force[2]<<std::endl;
                        //register force on the pulled vertex
                        active_instance->f_ext_[3*(active_instance->pulled_vertex_)+0]+=external_force[0];
                        active_instance->f_ext_[3*(active_instance->pulled_vertex_)+1]+=external_force[1];
                        active_instance->f_ext_[3*(active_instance->pulled_vertex_)+2]+=external_force[2];
                        //distributing force over the neighboring vertices
                        set<int> affected_vertices;
                        set<int> last_layer_vertices;
                        affected_vertices.insert(active_instance->pulled_vertex_);
                        last_layer_vertices.insert(active_instance->pulled_vertex_);
                        for(int j=1;j<active_instance->force_neighbor_size_;++j)
                        {
                            //linear kernel
                            double force_magnitude=1.0*(active_instance->force_neighbor_size_-j)/active_instance->force_neighbor_size_;
                            set<int> new_affected_vertices;
                            for(set<int>::iterator iter=last_layer_vertices.begin();iter!=last_layer_vertices.end();++iter)
                            {
                                //traverse all neighbors and check if they were already previously inserted
                                int vtx=*iter;
                                int deg=active_instance->mesh_graph_->GetNumNeighbors(vtx);
                                for(int k=0;k<deg;++k)
                                {
                                    int vtx_neighbor=active_instance->mesh_graph_->GetNeighbor(vtx,k);
                                    if(affected_vertices.find(vtx_neighbor)==affected_vertices.end())
                                        new_affected_vertices.insert(vtx_neighbor);
                                }
                            }
                            last_layer_vertices.clear();
                            for(set<int>::iterator iter=new_affected_vertices.begin();iter!=new_affected_vertices.end();iter++)
                            {
                                //apply forces
                                active_instance->f_ext_[3*(*iter)+0]+=force_magnitude*external_force[0];
                                active_instance->f_ext_[3*(*iter)+1]+=force_magnitude*external_force[1];
                                active_instance->f_ext_[3*(*iter)+2]+=force_magnitude*external_force[2];
                                //generate new layers
                                last_layer_vertices.insert(*iter);
                                affected_vertices.insert(*iter);
                            }
                        }
                    }
                }

                //plane--
                // if(active_instance->plane_num_>0)
            	// {
            	// 	active_instance->planes_->resolveContact(active_instance->render_surface_mesh_->GetMesh(),active_instance->f_col_,active_instance->collide_vert_num_);
                //     for(int i=0;i<3*active_instance->simulation_vertices_num_;++i)
                //     {
                //         active_instance->f_ext_[i]+=active_instance->f_col_[i];
                //     }
            	// }
                // std::cout<<"a\n";
                active_instance->simulator_->setExternalForces(active_instance->f_ext_);
                // std::cout<<"b\n";
                active_instance->simulator_->advanceStep();
                // std::cout<<"c\n";
            // }
            step_counter.StopCounter();
            active_instance->step_simulation_time_+=step_counter.GetElapsedTime();
            //time computing
            // if((active_instance->time_step_counter_)%active_instance->timer_sample_interval_==0&&active_instance->time_step_counter_>=active_instance->timer_sample_interval_)
            // {
            //     printf("Time in %d steps: %f, average step time: %f\n",active_instance->timer_sample_interval_,active_instance->step_simulation_time_,active_instance->step_simulation_time_/active_instance->timer_sample_interval_);
            //     active_instance->step_simulation_time_ = 0.0;
            //     getchar();
            // }
        }
        // std::cout<<"1\n";
        // getchar();
    	// std::cout<<"15\n";
        // getchar();
        if(active_instance->render_mesh_type_==VISUAL_MESH)
        {
            memset(active_instance->u_render_surface_,0.0,sizeof(double)*3*(active_instance->visual_mesh_->Getn()));
            VolumetricMesh::interpolate(active_instance->u_,active_instance->u_render_surface_,active_instance->visual_mesh_->Getn(),
                                        active_instance->object_interpolation_element_vertices_num_,active_instance->object_interpolation_vertices_,
                                        active_instance->object_interpolation_weights_);

            active_instance->visual_mesh_->SetVertexDeformations(active_instance->u_render_surface_);
            active_instance->render_surface_mesh_->SetVertexDeformations(active_instance->u_);
        }
        // active_instance->simulator_->setu(active_instance->u_);
    	// std::cout<<"16\n";
        // getchar();
        each_frame_performance_counter.StopCounter();
        each_frame_time=each_frame_performance_counter.GetElapsedTime();
        active_instance->total_simulation_time_+=each_frame_time;
        // if(active_instance->time_step_counter_==(active_instance->total_steps_-1))
        // {
        //     std::cout<<"Total simulation_time(write file time excluded) is "<<active_instance->total_simulation_time_<<" s;";
        //     std::cout<<"Average: "<<(active_instance->total_simulation_time_)/(active_instance->total_frames_)<<" s/frame.\n";
        // }
        if(each_frame_time>0)
            active_instance->fps_=1.0/each_frame_time;
        else
            active_instance->fps_=0.0;
        if(!active_instance->pause_simulation_)
            active_instance->time_step_counter_++;

        //save files per 4 steps
        if((active_instance->time_step_counter_%4==0)&&(active_instance->time_step_counter_>0)&&(active_instance->enable_save_objmesh_)&&(!active_instance->pause_simulation_))
        {
            active_instance->saveCurrentObjmesh(0);
        }
        // if((active_instance->time_step_counter_%4==0)&&(active_instance->time_step_counter_>0)&&(active_instance->save_tet_mesh_)&&(!active_instance->pause_simulation_))
        // {
        //     active_instance->saveCurrentTetmesh(0);
        // }
        // active_instance->simulator_->setu(active_instance->u_);
        glutPostRedisplay();
    }

}

void OpenGLDriver::reshapeFunction(int width, int height)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    glViewport(0,0,width,height);
    active_instance->window_width_ = width;
    active_instance->window_height_ = height;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f,1.0f*width/height,active_instance->znear_,active_instance->zfar_);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void OpenGLDriver::keyboardFunction(unsigned char key, int x, int y)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    switch(key)
    {
    case 27: //ESC, exit
        exitApplication(0);
        break;
    case 'z': //single step
        active_instance->simulator_->advanceStep();
        break;
    case 'a': //render axis
        active_instance->render_axis_ = !(active_instance->render_axis_);
        break;
    case 'w': //render wireframe
        active_instance->render_wireframe_ = !(active_instance->render_wireframe_);
        break;
    case 'b': //render fixed vertices
        active_instance->render_fixed_vertices_ = !(active_instance->render_fixed_vertices_);
        break;
    case 'f':
        active_instance->render_ext_force_ = !(active_instance->render_ext_force_);
        break;
    case 'g': //switch on/off gravity
        active_instance->add_gravity_=!(active_instance->add_gravity_);
        addGravitySwitch(active_instance->add_gravity_);
        // active_instance->simulator_->setGravity(active_instance->add_gravity_,active_instance->gravity_);
        break;
    case 'r': //reset camera
        active_instance->camera_->Reset();
        break;
    case 's':
        // active_instance->save_tet_mesh_ = !(active_instance->save_tet_mesh_ );
        active_instance->enable_save_objmesh_ = !(active_instance->enable_save_objmesh_);
        active_instance->save_tet_mesh_ = !(active_instance->save_tet_mesh_ );
        // active_instance->saveCurrentObjmesh(0);
        break;
    case 't':
        active_instance->save_tet_mesh_ = !(active_instance->save_tet_mesh_ );
        break;
    // case 'c':
    //     active_instance->save_reconstruct_mesh_ = !(active_instance->save_tet_mesh_ );
    //     active_instance->enable_save_objmesh_ = !(active_instance->enable_save_objmesh_);
    //     break;
    case 32: //space button, pause simulation
        active_instance->pause_simulation_ = !(active_instance->pause_simulation_);
        break;
    case 'v': //render vertices
        active_instance->render_vertices_ = !(active_instance->render_vertices_);
        break;
    case 'V':
        active_instance->render_velocity_ = !(active_instance->render_velocity_);
        break;
    case 'i':
        active_instance->render_velocity_scale_ *=2.0;
        break;
    case 'd':
        active_instance->render_velocity_scale_ *= 0.5;
        break;
    // case 'u':
    //     active_instance->render_dis_ = !(active_instance->render_dis_);
    //     break;
    case 'e': //render eigenfunctions
        active_instance->render_eigenfunction_ = !(active_instance->render_eigenfunction_);
        std::string static_text_content("Rendering eigenfunctions enabled: ");
        if((active_instance->render_mesh_type_==OBJECT_EIGEN_MESH)&&(active_instance->isload_object_eigen_))
        {
            static_text_content+="Yes";
            active_instance->render_volumetric_mesh_->RenderVertexColorMap(active_instance->simulation_mesh_,
                            active_instance->simulator_->objectEigenFunctions()[active_instance->current_render_eigen_idx_-1]);
        }
        else if((active_instance->render_mesh_type_==EXAMPLE_MESH)&&(active_instance->isload_example_eigen_))
        {
            active_instance->render_volumetric_mesh_->RenderVertexColorMap(active_instance->example_mesh_[active_instance->current_example_index_-1],
                            active_instance->simulator_->exampleEigenFunctions()[active_instance->current_example_index_-1][active_instance->current_render_eigen_idx_-1]);
        }
        else
        {
            static_text_content+="No";
            active_instance->render_eigenfunction_=false;
            std::cout<<"Error: eigenFunctions unloaded.\n";
        }
        active_instance->glui_rendering_eigenfunctions_enabled_->set_name(static_text_content.c_str());
        break;
    }
}

void OpenGLDriver::specialFunction(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_LEFT:
        break;
    case GLUT_KEY_RIGHT:
        break;
    case GLUT_KEY_DOWN:
        break;
    case GLUT_KEY_UP:
        break;
    case GLUT_KEY_PAGE_UP:
        break;
    case GLUT_KEY_PAGE_DOWN:
        break;
    case GLUT_KEY_HOME:
        break;
    case GLUT_KEY_END:
        break;
    case GLUT_KEY_INSERT:
        break;
    default:
        break;
    }
}

void OpenGLDriver::motionFunction(int x, int y)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    int mouse_delta_x = x - active_instance->mouse_pos_[0];
    int mouse_delta_y = y - active_instance->mouse_pos_[1];
    active_instance->mouse_pos_[0] = x;
    active_instance->mouse_pos_[1] = y;
    SphericalCamera *camera = active_instance->camera_;
    if(active_instance->left_button_down_)
    {}
    if(active_instance->right_button_down_) //camera rotation
    {
        const double factor = 0.2;
        camera->MoveRight(factor*mouse_delta_x);
        camera->MoveUp(factor*mouse_delta_y);
    }
    if(active_instance->middle_button_down_) //zoom in/out
    {
        const double factor = 0.1;
        camera->ZoomIn(active_instance->camera_radius_*factor*mouse_delta_y);
    }
}

void OpenGLDriver::mouseFunction(int button, int state, int x, int y)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    switch (button)
     {
        case GLUT_LEFT_BUTTON:
            active_instance->left_button_down_=(state==GLUT_DOWN);
            active_instance->shift_pressed_=(glutGetModifiers()==GLUT_ACTIVE_SHIFT);
            active_instance->alt_pressed_=(glutGetModifiers()==GLUT_ACTIVE_ALT);
            active_instance->ctrl_pressed_=(glutGetModifiers()==GLUT_ACTIVE_CTRL);
            if(active_instance->left_button_down_&&(!active_instance->shift_pressed_)&&(!active_instance->ctrl_pressed_)) //used for pulled vertex, apply force
            {
                //apply force to vertex
                GLdouble model[16];
                glGetDoublev(GL_MODELVIEW_MATRIX,model);
                GLdouble proj[16];
                glGetDoublev(GL_PROJECTION_MATRIX,proj);
                GLint view[4];
                glGetIntegerv(GL_VIEWPORT,view);
                int win_x=x;
                int win_y=view[3]-1-y;
                float z_value;
                glReadPixels(win_x,win_y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&z_value);
                GLubyte stencil_value;
                glReadPixels(win_x,win_y,1,1,GL_STENCIL_INDEX,GL_UNSIGNED_BYTE,&stencil_value);
                GLdouble world_x,world_y,world_z;
                gluUnProject(win_x,win_y,z_value,model,proj,view,&world_x,&world_y,&world_z);
                if(stencil_value==1)
                {
                    active_instance->drag_start_x_=x;
                    active_instance->drag_start_y_=y;
                    Vec3d pos(world_x,world_y,world_z);
                    //the pulled vertex is on the exterior surface of the volumetric mesh
                    //virtual int GetClosestVertex(Vec3d & queryPos, double * distance=NULL, double * auxVertexBuffer=NULL);
                    //get pulled_vertex_ index for simulation mesh
                    if(active_instance->isrender_surface_mesh_)
                    {
                        // if(active_instance->simulation_mode_==REDUCEDSPACE)                    //
                        //     active_instance->pulled_vertex_=active_instance->render_reduced_surface_mesh_->GetClosestVertex(pos);
                        // else
                            active_instance->pulled_vertex_=active_instance->render_surface_mesh_->GetClosestVertex(pos);
                    }
                    if(active_instance->isrender_volumetric_mesh_)
                        active_instance->pulled_vertex_=active_instance->simulation_mesh_->getClosestVertex(pos);
                    std::cout<<"Clicked on vertex "<<active_instance->pulled_vertex_<<" (0-indexed)\n";
                }
                else
                    std::cout<<"Clicked on empty stencil: "<<stencil_value<<".\n";
            }
            if(!active_instance->left_button_down_)
                active_instance->pulled_vertex_=-1;
            break;
        case GLUT_MIDDLE_BUTTON:
            active_instance->middle_button_down_=(state==GLUT_DOWN);
            active_instance->shift_pressed_=(glutGetModifiers()==GLUT_ACTIVE_SHIFT);
            active_instance->alt_pressed_=(glutGetModifiers()==GLUT_ACTIVE_ALT);
            active_instance->ctrl_pressed_=(glutGetModifiers()==GLUT_ACTIVE_CTRL);
            break;
        case GLUT_RIGHT_BUTTON:
            active_instance->right_button_down_=(state==GLUT_DOWN);
            active_instance->shift_pressed_=(glutGetModifiers()==GLUT_ACTIVE_SHIFT);
            active_instance->alt_pressed_=(glutGetModifiers()==GLUT_ACTIVE_ALT);
            active_instance->ctrl_pressed_=(glutGetModifiers()==GLUT_ACTIVE_CTRL);
            break;
    }
    active_instance->mouse_pos_[0]=x;
    active_instance->mouse_pos_[1]=y;
}
void OpenGLDriver::addGravitySwitch(bool add_gravity)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    active_instance->simulator_->setGravity(active_instance->add_gravity_,active_instance->gravity_);

    if(active_instance->add_gravity_)
        std::cout<<"Gravity switch on.\n";
    else
        std::cout<<"Gravity switch off.\n";
}
void OpenGLDriver::updateRenderMesh(int code)
{
    // std::cout<<"updateRenderMesh function:\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(active_instance->render_mesh_type_==VISUAL_MESH)
    {
        active_instance->isrender_volumetric_mesh_=false;
        active_instance->isrender_surface_mesh_=true;
        active_instance->render_eigen_index_spinner_->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);
    }
    else if(active_instance->render_mesh_type_==OBJECT_EIGEN_MESH)
    {
        if(active_instance->isload_object_eigen_)
        {
            active_instance->render_eigen_index_spinner_->set_int_limits(1,active_instance->interpolate_eigenfunction_num_,
                                                            GLUI_LIMIT_CLAMP);
        }
        else
            active_instance->render_eigen_index_spinner_->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);
        active_instance->isrender_volumetric_mesh_=true;
        active_instance->isrender_surface_mesh_=false;
    }
    else if(active_instance->render_mesh_type_==EXAMPLE_MESH)
    {
        if(active_instance->isload_example_eigen_)
            active_instance->render_eigen_index_spinner_->set_int_limits(1,active_instance->interpolate_eigenfunction_num_,GLUI_LIMIT_CLAMP);
        else
            active_instance->render_eigen_index_spinner_->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);
        active_instance->isrender_volumetric_mesh_=true;
        active_instance->isrender_surface_mesh_=false;
    }
    else
    {
    }
}
void OpenGLDriver::updateCurrentExample(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    active_instance->current_example_mesh_=active_instance->example_mesh_[active_instance->current_example_index_-1];
    updateRenderMesh(0);
}
void OpenGLDriver::changeCurrentEigenIndex(int code)
{
    OpenGLDriver* active_instance=OpenGLDriver::activeInstance();
    assert(active_instance);
    active_instance->updateRenderMesh(code);
}
void OpenGLDriver::changeSimulationMode(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!(active_instance->isload_object_eigen_&&active_instance->isload_example_eigen_))
    {
        std::cout<<"Error:eigenfunctions are unloaded.\n";
        active_instance->enable_example_simulation_=false;
        return;
    }
    active_instance->enable_example_simulation_=!(active_instance->enable_example_simulation_);
    active_instance->simulator_->setEnableExampleBasedSimulation(active_instance->enable_example_simulation_);
    if(active_instance->enable_example_simulation_)
        active_instance->change_simulation_mode_button_->set_name("Disable example-based simulation");
    else
    {
        active_instance->change_simulation_mode_button_->set_name("Enable example-based simulation");
    }
}

// void OpenGLDriver::enableReducedSimulation(int code)
// {
//     OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
//     assert(active_instance);
//     active_instance->enable_reduced_simulation_=!(active_instance->enable_reduced_simulation_);
//     if(active_instance->enable_example_simulation_)
//         active_instance->reduced_simulation_button_->set_name("Disable reduced simulation");
//     else
//     {
//         active_instance->reduced_simulation_button_->set_name("Enable reduced simulation");
//     }
// }
void OpenGLDriver::loadMassmatrix(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadMassmatrix(active_instance->mass_matrix_file_name_))
    {
        std::cout<<"Error:load mass matrix files failed.\n";
        exit(0);
    }
    std::cout<<"Load object mass matrix succeed.\n";
}
void OpenGLDriver::loadInertiaTensor(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadInertiaTensor(active_instance->inertia_tensor_file_name_))
    {
        std::cout<<"Error:load inertia tensor files failed.\n";
        exit(0);
    }
    std::cout<<"Load object inertia tensor succeed.\n";
}
void OpenGLDriver::loadObjectEigenfunctions(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadObjectEigenfunctions(active_instance->object_eigen_file_name_))
    {
        std::cout<<"Error: load object eigenfunctions failed.\n";
        exit(0);
    }
    std::cout<<"Load object eigenfunctions succeed!\n";

    active_instance->isload_object_eigen_=true;
    active_instance->glui_object_surface_eigenfunctions_loaded_->set_name("Simulation object surface eigenfunctions loaded:Yes");
    if(active_instance->render_mesh_type_==OBJECT_EIGEN_MESH)
        active_instance->render_eigen_index_spinner_->set_int_limits(1,active_instance->interpolate_eigenfunction_num_,GLUI_LIMIT_CLAMP);

    //get initial object eigencoefs
    // for(int i=0;i<active_instance->reconstruct_eigenfunction_num_;++i)
    //     for(int j=0;j<3;++j)
    //         active_instance->initial_object_eigencoefs_[i][j]=active_instance->simulator_->objectEigencoefs()[i][j];
    //
    // for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
    //     for(int j=0;j<3;++j)
    //         active_instance->deformed_object_eigencoefs_[i][j]=active_instance->initial_object_eigencoefs_[i][j];
    // active_instance->example_guided_reduced_force_ = new double[3*active_instance->interpolate_eigenfunction_num_];

}

void OpenGLDriver::saveObjectEigenfunctions(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->saveObjectEigenfunctions(active_instance->output_object_eigen_file_name_))
    {
        std::cout<<"Error: load example eigenfunctions failed.\n";
        exit(0);
    }
}

void OpenGLDriver::saveObjectEigenfunctionsColor(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->saveObjectEigenfunctions(active_instance->object_eigen_file_name_))
    {
        std::cout<<"Error: load example eigenfunctions failed.\n";
        exit(0);
    }
    std::string output_file_name;
    std::stringstream adaptor;
    std::string eigen_idx_str;
    for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
    {
        adaptor.str("");
        adaptor.clear();
        adaptor<<i+1;
        adaptor>>eigen_idx_str;
        output_file_name=eigen_idx_str+".txt";//"new-armadillo-left.txt";
        active_instance->render_volumetric_mesh_->SaveVertexColorMap(active_instance->simulation_mesh_,
                        active_instance->simulator_->objectEigenFunctions()[i],
                        "armadillo_selected.bou",output_file_name.c_str());
    }
}

void OpenGLDriver::loadExampleEigenfunctions(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadExampleEigenFunctions(active_instance->example_eigen_file_name_prefix_))
    {
        std::cout<<"Error: load example eigenfunctions failed.\n";
        exit(0);
    }

    active_instance->glui_current_example_eigenfunctions_loaded_->set_name("Eigenfunctions for current example loaded:Yes");
    if(active_instance->render_mesh_type_==EXAMPLE_MESH)
        active_instance->render_eigen_index_spinner_->set_int_limits(1,active_instance->interpolate_eigenfunction_num_,GLUI_LIMIT_CLAMP);
    active_instance->isload_example_eigen_=true;

    //get example eigencoefs
    //  for(int i=0;i<active_instance->example_num_;++i)
    //      for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
    //      {
    //         active_instance->example_eigencoefs_[i][j]=active_instance->simulator_->exampleEigencoefs()[i][j];
    //      }
    std::cout<<"Load example eigenfunctions succeed!\n";
}

void OpenGLDriver::saveExampleEigenfunctions(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->saveExampleEigenfunctions(active_instance->output_eigen_file_name_prefix_))
    {
        std::cout<<"Error: failed to save example eigenfunctions.\n";
        exit(0);
    }

}

void OpenGLDriver::saveExampleEigenfunctionsColor(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadExampleEigenFunctions(active_instance->example_eigen_file_name_prefix_))
    {
        std::cout<<"Error: load example eigenfunctions failed.\n";
        exit(0);
    }

    std::string output_file_name;
    std::stringstream adaptor;
    std::string eigen_idx_str;
    std::stringstream adaptor1;
    std::string example_idx_str;
    for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
    {
        adaptor.str("");
        adaptor.clear();
        adaptor<<i+1;
        adaptor>>eigen_idx_str;
        for(int j=0;j<active_instance->example_num_;++j)
        {
            adaptor1.str("");
            adaptor1.clear();
            adaptor1<<j;
            adaptor1>>example_idx_str;
            output_file_name="ex"+example_idx_str+"-"+eigen_idx_str+".txt";
            active_instance->render_volumetric_mesh_->SaveVertexColorMap(active_instance->example_mesh_[j],
                            active_instance->simulator_->exampleEigenFunctions()[j][i],
                            "armadillo_selected.bou",output_file_name.c_str());
        }
    }
}
void OpenGLDriver::saveCurrentObjmesh(int code)
{
    // PerformanceCounter each_frame_counter;
    // each_frame_counter.StartCounter();
    OpenGLDriver *active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    // if(((active_instance->time_step_counter_+1)%(int)(1.0/(active_instance->frame_rate_*active_instance->time_step_))==0)&&(active_instance->time_step_counter_>0))
    // {
        std::stringstream adaptor;
        std::string output_object_file_name,output_file_index_str;
        std::string output_file_name_base(active_instance->output_objmesh_file_name_base_);

        cout<<"frame "<<(active_instance->time_step_counter_+1)/(int)(1.0/(active_instance->frame_rate_*active_instance->time_step_))<<" begins \n";
        adaptor<<active_instance->output_file_index_++;
        adaptor>>output_file_index_str;
        output_object_file_name=output_file_name_base+output_file_index_str+".obj";
        if(active_instance->simulation_mode_==FULLSPACE)
        {
            ObjMesh *mesh=active_instance->visual_mesh_->GetMesh();
            if(active_instance->enable_save_surfacemtl_)
                mesh->save(output_object_file_name,1);
            else
                mesh->save(output_object_file_name,0);
            //  delete mesh;
        }
        else
        {
            // ObjMesh *mesh=active_instance->render_reduced_surface_mesh_->GetMesh();
            ObjMesh *mesh=active_instance->visual_mesh_->GetMesh();
            if(active_instance->enable_save_surfacemtl_)
                mesh->save(output_object_file_name,1);
            else
                mesh->save(output_object_file_name,0);
            // delete mesh;
        }

        std::cout<<output_object_file_name<<" saved.\n";
    // }
    // each_frame_counter.StopCounter();
    // active_instance->total_write_file_time_+=each_frame_counter.GetElapsedTime();
    // if(active_instance->time_step_counter_==(active_instance->total_steps_-1))
    // {
    //     std::cout<<"Total write file time: "<<active_instance->total_write_file_time_<<" s.\n";
    // }
}

void OpenGLDriver::saveCurrentTetmesh(int code)
{
    //std::cout<<"aaaaaaaaaaaaaaaa\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    std::stringstream adaptor;
    std::string output_tet_file_name,output_file_index_str;
    std::string output_file_name_base(active_instance->output_objmesh_file_name_base_);
    cout<<"frame "<<(active_instance->time_step_counter_+1)/(int)(1.0/(active_instance->frame_rate_*active_instance->time_step_))<<" begins \n";
    adaptor<<active_instance->output_tet_file_index_++;
    adaptor>>output_file_index_str;
    output_tet_file_name=output_file_name_base+output_file_index_str+".smesh";
    std::ofstream output_file(output_tet_file_name.c_str());
	if(!output_file)
	{
		std::cout<<"Error: failed to open "<<output_tet_file_name<<".\n";
		return;
	}
	output_file<<"*VERTICES"<<std::endl;
	output_file<<active_instance->simulation_mesh_->getNumVertices()<<" 3 0 0"<<std::endl;
	for(unsigned int i=0;i<active_instance->simulation_mesh_->getNumVertices();++i)
	{
        Vec3d new_dis;
        new_dis[0]=(*active_instance->simulation_mesh_->getVertex(i))[0]+active_instance->u_[3*i];
        new_dis[1]=(*active_instance->simulation_mesh_->getVertex(i))[1]+active_instance->u_[3*i+1];
        new_dis[2]=(*active_instance->simulation_mesh_->getVertex(i))[2]+active_instance->u_[3*i+2];
		output_file<<i+1<<" "<<new_dis[0]<<" ";
		output_file<<new_dis[1]<<" ";
		output_file<<new_dis[2]<<std::endl;
	}
	output_file<<"*ELEMENTS"<<std::endl;
	output_file<<"TET"<<std::endl;
	output_file<<active_instance->simulation_mesh_->getNumElements()<<" 4 0"<<std::endl;
	for(unsigned int i=0;i<active_instance->simulation_mesh_->getNumElements();++i)
	{
		output_file<<i+1;
		for(unsigned int j=0;j<active_instance->simulation_mesh_->getNumElementVertices();++j)
		{
			unsigned int global_idx=active_instance->simulation_mesh_->getVertexIndex(i,j);
			output_file<<" "<<global_idx+1;
		}
		output_file<<std::endl;
	}
	output_file.close();
}
void OpenGLDriver::loadReducedBasis(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadReducedBasis(active_instance->reduced_basis_file_name_))
    {
        std::cout<<"Error: load reduced basis failed.\n";
        return;
    }
    std::cout<<"Load reduced basis succeed!\n";
}

void OpenGLDriver::loadObjectCubicaData(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadObjectCubicaData(active_instance->object_cubica_file_name_))
    {
        std::cout<<"Error: load cubica data for reduced simulation failed.\n";
        return;
    }
}
void OpenGLDriver::loadLBObjectCubicaData(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadObjectLBCubicaData(active_instance->object_LB_cubica_file_name_))
    {
        std::cout<<"Error: load cubica data for LB subspace failed.\n";
        return;
    }
}
void OpenGLDriver::loadCorrespondenceData(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadCorrespondenceData(active_instance->corresponding_file_name_))
    {
        std::cout<<"Error: load correspondence data failed.\n";
        return;
    }
    active_instance->isload_correspondence_data_=true;
}

void OpenGLDriver::registerEigenfunctions(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->isload_correspondence_data_)
    {
        std::cout<<"Error: the corresponding data is unloaded.\n";
        return;
    }
    if(!active_instance->simulator_->registerEigenfunctions())
    {
        std::cout<<"Error: register eigenfunctions failed.\n";
        return;
    }
}

void OpenGLDriver::resetDeformation(int code)
{
    // OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    // assert(active_instance);
    // active_instance->integrator_base_->ResetToRest();
    // active_instance->integrator_base_->SetState(active_instance->u_initial_,active_instance->vel_initial_);
    // //set the displacement of volumetric surface
    // memcpy(active_instance->u_,active_instance->integrator_base_->Getq(),sizeof(double)*3*active_instance->simulation_vertices_num_);
    //
    // active_instance->visual_mesh_->SetVertexDeformations(active_instance->u_);
    // //interpolate deformation from volumetric mesh to rendering triangel mesh
    // VolumetricMesh::interpolate(active_instance->u_,active_instance->u_render_surface_,active_instance->visual_mesh_->Getn(),
    //                             active_instance->object_interpolation_element_vertices_num_,
    //                             active_instance->object_interpolation_vertices_,active_instance->object_interpolation_weights_);
    // active_instance->visual_mesh_->SetVertexDeformations(active_instance->u_render_surface_);
    // active_instance->time_step_counter_=0;
    // //active_instance->output_file_index_=1;
    //
    // //this stuff may have been changed in projectOnExampleManifold, reset them
    // active_instance->last_initial_weight_=1.5;
    // active_instance->integrator_base_sparse_->SetDampingMassCoef(active_instance->damping_mass_coef_);
    // cout<<"Deformation reset completed.\n";
}
void OpenGLDriver::exitApplication(int code)
{
    //TO DO: release memories
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    // assert(active_instance);
    // if(active_instance->simulator_)
    //     delete active_instance->simulator_;
    // if(active_instance->camera_)
    //     delete active_instance->camera_;
    // if(active_instance->lighting_)
    //     delete active_instance->lighting_;
    // if(active_instance->simulation_mesh_)
    //     delete active_instance->simulation_mesh_;
    // // if(active_instance->tet_mesh_)
    // //     delete active_instance->tet_mesh_;
    // if(active_instance->visual_mesh_)
    //     delete active_instance->visual_mesh_;
    // for(int i=0;i<active_instance->example_num_;++i)
    // {
    //     if(active_instance->example_mesh_[i])
    //         delete active_instance->example_mesh_[i];
    //     if(active_instance->example_eigencoefs_[i])
    //         delete active_instance->example_eigencoefs_[i];
    // }
    // if(active_instance->initial_object_configurations_)
    //     delete[] active_instance->initial_object_configurations_;
    // //delete[] deformed_object_configurations_;
    // if(active_instance->temp_deformed_object_dis_)
    //     delete[] active_instance->temp_deformed_object_dis_;
    // if(active_instance->example_guided_reduced_force_)
    //     delete[] active_instance->example_guided_reduced_force_;
    // if(active_instance->current_example_mesh_)
    //     delete active_instance->current_example_mesh_;
    // if(active_instance->render_volumetric_mesh_)
    //     delete active_instance->render_volumetric_mesh_;
    // if(active_instance->visual_mesh_)
    //     delete active_instance->visual_mesh_;
    // if(active_instance->u_render_surface_)//
    //     delete active_instance->u_render_surface_;
    // if(active_instance->mesh_graph_)
    //     delete active_instance->mesh_graph_;
    // if(active_instance->mass_matrix_)
    //     delete active_instance->mass_matrix_;
    // if(active_instance->laplacian_matrix_)//
    //     delete active_instance->laplacian_matrix_;
    // if(active_instance->laplacian_damping_matrix_)
    //     delete active_instance->laplacian_damping_matrix_;
    // if(active_instance->isotropic_material_)
    //     delete active_instance->isotropic_material_;
    // if(active_instance->example_isotropic_material_)
    //     delete active_instance->example_isotropic_material_;
    // if(active_instance->isotropic_hyperelastic_fem_)
    //     delete active_instance->isotropic_hyperelastic_fem_;
    // if(active_instance->force_model_)
    //     delete active_instance->force_model_;
    // if(active_instance->integrator_base_)
    //     delete active_instance->integrator_base_;
    // if(active_instance->integrator_base_sparse_)
    //     delete active_instance->integrator_base_sparse_;
    //
    // if(active_instance->planes_) delete active_instance->planes_;
    // if(active_instance->u_)
    //     delete active_instance->u_;
    // if(active_instance->vel_)
    //     delete active_instance->vel_;
    // if(active_instance->u_initial_)
    //     delete active_instance->u_initial_;
    // if(active_instance->vel_initial_)
    //     delete active_instance->vel_initial_;
    // if(active_instance->f_ext_)
    //     delete active_instance->f_ext_;
    // if(active_instance->f_col_)
    //     delete active_instance->f_col_;
    // if(active_instance->fixed_vertices_)
    //     delete active_instance->fixed_vertices_;
    // if(active_instance->fixed_dofs_)
    //     delete active_instance->fixed_dofs_;
    // if(active_instance->force_loads_)
    //     delete active_instance->force_loads_;
    // if(active_instance->object_interpolation_vertices_)
    //     delete active_instance->object_interpolation_vertices_;
    // if(active_instance->object_interpolation_weights_)
    //     delete active_instance->object_interpolation_weights_;
    // for(int i=0;i<active_instance->extra_objects_num_;++i)
    //     if(active_instance->extra_objects_[i])
    //         delete active_instance->extra_objects_[i];
    // //reduced Simulation
    //
    // if(active_instance->modal_matrix_)
    //     delete[] active_instance->modal_matrix_;
    // if(active_instance->render_reduced_surface_mesh_)
    //     delete[] active_instance->render_reduced_surface_mesh_;
    // if(active_instance->reduced_force_model_)
    //     delete[] active_instance->reduced_force_model_;
    // if(active_instance->reduced_neoHookean_force_model_)
    //     delete[] active_instance->reduced_neoHookean_force_model_;
    // if(active_instance->reduced_mass_matrix_)
    //     delete[] active_instance->reduced_mass_matrix_;
    // if(active_instance->reduced_stiffness_matrix_)
    //     delete[] active_instance->reduced_stiffness_matrix_;
    // if(active_instance->reduced_f_ext_)
    //     delete[] active_instance->reduced_f_ext_;
    // if(active_instance->q_)
    //     delete[] active_instance->q_;
    // if(active_instance->fq_)
    //     delete[] active_instance->fq_;
    // if(active_instance->fqBase_)
    //     delete[] active_instance->fqBase_;
    // if(active_instance->U_)
    //     delete[] active_instance->U_;

    exit(0);
}

void OpenGLDriver::drawAxis(double axis_length) const
{
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glBegin(GL_LINES);
    for(unsigned int i = 0; i < 3; ++i)
    {
        float color[3] = {0,0,0};
        color[i] = 1.0f;
        glColor3fv(color);
        float vertex[3] = {0,0,0};
        vertex[i] = axis_length;
        glVertex3fv(vertex);
        glVertex3f(0,0,0);
    }
    // glLineWidth(1.0);
    // drawAxes(10.0);
    glEnd();
}
//draw color table
void OpenGLDriver::drawIndexColorTable() const
{
    glPushMatrix();
    //draw color table, size 20*128
    GLfloat refband[128][20][3];
    for(int idx = 0; idx < 128; idx += 2)
    {//color grades
        for(int kk = 0; kk < 3; kk++)
        {//rgb
            for(int jj = 0; jj < 20; jj++)
            {//width
                refband[idx][jj][kk] = ColorTable::Color_Table[idx>>1][kk];
                refband[idx+1][jj][kk] = ColorTable::Color_Table[idx>>1][kk];
            }
        }
    }
    glDrawPixels(20, 128, GL_RGB, GL_FLOAT, refband);
    glPopMatrix();
}

}  //namespace RTLB
