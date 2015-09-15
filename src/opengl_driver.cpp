/*
 * @file: opengl_driver.cpp
 * @brief: UI interface for simulation
 * @author: Fei Zhu
 *
 */
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <sstream>
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
    std::cout<<"0";
    if(active_instance_)
        delete active_instance_;
    active_instance_ = this;
    //TO DO: init everything and enter mainloop
    this->simulator_=new RTLB::RealTimeExampleBasedDeformer();

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
    //OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
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
    config_file_.addOptionOptional("objectCubicaFilename",object_cubica_file_name_,"none");
    //solver and materials config
    config_file_.addOptionOptional("solver",solver_method_,"UNKNOWN");
    config_file_.addOptionOptional("deformableModel",deformable_model_,"StVK");
    config_file_.addOptionOptional("invertibleMaterial",invertible_material_,"none");
    config_file_.addOptionOptional("principalStretchThreshold",&principal_stretch_threshold_,principal_stretch_threshold_);

    //simulation object
    config_file_.addOptionOptional("simulationMeshFilename",simulation_mesh_file_name_,"none");
    config_file_.addOptionOptional("objectRenderSurfaceMeshFilename",visual_mesh_file_name_,"none");
    config_file_.addOptionOptional("objectInterpolationFilename",object_interpolation_file_name_,"none");

    //examples
    config_file_.addOptionOptional("exampleFilenameBase",example_file_name_prefix_,"none");
    //config_file_.addOptionOptional("exampleSurfaceMeshFilenameBase",example_surface_mesh_file_name_,"none");
    //config_file_.addOptionOptional("exampleInterpolationFilenameBase",example_interpolation_file_name_,"none");

    //eigen configuration
    config_file_.addOption("exampleNum",&example_num_);
    config_file_.addOption("exampleEigunfunctionNum",&example_eigenfunction_num_);
    config_file_.addOption("objectEigenfunctionNum",&object_eigenfunction_num_);
    config_file_.addOptionOptional("exampleEigenFunctionFilenameBase",example_eigen_file_name_,"none");
    config_file_.addOptionOptional("objectEigenFunctionFilename",object_eigen_file_name_,"none");
    config_file_.addOptionOptional("correspondingFunctionFilename",corresponding_file_name_,"none");
    //output
    //config_file_.addOptionOptional("outputFilenameBase",output_file_name_,"output");

    //Simulation
    config_file_.addOptionOptional("forceNeighborSize",&force_neighbor_size_,force_neighbor_size_);
    config_file_.addOptionOptional("enableEigenWeightControl",&enable_eigen_weight_control_,enable_eigen_weight_control_);
    config_file_.addOptionOptional("timestep",&time_step_,time_step_);
    config_file_.addOptionOptional("massMatrixFilename",mass_matrix_file_name_,"none");
    config_file_.addOptionOptional("dampingMassCoef",&damping_mass_coef_,damping_mass_coef_);
    config_file_.addOptionOptional("dampingStiffnessCoef",&damping_stiffness_coef_,damping_stiffness_coef_);
    config_file_.addOptionOptional("dampingLaplacianCoef",&damping_laplacian_coef_,damping_laplacian_coef_);
    config_file_.addOptionOptional("newmarkBeta",&newmark_beta_,newmark_beta_);
    config_file_.addOptionOptional("newmarkGamma",&newmark_gamma_,newmark_gamma_);
    config_file_.addOptionOptional("integratorEpsilon",&integrator_epsilon_,integrator_epsilon_);
    config_file_.addOptionOptional("deformableObjectCompliance",&deformable_object_compliance_,deformable_object_compliance_);
    config_file_.addOptionOptional("exampleStiffnessScale",&example_stiffness_scale_,example_stiffness_scale_);
    config_file_.addOptionOptional("maxIterations",&max_iterations_,max_iterations_);
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
    config_file_.addOptionOptional("exampleAffectedVerticesFilenameBase",example_affected_vertices_file_name_,"none");

    //parse the configuration file
    if(config_file_.parseOptions((char*)config_file_name.c_str())!=0)
    {
        std::cout<<"Error: parsing options failed.\n";
        exit(0);
    }
    //print the variables that were just parsed
    config_file_.printOptions();
    //set the solver based on config file input
    if(strcmp(solver_method_,"implicitNewmark")==0)
        solver_type_=IMPLICITNEWMARK;
    if(strcmp(solver_method_,"implicitBackwardEuler")==0)
        solver_type_=IMPLICITBACKWARDEULER;
    if(strcmp(solver_method_,"Euler")==0)
        solver_type_=EULER;
    if(strcmp(solver_method_,"sympleticEuler")==0)
        solver_type_=SYMPLECTICEULER;
    if(strcmp(solver_method_,"centralDifferences")==0)
        solver_type_=CENTRALDIFFERENCES;
    if(solver_type_==UNKNOWN)
    {
        std::cout<<"Error:unknown implicit solver specified."<<std::endl;
        exit(0);
    }
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
    if(example_num_ > 0)
    {
        adaptor.clear();
        adaptor<<example_eigenfunction_num_;
        adaptor>>eigenfunction_num_str;
        static_text_content.clear();
        static_text_content = std::string("Number of eigenfunctions for examples: ");
        static_text_content += eigenfunction_num_str;
        glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
        adaptor.clear();
        glui_current_example_eigenfunctions_loaded_=glui_->add_statictext_to_panel(eigen_panel,"Eigenfunctions for current example loaded:No");
        glui_->add_button_to_panel(eigen_panel,"Load Eigenfunctions for All Examples",0,loadExampleEigenfunctions);
        glui_->add_separator_to_panel(eigen_panel);
    }
    adaptor.clear();

    std::string object_eigenfunction_num_str;
    adaptor<<object_eigenfunction_num_;
    adaptor>>object_eigenfunction_num_str;
    static_text_content.clear();
    static_text_content="Number of eigenfunctions for simulation object:";
    static_text_content+=object_eigenfunction_num_str;
    glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    static_text_content.clear();
    static_text_content="Simulation object surface eigenfunctions loaded:No";
    glui_object_surface_eigenfunctions_loaded_=glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    glui_->add_button_to_panel(eigen_panel,"Load Eigenfunctions for Simulation Object Surface",1,loadObjectEigenfunctions);
    glui_->add_separator_to_panel(eigen_panel);
    static_text_content.clear();
    static_text_content="Rendering eigenfunctions enable:No";
    glui_rendering_eigenfunctions_enabled_=glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    render_eigen_index_spinner_=glui_->add_spinner_to_panel(eigen_panel,"Current eigenfunction index: ",GLUI_SPINNER_INT,&current_render_eigen_idx_,0,changeCurrentEigenIndex);
    render_eigen_index_spinner_->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);

    //coupled quasi-harmonic bases
    if(example_num_>0)
    {
        std::cout<<simulator_->correspondingFunctionNum()<<std::endl;
        unsigned int coupled_eigenfunction_num=simulator_->correspondingFunctionNum();
        GLUI_Panel *coupled_panel=glui_->add_panel("Coupled quasi-harmonics",GLUI_PANEL_EMBOSSED);
        coupled_panel->set_alignment(GLUI_ALIGN_LEFT);
        GLUI_Spinner *coupled_base_num_spinner=glui_->add_spinner_to_panel(coupled_panel,"Coupled Eigenfunction Number:",GLUI_SPINNER_INT,&coupled_eigenfunction_num);
        coupled_base_num_spinner->set_int_limits(1,(object_eigenfunction_num_<example_eigenfunction_num_?object_eigenfunction_num_:example_eigenfunction_num_));
        glui_->add_button_to_panel(coupled_panel,"Load All Corresponding Functions",0,loadCorrespondenceData);
        glui_->add_button_to_panel(coupled_panel,"register Eigenfunctions",0,registerEigenfunctions);
    }
    //save eigenfunctions
    GLUI_Panel *save_panel=glui_->add_panel("Save Eigenfunctions",GLUI_PANEL_EMBOSSED);

    if(example_num_>0)
    {
        GLUI_Spinner *saved_base_num_spinner=glui_->add_spinner_to_panel(save_panel,"Saved Eigenfunction Number: ",GLUI_SPINNER_INT,&save_eigenfunction_num_);
        saved_base_num_spinner->set_int_limits(1,(object_eigenfunction_num_<example_eigenfunction_num_?object_eigenfunction_num_:example_eigenfunction_num_));
        //glui_->add_button_to_panel(save_panel,"Save Eigenfunctions for Current Example",0,/*saveCurrentEigenfunctions*/saveExampleEigenfunctions);
        glui_->add_button_to_panel(save_panel,"Save Eigenfunctions for All Examples",0,saveExampleEigenfunctions);
    }
    glui_->add_button_to_panel(save_panel,"Save Eigenfunctions for Simulation Object Surface",1,saveObjectEigenfunctions);

    //simulation
    GLUI_Panel *sim_panel=glui_->add_panel("Simulation",GLUI_PANEL_EMBOSSED);
    //add_rollout("Simulation",true);
    sim_panel->set_alignment(GLUI_ALIGN_LEFT);
    glui_->add_button_to_panel(sim_panel,"Reset",0,/*resetDeformation*/changeSimulationMode);
    //glui_->add_button_to_panel(sim_panel,"Save Current object Surface",0,/*saveCurrentObjectSurface*/changeSimulationMode);
    if(example_num_>0)
        change_simulation_mode_button_=glui_->add_button_to_panel(sim_panel,"Enable Example-based Simulation",0,changeSimulationMode);


    //exit button
    glui_->add_button("Exit",0,exitApplication);
    glui_->sync_live();
    glui_->set_main_gfx_window(window_id_);
}

void OpenGLDriver::initSimulation()
{
    //TO DO: setup simulation
    std::cout<<lighting_config_file_name_<<"\n";
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
        }
        else
        {
            std::cout<<"Error: configuration file for planes in scene not specified.\n";
            exit(1);
        }
    }
    //set deformable material type
    if(strcmp(simulation_mesh_file_name_,"none")!=0)
    {
        if(strcmp(deformable_model_,"StVK")==0)
            deformable_object_type_=STVK;
        if(strcmp(deformable_model_,"CLFEM")==0)
            deformable_object_type_=COROTLINFEM;
        if(strcmp(deformable_model_,"LinearFEM")==0)
            deformable_object_type_=LINFEM;
        if(strcmp(deformable_model_,"InvertibleFEM")==0)
            deformable_object_type_=INVERTIBLEFEM;
    }
    std::cout<<deformable_model_<<":"<<deformable_object_type_<<"\n";
    if(deformable_object_type_==UNSPECIFIED)
    {
        std::cout<<"Error: no deformable model specified.\n";
        exit(1);
    }
    std::cout<<deformable_object_type_<<"\n";
    //load mesh
    if((deformable_object_type_==STVK)||(deformable_object_type_==COROTLINFEM)||(deformable_object_type_==LINFEM)||(deformable_object_type_==INVERTIBLEFEM))
    {
        std::cout<<"Loading object volumetric mesh from file "<<simulation_mesh_file_name_<<".\n";
        simulation_mesh_=VolumetricMeshLoader::load(simulation_mesh_file_name_);
        if(simulation_mesh_==NULL)
        {
            std::cout<<"Error: unable to load the object simulation mesh from "<<simulation_mesh_file_name_<<".\n";
            exit(1);
        }
        simulation_vertices_num_=simulation_mesh_->getNumVertices();
        std::cout<<"simulation mesh vertices number: "<<simulation_vertices_num_<<".\n";
        std::cout<<"simulation mesh elements number: "<<simulation_mesh_->getNumElements()<<".\n";
        mesh_graph_=GenerateMeshGraph::Generate(simulation_mesh_);
        object_volume_=simulation_mesh_->getVolume();
        object_vertex_volume_=new double[simulation_mesh_->getNumVertices()];
        for(unsigned int ele_idx=0;ele_idx<simulation_mesh_->getNumElements();++ele_idx)
        {
            for(unsigned int i=0;i<simulation_mesh_->getNumElementVertices();++i)
            {
                unsigned int global_idx=simulation_mesh_->getVertexIndex(ele_idx,i);
                object_vertex_volume_[global_idx]+=simulation_mesh_->getElementVolume(ele_idx);
            }
        }
        //load mass matrix
        if(strcmp(mass_matrix_file_name_,"none")==0)
        {
            std::cout<<"Error: mass matrix for the deformable model not specified.\n";
            exit(1);
        }
        std::cout<<"Loading the mass matrix from file "<<mass_matrix_file_name_<<".\n";
        //get the mass matrix
        SparseMatrixOutline *mass_matrix_outline;
        try
        {
            //3 is expansion flag to indicate this is a mass matrix and does 3x3 identify block expansion
            mass_matrix_outline=new SparseMatrixOutline(mass_matrix_file_name_,3);
        }
        catch(int exception_code)
        {
            std::cout<<"Error: loading mass matrix failed.\n";
            exit(1);
        }
        mass_matrix_=new SparseMatrix(mass_matrix_outline);
        delete(mass_matrix_outline);
        if(deformable_object_type_==STVK||deformable_object_type_==LINFEM)//LINFEM constructed from StVKInternalforces
        {
            unsigned int loading_flag=0;
            precomputed_integrals_=StVKElementABCDLoader::load(simulation_mesh_,loading_flag);
            if(precomputed_integrals_==NULL)
            {
                std::cout<<"Error: unable to load the StVK integrals.\n";
                exit(1);
            }
            std::cout<<"Generating internal forces and stiffness matrix models.\n";
            stvk_internal_force_=new StVKInternalForces(simulation_mesh_,precomputed_integrals_,add_gravity_,gravity_);
            stvk_stiffness_matrix_=new StVKStiffnessMatrix(stvk_internal_force_);
        }
    }
    else
    {
        std::cout<<"Error: unsupported material type.\n";
        exit(1);
    }
    int scale_rows=1;
    mesh_graph_->GetLaplacian(&laplacian_matrix_,scale_rows);
    mesh_graph_->GetLaplacian(&laplacian_damping_matrix_,scale_rows);
    laplacian_damping_matrix_->ScalarMultiply(damping_laplacian_coef_);

    //init surface mesh of the volumetric Mesh
    if(strcmp(visual_mesh_file_name_,"none")==0)
    {
        std::cout<"Error: the surface mesh of the volumetric mesh was not specified.\n";
        exit(1);
    }
    visual_mesh_=new SceneObjectDeformable(visual_mesh_file_name_);
    //init object render surface Mesh
    isrender_volumetric_mesh_=false;
    isrender_surface_mesh_=true;
    render_surface_mesh_=visual_mesh_;
    if(enable_textures_)
        render_surface_mesh_->SetUpTextures(SceneObject::MODULATE,SceneObject::NOMIPMAP);
    render_surface_mesh_->ResetDeformationToRest();
    render_surface_mesh_->BuildNeighboringStructure();
    render_surface_mesh_->BuildNormals();
    u_render_surface_=new double[3*visual_mesh_->Getn()];//---------------to check
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
    // //load the fixed vertices
    // //1-indexed notation
    if(strcmp(fixed_vertex_file_name_,"none")==0)
    {
        fixed_vertices_num_=0;
        fixed_vertices_=NULL;
    }
    else
    {
        if(LoadList::load(fixed_vertex_file_name_,&fixed_vertices_num_,&fixed_vertices_)!=0)
        {
            std::cout<<"Error: reading fixed vertices failed.\n";
            exit(1);
        }
        LoadList::sort(fixed_vertices_num_,fixed_vertices_);
        std::cout<<"Loaded "<<fixed_vertices_num_<<" fixed vertices. They are: \n";
        LoadList::print(fixed_vertices_num_,fixed_vertices_);
        //creat 0-indexed fixed DOFs
        fixed_dofs_num_=3*fixed_vertices_num_;
        fixed_dofs_=new int [fixed_dofs_num_];
        for(unsigned int i=0;i<fixed_dofs_num_;++i)
        {
            fixed_dofs_[3*i+0]=3*fixed_vertices_[i]-3;
            fixed_dofs_[3*i+1]=3*fixed_vertices_[i]-2;
            fixed_dofs_[3*i+2]=3*fixed_vertices_[i]-1;
        }
        for(unsigned int i=0;i<fixed_vertices_num_;++i)
            fixed_vertices_[i]--;
        std::cout<<"Boundary vertices processed.\n";
    }

    //allocate space for deformation and force vectors
    u_=new double[3*simulation_vertices_num_];
    vel_=new double[3*simulation_vertices_num_];
    f_ext_=new double[3*simulation_vertices_num_];
    f_col_=new double[3*simulation_vertices_num_];

    //load initial position
    if(strcmp(initial_position_file_name_,"none")!=0)
    {
        int m,n;
        int status=ReadMatrixFromDiskTextFile(initial_position_file_name_,&m,&n,&u_initial_);
        if(status)
            exit(1);
        if((m!=3*simulation_vertices_num_)||(n!=1))
        {
            std::cout<<"Error: initial position matrix size mismatch.\n";
            exit(1);
        }
    }
    else
        u_initial_=new double[3*simulation_vertices_num_];

    //load initial velocity
    if(strcmp(initial_velocity_file_name_,"none")!=0)
    {
        int m,n;
        int status=ReadMatrixFromDiskTextFile(initial_velocity_file_name_,&m,&n,&vel_initial_);
        if(status)
            exit(1);
        if((m!=3*simulation_vertices_num_)||(n!=1))
        {
            std::cout<<"Error: initial velocity matrix size mismatch.\n";
            exit(1);
        }
    }
    else
        vel_initial_=new double[3*simulation_vertices_num_];

    //load force load
    if(strcmp(force_loads_file_name_,"none")!=0)
    {
        int m;
        int status=ReadMatrixFromDiskTextFile(force_loads_file_name_,&m,&force_loads_num_,&force_loads_);
        if(status)
            exit(1);
        if(m!=3*simulation_vertices_num_)
        {
            std::cout<<"Error: mismatch in the dimension of the force load matrix.\n";
            exit(1);
        }
    }
    else
        u_initial_=new double[3*simulation_vertices_num_];

    //create force models, to be used by the integrator
    std::cout<<"Creating force model:\n";
    if(deformable_object_type_==STVK)
        force_model_=new StVKForceModel(stvk_internal_force_,stvk_stiffness_matrix_);
    if(deformable_object_type_==COROTLINFEM)
    {
        tet_mesh_=dynamic_cast<TetMesh*>(simulation_mesh_);
        if(tet_mesh_=NULL)
        {
            std::cout<<"Error: the input mesh is not a tet mesh (CLFEM deformable model).\n";
            exit(1);
        }
        corotational_linear_fem_=new CorotationalLinearFEM(tet_mesh_);
        force_model_=new CorotationalLinearFEMForceModel(corotational_linear_fem_,corotational_linearFEM_warp_);
    }
    if(deformable_object_type_==LINFEM)
        force_model_=new LinearFEMForceModel(stvk_internal_force_);
    if(deformable_object_type_==INVERTIBLEFEM)
    {
        tet_mesh_=dynamic_cast<TetMesh*>(simulation_mesh_);
        if(tet_mesh_==NULL)
        {
            std::cout<<"Error: the input mesh is not a tet mesh (CLFEM deformable model).\n";
            exit(1);
        }
        //create the invertible material model
        if(strcmp(invertible_material_,"StVK")==0)
            invertible_material_type_=INV_STVK;
        if(strcmp(invertible_material_,"neoHookean")==0)
            invertible_material_type_=INV_NEOHOOKEAN;
        switch(invertible_material_type_)
        {
            case INV_STVK:
                isotropic_material_=new StVKIsotropicMaterial(tet_mesh_);
                std::cout<<"Invertible material: StVK.\n";
                break;
            case INV_NEOHOOKEAN:
                isotropic_material_=new NeoHookeanIsotropicMaterial(tet_mesh_);
                std::cout<<"Invertible material: neo-Hookean.\n";
                break;
            default:
                std::cout<<"Error: invalid invertible material type.\n";
                exit(1);
                break;
        }
        isotropic_hyperelastic_fem_=new IsotropicHyperelasticFEM(tet_mesh_,isotropic_material_,principal_stretch_threshold_,add_gravity_,gravity_);
        force_model_=new IsotropicHyperelasticFEMForceModel(isotropic_hyperelastic_fem_);
    }
    //initialize the integrator
    std::cout<<"Initializing the integrator, n= "<<simulation_vertices_num_<<".\n";
    std::cout<<"Solver type: "<<solver_type_<<".\n";
    integrator_base_sparse_=NULL;
    if(solver_type_==IMPLICITNEWMARK)
    {
        integrator_base_sparse_=new ImplicitNewmarkSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,
                                                        positive_definite_,fixed_dofs_num_,fixed_dofs_,
                                                        damping_mass_coef_,damping_stiffness_coef_,max_iterations_,
                                                        integrator_epsilon_,newmark_beta_,newmark_gamma_,solver_threads_num_);
    }
    else if(solver_type_==IMPLICITBACKWARDEULER)
    {
        integrator_base_sparse_=new ImplicitBackwardEulerSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,
                                                        positive_definite_,fixed_dofs_num_,fixed_dofs_,
                                                        damping_mass_coef_,damping_stiffness_coef_,max_iterations_,
                                                        integrator_epsilon_,solver_threads_num_);
    }
    else if(solver_type_==EULER)
    {
        int symplectic=0;
        integrator_base_sparse_=new EulerSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,symplectic,
                                                fixed_dofs_num_,fixed_dofs_,damping_mass_coef_);
    }
    else if(solver_type_==SYMPLECTICEULER)
    {
        int symplectic=1;
        integrator_base_sparse_=new EulerSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,symplectic,
                                                fixed_dofs_num_,fixed_dofs_,damping_mass_coef_);
    }
    else if(solver_type_==CENTRALDIFFERENCES)
    {
        integrator_base_sparse_=new CentralDifferencesSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,
                                                            fixed_dofs_num_,fixed_dofs_,damping_mass_coef_,damping_stiffness_coef_,
                                                            central_difference_tangential_damping_update_mode_,solver_threads_num_);
    }
    else
    {
    }

    integrator_base_=integrator_base_sparse_;
    if(integrator_base_==NULL)
    {
        std::cout<<"Error: failed to initialize numerical integrator.\n";
        exit(1);
    }

    //set integration parameters
    integrator_base_sparse_->SetDampingMatrix(laplacian_damping_matrix_);
    integrator_base_->ResetToRest();
    integrator_base_->SetState(u_initial_,vel_initial_);
    integrator_base_->SetTimestep(time_step_);

    //load example volumetric meshes
    if(example_num_>0)
    {
        example_mesh_=new VolumetricMesh *[example_num_];
    	example_mesh_=new VolumetricMesh *[example_num_];
    	example_volume_=new double[example_num_];
    	example_vertex_volume_=new double *[example_num_];
        std::stringstream stream;
        for(unsigned int i=0;i<example_num_;++i)
        {
            std::string example_idx_str;
            stream.str("");
            stream.clear();
            stream<<i;
            stream>>example_idx_str;
            std::string example_file_name=example_file_name_prefix_+example_idx_str+std::string(".veg");
            example_mesh_[i]=VolumetricMeshLoader::load(example_file_name.c_str());
            if(example_mesh_[i]==NULL)
            {
                std::cout<<"Error: unable to load the example mesh from "<<example_file_name<<std::endl;
                exit(1);
            }
            example_volume_[i]=example_mesh_[i]->getVolume();
    		std::cout<<i<<":"<<example_volume_[i]<<std::endl;
    		example_vertex_volume_[i]=new double[example_mesh_[i]->getNumVertices()];
    		for(unsigned int ele_idx=0;ele_idx<example_mesh_[i]->getNumElements();++ele_idx)
    		{
    			for(unsigned int j=0;j<example_mesh_[i]->getNumElementVertices();++j)
    			{
    				unsigned int global_idx=example_mesh_[i]->getVertexIndex(ele_idx,j);
    				example_vertex_volume_[i][global_idx]+=example_mesh_[i]->getElementVolume(ele_idx);
    			}
    		}
        }
        current_example_index_=1;
        current_example_mesh_=example_mesh_[0];
    }
    //load extra objects
    if(extra_objects_num_>0)
    {
        extra_objects_=new SceneObjectDeformable*[extra_objects_num_];
        for(unsigned int i=1;i<extra_objects_num_;++i)
        {
            std::string current_extra_obj_mesh_name(extra_objects_file_name_base_);
            std::stringstream stream;
            stream<<i;
            std::string current_extra_obj_index_str;
            stream>>current_extra_obj_index_str;
            current_extra_obj_mesh_name+=current_extra_obj_index_str;
            current_extra_obj_mesh_name+=".obj";
            extra_objects_[i-1]=new SceneObjectDeformable(current_extra_obj_mesh_name.c_str());
        }
    }
    std::cout<<"init Simulation finish.\n";
}

void OpenGLDriver::initCamera()
{
    std::cout<<"initCamera function:\n";
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
            active_instance->render_surface_mesh_->SetLighting(active_instance->lighting_);
        glEnable(GL_LIGHTING);
        glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE);
        glStencilFunc(GL_ALWAYS,0,~(0u));
    }
    //render eigenfunction
    if(active_instance->render_eigenfunction_)
    {
        //active_instance->render_volumetric_mesh_->RenderVertexColorMap(active_instance->simulation_mesh_,);
        active_instance->drawIndexColorTable();//draw color table at left bottom corner of the window
    }

    //render extra objects
    if(active_instance->extra_objects_num_>0)//render the extra objects in sceneObjectDeformable
    {
        for(int i=0;i<active_instance->extra_objects_num_;++i)
            active_instance->extra_objects_[i]->Render();
    }
    //render different mode
    if(active_instance->isrender_surface_mesh_&&active_instance->render_mesh_type_==VISUAL_MESH)
    {//object surface mesh
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0,1.0);
        active_instance->render_surface_mesh_->Render();
        if(active_instance->render_vertices_)
        {
            glDisable(GL_LIGHTING);
            glColor3f(0.5,0.0,0.0);
            glPointSize(8.0);
            active_instance->render_surface_mesh_->RenderVertices();
            glEnable(GL_LIGHTING);
        }
        if(active_instance->render_wireframe_)
        {
            glDisable(GL_LIGHTING);
            glColor3f(0.0,0.0,0.0);
            active_instance->render_surface_mesh_->RenderVertices();
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
            glDisable(GL_LIGHTING);
            glColor3f(0.0,0.5,0.0);
            active_instance->render_volumetric_mesh_->Render(active_instance->simulation_mesh_);
            if(active_instance->render_vertices_)
            {
                glDisable(GL_LIGHTING);
                glColor3f(0.5,0.0,0.0);
                glPointSize(8.0);
                active_instance->render_volumetric_mesh_->RenderVertices(active_instance->simulation_mesh_);
                glEnable(GL_LIGHTING);
            }
            if(active_instance->render_wireframe_)
            {
                glDisable(GL_LIGHTING);
                glColor3f(0.0,0.0,0.0);
                active_instance->render_volumetric_mesh_->RenderWireframe(active_instance->simulation_mesh_);
                glEnable(GL_LIGHTING);
            }
            glDisable(GL_BLEND);
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
            glDisable(GL_LIGHTING);
            glColor3f(0.0,0.5,0.0);
            active_instance->render_volumetric_mesh_->Render(active_instance->example_mesh_[active_instance->current_example_index_-1]);
            if(active_instance->render_vertices_)
            {
                glDisable(GL_LIGHTING);
                glColor3f(0.5,0.0,0.0);
                glPointSize(8.0);
                active_instance->render_volumetric_mesh_->RenderVertices(active_instance->example_mesh_[active_instance->current_example_index_-1]);
                glEnable(GL_LIGHTING);
            }
            if(active_instance->render_wireframe_)
            {
                glDisable(GL_LIGHTING);
                glColor3f(0.0,0.0,0.0);
                active_instance->render_volumetric_mesh_->RenderWireframe(active_instance->example_mesh_[active_instance->current_example_index_-1]);
                glEnable(GL_LIGHTING);
            }
            glDisable(GL_BLEND);
        }
        active_instance->isrender_surface_mesh_=false;
    }
    //render planes
    if(active_instance->plane_num_>0)
    {
        active_instance->planes_->render();
    }
    glDisable(GL_LIGHTING);
    glStencilFunc(GL_ALWAYS,1,~(0u));
    glColor3f(1.0,0.0,0.0);
    glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP);
    //render axis
    if(active_instance->render_axis_)
    {
        glLineWidth(1.0);
        active_instance->drawAxis(1.0);
    }
    //render the currently pulled vertex
    if(active_instance->pulled_vertex_!=-1)
    {
        glColor3f(0.0,1.0,0.0);
        double pulled_vertex_pos[3];
        active_instance->visual_mesh_->GetSingleVertexPositionFromBuffer(active_instance->pulled_vertex_,&pulled_vertex_pos[0],&pulled_vertex_pos[1],&pulled_vertex_pos[2]);
        glEnable(GL_POLYGON_OFFSET_POINT);
        glPolygonOffset(-1.0,-1.0);
        glBegin(GL_POINTS);
        glVertex3f(pulled_vertex_pos[0],pulled_vertex_pos[1],pulled_vertex_pos[2]);
        glEnd();
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
    //render fixed vertices
    if(active_instance->render_fixed_vertices_)
    {
        for(int i=0;i<active_instance->fixed_vertices_num_;++i)
        {
            glColor3f(1,0,0);
            double fixed_vertex_pos[3];
            active_instance->visual_mesh_->GetSingleVertexPositionFromBuffer(active_instance->fixed_vertices_[i],&fixed_vertex_pos[0],&fixed_vertex_pos[1],&fixed_vertex_pos[2]);
            glEnable(GL_POLYGON_OFFSET_POINT);
            glPolygonOffset(-1.0,-1.0);
            glPointSize(12.0);
            glBegin(GL_POINTS);
            glVertex3f(fixed_vertex_pos[0],fixed_vertex_pos[1],fixed_vertex_pos[2]);
            glEnd();
            glDisable(GL_POLYGON_OFFSET_FILL);
        }
    }
    //render example guided surface deformation---to do
    glutSwapBuffers();
}

void OpenGLDriver::idleFunction()
{
    //std::cout<<"idleFunction:\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    glutSetWindow(active_instance->window_id_);
    if(!active_instance->pause_simulation_)
    {
        if(active_instance->left_button_down_)
        {
            if(active_instance->pulled_vertex_!=-1)
            {
                double force_x=active_instance->mouse_pos_[0]-active_instance->drag_start_x_;
                double force_y=(-1.0)*(active_instance->mouse_pos_[1]-active_instance->drag_start_y_);
                double external_force[3];
                active_instance->camera_->CameraVector2WorldVector_OrientationOnly3D(force_x,force_y,0,external_force);
                for(int i=0;i<3;++i)
                {
                    external_force[i]*=active_instance->deformable_object_compliance_;
                }
                std::cout<<active_instance->pulled_vertex_<<" fx: "<<force_x<<",fy: "<<force_y<<" | "<<external_force[0]<<",";
                std::cout<<external_force[1]<<","<<external_force[2]<<std::endl;
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
        //apply any scripted force loads
        if(active_instance->time_step_counter_<active_instance->force_loads_num_)
        {
            std::cout<<"External forces read from the text input file.\n";
            for(int i=0;i<3*active_instance->simulation_vertices_num_;++i)
                active_instance->f_ext_[i]+=active_instance->force_loads_[ELT(3*active_instance->simulation_vertices_num_,i,active_instance->time_step_counter_)];
        }
        //apply the force loads caused by the examples---to do

        //apply the penalty collision forces with planes in scene
        if(active_instance->plane_num_>0)
        {
            active_instance->planes_->resolveContact(active_instance->visual_mesh_->GetMesh(),active_instance->f_col_);
            for(int i=0;i<active_instance->simulation_vertices_num_;++i)
                active_instance->f_ext_[i]+=active_instance->f_col_[i];
        }
        //set forces to the integrator
        active_instance->integrator_base_sparse_->SetExternalForces(active_instance->f_ext_);
        //time step the dynamics
        int code=active_instance->integrator_base_->DoTimestep();
        std::cout<<".";
        ++active_instance->time_step_counter_;
        memcpy(active_instance->u_,active_instance->integrator_base_->Getq(),sizeof(double)*3*active_instance->simulation_vertices_num_);

    }
     active_instance->visual_mesh_->SetVertexDeformations(active_instance->u_);//set the displacement of volumetric Surface
    //interpolate deformations from volumetric mesh to object surface mesh, update its configuration
    VolumetricMesh::interpolate(active_instance->u_,active_instance->u_render_surface_,active_instance->visual_mesh_->Getn(),
                                active_instance->object_interpolation_element_vertices_num_,active_instance->object_interpolation_vertices_,
                                active_instance->object_interpolation_weights_);
    // if(active_instance->render_mesh_type_==VISUAL_MESH)
    //     active_instance->render_surface_mesh_->SetVertexDeformations(active_instance->u_render_surface_);//not done yet
    // else if(active_instance->render_mesh_type_==OBJECT_EIGEN_MESH)
    // {
    // //    active_instance->render_volumetric_mesh_->RenderWireframeDeformation(active_instance->simulation_mesh_,active_instance_->u_);
    // }
    // else
    // {
    //     //active_instance->render_volumetric_mesh_->RenderWireframeDeformation(active_instance->example_mesh_[active_instance->current_example_index_-1],active_instance_->u_);
    //     //active_instance->render_volumetric_mesh_->Render(active_instance->simulation_mesh_);
    // }
    //save object surface mesh to files--not done yet
    glutPostRedisplay();
}

void OpenGLDriver::reshapeFunction(int width, int height)
{
    std::cout<<"reshapeFunction:\n";
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
    std::cout<<"keyboardFunction:\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    switch(key)
    {
    case 27: //ESC, exit
        exitApplication(0);
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
    case 'r': //reset camera
        active_instance->camera_->Reset();
        break;
    case 32: //space button, pause simulation
        active_instance->pause_simulation_ = !(active_instance->pause_simulation_);
        break;
    case 'v': //render vertices
        active_instance->render_vertices_ = !(active_instance->render_vertices_);
        break;
    case 'e': //render eigenfunctions
        active_instance->render_eigenfunction_ = !(active_instance->render_eigenfunction_);
        std::string static_text_content("Rendering eigenfunctions enabled: ");
        if(active_instance->render_eigenfunction_)
        {
            static_text_content+="Yes";
            active_instance->setEigenfunctionColors();
        }
        else
        {
            static_text_content="No";
        }
        active_instance->glui_rendering_eigenfunctions_enabled_->set_name(static_text_content.c_str());
        break;
    }
}

void OpenGLDriver::specialFunction(int key, int x, int y)
{
    std::cout<<"specialFunction:\n";
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
    std::cout<<"motionFunction:\n";
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
    std::cout<<"mouseFunction:\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    switch (button)
     {
        case GLUT_LEFT_BUTTON:
            active_instance->left_button_down_=(state==GLUT_DOWN);
            active_instance->shift_pressed_=(glutGetModifiers()==GLUT_ACTIVE_SHIFT);
            active_instance->alt_pressed_=(glutGetModifiers()==GLUT_ACTIVE_ALT);
            active_instance->ctrl_pressed_=(glutGetModifiers()==GLUT_ACTIVE_CTRL);
            if(active_instance->left_button_down_&&!active_instance->shift_pressed_&&active_instance->ctrl_pressed_) //used for pulled vertex, apply force
            {
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
                    active_instance->pulled_vertex_=active_instance->visual_mesh_->GetClosestVertex(pos);
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

void OpenGLDriver::updateRenderMesh(int code)
{
    std::cout<<"updateRenderMesh function:\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);

    std::cout<<"active_instance:"<<active_instance->render_mesh_type_<<"\n";
    if(active_instance->render_mesh_type_==VISUAL_MESH)
    {
        active_instance->render_surface_mesh_=active_instance->visual_mesh_;
        active_instance->render_surface_mesh_->Render();
        active_instance->isrender_volumetric_mesh_=false;
        active_instance->isrender_surface_mesh_=true;
        active_instance->render_eigen_index_spinner_->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);
     }
    else if(active_instance->render_mesh_type_==OBJECT_EIGEN_MESH)
    {
        if(active_instance->simulation_mesh_==NULL)
        {
            std::cout<<"Error: object simulation mesh is null.\n";
            exit(1);
        }
        else
        {
            //active_instance->render_volumetric_mesh_->Render(active_instance->simulation_mesh_);
            active_instance->render_volumetric_mesh_->RenderWireframe(active_instance->simulation_mesh_);
            if(active_instance->isrender_object_eigen_)
                active_instance->render_eigen_index_spinner_->set_int_limits(1,active_instance->object_eigenfunction_num_,GLUI_LIMIT_CLAMP);
            else
                active_instance->render_eigen_index_spinner_->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);
            active_instance->isrender_volumetric_mesh_=true;
            active_instance->isrender_surface_mesh_=false;
        }
    }
    else if(active_instance->render_mesh_type_==EXAMPLE_MESH)
    {
        if(active_instance->example_mesh_[active_instance->current_example_index_-1]==NULL)
        {
            std::cout<<"Error: example volumetric mesh is null.\n";
            exit(1);
        }
        else
        {
            active_instance->render_volumetric_mesh_->Render(active_instance->example_mesh_[active_instance->current_example_index_-1]);
            if(active_instance->isrender_example_eigen_)
                active_instance->render_eigen_index_spinner_->set_int_limits(1,active_instance->example_eigenfunction_num_,GLUI_LIMIT_CLAMP);
            else
                active_instance->render_eigen_index_spinner_->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);
            active_instance->isrender_volumetric_mesh_=true;
            active_instance->isrender_surface_mesh_=false;
        }
    }
    else
    {
    }
}
void OpenGLDriver::updateCurrentExample(int code)
{
    std::cout<<"updateCurrentExample function:\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    active_instance->current_example_mesh_=active_instance->example_mesh_[active_instance->current_example_index_-1];
    if(!code)
    {
        //eigenfunctions--to do later

    }
    updateRenderMesh(0);
}
void OpenGLDriver::changeCurrentEigenIndex(int code)
{
    OpenGLDriver* active_instance=OpenGLDriver::activeInstance();
    assert(active_instance);
    if(active_instance->render_eigenfunction_)
        active_instance->setEigenfunctionColors();
}
void OpenGLDriver::changeSimulationMode(int code)
{
    //TO DO
    std::cout<<"changeSimulationMode function:\n";
}

void OpenGLDriver::loadObjectEigenfunctions(int code)
{
    std::cout<<"loadObjectEigenfunctions:\n";
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadObjectEigenfunctions(active_instance->object_eigen_file_name_))
    {
        std::cout<<"Error: load object eigenfunctions failed.\n";
        return;
    }
    std::cout<<"Load object eigenfunctions succeed!\n";
    active_instance->isrender_object_eigen_=true;
    active_instance->glui_object_surface_eigenfunctions_loaded_->set_name("Simulation object surface eigenfunctions loaded:Yes");
    if(active_instance->render_mesh_type_==OBJECT_EIGEN_MESH)
        active_instance->render_eigen_index_spinner_->set_int_limits(1,active_instance->object_eigenfunction_num_,GLUI_LIMIT_CLAMP);
}

void OpenGLDriver::saveObjectEigenfunctions(int code)
{
    //TO DO
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->saveObjectEigenfunctions(active_instance->example_file_name_prefix_))
    {
        std::cout<<"Error: load example eigenfunctions failed.\n";
        return;
    }
}

void OpenGLDriver::loadExampleEigenfunctions(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadExampleEigenFunctions(active_instance->example_file_name_prefix_))
    {
        std::cout<<"Error: load example eigenfunctions failed.\n";
        return;
    }
    std::cout<<"Load example eigenfunctions succeed!\n";
    active_instance->glui_current_example_eigenfunctions_loaded_->set_name("Eigenfunctions for current example loaded:Yes");
    active_instance->isrender_example_eigen_=true;
    if(active_instance->render_mesh_type_==EXAMPLE_MESH)
        active_instance->render_eigen_index_spinner_->set_int_limits(1,active_instance->example_eigenfunction_num_,GLUI_LIMIT_CLAMP);
}

void OpenGLDriver::saveExampleEigenfunctions(int code)
{
    //TO DO
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->saveExampleEigenfunctions(active_instance->example_file_name_prefix_));
    {
        std::cout<<"Error: failed to save example eigenfunctions.\n";
        return;
    }

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
}

void OpenGLDriver::loadObjectCubicaData(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadObjectCubicaData(active_instance->object_cubica_file_name_))
    {
        std::cout<<"Error: load cubica data failed.\n";
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
}

void OpenGLDriver::registerEigenfunctions(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->registerEigenfunctions())
    {
        std::cout<<"Error: register eigenfunctions failed.\n";
        return;
    }
}

void OpenGLDriver::exitApplication(int code)
{
    //TO DO: release memories
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(active_instance->simulator_)
        delete active_instance->simulator_;
    if(active_instance->camera_)
        delete active_instance->camera_;
    if(active_instance->lighting_)
        delete active_instance->lighting_;
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
void OpenGLDriver::setEigenfunctionColors()
{
    //not done yet
    double *eigenfunction=NULL;
    int eigenfunction_dim=0;
    //isrender_volumetric_mesh_=false;
    if(render_mesh_type_==OBJECT_EIGEN_MESH)
    {
        if(!isrender_object_eigen_)
        {
            std::cout<<"Error: object eigenfuncion unloaded.\n";
            return;
        }
        eigenfunction=simulator_->objectEigenFunctions()[current_render_eigen_idx_-1];

        render_volumetric_mesh_->RenderVertexColorMap(simulation_mesh_,eigenfunction);
        
        //render_volumetric_mesh_->RenderSolidAndWireframe(simulation_mesh_);
        //eigenfunction_dim=simulation_mesh_->getNumVertices();
    }
    else if(render_mesh_type_==EXAMPLE_MESH)
    {
        if(!isrender_example_eigen_)
        {
            std::cout<<"Error: example eigenfunctions unloaded.\n";
            return;
        }
        eigenfunction=simulator_->exampleEigenFunctions()[current_example_index_-1][current_render_eigen_idx_-1];
        //eigenfunction_dim=example_mesh_[current_example_index_-1]->getNumVertices();
    }
    else
    {
        std::cout<<"Error: current surface is for rendering only. No eigenfunctions loaded.\n";
        return;
    }
    //std::cout<<eigenfunction_dim<<"--------------\n";
    //vector<Vec3d> custom_colors(eigenfunction_dim);
    // double min_value=eigenfunction[0],max_value=eigenfunction[0];
    // std::cout<<min_value<<","<<max_value<<"\n";
    // for(int i=1;i<eigenfunction_dim;++i)//plot colors relatively
    // {
    //     if(eigenfunction[i]>max_value)
    //         max_value=eigenfunction[i];
    //     if(eigenfunction[i]<min_value)
    //         min_value=eigenfunction[i];
    // }
    // std::cout<<min_value<<","<<max_value<<"\n";
    // for(int i=0;i<eigenfunction_dim;++i)
    // {
    //     double numerator=(eigenfunction[i]-min_value)*63;
    //     double denominator=max_value-min_value;
    //     int color_table_index=0;
    //     if(denominator>epsilon_)
    //     {
    //         while(denominator<1.0)//for numerical reasons, avoid extremely small denominator
    //         {
    //             numerator*=100.0;
    //             denominator*=100.0;
    //         }
    //         color_table_index=numerator/denominator;
    //     }
    //     Vec3d color(ColorTable::Color_Table[color_table_index]);
    //     custom_colors[i]=color;
    // }

    //isrender_volumetric_mesh_=false;
    // if(render_mesh_type_==OBJECT_EIGEN_MESH)
    // {
    //     std::cout<<"aaaaaaaaaa\n";
    //     render_volumetric_mesh_->RenderVertexColorMap(simulation_mesh_,eigenfunction,simulator_->objectEigenValues());
    //     std::cout<<"Bbbbb\n";
    // }
    // else if(render_mesh_type_==EXAMPLE_MESH)
    // {
    //     render_volumetric_mesh_->RenderVertexColorMap(example_mesh_[current_example_index_-1],eigenfunction,simulator_->exampleEigenValues()[current_example_index_-1]);
    // }
    // else
    // {}

}

}  //namespace RTLB
