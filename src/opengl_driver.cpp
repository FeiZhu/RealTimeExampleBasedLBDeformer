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
#include "GL/freeglut.h"
#include "GL/glui.h"
#include "camera.h"
#include "lighting.h"
#include "real_time_example_based_deformer.h"
#include "opengl_driver.h"


namespace RTLB{

OpenGLDriver* OpenGLDriver::active_instance_ = NULL;

OpenGLDriver::OpenGLDriver(const std::string &config_file_name)
{
    if(active_instance_)
        delete active_instance_;
    active_instance_ = this;
    //TO DO: init everything and enter mainloop
    this->simulator_=new RTLB::RealTimeExampleBasedDeformer();
    initGLUT();
    initGLUI();
    initGraphics();
    initConfigurations(config_file_name);
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

    //a=config_file.addOptionOptional("windowWidth",&window_width_,window_width_);
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
    config_file_.addOptionOptional("planesFilename",plane_file_name_,plane_file_name_);
    config_file_.addOptionOptional("planeNumber",&plane_num_,plane_num_);
    //cubica FILE
    config_file_.addOptionOptional("objectCubicaFilename",object_cubica_file_name_,"none");
    //solver and materials config
    config_file_.addOptionOptional("solver",solver_method_,"implicitBackwardEuler");
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
    config_file_.addOptionOptional("exampleEigenFunctionFilenameBase",example_eigen_file_name_,"eigen_example");
    config_file_.addOptionOptional("objectEigenFunctionFilename",object_eigen_file_name_,"eigen_object");
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
    if(solver_type_=UNKNOWN)
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


    //temp setting
    example_num_=2;
    //temp setting


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
        example_index_spinner->set_int_limits(0,example_num_-1,GLUI_LIMIT_CLAMP);
    }
    //eigenfunction information
    GLUI_Panel *eigen_panel = glui_->add_panel("Laplace-Beltrami Eigenfunctions",GLUI_PANEL_EMBOSSED);
    eigen_panel->set_alignment(GLUI_ALIGN_LEFT);
    std::string eigenfunction_num_str;
    if(example_num_ > 0)
    {
        adaptor.clear();
        example_eigenfunction_num_=simulator_->exampleEigenfunctionNum();
        adaptor<<example_eigenfunction_num_;
        adaptor>>eigenfunction_num_str;
        static_text_content.clear();
        static_text_content = std::string("Number of eigenfunctions for examples: ");
        static_text_content += eigenfunction_num_str;
        glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
        adaptor.clear();
        adaptor<<example_with_eigen_num_;
        std::string example_with_eigen_num_str;
        adaptor>>example_with_eigen_num_str;
        static_text_content.clear();
        static_text_content = std::string("Examples with eigenfunctions loaded: ");
        static_text_content += example_with_eigen_num_str;
        glui_example_with_eigenfunctions_loaded=glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
        glui_current_example_eigenfunctions_loaded=glui_->add_statictext_to_panel(eigen_panel,"Eigenfunctions for current example loaded:No");
        //glui_->add_button_to_panel(eigen_panel,"Load Eigenfunctions for Current Example",0,/*loadCurrentEigenFunctions*/loadExampleEigenfunctions);
        glui_->add_button_to_panel(eigen_panel,"Load Eigenfunctions for All Examples",0,loadExampleEigenfunctions);
        glui_->add_separator_to_panel(eigen_panel);
    }
    adaptor.clear();

    std::string object_eigenfunction_num_str;
    object_eigenfunction_num_=simulator_->objectEigenfunctionNum();
    adaptor<<object_eigenfunction_num_;
    adaptor>>object_eigenfunction_num_str;
    static_text_content.clear();
    static_text_content="Number of eigenfunctions for simulation object:";
    static_text_content+=object_eigenfunction_num_str;
    glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    static_text_content.clear();
    static_text_content="Simulation object surface eigenfunctions loaded:No";
    glui_object_surface_eigenfunctions_loaded=glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    glui_->add_button_to_panel(eigen_panel,"Load Eigenfunctions for SImulation Object Surface",1,/*loadCurrentEigenFunctions*/loadExampleEigenfunctions);
    glui_->add_separator_to_panel(eigen_panel);
    static_text_content.clear();
    static_text_content="Rendering eigenfunctions enable:No";
    glui_rendering_eigenfunctions_enabled=glui_->add_statictext_to_panel(eigen_panel,static_text_content.c_str());
    //render_eigen_index_spinner=glui_->add_spinner_to_panel(eigen_panel,"Current eigenfunction index: ",GLUI_SPINNER_INT,&current_render_eigenfunction,0,changeCurrentEigenIndex);
    //render_eigen_index_spinner->set_int_limits(1,object_eigenfunction_num_,GLUI_LIMIT_CLAMP);

    //coupled quasi-harmonic bases
    if(example_num_>0)
    {
        std::cout<<simulator_->correspondingFunctionNum()<<std::endl;
        unsigned int coupled_eigenfunction_num=simulator_->correspondingFunctionNum();
        GLUI_Panel *coupled_panel=glui_->add_panel("Coupled quasi-harmonics",GLUI_PANEL_EMBOSSED);
        //add_rollout_to_panel(eigen_panel,"Coupled quasi-harmonics",true);
        coupled_panel->set_alignment(GLUI_ALIGN_LEFT);
        GLUI_Spinner *coupled_base_num_spinner=glui_->add_spinner_to_panel(coupled_panel,"Coupled Eigenfunction Number:",GLUI_SPINNER_INT,&coupled_eigenfunction_num);
        coupled_base_num_spinner->set_int_limits(1,(object_eigenfunction_num_<example_eigenfunction_num_?object_eigenfunction_num_:example_eigenfunction_num_));
        glui_->add_button_to_panel(coupled_panel,"Load All Corresponding Functions",0,loadCorrespondenceData);
        glui_->add_button_to_panel(coupled_panel,"register Eigenfunctions",0,registerEigenfunctions);
    }
    //save eigenfunctions
    GLUI_Panel *save_panel=glui_->add_panel("Save Eigenfunctions",GLUI_PANEL_EMBOSSED);
    //glui_->add_rollout_to_panel(eigen_panel,"Save Eigenfunctions",true);
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
        change_simulation_mode_button=glui_->add_button_to_panel(sim_panel,"Enable Example-based Simulation",0,changeSimulationMode);


    //exit button
    glui_->add_button("Exit",0,exitApplication);
    glui_->sync_live();
    glui_->set_main_gfx_window(window_id_);
}

void OpenGLDriver::initSimulation()
{
    //TO DO: setup simulation
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
    if(deformable_object_type_==UNSPECIFIED)
    {
        std::cout<<"Error: no deformable model specified.\n";
        exit(1);
    }
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
        //mesh_graph_=GenerateMeshGraph::Generate(simulation_mesh_);

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
//    mesh_graph_->GetLaplacian(&laplacian_matrix_,scale_rows);
//    mesh_graph_->GetLaplacian(&laplacian_damping_matrix_,scale_rows);
    laplacian_damping_matrix_->ScalarMultiply(damping_laplacian_coef_);

    //init surface mesh of the volumetric Mesh
    if(strcmp(visual_mesh_file_name_,"none")==0)
    {
        std::cout<"Error: the surface mesh of the volumetric mesh was not specified.\n";
        exit(1);
    }
    visual_mesh_=new SceneObjectDeformable(visual_mesh_file_name_);

    //load example volumetric meshes
    if(example_num_>0)
    {
        example_mesh_=new VolumetricMesh *[example_num_];
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
        }
    }
    // //init eigenfunctions
    //
    // //load interpolation structure for object
    //
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
    }
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
    //
    // //allocate space for deformation and force vectors
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

    // //create force models, to be used by the integrator
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
        if(tet_mesh_=NULL)
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
    // //initialize the integrator
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
    //
    // //set integration parameters
    // integrator_base_sparse_->SetDampingMatrix(laplacian_damping_matrix_);
    // integrator_base_->ResetToRest();
    // integrator_base_->SetState(u_initial_,vel_initial_);
    // integrator_base_->setTimeStep(time_step_);

    //allocate storage for example
}

void OpenGLDriver::initCamera()
{
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
    // clear to white
    glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);
    // clear to light blue
    //glClearColor(233.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);
    // clear to gray
    //glClearColor(196.0 / 256, 196.0 / 256, 196.0 / 256, 0.0);
    // clear to brown
    //glClearColor(255.0 / 256, 156.0 / 256, 17.0 / 256, 0.0);
    // clear to medical blue
    //glClearColor(148.0 / 256, 199.0 / 256, 211.0 / 256, 0.0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_LINE_SMOOTH);

    reshapeFunction(window_width_,window_height_);
    //init Lighting
    try
    {
        lighting_ = new Lighting(lighting_config_file_name_);
    }
    catch(int exception_code)
    {
        std::cout<<"Error: "<<exception_code<<" reading lighting information from "<<lighting_config_file_name_<<".\n";
        exit(0);
    }

    //init camera
    initCamera();

    //printf ("Graphics initialization complete.\n")
}

void OpenGLDriver::displayFunction()
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    active_instance->camera_->Look();
    //to do
}

void OpenGLDriver::idleFunction()
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    glutSetWindow(active_instance->window_id_);
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

}

void OpenGLDriver::updateRenderMesh(int code)
{
    //TO DO
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(active_instance->render_mesh_type_==VISUAL_MESH)
    {
    //    active_instance->rendering_mesh_=visual_mesh_;
    //    active_instance->render_eigen_index_spinner->set_int_limits(1,active_instance->object_eigenfunction_num_,GLUI_LIMIT_CLAMP);
    }
    else if(active_instance->render_mesh_type_==OBJECT_EIGEN_MESH)
    {
    //    active_instance->rendering_mesh_=simulation_mesh_;
    //    active_instance->render_eigen_index_spinner->set_int_limits(1,active_instance->example_eigenfunction_num_,GLUI_LIMIT_CLAMP);

    }
    else
    {
    //    active_instance->rendering_mesh_=example_mesh_[current_example_index_];
    //    active_instance->render_eigen_index_spinner->set_int_limits(-1,-1,GLUI_LIMIT_CLAMP);

    }

}

void OpenGLDriver::updateCurrentExample(int code)
{
    //TO DO
    /*if(example_mesh_[current_example_index_-1]==NULL)
    {

    }
    else
    {
        current_example_mesh_=example_mesh_[current_example_index_-1];
    }
    if(!code)
    {
        if(example_eigenfunctions_loaded)
    }*/
}

void OpenGLDriver::changeSimulationMode(int code)
{
    //TO DO

}

void OpenGLDriver::loadObjectEigenfunctions(int code)
{
    OpenGLDriver* active_instance = OpenGLDriver::activeInstance();
    assert(active_instance);
    if(!active_instance->simulator_->loadObjectEigenfunctions(active_instance->object_eigen_file_name_))
    {
        std::cout<<"Error: load object eigenfunctions failed.\n";
        return;
    }
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

}  //namespace RTLB
