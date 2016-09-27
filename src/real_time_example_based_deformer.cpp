/*
 * @file: real_time_example_based_deformer.cpp
 * @brief: simulation driver for real-time simulation of example-based materials in
 *         Laplace-Beltrami shape space
 * @author: Fei Zhu,Mirror
 *
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include "coupled_quasi_harmonics.h"
#include "planes.h"
#include "real_time_example_based_deformer.h"
#include "matrixProjection.h"
using alglib::real_1d_array;
using alglib::real_2d_array;
using alglib::integer_1d_array;
using alglib::minbleicstate;
using alglib::minbleicreport;
using alglib::ae_int_t;
using alglib::minlmstate;
using alglib::minlmreport;
using alglib::mincgstate;
using alglib::mincgreport;

namespace RTLB{
RealTimeExampleBasedDeformer* RealTimeExampleBasedDeformer::active_instance_ = NULL;

RealTimeExampleBasedDeformer::RealTimeExampleBasedDeformer()
{
	if(active_instance_)
        delete active_instance_;
    active_instance_ = this;
}

RealTimeExampleBasedDeformer::~RealTimeExampleBasedDeformer()
{
	//delete tetrahedral meshes
	// for(int i=0;i<3*simulation_mesh_->getNumVertices();++i)
	// 	if(mass_[i])
	// 		delete[] mass_[i];
	// if(mass_)
	// 	delete[] mass_;
	if(simulation_mesh_)
		delete simulation_mesh_;
	for(unsigned int i=0;i<example_num_;++i)
		if(examples_[i])
			delete examples_[i];
	if(visual_mesh_)
		delete visual_mesh_;
	//delete simulation data
	delete[] displacement_;
	delete[] velocity_;
	delete[] external_force_;
	//delete reduced simulation data
	for(unsigned int i=0;i<r_;++i)
		if(reduced_basis_[i])
			delete[] reduced_basis_[i];
	delete[] reduced_basis_;
	delete[] reduced_basis_values_;
	//delete eigenfunction data and eigencoef
	delete[] object_eigencoefs_;
	delete[] object_vertex_volume_;
	delete[] object_eigenvalues_;
	delete[] object_eigencoefs_;
	for(unsigned int i=0;i<reconstruct_eigenfunction_num_;++i)
		if(object_eigenfunctions_[i])
			delete[] object_eigenfunctions_[i];
	for(unsigned int i=0;i<example_num_;++i)
	{
		if(example_eigenfunctions_[i])
		{
			for(unsigned int j=0;j<reconstruct_eigenfunction_num_;++j)
				delete[] example_eigenfunctions_[i][j];
			delete[] example_eigenfunctions_[i];
		}
		delete[] example_eigenvalues_[i];
		delete[] example_eigencoefs_[i];
		delete[] example_vertex_volume_[i];
		delete[] example_eigencoefs_[i];
	}
	delete[] example_eigencoefs_;
	delete[] example_eigenfunctions_;
	delete[] example_eigenvalues_;
	delete[] example_volume_;
	delete[] target_eigencoefs_;
	//delete cubica data
	delete[] object_cubica_weights_;
	delete[] object_cubica_elements_;
	delete[] object_LB_cubica_weights_;
	delete[] object_LB_cubica_elements_;
	//delete registration information
	for(unsigned int i=0;i<corresponding_function_num_;++i)
		if(object_corresponding_functions_[i])
			delete[] object_corresponding_functions_[i];
	for(unsigned int i=0;i<example_num_;++i)
	{
		if(example_corresponding_functions_[i])
		{
			for(unsigned int j=0;j<interpolate_eigenfunction_num_;++j)
				delete[] example_corresponding_functions_[i][j];
			delete[] example_corresponding_functions_[i];
		}
	}
	if(planes_==NULL)
		delete[] planes_;
	delete[] fixed_vertices_;
	//delete[] reduced_force_;

	for(int i=0;i<object_cubica_ele_num_;++i)
		delete[] restpos_[i];
	delete[] restpos_;
	for(int i=0;i<object_LB_cubica_ele_num_;++i)
		delete[] LB_restpos_[i];
	delete[] LB_restpos_;
	for(int i=0;i<object_LB_cubica_ele_num_;++i)
	{
		if(cubica_LB_tetsubBasis_[i])
		{
			for(int j=0;j<4;++j)
			{
				if(cubica_LB_tetsubBasis_[i][j])
					delete[] cubica_LB_tetsubBasis_[i][j];
			}
			cubica_LB_tetsubBasis_[i];
		}
	}
	// if(full_drag_force_)
	// 	delete[] full_drag_force_;
	// if(reduced_drag_force_)
	// 	delete[] reduced_drag_force_;
	//need to be checked
}

RealTimeExampleBasedDeformer* RealTimeExampleBasedDeformer::activeInstance()
{
    return active_instance_;
}
void RealTimeExampleBasedDeformer::setSimulationType(const std::string &simulation_type)
{
	simulation_type_=simulation_type.c_str();
	// std::cout<<simulation_type<<"...\n";
	if(strcmp(simulation_type_.c_str(),"fullspace")==0)
    {
        simulation_mode_=FULLSPACE;
    }
    if(strcmp(simulation_type_.c_str(),"reducedspace")==0)
    {
        simulation_mode_=REDUCEDSPACE;
    }
    if(strcmp(simulation_type_.c_str(),"UNKNOWN")==0)
    {
        std::cout<<"Error:unknown simulation mode specified."<<std::endl;
        exit(0);
    }
	if(simulation_mode_==REDUCEDSPACE)
		material_mode_=REDUCED_STVK;
	else
		material_mode_=INV_STVK;
}
void RealTimeExampleBasedDeformer::setMaterialType(const std::string &material_type)
{
	// material_type_=material_type.c_str();
	// if(strcmp(material_type_.c_str(),"reducedStVK")==0)
    // {
    //     material_mode_=REDUCED_STVK;
    // }
    // if(strcmp(material_type_.c_str(),"reducedNeoHookean")==0)
    // {
    //     material_mode_=REDUCED_NEOHOOKEAN;
    // }
	// if(strcmp(material_type_.c_str(),"StVK")==0)
	// {
	// 	material_mode_=INV_STVK;
	// }
	// if(strcmp(material_type_.c_str(),"NeoHookean")==0)
	// {
	// 	material_mode_=INV_NEOHOOKEAN;
	// }
	// std::cout<<material_mode_;
	if(simulation_mode_==REDUCEDSPACE)
		material_mode_=REDUCED_STVK;
	else
		material_mode_=INV_STVK;
	// getchar();
}
void RealTimeExampleBasedDeformer::setSolverType(const std::string &solver)
{
	solver_method_=solver;
	if(strcmp(solver_method_.c_str(),"implicitNewmark")==0)
        solver_type_=IMPLICITNEWMARK;
    if(strcmp(solver_method_.c_str(),"implicitBackwardEuler")==0)
        solver_type_=IMPLICITBACKWARDEULER;
    if(strcmp(solver_method_.c_str(),"Euler")==0)
        solver_type_=EULER;
    if(strcmp(solver_method_.c_str(),"centralDifferences")==0)
        solver_type_=CENTRALDIFFERENCES;
    if(strcmp(solver_method_.c_str(),"reducedCentralDifferences")==0)
        solver_type_=REDUCEDCENTRALDIFFERENCES;
    if(strcmp(solver_method_.c_str(),"reducedImplicitNewmark")==0)
        solver_type_=REDUCEDIMPLICITNEWMARK;
    if(strcmp(solver_method_.c_str(),"reducedImplicitBackwardEuler")==0)
        solver_type_=REDUCEDIMPLICITBACKWARDEULER;
    if(solver_type_==UNKNOWN)
    {
        std::cout<<"Error:unknown implicit solver specified."<<std::endl;
        exit(0);
    }
}
void RealTimeExampleBasedDeformer::setu(double *u)
{
	memcpy(u_,u,sizeof(double)*3*simulation_vertices_num_);
}
// void RealTimeExampleBasedDeformer::setInitialVel(double *vel)
// {
// 	memcpy(vel_initial_,vel,sizeof(double)*3*simulation_vertices_num_);
// }
// void RealTimeExampleBasedDeformer::setInitialDis(double *u)
// {
// 	memcpy(u_initial_,u,sizeof(double)*3*simulation_vertices_num_);
// }
// void RealTimeExampleBasedDeformer::setInitialForceLoad(double *force)
// {
// 	memcpy(f_ext_,force,sizeof(double)*3*simulation_vertices_num_);
// }
void RealTimeExampleBasedDeformer::setVelAfterCollision(double *vel)
{
	memcpy(collide_vel_,vel,sizeof(double)*3*simulation_vertices_num_);
}
void RealTimeExampleBasedDeformer::setupSimulation()
{
	//full space
	if(!simulation_mesh_)
	{
		std::cout<<"simulation mesh is null.\n";
		exit(0);
	}
	total_steps_=(int)((1.0/time_step_)/frame_rate_)*total_frames_;
    mesh_graph_=GenerateMeshGraph::Generate(simulation_mesh_);
    int scale_rows=1;
    mesh_graph_->GetLaplacian(&laplacian_matrix_,scale_rows);
    mesh_graph_->GetLaplacian(&laplacian_damping_matrix_,scale_rows);
    laplacian_damping_matrix_->ScalarMultiply(damping_laplacian_coef_);
	// simulation_vertices_num_=simulation_mesh_->getNumVertices();

	// full_drag_force_=new double[3*simulation_vertices_num_];
	// memset(full_drag_force_,0.0,sizeof(double)*3*simulation_vertices_num_);

	u_=new double[3*simulation_vertices_num_];
	memset(u_,0.0,sizeof(double)*3*simulation_vertices_num_);
	vel_=new double[3*simulation_vertices_num_];
	memset(vel_,0.0,sizeof(double)*3*simulation_vertices_num_);
	collide_vel_=new double[3*simulation_vertices_num_];
	memset(collide_vel_,0.0,sizeof(double)*3*simulation_vertices_num_);
	u_initial_=new double[3*simulation_vertices_num_];
	memset(u_initial_,0.0,sizeof(double)*3*simulation_vertices_num_);
	vel_initial_=new double[3*simulation_vertices_num_];
	memset(vel_initial_,0.0,sizeof(double)*3*simulation_vertices_num_);
	f_ext_=new double[3*simulation_vertices_num_];
	memset(f_ext_,0.0,sizeof(double)*3*simulation_vertices_num_);
	example_f_=new double[3*simulation_vertices_num_];
	memset(example_f_,0.0,sizeof(double)*3*simulation_vertices_num_);
	//reduced space
	// reduced_drag_force_=new double[r_];
	// memset(reduced_drag_force_,0.0,sizeof(double)*r_);
	q_=new double[r_];
	memset(q_,0.0,sizeof(double)*r_);
	temp_q_=new double[r_];
	memset(temp_q_,0.0,sizeof(double)*r_);
	temp_grad_=new double[r_];
	memset(temp_grad_,0.0,sizeof(double)*r_);
	qvel_=new double[r_];
	memset(qvel_,0.0,sizeof(double)*r_);
	qaccel_=new double[r_];
	memset(qaccel_,0.0,sizeof(double)*r_);
	fq_=new double[r_];
	memset(fq_,0.0,sizeof(double)*r_);
	fqBase_=new double[r_];
	memset(fqBase_,0.0,sizeof(double)*r_);
	fq_ext_=new double[r_];
	memset(fq_ext_,0.0,sizeof(double)*r_);

	restpos_ = new double*[object_cubica_ele_num_];
	for(int i=0;i<object_cubica_ele_num_;++i)
	{
		restpos_[i] = new double[12];//3n*1
		int ele=object_cubica_elements_[i];
		for(int j=0;j<4;++j)
		{
			int global_idx=simulation_mesh_->getVertexIndex(ele,j);
			restpos_[i][3*j]=(*simulation_mesh_->getVertex(global_idx))[0];
			restpos_[i][3*j+1]=(*simulation_mesh_->getVertex(global_idx))[1];
			restpos_[i][3*j+2]=(*simulation_mesh_->getVertex(global_idx))[2];
		}
	}
	if(isload_LB_cubica_)
	{
		LB_restpos_ = new double*[object_LB_cubica_ele_num_];
		for(int i=0;i<object_LB_cubica_ele_num_;++i)
		{
			LB_restpos_[i] = new double[12];//3n*1
			int ele=object_LB_cubica_elements_[i];
			for(int j=0;j<4;++j)
			{
				int global_idx=simulation_mesh_->getVertexIndex(ele,j);
				LB_restpos_[i][3*j]=(*simulation_mesh_->getVertex(global_idx))[0];
				LB_restpos_[i][3*j+1]=(*simulation_mesh_->getVertex(global_idx))[1];
				LB_restpos_[i][3*j+2]=(*simulation_mesh_->getVertex(global_idx))[2];
			}
		}
		deformed_=new double[12];
	}

	// double *mass_=new double[3*simulation_vertices_num_];
	// memset(mass_,0.0,sizeof(double)*3*simulation_vertices_num_);
	//temp compute mass for each vertex
	// for(int el=0;el<simulation_mesh_->getNumElements();++el)
	// {
	// 	for(int i=0;i<4;++i)
	// 	{
	// 		int vertID=simulation_mesh_->getVertexIndex(el,i);
	// 		mass_[vertID]+=0.25*simulation_mesh_->getElementVolume(el)*simulation_mesh_->getElementDensity(el)/3.0;
	// 	}
	// }
	// std::cout<<"mass:\n";
	//temp compute mass for each vertex-end
    //create force models, to be used by the integrator,deformable_object_type_=INVERTIBLEFEM,neohookean material
    std::cout<<"Creating force model:\n";
	tet_mesh_=dynamic_cast<TetMesh*>(simulation_mesh_);
	if(tet_mesh_==NULL)
	{
		std::cout<<"Error: the input mesh is not a tet mesh (CLFEM deformable model).\n";
		exit(1);
	}
	if(simulation_mode_==FULLSPACE)
	{
		if(material_mode_==INV_STVK)
		{
			isotropic_material_=new StVKIsotropicMaterial(tet_mesh_);
			cout<<"Invertible material:StVK.\n";
		}
		else if(material_mode_==INV_NEOHOOKEAN)
		{
			isotropic_material_=new NeoHookeanIsotropicMaterial(tet_mesh_);
			std::cout<<"Invertible material: neo-Hookean.\n";
		}
		isotropic_hyperelastic_fem_=new IsotropicHyperelasticFEM(tet_mesh_,isotropic_material_,principal_stretch_threshold_,add_gravity_,gravity_);
		force_model_=new IsotropicHyperelasticFEMForceModel(isotropic_hyperelastic_fem_);
	}
	else
	{
		//reduced neohookean
		if(material_mode_==REDUCED_NEOHOOKEAN)
		{
			reduced_neoHookean_force_model_ = new ReducedNeoHookeanForceModel(r_,simulation_mesh_,U_,object_cubica_ele_num_,
									object_cubica_weights_,object_cubica_elements_,restpos_,add_gravity_,gravity_);
			reduced_force_model_=reduced_neoHookean_force_model_;
		}
		//reduced stvk cubature
		if(material_mode_==REDUCED_STVK)
		{
			reduced_stvk_cubature_force_model_ = new ReducedStVKCubatureForceModel(r_,simulation_mesh_,U_,object_cubica_ele_num_,
										object_cubica_weights_,object_cubica_elements_,restpos_,add_gravity_,gravity_);
			reduced_force_model_=reduced_stvk_cubature_force_model_;
		}

	}
    //initialize the integrator
    std::cout<<"Initializing the integrator, n= "<<simulation_vertices_num_<<".\n";
    std::cout<<"Solver type: "<<solver_type_<<".\n";
    integrator_base_sparse_=NULL;
    integrator_base_dense_=NULL;

    if(solver_type_==IMPLICITNEWMARK)
    {
        implicit_newmark_sparse_=new ImplicitNewmarkSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,
                                                        positive_definite_,fixed_dofs_num_,fixed_dofs_,
                                                        damping_mass_coef_,damping_stiffness_coef_,max_iterations_,
                                                        integrator_epsilon_,newmark_beta_,newmark_gamma_,solver_threads_num_);

        integrator_base_sparse_=implicit_newmark_sparse_;

    }
    else if(solver_type_==IMPLICITBACKWARDEULER)
    {
        implicit_newmark_sparse_=new ImplicitBackwardEulerSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,
                                                        positive_definite_,fixed_dofs_num_,fixed_dofs_,
                                                        damping_mass_coef_,damping_stiffness_coef_,max_iterations_,
                                                        integrator_epsilon_,solver_threads_num_);
        integrator_base_sparse_=implicit_newmark_sparse_;
    }
    else if(solver_type_==EULER)
    {
        int symplectic=0;
        integrator_base_sparse_=new EulerSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,symplectic,
                                                fixed_dofs_num_,fixed_dofs_,damping_mass_coef_);
    }
    // else if(solver_type_==CENTRALDIFFERENCES)
    // {
    //     integrator_base_sparse_=new CentralDifferencesSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,
    //                                                         fixed_dofs_num_,fixed_dofs_,damping_mass_coef_,damping_stiffness_coef_,
    //                                                         central_difference_tangential_damping_update_mode_,solver_threads_num_);
    // }
    // else if(solver_type_==REDUCEDCENTRALDIFFERENCES)
    // {
    //     central_differences_dense_ = new CentralDifferencesDense(r_,time_step_,reduced_mass_matrix_,reduced_force_model_,damping_mass_coef_,
    //                                                         damping_stiffness_coef_,central_difference_tangential_damping_update_mode_);
    //     integrator_base_dense_=central_differences_dense_;
    //     simulation_mode_=REDUCEDSPACE;
    // }
    else if(solver_type_==REDUCEDIMPLICITNEWMARK)
    {
		std::cout<<"REDUCEDIMPLICITNEWMARK:\n";
        implicit_newmark_dense_ = new ImplicitNewmarkDense(r_,time_step_,reduced_mass_matrix_,reduced_force_model_,
                                                        ImplicitNewmarkDense::positiveDefiniteMatrixSolver,damping_mass_coef_,damping_stiffness_coef_,
                                                        max_iterations_,integrator_epsilon_,newmark_beta_,newmark_gamma_);
        integrator_base_dense_=implicit_newmark_dense_;
        simulation_mode_=REDUCEDSPACE;
		std::cout<<"REDUCEDIMPLICITNEWMARK-end\n";
    }
    else if(solver_type_==REDUCEDIMPLICITBACKWARDEULER)
    {
        implicit_backward_euler_dense_ = new ImplicitBackwardEulerDense(r_,time_step_,reduced_mass_matrix_,reduced_force_model_,
                                                                ImplicitBackwardEulerDense::positiveDefiniteMatrixSolver,damping_mass_coef_,
                                                                damping_stiffness_coef_,max_iterations_,integrator_epsilon_);
        integrator_base_dense_=implicit_backward_euler_dense_;
        simulation_mode_=REDUCEDSPACE;
    }
    else
    {
    }
	if(simulation_mode_==FULLSPACE)
	{
		//load initial velocity
		std::cout<<"initial_velocity_file_name_:"<<initial_velocity_file_name_<<"\n";
		if(initial_velocity_file_name_!="none")
		{
			int m1,n1;
			int status=ReadMatrixFromDiskTextFile(initial_velocity_file_name_.c_str(),&m1,&n1,&vel_initial_);
			if(status)
				exit(1);
			if((m1!=3*simulation_vertices_num_)||(n1!=1))
			{
				cout<<"Error:initial velocity matrix size mismatch.\n";
				exit(1);
			}
			cout<<"Initial velocity loaded succeed.\n";
		}
		//load initial velocity
		// if(strcmp(initial_position_file_name_,"none")!=0)
		// {
		// 	int m1,n1;
		// 	int status=ReadMatrixFromDiskTextFile(initial_position_file_name_,&m1,&n1,&u_initial_);
		// 	if(status)
		// 		exit(1);
		// 	if((m1!=3*simulation_vertices_num_)||(n1!=1))
		// 	{
		// 		cout<<"Error:initial position matrix size mismatch.\n";
		// 		exit(1);
		// 	}
		// 	cout<<"Initial position loaded succeed.\n";
		// }
	}

	//load initial force
	std::cout<<"force_loads_file_name_:"<<force_loads_file_name_<<"\n";
	if(force_loads_file_name_!="none")
	{
		int m1;
		int status=ReadMatrixFromDiskTextFile(force_loads_file_name_.c_str(),&m1,&force_loads_num_,&force_loads_);
		if(status)
			exit(1);
		if(m1!=3*simulation_vertices_num_)
		{
			cout<<"Mismatch in the dimension of the force load matrix.\n";
			exit(1);
		}
		cout<<"Initial force loaded succeed.\n";
	}

	//load local affected vertices for local examples-based simulation
	std::cout<<"object_affected_vertices_file_name_:"<<object_affected_vertices_file_name_<<"\n";
	if(object_affected_vertices_file_name_=="none")
    {
        object_affected_vertices_num_=0;
        object_affected_vertices_=NULL;
    }
    else
    {
		// std::cout<<"...3\n";
        if(LoadList::load(object_affected_vertices_file_name_.c_str(),&object_affected_vertices_num_,&object_affected_vertices_)!=0)
        {
            cout<<"Error reading object vertices affected by local examples.\n";
            exit(1);
        }
        LoadList::sort(object_affected_vertices_num_,object_affected_vertices_);
        cout<<"Loaded "<<object_affected_vertices_num_<<" object vertices affected by local examples. They are: \n";
        LoadList::print(object_affected_vertices_num_,object_affected_vertices_);
    }
	if(example_num_>0)
	{
		example_affected_vertices_num_=new int[example_num_];
        example_affected_vertices_=new int*[example_num_];
        bool empty_file_name;
        if(example_affected_vertices_file_base_=="none")
            empty_file_name=true;
        else
            empty_file_name=false;
        for(int i=0;i<example_num_;++i)
        {
            if(empty_file_name)
            {
                example_affected_vertices_num_[i]=0;
                example_affected_vertices_[i]=NULL;
                continue;
            }
            stringstream adaptor;
            string example_idx_str;
            adaptor<<i;
            adaptor>>example_idx_str;
            string cur_example_affected_vertices_file_name=example_affected_vertices_file_base_+example_idx_str+".bou";
            if(LoadList::load(cur_example_affected_vertices_file_name.c_str(),&example_affected_vertices_num_[i],&example_affected_vertices_[i])!=0)
            {
                cout<<"Error reading vertices belong to local region on example"<<example_idx_str<<".\n";
                exit(1);
            }
            LoadList::sort(example_affected_vertices_num_[i],example_affected_vertices_[i]);
            cout<<"Loaded "<<example_affected_vertices_num_[i]<<" vertices belong to local region on example"<<example_idx_str<<". They are:\n";
            LoadList::print(example_affected_vertices_num_[i],example_affected_vertices_[i]);
        }
	}
    //set integration parameters
    if(simulation_mode_==REDUCEDSPACE)
    {
        integrator_base_=integrator_base_dense_;
	    if(integrator_base_==NULL)
        {
            std::cout<<"Error: failed to initialize reduced numerical integrator.\n";
            exit(1);
        }
	    integrator_base_->SetTimestep(time_step_);
		integrator_base_->SetState(q_,qvel_);
		preAllocateLocalFrameCorrespondingVertices();
		if(!with_constrains_)
			rigidBodyPreComputation();
    }
    else
    {
        integrator_base_=integrator_base_sparse_;
        if(integrator_base_==NULL)
        {
            std::cout<<"Error: failed to initialize numerical integrator.\n";
            exit(1);
        }
        integrator_base_sparse_->SetDampingMatrix(laplacian_damping_matrix_);
        integrator_base_->ResetToRest();
        integrator_base_->SetState(u_initial_,vel_initial_);
        integrator_base_->SetTimestep(time_step_);
    }
	//load initial tet mesh
	std::cout<<"load initial tet mesh: "<<initial_tetmesh_file_name_<<"\n";
	if(initial_tetmesh_file_name_!="none")
	{
		VolumetricMesh *init_simulation_mesh=VolumetricMeshLoader::load(initial_tetmesh_file_name_.c_str());
		for(unsigned int i=0;i<simulation_vertices_num_;++i)
		{
			Vec3d vert_before=*init_simulation_mesh->getVertex(i);
			Vec3d vert_after=*simulation_mesh_->getVertex(i);
			u_[3*i]=vert_before[0]-vert_after[0];
			u_[3*i+1]=vert_before[1]-vert_after[1];
			u_[3*i+2]=vert_before[2]-vert_after[2];
		}
		std::cout<<"init-q:\n";
		for(int i=0;i<r_;++i)
			std::cout<<q_[i]<<",";
		if(simulation_mode_==REDUCEDSPACE)
	    {
			modal_matrix_->ProjectVector(u_,q_);

			integrator_base_dense_->SetState(q_);
		}
		else
		{
			//memset(vel_initial_,0.0,sizeof(double)*simulation_vertices_num_*3);
			integrator_base_->SetState(u_,vel_initial_);
		}
		delete init_simulation_mesh;
		std::cout<<"u_:\n";
		for(int i=0;i<30;++i)
		{
			std::cout<<u_[i]<<",";
		}
		// getchar();
	}
	std::cout<<"init example simulation\n";
	// getchar();
	//initial for example simulation
	target_eigencoefs_=new Vec3d[interpolate_eigenfunction_num_];
	memset(target_eigencoefs_,0.0,sizeof(Vec3d)*interpolate_eigenfunction_num_);
	energy_grad_=new Vec3d[interpolate_eigenfunction_num_];
	memset(energy_grad_,0.0,sizeof(Vec3d)*interpolate_eigenfunction_num_);
	// target_eigencoefs_diff_=new Vec3d[interpolate_eigenfunction_num_];
	// memset(target_eigencoefs_diff_,0.0,sizeof(Vec3d)*interpolate_eigenfunction_num_);
	object_current_eigencoefs_ = new Vec3d[interpolate_eigenfunction_num_];
	memset(object_current_eigencoefs_,0.0,sizeof(Vec3d)*interpolate_eigenfunction_num_);
	example_guided_deformation_ = new double[3*simulation_vertices_num_];
	memset(example_guided_deformation_,0.0,sizeof(double)*3*simulation_vertices_num_);
	temp_example_guided_deformation_ = new double[3*simulation_vertices_num_];
	memset(temp_example_guided_deformation_,0.0,sizeof(double)*3*simulation_vertices_num_);
	if(simulation_mode_==REDUCEDSPACE)
	{

		target_deformation_ = new double[3*simulation_vertices_num_];
		memset(target_deformation_,0.0,sizeof(double)*3*simulation_vertices_num_);
		// example_based_LB_forces_ = new double[3*interpolate_eigenfunction_num_];
		// memset(example_based_LB_forces_,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
		example_based_forces_ = new double[3*simulation_vertices_num_];
		memset(example_based_forces_,0.0,sizeof(double)*3*simulation_vertices_num_);
		example_based_q_ = new double[r_];
		memset(example_based_q_,0.0,sizeof(double)*r_);
		example_based_fq_ = new double[r_];
		memset(example_based_fq_,0.0,sizeof(double)*r_);
	}
	std::cout<<"initSimulation-end.\n";

	// std::cout<<"reduced:testEnergyGradients:\n";
	// reduced_neoHookean_force_model_->testEnergyGradients();
	// getchar();
	// std::cout<<"testInternalForceGradients:\n";
	// reduced_neoHookean_force_model_->testObjectiveGradients();
	// std::cout<<"testInternalForceGradients-end:\n";
	// getchar();
}
void RealTimeExampleBasedDeformer::rigidBodyPreComputation()
{
	if(!simulation_mesh_)
	{
		std::cout<<"simulation mesh unloaded!\n";
		exit(0);
	}
	//compute mass for each vertex
	double *row_mass=new double[3*simulation_mesh_->getNumVertices()];
	memset(row_mass,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	vert_mass_=new double[simulation_mesh_->getNumVertices()];
	memset(vert_mass_,0.0,sizeof(double)*simulation_mesh_->getNumVertices());
	// vert_f_ext_=new double[3*simulation_mesh_->getNumVertices()];
	// memset(vert_f_ext_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	mass_matrix_->SumRowEntries(row_mass);
	for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	{
		vert_mass_[i]+=row_mass[3*i]+row_mass[3*i+1]+row_mass[3*i+2];
		total_mass_+=vert_mass_[i];
	}
	std::cout<<"total_mass_:"<<total_mass_<<"\n";

	delete[] row_mass;
	for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	{
		rigid_center_[0]+=vert_mass_[i]*(*simulation_mesh_->getVertex(i))[0]/total_mass_;
		rigid_center_[1]+=vert_mass_[i]*(*simulation_mesh_->getVertex(i))[1]/total_mass_;
		rigid_center_[2]+=vert_mass_[i]*(*simulation_mesh_->getVertex(i))[2]/total_mass_;
	}
	// std::cout<<rigid_center_[0]<<","<<rigid_center_[1]<<","<<rigid_center_[2]<<"\n";
	// rigid_center_[0]=1.35967;
	// rigid_center_[1]=-3.72632;
	// rigid_center_[2]=1.24794;
	// initial_inertia_tensor_[0]=2873.9;
	// initial_inertia_tensor_[1]=initial_inertia_tensor_[3]=559.417;
	// initial_inertia_tensor_[2]=initial_inertia_tensor_[6]=-103.085;
	// initial_inertia_tensor_[4]=1371.21;
	// initial_inertia_tensor_[5]=initial_inertia_tensor_[7]=-167.487;
	// initial_inertia_tensor_[8]=3371.17;
	// initial_inertia_tensor_[0]=initial_inertia_tensor_[4]=initial_inertia_tensor_[8]=1.0;
	// initial_inertia_tensor_[1]=initial_inertia_tensor_[2]=initial_inertia_tensor_[3]=initial_inertia_tensor_[5]=0.0;
	// initial_inertia_tensor_[6]=initial_inertia_tensor_[7]=0.0;
	rigid_=new RigidBody_GeneralTensor(total_mass_,initial_inertia_tensor_);
	// rigid_=new RigidBody(total_mass_,initial_inertia_tensor_[0],initial_inertia_tensor_[4],initial_inertia_tensor_[8]);
	rigid_->SetPosition(rigid_center_[0],rigid_center_[1],rigid_center_[2]);
	rigid_->SetLinearDamping(0.2);
	rigid_->SetRotationalDamping(0.3);
	rigid_->SetVelocity(initial_rigid_vel_[0],initial_rigid_vel_[1],initial_rigid_vel_[2]);
	// rigid_->SetAngularVelocity(0.0,0.0,0.3);
	linear_velocity_[0]=initial_rigid_vel_[0];
	linear_velocity_[1]=initial_rigid_vel_[1];
	linear_velocity_[2]=initial_rigid_vel_[2];

	local_reference_=new double[3*simulation_mesh_->getNumVertices()];
	memset(local_reference_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	local_u_=new double[3*simulation_mesh_->getNumVertices()];
	memset(local_u_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	global_u_=new double[3*simulation_mesh_->getNumVertices()];
	memset(global_u_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	Uqvel_=new double[3*simulation_mesh_->getNumVertices()];
	memset(Uqvel_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	Uqaccel_=new double[3*simulation_mesh_->getNumVertices()];
	memset(Uqaccel_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	{
		for(int j=0;j<3;++j)
			local_reference_[3*i+j]=(*simulation_mesh_->getVertex(i))[j]-rigid_center_[j];
	}
}
void RealTimeExampleBasedDeformer::advanceStep()
{
	if(simulation_mode_==FULLSPACE)
		fullspaceSimulation();
	else
	{
		if(with_constrains_)
			reducedspaceSimulationWithConstraints();
		else
			reducedspaceSimulationWithoutConstraints();
	}
}
void RealTimeExampleBasedDeformer::setExternalForces(double *ext_forces)
{
	memcpy(f_ext_,ext_forces,sizeof(double)*3*simulation_mesh_->getNumVertices());
	// memset(f_ext_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
}

void RealTimeExampleBasedDeformer::setReducedExternalForces(double *ext_forces)
{
	memcpy(fq_,ext_forces,sizeof(double)*r_);
}
void RealTimeExampleBasedDeformer::fullspaceSimulation()
{
	// std::cout<<"fullspaceSimulation:\n";
	// memcpy(f_ext_,full_drag_force_,sizeof(double)*3*simulation_vertices_num_);
    static int count=1;
	static double temp_coef=1.0;
	if(time_step_counter_<force_loads_num_)
	// if(count_num<force_loads_num_)
	{
		// std::cout<<"External forces read from the text input file.\n";
		for(int i=0;i<3*simulation_vertices_num_;++i)
			f_ext_[i]+=force_loads_[ELT(3*simulation_vertices_num_,i,time_step_counter_)];
		// count_num++;

	}
	//apply the force loads caused by the examples
	// static int count=1;
	if(enable_example_simulation_)
	{
		PerformanceCounter project_counter;
		project_counter.StartCounter();
		projectOnEigenFunctions(simulation_mesh_,u_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
								interpolate_eigenfunction_num_,object_current_eigencoefs_,
								object_affected_vertices_num_,object_affected_vertices_);
		project_counter.StopCounter();
		total_projection_time_+=project_counter.GetElapsedTime();
		// if(count==1)
		// {
		// 	for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 		for(int j=0;j<3;++j)
		// 		{
		// 			object_eigencoefs_[i][j]=example_eigencoefs_[0][i][j];
		// 			// object_current_eigencoefs_[i][j]=example_eigencoefs_[0][i][j];
		// 		}
		// 	count++;
		// }
		// projectReducedSubspaceOnLBSubspace(fq_,object_surface_eigencoefs_);
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			object_current_eigencoefs_[i]=object_current_eigencoefs_[i]-object_eigencoefs_[i]+example_eigencoefs_[0][i];
		memcpy(target_eigencoefs_,object_current_eigencoefs_,sizeof(Vec3d)*interpolate_eigenfunction_num_);
		PerformanceCounter target_counter;
		target_counter.StartCounter();
		//find target shape for simulation object
		projectOnExampleManifold(object_current_eigencoefs_,target_eigencoefs_);
		target_counter.StopCounter();
		total_target_time_+=target_counter.GetElapsedTime();

		// for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 	target_eigencoefs_[i]=example_eigencoefs_[1][i];
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			target_eigencoefs_[i]=object_current_eigencoefs_[i]-target_eigencoefs_[i];
		// std::cout<<"-----------------before-reconstruct:\n";
		PerformanceCounter reconstruct_counter;
		reconstruct_counter.StartCounter();
		reconstructFromEigenCoefs(target_eigencoefs_,example_guided_deformation_);
		reconstruct_counter.StopCounter();
		total_reconstruction_time_+=reconstruct_counter.GetElapsedTime();
		// std::cout<<"-----------------after-reconstruct:\n";

		if(object_affected_vertices_num_>0)
		{
			memcpy(temp_example_guided_deformation_,example_guided_deformation_,sizeof(double)*3*simulation_vertices_num_);
			memset(example_guided_deformation_,0.0,sizeof(double)*3*simulation_vertices_num_);
			for(int i=0;i<object_affected_vertices_num_;++i)
			{
				int idx=object_affected_vertices_[i]-1;
				example_guided_deformation_[3*idx+0]=temp_example_guided_deformation_[3*idx+0];
				example_guided_deformation_[3*idx+1]=temp_example_guided_deformation_[3*idx+1];
				example_guided_deformation_[3*idx+2]=temp_example_guided_deformation_[3*idx+2];
			}
		}

		for(int i=0;i<3*simulation_vertices_num_;++i)
			example_f_[i]=temp_coef*example_stiffness_scale_*example_guided_deformation_[i];
		// temp_coef=temp_coef+0.01;
		//the internal forces are returned with the sign corresponding to f_int(x) on the left side of the equation M * x'' + f_int(x) = f_ext
		//so they must to be subtracted from external force
		for(int i=0;i<3*simulation_vertices_num_;++i)
			f_ext_[i]-=example_f_[i];
		// for(int i=0;i<3*simulation_vertices_num_;++i)
		// 	f_ext_[i]+=(-1.0)*example_stiffness_scale_*example_guided_deformation_[i];
		//vert_idx:0
		// std::cout<<"0:"<<f_ext_[0]<<","<<f_ext_[1]<<","<<f_ext_[2]<<"\n";
		// std::cout<<"1:"<<f_ext_[3*1+0]<<","<<f_ext_[3*1+1]<<","<<f_ext_[3*1+2]<<"\n";
		// std::cout<<"4:"<<f_ext_[3*4+0]<<","<<f_ext_[3*4+1]<<","<<f_ext_[3*4+2]<<"\n";
		// std::cout<<"5:"<<f_ext_[3*5+0]<<","<<f_ext_[3*5+1]<<","<<f_ext_[3*5+2]<<"\n";
		// getchar();
		// if(time_step_counter_%timer_sample_interval_==0&&time_step_counter_>=timer_sample_interval_)
        // {
		// 	// std::cout<<timer_sample_interval_<<"------\n";
		// 	std::cout<<"Total project solver time: "<<total_projection_time_/timer_sample_interval_<<" s.\n";
        //     std::cout<<"Total target solver time: "<<total_target_time_/timer_sample_interval_<<" s.\n";
		// 	std::cout<<"Total reconstruct time: "<<total_reconstruction_time_/timer_sample_interval_<<" s.\n";
        //     total_projection_time_=total_target_time_=total_reconstruction_time_= 0;
        //     // getchar();
        // }

	}
	// else
	// {
	// 	for(int i=0;i<3*simulation_vertices_num_;++i)
	// 		f_ext_[i]*=damping_example_stiffness_;
	// }
	//set forces to the integrator
	// std::cout<<"c\n";
	integrator_base_sparse_->SetExternalForces(f_ext_);
	// std::cout<<"d\n";
	//time step the dynamics
	// PerformanceCounter step_counter;
	// step_counter.StartCounter();
	if(time_step_counter_ < total_steps_)
	{
		int code=integrator_base_->DoTimestep();
		std::cout<<".";fflush(NULL);
		++time_step_counter_;
	}
	// std::cout<<"e\n";
	// step_counter.StopCounter();
	// total_time_+=step_counter.GetElapsedTime();
	//
	// if(time_step_counter_%timer_sample_interval_==0&&time_step_counter_>=timer_sample_interval_)
	// {
	// 	std::cout<<"step computation time: "<<total_time_/timer_sample_interval_<<" s.\n";
	// 	total_time_= 0;
	// }

	memcpy(u_,integrator_base_->Getq(),sizeof(double)*3*(simulation_vertices_num_));
	// for(int i=0;i<30;++i)
	// std::cout<<u_[i]<<",";
	if(time_step_counter_==(total_steps_-1))
	{
		std::cout<<"Total target configuration solver time: "<<total_target_time_<<" s.\n";
	}
}
void RealTimeExampleBasedDeformer::reducedspaceSimulationWithConstraints()
{
	memset(fq_,0.0,sizeof(double)*r_);
	static int count=1;
	// reduced_stvk_cubature_force_model_->testEnergyGradients();
	// testEnergyGradients();
	// std::cout<<"2\n";

	if(enable_example_simulation_)
	{

		PerformanceCounter project_counter;
		project_counter.StartCounter();
		projectOnEigenFunctions(simulation_mesh_,u_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
								interpolate_eigenfunction_num_,object_current_eigencoefs_,
								object_affected_vertices_num_,object_affected_vertices_);
		project_counter.StopCounter();
		total_projection_time_+=project_counter.GetElapsedTime();
		// if(count==1)
		// {
		// 	for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 		for(int j=0;j<3;++j)
		// 		{
		// 			object_eigencoefs_[i][j]=example_eigencoefs_[0][i][j];
		// 			// object_current_eigencoefs_[i][j]=example_eigencoefs_[0][i][j];
		// 		}
		// 	count++;
		// }
		// std::cout<<"object-coefs:\n";
		// for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 	std::cout<<object_current_eigencoefs_[i]<"\n";
		// std::cout<<"\n";
		// //
		// std::cout<<"ex0-coefs:\n";
		// for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 	std::cout<<example_eigencoefs_[0][i]<"\n";
		// std::cout<<"\n";
		//
		// std::cout<<"ex1-coefs:\n";
		// for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 	std::cout<<example_eigencoefs_[1][i]<"\n";
		// std::cout<<"\n";

		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			object_current_eigencoefs_[i]=object_current_eigencoefs_[i]-object_eigencoefs_[i]+example_eigencoefs_[0][i];
		memcpy(target_eigencoefs_,object_current_eigencoefs_,sizeof(Vec3d)*interpolate_eigenfunction_num_);
		// PerformanceCounter target_counter;
		// target_counter.StartCounter();
		PerformanceCounter target_counter;
		target_counter.StartCounter();
		projectOnExampleManifold(object_current_eigencoefs_,target_eigencoefs_);

		// std::cout<<"target-coefs:\n";
		// for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 	std::cout<<target_eigencoefs_[i]<"\n";
		// std::cout<<"\n";
// getchar();
		target_counter.StopCounter();
		total_target_time_+=target_counter.GetElapsedTime();
		// target_counter.StopCounter();
		// total_target_time_+=target_counter.GetElapsedTime();
		// std::cout<<"target_counter:"<<target_counter.GetElapsedTime()<<"\n";
		//--------------------------------------------------------------------------------------------
		//----------------used for reduced elastic force-reconstruct target configuration first--------
		// reconstructFromEigenCoefs(target_eigencoefs_,target_reconstruction_);
		//---------------------------------------------------------------------------------------------
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			target_eigencoefs_[i]=object_current_eigencoefs_[i]-target_eigencoefs_[i];
		PerformanceCounter reconstruct_counter;
		reconstruct_counter.StartCounter();
		reconstructFromEigenCoefs(target_eigencoefs_,example_guided_deformation_);
		reconstruct_counter.StopCounter();
		total_reconstruction_time_+=reconstruct_counter.GetElapsedTime();
		// reconstruct_counter.StopCounter();
		// std::cout<<"reconstruct_counter:"<<reconstruct_counter.GetElapsedTime()<<"\n";
		if(object_affected_vertices_num_>0)
		{
			memcpy(temp_example_guided_deformation_,example_guided_deformation_,sizeof(double)*3*simulation_vertices_num_);
			memset(example_guided_deformation_,0.0,sizeof(double)*3*simulation_vertices_num_);
			for(int i=0;i<object_affected_vertices_num_;++i)
			{
				int idx=object_affected_vertices_[i]-1;
				example_guided_deformation_[3*idx+0]=temp_example_guided_deformation_[3*idx+0];
				example_guided_deformation_[3*idx+1]=temp_example_guided_deformation_[3*idx+1];
				example_guided_deformation_[3*idx+2]=temp_example_guided_deformation_[3*idx+2];
			}
		}
		//method 1-end
		//method 2
		// reconstructFromEigenCoefs(target_eigencoefs_,2);//reference+target_deformation=target configuration
		// if(object_affected_vertices_num_>0)
		// {
		// 	// memcpy(temp_example_guided_deformation_,example_guided_deformation_,sizeof(double)*3*simulation_vertices_num_);
		// 	memset(example_guided_deformation_,0.0,sizeof(double)*3*simulation_vertices_num_);
		// 	for(int i=0;i<object_affected_vertices_num_;++i)
		// 	{
		// 		int idx=object_affected_vertices_[i]-1;
		// 		example_guided_deformation_[3*idx+0]=(*simulation_mesh_->getVertex(idx))[0]+u_[3*idx+0]-target_reconstruction_[3*idx+0];
		// 		example_guided_deformation_[3*idx+1]=(*simulation_mesh_->getVertex(idx))[1]+u_[3*idx+1]-target_reconstruction_[3*idx+1];
		// 		example_guided_deformation_[3*idx+2]=(*simulation_mesh_->getVertex(idx))[2]+u_[3*idx+2]-target_reconstruction_[3*idx+2];
		// 	}
		// }
		// else
		// {
		// 	// getchar();
		// 	for(int i=0;i<simulation_mesh_->getNumVertices();++i)
		// 	{
		// 		example_guided_deformation_[3*i+0]=(*simulation_mesh_->getVertex(i))[0]+u_[3*i+0]-target_reconstruction_[3*i+0];
		// 		example_guided_deformation_[3*i+1]=(*simulation_mesh_->getVertex(i))[1]+u_[3*i+1]-target_reconstruction_[3*i+1];
		// 		example_guided_deformation_[3*i+2]=(*simulation_mesh_->getVertex(i))[2]+u_[3*i+2]-target_reconstruction_[3*i+2];
		// 		// target_deformation_[3*i+0]=target_reconstruction_[3*i+0]-(*simulation_mesh_->getVertex(i))[0];
		// 		// target_deformation_[3*i+1]=target_reconstruction_[3*i+1]-(*simulation_mesh_->getVertex(i))[1];
		// 		// target_deformation_[3*i+2]=target_reconstruction_[3*i+2]-(*simulation_mesh_->getVertex(i))[2];
		// 	}
		// }
		//method 2-end
		// linear force
		for(int i=0;i<3*simulation_vertices_num_;++i)
		{
			f_ext_[i]=(-1.0)*example_stiffness_scale_*example_guided_deformation_[i];
			// if(f_ext_[i]>1.0e-6)
			// 	std::cout<<f_ext_[i]<<",";
		}
// getchar();
		//----------------------------------------------------------------------------
		//--------------reduced stvk elastic force-------------------------------------
		// modal_matrix_->ProjectVector(example_guided_deformation_,example_based_q_);
		//
		// if(material_mode_==REDUCED_STVK)
		// 	// reduced_stvk_cubature_force_model_->computeReducedElasticInternalForce(example_guided_deformation_,fq_,target_reconstruction_);//target_deformation_);
		// 	reduced_stvk_cubature_force_model_->computeReducedInternalForce(example_based_q_,fq_);
		//-----------------------------------------------------------------------------
		// if(material_mode_==REDUCED_NEOHOOKEAN)
		// 	reduced_neoHookean_force_model_->computeReducedElasticInternalForce(example_guided_deformation_,example_based_fq_,target_reconstruction_);//target_deformation_);

			    // double *global_forces=new double[3*simulation_mesh_->getNumVertices()];
			    // memset(global_forces,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
		// reconsturctFromEigenCoefs(target_eigencoefs_diff_,target_deformation_);
		// modal_matrix_->ProjectVector(target_deformation_,target_q_);
		// reduced_neoHookean_force_model_->computeReducedInternalForce(target_q_,example_based_fq_);
		//....test
		// for(int i=0;i<3*simulation_vertices_num_;++i)
		// 	f_ext_[i]=-0.0001*global_forces[i];
		// delete[] global_forces;


		// if(time_step_counter_%timer_sample_interval_==0&&time_step_counter_>=timer_sample_interval_)
        // {
		// 	// std::cout<<timer_sample_interval_<<"------\n";
		// 	std::cout<<"Total project solver time: "<<total_projection_time_/timer_sample_interval_<<" s.\n";
        //     std::cout<<"Total target solver time: "<<total_target_time_/timer_sample_interval_<<" s.\n";
		// 	std::cout<<"Total reconstruct time: "<<total_reconstruction_time_/timer_sample_interval_<<" s.\n";
        //     total_projection_time_=total_target_time_=total_reconstruction_time_= 0;
        //     // getchar();
        // }
	}
	//linear-force
	modal_matrix_->ProjectVector(f_ext_,fq_);
	// for(int i=0;i<r_;++i)
	// 	fq_[i]*=-2.0;
	integrator_base_dense_->SetExternalForces(fq_);
	// getchar();

	// for(int i=0;i<r_;++i)
	// std::cout<<q_[i]<<",";
	// std::cout<<"\n";
	// std::cout<<"3\n";
	// modal_matrix_->ProjectVector(u_,q_);
	// integrator_base_dense_->SetState(q_);
	// std::cout<<"1:\n";
	// for(int i=0;i<30;++i)
	// {
	// 	std::cout<<u_[i]<<",";
	// }
	// modal_matrix_->ProjectVector(u_,q_);
	// std::cout<<"before:q_\n";
	// for(int i=0;i<r_;++i)
	// std::cout<<q_[i]<<",";
	// std::cout<<"\n";
	if(time_step_counter_ < total_steps_)
	{
		int code=integrator_base_dense_->DoTimestep();
		std::cout<<"."<<code;fflush(NULL);
		++time_step_counter_;
	}
	// std::cout<<"2:\n";
	// for(int i=0;i<30;++i)
	// {
	// 	std::cout<<u_[i]<<",";
	// }
	memcpy(q_,integrator_base_dense_->Getq(),sizeof(double)*r_);
	// std::cout<<"after:q_\n";
	// for(int i=0;i<r_;++i)
	// std::cout<<q_[i]<<",";
	// std::cout<<"\n";
	// std::cout<<"3:\n";
	// for(int i=0;i<30;++i)
	// {
	// 	std::cout<<u_[i]<<",";
	// }
	// std::cout<<"4\n";
	// for(int i=0;i<r_;++i)
	// std::cout<<q_[i]<<",";
	// std::cout<<"\n";
	// // getchar();
	// memset(u_,0.0,sizeof(double)*3*simulation_vertices_num_);
	// for(int i=0;i<3*simulation_vertices_num_;++i)
	// {
	// 	for(int j=0;j<r_;++j)
	// 	{
	// 		u_[i]+=U_[3*simulation_vertices_num_*j+i]*q_[j];
	// 	}
	// }
	// std::cout<<"4:\n";
	// for(int i=0;i<30;++i)
	// {
	// 	std::cout<<u_[i]<<",";
	// }
	modal_matrix_->AssembleVector(q_,u_);
	// std::cout<<"5:\n";
	// for(int i=0;i<30;++i)
	// 	std::cout<<u_[i]<<",";
	// modal_matrix_->ProjectVector(u_,q_);
	// std::cout<<"after-after:q_\n";
	// for(int i=0;i<r_;++i)
	// std::cout<<q_[i]<<",";
	// std::cout<<"\n";
	// setu(u_);
	// std::cout<<"5\n";
	if(time_step_counter_==(total_steps_-1))
	{
		std::cout<<"Total target configuration solver time: "<<total_target_time_<<" s.\n";
	}
	// getchar();
}
void RealTimeExampleBasedDeformer::reducedspaceSimulationWithoutConstraints()
{
	memset(fq_,0.0,sizeof(double)*r_);
	// static int count1=1;
	// static double temp_coef=1.005;
	// std::cout<<"a\n";
	// getchar();
	rigid_->ResetWrenches();
	rigid_external_force_[0]=rigid_external_force_[1]=rigid_external_force_[2]=0.0;
	memset(example_f_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	// memset(f_ext_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	// //compute collision force
	// if(collide_vert_num_>0)
	// {
	// 	for(int i=0;i<simulation_vertices_num_;++i)
	// 	{
	// 		for(int j=0;j<3;++j)
	// 			f_ext_[3*i+j]+=f_ext_[3*i+j]*vert_mass_[i]*2.0;
	// 	}
	// }

    // static int count_num=1;
	if(time_step_counter_<force_loads_num_)
	{
		// std::cout<<"External forces read from the text input file.\n";
		for(int i=0;i<3*simulation_vertices_num_;++i)
			f_ext_[i]+=force_loads_[ELT(3*simulation_vertices_num_,i,time_step_counter_)];

	}
	//transfer u_ in global to local frame
	rigid_->GetRotation(R_);
	Mat3d R_matrix(R_);
	Mat3d invR=inv(R_matrix);
	rigid_->GetPosition(&t_[0],&t_[1],&t_[2]);
	// for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	// {
	// 	Vec3d global_x(0.0);
	// 	global_x[0]=(*simulation_mesh_->getVertex(i))[0]+u_[3*i];
	// 	global_x[1]=(*simulation_mesh_->getVertex(i))[1]+u_[3*i+1];
	// 	global_x[2]=(*simulation_mesh_->getVertex(i))[2]+u_[3*i+2];
	// 	Vec3d local_x=invR*(global_x-t_);
	// 	local_u_[3*i]=local_x[0]-local_reference_[3*i];
	// 	local_u_[3*i+1]=local_x[1]-local_reference_[3*i+1];
	// 	local_u_[3*i+2]=local_x[2]-local_reference_[3*i+2];
	//
	// }
	if(enable_example_simulation_)
	{
		PerformanceCounter project_counter;
		project_counter.StartCounter();
		projectOnEigenFunctions1(local_u_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
								interpolate_eigenfunction_num_,object_current_eigencoefs_);
		project_counter.StopCounter();
		total_projection_time_+=project_counter.GetElapsedTime();
		// projectOnEigenFunctions(simulation_mesh_,local_u_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
		// 						interpolate_eigenfunction_num_,object_current_eigencoefs_);
		// projectOnEigenFunctions(simulation_mesh_,u_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
		// 						interpolate_eigenfunction_num_,object_current_eigencoefs_);
		// for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 	std::cout<<object_current_eigencoefs_[i]<<",";
		// if(count1==1)
		// {
		// 	for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 		for(int j=0;j<3;++j)
		// 		{
		// 			object_eigencoefs_[i][j]=example_eigencoefs_[0][i][j];
		// 			object_current_eigencoefs_[i][j]=example_eigencoefs_[0][i][j];
		// 		}
		// 	count1++;
		// }
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			object_current_eigencoefs_[i]=object_current_eigencoefs_[i]-object_eigencoefs_[i]+example_eigencoefs_[0][i];
		memcpy(target_eigencoefs_,object_current_eigencoefs_,sizeof(Vec3d)*interpolate_eigenfunction_num_);
		// PerformanceCounter target_counter;
		// target_counter.StartCounter();
		PerformanceCounter target_counter;
		target_counter.StartCounter();
		projectOnExampleManifold(object_current_eigencoefs_,target_eigencoefs_);

		target_counter.StopCounter();
		total_target_time_+=target_counter.GetElapsedTime();
		// target_counter.StopCounter();
		// total_target_time_+=target_counter.GetElapsedTime();
		//method 1:reconstruct by dis_coefs-begin:
		//--------------------------------------------------------------------------------------------
		//----------------used for reduced elastic force-reconstruct target configuration first--------
		// reconstructFromEigenCoefs(target_eigencoefs_,target_reconstruction_);
		//---------------------------------------------------------------------------------------------
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			target_eigencoefs_[i]=object_current_eigencoefs_[i]-target_eigencoefs_[i];
		// for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 	std::cout<<target_eigencoefs_[i]<<"--";
		// std::cout<<"\n";
		PerformanceCounter reconstruct_counter;
		reconstruct_counter.StartCounter();
		reconstructFromEigenCoefs(target_eigencoefs_,example_guided_deformation_);
		reconstruct_counter.StopCounter();
		total_reconstruction_time_+=reconstruct_counter.GetElapsedTime();
		//method 1:reconstruct by dis_coefs-end.
		// method 2: construct local frame-begin:
		// reconstructFromEigenCoefs(target_eigencoefs_,2);//reference+target_deformation=target configuration
		// for(int i=0;i<simulation_mesh_->getNumVertices();++i)
		// {
		// 	example_guided_deformation_[3*i+0]=(*simulation_mesh_->getVertex(i))[0]+u_[3*i+0]-target_reconstruction_[3*i+0];
		// 	example_guided_deformation_[3*i+1]=(*simulation_mesh_->getVertex(i))[1]+u_[3*i+1]-target_reconstruction_[3*i+1];
		// 	example_guided_deformation_[3*i+2]=(*simulation_mesh_->getVertex(i))[2]+u_[3*i+2]-target_reconstruction_[3*i+2];
		// 	// example_guided_deformation_[3*i+0]=(*simulation_mesh_->getVertex(i))[0]+u_[3*i+0]-example_guided_deformation_[3*i+0];
		// 	// example_guided_deformation_[3*i+1]=(*simulation_mesh_->getVertex(i))[1]+u_[3*i+1]-example_guided_deformation_[3*i+1];
		// 	// example_guided_deformation_[3*i+2]=(*simulation_mesh_->getVertex(i))[2]+u_[3*i+2]-example_guided_deformation_[3*i+2];
		// 	// example_guided_deformation_[3*i+0]=local_u_[3*i]+local_reference_[3*i]-example_guided_deformation_[3*i+0];
		// 	// example_guided_deformation_[3*i+1]=local_u_[3*i+1]+local_reference_[3*i+1]-example_guided_deformation_[3*i+1];
		// 	// example_guided_deformation_[3*i+2]=local_u_[3*i+2]+local_reference_[3*i+2]-example_guided_deformation_[3*i+2];
		// 	// example_guided_deformation_[3*i+0]=u_[3*i+0]-target_reconstruction_[3*i+0];
		// 	// example_guided_deformation_[3*i+1]=u_[3*i+1]-target_reconstruction_[3*i+1];
		// 	// example_guided_deformation_[3*i+2]=u_[3*i+2]-target_reconstruction_[3*i+2];
		// }
		//method 2: construct local frame-end.
		// linear force
		for(int i=0;i<3*simulation_vertices_num_;++i)
			example_f_[i]=(-1.0)*example_stiffness_scale_*example_guided_deformation_[i];
		//reduced elastic force
		// for(int i=0;i<simulation_mesh_->getNumVertices();++i)
		// {
		// 	target_reconstruction_[3*i+0]=(*simulation_mesh_->getVertex(i))[0]+u_[3*i+0]-example_guided_deformation_[3*i+0];
		// 	target_reconstruction_[3*i+1]=(*simulation_mesh_->getVertex(i))[1]+u_[3*i+1]-example_guided_deformation_[3*i+1];
		// 	target_reconstruction_[3*i+2]=(*simulation_mesh_->getVertex(i))[2]+u_[3*i+2]-example_guided_deformation_[3*i+2];
		// }
		//----------------------------------------------------------------------------
		//--------------reduced stvk elastic force-------------------------------------
		// modal_matrix_->ProjectVector(example_guided_deformation_,example_based_q_);
		//
		// if(material_mode_==REDUCED_STVK)
		// 	// reduced_stvk_cubature_force_model_->computeReducedElasticInternalForce(example_guided_deformation_,fq_,target_reconstruction_);//target_deformation_);
		// 	reduced_stvk_cubature_force_model_->computeReducedInternalForce(example_based_q_,fq_);
		// for(int i=0;i<r_;++i)
		// 	fq_[i]*=-example_stiffness_scale_;
		// modal_matrix_->AssembleVector(fq_,example_f_);
		//-----------------------------------------------------------------------------
		// if(material_mode_==REDUCED_STVK)
		// 	reduced_stvk_cubature_force_model_->computeReducedElasticInternalForce(example_guided_deformation_,example_based_fq_,target_reconstruction_);//target_deformation_);
		// // if(material_mode_==REDUCED_NEOHOOKEAN)
		// // 	reduced_neoHookean_force_model_->computeReducedElasticInternalForce(example_based_q_,example_based_fq_,example_guided_deformation_);//target_deformation_);
		// for(int i=0;i<r_;++i)
		// 	fq_[i]=(-1.0)*example_based_fq_[i];
		// modal_matrix_->AssembleVector(fq_,example_f_);

		// if(time_step_counter_%timer_sample_interval_==0&&time_step_counter_>=timer_sample_interval_)
        // {
		// 	// std::cout<<timer_sample_interval_<<"------\n";
		// 	std::cout<<"Total project solver time: "<<total_projection_time_/timer_sample_interval_<<" s.\n";
        //     std::cout<<"Total target solver time: "<<total_target_time_/timer_sample_interval_<<" s.\n";
		// 	std::cout<<"Total reconstruct time: "<<total_reconstruction_time_/timer_sample_interval_<<" s.\n";
        //     total_projection_time_=total_target_time_=total_reconstruction_time_= 0;
        //     // getchar();
        // }
	}
	else
	{
		for(int i=0;i<3*simulation_vertices_num_;++i)
			example_f_[i]*=damping_example_stiffness_;
	}
	// memcpy(fq_,example_f_,sizeof(double)*r_);
	// memcpy(qaccel_,integrator_base_->Getqaccel(),sizeof(double)*r_);
	// modal_matrix_->AssembleVector(qaccel_,Uqaccel_);
	//compute non-inertial force

	double total_torquex,total_torquey,total_torquez;
	double alpha=1.0e-2;
	#pragma omp parallel for
	for(int i=0;i<simulation_vertices_num_;++i)
	{
		double torquex,torquey,torquez;
		Vec3d pos;
		pos[0]=(*simulation_mesh_->getVertex(i))[0]+u_[3*i+0];
		pos[1]=(*simulation_mesh_->getVertex(i))[1]+u_[3*i+1];
		pos[2]=(*simulation_mesh_->getVertex(i))[2]+u_[3*i+2];

		Vec3d vert_f;

		for(int j=0;j<3;++j)
		{
			// example_f_[3*i+j]-=vert_mass_[i]*Uqaccel_[j];
			if(enable_example_simulation_)
				// rigid_external_force_[j]+=example_f_[3*i+j];
				vert_f[j]=example_f_[3*i+j];
			rigid_external_force_[j]+=f_ext_[3*i+j];//-vert_mass_[i]*Uqaccel_[j];
		}
		if(enable_example_simulation_)
		{
			Vec3d temp=R_matrix*vert_f;
			for(int j=0;j<3;++j)
				rigid_external_force_[j]+=temp[j];
		}

		if(enable_example_simulation_)
		{
			rigid_->ComputeTorque(pos[0],pos[1],pos[2],example_f_[3*i+0],example_f_[3*i+1],
								example_f_[3*i+2],&torquex,&torquey,&torquez);
			total_torquex+=torquex*torque_coef_;
			total_torquey+=torquey*torque_coef_;
			total_torquez+=torquez*torque_coef_;
		}
		rigid_->ComputeTorque(pos[0],pos[1],pos[2],f_ext_[3*i+0],f_ext_[3*i+1],
								f_ext_[3*i+2],&torquex,&torquey,&torquez);
		// rigid_->ComputeTorque(pos[0],pos[1],pos[2],rigid_external_force_[0],rigid_external_force_[1],
		// 						rigid_external_force_[2],&torquex,&torquey,&torquez);
		total_torquex+=torquex*alpha;
		total_torquey+=torquey*alpha;
		total_torquez+=torquez*alpha;
	}
	for(int i=0;i<3;++i)
	{
		if(collide_vert_num_>col_limited_num_)
		{
			// rigid_external_force_[i]=rigid_external_force_[i]*100.0/collide_vert_num_;
			total_torquex=total_torquex*col_limited_num_/collide_vert_num_;
			total_torquey=total_torquey*col_limited_num_/collide_vert_num_;
			total_torquez=total_torquez*col_limited_num_/collide_vert_num_;
		}
	}
	rigid_->SetExternalForce(rigid_external_force_[0],rigid_external_force_[1],rigid_external_force_[2]);
	if(enable_example_simulation_)
		rigid_->SetExternalTorque(total_torquex,total_torquey,total_torquez);
	// static int count_force=1;

	// rigid_->GetVelocity(&linear_velocity_[0],&linear_velocity_[1],&linear_velocity_[2]);
	// std::cout<<"before-lin-vel:"<<linear_velocity_[0]<<","<<linear_velocity_[1]<<","<<linear_velocity_[2]<<"\n";
	// rigid_->GetAngularVelocity(&angular_velocity_[0],&angular_velocity_[1],&angular_velocity_[2]);
	// std::cout<<"before-ang-vel:"<<angular_velocity_[0]<<","<<angular_velocity_[1]<<","<<angular_velocity_[2]<<"\n";
	// rigid_->GetAngularMomentum(&am_[0],&am_[1],&am_[2]);
	// std::cout<<"before-ang-momentum:"<<am_[0]<<","<<am_[1]<<","<<am_[2]<<"\n";
	// std::cout<<"R_matrix:"<<R_matrix<<"\n";
	if(add_gravity_)
		rigid_->AddExternalForce(0.0,(-1.0)*total_mass_*gravity_,0.0);
	rigid_->EulerStep(time_step_);
	rigid_->GetAngularVelocity(&new_angular_velocity_[0],&new_angular_velocity_[1],&new_angular_velocity_[2]);
	rigid_->GetVelocity(&new_linear_velocity_[0],&new_linear_velocity_[1],&new_linear_velocity_[2]);
	// std::cout<<"after-lin-vel:"<<new_linear_velocity_[0]<<","<<new_linear_velocity_[1]<<","<<new_linear_velocity_[2]<<"\n";
	// std::cout<<"after-ang-vel:"<<new_angular_velocity_[0]<<","<<new_angular_velocity_[1]<<","<<new_angular_velocity_[2]<<"\n";
	// rigid_->GetAngularMomentum(&am_[0],&am_[1],&am_[2]);
	// std::cout<<"after-ang-momentum:"<<am_[0]<<","<<am_[1]<<","<<am_[2]<<"\n";
	// getchar();
	rigid_->GetRotation(R_);
	Mat3d new_R_matrix(R_);
	// std::cout<<"new_R_matrix:"<<new_R_matrix<<"\n";
	Mat3d new_invR=inv(new_R_matrix);
	rigid_->GetPosition(&t_[0],&t_[1],&t_[2]);
	Vec3d linear_vel_grad=(new_linear_velocity_-linear_velocity_)/time_step_;
	Vec3d angular_vel_grad=(new_angular_velocity_-angular_velocity_)/time_step_;
	angular_velocity_=new_angular_velocity_;
	linear_velocity_=new_linear_velocity_;
	if(add_gravity_)
	{
		for(int i=0;i<simulation_mesh_->getNumVertices();++i)
		{
			f_ext_[3*i+1]+=(-1.0)*vert_mass_[i]*gravity_;
		}
	}

	memcpy(qvel_,integrator_base_->Getqvel(),sizeof(double)*r_);
	modal_matrix_->AssembleVector(qvel_,Uqvel_);
	static int count=1;
	#pragma omp parallel for
	for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	{
		// Vec3d vert_global_dis;
		Vec3d vert_global_pos,vert_local_pos,vert_global_u,vert_local_u;
		Vec3d vert_f;
		Vec3d local_dis_grad;
		Vec3d Uqaccel;
		for(int j=0;j<3;++j)
		{
			vert_global_pos[j]=(*simulation_mesh_->getVertex(i))[j]+u_[3*i+j];
			vert_global_u[j]=u_[3*i+j];
			vert_f[j]=f_ext_[3*i+j];
			local_dis_grad[j]=Uqvel_[3*i+j];//-linear_velocity_[j];
		}
		// gv[0]+=f_ext_[3*i+0]-vert_mass_[i]*Uqaccel_[3*i+0];
		// gv[1]+=f_ext_[3*i+1]-vert_mass_[i]*Uqaccel_[3*i+1];
		// gv[2]+=f_ext_[3*i+2]-vert_mass_[i]*Uqaccel_[3*i+2];
		vert_local_pos=new_R_matrix*(vert_global_pos-t_);
		// vert_local_u=new_invR*vert_global_u;
		// for(int j=0;j<3;++j)
		// {
		// 	local_u_[3*i+j]=vert_local_u[j];
		// }
		//coriolis force-not done yet
		Vec3d temp1=cross(angular_vel_grad,vert_local_pos);//omega x local_frame_position
		// local_dis_grad-=temp1;
		Vec3d temp2=cross(angular_velocity_,local_dis_grad);//omega x dot(U_i,qvel)
		Vec3d f_cor=(-2.0)*vert_mass_[i]*temp2;
		//inertial force:done
		Vec3d f_ine=(-1.0)*vert_mass_[i]*linear_vel_grad;
		//euler force:
		Vec3d f_eul=(-1.0)*vert_mass_[i]*temp1;
		//centrifugal force:done
		Vec3d wxq=cross(new_angular_velocity_,vert_local_pos);//omega x local_frame_pos
		Vec3d temp4=cross(new_angular_velocity_,wxq); //omega x (omega x local_frame_pos)
		Vec3d f_cen=(-1.0)*vert_mass_[i]*temp4;
		Vec3d f_fic=f_cor+f_ine+f_eul+f_cen;
		Vec3d temp5=new_invR*vert_f; //R^T*f_ext
		//example_f is computed in the non-inertia frame
		f_ext_[3*i]=temp5[0]+f_fic[0]+example_f_[3*i];
		f_ext_[3*i+1]=temp5[1]+f_fic[1]+example_f_[3*i+1];
		f_ext_[3*i+2]=temp5[2]+f_fic[2]+example_f_[3*i+2];
		// Vec3d temp=cross(angular_velocity_,local_x);
		// Uqvel_[3*i+0]=collide_vel_[3*i+0]-linear_velocity_[0]-wxq[0];
		// Uqvel_[3*i+1]=collide_vel_[3*i+1]-linear_velocity_[1]-wxq[1];
		// Uqvel_[3*i+2]=collide_vel_[3*i+2]-linear_velocity_[2]-wxq[2];
		// if(i==0)
		// {
		// 	std::cout<<"f_cor:"<<f_cor<<"\n";
		// 	std::cout<<"f_ine:"<<f_ine<<"\n";
		// 	std::cout<<"f_eul:"<<f_eul<<"\n";
		// 	std::cout<<"f_cen:"<<f_cen<<"\n";
		// 	std::cout<<"f_fic:"<<f_fic<<"\n";
		// }
	}
	modal_matrix_->ProjectVector(f_ext_,fq_);
	integrator_base_dense_->SetExternalForces(fq_);
	// int a=integrator_base_dense_->SetState(q_,qvel_);
	if(time_step_counter_ < total_steps_)
	{
		// PerformanceCounter counter2;
		// counter2.StartCounter();
		int code=integrator_base_dense_->DoTimestep();
			std::cout<<".";fflush(NULL);
		// counter2.StopCounter();
		// std::cout<<"integrator DoTimestep:"<<counter2.GetElapsedTime()<<"\n";
		++time_step_counter_;
	}
	//local q_n+1 in non-inertial frame, project it in inertial frame
	memcpy(q_,integrator_base_dense_->Getq(),sizeof(double)*r_);
	modal_matrix_->AssembleVector(q_,local_u_);
	// memcpy(qvel_,integrator_base_dense_->Getqvel(),sizeof(double)*r_);
	// modal_matrix_->AssembleVector(qvel_,Uqvel_);

	// //q_ is in non-inertial frame at time_step n;
	// //compute pos in non-inertial frame, project the position on inertial frame--global_x
	// //compute the displacement u_ in inertial frame, u_=current_defomration_configuration-global_reference_time_0 for rendering

	//convert local_u_ from non-inertia frame at time n to non-inertia frame at time n+1, store still in local_u_
	//get u_ in inertia frame, u_ contains both the deformation displacement in non=inertia and rigid displacement in inertia,
	//we adopt project position from non-inertia position to inertia world space, then compute u_;
	// modal_matrix_->AssembleVector(qvel_,Uqvel_);
	#pragma omp parallel for
	for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	{
		Vec3d local_x(0.0),local_v(0.0);
		local_x[0]=local_u_[3*i]+local_reference_[3*i];
		local_x[1]=local_u_[3*i+1]+local_reference_[3*i+1];
		local_x[2]=local_u_[3*i+2]+local_reference_[3*i+2];
		Vec3d global_x=new_R_matrix*local_x+t_;//R_n
		// //u_ used for rendering,generated by inertia force and non-inertial force
		u_[3*i]=global_x[0]-(*simulation_mesh_->getVertex(i))[0];
		u_[3*i+1]=global_x[1]-(*simulation_mesh_->getVertex(i))[1];
		u_[3*i+2]=global_x[2]-(*simulation_mesh_->getVertex(i))[2];
		// Vec3d temp=cross(new_angular_velocity_,local_x);//omega x dot(U_i,qvel)
		//
		// for(int j=0;j<3;++j)
		// 	local_v[j]=Uqvel_[3*i+j];
		// for(int j=0;j<3;++j)
		// 	vel_[3*i+j]=new_linear_velocity_[j]+temp[j]+Uqvel_[3*i+j];
	}
	if(time_step_counter_==(total_steps_-1))
	{
		std::cout<<"Total target configuration solver time: "<<total_target_time_<<" s.\n";
	}
}

void RealTimeExampleBasedDeformer::reconstructFromEigenCoefs(Vec3d *target_eigencoefs,double *vert_pos)
{
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// std::cout<<reconstruct_eigenfunction_num_<<"\n";
	double scale_factor;
	memset(vert_pos,0.0,sizeof(double)*3*simulation_vertices_num_);
	// memset(vert_pos,0.0,sizeof(double)*3*examples_[1]->getNumVertices());
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	target_eigencoefs[i]=example_eigencoefs_[2][i];
	// 	std::cout<<target_eigencoefs_[i]<<",";
	// 	// std::cout<<example_eigencoefs_[0][i]<<",";
	// }
	// std::cout<<"\n";
	// std::cout<<"ex0:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
		// target_eigencoefs[i]=object_eigencoefs_[i];
		// std::cout<<target_eigencoefs_[i]<<",";
		// std::cout<<example_eigencoefs_[0][i]<<",";
	// }
	// std::cout<<"\n";
	// std::cout<<"ex1:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	target_eigencoefs[i]=example_eigencoefs_[1][i];
	// 	// std::cout<<target_eigencoefs_[i]<<",";
	// 	std::cout<<example_eigencoefs_[1][i]<<",";
	// }
	// std::cout<<"\n";
	for(unsigned int vert_idx=0;vert_idx<simulation_vertices_num_;++vert_idx)
	{
		for(unsigned int i=0;i<interpolate_eigenfunction_num_;++i)
		{
			scale_factor=1.0/object_eigenvalues_[i];
			vert_pos[3*vert_idx]+=target_eigencoefs[i][0]*object_eigenfunctions_[i][vert_idx]*scale_factor;
			vert_pos[3*vert_idx+1]+=target_eigencoefs[i][1]*object_eigenfunctions_[i][vert_idx]*scale_factor;
			vert_pos[3*vert_idx+2]+=target_eigencoefs[i][2]*object_eigenfunctions_[i][vert_idx]*scale_factor;

		}
		// for(unsigned int i=interpolate_eigenfunction_num_;i<reconstruct_eigenfunction_num_;++i)
		// {
		// 	scale_factor=1.0/object_eigenvalues_[i];
		// 	vert_pos[3*vert_idx]+=object_eigencoefs_[i][0]*object_eigenfunctions_[i][vert_idx]*scale_factor;
		// 	vert_pos[3*vert_idx+1]+=object_eigencoefs_[i][1]*object_eigenfunctions_[i][vert_idx]*scale_factor;
		// 	vert_pos[3*vert_idx+2]+=object_eigencoefs_[i][2]*object_eigenfunctions_[i][vert_idx]*scale_factor;
		// }
	}

	// rigid_->GetRotation(R_);
	// Mat3d new_R_matrix(R_);
	// for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	// {
	// 	Vec3d local_x(0.0);
	// 	local_x[0]=vert_pos[3*i];
	// 	local_x[1]=vert_pos[3*i+1];
	// 	local_x[2]=vert_pos[3*i+2];
	// 	Vec3d global_x=new_R_matrix*local_x;//R_n
	// 	// //u_ used for rendering,generated by inertia force and non-inertial force
	//
	// 	vert_pos[3*i]=local_x[0];
	// 	vert_pos[3*i+1]=local_x[1];
	// 	vert_pos[3*i+2]=local_x[2];
	// }
	// saveReconstructMesh(vert_pos);
	// getchar();
}
void RealTimeExampleBasedDeformer::reconstructFromEigenCoefs(Vec3d *target_eigencoefs,int flag)
{
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// std::cout<<reconstruct_eigenfunction_num_<<"\n";
	double scale_factor;
	double *vert_pos=new double[3*simulation_vertices_num_];
	memset(vert_pos,0.0,sizeof(double)*3*simulation_vertices_num_);
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	//
	// 	std::cout<<target_eigencoefs_[i]<<",";
	// }
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	target_eigencoefs[i]=example_eigencoefs_[0][i];
	// 	std::cout<<target_eigencoefs_[i]<<",";
	// }
	for(unsigned int vert_idx=0;vert_idx<simulation_vertices_num_;++vert_idx)
	{
		for(unsigned int i=0;i<interpolate_eigenfunction_num_;++i)
		{
			scale_factor=1.0/object_eigenvalues_[i];
			vert_pos[3*vert_idx]+=target_eigencoefs[i][0]*object_eigenfunctions_[i][vert_idx]*scale_factor;
			vert_pos[3*vert_idx+1]+=target_eigencoefs[i][1]*object_eigenfunctions_[i][vert_idx]*scale_factor;
			vert_pos[3*vert_idx+2]+=target_eigencoefs[i][2]*object_eigenfunctions_[i][vert_idx]*scale_factor;
		}
		// for(unsigned int i=interpolate_eigenfunction_num_;i<reconstruct_eigenfunction_num_;++i)
		// {
		// 	scale_factor=1.0/object_eigenvalues_[i];
		// 	vert_pos[3*vert_idx]+=object_eigencoefs_[i][0]*object_eigenfunctions_[i][vert_idx]*scale_factor;
		// 	vert_pos[3*vert_idx+1]+=object_eigencoefs_[i][1]*object_eigenfunctions_[i][vert_idx]*scale_factor;
		// 	vert_pos[3*vert_idx+2]+=object_eigencoefs_[i][2]*object_eigenfunctions_[i][vert_idx]*scale_factor;
		// }
	}
	std::cout<<"before generate new detail vector:\n";
	if(flag==1)
	{//initial:generate local detail vector
		generateLocalDetailVector(vert_pos);
		// std::cout<<"flag==1\n";
	}
	else if(flag==2)
	{
		generateNewDetailVector(vert_pos);
		// std::cout<<"flag==2\n";
		// saveReconstructMesh(vert_pos);
		// getchar();
	}
	delete[] vert_pos;
}
void RealTimeExampleBasedDeformer::saveReconstructMesh1(double *new_pos)
{
	static int id=0;
	++id;
	std::string id_str;
	std::stringstream stream;
	stream<<id;
	stream>>id_str;
	std::string file_name="target/a"+id_str+".smesh";
	// std::string file_name="obj_obj-detail.smesh";
	std::ofstream output_file(file_name.c_str());
	if(!output_file)
	{
		std::cout<<"Error: failed to open "<<file_name<<".\n";
		return;
	}
	output_file<<"*VERTICES"<<std::endl;
	output_file<<examples_[1]->getNumVertices()<<" 3 0 0"<<std::endl;
	for(unsigned int i=0;i<examples_[1]->getNumVertices();++i)
	{
		Vec3d pos;
		// pos[0]=(*simulation_mesh_->getVertex(i))[0]-new_pos[3*i+0];
		// pos[1]=(*simulation_mesh_->getVertex(i))[1]-new_pos[3*i+1];
		// pos[2]=(*simulation_mesh_->getVertex(i))[2]-new_pos[3*i+2];
		// pos[0]=(*examples_[1]->getVertex(i))[0];
		// pos[1]=(*examples_[1]->getVertex(i))[1];
		// pos[2]=(*examples_[1]->getVertex(i))[2];
		// pos[0]=target_reconstruction_[3*i+0];
		// pos[1]=target_reconstruction_[3*i+1];
		// pos[2]=target_reconstruction_[3*i+2];
		// pos[0]=(*simulation_mesh_->getVertex(i))[0]+new_pos[3*i+0];
		// pos[1]=(*simulation_mesh_->getVertex(i))[1]+new_pos[3*i+1];
		// pos[2]=(*simulation_mesh_->getVertex(i))[2]+new_pos[3*i+2];
		pos[0]=new_pos[3*i+0];
		pos[1]=new_pos[3*i+1];
		pos[2]=new_pos[3*i+2];

		// pos[0]=(*simulation_mesh_->getVertex(i))[0]+u_[3*i+0];
		// pos[1]=(*simulation_mesh_->getVertex(i))[1]+u_[3*i+1];
		// pos[2]=(*simulation_mesh_->getVertex(i))[2]+u_[3*i+2];
		output_file<<i+1<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<std::endl;
	}
	output_file<<"*ELEMENTS"<<std::endl;
	output_file<<"TET"<<std::endl;
	output_file<<examples_[1]->getNumElements()<<" 4 0"<<std::endl;
	for(unsigned int i=0;i<examples_[1]->getNumElements();++i)
	{
		output_file<<i+1;
		for(int j=0;j<4;++j)
		{
			int vertID=examples_[1]->getVertexIndex(i,j)+1;
			output_file<<" "<<vertID;
		}
		output_file<<std::endl;
	}
}
void RealTimeExampleBasedDeformer::saveReconstructMesh(double *new_pos)
{
	static int id=0;
	++id;
	std::string id_str;
	std::stringstream stream;
	stream<<id;
	stream>>id_str;
	std::string file_name="target/a"+id_str+".smesh";
	// std::string file_name="obj_obj-detail.smesh";
	std::ofstream output_file(file_name.c_str());
	if(!output_file)
	{
		std::cout<<"Error: failed to open "<<file_name<<".\n";
		return;
	}
	output_file<<"*VERTICES"<<std::endl;
	output_file<<simulation_mesh_->getNumVertices()<<" 3 0 0"<<std::endl;
	for(unsigned int i=0;i<simulation_mesh_->getNumVertices();++i)
	{
		Vec3d pos;
		// pos[0]=(*simulation_mesh_->getVertex(i))[0]-new_pos[3*i+0];
		// pos[1]=(*simulation_mesh_->getVertex(i))[1]-new_pos[3*i+1];
		// pos[2]=(*simulation_mesh_->getVertex(i))[2]-new_pos[3*i+2];
		// pos[0]=(*examples_[1]->getVertex(i))[0];
		// pos[1]=(*examples_[1]->getVertex(i))[1];
		// pos[2]=(*examples_[1]->getVertex(i))[2];
		// pos[0]=target_reconstruction_[3*i+0];
		// pos[1]=target_reconstruction_[3*i+1];
		// pos[2]=target_reconstruction_[3*i+2];
		// pos[0]=(*simulation_mesh_->getVertex(i))[0]+new_pos[3*i+0];
		// pos[1]=(*simulation_mesh_->getVertex(i))[1]+new_pos[3*i+1];
		// pos[2]=(*simulation_mesh_->getVertex(i))[2]+new_pos[3*i+2];
		pos[0]=new_pos[3*i+0];
		pos[1]=new_pos[3*i+1];
		pos[2]=new_pos[3*i+2];

		// pos[0]=(*simulation_mesh_->getVertex(i))[0]+u_[3*i+0];
		// pos[1]=(*simulation_mesh_->getVertex(i))[1]+u_[3*i+1];
		// pos[2]=(*simulation_mesh_->getVertex(i))[2]+u_[3*i+2];
		output_file<<i+1<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<std::endl;
	}
	output_file<<"*ELEMENTS"<<std::endl;
	output_file<<"TET"<<std::endl;
	output_file<<simulation_mesh_->getNumElements()<<" 4 0"<<std::endl;
	for(unsigned int i=0;i<simulation_mesh_->getNumElements();++i)
	{
		output_file<<i+1;
		for(int j=0;j<4;++j)
		{
			int vertID=simulation_mesh_->getVertexIndex(i,j)+1;
			output_file<<" "<<vertID;
		}
		output_file<<std::endl;
	}
}

//initial detail vector generation:
void RealTimeExampleBasedDeformer::generateLocalDetailVector(const double *vert_pos)
{
	if(!simulation_mesh_)
	{
		std::cout<<"simulation mesh unloaded!\n";
		exit(0);
	}
	//construct local frame for each vertex
	for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	{//global=R*local+pos
		Vec3d pos1,pos2,pos3;
		pos1[0]=vert_pos[3*i+0];
		pos1[1]=vert_pos[3*i+1];
		pos1[2]=vert_pos[3*i+2];
		pos2[0]=vert_pos[3*vert_vertex1[i]+0];
		pos2[1]=vert_pos[3*vert_vertex1[i]+1];
		pos2[2]=vert_pos[3*vert_vertex1[i]+2];
		pos3[0]=vert_pos[3*vert_vertex2[i]+0];
		pos3[1]=vert_pos[3*vert_vertex2[i]+1];
		pos3[2]=vert_pos[3*vert_vertex2[i]+2];
		Vec3d local_frame1=pos2-pos1;
		Vec3d local_frame2=cross(pos2-pos1,pos3-pos1);
		Vec3d local_frame3=cross(local_frame1,local_frame2);
		local_frame1.normalize();
		local_frame2.normalize();
		local_frame3.normalize();
		if(i==0)
			std::cout<<"local_frame1"<<local_frame1<<"\n";
		if(i==0)
			std::cout<<"local_frame2"<<local_frame2<<"\n";
		if(i==0)
			std::cout<<"local_frame3"<<local_frame3<<"\n";
		//compute R from global to local frame
		Mat3d R(1.0);
		for(int j=0;j<3;++j)
		{
			R[j][0]=local_frame1[j];
			R[j][1]=local_frame2[j];
			R[j][2]=local_frame3[j];
		}
		Vec3d global_vert_detail_vector;
		for(int j=0;j<3;++j)
		{
			global_vert_detail_vector[j]=(*simulation_mesh_->getVertex(i))[j]-vert_pos[3*i+j];
			// global_vert_detail_vector[j]=examples_deformation0_[3*i+j]+(*simulation_mesh_->getVertex(i))[j]-vert_pos[3*i+j];
		}
		Mat3d invR=inv(R);
		//global detail vector to local detail vector
		for(int j=0;j<3;++j)
		{
			local_detail_vector_[3*i+j]=invR[j][0]*(global_vert_detail_vector[0])+invR[j][1]*(global_vert_detail_vector[1])
											+invR[j][2]*(global_vert_detail_vector[2]);

		}
		//project local detail vector to non-inertia frame
		// if(!with_constrains_)
		// {
		// 	rigid_->get
		// }
	}
}
void RealTimeExampleBasedDeformer::generateNewDetailVector(const double *vert_pos)
{
	if(!simulation_mesh_)
	{
		std::cout<<"simulation mesh unloaded!\n";
		exit(0);
	}
	//construct local frame for each vertex
	for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	{
		Vec3d pos1,pos2,pos3;
		pos1[0]=vert_pos[3*i+0];
		pos1[1]=vert_pos[3*i+1];
		pos1[2]=vert_pos[3*i+2];
		pos2[0]=vert_pos[3*vert_vertex1[i]+0];
		pos2[1]=vert_pos[3*vert_vertex1[i]+1];
		pos2[2]=vert_pos[3*vert_vertex1[i]+2];
		pos3[0]=vert_pos[3*vert_vertex2[i]+0];
		pos3[1]=vert_pos[3*vert_vertex2[i]+1];
		pos3[2]=vert_pos[3*vert_vertex2[i]+2];
		Vec3d local_frame1=pos2-pos1;
		Vec3d local_frame2=cross(pos2-pos1,pos3-pos1);
		Vec3d local_frame3=cross(local_frame1,local_frame2);
		local_frame1.normalize();
		local_frame2.normalize();
		local_frame3.normalize();
		if(i==0)
			std::cout<<"local_frame1"<<local_frame1<<"\n";
		if(i==0)
			std::cout<<"local_frame2"<<local_frame2<<"\n";
		if(i==0)
			std::cout<<"local_frame3"<<local_frame3<<"\n";
		//compute R from global to local frame
		Mat3d R(1.0);
		for(int j=0;j<3;++j)
		{
			R[j][0]=local_frame1[j];
			R[j][1]=local_frame2[j];
			R[j][2]=local_frame3[j];

		}
		if(i==0)
			std::cout<<R<<"\n";
		for(int j=0;j<3;++j)
		{
			target_reconstruction_[3*i+j]=R[j][0]*local_detail_vector_[3*i+0]+R[j][1]*local_detail_vector_[3*i+1]
											+R[j][2]*local_detail_vector_[3*i+2]+pos1[j];
		}

		Mat3d invR=inv(R);
		// std::cout<<"invR:"<<invR<<"\n";
		//update local detail vector for current step,global detail vector to local detail vector
		for(int j=0;j<3;++j)
		{
			local_detail_vector_[3*i+j]=invR[j][0]*(target_reconstruction_[3*i]-vert_pos[3*i])
										+invR[j][1]*(target_reconstruction_[3*i+1]-vert_pos[3*i+1])
										+invR[j][2]*(target_reconstruction_[3*i+2]-vert_pos[3*i+2]);
			// std::cout<<local_detail_vector_[3*i+j]<<",";

		}
	}
	// saveReconstructMesh();
	// getchar();
}

void RealTimeExampleBasedDeformer::projectOnEigenFunctions1(double *displacement, double *vertex_volume,
                                                           double **eigenfunctions, double *eigenvalues,
														   unsigned int eigenfunction_num,Vec3d *eigencoefs)
{
	double scale_factor;
	// std::cout<<"---------show-----------\n";
	// std::cout<<"eigenfun_num:"<<eigenfunction_num<<"\n";
	for(unsigned int i=0;i<eigenfunction_num;++i)
	{
		for(unsigned int j=0;j<3;++j)
		{
			eigencoefs[i][j]=0.0;
		}
		scale_factor=eigenvalues[i];
		for(unsigned int vert_idx=0;vert_idx<simulation_mesh_->getNumVertices();++vert_idx)
		{
			Vec3d vert_pos;
			for(unsigned int dim=0;dim<3;++dim)
			{
				vert_pos[dim]=local_reference_[3*vert_idx+dim]+displacement[3*vert_idx+dim];
				eigencoefs[i][dim]+=vert_pos[dim]*eigenfunctions[i][vert_idx]*vertex_volume[vert_idx]*scale_factor;
			}
		}
	}
}
//the result eigencoefs is a shape not a displacement on LB space
void RealTimeExampleBasedDeformer::projectOnEigenFunctions(VolumetricMesh *mesh, double *displacement, double *vertex_volume,
                                                           double **eigenfunctions, double *eigenvalues, unsigned int eigenfunction_num,
                                                           Vec3d *eigencoefs,int affected_vertices_num,int *affected_vertices)
{
	//local example mode is not written yet
    int vert_num=mesh->getNumVertices();
	double scale_factor;
	// std::cout<<"---------show-----------\n";
	// std::cout<<"eigenfun_num:"<<eigenfunction_num<<"\n";
	for(unsigned int i=0;i<eigenfunction_num;++i)
	{
		for(unsigned int j=0;j<3;++j)
		{
			eigencoefs[i][j]=0.0;
		}
		scale_factor=eigenvalues[i];
		//local example mode, only project local regions onto eigenfunctions
		if(affected_vertices_num>0)
		{
			for(unsigned int j=0;j<affected_vertices_num;++j)
			{
				int vert_idx=affected_vertices[j]-1;
				Vec3d vert_pos;
				for(unsigned int dim=0;dim<3;++dim)
				{
					vert_pos[dim]=(*mesh->getVertex(vert_idx))[dim]+displacement[3*vert_idx+dim];
					eigencoefs[i][dim]+=vert_pos[dim]*eigenfunctions[i][vert_idx]*vertex_volume[vert_idx]*scale_factor;
				}
			}
			continue;
		}
		for(unsigned int vert_idx=0;vert_idx<vert_num;++vert_idx)
		{
			Vec3d vert_pos;
			for(unsigned int dim=0;dim<3;++dim)
			{
				vert_pos[dim]=(*mesh->getVertex(vert_idx))[dim]+displacement[3*vert_idx+dim];
				eigencoefs[i][dim]+=vert_pos[dim]*eigenfunctions[i][vert_idx]*vertex_volume[vert_idx]*scale_factor;
			}
		}
	}
}

void RealTimeExampleBasedDeformer::projectOnExampleManifold(Vec3d *input_eigencoefs, Vec3d *projected_eigencoefs)
{
	//suppose the target configuration is: w(0)*example_eigencoefs_(0)+w(1)*example_eigencoefs_(1)+...+w(example_num_)*example_eigencoefs_(exmaple_num_-1)
	//w(0)+w(1)+...+w(example_num_-1)=1, 0<=w(i)<=1, i=0,1,...,example_num_-1
	//compute the weights by minimizing 1/2*dist(target,input)*dist(target_input), least square minimization with constraints
	//solve this minimization using optimization solver from ALGLIB
	// PerformanceCounter counter1;
	// counter1.StartCounter();
	double *temp_buffer = new double[example_num_+1];
	memset(temp_buffer,0.0,sizeof(double)*(example_num_+1));
	real_1d_array target_weights;
	target_weights.setcontent(example_num_,temp_buffer); //initial guess
	//prepare boundarty constraints:0<=w(i)<=1, i=0...example_num-1;
	real_1d_array lower_boundary;
	lower_boundary.setcontent(example_num_,temp_buffer);
	real_1d_array upper_boundary;
	for(int i=0;i<example_num_+1;++i)
		temp_buffer[i]=1.0;

	// std::cout<<"....target_weights:"<<target_weights[0]<<"......"<<target_weights[1]<<"\n";
	upper_boundary.setcontent(example_num_,temp_buffer);
	//prepare the equality constraints: w(0)+w(1)+...+w(example_num_-1)=1
	real_2d_array equality_constraint;
	equality_constraint.setcontent(1,example_num_+1,temp_buffer);
	integer_1d_array constraint_type="[0]";
	//stop conditions
	double epsg=integrator_epsilon_;
	double epsf=integrator_epsilon_;
	double epsx=integrator_epsilon_;
	ae_int_t max_iterations= 1000;
	//start solving
	minbleicstate state;
	minbleicreport report;
	minbleiccreate(target_weights,state);
	minbleicsetbc(state,lower_boundary,upper_boundary);
	minbleicsetlc(state,equality_constraint,constraint_type);
	minbleicsetcond(state,epsg,epsf,epsx,max_iterations);
	minbleicoptimize(state,evaluateObjectiveAndGradient1,NULL,(void*)input_eigencoefs);//the last parameter is passed to evaluateObjectiveAndGradient
	minbleicresults(state,target_weights,report);
	//print some information about the solve if needed
	cout<<"Step1 converged in "<<report.iterationscount<<" iterations. Terminatin type of the solver: "<<report.terminationtype<<"\n";

	//bias the weights towards example other than the rest pose if necessary
	//so that the examples take effect even when natrual deformation is far from Example
	//only trigger this while the object is defroming away from initial configuration
	//we want the object to recover, thus no bias is needed when deforming back to rest configuration
	//note: all the parameters here need tuning for different demos

	bool is_deform_away=(last_initial_weight_-target_weights[0]>=1.0e-5);
	double lower_bound=0.05,upper_bound=1.01;
	bool initial_weight_in_range = (target_weights[0]>lower_bound)&&(target_weights[0]<upper_bound);
	enable_eigen_weight_control_=true;
	// std::cout<<"\n";
	std::cout<<"....last_initial_weight_:"<<last_initial_weight_<<"......target_weights[0]:"<<target_weights[0]<<"\n";
	std::cout<<"....is_deform_away:"<<is_deform_away<<"......initial_weight_in_range:"<<initial_weight_in_range<<"\n";
	if(enable_eigen_weight_control_&&example_num_>1
		// &&initial_weight_in_range
    	)
	{
		// double bias_factor=0.8; //different biases are used for different simulationMesh
		//example 1 is by default the rest pose
		double other_weight_delta=(1-example_bias_)*target_weights[0]/(example_num_-1);
		target_weights[0]*=example_bias_;
		for(int i=1; i<example_num_;++i)
			target_weights[i]+=other_weight_delta;
		last_initial_weight_=target_weights[0];
	}

	//for some trick-begin
	// if(target_weights[0]<0.1)
	// 	target_weights[0]=0.1;
	//for some trick-end
	//pure_example_linear_interpolation=true;
	//compute the target coefficients with the weights
	for(int i=0;i<interpolate_eigenfunction_num_;++i)
	{
		projected_eigencoefs[i][0]=projected_eigencoefs[i][1]=projected_eigencoefs[i][2]=0.0;
		for(int j=0;j<example_num_;++j)
		{
			projected_eigencoefs[i][0]+=target_weights[j]*example_eigencoefs_[j][i][0];
			projected_eigencoefs[i][1]+=target_weights[j]*example_eigencoefs_[j][i][1];
			projected_eigencoefs[i][2]+=target_weights[j]*example_eigencoefs_[j][i][2];
		}
	}

	// std::cout<<"ex0-eigencoefs:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<example_eigencoefs_[0][i]<<",";
	// std::cout<<"\n";
	//
	// std::cout<<"ex1-eigencoefs:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<example_eigencoefs_[1][i]<<",";
	// std::cout<<"\n";
	//
	// std::cout<<"ex0-eigencoefs:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<projected_eigencoefs[i]<<",";
	// std::cout<<"\n";

	free(temp_buffer);
	if(simulation_mode_==FULLSPACE)
		pure_example_linear_interpolation_=true;
	else
		pure_example_linear_interpolation_=false;

	// counter1.StopCounter();
	// std::cout<<"step 1:"<<counter1.GetElapsedTime()<<"\n";
	// counter1.StopCounter();
	// std::cout<<"elapseTime:"<<counter1.GetElapsedTime()<<"\n";
	// target_weights[0]=0.5;
	// target_weights[1]=0.5;
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	for(int j=0;j<3;++j)
	// 		projected_eigencoefs[i][j]=example_eigencoefs_[1][i][j];
	// PerformanceCounter counter2;
	// counter2.StartCounter();
	//alternative step2

	// PerformanceCounter counter2;
	// counter2.StartCounter();
	// std::cout<<"---a--\n";
	if(!pure_example_linear_interpolation_)
	{
		//compute for evaluateObjectiveAndGradient4-begin:
		// double **temp_coef=new double*[example_num_];
		// for(unsigned int idx=0;idx<example_num_;++idx)
		// {
		// 	temp_coef=new double [interpolate_eigenfunction_num_];
		// 	for(unsigned int i=0;i<interpolate_eigenfunction_num_;++i)
		// 		temp_coef[idx][i]=example_eigencoefs_[idx][i]-object_eigencoefs_[i];
		// 	ex_E_mat_
		// }
		//end
		double *input_data=new double[example_num_];
		for(int i=0;i<example_num_;++i)
			input_data[i]=target_weights[i];
		temp_buffer=new double[3*interpolate_eigenfunction_num_];
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
		{
			for(int j=0;j<3;++j)//current configuration or linear interpolation can be used as initial gueses
				temp_buffer[3*i+j]=projected_eigencoefs[i][j]; //input_eigen_coefs[i][j]
		}
		real_1d_array new_eigencoefs;
		mincgstate new_state;
		mincgreport new_report;
		max_iterations=1000;
		epsg*=100.0;
		epsf*=100.0;
		epsx*=100.0;
		new_eigencoefs.setcontent(3*interpolate_eigenfunction_num_,temp_buffer);
		mincgcreate(new_eigencoefs,new_state);
		mincgsetcond(new_state,epsg,epsf,epsx,max_iterations);
		// std::cout<<"----1----\n";
		mincgoptimize(new_state,evaluateObjectiveAndGradient2,NULL,(void*)input_data);
		// std::cout<<"----2----\n";
		// getchar();
		mincgresults(new_state,new_eigencoefs,new_report);
		//print some information about the solver if needed
		cout<<"Step2 converged in "<<new_report.iterationscount<<" iterations. Termination type of the solver: "<<new_report.terminationtype<<".\n";

		delete[] temp_buffer;
		delete[] input_data;
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			for(int j=0;j<3;++j)
				projected_eigencoefs[i][j]=new_eigencoefs[3*i+j];
				// projected_eigencoefs[i][j]=example_eigencoefs_[1][i][j];
	}
	// counter2.StopCounter();
	// std::cout<<"step 2:"<<counter2.GetElapsedTime()<<"\n";

	cout<<"The target_weights in example manifold is:";
	for(int i=0;i<example_num_;++i)
		cout<<target_weights[i]<<",";
	cout<<"\n";

	// counter2.StopCounter();
	// std::cout<<"step2:elapseTime:"<<counter2.GetElapsedTime()<<"\n";
		// counter2.StopCounter();
	// std::cout<<"target:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<projected_eigencoefs[i]<<"\n";
	// std::cout<<"\n";
	// getchar();

}

//helper function, used in projectOnExampleManifold()
//evaluate the value of objective function and its gradient, ALGLIB needs it
//'ptr' is passed as input_eigencoefs in projectOnExampleManifold()
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient1(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
{
	// std::cout<<"--grad1-begin\n";
	RealTimeExampleBasedDeformer* active_instance=RealTimeExampleBasedDeformer::activeInstance();
	assert(active_instance);
	Vec3d *input_eigencoefs=(Vec3d*)ptr;
	Vec3d *temp_target_eigencoefs=new Vec3d[active_instance->interpolate_eigenfunction_num_];
	//evaluate the objective function
	func=0.0;
	for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
	{
		temp_target_eigencoefs[i][0]=temp_target_eigencoefs[i][1]=temp_target_eigencoefs[i][2]=0.0;
		for(int j=0;j<active_instance->example_num_;++j)
		{
			temp_target_eigencoefs[i][0]+=x[j]*active_instance->example_eigencoefs_[j][i][0];
			temp_target_eigencoefs[i][1]+=x[j]*active_instance->example_eigencoefs_[j][i][1];
			temp_target_eigencoefs[i][2]+=x[j]*active_instance->example_eigencoefs_[j][i][2];
		}
		func+=(temp_target_eigencoefs[i][0]-input_eigencoefs[i][0])*(temp_target_eigencoefs[i][0]-input_eigencoefs[i][0]);
		func+=(temp_target_eigencoefs[i][1]-input_eigencoefs[i][1])*(temp_target_eigencoefs[i][1]-input_eigencoefs[i][1]);
		func+=(temp_target_eigencoefs[i][2]-input_eigencoefs[i][2])*(temp_target_eigencoefs[i][2]-input_eigencoefs[i][2]);
		// std::cout<<"temp_target_eigencoefs:"<<temp_target_eigencoefs[i][0]<<","<<temp_target_eigencoefs[i][1]<<","<<temp_target_eigencoefs[i][2]<<"\n";
		// std::cout<<"input_eigencoefs:"<<input_eigencoefs[i][0]<<","<<input_eigencoefs[i][1]<<","<<input_eigencoefs[i][2]<<"\n";
	}
	func*=0.5;
	// std::cout<<"---------------func="<<func<<"\n";
	//analytically evaluate the gradient of the objective function
	for(int i=0;i<active_instance->example_num_;++i)
	{
		grad[i]=0.0;
		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			grad[i]+=(temp_target_eigencoefs[j][0]-input_eigencoefs[j][0])*active_instance->example_eigencoefs_[i][j][0];
			grad[i]+=(temp_target_eigencoefs[j][1]-input_eigencoefs[j][1])*active_instance->example_eigencoefs_[i][j][1];
			grad[i]+=(temp_target_eigencoefs[j][2]-input_eigencoefs[j][2])*active_instance->example_eigencoefs_[i][j][2];
		}
	}
	delete[] temp_target_eigencoefs;
	// std::cout<<"grad1--end.\n";
}
//step2, minimization of wieghted deformation energy to the examples
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient2(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr)
{
	// std::cout<<"---3---\n";
	PerformanceCounter counter;
	counter.StartCounter();
	RealTimeExampleBasedDeformer* active_instance=RealTimeExampleBasedDeformer::activeInstance();
	assert(active_instance);
	func=0.0;
	for(int j=0;j<3*active_instance->interpolate_eigenfunction_num_;++j)
		grad[j]=0.0;
	double *target_weights=(double*)ptr;
	Vec3d *temp_eigencoefs=new Vec3d[active_instance->interpolate_eigenfunction_num_];
	double energy=0.0;
	//temp code
	// double *temp_buffer = (double*)ptr;
	// double *temp_buffer =new double[active_instance->example_num_+1];
	// memset(temp_buffer,0.0,sizeof(double)*(active_instance->example_num_+1));
	//
	// real_1d_array target_weights;
	// target_weights.setcontent(active_instance->example_num_,temp_buffer); //initial guess
	// //prepare boundarty constraints:0<=w(i)<=1, i=0...example_num-1;
	// real_1d_array lower_boundary;
	// lower_boundary.setcontent(active_instance->example_num_,temp_buffer);
	// real_1d_array upper_boundary;
	// for(int i=0;i<active_instance->example_num_+1;++i)
	// 	temp_buffer[i]=1.0;
	//
	// // std::cout<<"....target_weights:"<<target_weights[0]<<"......"<<target_weights[1]<<"\n";
	// upper_boundary.setcontent(active_instance->example_num_,temp_buffer);
	// //prepare the equality constraints: w(0)+w(1)+...+w(example_num_-1)=1
	// real_2d_array equality_constraint;
	// equality_constraint.setcontent(1,active_instance->example_num_+1,temp_buffer);
	// integer_1d_array constraint_type="[0]";
	// //stop conditions
	// double epsg=active_instance->integrator_epsilon_;
	// double epsf=active_instance->integrator_epsilon_;
	// double epsx=active_instance->integrator_epsilon_;
	// ae_int_t max_iterations= 1000;
	// //start solving
	// Vec3d *input_eigencoefs;
	// for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
	// 	for(int j=0;j<3;++j)
	// 	{
	// 		input_eigencoefs[i][j]=x[3*i+j];
	// 		std::cout<<input_eigencoefs[i][j]<<",";
	// 	}
	// std::cout<<"\n";
	// // getchar();
	// minbleicstate state;
	// minbleicreport report;
	// minbleiccreate(target_weights,state);
	// minbleicsetbc(state,lower_boundary,upper_boundary);
	// minbleicsetlc(state,equality_constraint,constraint_type);
	// minbleicsetcond(state,epsg,epsf,epsx,max_iterations);
	// minbleicoptimize(state,evaluateObjectiveAndGradient1,NULL,(void*)input_eigencoefs);//the last parameter is passed to evaluateObjectiveAndGradient
	// minbleicresults(state,target_weights,report);
	// std::cout<<"target_weights:"<<target_weights[0]<<"......"<<target_weights[1]<<"......"<<target_weights[2]<<"\n";
	//temp code
	for(int i=0;i<active_instance->example_num_;++i)
	{
		// std::cout<<"----a---\n";
	// 	//compute displacement in LB space for each example
		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
			{
				temp_eigencoefs[j][k]=x[3*j+k]-active_instance->example_eigencoefs_[i][j][k];
			}
		}
		// std::cout<<"----b---\n";
		//method 1: compute energy and energy gradient using LB basis and its cubature optimization ELEMENTS
		// active_instance->computeReducedEnergy(temp_eigencoefs,energy);
		// func+=target_weights[i]*energy;
		// double *energy_grad=new double[3*active_instance->interpolate_eigenfunction_num_];
		// memset(energy_grad,0.0,sizeof(double)*3*active_instance->interpolate_eigenfunction_num_);
		// active_instance->computeReducedInternalForce(temp_eigencoefs,energy_grad);
		//
		// for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		// {
		// 	for(int k=0;k<3;++k)
		// 		grad[3*j+k]+=target_weights[i]*energy_grad[3*j+k];
		// }
		//method 1-end

		//method 2: compute energy and internal force using modal basis and cubature optimization elements
		// active_instance->reconstructFromEigenCoefs(temp_eigencoefs,active_instance->example_guided_deformation_);
		// active_instance->modal_matrix_->ProjectVector(active_instance->example_guided_deformation_,active_instance->temp_q_);
		// active_instance->reduced_stvk_cubature_force_model_->computeReducedEnergy(active_instance->temp_q_,energy);
		// func+=target_weights[i]*energy;
		// active_instance->reduced_stvk_cubature_force_model_->computeReducedInternalForce(active_instance->temp_q_,active_instance->temp_grad_);
		// active_instance->modal_matrix_->AssembleVector(active_instance->temp_grad_,active_instance->example_guided_deformation_);

		// active_instance->projectOnEigenFunctions1(active_instance->example_guided_deformation_,
		// 						active_instance->object_vertex_volume_,active_instance->object_eigenfunctions_,
		// 						active_instance->object_eigenvalues_,
		// 						active_instance->interpolate_eigenfunction_num_,active_instance->energy_grad_);

		// for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		// {
		// 	for(int k=0;k<3;++k)
		// 		grad[3*j+k]+=target_weights[i]*active_instance->energy_grad_[j][k];
		// }
		//method 2-end

		//method 3: improve method 2, LB space to reduced space
		// static double counter1_time=0.0,counter2_time=0.0,counter3_time=0.0;		//
		// PerformanceCounter counter1;
		// counter1.StartCounter();
		memset(active_instance->temp_q_,0.0,sizeof(double)*active_instance->r_);
		//test:LB TO REDUCED BASIS
		for(unsigned int j=0;j<active_instance->r_;++j)
		{
			for(unsigned int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
			{
				for(unsigned int k=0;k<3;++k)
					active_instance->temp_q_[j]+=active_instance->LB_to_reducespace_[j][3*i+k]*temp_eigencoefs[i][k];
			}
		}
		// std::cout<<"----c---\n";
		// counter1.StopCounter();
		// counter1_time=counter1.GetElapsedTime();
		// std::cout<<"project q:"<<counter1_time<<"\n";

		// PerformanceCounter counter2;
		// counter2.StartCounter();
		active_instance->reduced_stvk_cubature_force_model_->computeReducedEnergy(active_instance->temp_q_,energy);
		func+=target_weights[i]*energy;
		active_instance->reduced_stvk_cubature_force_model_->computeReducedInternalForce(active_instance->temp_q_,active_instance->temp_grad_);
		// counter2.StopCounter();
		// counter2_time+=counter2.GetElapsedTime();
		// std::cout<<"compute energy and grad:"<<counter2_time<<"\n";

		// PerformanceCounter counter3;
		// counter3.StartCounter();
		for(unsigned int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			active_instance->energy_grad_[j][0]=active_instance->energy_grad_[j][1]=active_instance->energy_grad_[j][2]=0.0;
			for(unsigned int i=0;i<active_instance->r_;++i)
			{
				for(unsigned int k=0;k<3;++k)
				{
					active_instance->energy_grad_[j][k]+=active_instance->temp_grad_[i]*active_instance->reducespace_to_LB_[i][3*j+k];
				}
			}
		}
		// std::cout<<"----d---\n";
		// counter3.StopCounter();
		// counter3_time+=counter3.GetElapsedTime();
		// std::cout<<"project to LB:"<<counter3_time<<"\n";

		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
				grad[3*j+k]+=target_weights[i]*active_instance->energy_grad_[j][k];
		}
		// std::cout<<"----e---\n";
		// getchar();
		//method 3-end
	}
	delete[] temp_eigencoefs;
	// std::cout<<"---4---\n";
	// getchar();
}
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient3(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr)
{
	// PerformanceCounter counter;
	// counter.StartCounter();
	RealTimeExampleBasedDeformer* active_instance=RealTimeExampleBasedDeformer::activeInstance();
	assert(active_instance);
	func=0.0;
	for(int j=0;j<3*active_instance->interpolate_eigenfunction_num_;++j)
		grad[j]=0.0;
	double *target_weights=(double*)ptr;
	Vec3d *temp_eigencoefs=new Vec3d[active_instance->interpolate_eigenfunction_num_];
	double energy=0.0;
	for(int i=0;i<active_instance->example_num_;++i)
	{
	// 	//compute displacement in LB space for each example
		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
			{
				temp_eigencoefs[j][k]=x[3*j+k]-active_instance->example_eigencoefs_[i][j][k];
				std::cout<<temp_eigencoefs[j][k]<<",";
				// temp_eigencoefs[j][k]=active_instance->object_eigencoefs_[j][k]-active_instance->example_eigencoefs_[i][j][k];
			}
		}
		//active_instance->computeF(temp_eigencoefs);

		//method 1: compute energy and energy gradient using LB basis and its cubature optimization ELEMENTS
		active_instance->computeEnergy(temp_eigencoefs,energy);
		std::cout<<"energy:"<<energy<<"\n";
		func+=target_weights[i]*energy;
		double *energy_grad=new double[3*active_instance->interpolate_eigenfunction_num_];
		memset(energy_grad,0.0,sizeof(double)*3*active_instance->interpolate_eigenfunction_num_);
		active_instance->computeEnergyGrad(temp_eigencoefs,energy_grad);
		//method 1-end

		std::cout<<"energy-grad:\n";
		for(int i=0;i<3*active_instance->interpolate_eigenfunction_num_;++i)
			std::cout<<energy_grad[i]<<",";
		std::cout<<"\n";
		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
				grad[3*j+k]+=energy_grad[3*j+k];
		}
		std::cout<<"grad:\n";
		for(int i=0;i<3*active_instance->interpolate_eigenfunction_num_;++i)
			std::cout<<grad[i]<<",";
		std::cout<<"\n";
		delete[] energy_grad;
	}
	delete[] temp_eigencoefs;
}
/*
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient4(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr)
{
}*/
void RealTimeExampleBasedDeformer::preAllocateLocalFrameCorrespondingVertices()
{
	std::cout<<"preAllocateLocalFrameCorrespondingVertices-begin.\n";
	if(!simulation_mesh_)
	{
		std::cout<<"simulation mesh unloaded!\n";
		exit(0);
	}
	for(int vertID=0;vertID<simulation_mesh_->getNumVertices();++vertID)
	{
		vert_vertex1[vertID]=-1;
		vert_vertex2[vertID]=-1;
	}
	double *vert_quality=new double[simulation_mesh_->getNumVertices()];
	memset(vert_quality,0.0,sizeof(double)*simulation_mesh_->getNumVertices());
	for(int eleID=0;eleID<simulation_mesh_->getNumElements();++eleID)
	{
		int *global_vert=new int[4];
		memset(global_vert,0,sizeof(int)*4);
		double *cell_length=new double[6];
		memset(cell_length,0.0,sizeof(double)*6);
		for(int i=0;i<4;++i)
		{
			global_vert[i]=simulation_mesh_->getVertexIndex(eleID,i);
		}
		for(int i=0;i<3;++i)
		{
			for(int j=i+1;j<4;++j)
			{
				Vec3d vert1=(*simulation_mesh_->getVertex(global_vert[i]));
				Vec3d vert2=(*simulation_mesh_->getVertex(global_vert[j]));
				if(i==0)
				{
					cell_length[j-1]=sqrt((vert1[0]-vert2[0])*(vert1[0]-vert2[0])+(vert1[1]-vert2[1])*(vert1[1]-vert2[1])+(vert1[2]-vert2[2])*(vert1[2]-vert2[2]));
					// std::cout<<"cell_length"<<j-1<<":"<<cell_length[j-1]<<"\n";
				}
				else
				{
					cell_length[i+j]=sqrt((vert1[0]-vert2[0])*(vert1[0]-vert2[0])+(vert1[1]-vert2[1])*(vert1[1]-vert2[1])+(vert1[2]-vert2[2])*(vert1[2]-vert2[2]));
					// std::cout<<"cell_length"<<i+j<<":"<<cell_length[i+j]<<"\n";
				}

			}
		}
		// for(int i=0;i<6;++i)
		// {
		// 	std::cout<<cell_length[i]<<"\n";
		// }
		// getchar();

		double l_rms_temp=0.0,l_harm=0.0,l_rms=0.0;
		for(int i=0;i<6;++i)
		{
			l_rms_temp+=cell_length[i]*cell_length[i];
			l_harm+=1.0/cell_length[i];
		}
		l_rms=(l_rms_temp/6.0)*(l_rms_temp/6.0);
		l_harm=6.0/l_harm;
		double tet_quality=6.0*sqrt(2.0)*(simulation_mesh_->getElementVolume(eleID))*l_harm/l_rms;
		for(int i=0;i<4;++i)
		{
			if(vert_quality[global_vert[i]]<tet_quality)
			{
				vert_vertex1[global_vert[i]]=global_vert[(i+1)%4];
				vert_vertex2[global_vert[i]]=global_vert[(i+2)%4];
				vert_quality[global_vert[i]]=tet_quality;
			}
			// if(vert_vertex1[global_vert[i]]==-1)
				// vert_vertex1[global_vert[i]]=global_vert[(i+1)%4];
			// if(vert_vertex2[global_vert[i]]==-1)
				// vert_vertex2[global_vert[i]]=global_vert[(i+2)%4];
		}
		// std::cout<<"ele:"<<eleID<<","<<simulation_mesh_->getElementVolume(eleID)<<"\n";
		delete[] cell_length;
		delete[] global_vert;
	}
	// for(int i=0;i<simulation_mesh_->getNumVertices();++i)
	// 	std::cout<<vert_quality[i]<<"\n";
	delete[] vert_quality;
	std::cout<<"preAllocateLocalFrameCorrespondingVertices-end.\n";
	// for(int vertID=0;vertID<simulation_mesh_->getNumVertices();++vertID)
	// {
	// 	std::cout<<"vertID:"<<vertID<<","<<vert_vertex1[vertID]<<","<<vert_vertex2[vertID]<<"\n";
	// 	if((vert_vertex1[vertID]==-1)||(vert_vertex2[vertID]==-1))
	// 		getchar();
	// }
	// getchar();
}
bool RealTimeExampleBasedDeformer::loadMassmatrix(const std::string &file_name)
{
	//get the mass matrix
    SparseMatrixOutline *mass_matrix_outline;
    try
    {
        //3 is expansion flag to indicate this is a mass matrix and does 3x3 identify block expansion
        mass_matrix_outline=new SparseMatrixOutline(file_name.c_str(),1);
    }
    catch(int exception_code)
    {
        std::cout<<"Error: loading mass matrix failed.\n";
        exit(1);
    }
    mass_matrix_=new SparseMatrix(mass_matrix_outline);
    delete(mass_matrix_outline);
	//compute for reduced simulation
	if(simulation_mode_==REDUCEDSPACE)
	{
		reduced_mass_matrix_=new double[r_*r_];
	    memset(reduced_mass_matrix_,0.0,sizeof(double)*r_*r_);
		// for(int i=0;i<r_;++i)
		// 	reduced_mass_matrix_[ELT(r_,i,i)] = 1.0;
	    //U_ is column major, the column number is r_
	    U_ = new double[3*simulation_vertices_num_*r_];
	    memset(U_,0.0,sizeof(double)*3*simulation_vertices_num_*r_);
	    for(int j=0;j<r_;++j)
	    {
	        for(int i=0;i<3*simulation_vertices_num_;++i)
	        {
	            U_[3*simulation_vertices_num_*j+i]=reduced_basis_[j][i];
	        }
	    }
		//reduced_mass_matrix_ is identity()
	    mass_matrix_->ConjugateMatrix(U_,r_,reduced_mass_matrix_);
		// for(int i=0; i<r_; i++)
	    // 	reduced_mass_matrix_[ELT(r_,i,i)] = 1.0;
		// for(int i=0;i<r_;++i)
		// {
		// 	for(int j=0;j<r_;++j)
		// 		std::cout<<reduced_mass_matrix_[i*r_+j]<<",";
		// 	std::cout<<"\n";
		// }
	    modal_matrix_ = new ModalMatrix(simulation_vertices_num_,r_,U_);
	}

	return true;
}
//
bool RealTimeExampleBasedDeformer::loadInertiaTensor(const std::string &file_name)
{
    std::fstream input_file(file_name.c_str());
	if(!input_file)
	{
		std::cout<<"Error: failed open "<<file_name<<" .\n";
        return false;
	}
	int num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		double temp_value;
		input_file>>temp_value;
		initial_inertia_tensor_[num]=temp_value;
		if(num>=9)
			break;
		num++;
	}
	input_file.close();
	return true;
}
bool RealTimeExampleBasedDeformer::loadSimulationMesh(const std::string &file_name)
{
    std::cout<<"Loading simulation mesh: "<<file_name<<std::endl;
	simulation_mesh_=VolumetricMeshLoader::load(file_name.c_str());
	if(simulation_mesh_==NULL)
	{
		std::cout<<"Error: unable to load the simulation mesh from "<<file_name<<std::endl;
		return false;
	}
	simulation_vertices_num_=simulation_mesh_->getNumVertices();
	object_volume_=simulation_mesh_->getVolume();
	object_vertex_volume_=new double[simulation_mesh_->getNumVertices()];
	memset(object_vertex_volume_,0.0,sizeof(double)*simulation_mesh_->getNumVertices());
	for(unsigned int ele_idx=0;ele_idx<simulation_mesh_->getNumElements();++ele_idx)
	{
		for(unsigned int i=0;i<simulation_mesh_->getNumElementVertices();++i)
		{
			unsigned int global_vert_idx=simulation_mesh_->getVertexIndex(ele_idx,i);
			object_vertex_volume_[global_vert_idx]+=simulation_mesh_->getElementVolume(ele_idx)/4.0;
		}
	}
	//material setting, we use homogeneous materials,get the first element material
	VolumetricMesh::Material * material = simulation_mesh_->getElementMaterial(0);
	VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
	if (eNuMaterial == NULL)
	{
		std::cout<<"Error: mesh does not consist of E, nu materials.\n";
		exit(0);
	}
	lamda_ = eNuMaterial->getLambda();
	mu_ = eNuMaterial->getMu();
    return true;
}

bool RealTimeExampleBasedDeformer::loadExamples(const std::string &file_name_prefix, unsigned int example_num)
{
   	std::cout<<"loading examples: "<<std::endl;
	example_num_=example_num;
	examples_=new VolumetricMesh *[example_num];
	example_volume_=new double[example_num];
	example_vertex_volume_=new double *[example_num];
    std::stringstream stream;
	// std::cout<<"\n";
	for(unsigned int i=0;i<example_num;++i)
	{
        std::string example_num_str="",file_name="";
		stream.str("");
		stream.clear();
		stream<<i;
		stream>>example_num_str;
		file_name=file_name_prefix+example_num_str+string(".veg");
		examples_[i]=VolumetricMeshLoader::load(file_name.c_str());
        if(examples_[i]==NULL)
        {
            std::cout<<"Error: unable to load the example mesh from "<<file_name<<std::endl;
            return false;
        }
		example_volume_[i]=examples_[i]->getVolume();
		std::cout<<i<<":"<<example_volume_[i]<<std::endl;
		example_vertex_volume_[i]=new double[examples_[i]->getNumVertices()];
		memset(example_vertex_volume_[i],0.0,sizeof(double)*examples_[i]->getNumVertices());
		for(unsigned int ele_idx=0;ele_idx<examples_[i]->getNumElements();++ele_idx)
		{
			for(unsigned int j=0;j<examples_[i]->getNumElementVertices();++j)
			{
				unsigned int global_idx=examples_[i]->getVertexIndex(ele_idx,j);

				example_vertex_volume_[i][global_idx]+=examples_[i]->getElementVolume(ele_idx)/4.0;
			}
		}
		// std::cout<<"ex"<<i<<":\n";
		// for(int j=0;j<examples_[i]->getNumVertices();++j)
		// 	std::cout<<example_vertex_volume_[i][j]<<",";
		// 	std::cout<<"\n";
	}
	std::cout<<"loading examples succeed.\n";
	//temp
	//load temp target smesh
	// std::string temp_file_name="3.veg";
	// temp_mesh_=VolumetricMeshLoader::load(temp_file_name.c_str());
	// if(temp_mesh_==NULL)
	// {
	// 	std::cout<<"Error: unable to load temp mesh from"<<temp_mesh_<<std::endl;
	// 	return false;
	// }
	// std::cout<<"Num:"<<temp_mesh_->getNumVertices()<<"\n";
	// examples_deformation0_=new double[3*temp_mesh_->getNumVertices()];
	// for(int i=0;i<temp_mesh_->getNumVertices();++i)
	// {
	// 	examples_deformation0_[3*i+0]=(*temp_mesh_->getVertex(i))[0]-(*simulation_mesh_->getVertex(i))[0];
	// 	examples_deformation0_[3*i+1]=(*temp_mesh_->getVertex(i))[1]-(*simulation_mesh_->getVertex(i))[1];
	// 	examples_deformation0_[3*i+2]=(*temp_mesh_->getVertex(i))[2]-(*simulation_mesh_->getVertex(i))[2];
	// }
	// std::string temp_file_name1="4.veg";
	// temp_mesh_=VolumetricMeshLoader::load(temp_file_name1.c_str());
	// if(temp_mesh_==NULL)
	// {
	// 	std::cout<<"Error: unable to load temp mesh from"<<temp_mesh_<<std::endl;
	// 	return false;
	// }
	// std::cout<<"Num:"<<temp_mesh_->getNumVertices()<<"\n";
	// examples_deformation_=new double[3*temp_mesh_->getNumVertices()];
	// for(int i=0;i<temp_mesh_->getNumVertices();++i)
	// {
	// 	examples_deformation_[3*i+0]=(*temp_mesh_->getVertex(i))[0]-(*simulation_mesh_->getVertex(i))[0];
	// 	examples_deformation_[3*i+1]=(*temp_mesh_->getVertex(i))[1]-(*simulation_mesh_->getVertex(i))[1];
	// 	examples_deformation_[3*i+2]=(*temp_mesh_->getVertex(i))[2]-(*simulation_mesh_->getVertex(i))[2];
	// }
	//temp
    return true;

}

bool RealTimeExampleBasedDeformer::loadVisualMesh(const std::string &file_name)
{
    visual_mesh_=new SceneObjectDeformable(file_name.c_str());
    if(visual_mesh_==NULL)
    {
        std::cout<<"Error: unable to load visual mesh from "<<file_name<<std::endl;
        return false;
    }
	std::cout<<"visual_mesh_:"<<visual_mesh_->GetNumVertices();
    return true;
}

//load reduced basis:first line is row_num=vertex_num,second line is col_num=basis num;
//basis are 3*vertex_num*basis_num
bool RealTimeExampleBasedDeformer::loadReducedBasis(const std::string &file_name)
{
	std::cout<<"file_name:"<<file_name<<"\n";
	std::cout<<"load reduced basis:\n";
	if(simulation_mesh_==NULL)
	{
		std::cout<<"The simulation tet mesh is null.\n";
        return false;
	}
    std::fstream input_file(file_name.c_str());
	if(!input_file)
	{
		std::cout<<"Error: failed open "<<file_name<<" .\n";
        return false;
	}
    std::string temp_str;

	while(std::getline(input_file,temp_str))
	{
		if(temp_str.compare(0,7,string("*values"))==0)
			break;
	}
	//read file, save as reduced_basis_[reduced_idx][3*vert_idx+dim]
	std::string reduced_basis_value_num_str;
	std::getline(input_file,reduced_basis_value_num_str);
	r_=atoi(reduced_basis_value_num_str.c_str());
	reduced_basis_values_=new double[r_];

	unsigned int num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		double temp_value;
		input_file>>temp_value;
		reduced_basis_values_[num]=temp_value;
		num++;
		if(num>r_-1)
			break;
	}
	while(std::getline(input_file,temp_str))
	{
		if(temp_str.compare(0,8,string("*vectors"))==0)
			break;
	}
	unsigned int reduced_row_num,reduced_col_num;
	std::getline(input_file,temp_str);
	reduced_row_num=atoi(temp_str.c_str());
	std::getline(input_file,temp_str);
	reduced_col_num=atoi(temp_str.c_str());
	std::cout<<reduced_row_num<<","<<reduced_col_num<<"\n";
	r_=reduced_col_num;
	reduced_basis_=new double *[reduced_col_num];
	// if(reduced_row_num!=simulation_mesh_->getNumVertices()*3)
	// {
	// 	std::cout<<"The input reduced basis is not correct!\n";
	// 	return false;
	// }
	for(unsigned int i=0;i<reduced_col_num;++i)
	{
		reduced_basis_[i]=new double[reduced_row_num];
	}
	for(unsigned int i=0;i<reduced_col_num;++i)
	{
		for(unsigned int j=0;j<reduced_row_num;++j)
			reduced_basis_[i][j]=0.0;
	}
	unsigned int line_num=0,str_num=0,total_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		++total_num;
		double temp_value;
		input_file>>temp_value;
		reduced_basis_[str_num][line_num]=temp_value;
		if(str_num>=reduced_col_num-1)
		{
			str_num=0;
			if(line_num<reduced_row_num-1)
				++line_num;
		}
		else
			++str_num;
		if(total_num>=reduced_row_num*reduced_col_num)
			break;
	}
	input_file.close();
	// std::cout<<reduced_basis_[reduced_col_num-1][reduced_row_num-1];
	//normalize eigenfunction with respect to the s-inner product: <f,f>=1 on the volumetric vertices
	// for(int i=0;i<r_;++i)
	// {
	// 	double sum=0.0;
	// 	for(int j=0;j<3*simulation_mesh_->getNumVertices();++j)
	// 		sum+=reduced_basis_[i][j]*reduced_basis_[i][j];
	// 	if(sum>epsilon_)
	// 		sum=sqrt(sum);
	// 	else
	// 		std::cout<<"object eigenfunction sum is 0\n";
	// 	for(int j=0;j<3*simulation_mesh_->getNumVertices();++j)
	// 		reduced_basis_[i][j]/=sum;
	// }
	isload_reduced_basis_ = true;
    return true;
}
//format:.eigencoef
//first line:*eigenValues; second line:eigenvalues_num; third line: eigen values
//*eigenVectors; eigenfunctions_row_num, eigenfunctions_col_num
bool RealTimeExampleBasedDeformer::loadObjectEigenfunctions(const std::string &file_name)
{
	std::cout<<"Load object eigenfunction:\n";
    std::fstream input_file(file_name.c_str());

	if(!input_file)
	{
		std::cout<<"Error: Cannot open file "<<file_name<<std::endl;
		return false;
	}
    std::string temp_str;
	//read file,save as object_eigenvalues_[eigen_idx]
	while(std::getline(input_file,temp_str))
	{
		if(temp_str.compare(0,12,string("*eigenValues"))==0)
			break;
	}
    std::string reconstruct_eigenfunction_num_str;
	std::getline(input_file,reconstruct_eigenfunction_num_str);
	reconstruct_eigenfunction_num_=atoi(reconstruct_eigenfunction_num_str.c_str());
	object_eigenvalues_=new double[reconstruct_eigenfunction_num_];

	unsigned int str_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		double temp_value;
		input_file>>temp_value;
		object_eigenvalues_[str_num]=temp_value;
		str_num++;
		if(str_num>=reconstruct_eigenfunction_num_)
			break;
	}
	std::cout<<std::endl;
	str_num=0;
	// std::cout<<"object eigenvalues:\n";
	// for(int i=0;i<reconstruct_eigenfunction_num_;++i)
	// 	std::cout<<object_eigenvalues_[i]<<",";
	//read eigenvectors, save as object_eigenfunctions_[eigen_idx][vert_idx]
	while(std::getline(input_file,temp_str))
	{
		if(temp_str.compare(0,13,string("*eigenVectors"))==0)
			break;
	}
	unsigned int eigenfunction_row_num=0,eigenfunction_col_num=0;
	string eigenfunction_row_num_str,eigenfunction_col_num_str;
	std::getline(input_file,eigenfunction_row_num_str);
	eigenfunction_row_num=atoi(eigenfunction_row_num_str.c_str());
	std::getline(input_file,eigenfunction_col_num_str);
	eigenfunction_col_num=atoi(eigenfunction_col_num_str.c_str());
	object_eigenfunctions_=new double *[eigenfunction_col_num];
	for(unsigned int i=0;i<eigenfunction_col_num;++i)
	{
		object_eigenfunctions_[i]=new double [eigenfunction_row_num];
		for(unsigned int j=0;j<eigenfunction_row_num;++j)
		{
			object_eigenfunctions_[i][j]=0.0;
		}
	}
	unsigned int line_num=0;
	str_num=0;
	unsigned int total_str_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		++total_str_num;
		double temp_value1;
		input_file>>temp_value1;
		object_eigenfunctions_[str_num][line_num]=temp_value1;
		if(str_num>=eigenfunction_col_num-1)
		{
			str_num=0;
			if(line_num<eigenfunction_row_num-1)
				++line_num;
		}
		else
			++str_num;
		if(total_str_num>=eigenfunction_row_num*eigenfunction_col_num)
			break;
	}
	input_file.close();
	//normalize eigenfunction with respect to the s-inner product: <f,f>=1 on the volumetric vertices
	for(int i=0;i<eigenfunction_col_num;++i)
	{
		double sum=0.0;
		for(int j=0;j<simulation_mesh_->getNumVertices();++j)
			sum+=object_eigenfunctions_[i][j]*object_eigenfunctions_[i][j]*object_vertex_volume_[j];
		if(sum>epsilon_)
			sum=sqrt(sum);
		else
			std::cout<<"object eigenfunction sum is 0\n";
		for(int j=0;j<simulation_mesh_->getNumVertices();++j)
			object_eigenfunctions_[i][j]/=sum;
	}
	//project on eigenfunctions
	double *dis=new double[3*simulation_mesh_->getNumVertices()];
	for(int i=0;i<3*simulation_mesh_->getNumVertices();++i)
		dis[i]=0.0;
	object_eigencoefs_=new Vec3d[eigenfunction_col_num];
    projectOnEigenFunctions(simulation_mesh_,dis,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
							eigenfunction_col_num,object_eigencoefs_,object_affected_vertices_num_,object_affected_vertices_);

	delete[] dis;
	if(simulation_mode_==REDUCEDSPACE)
    {
		//temp compute initial object eigen skeleton
		// initial_eigen_skeleton_=new double[3*simulation_mesh_->getNumVertices()];
		// local_detail_vector_=new double[3*simulation_mesh_->getNumVertices()];
		// memset(local_detail_vector_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
		// reconstructFromEigenCoefs(object_eigencoefs_,1);
		target_reconstruction_=new double[3*simulation_mesh_->getNumVertices()];
		memset(target_reconstruction_,0.0,sizeof(double)*3*simulation_mesh_->getNumVertices());
	}
		// std::cout<<"loadObjectEigenfunctions---end.\n";
	if(isload_LB_cubica_)
	{
		cubica_LB_tetsubBasis_=new double**[object_LB_cubica_ele_num_];
		for(int i=0;i<object_LB_cubica_ele_num_;++i)
		{
			cubica_LB_tetsubBasis_[i]=new double*[4];
			for(int j=0;j<4;++j)
				cubica_LB_tetsubBasis_[i][j]=new double[interpolate_eigenfunction_num_];
		}
		for(int i=0;i<object_LB_cubica_ele_num_;++i)
		{
			int ele=object_LB_cubica_elements_[i];
			for(int j=0;j<4;++j)
	        {
	            int vertID=simulation_mesh_->getVertexIndex(ele,j);
	            for(int k=0;k<interpolate_eigenfunction_num_;++k)
					cubica_LB_tetsubBasis_[i][j][k]=object_eigenfunctions_[k][vertID];
	        }
		}
	}
	// std::cout<<"object_eigencoefs:\n";
	// for(int j=0;j<45;++j)
	// 	std::cout<<object_eigencoefs_[j]<<",";
	//compute LB to reduced basis matrix:
	LB_to_reducespace_=new double*[r_];
	for(unsigned int i=0;i<r_;++i)
	{
		LB_to_reducespace_[i]=new double[3*interpolate_eigenfunction_num_];
		for(unsigned int j=0;j<3*interpolate_eigenfunction_num_;++j)
			LB_to_reducespace_[i][j]=0.0;
	}
	//reduced_basis_:[basis_num][vertex_num];-----rx3n
	//eigenfunction:[eigen_num][vertex_num];------mxn
	for(unsigned int i=0;i<r_;++i)
	{
		for(unsigned int j=0;j<interpolate_eigenfunction_num_;++j)
		{
			for(unsigned int k=0;k<simulation_mesh_->getNumVertices();++k)
			{
				for(unsigned int l=0;l<3;++l)
					LB_to_reducespace_[i][3*j+l]+=reduced_basis_[i][3*k+l]*object_eigenfunctions_[j][k]/object_eigenvalues_[j];
			}
		}
	}
	//compute reduced space to LB space

	reducespace_to_LB_=new double*[r_];
	for(unsigned int i=0;i<r_;++i)
	{
		reducespace_to_LB_[i]=new double[3*interpolate_eigenfunction_num_];
		for(unsigned int j=0;j<3*interpolate_eigenfunction_num_;++j)
			reducespace_to_LB_[i][j]=0.0;
	}
	//reduced_basis_:[basis_num][vertex_num];-----rx3n
	//eigenfunction:[eigen_num][vertex_num];------mxn
	for(unsigned int i=0;i<r_;++i)
	{
		for(unsigned int j=0;j<interpolate_eigenfunction_num_;++j)
		{
			for(unsigned int k=0;k<simulation_mesh_->getNumVertices();++k)
			{
				for(unsigned int l=0;l<3;++l)
					reducespace_to_LB_[i][3*j+l]+=reduced_basis_[i][3*k+l]*object_eigenfunctions_[j][k]*object_vertex_volume_[k]*object_eigenvalues_[j];
			}
		}
	}
	// std::cout<<"aaaaa\n";
	// 	//temp--4 reconstruct
	// 	temp_eigencoefs_=new Vec3d[reconstruct_eigenfunction_num_];
	//     projectOnEigenFunctions(simulation_mesh_,examples_deformation_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
	// 							reconstruct_eigenfunction_num_,temp_eigencoefs_);
	// std::cout<<"Bbbb\n";
	//3-reconstruct
	// temp_eigencoefs0_=new Vec3d[reconstruct_eigenfunction_num_];
    // projectOnEigenFunctions(simulation_mesh_,examples_deformation0_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
	// 						reconstruct_eigenfunction_num_,temp_eigencoefs0_);
	// reconstructFromEigenCoefs(temp_eigencoefs0_,1);


	std::cout<<"object_eigencoefs:\n";
		for(int j=0;j<interpolate_eigenfunction_num_;++j)
			std::cout<<object_eigencoefs_[j]<<",";
	// getchar();
	std::cout<<"load object eigenfunctions done.\n";
    return true;
}

bool RealTimeExampleBasedDeformer::loadExampleEigenFunctions(const std::string &file_name_prefix)
{
   	std::cout<<"Load example eigenfunctions:\n";
	if(example_num_==0)
	{
		std::cout<<"eigenfunction num for examples is zero.\n";
		return false;
	}
	example_eigenfunctions_ = new double **[example_num_];
	example_eigenvalues_ = new double *[example_num_];
	example_eigencoefs_=new Vec3d*[example_num_];
	object_example_eigencoefs_=new Vec3d*[example_num_];
	// for(int i=0;i<example_num_;++i)
	// 	object_example_eigencoefs_[i]=new Vec3d[reconstruct_eigenfunction_num_];

	for(unsigned int ex_num=0;ex_num<example_num_;++ex_num)
	{
        std::string file_num_str,file_name;
        std::stringstream adaptor;
		adaptor.str("");
		adaptor.clear();
		adaptor<<ex_num;
		adaptor>>file_num_str;
		//read file, file format is .eigen,
		file_name=file_name_prefix+file_num_str+std::string(".eigen");
        std::fstream input_file(file_name.c_str());
		if(!input_file)
		{
			std::cout<<"Error: Cannot open file "<<file_name<<std::endl;
			return false;
		}
        std::string temp_str;
		//read eigenvalues first,save as example_eigenvalues_[ex_idx][eigen_idx]
		while(std::getline(input_file,temp_str))
		{
			if(temp_str.compare(0,12,std::string("*eigenValues"))==0)
				break;
		}
        std::string reconstruct_eigenfunction_num_str;
		std::getline(input_file,reconstruct_eigenfunction_num_str);
		reconstruct_eigenfunction_num_=atoi(reconstruct_eigenfunction_num_str.c_str());
		example_eigenvalues_[ex_num]=new double[reconstruct_eigenfunction_num_];
		for(unsigned int j=0;j<reconstruct_eigenfunction_num_;++j)
		{
			example_eigenvalues_[ex_num][j]=0.0;
		}
		unsigned int str_num=0;
		while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
		{
			double temp_value;
			input_file>>temp_value;
			example_eigenvalues_[ex_num][str_num]=temp_value;
			str_num++;
			if(str_num>=reconstruct_eigenfunction_num_)
				break;
		}
		//read eigenvectors, save as example_eigenfunctions_[ex_idx][eigen_idx][vert_idx]
		while(std::getline(input_file,temp_str))
		{
			if(temp_str.compare(0,13,std::string("*eigenVectors"))==0)
				break;
		}
		unsigned int eigenfunction_row_num=0,eigenfunction_col_num=0;
        std::string eigenfunction_row_num_str,eigenfunction_col_num_str;
		std::getline(input_file,eigenfunction_row_num_str);
		std::getline(input_file,eigenfunction_col_num_str);
		eigenfunction_row_num=atoi(eigenfunction_row_num_str.c_str());
		eigenfunction_col_num=atoi(eigenfunction_col_num_str.c_str());
		example_eigenfunctions_[ex_num]=new double *[eigenfunction_col_num];
		for(unsigned int i=0;i<eigenfunction_col_num;++i)
		{
			example_eigenfunctions_[ex_num][i]=new double [eigenfunction_row_num];
		}
		for(unsigned int j=0;j<eigenfunction_col_num;++j)
			for(unsigned int k=0;k<eigenfunction_row_num;++k)
				example_eigenfunctions_[ex_num][j][k]=0.0;
		unsigned int line_num=0;
		str_num=0;
		unsigned int total_str_num=0;
		while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
		{
			++total_str_num;
			double temp_value1=0.0;
			input_file>>temp_value1;
			example_eigenfunctions_[ex_num][str_num][line_num]=temp_value1;
			if(str_num>=eigenfunction_col_num-1)
			{
				str_num=0;
				if(line_num<eigenfunction_row_num-1)
					++line_num;
			}
			else
				++str_num;
			if(total_str_num>=eigenfunction_row_num*eigenfunction_col_num)
				break;
		}
		input_file.close();
		//normalize example[i] eigenfunction with respect to the s-inner product: <f,f>=1 on the volumetric vertcies
		for(int i=0;i<eigenfunction_col_num;++i)
		{
			double sum=0.0;
			for(int j=0;j<examples_[ex_num]->getNumVertices();++j)
			{
				sum+=example_eigenfunctions_[ex_num][i][j]*example_eigenfunctions_[ex_num][i][j]*example_vertex_volume_[ex_num][j];
			}
			if(sum>epsilon_)
				sum=sqrt(sum);
			else
				std::cout<<"object eigenfunction sum is 0\n";
			for(int j=0;j<examples_[ex_num]->getNumVertices();++j)
				example_eigenfunctions_[ex_num][i][j]/=sum;
		}
		example_eigencoefs_[ex_num]=new Vec3d[eigenfunction_col_num];
        double *dis=new double[3*examples_[ex_num]->getNumVertices()];
        for(int j=0;j<3*examples_[ex_num]->getNumVertices();++j)
            dis[j]=0.0;
		// if(ex_num==2)
        // 	for(int j=0;j<examples_[ex_num]->getNumVertices();++j)
        //     	dis[3*j]=-1.0;

		// if((!with_constrains_)&&(simulation_mode_==REDUCEDSPACE))
		// {
		// 	rigid_->GetPosition(&t_[0],&t_[1],&t_[2]);
		// 	for(int j=0;j<examples_[ex_num]->getNumVertices();++j)
		// 	{
		// 		for(int k=0;k<3;++k)
		// 			dis[3*j+k]=t_[k];
		// 	}
		// }
        projectOnEigenFunctions(examples_[ex_num],dis,example_vertex_volume_[ex_num],
                                example_eigenfunctions_[ex_num],example_eigenvalues_[ex_num],
                                eigenfunction_col_num,example_eigencoefs_[ex_num],
								example_affected_vertices_num_[ex_num],example_affected_vertices_[ex_num]);

		// reconstructFromEigenCoefs(example_eigencoefs_[ex_num],dis);
		// for(int j=0;j<20;++j)
		// 	if(ex_num==2)
		// 		std::cout<<dis[3*j]<<","<<dis[3*j+1]<<","<<dis[3*j+2]<<"\n";
		// 	else
		// 		std::cout<<dis[3*j]<<","<<dis[3*j+1]<<","<<dis[3*j+2]<<"\n";
		// getchar();
		//object projectOnExampleEigunFunctions
		// projectOnEigenFunctions(examples_[ex_num],dis,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
		// 						reconstruct_eigenfunction_num_,object_example_eigencoefs_[ex_num]);
		// std::cout<<"a\n";
		// if(ex_num==1)
		// {
		// // 	Vec3d *temp_eigencoefs=new Vec3d[interpolate_eigenfunction_num_];
		// // 	for(unsigned int i=0;i<interpolate_eigenfunction_num_;++i)
		// // 		for(unsigned int j=0;j<3;++j)
		// // 			temp_eigencoefs[i][j]=example_eigencoefs_[ex_num][i][j];//object_eigencoefs_[i][j];
		// //
		// // 	reconstructFromEigenCoefs(temp_eigencoefs,example_guided_deformation_);
		// 	reconstructFromEigenCoefs(object_eigencoefs_,2);
		// }

		delete[] dis;
		std::cout<<"example-example_eigencoefs_ for ex "<<ex_num<<"\n";
		for(int i=0;i<3;++i)
		{
			for(int j=0;j<10;++j)
				std::cout<<example_eigencoefs_[ex_num][i][j]<<",";
			std::cout<<"\n";
		}
	}
	// double *dis;
	// memset(dis,0.0,sizeof(double)*3*examples_[1]->getNumVertices());
	// reconstructFromEigenCoefs(example_eigencoefs_[1],example_guided_deformation_);
	// delete[] dis;
	// double *dis=new double[3*simulation_vertices_num_];
	// for(int j=0;j<3*simulation_vertices_num_;++j)
	// 	dis[j]=0.0;
	// reconstructFromEigenCoefs(example_eigencoefs_[1],dis);
	// 		delete[] dis;
	//Presume system have 2 examples
	// example_length_[0]=example_length_[1]=example_length_[2]=1.0;
	// Vec3d *example_diff=new Vec3d[interpolate_eigenfunction_num_];
	// std::cout<<"\n";
	// std::cout<<"object_example_eigencoefs_0:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	std::cout<<object_example_eigencoefs_[0][i]<<",";
	// }
	// std::cout<<"\n";
	//
	// std::cout<<"\n";
	// std::cout<<"object_example_eigencoefs_1:\n";
	// Vec3d example_example_length;
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	std::cout<<object_example_eigencoefs_[1][i]<<",";
	// }
	// std::cout<<"\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	for(int j=0;j<3;++j)
	// 	{
	// 		example_diff[i][j]=object_example_eigencoefs_[1][i][j]-object_example_eigencoefs_[0][i][j];
	// 		example_length_[j]+=example_diff[i][j]*example_diff[i][j];
	// 	}
	// }
	// for(int j=0;j<3;++j)
	// 	example_length_[j]=sqrt(example_length_[j]);//b1-a1
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	for(int j=0;j<3;++j)
	// 	{
	// 		example_diff[i][j]=example_eigencoefs_[1][i][j]-example_eigencoefs_[0][i][j];
	// 		example_example_length[j]+=example_diff[i][j]*example_diff[i][j];
	// 	}
	// }
	// for(int j=0;j<3;++j)
	// 	example_example_length[j]=sqrt(example_example_length[j]);//b0-a0
	// for(int j=0;j<3;++j)
	// 	example_length_[j]=example_example_length[j]/example_length_[j];
	//
	// delete[] example_diff;
	// std::cout<<"loadexampleEigenfunctions---end.\n";
    return true;
}

bool RealTimeExampleBasedDeformer::loadPlanesInScene(const std::string &file_name, unsigned int plane_num)
{
    if(file_name=="")
	{
		std::cout<<"Error: cannot read "<<file_name<<".\n";
		return false;
	}
	std::cout<<"file_name:"<<file_name<<"\n";
	std::cout<<"plane_num:"<<plane_num<<"\n";
	plane_num_=plane_num;
	planes_=new Planes(file_name.c_str(),plane_num);
    if(planes_ == NULL)
    {
        std::cout<<"Error: failed to load planes from "<<file_name<<".\n";
        return false;
    }
    return true;
}

bool RealTimeExampleBasedDeformer::loadFixedVertices(const std::string &file_name)
{
	if(file_name=="")
	{
		std::cout<<"Error: cannot read "<<file_name<<".\n";
		return false;
	}
	if (LoadList::load(file_name.c_str(), &fixed_vertex_num_,&fixed_vertices_) != 0)
	{
		printf("Error reading fixed vertices.\n");
		exit(1);
	}
	LoadList::sort(fixed_vertex_num_, fixed_vertices_);

	printf("Loaded %d fixed vertices. They are:\n",fixed_vertex_num_);
	LoadList::print(fixed_vertex_num_,fixed_vertices_);
	// create 0-indexed fixed DOFs
	fixed_dofs_num_ = 3 * fixed_vertex_num_;
	fixed_dofs_ = new int[fixed_dofs_num_];
	for(int i=0; i<fixed_vertex_num_; i++)
	{
		fixed_dofs_[3*i+0] = 3*fixed_vertices_[i]-3;
		fixed_dofs_[3*i+1] = 3*fixed_vertices_[i]-2;
		fixed_dofs_[3*i+2] = 3*fixed_vertices_[i]-1;
	}
	for(int i=0; i<fixed_vertex_num_; i++)
		fixed_vertices_[i]--;
	printf("Boundary vertices processed.\n");
	std::cout<<"Load fixed vertices file done.\n";
	return true;
}
//tetID : 1-indexed,load cubica data for reduced simulation
bool RealTimeExampleBasedDeformer::loadObjectCubicaData(const std::string &file_name)
{
	std::cout<<"Load object cubica data...\n";
	std::fstream input_file(file_name.c_str());
	if(!input_file)
	{
		std::cout<<"Error: failed to open file "<<file_name<<".\n";
		return false;
	}
	string temp_str;
	std::getline(input_file,temp_str);
	object_cubica_ele_num_=atoi(temp_str.c_str());
	object_cubica_elements_ = new unsigned int[object_cubica_ele_num_];
	object_cubica_weights_ = new double [object_cubica_ele_num_];
	while(std::getline(input_file,temp_str))
	{
		if(temp_str.compare(0,6,string("*tetID"))==0)
			break;
	}
	unsigned int str_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		unsigned int temp_value;
		input_file>>temp_value;
		object_cubica_elements_[str_num]=temp_value-1;
		str_num++;
		if(str_num>=object_cubica_ele_num_)
			break;
	}
	while(std::getline(input_file,temp_str))
	{
		if(temp_str.compare(0,10,string("*tetWeight"))==0)
			break;
	}
	str_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		double temp_value;
		input_file>>temp_value;
		object_cubica_weights_[str_num]=temp_value;
		++str_num;
		if(str_num>=object_cubica_ele_num_)
			break;
	}
	// std::cout<<"Load object cubica data succeed.\n";

	//test all cubica elements:
	// object_cubica_ele_num_=simulation_mesh_->getNumElements();
	// object_cubica_elements_=new unsigned int[object_cubica_ele_num_];
	// object_cubica_weights_=new double[object_cubica_ele_num_];
	// for(int i=0;i<object_cubica_ele_num_;++i)
	// {
	// 	object_cubica_elements_[i]=i;
	// 	double ele_volume=simulation_mesh_->getElementVolume(i);
	// 	object_cubica_weights_[i]=1.0*ele_volume;
	// }
	isload_cubica_ = true;
	return true;
}
//load cubica data for LB subspace : eigenfunction as input basis
bool RealTimeExampleBasedDeformer::loadObjectLBCubicaData(const std::string &file_name)
{
	std::cout<<"Load object cubica data..."<<file_name<<"\n";
	std::fstream input_file(file_name.c_str());
	if(!input_file)
	{
		std::cout<<"Error: failed to open file "<<file_name<<".\n";
		return false;
	}
	string temp_str;
	std::getline(input_file,temp_str);
	object_LB_cubica_ele_num_=atoi(temp_str.c_str());
	object_LB_cubica_elements_ = new unsigned int[object_LB_cubica_ele_num_];
	object_LB_cubica_weights_ = new double [object_LB_cubica_ele_num_];
	while(std::getline(input_file,temp_str))
	{
		if(temp_str.compare(0,6,string("*tetID"))==0)
			break;
	}
	unsigned int str_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		unsigned int temp_value;
		input_file>>temp_value;
		object_LB_cubica_elements_[str_num]=temp_value-1;
		str_num++;
		if(str_num>=object_LB_cubica_ele_num_)
			break;
	}
	while(std::getline(input_file,temp_str))
	{
		if(temp_str.compare(0,10,string("*tetWeight"))==0)
			break;
	}
	str_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		double temp_value;
		input_file>>temp_value;
		object_LB_cubica_weights_[str_num]=temp_value;
		++str_num;
		if(str_num>=object_LB_cubica_ele_num_)
			break;
	}
	isload_LB_cubica_ = true;
	return true;
}
bool RealTimeExampleBasedDeformer::saveSimulationMesh(const std::string &file_name) const
{
    if(simulation_mesh_==NULL)
	{
		std::cout<<"Error: simulation mesh is null.\n";
		return false;
	}
    std::ofstream output_file(file_name.c_str());
	if(!output_file)
	{
		std::cout<<"Error: failed to open "<<file_name<<".\n";
		return false;
	}
	output_file<<"*VERTICES"<<std::endl;
	output_file<<simulation_mesh_->getNumVertices()<<" 3 0 0"<<std::endl;
	for(unsigned int i=0;i<simulation_mesh_->getNumVertices();++i)
	{
		output_file<<i+1<<" "<<(*simulation_mesh_->getVertex(i))[0]<<" ";
		output_file<<(*simulation_mesh_->getVertex(i))[1]<<" ";
		output_file<<(*simulation_mesh_->getVertex(i))[2]<<std::endl;
	}
	output_file<<"*ELEMENTS"<<std::endl;
	output_file<<"TET"<<std::endl;
	output_file<<simulation_mesh_->getNumElements()<<" 4 0"<<std::endl;
	for(unsigned int i=0;i<simulation_mesh_->getNumElements();++i)
	{
		output_file<<i+1;
		for(unsigned int j=0;j<simulation_mesh_->getNumElementVertices();++j)
		{
			unsigned int global_idx=simulation_mesh_->getVertexIndex(i,j);
			output_file<<" "<<global_idx+1;
		}
		output_file<<std::endl;
	}
	output_file.close();
    return true;
}
bool RealTimeExampleBasedDeformer::saveExamples(const std::string &file_name_prefix) const
{
	std::string num_str;
	std::stringstream stream;
	for(unsigned int i=0;i<example_num_;++i)
	{
		if(examples_[i]==NULL)
		{
			std::cout<<"Error: example mesh "<<i<<" is null.\n";
		    return false;
		}
		stream.str("");
		stream.clear();
		stream<<i;
		stream>>num_str;
		std::string file_name=file_name_prefix+num_str+".smesh";
        std::ofstream output_file(file_name.c_str());
		if(!output_file)
		{
			std::cout<<"Error: failed to open "<<file_name<<".\n";
			return false;
		}
		output_file<<"*VERTICES"<<std::endl;
		output_file<<examples_[i]->getNumVertices()<<" 3 0 0"<<std::endl;
		for(unsigned int j=0;j<examples_[i]->getNumVertices();++j)
		{
			output_file<<j+1<<" "<<(*examples_[i]->getVertex(j))[0]<<" ";
			output_file<<(*examples_[i]->getVertex(j))[1]<<" ";
			output_file<<(*examples_[i]->getVertex(j))[2]<<std::endl;
		}
		output_file<<"*ELEMENTS"<<std::endl;
		output_file<<"TET"<<std::endl;
		output_file<<examples_[i]->getNumElements()<<" 4 0"<<std::endl;
		for(unsigned int j=0;j<examples_[i]->getNumElements();++j)
		{
			output_file<<j+1;
			for(unsigned int k=0;k<examples_[i]->getNumElementVertices();++k)
			{
				unsigned int global_idx=examples_[i]->getVertexIndex(j,k);
				output_file<<" "<<global_idx+1;
			}
			output_file<<std::endl;
		}
		output_file.close();
	}
    return true;
}

bool RealTimeExampleBasedDeformer::saveVisualMesh(const std::string &file_name) const
{
    //need to be modified later
	if(visual_mesh_==NULL)
	{
		std::cout<<"Error: the visual mesh is null.\n";
		return false;
	}
	ObjMesh *mesh=visual_mesh_->GetMesh();
	mesh->save(file_name.c_str(),0);
	std::cout<<file_name<<"saved.\n";
    return true;
}

bool RealTimeExampleBasedDeformer::saveObjectEigenfunctions(const std::string &file_name) const
{
	if(simulation_mesh_==NULL)
	{
		std::cout<<"The simulation mesh is null.\n";
		return false;
	}
	if((object_eigenfunctions_==NULL)||(object_eigenvalues_==NULL))
	{
		std::cout<<"Eigenfunctions or eigenvalues for object is null.\n";
		return false;
	}
    std::ofstream output_file(file_name.c_str());
	std::cout<<file_name<<std::endl;
	if(!output_file)
	{
		std::cout<<"Error:unable to open file "<<file_name<<".\n";
		return false;
	}
	//write eigenvalues first
	output_file<<"*eigenValues"<<std::endl;
	output_file<<reconstruct_eigenfunction_num_<<std::endl;
	for(unsigned int i=0;i<reconstruct_eigenfunction_num_;++i)
	{
		output_file<<object_eigenvalues_[i]<<" ";
	}
	output_file<<std::endl;
	//write eigenvectors
	output_file<<"*eigenVectors"<<std::endl;
	unsigned int vert_num=simulation_mesh_->getNumVertices();
	output_file<<vert_num<<std::endl;
	output_file<<reconstruct_eigenfunction_num_<<std::endl;
	for(unsigned int i=0;i<vert_num;++i)
	{
		for(unsigned int j=0;j<reconstruct_eigenfunction_num_;++j)
		{
			output_file<<object_eigenfunctions_[j][i]<<" ";
		}
		output_file<<std::endl;
	}
	output_file.close();
	std::cout<<"Simulation object eigen function saved.\n";
    return true;
}

bool RealTimeExampleBasedDeformer::saveExampleEigenfunctions(const std::string &file_name_prefix) const
{
	std::string num_str;
	std::stringstream stream;
	for(unsigned int i=0;i<example_num_;++i)
	{
		if(examples_[i]==NULL)
		{
			std::cout<<"Error: example mesh "<<i<<" is null.\n";
			return false;
		}
		if((example_eigenfunctions_[i]==NULL)||(example_eigenvalues_[i]==NULL))
		{
			std::cout<<"Error: eigenfunction or eigenvalue for example "<<i<<" is null.\n";
			return false;
		}
		stream.str("");
		stream.clear();
		stream<<i;
		stream>>num_str;
		std::string file_name=file_name_prefix+num_str+std::string(".eigen");
		std::cout<<file_name<<std::endl;
        std::ofstream output_file(file_name.c_str());
		if(!output_file)
		{
			std::cout<<"Error:unable to open file "<<file_name<<".\n";
			return false;
		}
		//write file,eigenvalues first
		output_file<<"*eigenValues"<<std::endl;
		output_file<<reconstruct_eigenfunction_num_<<std::endl;
		for(unsigned int j=0;j<reconstruct_eigenfunction_num_;++j)
		{
			output_file<<example_eigenvalues_[i][j]<<" ";
		}
		output_file<<std::endl;
		//write eigenvectors
		output_file<<"*eigenVectors"<<std::endl;
		unsigned int vert_num=examples_[i]->getNumVertices();
		output_file<<vert_num<<std::endl;
		output_file<<reconstruct_eigenfunction_num_<<std::endl;
		for(unsigned int j=0;j<vert_num;++j)
		{
			for(unsigned int k=0;k<reconstruct_eigenfunction_num_;++k)
			{
				output_file<<example_eigenfunctions_[i][k][j]<<" ";
			}
			output_file<<std::endl;
		}
		output_file.close();
		std::cout<<"Example "<<i<<" eigen function saved.\n";
	}
	return true;
}

 VolumetricMesh* RealTimeExampleBasedDeformer::exampleMesh(unsigned int example_idx) const
{
	assert(example_idx < example_num_);
	return examples_[example_idx];
}

//file format of corresponding function file:
//first line is an integer p indicating the number of corresponding functions (regions on the object)
//the following are (1+example_pos_num)*p  lines, every (1+example_num_) consecutive lines are a group
//meaning corresponding region on the object and examples
//each line is a list of integers separated by comma, meanning the point in this region
bool RealTimeExampleBasedDeformer::loadCorrespondenceData(const std::string &file_name)
{
	std::cout<<"Load input correspondence data:\n";
	if(simulation_mesh_==NULL)
	{
		std::cout<<"Error: the simulation mesh is null.\n";
		return false;
	}
	std::cout<<simulation_mesh_->getNumVertices()<<"\n";
	if(example_num_==0)
	{
		std::cout<<"Error: the input example num is zero.\n";
		return false;
	}
	else
	{
		for(unsigned int i=0;i<example_num_;++i)
		{
			if(examples_[i]==NULL)
			{
				std::cout<<"Error: example "<<i<<" tet mesh is null.\n";
				return false;
			}
		}
	}
    std::fstream input_file(file_name.c_str());
	if(!input_file)
	{
		std::cout<<"Error: failed to open "<<file_name<<std::endl;
		return false;
	}
	//read corresponding function number from file
	std::string line_str;
	std::getline(input_file,line_str);
	corresponding_function_num_=atoi(line_str.c_str());
	//allocate space for corresponding functions, each col is a function
	//save as object_corresponding_functions_[vert_idx][cor_region_idx];
	//save as example_corresponding_functions_[ex_idx][vert_idx][cor_region_idx];
	int object_eigen_vert_num=simulation_mesh_->getNumVertices();

	object_corresponding_functions_= new double *[object_eigen_vert_num];
	for(unsigned int i=0;i<object_eigen_vert_num;++i)
	{
		object_corresponding_functions_[i]=new double[corresponding_function_num_];
		for(unsigned int j=0;j<corresponding_function_num_;++j)
			object_corresponding_functions_[i][j]=0.0;
	}
	example_corresponding_functions_=new double **[example_num_];
	for(unsigned int i=0;i<example_num_;++i)
	{
		int ex_eigen_vert_num=examples_[i]->getNumVertices();
		example_corresponding_functions_[i]=new double *[ex_eigen_vert_num];
		for(unsigned int j=0;j<ex_eigen_vert_num;++j)
		{
			example_corresponding_functions_[i][j]=new double[corresponding_function_num_];
			for(unsigned int k=0;k<corresponding_function_num_;++k)
				example_corresponding_functions_[i][j][k]=0.0;
		}
	}
	//read data from file
	for(unsigned int i=0;i<corresponding_function_num_;++i)
	{
		//points on the object
		std::getline(input_file,line_str);
		int j=0;
		int points_in_region=0;
		while(j<line_str.size())
		{
			int k=j;
			while(line_str[k]!=';')
				++k;
			std::string idx_str=line_str.substr(j,k-j);
			int idx=atoi(idx_str.c_str());
			object_corresponding_functions_[idx-1][i]=1;
			j=k+1;
			++points_in_region;
		}
		if(points_in_region>1)//corresponding functions are region based
			is_region_based_correspondence_=true;
		for(unsigned int ex_idx=0;ex_idx<example_num_;++ex_idx)
		{
			//points on the examples
			std::getline(input_file,line_str);
			j=0;
			while(j<line_str.size())
			{
				int k=j;
				while(line_str[k]!=';')
					++k;
				std::string idx_str=line_str.substr(j,k-j);
				int idx=atoi(idx_str.c_str());
				example_corresponding_functions_[ex_idx][idx-1][i]=1;
				j=k+1;
			}
		}
	}
	input_file.close();
	std::cout<<"Loaded "<<corresponding_function_num_<<" corresponding functions from: "<<file_name<<std::endl;
    return true;
}

bool RealTimeExampleBasedDeformer::registerEigenfunctions()
{
	if(object_eigenfunctions_==NULL)
	{
		std::cout<<"Error: object_eigenfunction is not loaded.\n";
		return false;
	}
	if(simulation_mesh_==NULL)
	{
		std::cout<<"Error: simulation volumetric mesh is not loaded.\n";
		return false;
	}
	if(corresponding_function_num_==0)
	{
		std::cout<<"Error: corresponding data is not loaded.\n";
	}
	if(example_num_==0)
	{
		std::cout<<"Error: example num is zero.\n";
		return false;
	}
	else
	{
		for(unsigned int i=0;i<example_num_;++i)
		{
			if(examples_[i]==NULL)
			{
				std::cout<<"Error: example "<<i<<" tet mesh is not loaded.\n";
				return false;
			}
			if(example_eigenfunctions_[i]==NULL)
			{
				std::cout<<"Error: eigenfunction of example "<<i<<" is not loaded.\n";
				return false;
			}
			if(example_corresponding_functions_[i]==NULL)
			{
				std::cout<<"Error: example "<<i<<" corresponding function is not loaded.\n";
				return false;
			}
		}
	}
	if(object_corresponding_functions_==NULL)
	{
		std::cout<<"Error: object corresponding functions are not loaded.\n";
		return false;
	}
	int row_b =simulation_mesh_->getNumVertices();
	//object_vertex_volume_=new double [row_b];
	int col_b = interpolate_eigenfunction_num_;
	//change to col order,each col is an eigenfunction
	double **obj_eigenfunction_col=new double *[row_b];
	double **new_obj_eigenfunction_col=new double *[row_b];
	for(unsigned int i=0;i<row_b;++i)
	{
		obj_eigenfunction_col[i]=new double[col_b];
		new_obj_eigenfunction_col[i]=new double[col_b];
	}
	//change to col order, each col is an eigenfunction
	for(int i=0;i<col_b;++i)
		for(int j=0;j<row_b;++j)
			obj_eigenfunction_col[j][i]=object_eigenfunctions_[i][j];
	std::cout<<"\n";
	for(int i=0;i<example_num_;++i)
	{
		double **G=object_corresponding_functions_,**F=example_corresponding_functions_[i];
		int p=corresponding_function_num_;
		int row_a=examples_[i]->getNumVertices();
		int col_a=interpolate_eigenfunction_num_;
		double **ex_eigenfunction_col=new double *[row_a];
		double **new_ex_eigenfunction_col=new double *[row_a];
		for(int j=0;j<row_a;++j)
		{
			ex_eigenfunction_col[j]=new double[col_a];
			new_ex_eigenfunction_col[j]=new double[col_a];
		}
		for(int j=0;j<col_a;++j)
			for(int k=0;k<row_a;++k)
				ex_eigenfunction_col[k][j]=example_eigenfunctions_[i][j][k];
		int col=col_a<col_b?col_a:col_b;
		//std::cout<<"col:"<<col_a<<",col_b:"<<col_b<<"\n";
		CoupledQuasiHarmonics quasi_base_generator(ex_eigenfunction_col,obj_eigenfunction_col,example_eigenvalues_[i],
													object_eigenvalues_,row_a,row_b,col,
													example_vertex_volume_[i],object_vertex_volume_,F,G,p);
		//couple example with object
		quasi_base_generator.setObjectiveType(CoupledQuasiHarmonics::ONLY_A);

		//set optimize package
		//quasi_base_generator.setPackage(CoupledQuasiHarmonics::NLOPT);
		quasi_base_generator.setPackage(CoupledQuasiHarmonics::OPTPP);

		//set different options for region-based corresponding and point-wise corresponding
		if(is_region_based_correspondence_)
		{
			quasi_base_generator.setInnerProductType(CoupledQuasiHarmonics::AREA_WEIGHTED);
			quasi_base_generator.enableScaleCorrespondingFunction();
		}
		else
		{
			quasi_base_generator.setInnerProductType(CoupledQuasiHarmonics::STANDARD);
			quasi_base_generator.disableScaleCorrespondingFunction();
		}
		//quasi_base_generator.setMu(100);
		std::cout<<"Coupling object and example "<<i+1<<std::endl;
		quasi_base_generator.getCoupledBases(/*coupled_eigenfunction_num_*/col,new_ex_eigenfunction_col,new_obj_eigenfunction_col);
		std::cout<<"Done.\n";

		//transform col-order eigen-bases to row-order so that each row is one eigenfunction
		for(int j=0;j<row_a;++j)
			for(int k=0;k<interpolate_eigenfunction_num_;++k)
				example_eigenfunctions_[i][k][j]=new_ex_eigenfunction_col[j][k];

		for(int j=0;j<row_a;++j)
		{
			if(ex_eigenfunction_col[j])
				delete[] ex_eigenfunction_col[j];
			if(new_ex_eigenfunction_col[j])
				delete[] new_ex_eigenfunction_col[j];
		}
		 delete[] ex_eigenfunction_col;
		 delete[] new_ex_eigenfunction_col;
	}
	for(int i=0;i<row_b;++i)
	{
		delete[] obj_eigenfunction_col[i];
		delete[] new_obj_eigenfunction_col[i];
	}
	delete[] obj_eigenfunction_col;
	delete[] new_obj_eigenfunction_col;
    return true;
}
void RealTimeExampleBasedDeformer::setGravity(bool add_gravity,double gravity)
{
	if(simulation_mode_==FULLSPACE)
	{
		isotropic_hyperelastic_fem_->SetGravity(add_gravity);
	}
	else
	{
		if(!with_constrains_)
		{
			if(!rigid_)
			{
				std::cout<<"rigid body doesn't defined!\n";
				exit(0);
			}
			add_gravity_=add_gravity;
			// if(add_gravity_)
			// 	rigid_->AddGravity(gravity_,0.0,1.0,0.0);
			// std::cout<<"----"<<add_gravity_<<"\n";
			// getchar();
			// rigid_->SetAngularVelocity(0.1,0,0.1);
			// rigid_->SetVelocity(-1.0,0.0,0.0);
			// if(!add_gravity_)
			// rigid_->ResetWrenches();
		}
	// 	if(material_mode_==REDUCED_NEOHOOKEAN)
	// 		reduced_neoHookean_force_model_->SetGravity(add_gravity,gravity,U_);
	// 	if(material_mode_==REDUCED_STVK)
	// 		reduced_stvk_cubature_force_model_->SetGravity(add_gravity,gravity,U_);
	}


}
void RealTimeExampleBasedDeformer::computeEnergy(const Vec3d *reduced_dis,double &energy) const
{
	energy=0.0;
	// double *temp=new double[3*interpolate_eigenfunction_num_];
	// memset(temp,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
	// for(unsigned int i=0;i<3*interpo)
	for(unsigned int i=0;i<interpolate_eigenfunction_num_;++i)
	{
		for(unsigned int j=0;j<3;++j)
		{
			energy+=reduced_dis[i][j]*reduced_dis[i][j];
			// temp[3*i+j]=reduced_dis[i][j];
		}
	}
	// for(unsigned int i=0;i<3*interpolate_eigenfunction_num_;++i)
	// 	energy+=temp[i]*temp[i];
	// delete[] temp;
}

void RealTimeExampleBasedDeformer::computeEnergyGrad(const Vec3d *reduced_dis,double *forces) const
{//r*1
	// PerformanceCounter p_counter;
	// p_counter.StartCounter();
	memset(forces,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
	for(unsigned int i=0;i<interpolate_eigenfunction_num_;++i)
	{
		for(unsigned int j=0;j<3;++j)
		{
			forces[3*i+j]=2.0*reduced_dis[i][j];
		}
	}
}
Mat3d RealTimeExampleBasedDeformer::computeDs(const double *reduced_dis) const
{//reduced dis is 12*1 for each element
	Mat3d Ds(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
	{
		Ds[0][i]=reduced_dis[3*i]-reduced_dis[9];
		Ds[1][i]=reduced_dis[3*i+1]-reduced_dis[10];
		Ds[2][i]=reduced_dis[3*i+2]-reduced_dis[11];
	}
	return Ds;
}
Mat3d RealTimeExampleBasedDeformer::computeDmInv(const int &ele) const
{//3x3
	int *global_idx=new int[4];
	Vec3d *vert_pos=new Vec3d[4];
	for(int j=0;j<4;++j)
	{
		global_idx[j]=simulation_mesh_->getVertexIndex(ele,j);
		vert_pos[j]=*simulation_mesh_->getVertex(global_idx[j]);
	}
	Mat3d Dm(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
	{
		Dm[i][0]=vert_pos[0][i]-vert_pos[3][i];
		Dm[i][1]=vert_pos[1][i]-vert_pos[3][i];
		Dm[i][2]=vert_pos[2][i]-vert_pos[3][i];
	}
	Mat3d DmInv=inv(Dm);
	delete[] global_idx;
	delete[] vert_pos;
	return DmInv;
}
Mat3d RealTimeExampleBasedDeformer::computeF(const int &cubica_idx,const Vec3d *reduced_dis) const
{//for all cubica elements
	int ele=object_LB_cubica_elements_[cubica_idx];
	//compute displacement x=X+Uq
	double *temp=new double[12];
	memset(temp,0.0,sizeof(double)*12);
	for(int i=0;i<4;++i)
	{
		for(int j=0;j<interpolate_eigenfunction_num_;++j)
		{
			temp[3*i+0]+=cubica_LB_tetsubBasis_[cubica_idx][i][j]*reduced_dis[j][0];
			temp[3*i+1]+=cubica_LB_tetsubBasis_[cubica_idx][i][j]*reduced_dis[j][1];
			temp[3*i+2]+=cubica_LB_tetsubBasis_[cubica_idx][i][j]*reduced_dis[j][2];
		}
	}
	memset(deformed_,0.0,sizeof(double)*12);
	for(int i=0;i<4;++i)
	{
		int vertID=simulation_mesh_->getVertexIndex(ele,i);
		deformed_[3*i]=LB_restpos_[cubica_idx][3*i]+temp[3*i];
		deformed_[3*i+1]=LB_restpos_[cubica_idx][3*i+1]+temp[3*i+1];
		deformed_[3*i+2]=LB_restpos_[cubica_idx][3*i+2]+temp[3*i+2];
	}
	Mat3d Ds=computeDs(deformed_);
	Mat3d F=Ds*computeDmInv(ele);
	delete[] temp;
	return F;
}
Mat3d RealTimeExampleBasedDeformer::computeF_gradient(const int &ele,const int &vert_idx,const int &vert_idx_dim) const
{
	//vert_idx denotes the j-th vertex, vert_idx_dim is the k-th coordinate of the j-th vertex
	//we get the result as dF_i/dx_j^k, which is the derivative force of vertex i to the vertex j on the coordinate k
	//identity vector definition
	if((vert_idx>4)||(vert_idx_dim>3))
	{
		std::cout<<"the vert_idx or the vert_idx_dim is out of range, they should be smaller than 3";
	}
	Mat3d result_matrix(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	vector<Vec3d> e_vector;
	e_vector.resize(3,Vec3d(0.0,0.0,0.0));
	e_vector[0][0]=1.0;
	e_vector[1][1]=1.0;
	e_vector[2][2]=1.0;
	Mat3d e_matrix(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	//for the j=1,2,3 we have dF/dx_j^k=e_k*e_j^T*Dm^-1
	for(unsigned int row_idx=0;row_idx<3;++row_idx)
	{
		for(unsigned int col_idx=0;col_idx<3;++col_idx)
		{
			e_matrix[row_idx][col_idx]=e_vector[vert_idx_dim][row_idx]*e_vector[vert_idx][col_idx];
		}
	}
	//compute dF/dx_4
	if(vert_idx==3)
	{
		Mat3d vert_cord_matrix(0.0);
		if(vert_idx_dim==0)
		{
			vert_cord_matrix[0][0]=vert_cord_matrix[0][1]=vert_cord_matrix[0][2]=-1.0;
		}
		else if(vert_idx_dim==1)
		{
			vert_cord_matrix[1][0]=vert_cord_matrix[1][1]=vert_cord_matrix[1][2]=-1.0;
		}
		else if(vert_idx_dim==2)
		{
			vert_cord_matrix[2][0]=vert_cord_matrix[2][1]=vert_cord_matrix[2][2]=-1.0;
		}
		else
		{
		}
		result_matrix=vert_cord_matrix*computeDmInv(ele);
	}
	else
	{
		//compute dF/dx_j^k,j=1,2,3;k=1,2,3
		result_matrix=e_matrix*computeDmInv(ele);
	}
	return result_matrix;
}

Mat3d RealTimeExampleBasedDeformer::firstPiolaKirchhoff(Mat3d &F) const
{//3*3
	//reduced P for neohookean model
	//Mat3d P=mu_*(F-trans(inv(F)))+lamda_*log(det(F))*trans(inv(F));
	// P for stvk model
    Mat3d I(1.0);
    Mat3d temp=trans(F)*F;
	Mat3d E=0.5*(temp-I);
    double trace_E=E[0][0]+E[1][1]+E[2][2];
    Mat3d P=F*(2*mu_*E+lamda_*trace_E*I);
	return P;
}
Mat3d RealTimeExampleBasedDeformer::computeP_gradient(const int &ele,const Mat3d &F,const int &vert_idx,
														const int &vert_idx_dim) const
{
	//P_gradient for neohookean model
	// Mat3d pFpx=computeF_gradient(ele,vert_idx,vert_idx_dim);
	// Mat3d temp=inv(F)*pFpx;
	// double trace=temp[0][0]+temp[1][1]+temp[2][2];
	// Mat3d tiF=trans(inv(F));
	// Mat3d pPpx=mu_*pFpx+(mu_-lamda_*log(det(F)))*tiF*trans(pFpx)*tiF+lamda_*trace*tiF;
	//P_gradient for stvk model
	Mat3d I(1.0);
	Mat3d pFpx=computeF_gradient(ele,vert_idx,vert_idx_dim);
    Mat3d temp=trans(F)*F;
	Mat3d E=0.5*(temp-I);
	double trace_E=E[0][0]+E[1][1]+E[2][2];
    Mat3d temp1=trans(pFpx)*F;
    Mat3d temp2=trans(F)*pFpx;
    Mat3d pEpF=0.5*(temp1+temp2);
    double trace_pE=pEpF[0][0]+pEpF[1][1]+pEpF[2][2];
    Mat3d temp3=pFpx*(2*mu_*E+lamda_*trace_E*I);
    Mat3d temp4=F*(2*mu_*pEpF+lamda_*trace_pE*I);
    Mat3d pPpF=temp3+temp4;
	return pPpF;
}

// void RealTimeExampleBasedDeformer::elasticTensor(const double &youngs_modulus, const double &possion_rate) const
// {
// 	elastic_tensor_(0,1) = elastic_tensor_(1,0) = elastic_tensor_(0,2) = elastic_tensor_(2,0) = possion_rate;
// 	elastic_tensor_(1,2) = elastic_tensor_(2,1) = possion_rate;
// 	elastic_tensor_(0,0) = elastic_tensor_(1,1) = elastic_tensor_(2,2) = 1.0-possion_rate;
// 	elastic_tensor_(3,3) = elastic_tensor_(4,4) = elastic_tensor_(5,5) = 0.5*(1.0-2.0*possion_rate);
// 	Scalar coef=youngs_modulus/((1.0+possion_rate)*(1.0-2.0*possion_rate));
// 	elastic_tensor_ *= coef;
// }
// void RealTimeExampleBasedDeformer::computeLinearEnergy(const Vec3d *reduced_dis,double &energy) const
// {
// 	energy=0.0;
// 	for(int cubica_idx=0;cubica_idx<object_LB_cubica_ele_num_;++cubica_idx)
// 	{
// 		Mat3d F=computeF(cubica_idx,reduced_dis);
//         Mat3d I(1.0);
//         Mat3d temp=trans(F)*F;
//     	Mat3d E=0.5*(temp-I);
// 		Matrix<double> E_mat(6,1,0.0);
// 		for(unsigned int i=0;i<3;++i)
// 		{
// 			E_mat(i,0)=E[i][i];
// 		}
// 		E_mat(3,0)=2.0*E[1][2];
// 		E_mat(4,0)=2.0*E[0][2];
// 		E_mat(5,0)=2.0*E[0][1];
// 		Matrix<double> temp=E_mat.MultiplyT(elastic_tensor_);
//         Matrix<double> temp1=temp*E_mat;
//         double element_energy=0.5*temp1;
// 		energy += object_LB_cubica_weights_[cubica_idx]*element_energy;
// 	}
// }
// void RealTimeExampleBasedDeformer::computeLinearEnergyGradient(const Vec3d *reduced_dis,double *energy_grad) const
// {
// 	memset(energy_grad,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
// 	for(int cubica_idx=0;cubica_idx<object_LB_cubica_ele_num_;++cubica_idx)
// 	{
// 		int ele=object_LB_cubica_elements_[cubica_idx];
// 		Mat3d F=computeF(cubica_idx,reduced_dis);
// 		Mat3d P=firstPiolaKirchhoff(F);
// 		Mat3d temp=P*trans(computeDmInv(ele));
//
// 		double *temp1=new double[12];
// 		memset(temp1,0.0,sizeof(double)*12);
// 		for(int j=0;j<4;++j)
// 		{
// 			for(int i=0;i<3;++i)
// 			{
// 				if(j==3)
// 					temp1[3*j+i]=(-1.0)*(temp[i][0]+temp[i][1]+temp[i][2]);
// 				else
// 					temp1[3*j+i]=temp[i][j];
// 			}
// 		}
// 		double *g=new double[3*interpolate_eigenfunction_num_];//3rx1
// 		memset(g,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
// 		for(int i=0;i<interpolate_eigenfunction_num_;++i)
// 		{
// 			for(int j=0;j<4;++j)
// 			{
// 				g[3*i]+=cubica_LB_tetsubBasis_[cubica_idx][j][i]*temp1[3*j];
// 				g[3*i+1]+=cubica_LB_tetsubBasis_[cubica_idx][j][i]*temp1[3*j+1];
// 				g[3*i+2]+=cubica_LB_tetsubBasis_[cubica_idx][j][i]*temp1[3*j+2];
// 			}
// 		}
// 		delete[] temp1;
// 		for(int i=0;i<3*interpolate_eigenfunction_num_;++i)
// 			energy_grad[i] += object_LB_cubica_weights_[cubica_idx]*g[i];
// 		delete[] g;
// 	}
// }

void RealTimeExampleBasedDeformer::computeReducedEnergy(const Vec3d *reduced_dis,double &energy) const
{
	energy=0.0;
	// PerformanceCounter counter1;
	// counter1.StartCounter();
	//energy for neohookean model
	// for(int cubica_idx=0;cubica_idx<object_LB_cubica_ele_num_;++cubica_idx)
	// {
	// 	Mat3d F=computeF(cubica_idx,reduced_dis);
	// 	Mat3d temp=trans(F)*F;
	// 	double trace_c=temp[0][0]+temp[1][1]+temp[2][2];
	// 	double lnJ=log(det(F));
	// 	double element_energy=0.5*mu_*(trace_c-3)-mu_*lnJ+0.5*lamda_*lnJ*lnJ;
	// 	energy += object_LB_cubica_weights_[cubica_idx]*element_energy;
	// }
	//energy for reduced stvk model
	for(int cubica_idx=0;cubica_idx<object_LB_cubica_ele_num_;++cubica_idx)
	{
		Mat3d F=computeF(cubica_idx,reduced_dis);
        Mat3d I(1.0);
        Mat3d temp=trans(F)*F;
    	Mat3d E=0.5*(temp-I);
        double trace_E=E[0][0]+E[1][1]+E[2][2];
        double doubleE=0.0;
        for(int i=0;i<3;++i)
            for(int j=0;j<3;++j)
                doubleE+=E[i][j]*E[i][j];
        double element_energy=0.5*lamda_*trace_E*trace_E+mu_*doubleE;
		energy += object_LB_cubica_weights_[cubica_idx]*element_energy;
	}
	// counter1.StopCounter();
	std::cout<<"energy---:"<<energy<<"\n";
}
void RealTimeExampleBasedDeformer::computeReducedInternalForce(const Vec3d *reduced_dis,double *forces) const
{//r*1
	// PerformanceCounter p_counter;
	// p_counter.StartCounter();
	memset(forces,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);

	// PerformanceCounter counter1;
	// counter1.StartCounter();
	for(int cubica_idx=0;cubica_idx<object_LB_cubica_ele_num_;++cubica_idx)
	{
		// PerformanceCounter counter4;
		// counter4.StartCounter();
		int ele=object_LB_cubica_elements_[cubica_idx];
		Mat3d F=computeF(cubica_idx,reduced_dis);
		Mat3d P=firstPiolaKirchhoff(F);
		Mat3d temp=P*trans(computeDmInv(ele));

		double *temp1=new double[12];
		memset(temp1,0.0,sizeof(double)*12);
		for(int j=0;j<4;++j)
		{
			for(int i=0;i<3;++i)
			{
				if(j==3)
					temp1[3*j+i]=(-1.0)*(temp[i][0]+temp[i][1]+temp[i][2]);
				else
					temp1[3*j+i]=temp[i][j];
			}
		}
		double *g=new double[3*interpolate_eigenfunction_num_];//3rx1
		memset(g,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
		// counter4.StopCounter();
		// std::cout<<"before-multi:"<<counter4.GetElapsedTime()<<"\n";
		// PerformanceCounter counter3;
		// counter3.StartCounter();
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
		{
			for(int j=0;j<4;++j)
			{
				g[3*i]+=cubica_LB_tetsubBasis_[cubica_idx][j][i]*temp1[3*j];
				g[3*i+1]+=cubica_LB_tetsubBasis_[cubica_idx][j][i]*temp1[3*j+1];
				g[3*i+2]+=cubica_LB_tetsubBasis_[cubica_idx][j][i]*temp1[3*j+2];
			}
		}
		delete[] temp1;
		// counter3.StopCounter();
		// std::cout<<"multi:"<<counter3.GetElapsedTime()<<"\n";
		// PerformanceCounter counter5;
		// counter5.StartCounter();
		for(int i=0;i<3*interpolate_eigenfunction_num_;++i)
			forces[i] += object_LB_cubica_weights_[cubica_idx]*g[i];
		// counter4.StopCounter();
		// std::cout<<"after multi:"<<counter4.GetElapsedTime()<<"\n";
		delete[] g;
	}
	// counter1.StopCounter();
	// std::cout<<"internal force:"<<counter1.GetElapsedTime()<<"\n";
	// p_counter.StopCounter();
	// std::cout<<"compute internal force:"<<p_counter.GetElapsedTime()<<"\n";
}
void RealTimeExampleBasedDeformer::testEnergyGradients()
{
	//test energy and energy gradient-internal force
	// std::cout<<"--a\n";
	real_1d_array x;
	int num=3*interpolate_eigenfunction_num_;
	double *temp_buffer=new double[num];
	srand((unsigned)time(0));
	int lowest=1,highest=10;
	int range=(highest-lowest)+1;
	// std::cout<<"--b\n";
	//randomly generate initial x in range[0.1,1]
	for(int i=0;i<num;++i)
		temp_buffer[i]=(lowest+rand()%range)/100.0;
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	for(int j=0;j<3;++j)
	// 		temp_buffer[3*i+j]=object_eigencoefs_[i][j];
	// std::cout<<"--c\n";
	x.setcontent(num,temp_buffer);
	double f_min,f_plus,f;
	real_1d_array temp_grad,grad;
	temp_grad.setcontent(num,temp_buffer);
	grad.setcontent(num,temp_buffer);
	double target_weights[2]={0.5,0.5};
	// std::cout<<"--d\n";
	evaluateObjectiveAndGradient2(x,f,grad,(void*)target_weights);
	// evaluateObjectiveAndGradient3(x,f,grad,(void*)target_weights);
	for(int i=0;i<num;++i)
	{
		temp_buffer[i]=(lowest+rand()%range)/1.0e8;
		x[i]+=temp_buffer[i];
	}
	// std::cout<<"--e\n";
	evaluateObjectiveAndGradient2(x,f_plus,temp_grad,(void*)target_weights);
	// evaluateObjectiveAndGradient3(x,f_plus,temp_grad,(void*)target_weights);
	for(int i=0;i<num;++i)
	{
		x[i]-=2*temp_buffer[i];
	}
	evaluateObjectiveAndGradient2(x,f_min,temp_grad,(void*)target_weights);
	// evaluateObjectiveAndGradient3(x,f_min,temp_grad,(void*)target_weights);
	double f_grad_times_dx=0.0;
	for(int i=0;i<num;++i)
		f_grad_times_dx+=grad[i]*2*temp_buffer[i];
		std::cout<<f_plus<<"...\n";
	std::cout<<"Objective, df, analytic: "<<setprecision(15)<<f_grad_times_dx<<", numerial: "<<setprecision(15)<<f_plus-f_min;
	std::cout<<", rel_error= "<<(f_plus-f_min-f_grad_times_dx)/(fabs(f_grad_times_dx)>1e-20?fabs(f_grad_times_dx):1e-20)<<"\n";
	delete[] temp_buffer;
}
// void RealTimeExampleBasedDeformer::testObjectiveGradients()
// {
// 	//test internal force and force gradient--stiffness matrix
// 	int num=3*interpolate_eigenfunction_num_;
// 	Vec3d *q=new Vec3d[interpolate_eigenfunction_num_];
// 	//double *x=new double[num];
// 	double *dx=new double[num];
// 	double *f=new double[num];
// 	double *f_plus=new double[num];
// 	double *f_min=new double[num];
// 	Matrix<double> K((int)interpolate_eigenfunction_num_,(int)interpolate_eigenfunction_num_);
// 	double *dK=new double[num];
// 	srand((unsigned)time(0));
// 	//randomly generate displacement in range [0.1,1]
// 	double x[24]={0.3,0.1,0.4,0.3,0.8,1,0.5,1,0.7,0.2,0.2,0.3,0.2,0.4,0.6,0.2,0.7,0.3,0.8,0.6,1,1,0.3,0.5};
// 	//double x[24]={0.4,0.4,0.4,0.3,0.3,0.3,0.5,0.5,0.5,0.2,0.2,0.2,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.5,0.5,0.5};
// 	for(int i=0;i<num;++i)
// 	{
// 		x[i]=(1+rand()%10)/10.0;
// 		std::cout<<x[i]<<",";
// 	}
// 	for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 	{
// 		for(int k=0;k<3;++k)
// 		{
// 			q[j][k]=x[3*j+k];
// 		}
// 	}
// 	computeReducedStiffnessMatrix(q,K);
// 	//
// 	std::cout<<"K:\n";
// 	for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 	{
// 		for(int k=0;k<interpolate_eigenfunction_num_;++k)
// 		{
// 			std::cout<<K(j,k)<<",";
// 		}
// 		std::cout<<"\n";
// 	}
// 	//perturb a little bit
// 	for(int i=0;i<num;++i)
// 	{
// 		dx[i]=(1+rand()%10)/1.0e7;
// 		x[i]+=dx[i];
// 	}
// 	for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 	{
// 		for(int k=0;k<3;++k)
// 		{
// 			q[j][k]=x[3*j+k];
// 		}
// 	}
// 	computeReducedInternalForce(q,f_plus);
// 	for(int i=0;i<num;++i)
// 	{
// 		//dx[i]*=1.0;
// 		x[i]-=2*dx[i];
// 	}
// 	for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 	{
// 		for(int k=0;k<3;++k)
// 		{
// 			q[j][k]=x[3*j+k];
// 		}
// 	}
// 	computeReducedInternalForce(q,f_min);
// 	Matrix<double> dis_matrix((int)interpolate_eigenfunction_num_,3);
// 	for(int i=0;i<interpolate_eigenfunction_num_;++i)
// 		for(int j=0;j<3;++j)
// 			dis_matrix(i,j)=2.0*dx[3*i+j];
// 	Matrix<double> temp_matrix=K*dis_matrix;
// 	for(int i=0;i<interpolate_eigenfunction_num_;++i)
// 	{
// 		for(int j=0;j<3;++j)
// 			dK[3*i+j]=temp_matrix(i,j);
// 	}
// 	double max_rel_error=1000;
// 	double df,dk;
// 	for(int i=0;i<3*interpolate_eigenfunction_num_;++i)
// 	{
// 		double rel_error = fabs((f_plus[i]-f_min[i]-dK[i])/dK[i]);
// 		//std::cout<<"f_plus:"<<f_plus[i]<<"f_min:"<<f_min[i]<<"dk:"<<dK[i]<<"\n";
// 		if(rel_error < max_rel_error)
// 		{
// 			df=f_plus[i]-f_min[i];
// 			dk = dK[i];
// 			max_rel_error = rel_error;
// 		}
// 	}
// 	std::cout<<"Max relative error:\n";
// 	std::cout<<"df: "<<setprecision(15)<<df<<"\n";
// 	std::cout<<"k*dk: "<<setprecision(15)<<dk<<"\n";
// 	std::cout<<"rel error: "<<max_rel_error<<"\n";
// 	//delete[] x;
// 	delete[] dx;
// 	delete[] f;
// 	delete[] f_plus;
// 	delete[] f_min;
// 	delete[] dK;
// }
int RealTimeExampleBasedDeformer::ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const
{
    // The code handles the following necessary special situations (see the code below) :

    //---------------------------------------------------------
    // 1. det(V) == -1
    //    - simply multiply a column of V by -1
    //---------------------------------------------------------
    // 2. An entry of Fhat is near zero
    //---------------------------------------------------------
    // 3. Tet is inverted.
    //    - check if det(U) == -1
    //    - If yes, then negate the minimal element of Fhat
    //      and the corresponding column of U
    //---------------------------------------------------------

    double modifiedSVD_singularValue_eps = 1e-8;

    // form F^T F and do eigendecomposition
    Mat3d normalEq = trans(F) * F;
    Vec3d eigenValues;
    Vec3d eigenVectors[3];

    // note that normalEq is changed after calling eigen_sym
    eigen_sym(normalEq, eigenValues, eigenVectors);

    V.set(eigenVectors[0][0], eigenVectors[1][0], eigenVectors[2][0],
          eigenVectors[0][1], eigenVectors[1][1], eigenVectors[2][1],
          eigenVectors[0][2], eigenVectors[1][2], eigenVectors[2][2]);
    /*
      printf("--- original V ---\n");
      V.print();
      printf("--- eigenValues ---\n");
      printf("%G %G %G\n", eigenValues[0], eigenValues[1], eigenValues[2]);
    */

    // Handle situation:
    // 1. det(V) == -1
    //    - simply multiply a column of V by -1
    if (det(V) < 0.0)
    {
        // convert V into a rotation (multiply column 1 by -1)
        V[0][0] *= -1.0;
        V[1][0] *= -1.0;
        V[2][0] *= -1.0;
    }

    Fhat[0] = (eigenValues[0] > 0.0) ? sqrt(eigenValues[0]) : 0.0;
    Fhat[1] = (eigenValues[1] > 0.0) ? sqrt(eigenValues[1]) : 0.0;
    Fhat[2] = (eigenValues[2] > 0.0) ? sqrt(eigenValues[2]) : 0.0;

    //printf("--- Fhat ---\n");
    //printf("%G %G %G\n", Fhat[0][0], Fhat[1][1], Fhat[2][2]);

    // compute inverse of singular values
    // also check if singular values are close to zero
    Vec3d FhatInverse;
    FhatInverse[0] = (Fhat[0] > modifiedSVD_singularValue_eps) ? (1.0 / Fhat[0]) : 0.0;
    FhatInverse[1] = (Fhat[1] > modifiedSVD_singularValue_eps) ? (1.0 / Fhat[1]) : 0.0;
    FhatInverse[2] = (Fhat[2] > modifiedSVD_singularValue_eps) ? (1.0 / Fhat[2]) : 0.0;

    // compute U using the formula:
    // U = F * V * diag(FhatInverse)
    U = F * V;
    U.multiplyDiagRight(FhatInverse);

    // In theory, U is now orthonormal, U^T U = U U^T = I .. it may be a rotation or a reflection, depending on F.
    // But in practice, if singular values are small or zero, it may not be orthonormal, so we need to fix it.
    // Handle situation:
    // 2. An entry of Fhat is near zero
    // ---------------------------------------------------------

    /*
      printf("--- FhatInverse ---\n");
      FhatInverse.print();
      printf(" --- U ---\n");
      U.print();
    */

    if ((Fhat[0] < modifiedSVD_singularValue_eps) && (Fhat[1] < modifiedSVD_singularValue_eps) && (Fhat[2] < modifiedSVD_singularValue_eps))
    {
        // extreme case, material has collapsed almost to a point
        // see [Irving 04], p. 4
        U.set(1.0, 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, 1.0);
    }
    else
    {
        int done = 0;
        for(int dim=0; dim<3; dim++)
        {
            int dimA = dim;
            int dimB = (dim + 1) % 3;
            int dimC = (dim + 2) % 3;
            if ((Fhat[dimB] < modifiedSVD_singularValue_eps) && (Fhat[dimC] < modifiedSVD_singularValue_eps))
            {
                // only the column dimA can be trusted, columns dimB and dimC correspond to tiny singular values
                Vec3d tmpVec1(U[0][dimA], U[1][dimA], U[2][dimA]); // column dimA
                Vec3d tmpVec2;
                FindOrthonormalVector(tmpVec1, tmpVec2);
                Vec3d tmpVec3 = norm(cross(tmpVec1, tmpVec2));
                U[0][dimB] = tmpVec2[0];
                U[1][dimB] = tmpVec2[1];
                U[2][dimB] = tmpVec2[2];
                U[0][dimC] = tmpVec3[0];
                U[1][dimC] = tmpVec3[1];
                U[2][dimC] = tmpVec3[2];
                if (det(U) < 0.0)
                {
                    U[0][dimB] *= -1.0;
                    U[1][dimB] *= -1.0;
                    U[2][dimB] *= -1.0;
                }
                done = 1;
                break; // out of for
            }
        }

        if (!done)
        {
            for(int dim=0; dim<3; dim++)
            {
                int dimA = dim;
                int dimB = (dim + 1) % 3;
                int dimC = (dim + 2) % 3;

                if (Fhat[dimA] < modifiedSVD_singularValue_eps)
                {
                    // columns dimB and dimC are both good, but column dimA corresponds to a tiny singular value
                    Vec3d tmpVec1(U[0][dimB], U[1][dimB], U[2][dimB]); // column dimB
                    Vec3d tmpVec2(U[0][dimC], U[1][dimC], U[2][dimC]); // column dimC
                    Vec3d tmpVec3 = norm(cross(tmpVec1, tmpVec2));
                    U[0][dimA] = tmpVec3[0];
                    U[1][dimA] = tmpVec3[1];
                    U[2][dimA] = tmpVec3[2];
                    if (det(U) < 0.0)
                    {
                        U[0][dimA] *= -1.0;
                        U[1][dimA] *= -1.0;
                        U[2][dimA] *= -1.0;
                    }
                    done = 1;
                    break; // out of for
                }
            }
        }

        if (!done)
        {
            // Handle situation:
            // 3. Tet is inverted.
            //    - check if det(U) == -1
            //    - If yes, then negate the minimal element of Fhat
            //      and the corresponding column of U

            double detU = det(U);
            if (detU < 0.0)
            {
                // tet is inverted
                // find smallest singular value (they are all non-negative)
                int smallestSingularValueIndex = 0;
                for(int dim=1; dim<3; dim++)
                    if (Fhat[dim] < Fhat[smallestSingularValueIndex])
                        smallestSingularValueIndex = dim;

                // negate smallest singular value
                Fhat[smallestSingularValueIndex] *= -1.0;
                U[0][smallestSingularValueIndex] *= -1.0;
                U[1][smallestSingularValueIndex] *= -1.0;
                U[2][smallestSingularValueIndex] *= -1.0;
            }
        }
    }

    /*
      printf("U = \n");
      U.print();
      printf("Fhat = \n");
      Fhat.print();
      printf("V = \n");
      V.print();
    */

    return 0;
}

void RealTimeExampleBasedDeformer::FindOrthonormalVector(Vec3d & v, Vec3d & result) const
{
    // find smallest abs component of v
    int smallestIndex = 0;
    for(int dim=1; dim<3; dim++)
        if (fabs(v[dim]) < fabs(v[smallestIndex]))
            smallestIndex = dim;

    Vec3d axis(0.0, 0.0, 0.0);
    axis[smallestIndex] = 1.0;

    // this cross-product will be non-zero (as long as v is not zero)
    result = norm(cross(v, axis));
}
Mat3d RealTimeExampleBasedDeformer::computeDeformationGradient(const Mat3d &init_matrix,const Mat3d &deformed_matrix/*Vec3d *init_pos,Vec3d *deform_pos*/)
{
	Mat3d result(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	result=deformed_matrix*inv(init_matrix);

    //handle inversion
    Mat3d U,V;
    Vec3d Fhat;
    ModifiedSVD(result,U,Fhat,V);
    //clamphat if below the principle stretch threshold
    double principle_threshold = 0.6;
    for(unsigned int i = 0; i < 3 ; ++i)
        if(Fhat[i] < principle_threshold)
            Fhat[i] = principle_threshold;
    Mat3d Fhat_mat;
    for(unsigned int i = 0; i < 3; ++i)
        for(unsigned int j = 0; j < 3; ++j)
            Fhat_mat[i][j] = (i==j)?Fhat[i]:0;
    result = U*Fhat_mat*trans(V);
    return result;
}
}  //namespace RTLB
