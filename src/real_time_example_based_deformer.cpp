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
	// getchar();
}
void RealTimeExampleBasedDeformer::setu(double *u)
{
	memcpy(u_,u,sizeof(double)*3*simulation_vertices_num_);
}
void RealTimeExampleBasedDeformer::setupSimulation()
{
	std::cout<<"setupSimulation begins\n";
	if(!simulation_mesh_)
	{
		std::cout<<"simulation mesh is null.\n";
		exit(0);
	}

    mesh_graph_=GenerateMeshGraph::Generate(simulation_mesh_);
    int scale_rows=1;
    mesh_graph_->GetLaplacian(&laplacian_matrix_,scale_rows);
    mesh_graph_->GetLaplacian(&laplacian_damping_matrix_,scale_rows);
    laplacian_damping_matrix_->ScalarMultiply(damping_laplacian_coef_);
	// simulation_vertices_num_=simulation_mesh_->getNumVertices();
	u_=new double[3*simulation_vertices_num_];
	memset(u_,0.0,sizeof(double)*3*simulation_vertices_num_);
	vel_=new double[3*simulation_vertices_num_];
	memset(vel_,0.0,sizeof(double)*3*simulation_vertices_num_);
	u_initial_=new double[3*simulation_vertices_num_];
	memset(u_initial_,0.0,sizeof(double)*3*simulation_vertices_num_);
	vel_initial_=new double[3*simulation_vertices_num_];
	memset(vel_initial_,0.0,sizeof(double)*3*simulation_vertices_num_);
	f_ext_=new double[3*simulation_vertices_num_];
	memset(f_ext_,0.0,sizeof(double)*3*simulation_vertices_num_);
	f_col_=new double[3*simulation_vertices_num_];
	memset(f_col_,0.0,sizeof(double)*3*simulation_vertices_num_);
	// full_drag_force_=new double[3*simulation_vertices_num_];
	// memset(full_drag_force_,0.0,sizeof(double)*3*simulation_vertices_num_);

	//reduced space
	// reduced_drag_force_=new double[r_];
	// memset(reduced_drag_force_,0.0,sizeof(double)*r_);
	std::cout<<"r:------------"<<r_<<"\n";
	q_=new double[r_];
	memset(q_,0.0,sizeof(double)*r_);
	qvel_=new double[r_];
	memset(qvel_,0.0,sizeof(double)*r_);
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
	// //temp compute mass for each vertex
	// for(int el=0;el<simulation_mesh_->getNumElements();++el)
	// {
	// 	for(int i=0;i<4;++i)
	// 	{
	// 		int vertID=simulation_mesh_->getVertexIndex(el,i);
	// 		mass_[vertID]+=0.25*simulation_mesh_->getElementVolume(el)*simulation_mesh_->getElementDensity(el)/3.0;
	// 	}
	// }
	// std::cout<<"mass:\n";
	// for(int i=0;i<20;++i)
	// 	std::cout<<i<<":"<<mass_[i]<<";\n";
	// 	getchar();
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
		isotropic_material_=new NeoHookeanIsotropicMaterial(tet_mesh_);
		std::cout<<"Invertible material: neo-Hookean.\n";
		isotropic_hyperelastic_fem_=new IsotropicHyperelasticFEM(tet_mesh_,isotropic_material_,principal_stretch_threshold_,add_gravity_,gravity_);
		force_model_=new IsotropicHyperelasticFEMForceModel(isotropic_hyperelastic_fem_);
	}
	else
	{
		reduced_neoHookean_force_model_ = new ReducedNeoHookeanForceModel(r_,simulation_mesh_,U_,object_cubica_ele_num_,
									object_cubica_weights_,object_cubica_elements_,restpos_,add_gravity_,gravity_);
		reduced_force_model_=reduced_neoHookean_force_model_;
		if(reduced_force_model_)
			std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!reduced force model is not null\n";
	}
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
    else if(solver_type_==CENTRALDIFFERENCES)
    {
        integrator_base_sparse_=new CentralDifferencesSparse(3*simulation_vertices_num_,time_step_,mass_matrix_,force_model_,
                                                            fixed_dofs_num_,fixed_dofs_,damping_mass_coef_,damping_stiffness_coef_,
                                                            central_difference_tangential_damping_update_mode_,solver_threads_num_);
    }
    else if(solver_type_==REDUCEDCENTRALDIFFERENCES)
    {
        central_differences_dense_ = new CentralDifferencesDense(r_,time_step_,reduced_mass_matrix_,reduced_neoHookean_force_model_,damping_mass_coef_,
                                                            damping_stiffness_coef_,central_difference_tangential_damping_update_mode_);
        integrator_base_dense_=central_differences_dense_;
    }
    else if(solver_type_==REDUCEDIMPLICITNEWMARK)
    {
        implicit_newmark_dense_ = new ImplicitNewmarkDense(r_,time_step_,reduced_mass_matrix_,reduced_force_model_,
                                                        ImplicitNewmarkDense::positiveDefiniteMatrixSolver,damping_mass_coef_,damping_stiffness_coef_,
                                                        max_iterations_,integrator_epsilon_,newmark_beta_,newmark_gamma_);
        integrator_base_dense_=implicit_newmark_dense_;
		std::cout<<"reducedImplicitNewmark\n";
		if(integrator_base_dense_)
			std::cout<<"...............integrator base dense is not null\n";
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
	std::cout<<"--------------------------\n";
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
		// std::cout<<"a\n";
		// for(int i=0;i<r_;++i)
		// {
		// 	std::cout<<q_[i]<<","<<qvel_[i]<<",--";
		// }
		// 	std::cout<<"\n";
		// if(integrator_base_dense_)
		// 	std::cout<<"...............integrator base is not null\n";
		// implicit_newmark_dense_->SetState(q_,qvel_);
		// std::cout<<"b\n";
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
	std::cout<<"init example simulation\n";
	//initial for example simulation
	target_eigencoefs_=new Vec3d[interpolate_eigenfunction_num_];
	memset(target_eigencoefs_,0.0,sizeof(Vec3d)*interpolate_eigenfunction_num_);
	target_eigencoefs_diff_=new Vec3d[interpolate_eigenfunction_num_];
	memset(target_eigencoefs_diff_,0.0,sizeof(Vec3d)*interpolate_eigenfunction_num_);
	object_current_eigencoefs_ = new Vec3d[interpolate_eigenfunction_num_];
	memset(object_current_eigencoefs_,0.0,sizeof(Vec3d)*interpolate_eigenfunction_num_);
	example_guided_deformation_ = new double[3*simulation_vertices_num_];
	memset(example_guided_deformation_,0.0,sizeof(double)*3*simulation_vertices_num_);
	target_deformation_ = new double[3*simulation_vertices_num_];
	memset(target_deformation_,0.0,sizeof(double)*3*simulation_vertices_num_);
	example_based_LB_forces_ = new double[3*interpolate_eigenfunction_num_];
	memset(example_based_LB_forces_,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
	example_based_forces_ = new double[3*simulation_vertices_num_];
	memset(example_based_forces_,0.0,sizeof(double)*3*simulation_vertices_num_);
	example_based_q_ = new double[r_];
	memset(example_based_q_,0.0,sizeof(double)*r_);
	example_based_fq_ = new double[r_];
	memset(example_based_fq_,0.0,sizeof(double)*r_);


	// std::cout<<"reduced:testEnergyGradients:\n";
	// reduced_neoHookean_force_model_->testEnergyGradients();
	// getchar();
	// std::cout<<"testInternalForceGradients:\n";
	// reduced_neoHookean_force_model_->testObjectiveGradients();
	// std::cout<<"testInternalForceGradients-end:\n";
	// getchar();
}

void RealTimeExampleBasedDeformer::advanceStep()
{
	if(simulation_mode_==FULLSPACE)
		fullspaceSimulation();
	else
		reducedspaceSimulation();
}
void RealTimeExampleBasedDeformer::setExternalForces(double *ext_forces)
{
	memcpy(f_ext_,ext_forces,sizeof(double)*3*simulation_mesh_->getNumVertices());
}

void RealTimeExampleBasedDeformer::setReducedExternalForces(double *ext_forces)
{
	memcpy(fq_,ext_forces,sizeof(double)*r_);
}
void RealTimeExampleBasedDeformer::fullspaceSimulation()
{
	// std::cout<<"fullspaceSimulation:\n";
	// memcpy(f_ext_,full_drag_force_,sizeof(double)*3*simulation_vertices_num_);
	// if(time_step_counter_<force_loads_num_)
	// {
	// 	std::cout<<"External forces read from the text input file.\n";
	// 	for(int i=0;i<3*simulation_vertices_num_;++i)
	// 		f_ext_[i]+=force_loads_[ELT(3*active_instance->simulation_vertices_num_,i,active_instance->time_step_counter_)];
	// }
	//apply the force loads caused by the examples
	if(enable_example_simulation_)
	{
		projectOnEigenFunctions(simulation_mesh_,u_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
								interpolate_eigenfunction_num_,object_current_eigencoefs_);
		// projectReducedSubspaceOnLBSubspace(fq_,object_surface_eigencoefs_);
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			object_current_eigencoefs_[i]=object_current_eigencoefs_[i]-object_eigencoefs_[i]+example_eigencoefs_[0][i];
		memcpy(target_eigencoefs_,object_current_eigencoefs_,sizeof(Vec3d)*interpolate_eigenfunction_num_);
		//find target shape for simulation object
		projectOnExampleManifold(object_current_eigencoefs_,target_eigencoefs_);
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			target_eigencoefs_[i]=object_current_eigencoefs_[i]-target_eigencoefs_[i];
		reconstructFromEigenCoefs(target_eigencoefs_,example_guided_deformation_);

		//explicit force;
		// VolumetricMesh::Material *material=simulation_mesh_->getMaterial(0);
		// VolumetricMesh::ENuMaterial *enu_material=downcastENuMaterial(material);
		// double new_E=1.0*enu_material->getE();
		for(int i=0;i<3*simulation_vertices_num_;++i)
			f_ext_[i]-=30.0*example_guided_deformation_[i];
	}
	//apply the penalty collision forces with planes in scene
	// if(plane_num_>0)
	// {
	// 	planes_->resolveContact(simulation_mesh_,f_col_);
	// 	//active_instance->planes_->resolveContact(resolveContact(const ObjMesh *mesh/*,double *forces*/,const double *vel,double *u_new,double *vel_new))
	//
	// 	for(int i=0;i<3*simulation_vertices_num_;++i)
	// 		f_ext_[i]+=f_col_[i];
	// }
	//set forces to the integrator
	integrator_base_sparse_->SetExternalForces(f_ext_);
	//time step the dynamics
	// if(time_step_counter_ < total_steps_)
	// {
		int code=integrator_base_->DoTimestep();
		std::cout<<".";fflush(NULL);
		++time_step_counter_;
	// }
	memcpy(u_,integrator_base_->Getq(),sizeof(double)*3*(simulation_vertices_num_));

}
void RealTimeExampleBasedDeformer::reducedspaceSimulation()
{
	//apply any scripted force loads for reduced space ---not done yet
	// if(time_step_counter_<force_loads_num_)
	// {
	// 	// std::cout<<"External forces read from the text input file.\n";
	// 	for(int i=0;i<3*simulation_vertices_num_;++i)
	// 		f_ext_[i]+=force_loads_[ELT(3*simulation_vertices_num_,i,time_step_counter_)];
	// 	ProjectVector(3*simulation_vertices_num_,r_,U_,reduced_force_loads_,f_ext_);
	// 	for(int i=0;i<active_instance->r_;++i)
	// 		fq_[i] += reduced_force_loads_[i];
	// }
	if(enable_example_simulation_)
	{
		//project 3nx1 reduced subspace on mx3 LB subspace
		//first: project rx1 reduced displacement to


		PerformanceCounter counter3;
		counter3.StartCounter();
		projectOnEigenFunctions(simulation_mesh_,u_,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
								interpolate_eigenfunction_num_,object_current_eigencoefs_);
		counter3.StopCounter();
		std::cout<<"project on eigenfunctions:"<<counter3.GetElapsedTime()<<"\n";

		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			object_current_eigencoefs_[i]=object_current_eigencoefs_[i]-object_eigencoefs_[i]+example_eigencoefs_[0][i];
		memcpy(target_eigencoefs_,object_current_eigencoefs_,sizeof(Vec3d)*interpolate_eigenfunction_num_);

		PerformanceCounter counter31;
		counter31.StartCounter();
		//find target shape for simulation object
		projectOnExampleManifold(object_current_eigencoefs_,target_eigencoefs_);
		counter31.StopCounter();
		std::cout<<"find target_eigencoefs:"<<counter31.GetElapsedTime()<<"\n";
		// for(int i=0;i<interpolate_eigenfunction_num_;++i)
		// 	std::cout<<target_eigencoefs_[i]<<",\n";
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
		{
			target_eigencoefs_diff_[i]=target_eigencoefs_[i]-object_eigencoefs_[i];
			target_eigencoefs_[i]=object_current_eigencoefs_[i]-target_eigencoefs_[i];

		}

		PerformanceCounter counter32;
		counter32.StartCounter();
		reconstructFromEigenCoefs(target_eigencoefs_diff_,target_deformation_,0);//reference+target_deformation=target configuration
		reconstructFromEigenCoefs(target_eigencoefs_,example_guided_deformation_,1);
		counter32.StopCounter();
		std::cout<<"reconstruct 2 shapes:"<<counter32.GetElapsedTime()<<"\n";

		PerformanceCounter counter33;
		counter33.StartCounter();
		modal_matrix_->ProjectVector(example_guided_deformation_,example_based_q_);
		counter33.StopCounter();
		std::cout<<"project on reduced simulation:"<<counter33.GetElapsedTime()<<"\n";


		PerformanceCounter counter34;
		counter34.StartCounter();
		reduced_neoHookean_force_model_->computeReducedElasticInternalForce(example_based_q_,example_based_fq_,target_deformation_);
		counter34.StopCounter();
		std::cout<<"compute reduced elastic internal forces:"<<counter34.GetElapsedTime()<<"\n";

		for(int i=0;i<r_;++i)
			fq_[i]-=0.05*example_based_fq_[i];
		//test.................
		// reconsturctFromEigenCoefs(target_eigencoefs_diff_,target_deformation_);
		// modal_matrix_->ProjectVector(target_deformation_,target_q_);
		// reduced_neoHookean_force_model_->computeReducedInternalForce(target_q_,example_based_fq_);
		//....test
		//linear force
		// for(int i=0;i<3*simulation_vertices_num_;++i)
		// 	f_ext_[i]=(-100.0)*example_guided_deformation_[i];
		// modal_matrix_->ProjectVector(f_ext_,example_based_fq_);
		// for(int i=0;i<r_;++i)
		// 	fq_[i]+=example_based_fq_[i];
	}
	//apply the penalty collision forces with planes in scene in reduced space
	// if(plane_num_>0)
	// {
	// 	planes_->resolveContact(simulation_mesh_,f_col_);
	// 	ProjectVector(3*simulation_vertices_num_,r_,U_,fq_ext_,f_col_);
	// 	for(int i=0;i<r_;++i)
	// 		fq_[i]+=0;
	// }
	std::cout<<"1\n";
	PerformanceCounter counter1;
	counter1.StartCounter();
	integrator_base_dense_->SetExternalForces(fq_);
std::cout<<"2\n";
	counter1.StopCounter();
	// std::cout<<"integrator f_ext:"<<counter1.GetElapsedTime()<<"\n";
	// if(time_step_counter_ < active_instance->total_steps_)
	// {
	PerformanceCounter counter2;
	counter2.StartCounter();
		int code=integrator_base_dense_->DoTimestep();
		std::cout<<".";fflush(NULL);
		std::cout<<"3\n";
	counter2.StopCounter();
	// std::cout<<"integrator DoTimestep:"<<counter2.GetElapsedTime()<<"\n";
		// ++time_step_counter_;
	// }
	std::cout<<",,,,,,,,,,,,,,,,,,,,,,,,\n";
	memcpy(q_,integrator_base_->Getq(),sizeof(double)*r_);
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
	for(unsigned int ele_idx=0;ele_idx<simulation_mesh_->getNumElements();++ele_idx)
	{
		for(unsigned int i=0;i<simulation_mesh_->getNumElementVertices();++i)
		{
			unsigned int global_vert_idx=simulation_mesh_->getVertexIndex(ele_idx,i);
			object_vertex_volume_[global_vert_idx]+=simulation_mesh_->getElementVolume(ele_idx);
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
		for(unsigned int ele_idx=0;ele_idx<examples_[i]->getNumElements();++ele_idx)
		{
			for(unsigned int j=0;j<examples_[i]->getNumElementVertices();++j)
			{
				unsigned int global_idx=examples_[i]->getVertexIndex(ele_idx,j);
				example_vertex_volume_[i][global_idx]+=examples_[i]->getElementVolume(ele_idx);
			}
		}
	}
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
	//read file, save as reduced_basis_[3*vert_idx+dim][reduced_idx]
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
	std::cout<<reduced_basis_[reduced_col_num-1][reduced_row_num-1];
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
	for(int i=0;i<reconstruct_eigenfunction_num_;++i)
		std::cout<<object_eigenvalues_[i]<<",";
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
	std::cout<<"eee\n";
	unsigned int line_num=0;
	str_num=0;
	unsigned int total_str_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		++total_str_num;
		double temp_value1;
		input_file>>temp_value1;
		object_eigenfunctions_[str_num][line_num]=temp_value1;
		if(str_num>=reconstruct_eigenfunction_num_-1)
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

	// for(int i=0;i<reconstruct_eigenfunction_num_;++i)
	// 	for(int j=0;j<10;++j)
	// 		std::cout<<object_eigenfunctions_[i][j]<<",";
	//normalize eigenfunction with respect to the s-inner product: <f,f>=1 on the volumetric vertices
	for(int i=0;i<reconstruct_eigenfunction_num_;++i)
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
	object_eigencoefs_=new Vec3d[reconstruct_eigenfunction_num_];
    projectOnEigenFunctions(simulation_mesh_,dis,object_vertex_volume_,object_eigenfunctions_,object_eigenvalues_,
							reconstruct_eigenfunction_num_,object_eigencoefs_);

	delete[] dis;
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
			double temp_value1;
			input_file>>temp_value1;
			example_eigenfunctions_[ex_num][str_num][line_num]=temp_value1;
			if(str_num>=reconstruct_eigenfunction_num_-1)
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
		for(int i=0;i<reconstruct_eigenfunction_num_;++i)
		{
			double sum=0.0;
			for(int j=0;j<examples_[ex_num]->getNumVertices();++j)
				sum+=example_eigenfunctions_[ex_num][i][j]*example_eigenfunctions_[ex_num][i][j]*example_vertex_volume_[ex_num][j];
			if(sum>epsilon_)
				sum=sqrt(sum);
			else
				std::cout<<"object eigenfunction sum is 0\n";
			for(int j=0;j<examples_[ex_num]->getNumVertices();++j)
				example_eigenfunctions_[ex_num][i][j]/=sum;
		}
		example_eigencoefs_[ex_num]=new Vec3d[reconstruct_eigenfunction_num_];
        double *dis=new double[3*examples_[ex_num]->getNumVertices()];
        for(int j=0;j<3*examples_[ex_num]->getNumVertices();++j)
            dis[j]=0.0;
        projectOnEigenFunctions(examples_[ex_num],dis,example_vertex_volume_[ex_num],
                                example_eigenfunctions_[ex_num],example_eigenvalues_[ex_num],
                                reconstruct_eigenfunction_num_,example_eigencoefs_[ex_num]);

        delete[] dis;
	}
	std::cout<<"loadexampleEigenfunctions---end.\n";
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
			object_corresponding_functions_[i][j]=0;
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
				example_corresponding_functions_[i][j][k]=0;
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
		isotropic_hyperelastic_fem_->SetGravity(add_gravity);
	else
		reduced_neoHookean_force_model_->SetGravity(add_gravity,gravity,U_);

}
//the result eigencoefs is a shape not a displacement on LB space
void RealTimeExampleBasedDeformer::projectOnEigenFunctions(VolumetricMesh *mesh, double *displacement, double *vertex_volume,
                                                           double **eigenfunctions, double *eigenvalues, unsigned int eigenfunction_num,
                                                           Vec3d *eigencoefs)
{
	//local example mode is not written yet
    int vert_num=mesh->getNumVertices();
	double scale_factor;
	for(unsigned int i=0;i<eigenfunction_num;++i)
	{
		for(unsigned int j=0;j<3;++j)
		{
			eigencoefs[i][j]=0.0;
		}
		scale_factor=eigenvalues[i];
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

//input:m*1 eigencoefs (m:eigenfunction_num)
//output:3n*1 Euclidean pos (n:volumetricMesh vertices num)
//reconstruct step just used for simulation object and contains two steps
//step1: reconstruct the lower frequencies from the interpolation shapes
//step2: contains the higher frequencies from the simulation object eigencoefs
// void RealTimeExampleBasedDeformer::reconstructFromEigenCoefs(double **eigenfunctions,double *eigenvalues,Vec3d *eigencoefs,Vec3d *target_eigencoefs,
//                                                             const int &eigenfunction_num,const int &input_reconstruct_eigenfunction_num,const int &vert_num, double *vert_pos)
// {
// 	double scale_factor;
// 	for(unsigned int vert_idx=0;vert_idx<vert_num;++vert_idx)
// 	{
// 		vert_pos[3*vert_idx]=vert_pos[3*vert_idx+1]=vert_pos[3*vert_idx+2]=0.0;
// 		for(unsigned int i=0;i<20;++i)
// 		{
// 			scale_factor=1.0/eigenvalues[i];
// 			vert_pos[3*vert_idx]+=target_eigencoefs[i][0]*eigenfunctions[i][vert_idx]*scale_factor;
// 			vert_pos[3*vert_idx+1]+=target_eigencoefs[i][1]*eigenfunctions[i][vert_idx]*scale_factor;
// 			vert_pos[3*vert_idx+2]+=target_eigencoefs[i][2]*eigenfunctions[i][vert_idx]*scale_factor;
//
// 		}
// 		//reconstruct the higher spectral pose eigencoefs from the simulation object itself
// 		// for(unsigned int i=eigenfunction_num;i<input_reconstruct_eigenfunction_num;++i)
// 		// {
// 		// 	scale_factor=1.0/eigenvalues[i];
// 		// 	vert_pos[3*vert_idx]+=eigencoefs[i][0]*eigenfunctions[i][vert_idx]*scale_factor;
// 		// 	vert_pos[3*vert_idx+1]+=eigencoefs[i][1]*eigenfunctions[i][vert_idx]*scale_factor;
// 		// 	vert_pos[3*vert_idx+2]+=eigencoefs[i][2]*eigenfunctions[i][vert_idx]*scale_factor;
// 		// }
// 	}
//
// }
void RealTimeExampleBasedDeformer::reconstructFromEigenCoefs(Vec3d *target_eigencoefs,double *vert_pos,int flag)
{
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// std::cout<<reconstruct_eigenfunction_num_<<"\n";
	double scale_factor;
	memset(vert_pos,0.0,sizeof(double)*3*simulation_vertices_num_);
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
	// if(flag==1)
	// 	saveReconstructMesh(vert_pos);
}
void RealTimeExampleBasedDeformer::saveReconstructMesh(double *new_pos)
{
	static int id=0;
	++id;
	std::string id_str;
	std::stringstream stream;
	stream<<id;
	stream>>id_str;
	std::string file_name="a"+id_str+".smesh";
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
		pos[0]=(*simulation_mesh_->getVertex(i))[0]+target_deformation_[3*i]+new_pos[3*i+0];
		pos[1]=(*simulation_mesh_->getVertex(i))[1]+target_deformation_[3*i+1]+new_pos[3*i+1];
		pos[2]=(*simulation_mesh_->getVertex(i))[2]+target_deformation_[3*i+2]+new_pos[3*i+2];
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

	//bool is_deform_away=(last_initial_weight_-target_weights[0]>=1.0e5);
	double lower_bound=0.5,upper_bound=0.999;
	//bool initial_weight_in_range = (target_weights[0]>lower_bound)&&(target_weights[0]<upper_bound);
	enable_eigen_weight_control_=true;

	// std::cout<<"....target_weights:"<<target_weights[0]<<"......"<<target_weights[1]<<"\n";
	if(enable_eigen_weight_control_&&example_num_>1)
	{
		double bias_factor=0.5; //different biases are used for different simulationMesh
		//example 1 is by default the rest pose
		double other_weight_delta=(1-bias_factor)*target_weights[0]/(example_num_-1);
		target_weights[0]*=bias_factor;
		for(int i=1; i<example_num_;++i)
			target_weights[i]+=other_weight_delta;
		//last_initial_weight_=target_weights[0];
	}
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
	free(temp_buffer);
	pure_example_linear_interpolation_=false;

	// counter1.StopCounter();
	// std::cout<<"elapseTime:"<<counter1.GetElapsedTime()<<"\n";
	// target_weights[0]=0.0;
	// target_weights[1]=1.0;
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	for(int j=0;j<3;++j)
	// 		projected_eigencoefs[i][j]=example_eigencoefs_[1][i][j];
	// PerformanceCounter counter2;
	// counter2.StartCounter();
	//alternative step2
	if(pure_example_linear_interpolation_)
	{
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
		//epsg=integrator_epsilon_*0.0001;
		new_eigencoefs.setcontent(3*interpolate_eigenfunction_num_,temp_buffer);
		mincgcreate(new_eigencoefs,new_state);
		mincgsetcond(new_state,epsg,epsf,epsx,max_iterations);
		mincgoptimize(new_state,evaluateObjectiveAndGradient2,NULL,(void*)input_data);
		mincgresults(new_state,new_eigencoefs,new_report);
		//print some information about the solver if needed
		cout<<"Step2 converged in "<<new_report.iterationscount<<" iterations. Termination type of the solver: "<<new_report.terminationtype<<".\n";

		delete[] temp_buffer;
		delete[] input_data;
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			for(int j=0;j<3;++j)
				projected_eigencoefs[i][j]=new_eigencoefs[3*i+j];
	}

	cout<<"The target_weights in example manifold is:";
	for(int i=0;i<example_num_;++i)
		cout<<target_weights[i]<<",";
	cout<<"\n";

	// counter2.StopCounter();
	// std::cout<<"step2:elapseTime:"<<counter2.GetElapsedTime()<<"\n";
		// counter2.StopCounter();
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<projected_eigencoefs[i]<<"\n";

}

//helper function, used in projectOnExampleManifold()
//evaluate the value of objective function and its gradient, ALGLIB needs it
//'ptr' is passed as input_eigencoefs in projectOnExampleManifold()
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient1(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
{
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
	}
	func*=0.5;
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
}
//step2, minimization of wieghted deformation energy to the examples
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient2(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr)
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
	for(int i=0;i<active_instance->example_num_;++i)
	{
	// 	//compute displacement in LB space for each example
		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
			{
				temp_eigencoefs[j][k]=x[3*j+k]-active_instance->example_eigencoefs_[i][j][k];
			}
		}
		//active_instance->computeF(temp_eigencoefs);
		double energy=0.0;
		active_instance->computeReducedEnergy(temp_eigencoefs,energy);
		func+=target_weights[i]*energy;
		double *energy_grad=new double[3*active_instance->interpolate_eigenfunction_num_];
		memset(energy_grad,0.0,sizeof(double)*3*active_instance->interpolate_eigenfunction_num_);
		active_instance->computeReducedInternalForce(temp_eigencoefs,energy_grad);
		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
				grad[3*j+k]+=target_weights[i]*energy_grad[3*j+k];
		}
		delete[] energy_grad;
	}
	delete[] temp_eigencoefs;
	// counter.StopCounter();
	// std::cout<<"evaluateObject-step2:"<<counter.GetElapsedTime()<<"\n";
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
	Mat3d P(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	P=mu_*(F-trans(inv(F)))+lamda_*log(det(F))*trans(inv(F));
	return P;
}
Mat3d RealTimeExampleBasedDeformer::computeP_gradient(const int &ele,const Mat3d &F,const int &vert_idx,
														const int &vert_idx_dim) const
{
	Mat3d pFpx=computeF_gradient(ele,vert_idx,vert_idx_dim);
	Mat3d temp=inv(F)*pFpx;
	double trace=temp[0][0]+temp[1][1]+temp[2][2];
	Mat3d tiF=trans(inv(F));
	Mat3d pPpx=mu_*pFpx+(mu_-lamda_*log(det(F)))*tiF*trans(pFpx)*tiF+lamda_*trace*tiF;
	return pPpx;
}
void RealTimeExampleBasedDeformer::computeReducedEnergy(const Vec3d *reduced_dis,double &energy) const
{
	energy=0.0;
	// PerformanceCounter counter1;
	// counter1.StartCounter();
	for(int cubica_idx=0;cubica_idx<object_LB_cubica_ele_num_;++cubica_idx)
	{
		Mat3d F=computeF(cubica_idx,reduced_dis);
		Mat3d temp=trans(F)*F;
		double trace_c=temp[0][0]+temp[1][1]+temp[2][2];
		double lnJ=log(det(F));
		double element_energy=0.5*mu_*(trace_c-3)-mu_*lnJ+0.5*lamda_*lnJ*lnJ;
		energy += object_LB_cubica_weights_[cubica_idx]*element_energy;
	}
	// counter1.StopCounter();
	// std::cout<<"energy:"<<counter1.GetElapsedTime()<<"\n";
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
// void RealTimeExampleBasedDeformer::computeReducedStiffnessMatrix(const Vec3d *reduced_dis,Matrix<double> &reduced_K) const
// {
// 	// std::cout<<".................begin............\n";
// 	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
// 	// 	for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 	// 		std::cout<<reduced_K(i,j)<<"..";
// 	// 		std::cout<<"....end\n";
// 	Matrix<double> K((int)interpolate_eigenfunction_num_,(int)interpolate_eigenfunction_num_);
// 	for(int cubica_idx=0;cubica_idx<object_LB_cubica_ele_num_;++cubica_idx)
// 	{
// 		int ele=object_LB_cubica_elements_[cubica_idx];
// 		Matrix<double> subU=tetSubBasis(ele);//12*r
// 		Matrix<double> ele_K(12,12);
// 		//Mat3d F=F_[cubica_idx];
// 		Mat3d F=computeF(cubica_idx,reduced_dis);
// 		Mat3d DmInv=computeDmInv(ele);
// 		for(int i=0;i<4;++i)
// 		{
// 			vector<Mat3d> g_derivative(3);//computes dg/dx_j^0,dg/dx_j^1,dg/dx_j^2
// 			g_derivative.clear();
// 			for(int j=0;j<3;++j)
// 			{
//                 g_derivative[j].set(0.0);
// 				g_derivative[j]=computeP_gradient(ele,F,i,j)*trans(DmInv);
// 			}
//
// 			//computes df0/dx_j,df1/dx_j,df2/dx_j,df3/dx_j
// 			for(int j=0;j<4;++j)
// 			{
// 				Matrix<double> f_derivative(3,3);
// 				for(int row=0;row<3;++row)
// 				{
// 					for(int col=0;col<3;++col)
// 					{
// 						if(j==3)
// 							f_derivative(row,col)=(-1.0)*(g_derivative[row][col][0]+g_derivative[row][col][1]+g_derivative[row][col][2]);
// 						else
// 							f_derivative(row,col)=g_derivative[row][col][j];
// 						ele_K(3*i+row,3*j+col)=f_derivative(row,col);
// 					}
// 				}
// 			}
// 		}
// 		// std::cout<<"```````````````K```````````````:\n";
// 		K+=object_LB_cubica_weights_[cubica_idx]*Transpose(subU)*ele_K*subU;
// 	}
// 	reduced_K=K;
// }
void RealTimeExampleBasedDeformer::testEnergyGradients()
{
	//test energy and energy gradient-internal force
	real_1d_array x;
	int num=3*interpolate_eigenfunction_num_;
	double *temp_buffer=new double[num];
	srand((unsigned)time(0));
	int lowest=1,highest=10;
	int range=(highest-lowest)+1;
	//randomly generate initial x in range[0.1,1]
	for(int i=0;i<num;++i)
		temp_buffer[i]=(lowest+rand()%range)/100.0;
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	for(int j=0;j<3;++j)
	// 		temp_buffer[3*i+j]=object_eigencoefs_[i][j];
	x.setcontent(num,temp_buffer);
	double f_min,f_plus,f;
	real_1d_array temp_grad,grad;
	temp_grad.setcontent(num,temp_buffer);
	grad.setcontent(num,temp_buffer);
	double target_weights[2]={0.5,0.5};
	evaluateObjectiveAndGradient2(x,f,grad,(void*)target_weights);
	for(int i=0;i<num;++i)
	{
		temp_buffer[i]=(lowest+rand()%range)/1.0e8;
		x[i]+=temp_buffer[i];
	}
	evaluateObjectiveAndGradient2(x,f_plus,temp_grad,(void*)target_weights);
	for(int i=0;i<num;++i)
	{
		x[i]-=2*temp_buffer[i];
	}
	evaluateObjectiveAndGradient2(x,f_min,temp_grad,(void*)target_weights);
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
