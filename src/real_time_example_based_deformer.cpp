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
#include "volumetricMesh.h"
#include "volumetricMeshLoader.h"
#include "volumetricMeshENuMaterial.h"
#include "tetMesh.h"
#include "sceneObjectDeformable.h"
#include "coupled_quasi_harmonics.h"
#include "planes.h"
#include "real_time_example_based_deformer.h"
#include "optimization.h"
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
	delete[] reduced_displacement_;
	delete[] reduced_velocity_;
	for(unsigned int i=0;i<reduced_basis_num_;++i)
		if(reduced_basis_[i])
			delete[] reduced_basis_[i];
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
	for(unsigned int i=0;i<example_num_;++i)
	{
		delete[] example_cubica_weights_[i];
		delete[] example_cubica_elements_[i];
	}
	delete[] example_cubica_elements_[example_num_];
	delete[] example_cubica_weights_;
	delete[] example_cubica_elements_;
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
}

RealTimeExampleBasedDeformer* RealTimeExampleBasedDeformer::activeInstance()
{
    return active_instance_;
}
void RealTimeExampleBasedDeformer::setupSimulation()
{
    //TO DO
}

void RealTimeExampleBasedDeformer::advanceStep()
{
    //TO DO
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
	//init element displacement matrix Dm
	int ele_num=simulation_mesh_->getNumElements();
	object_init_element_dis_matrix_ = new Mat3d[ele_num];
	for(unsigned int ele=0;ele<ele_num;++ele)
	{
		int *global_idx=new int[ele_num];
		Vec3d *vert_pos = new Vec3d[ele_num];
		for(unsigned int i=0;i<simulation_mesh_->getNumElementVertices();++i)
		{
			global_idx[i] = simulation_mesh_->getVertexIndex(ele,i);
			vert_pos[i] = *simulation_mesh_->getVertex(global_idx[i]);
		}
		for(unsigned int dim=0; dim<3; ++dim)
		{
			object_init_element_dis_matrix_[ele][dim][0]=vert_pos[0][dim]-vert_pos[3][dim];
			object_init_element_dis_matrix_[ele][dim][1]=vert_pos[1][dim]-vert_pos[3][dim];
			object_init_element_dis_matrix_[ele][dim][2]=vert_pos[2][dim]-vert_pos[3][dim];
		}
		delete[] global_idx;
		delete[] vert_pos;
	}
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
	//init example element displacement matrix Dm
	example_init_element_dis_matrix_ = new Mat3d*[example_num_];
	for(unsigned int ex_idx=0;ex_idx<example_num_;++ex_idx)
	{
		int ele_num=examples_[ex_idx]->getNumElements();
		example_init_element_dis_matrix_[ex_idx] = new Mat3d[ele_num];
		for(unsigned int ele=0;ele<ele_num;++ele)
		{
			int *vert_idx=new int[ele_num];
			Vec3d *vert_pos = new Vec3d[ele_num];
			for(unsigned int i=0;i<examples_[ex_idx]->getNumElementVertices();++i)
			{
				vert_idx[i] = examples_[ex_idx]->getVertexIndex(ele,i);
				vert_pos[i] = *examples_[ex_idx]->getVertex(vert_idx[i]);
			}
			for(unsigned int dim=0; dim<3; ++dim)
			{
				example_init_element_dis_matrix_[ex_idx][ele][dim][0]=vert_pos[0][dim]-vert_pos[3][dim];
				example_init_element_dis_matrix_[ex_idx][ele][dim][1]=vert_pos[1][dim]-vert_pos[3][dim];
				example_init_element_dis_matrix_[ex_idx][ele][dim][2]=vert_pos[2][dim]-vert_pos[3][dim];
			}
			delete[] vert_idx;
			delete[] vert_pos;
		}
	}
	ex_dis_=new double*[example_num_];
	for(int i=0;i<example_num_;++i)
		ex_dis_[i]=new double[3*examples_[i]->getNumVertices()];
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
	//read file, save as reduced_basis_[3*vert_idx+dim][reduced_idx]
	unsigned int reduced_row_num,reduced_col_num;
	std::getline(input_file,temp_str);
	reduced_row_num=atoi(temp_str.c_str());
	std::getline(input_file,temp_str);
	reduced_col_num=atoi(temp_str.c_str());
	reduced_basis_num_=reduced_col_num;
	reduced_basis_=new double *[reduced_row_num];
	std::cout<<reduced_row_num<<std::endl;
	std::cout<<reduced_col_num<<std::endl;
	if(reduced_row_num*3!=simulation_mesh_->getNumVertices())
	{
		std::cout<<"The input reduced basis is not correct!\n";
		return false;
	}
	for(unsigned int i=0;i<reduced_row_num;++i)
	{
		reduced_basis_[i]=new double[reduced_col_num];
	}
	unsigned int line_num=0,str_num=0,total_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		++total_num;
		double temp_value;
		input_file>>temp_value;
		reduced_basis_[line_num][str_num]=temp_value;
		if(str_num>=reduced_col_num-1)
		{
			str_num=0;
			if(line_num<=reduced_row_num-1)
				++line_num;
		}
		else
			++str_num;
		if(total_num>=reduced_row_num*reduced_col_num)
			break;
	}
	input_file.close();
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
	// std::cout<<"-------------\n";
	// for(int i=0;i<reconstruct_eigenfunction_num_;++i)
	// 	std::cout<<object_eigencoefs_[i]<<"\n";
	delete[] dis;
	//temp store the last 3 Eigenfunctions
	// std::ofstream wfile("gorilla_l3.eigen");
	// if(!wfile)
	// {
	// 	std::cout<<"error:\n";
	// }
	// wfile<<"*eigenValues"<<std::endl;
	// wfile<<"3"<<std::endl;
	// for(int i=4;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	if(i==interpolate_eigenfunction_num_-1)
	// 		wfile<<object_eigenvalues_[i];
	// 	else
	// 		wfile<<object_eigenvalues_[i]<<" ";
	// }
	// wfile<<std::endl;
	// wfile<<"*eigenVectors:"<<std::endl;
	// wfile<<simulation_mesh_->getNumVertices()<<std::endl;
	// wfile<<"3"<<std::endl;
	// for(int j=0;j<simulation_mesh_->getNumVertices();++j)
	// {
	// 	for(int i=4;i<7;++i)
	// 	{
	// 		wfile<<object_eigenfunctions_[i][j]<<" ";
	// 	}
	// 	wfile<<std::endl;
	// }
	// wfile.close();
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
	// for(int i=0;i<example_num_;++i)
	// {
	// 	std::cout<<"example"<<i<<"\n";
	// 	for(int j=0;j<reconstruct_eigenfunction_num_;++j)
	// 		std::cout<<example_eigencoefs_[i][j]<<"\n";
	// }

	// std::ofstream wfile("new1_ex_l3_0.eigen");
	// if(!wfile)
	// {
	// 	std::cout<<"error:\n";
	// }
	// wfile<<"*eigenValues"<<std::endl;
	// wfile<<"3"<<std::endl;
	// for(int i=4;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	if(i==interpolate_eigenfunction_num_-1)
	// 		wfile<<example_eigenvalues_[0][i];
	// 	else
	// 		wfile<<example_eigenvalues_[0][i]<<" ";
	// }
	// wfile<<std::endl;
	// wfile<<"*eigenVectors:"<<std::endl;
	// wfile<<examples_[0]->getNumVertices()<<std::endl;
	// wfile<<"3"<<std::endl;
	// for(int j=0;j<examples_[0]->getNumVertices();++j)
	// {
	// 	for(int i=4;i<7;++i)
	// 	{
	// 		wfile<<example_eigenfunctions_[0][i][j]<<" ";
	// 		std::cout<<example_eigenfunctions_[0][i][j]<<" ";
	// 	}
	// 	wfile<<std::endl;
	// }
	// wfile.close();
    return true;
}

bool RealTimeExampleBasedDeformer::loadPlanesInScene(const std::string &file_name, unsigned int plane_num)
{
    if(file_name=="")
	{
		std::cout<<"Error: cannot read "<<file_name<<".\n";
		return false;
	}
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
	std::fstream input_file(file_name);
	if(!input_file)
	{
		std::cout<<"Error: failed to open "<<file_name<<" .\n";
		return false;
	}
	std::cout<<"Load fixed vertices file:\n";
	string temp_str;
	std::vector<unsigned int> fixed_vertices_vector;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		input_file>>temp_str;
		fixed_vertices_vector.push_back(atoi(temp_str.c_str()));
	}
	fixed_vertex_num_=fixed_vertices_vector.size();
	fixed_vertices_=new unsigned int[fixed_vertices_vector.size()];
	for(unsigned int i=0;i<fixed_vertices_vector.size();++i)
	{
		fixed_vertices_[i]=fixed_vertices_vector[i];
	}
	input_file.close();
	std::cout<<"Load fixed vertices file done.\n";
	return true;
}
//tetID : 0-indexed
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
		object_cubica_elements_[str_num]=temp_value;
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
	std::cout<<"Load object cubica data succeed.\n";
	//isload_object_cubica_=true;
	return true;
}
bool RealTimeExampleBasedDeformer::loadExampleCubicaData(const std::string &file_name_prefix)
{
	std::cout<<"Load cubica data for examples...\n";
	if(example_num_==0)
	{
		std::cout<<"examle cubica file unloaded.\n";
		return false;
	}
	example_cubica_ele_num_ = new unsigned int[example_num_];
	example_cubica_elements_ = new unsigned int*[example_num_];
	example_cubica_weights_ = new double*[example_num_];
	for(unsigned int i=0;i<example_num_;++i)
	{
        std::string file_num_str,file_name;
        std::stringstream adaptor;
		adaptor.str("");
		adaptor.clear();
		adaptor<<i;
		adaptor>>file_num_str;
		//read file, file format is .eigen,
		file_name=file_name_prefix+file_num_str+std::string(".cubature");
        std::fstream input_file(file_name.c_str());
		if(!input_file)
		{
			std::cout<<"Error: Cannot open file "<<file_name<<std::endl;
			return false;
		}
		string temp_str;
		std::getline(input_file,temp_str);
		example_cubica_ele_num_[i]=atoi(temp_str.c_str());
		example_cubica_elements_[i] = new unsigned int[example_cubica_ele_num_[i]];
		example_cubica_weights_[i] = new double [example_cubica_ele_num_[i]];
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
			example_cubica_elements_[i][str_num]=temp_value;
			str_num++;
			if(str_num>=example_cubica_ele_num_[i])
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
			example_cubica_weights_[i][str_num]=temp_value;
			++str_num;
			if(str_num>=example_cubica_ele_num_[i])
				break;
		}
	}
	std::cout<<"Load example cubica data succeed.\n";
	//isload_example_cubica_=true;
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
	std::cout<<"0\n";
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
		std::cout<<"1\n";
		if(points_in_region>1)//corresponding functions are region based
			is_region_based_correspondence_=true;
		for(unsigned int ex_idx=0;ex_idx<example_num_;++ex_idx)
		{
			//points on the examples
			std::getline(input_file,line_str);
			//std::cout<<line_str<<"\n";
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
	std::cout<<"2\n";
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
	//std::cout<<"col_b:"<<col_b<<",row_b:"<<row_b<<"\n";
	//change to col order, each col is an eigenfunction
	for(int i=0;i<col_b;++i)
		for(int j=0;j<row_b;++j)
			obj_eigenfunction_col[j][i]=object_eigenfunctions_[i][j];
	//std::cout<<"example_num:"<<example_num_<<"\n";
//	for(int j=0;j<7;++j)
	//std::cout<<example_eigenfunctions_[0][j][2615]<<",";
	std::cout<<"\n";
	for(int i=0;i<example_num_;++i)
	{
	//	std::cout<<"example_eigenfunctions_(0,6,2615):"<<example_eigenfunctions_[0][6][2615]<<"\n";
		//std::cout<<"~~~~~~~~\n";
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
		std::cout<<"col_a:"<<col_a<<",row_a:"<<row_a<<"\n";
		for(int j=0;j<col_a;++j)
			for(int k=0;k<row_a;++k)
				ex_eigenfunction_col[k][j]=example_eigenfunctions_[i][j][k];

		std::cout<<"\n";
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
		//std::cout<<"example_num_:"<<example_num_<<"row:"<<row_a<<"corresponding_num:"<<corresponding_function_num_<<"\n";
	//	std::cout<<"row_a:"<<row_a<<",corresponding_num:"<<corresponding_function_num_<<"\n";
		//std::cout<<"example_eigenfunctions_(0,6,2615):"<<example_eigenfunctions_[0][6][2615]<<"\n";
		//std::cout<<"new:"<<new_ex_eigenfunction_col[2615][6];
		for(int j=0;j<row_a;++j)
			for(int k=0;k<interpolate_eigenfunction_num_;++k)
				example_eigenfunctions_[i][k][j]=new_ex_eigenfunction_col[j][k];
		std::cout<<"----example_eigenfunctions changed.\n";

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
std::cout<<"e\n";
	delete[] obj_eigenfunction_col;
	delete[] new_obj_eigenfunction_col;
    return true;
}

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
			//std::cout<<mesh->getVertex(vert_idx)<<std::endl;
			//std::cout<<mesh->getVertex(vert_idx)[vert_idx]<<std::endl;
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
void RealTimeExampleBasedDeformer::reconstructFromEigenCoefs(double **eigenfunctions,double *eigenvalues,Vec3d *eigencoefs,Vec3d *target_eigencoefs,
                                                            const int &eigenfunction_num,const int &input_reconstruct_eigenfunction_num,const int &vert_num, double *vert_pos)
{
	//std::cout<<"before reconstruct:\n";
	// for(int i=0;i<eigenfunction_num;++i)
	// 	std::cout<<target_eigencoefs[i]<<",";
	// std::cout<<std::endl;
    double scale_factor;
	for(unsigned int vert_idx=0;vert_idx<vert_num;++vert_idx)
	{
		vert_pos[3*vert_idx]=vert_pos[3*vert_idx+1]=vert_pos[3*vert_idx+2]=0.0;
		for(unsigned int i=0;i<eigenfunction_num;++i)
		{
			scale_factor=1.0/eigenvalues[i];
			vert_pos[3*vert_idx]+=target_eigencoefs[i][0]*eigenfunctions[i][vert_idx]*scale_factor;
			vert_pos[3*vert_idx+1]+=target_eigencoefs[i][1]*eigenfunctions[i][vert_idx]*scale_factor;
			vert_pos[3*vert_idx+2]+=target_eigencoefs[i][2]*eigenfunctions[i][vert_idx]*scale_factor;

		}
		//reconstruct the higher spectral pose eigencoefs from the simulation object itself
		for(unsigned int i=eigenfunction_num;i<input_reconstruct_eigenfunction_num;++i)
		{
			//std::cout<<i<<"\n";
			scale_factor=1.0/eigenvalues[i];
			vert_pos[3*vert_idx]+=eigencoefs[i][0]*eigenfunctions[i][vert_idx]*scale_factor;
			vert_pos[3*vert_idx+1]+=eigencoefs[i][1]*eigenfunctions[i][vert_idx]*scale_factor;
			vert_pos[3*vert_idx+2]+=eigencoefs[i][2]*eigenfunctions[i][vert_idx]*scale_factor;
		}
		// if(vert_idx==0)
		// std::cout<<"vert_pos0:"<<vert_pos[3*vert_idx]<<","<<vert_pos[3*vert_idx+1]<<","<<vert_pos[3*vert_idx+2]<<"\n";
	}

}

void RealTimeExampleBasedDeformer::projectOnExampleManifold(Vec3d *input_eigencoefs, Vec3d *projected_eigencoefs)
{
	//suppose the target configuration is: w(0)*example_eigencoefs_(0)+w(1)*example_eigencoefs_(1)+...+w(example_num_)*example_eigencoefs_(exmaple_num_-1)
	//w(0)+w(1)+...+w(example_num_-1)=1, 0<=w(i)<=1, i=0,1,...,example_num_-1
	//compute the weights by minimizing 1/2*dist(target,input)*dist(target_input), least square minimization with constraints
	//solve this minimization using optimization solver from ALGLIB

	// std::cout<<"projectOnExampleManifold:-a\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<input_eigencoefs[i]<<"\n";

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
	upper_boundary.setcontent(example_num_,temp_buffer);
	//prepare the equality constraints: w(0)+w(1)+...+w(example_num_-1)=1
	real_2d_array equality_constraint;
	equality_constraint.setcontent(1,example_num_+1,temp_buffer);
	integer_1d_array constraint_type="[0]";
	//stop conditions
	double epsg=integrator_epsilon_*0.1;
	double epsf=integrator_epsilon_*0.1;
	double epsx=integrator_epsilon_*0.1;
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
	//
	std::cout<<"eeeeeeeeeeeeeeeeeeeeenalbe:"<<enable_eigen_weight_control_<<"\n";
	enable_eigen_weight_control_=true;
	if(enable_eigen_weight_control_&&example_num_>1)
	{
		double bias_factor=0.85; //different biases are used for different simulationMesh
		//example 1 is by default the rest pose
		double other_weight_delta=(1-bias_factor)*target_weights[0]/(example_num_-1);
		target_weights[0]*=bias_factor;
		for(int i=1; i<example_num_;++i)
			target_weights[i]+=other_weight_delta;
		//last_initial_weight_=target_weights[0];
	}
	//pure_example_linear_interpolation=true;
	//compute the target coefficients with the weights
	//  target_weights[0]=0.2;
	//  target_weights[1]=0.8;
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
	//testObjectiveGradients();
	//alternative step2
	if(!pure_example_linear_interpolation_)
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
		max_iterations=100;
		new_eigencoefs.setcontent(3*interpolate_eigenfunction_num_,temp_buffer);
		mincgcreate(new_eigencoefs,new_state);
		mincgsetcond(new_state,epsg,epsf,epsx,max_iterations);
		mincgoptimize(new_state,evaluateObjectiveAndGradient2,NULL,(void*)input_data);
		mincgresults(new_state,new_eigencoefs,new_report);
		//print some information about the solver if needed
		cout<<"Step2 converged in "<<new_report.iterationscount<<" iterations. Terminatin type of the solver: "<<new_report.terminationtype<<".\n";

		delete[] temp_buffer;
		delete[] input_data;
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
			for(int j=0;j<3;++j)
				projected_eigencoefs[i][j]=new_eigencoefs[3*i+j];
	}

	cout<<"The target_weights in example manifold is:\n";
	for(int i=0;i<example_num_;++i)
		cout<<target_weights[i]<<",";
	cout<<"\n";
	std::cout<<"target_eigencoefs:\n";
	for(int i=0;i<interpolate_eigenfunction_num_;++i)
		std::cout<<projected_eigencoefs[i]<<"\n";
	std::cout<<"----------------------------------------------------weight compute end.\n";
//	getchar();
}

//helper function, used in projectOnExampleManifold()
//evaluate the value of objective function and its gradient, ALGLIB needs it
//'ptr' is passed as input_eigencoefs in projectOnExampleManifold()
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient1(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
{
	//std::cout<<"evaluateObjectiveAndGradient1---start:\n";

	RealTimeExampleBasedDeformer* active_instance=RealTimeExampleBasedDeformer::activeInstance();
	assert(active_instance);
	Vec3d *input_eigencoefs=(Vec3d*)ptr;
	Vec3d *temp_target_eigencoefs=new Vec3d[active_instance->interpolate_eigenfunction_num_];
	//evaluate the objective function
	func=0.0;
	// std::cout<<"input_eigencoefs:\n";
	// for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
	// 	std::cout<<input_eigencoefs[j]<<"\n";
	// for(int i=0;i<active_instance->example_num_;++i)
	// 	std::cout<<"x:"<<x[i]<<",";
	// std::cout<<"example_eigencoefs:\n";
	// for(int j=0;j<active_instance->example_num_;++j)
	// 	for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
	// 		std::cout<<active_instance->example_eigencoefs_[j][i]<<",";
	// std::cout<<"\n";
	// std::cout<<"temp_target_eigencoefs:\n";
	//
	// std::cout<<"~~~~~~~~~~~~~~:"<<x[0]<<","<<x[1]<<"\n";
	for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
	{
		temp_target_eigencoefs[i][0]=temp_target_eigencoefs[i][1]=temp_target_eigencoefs[i][2]=0.0;
		for(int j=0;j<active_instance->example_num_;++j)
		{
			temp_target_eigencoefs[i][0]+=x[j]*active_instance->example_eigencoefs_[j][i][0];
			temp_target_eigencoefs[i][1]+=x[j]*active_instance->example_eigencoefs_[j][i][1];
			temp_target_eigencoefs[i][2]+=x[j]*active_instance->example_eigencoefs_[j][i][2];

		//	std::cout<<"example"<<j<<":"<<active_instance->example_eigencoefs_[j][i]<<"\n";
		}
	//	std::cout<<temp_target_eigencoefs[i]<<"\n";
		func+=(temp_target_eigencoefs[i][0]-input_eigencoefs[i][0])*(temp_target_eigencoefs[i][0]-input_eigencoefs[i][0]);
		func+=(temp_target_eigencoefs[i][1]-input_eigencoefs[i][1])*(temp_target_eigencoefs[i][1]-input_eigencoefs[i][1]);
		func+=(temp_target_eigencoefs[i][2]-input_eigencoefs[i][2])*(temp_target_eigencoefs[i][2]-input_eigencoefs[i][2]);
	//	std::cout<<"temp_target_eigencoefs:"<<temp_target_eigencoefs[i][0]<<","<<temp_target_eigencoefs[i][1]<<","<<temp_target_eigencoefs[i][2]<<"\n";
	//	std::cout<<"input_eigencoefs:"<<input_eigencoefs[i][0]<<","<<input_eigencoefs[i][1]<<","<<input_eigencoefs[i][2]<<"\n";
	}
	func*=0.5;
	// std::cout<<"func:"<<func<<",\n";
	// std::cout<<"grad:\n";
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
		//std::cout<<i<<":"<<grad[i]<<",";
	}
	delete[] temp_target_eigencoefs;
}
void RealTimeExampleBasedDeformer::testObjectiveGradients()
{
	real_1d_array x;

	int num=interpolate_eigenfunction_num_*3;
	double *temp_buffer=new double[num];
	for(int i=0;i<interpolate_eigenfunction_num_;++i)
		for(int j=0;j<3;++j)
			temp_buffer[3*i+j]=example_eigencoefs_[0][i][j];
	x.setcontent(num,temp_buffer);
	for(int i=0;i<num;++i)
		temp_buffer[i]=0.0;
	double target_weights[2]={1.0,0.0};
	double f=0.0;
	real_1d_array grad;
	grad.setcontent(num,temp_buffer);
	evaluateObjectiveAndGradient2(x,f,grad,(void*)target_weights);
	delete[] temp_buffer;
}
//step2, minimization of wieghted deformation energy to the examples
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient2(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr)
{
	std::cout<<"evaluateObjectiveAndGradient2---start:\n";
	RealTimeExampleBasedDeformer* active_instance=RealTimeExampleBasedDeformer::activeInstance();
	assert(active_instance);
	// if(!active_instance->isload_cubica_file)
	// {
	// 	std::cout<<"Error: cubica file unloaded.\n";
	// 	exit(0);
	// }
	for(int i=0;i<3*active_instance->interpolate_eigenfunction_num_;++i)
	{
		std::cout<<"x:"<<x[i]<<",";
	}
	double *target_weights=(double*)ptr;
	double total_energy=0.0;
	Vec3d *temp_eigencoefs=new Vec3d[active_instance->interpolate_eigenfunction_num_];
//	std::cout<<"temp_eigencoefs:\n";
	// for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
	// {
	// 	for(int j=0;j<3;++j)
	// 	{
	// 		temp_eigencoefs[i][j]=x[3*i+j];
	// 		//temp_eigencoefs[i][j]=active_instance->example_eigencoefs_[0][i][j];
	// 	//	std::cout<<temp_eigencoefs[i][j]<<",";
	//
	// 	}
	//  //	std::cout<<"\n";
	// }

	for(int j=0;j<3*active_instance->interpolate_eigenfunction_num_;++j)
		grad[j]=0.0;
	for(int i=0;i<active_instance->example_num_;++i)
	{
	//	std::cout<<"example"<<i<<std::endl;
		//compute displacement in LB space for each example
		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
			{
				temp_eigencoefs[j][k]=x[3*j+k]-active_instance->example_eigencoefs_[i][j][k];
			//	std::cout<<temp_eigencoefs[i][j]<<",";
			}
		//	std::cout<<"\n";
		}
		double energy=0.0;
		double *energy_grad=new double[3*active_instance->interpolate_eigenfunction_num_];
		//convert x to deformed_pos of example[i],reconstruct shape from LB space to Euclidean space
		//double *deformed_pos=new double[3*active_instance->simulation_mesh_->getNumVertices()];
		double *displacement=new double[3*active_instance->simulation_mesh_->getNumVertices()];
		//memset(active_instance->ex_dis_[i],0.0,sizeof(double)*(3*active_instance->simulation_mesh_->getNumVertices()));
		memset(displacement,0.0,sizeof(double)*(3*active_instance->simulation_mesh_->getNumVertices()));
		//memset(deformed_pos,0.0,sizeof(double)*(3*active_instance->simulation_mesh_->getNumVertices()));
		std::cout<<"example"<<i<<":\n";
		//reconstruct full space configuration for x each iteration step
		// active_instance->reconstructFromEigenCoefs(active_instance->exampleEigenFunctions()[i],active_instance->exampleEigenValues()[i],
		// 							active_instance->exampleEigencoefs()[i],temp_eigencoefs,
		// 							active_instance->interpolate_eigenfunction_num_,active_instance->reconstruct_eigenfunction_num_,
		// 							active_instance->examples_[i]->getNumVertices(),deformed_pos);
		active_instance->reconstructFromEigenCoefs(active_instance->objectEigenFunctions(),active_instance->objectEigenValues(),
									active_instance->objectEigencoefs(),temp_eigencoefs,
									active_instance->interpolate_eigenfunction_num_,active_instance->reconstruct_eigenfunction_num_,
									active_instance->simulation_mesh_->getNumVertices(),displacement);
		//double scale_factor;
		// for(int vert_idx=0;vert_idx<active_instance->simulation_mesh_->getNumVertices();++vert_idx)
		// {
		// 	Vec3d init_pos=*active_instance->simulation_mesh_->getVertex(vert_idx);
		// 	for(int dim=0;dim<3;++dim)
		// 	{
		// 		displacement[3*vert_idx+dim]=deformed_pos[3*vert_idx+dim]/*-init_pos[dim]*/;
		// 		//active_instance->ex_dis_[i][3*vert_idx+dim]=displacement[3*vert_idx+dim];
		// 	}
		//
		// }
		//getchar();
		// active_instance->computeReducedEnergyAndGradient(active_instance->examples_[i],NULL,displacement,active_instance->example_cubica_ele_num_[i],
	 // 		active_instance->example_cubica_elements_[i],active_instance->example_cubica_weights_[i],1,i,energy,energy_grad);
		active_instance->computeReducedEnergyAndGradient(active_instance->simulation_mesh_,NULL,displacement,active_instance->object_cubica_ele_num_,
 				active_instance->object_cubica_elements_,active_instance->object_cubica_weights_,0,0,energy,energy_grad);
		func+=target_weights[i]*energy;


		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
				grad[3*j+k]+=target_weights[i]*energy_grad[3*j+k];
		}


		delete[] displacement;
		//delete[] deformed_pos;
		delete[] energy_grad;
		//delete[] active_instance_->ex_dis_;
	}
	delete[] temp_eigencoefs;



	//compute gradient force
}


//compute for E(deformed,source)
//1: mesh is the source mesh, first get the mesh material of source
//2: souce pos could be Null, then get the source pos from the given mesh
//3: get the cubica_elements from source mesh
//4: get the vertices position (deformed pos and source pos) for the corresponding elements
void RealTimeExampleBasedDeformer::computeReducedEnergyAndGradient(VolumetricMesh *mesh,const double *init_pos,const double *displacement, const unsigned int cubica_num,
									const unsigned int *cubica_elements, double *cubica_weights,const unsigned int example_flag,
									const unsigned int dis_ex_idx,double &energy,double *grad)
{
	//std::cout<<"computeReducedEnergyAndGradient---start:\n";
	// for(int i=0;i<3*mesh->getNumVertices();++i)
	// 	std::cout<<displacement[i]<<",";
	// std::cout<<std::endl;
	//clear grad vector
	energy=0.0;
	for(int i=0;i<3*interpolate_eigenfunction_num_;++i)
		grad[i]=0.0;
	//material setting, we use homogeneous materials,get the first element material
	VolumetricMesh::Material * material = mesh->getElementMaterial(0);
	VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
	if (eNuMaterial == NULL)
	{
		std::cout<<"Error: mesh does not consist of E, nu materials.\n";
		exit(0);
	}
	double lambda = eNuMaterial->getLambda();
	double mu = eNuMaterial->getMu();
	TetMesh *tet_mesh=dynamic_cast<TetMesh*>(mesh);

	//used for the method1:compute g on full space and redesign a 3n*1 vector
	//the four vertices contains displacement, the other points displacement is zero
//	double *new_dis=new double[3*mesh->getNumVertices()];
	double *g=new double[3*interpolate_eigenfunction_num_];
	for(int i=0;i<3*interpolate_eigenfunction_num_;++i)
		g[i]=0.0;
	//double total_weight=0.0;
	// for(int cubica_idx=0;cubica_idx<cubica_num;++cubica_idx)
	// 	total_weight+=cubica_weights[cubica_idx];
	// std::cout<<"total_weight:"<<total_weight<<std::endl;
	//generate new dis and energy value
	//for(int cubica_idx=0;cubica_idx<cubica_num;++cubica_idx)
	for(int ele=0;ele<mesh->getNumElements();++ele)
	{
		//int ele=cubica_elements[cubica_idx];
		//std::cout<<"eleeeeeeeeeeeeeeeeeeee:"<<ele<<"\n";
		double ele_volume=tet_mesh->getTetVolume(tet_mesh->getVertex(ele,0),tet_mesh->getVertex(ele,1),
								tet_mesh->getVertex(ele,2),tet_mesh->getVertex(ele,3));
	//	std::cout<<ele_volume<<"a\n";
	//	memset(new_dis,0.0,sizeof(double)*(3*mesh->getNumVertices()));
		//get global_idx and vert_pos for each cubica element
		int *global_idx=new int[mesh->getNumElementVertices()];
		Vec3d *vert_pos=new Vec3d[mesh->getNumElementVertices()];
		for(int j=0;j<mesh->getNumElementVertices();++j)
		{
			global_idx[j]=mesh->getVertexIndex(ele,j);
			vert_pos[j]=*mesh->getVertex(global_idx[j]);
		}
		Vec3d *init_ele_dis=new Vec3d[mesh->getNumElementVertices()];
		Vec3d *deformed_ele_dis=new Vec3d[mesh->getNumElementVertices()];
		for(int j=0;j<mesh->getNumElementVertices();++j)
		{
			for(int k=0;k<3;++k)
			{
				if(init_pos==NULL)
					init_ele_dis[j][k]=vert_pos[j][k];
				else
					init_ele_dis[j][k]=init_pos[3*global_idx[j]+k];
				deformed_ele_dis[j][k]=init_ele_dis[j][k]+displacement[3*global_idx[j]+k];
				//std::cout<<"dis:"<<displacement[3*global_idx[j]+k]<<",";
			//	new_dis[3*global_idx[j]+k]=displacement[3*global_idx[j]+k];//generate n Vec3d for all vertices
			}
		}
		//std::cout<<std::endl;
		// for(int j=0;j<mesh->getNumElementVertices();++j)
		// 	std::cout<<"init_ele_dis:"<<init_ele_dis[j]<<std::endl;
		//get init matrix for one ele
		Mat3d temp_init_matrix(0.0);
		//compute init displacement matrix Dm, if input parameter init_pos is not null, the init_pos is used for the current configuration
		if(init_pos==NULL)
		{
			//if(example_flag==0)
				temp_init_matrix=object_init_element_dis_matrix_[ele];
			// else
			// 	temp_init_matrix=example_init_element_dis_matrix_[dis_ex_idx][ele];
		}
		else
		{
			for(int i=0;i<3;++i)
			{
				temp_init_matrix[i][0]=init_ele_dis[0][i]-init_ele_dis[3][i];
				temp_init_matrix[i][1]=init_ele_dis[1][i]-init_ele_dis[3][i];
				temp_init_matrix[i][2]=init_ele_dis[2][i]-init_ele_dis[3][i];
			}
		}
		//compute deformed_displacement_matrix
		Vec3d *vec=new Vec3d[3];
		for(int i=0;i<3;++i)
		{
			vec[i][0]=deformed_ele_dis[0][i]-deformed_ele_dis[3][i];
			vec[i][1]=deformed_ele_dis[1][i]-deformed_ele_dis[3][i];
			vec[i][2]=deformed_ele_dis[2][i]-deformed_ele_dis[3][i];
		}
		Mat3d deformed_ele_dis_matrix(vec[0],vec[1],vec[2]);
		//compute energy
		Mat3d F=computeDeformationGradient(temp_init_matrix,deformed_ele_dis_matrix);
		//  std::cout<<"temp_init_matrix:"<<temp_init_matrix<<std::endl;
		//  std::cout<<"deformed_ele_dis_matrix:"<<deformed_ele_dis_matrix<<std::endl;
		//  std::cout<<"F:"<<F<<std::endl;
		// std::cout<<"det(temp_init_matrix):"<<det(temp_init_matrix)<<std::endl;
		// std::cout<<"det(deformed_ele_dis_matrix):"<<det(deformed_ele_dis_matrix)<<std::endl;
		//std::cout<<"det(F):"<<det(F)<<std::endl;
		Mat3d temp=trans(F)*F;
		double trace_c=temp[0][0]+temp[1][1]+temp[2][2];
		double lnJ=log(det(F));
		double element_energy=0.5*mu*(trace_c-3)-mu*lnJ+0.5*lambda*lnJ*lnJ;
		//std::cout<<"ele_energy:"<<element_energy<<std::endl;
	//	std::cout<<"ele_energy:"<<element_energy<<std::endl;
	//	std::cout<<"cubica_weights:"<<cubica_weights[cubica_idx]<<std::endl;

	//	energy += cubica_weights[cubica_idx]*element_energy;
		energy +=element_energy;

		delete[] global_idx;
		delete[] vert_pos;
		delete[] init_ele_dis;
		delete[] deformed_ele_dis;
		delete[] vec;
	}
		//compute energy gradient
		computeForceOnReducedSubSpace(mesh,init_pos,displacement,example_flag,dis_ex_idx,g);
		for(int ele=0;ele<mesh->getNumElements();++ele)
		 {
			for(int eigen_idx=0;eigen_idx<interpolate_eigenfunction_num_;++eigen_idx)
			{
				// grad[3*eigen_idx]+=cubica_weights[cubica_idx]*g[3*eigen_idx];
				// grad[3*eigen_idx+1]+=cubica_weights[cubica_idx]*g[3*eigen_idx+1];
				// grad[3*eigen_idx+2]+=cubica_weights[cubica_idx]*g[3*eigen_idx+2];
				grad[3*eigen_idx]=g[3*eigen_idx];
				grad[3*eigen_idx+1]=g[3*eigen_idx+1];
				grad[3*eigen_idx+2]=g[3*eigen_idx+2];
			}
		}
	//}
		std::cout<<"energy:"<<energy<<"c\n";
	//delete[] new_dis;
	delete[] g;
}
void RealTimeExampleBasedDeformer::computeForceOnReducedSubSpace(VolumetricMesh *mesh,const double *init_pos,const double *displacement,const unsigned int example_flag,
									const unsigned int dis_ex_idx,double *g)
{
//	std::cout<<"computeForceOnReducedSubSpace:\n";
	//clear g
	for(int i=0;i<3*interpolate_eigenfunction_num_;++i)
		g[i]=0.0;
	VolumetricMesh::Material * material = mesh->getElementMaterial(0);
	VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
	if (eNuMaterial == NULL)
	{
		std::cout<<"Error: mesh does not consist of E, nu materials.\n";
		exit(0);
	}
	double lambda = eNuMaterial->getLambda();
	double mu = eNuMaterial->getMu();
	TetMesh *tet_mesh=dynamic_cast<TetMesh*>(mesh);
	double *forces=new double[3*mesh->getNumVertices()];
	for(int i=0;i<3*mesh->getNumVertices();++i)
		forces[i]=0.0;
	for(int ele=0;ele<mesh->getNumElements();++ele)
	{
		int *global_idx=new int[mesh->getNumElementVertices()];
		Vec3d *vert_pos=new Vec3d[mesh->getNumElementVertices()];
		double ele_volume=tet_mesh->getTetVolume(tet_mesh->getVertex(ele,0),tet_mesh->getVertex(ele,1),
								tet_mesh->getVertex(ele,2),tet_mesh->getVertex(ele,3));
		for(int i=0;i<mesh->getNumElementVertices();++i)
		{
			global_idx[i]=mesh->getVertexIndex(ele,i);
			vert_pos[i]=*mesh->getVertex(global_idx[i]);
		}
		Vec3d *init_ele_dis=new Vec3d[mesh->getNumElementVertices()];
		Vec3d *deformed_ele_dis=new Vec3d[mesh->getNumElementVertices()];
		//get deformed_ele_displacement vector to compute for F
		Mat3d temp_init_matrix(0.0)	;
		//compute init displacement matrix Dm, if input parameter init_pos is not null, the init_pos is used for the current configuration
		if(init_pos==NULL)
		{
			if(example_flag==0)
				temp_init_matrix=object_init_element_dis_matrix_[ele];
			else
				temp_init_matrix=example_init_element_dis_matrix_[dis_ex_idx][ele];
			for(int i=0;i<mesh->getNumElementVertices();++i)
			{
				for(int j=0;j<3;++j)
					deformed_ele_dis[i][j]=vert_pos[i][j]+displacement[3*global_idx[i]+j];
			}
		}
		else
		{
			for(int i=0;i<3;++i)
			{
				temp_init_matrix[i][0]=init_pos[3*global_idx[0]+i]-init_pos[3*global_idx[3]+i];
				temp_init_matrix[i][1]=init_pos[3*global_idx[1]+i]-init_pos[3*global_idx[3]+i];
				temp_init_matrix[i][2]=init_pos[3*global_idx[2]+i]-init_pos[3*global_idx[3]+i];
			}
			for(int i=0;i<mesh->getNumElementVertices();++i)
			{
				for(int j=0;j<3;++j)
					deformed_ele_dis[i][j]=init_pos[3*global_idx[i]+j]+displacement[3*global_idx[i]+j];
			}
		}
		Vec3d *vec=new Vec3d[3];
		//compute deformed_displacement_matrix Ds
		for(int i=0;i<3;++i)
		{
			vec[i][0]=deformed_ele_dis[0][i]-deformed_ele_dis[3][i];
			vec[i][1]=deformed_ele_dis[1][i]-deformed_ele_dis[3][i];
			vec[i][2]=deformed_ele_dis[2][i]-deformed_ele_dis[3][i];
		}
		Mat3d deformed_ele_dis_matrix(vec[0],vec[1],vec[2]);
		Mat3d F=computeDeformationGradient(temp_init_matrix,deformed_ele_dis_matrix);
		if(ele==24320)
		{
	//		std::cout<<"object_init_element_dis_matrix_[ele]:"<<object_init_element_dis_matrix_[ele]<<"\n";
		// 	// std::cout<<"init_pos:\n";
		// 	// for(int i=0;i<mesh->getNumElementVertices();++i)
		// 	// {
		// 	// 	for(int j=0;j<3;++j)
		// 	// 	{
		// 	// 		std::cout<<init_pos[3*global_idx[i]+j]<<",";
		// 	// 	}
		// 	//
		// 	// }
		//
		// 		std::cout<<"displacement:\n";
		// 	for(int i=0;i<mesh->getNumElementVertices();++i)
		// 	{
		// 		for(int j=0;j<3;++j)
		// 		{
		// 			std::cout<<displacement[3*global_idx[i]+j]<<",";
		// 		}
		// 		//std::cout<<displacement[3*global_idx[i]+j+1]<<"extra,\n";
		// 	}
		//
		// 	std::cout<<"reduced-F:"<<det(F)<<"\n";
		// 	std::cout<<"temp_init_matrix:"<<temp_init_matrix<<"\n";
		// 	std::cout<<"deformed_ele_dis_matrix:"<<deformed_ele_dis_matrix<<"\n";
		}

		double lnJ=log(det(F));
		Mat3d P=mu*(F-trans(inv(F)))+lambda*lnJ*trans(inv(F));
		Mat3d H=ele_volume*P*inv(trans(temp_init_matrix));
		for(int i=0;i<mesh->getNumElementVertices();++i)
		 {
			for(int j=0;j<3;++j)
			{
				if(i==3)
					forces[3*global_idx[i]+j]+=(-1.0)*(H[j][0]+H[j][1]+H[j][2]);
				else
					forces[3*global_idx[i]+j]+=H[j][i];
			}
		}
		delete[] vec;
		delete[] deformed_ele_dis;
		delete[] global_idx;
		delete[] vert_pos;
	}
	// std::cout<<"compute g:\n";
	// for(int i=0;i<mesh->getNumVertices();++i)
	// 	if((fabs(forces[3*i])>1.0e-7)||(fabs(forces[3*i+1])>1.0e-7)||(fabs(forces[3*i+2])>1.0e-7))
	// 		std::cout<<i<<":"<<forces[3*i]<<",\n";
	//getchar();
	//project on reduced subspace
	if(example_flag==1)
	{
		for(int eigen_idx=0;eigen_idx<interpolate_eigenfunction_num_;++eigen_idx)
			for(int vert_idx=0;vert_idx<mesh->getNumVertices();++vert_idx)
			{
				g[3*eigen_idx]+=forces[3*vert_idx]*example_eigenfunctions_[dis_ex_idx][eigen_idx][vert_idx];
				g[3*eigen_idx+1]+=forces[3*vert_idx+1]*example_eigenfunctions_[dis_ex_idx][eigen_idx][vert_idx];
				g[3*eigen_idx+2]+=forces[3*vert_idx+2]*example_eigenfunctions_[dis_ex_idx][eigen_idx][vert_idx];
			}
	}
	else
	{
		for(int eigen_idx=0;eigen_idx<interpolate_eigenfunction_num_;++eigen_idx)
			for(int vert_idx=0;vert_idx<mesh->getNumVertices();++vert_idx)
			{
				g[3*eigen_idx]+=forces[3*vert_idx]*object_eigenfunctions_[eigen_idx][vert_idx];
				g[3*eigen_idx+1]+=forces[3*vert_idx+1]*object_eigenfunctions_[eigen_idx][vert_idx];
				g[3*eigen_idx+2]+=forces[3*vert_idx+2]*object_eigenfunctions_[eigen_idx][vert_idx];
				// if((vert_idx==1066)/*||(vert_idx==857)||(vert_idx==1073)||(vert_idx==1527)*/)
				// {
				// 	std::cout<<"idx:"<<vert_idx<<"...\n";
				// 	std::cout<<"fx:"<<eigen_idx<<":"<<forces[3*vert_idx]<<","<<g[3*eigen_idx]<<"...\n";
				// 	std::cout<<"fy:"<<eigen_idx<<":"<<forces[3*vert_idx+1]<<","<<g[3*eigen_idx+1]<<"...\n";
				// 	std::cout<<"fz:"<<eigen_idx<<":"<<forces[3*vert_idx+2]<<","<<g[3*eigen_idx+2]<<"...\n";
				//
				// // 	std::cout<<"g:"<<g[3*eigen_idx+1]<<"...\n";
				// // 	std::cout<<"g:"<<g[3*eigen_idx+2]<<"...\n";
				// }
			}
		// if((vert_idx==1066)&&(vert_idx==857)&&(vert_idx==1073)&&(vert_idx==1527))
		// {
		// 	std::cout<<forces[3]
		// }
		// std::cout<<"1066"<<forces[3*1066]<<","<<forces[3*1066+1]<<","<<forces[3*1066+2]<<"\n";
		// std::cout<<"857"<<forces[3*857]<<","<<forces[3*857+1]<<","<<forces[3*857+2]<<"\n";
	}

	//  for(int i=0;i<21;++i)
	//  	std::cout<<g[i]<<",";
	delete[] forces;
}

Mat3d RealTimeExampleBasedDeformer::computeDeformationGradient(const Mat3d &init_matrix,const Mat3d &deformed_matrix/*Vec3d *init_pos,Vec3d *deform_pos*/)
{
	// Mat3d F=deformed_matrix*inv(init_matrix);
    // return F;
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
void RealTimeExampleBasedDeformer::generateE() const
{
	//this matrix is later used to compute all the deformation gradients F from q in one shot.
	//the 3x3 F matrix is flattened in to a 9-vector

}
}  //namespace RTLB
