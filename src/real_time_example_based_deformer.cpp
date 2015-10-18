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
	// for(unsigned int i=0;i<example_num_;++i)
	// {
	// 	delete[] example_cubica_weights_[i];
	// 	delete[] example_cubica_elements_[i];
	// }
	// delete[] example_cubica_elements_[example_num_];
	// delete[] example_cubica_weights_;
	// delete[] example_cubica_elements_;
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
	delete[] F_;
	for(int i=0;i<object_cubica_ele_num_;++i)
		delete[] restpos_[i];
	delete[] restpos_;
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
	// std::cout<<reduced_row_num<<std::endl;
	// std::cout<<reduced_col_num<<std::endl;
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

	delete[] dis;

	isPreComputeReducedData_=true;
	if(isPreComputeReducedData_)
	{
		preComputeForReducedSimulation();
		isPreComputeReducedData_=false;
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

// bool RealTimeExampleBasedDeformer::loadPlanesInScene(const std::string &file_name, unsigned int plane_num)
// {
//     if(file_name=="")
// 	{
// 		std::cout<<"Error: cannot read "<<file_name<<".\n";
// 		return false;
// 	}
// 	planes_=new Planes(file_name.c_str(),plane_num);
//     if(planes_ == NULL)
//     {
//         std::cout<<"Error: failed to load planes from "<<file_name<<".\n";
//         return false;
//     }
//     return true;
// }

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
	for(int i=0;i<eigenfunction_num;++i)
		std::cout<<target_eigencoefs[i]<<",";
	std::cout<<std::endl;
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
	//std::cout<<"eeeeeeeeeeeeeeeeeeeeenalbe:"<<enable_eigen_weight_control_<<"\n";
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
	//std::cout<<"----------------------------------------------------weight compute end.\n";
//	getchar();
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
void RealTimeExampleBasedDeformer::testObjectiveGradients()
{
	real_1d_array x;

	// int num=interpolate_eigenfunction_num_*3;
	// double *temp_buffer=new double[num];
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	for(int j=0;j<3;++j)
	// 		temp_buffer[3*i+j]=example_eigencoefs_[0][i][j];
	// x.setcontent(num,temp_buffer);
	// for(int i=0;i<num;++i)
	// 	temp_buffer[i]=0.0;
	// double target_weights[2]={1.0,0.0};
	// double f=0.0;
	// real_1d_array grad;
	// grad.setcontent(num,temp_buffer);
	// evaluateObjectiveAndGradient2(x,f,grad,(void*)target_weights);
	// delete[] temp_buffer;



	int num=3*interpolate_eigenfunction_num_;
	double *temp_buffer=new double[num];
	srand((unsigned)time(0));
	int lowest=1,highest=10;
	int range=(highest-lowest)+1;
	//randomly generate initial x in range[0.1,1]
	//for(int i=0;i<num;++i)
		//temp_buffer[i]=(lowest+rand()%range)/10.0;
	for(int i=0;i<interpolate_eigenfunction_num_;++i)
		for(int j=0;j<3;++j)
			temp_buffer[3*i+j]=object_eigencoefs_[i][j];
	x.setcontent(num,temp_buffer);
	double f_min,f_plus,f;
	real_1d_array temp_grad,grad;
	temp_grad.setcontent(num,temp_buffer);
	grad.setcontent(num,temp_buffer);
	double target_weights[2]={0.5,0.5};
	evaluateObjectiveAndGradient2(x,f,grad,(void*)target_weights);
	for(int i=0;i<num;++i)
	{
		temp_buffer[i]=(lowest+rand()%range)/1.0e7;
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
	std::cout<<"Objective, df, analytic: "<<setprecision(15)<<f_grad_times_dx<<", numerial: "<<setprecision(15)<<f_plus-f_min;
	std::cout<<", rel_error= "<<(f_plus-f_min-f_grad_times_dx)/(fabs(f_grad_times_dx)>1e-20?fabs(f_grad_times_dx):1e-20)<<"\n";
	delete[] temp_buffer;
}
//step2, minimization of wieghted deformation energy to the examples
void RealTimeExampleBasedDeformer::evaluateObjectiveAndGradient2(const real_1d_array &x,double &func, real_1d_array &grad, void *ptr)
{
	RealTimeExampleBasedDeformer* active_instance=RealTimeExampleBasedDeformer::activeInstance();
	assert(active_instance);
	for(int i=0;i<3*active_instance->interpolate_eigenfunction_num_;++i)
	{
		std::cout<<"x:"<<x[i]<<",";
	}

	// for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
	// {
	// 	//for(int j=0;j<3;++j)
	// 	//std::cout<<"x:"<<x[i]<<",";
	// 	//x[3*i+j]=active_instance->object_eigencoefs_[i][j];
	// 	std::cout<<"e:"<<active_instance->object_eigencoefs_[i][0]<<","<<active_instance->object_eigencoefs_[i][1]<<","<<active_instance->object_eigencoefs_[i][2]<<",\n";
	// }
	// for(int i=0;i<active_instance->interpolate_eigenfunction_num_;++i)
	// {
	// 	//for(int j=0;j<3;++j)
	// 	//std::cout<<"x:"<<x[i]<<",";
	// 	//x[3*i+j]=active_instance->object_eigencoefs_[i][j];
	// 	std::cout<<"e:"<<active_instance->example_eigencoefs_[0][i][0]<<","<<active_instance->example_eigencoefs_[0][i][1]<<","<<active_instance->example_eigencoefs_[0][i][2]<<",\n";
	// }
	func=0.0;
	for(int j=0;j<3*active_instance->interpolate_eigenfunction_num_;++j)
		grad[j]=0.0;
	double *target_weights=(double*)ptr;
	Vec3d *temp_eigencoefs=new Vec3d[active_instance->interpolate_eigenfunction_num_];
	for(int i=0;i<active_instance->example_num_;++i)
	{
		std::cout<<"example"<<i<<std::endl;
		std::cout<<"temp_eigencoefs:"<<std::endl;
	// 	//compute displacement in LB space for each example
		for(int j=0;j<active_instance->interpolate_eigenfunction_num_;++j)
		{
			for(int k=0;k<3;++k)
			{
		//	std::cout<<active_instance->example_eigencoefs_[i][j][k]<<",";
				temp_eigencoefs[j][k]=x[3*j+k]-active_instance->example_eigencoefs_[i][j][k];

			}
			std::cout<<temp_eigencoefs[j][0]<<","<<temp_eigencoefs[j][1]<<","<<temp_eigencoefs[j][2]<<",\n";
		}
		active_instance->computeF(temp_eigencoefs);
		std::cout<<"~~~~~~~~~~~~~~~~~~"<<active_instance->F_[0]<<std::endl;
		double energy=0.0;
		active_instance->computeReducedEnergy(temp_eigencoefs,energy);
		func+=target_weights[i]*energy;
		std::cout<<"func:energy done, value is "<<func<<" .\n";

	//	getchar();
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
	//compute gradient force
}
void RealTimeExampleBasedDeformer::preComputeForReducedSimulation()
{
	//precompute for cubica simulation
	//generate cubica reduced E and H
//	generateE();
//	generateH();
	//compute restposition for all cubica elements
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
	F_ = new Mat3d[object_cubica_ele_num_];
}
//compute basis matrix Du
Matrix RealTimeExampleBasedDeformer::vertexSubBasis(const int &vert_idx) const
{//vertex_subBasis is 1*m
	Matrix vertex_subBasis(1,interpolate_eigenfunction_num_);
//	std::cout<<"interpolate_eigenfunction_num_:"<<interpolate_eigenfunction_num_<<"\n";
//	std::cout<<"object_eigenfunctions_[i][vert_idx]:"<<object_eigenfunctions_[0][vert_idx]<<",";
	for(int i=0;i<interpolate_eigenfunction_num_;++i)
	{
		vertex_subBasis(1,i+1)=object_eigenfunctions_[i][vert_idx];
//		std::cout<<"vertex_subBasis:"<<vertex_subBasis(1,i+1)<<",";
	}
	return vertex_subBasis;
}
Matrix RealTimeExampleBasedDeformer::tetSubBasis(const int &ele) const
{//tet_subBasis is 12*m
	Matrix tet_subBasis(12,interpolate_eigenfunction_num_);
	int *global_idx=new int[4];
	for(int i=0;i<4;++i)
		global_idx[i]=simulation_mesh_->getVertexIndex(ele,i);
	for(int j=0;j<4;++j)
	{
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
		{
			tet_subBasis(3*j+1,i+1)=object_eigenfunctions_[i][global_idx[j]];
			tet_subBasis(3*j+2,i+1)=object_eigenfunctions_[i][global_idx[j]];
			tet_subBasis(3*j+3,i+1)=object_eigenfunctions_[i][global_idx[j]];
		//	std::cout<<tet_subBasis(3*j+1,i+1)<<","<<tet_subBasis(3*j+2,i+1)<<","<<tet_subBasis(3*j+3,i+1)<<"\n";
		}
	}
	delete[] global_idx;
	return tet_subBasis;
}
void RealTimeExampleBasedDeformer::test()
{
	//---------------1---------------------test computeDs:
	// double *dis=new double[12];
	// memset(dis,0.0,sizeof(double)*12);
	// for(int i=1;i<=12;++i)
	// {
	// 	dis[i-1]=1.0*i;
	// 	std::cout<<dis[i]<<",";
	// }
	// dis[2]=3.3;dis[4]=5.1;dis[5]=6.2;dis[8]=9.5;
	// dis[9]=5.0;dis[10]=7.0;dis[11]=8.0;
	// Mat3d F=computeDs(dis);
	// delete[] dis;
	//--------------2------------------test computeDmInv
	//Mat3d DmInv=computeDmInv(0);
	//--------------3------------------test computeF
	Vec3d *reduced_dis=new Vec3d[2];
	reduced_dis[0][0]=1.0;reduced_dis[0][1]=2.0;reduced_dis[0][2]=3.0;
	reduced_dis[1][0]=2.1;reduced_dis[1][1]=3.2;reduced_dis[1][2]=4.5;
	computeF(reduced_dis);
	delete[] reduced_dis;


}
Mat3d RealTimeExampleBasedDeformer::computeDs(const double *reduced_dis) const
{//reduced dis is 12*1 for each element
	//std::cout<<"reduced_dis:\n";
	// for(int i=0;i<12;++i)
	// 	std::cout<<reduced_dis[i]<<",";
	Mat3d Ds(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
	{
		Ds[0][i]=reduced_dis[3*i]-reduced_dis[9];
		Ds[1][i]=reduced_dis[3*i+1]-reduced_dis[10];
		Ds[2][i]=reduced_dis[3*i+2]-reduced_dis[11];
	}
	//std::cout<<Ds<<"\n";
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
	//	std::cout<<"vert_pos:"<<vert_pos[j]<<",";
	}
	//std::cout<<"\n";
	Mat3d Dm(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
	{

		Dm[i][0]=vert_pos[0][i]-vert_pos[3][i];
		Dm[i][1]=vert_pos[1][i]-vert_pos[3][i];
		Dm[i][2]=vert_pos[2][i]-vert_pos[3][i];
		// std::cout<<"DmInv("<<i+1<<",1):"<<Dm(i+1,1)<<",";
		// std::cout<<"DmInv("<<i+1<<",2):"<<Dm(i+1,2)<<",";
		// std::cout<<"DmInv("<<i+1<<",3):"<<Dm(i+1,3)<<",";
	}
//	std::cout<<"Dm"<<Dm<<"\n";
	Mat3d DmInv=inv(Dm);
//	std::cout<<"DmInv"<<DmInv<<"\n";
	delete[] global_idx;
	delete[] vert_pos;
//	std::cout<<"computeDmInv-Done\n";
	return DmInv;
}
void RealTimeExampleBasedDeformer::computeF(const Vec3d *reduced_dis) const
{//for all cubica elements
	Matrix reduced_dis_matrix(interpolate_eigenfunction_num_,3);//rx3
	// std::cout<<"reduced_dis:"<<"\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<reduced_dis[i][0]<<","<<reduced_dis[i][1]<<","<<reduced_dis[i][2]<<"\n";
	for(int i=0;i<interpolate_eigenfunction_num_;++i)
	{
		reduced_dis_matrix(i+1,1)=reduced_dis[i][0];
		reduced_dis_matrix(i+1,2)=reduced_dis[i][1];
		reduced_dis_matrix(i+1,3)=reduced_dis[i][2];
	}
	for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
	{
		int ele=object_cubica_elements_[cubica_idx];
		//compute displacement x=X+Uq
		double *deformed=new double[12];
		memset(deformed,0.0,sizeof(double)*12);
		for(int j=0;j<4;++j)
		{
			int vertID=simulation_mesh_->getVertexIndex(ele,j);
			Matrix subU=vertexSubBasis(vertID);//1xr
			// for(int i=0;i<interpolate_eigenfunction_num_;++i)
			// std::cout<<"subu:"<<subU(1,i+1)<<","<<"\n";
			Matrix temp=subU*reduced_dis_matrix;//1xr,rx3->1x3
			// for(int i=0;i<3;++i)
			// std::cout<<"temp:"<<temp(1,i+1)<<","<<"\n";
			deformed[3*j]=restpos_[cubica_idx][3*j]+temp(1,1);
			deformed[3*j+1]=restpos_[cubica_idx][3*j+1]+temp(1,2);
			deformed[3*j+2]=restpos_[cubica_idx][3*j+2]+temp(1,3);
		}

		Mat3d Ds=computeDs(deformed);
//		std::cout<<"Ds:"<<Ds<<"\n";
//		std::cout<<"Dm"<<computeDmInv(ele)<<"\n";
		F_[cubica_idx]=Ds*computeDmInv(ele);//F for each ele
	//	std::cout<<"F:"<<F_[cubica_idx]<<"\n";
//	std::cout<<"detF:"<<det(F_[cubica_idx])<<"~~~~~~~~~\n";
		delete[] deformed;
	}
}
Mat3d RealTimeExampleBasedDeformer::firstPiolaKirchhoff(Mat3d &F) const
{//3*3
	Mat3d P(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	P=mu_*(F-trans(inv(F)))+lamda_*log(det(F))*trans(inv(F));
	return P;
}


void RealTimeExampleBasedDeformer::computeReducedInternalForce(const Vec3d *reduced_dis,double *forces)
{//r*1
	// memset(force,0.0,sizeof(double)*interpolate_eigenfunction_num_);
	// //memset(reduced_dis,0.0,sizeof(double)*interpolate_eigenfunction_num_);
	// Matrix g(9*object_cubica_ele_num_,1);//9n*1
	// double *reduced_F=new double[9*object_cubica_ele_num_];
	// memset(reduced_F,0.0,sizeof(Matrix)*9*object_cubica_ele_num_);
	// computeReducedF(reduced_dis,reduced_F);//generate: reduced_F
	// for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
	// {
	// 	int ele=object_cubica_elements_[cubica_idx];
	// 	//computation of F not done yet
	// 	Mat3d F;
	// 	for(int i=0;i<3;++i)
	// 		for(int j=0;j<3;++j)
	// 			F[i][j]=reduced_F[9*cubica_idx+3*i+j];
	// 	//forceDensity is the flatten P(P is the first PiolaKirchhoff)
	// 	Mat3d P=firstPiolaKirchhoff(F);
	// 	double *flat_P=new double[9];
	// 	flatten(P,flat_P);//9*1
	//
	// 	for(int i=0; i<9; ++i)
	// 		g(9*cubica_idx+i,1)=flat_P[i];
	// 	delete[] flat_P;
	// }
	// const Matrix temp_matrix=H_.t()*g;//H_:9n*r,g:9n*1->get temp_matrix:r*1
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// {
	// 	force[3*i]=temp_matrix(i,1);
	// 	force[3*i+1]=temp_matrix(i,1);
	// 	force[3*i+2]=temp_matrix(i,1);
	// }
	//
	// delete[] reduced_F;
//	memset(forces,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
	//computeF(reduced_dis);
	for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
	{
		int ele=object_cubica_elements_[cubica_idx];
		Mat3d P=firstPiolaKirchhoff(F_[cubica_idx]);
		Mat3d temp=P*trans(computeDmInv(ele));
		//compute P*trans(DmInv) for all vertices of each element
		//temp_matrix is 3x4, each column represents a vertex, each row denotes a direction
		Matrix temp_matrix(3,4);
		for(int i=1;i<=3;++i)
		{
			for(int j=1;j<=4;++j)
			{
				if(j==4)
					temp_matrix(i,j)=(-1.0)*(temp[i-1][0]+temp[i-1][1]+temp[i-1][2]);
				else
					temp_matrix(i,j)=temp[i-1][j-1];
			}
		}
		Matrix subU=tetSubBasis(ele);//4xr
		double *g=new double[3*interpolate_eigenfunction_num_];//3rx1
		memset(g,0.0,sizeof(double)*3*interpolate_eigenfunction_num_);
		//for eache ele:g:nxr; g_3i=-sum(U_vert x^i*(P*Dm^-T)_vertex x^k),k=x,y,z direction
		for(int i=0;i<interpolate_eigenfunction_num_;++i)
		{
			for(int j=0;j<4;++j)
			{
				g[3*i]+=subU(j+1,i+1)*temp_matrix(1,j+1);
				g[3*i+1]+=subU(j+1,i+1)*temp_matrix(2,j+1);
				g[3*i+2]+=subU(j+1,i+1)*temp_matrix(3,j+1);
			}
		}
		for(int i=0;i<3*interpolate_eigenfunction_num_;++i)
			forces[i] += g[i];
		delete[] g;
	}
}
//compute for E(deformed,source)
//1: mesh is the source mesh, first get the mesh material of source
//2: souce pos could be Null, then get the source pos from the given mesh
//3: get the cubica_elements from source mesh
//4: get the vertices position (deformed pos and source pos) for the corresponding elements
void RealTimeExampleBasedDeformer::computeReducedEnergy(const Vec3d *reduced_dis,double &energy)
{
//	std::cout<<"computeReducedEnergy:\n";
	energy=0.0;
//	std::cout<<"---reduced_dis:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<reduced_dis[i][0]<<","<<reduced_dis[i][1]<<","<<reduced_dis[i][2]<<",";
//	computeF(reduced_dis);//compute F for all cubica elements
	//std::cout<<"computeF done\n";
	for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
	{
		int ele=object_cubica_elements_[cubica_idx];
		//std::cout<<"F:\n";
		//std::cout<<F_[cubica_idx]<<",\n";

		Mat3d temp=trans(F_[cubica_idx])*F_[cubica_idx];
		double trace_c=temp[0][0]+temp[1][1]+temp[2][2];
		double lnJ=log(det(F_[cubica_idx]));
		double element_energy=0.5*mu_*(trace_c-3)-mu_*lnJ+0.5*lamda_*lnJ*lnJ;
		energy += object_cubica_weights_[cubica_idx]*element_energy;
	//	std::cout<<"energy:"<<energy<<",";
	//	getchar();
	}
	// std::cout<<"energy:"<<energy<<"c\n";
	// std::cout<<"computeReducedEnergy:-end\n";
	// std::cout<<"computeReducedEnergy:\n";
	// energy=0.0;
	// double *reduced_F=new double[9*object_cubica_ele_num_];
	// memset(reduced_F,0.0,sizeof(Matrix)*9*object_cubica_ele_num_);
	// std::cout<<"---reduced_dis:\n";
	// for(int i=0;i<interpolate_eigenfunction_num_;++i)
	// 	std::cout<<reduced_dis[i][0]<<","<<reduced_dis[i][1]<<","<<reduced_dis[i][2]<<",";
	// computeReducedF(reduced_dis,reduced_F);//generate: reduced_F
	// std::cout<<"reduced_F:\n";
	// for(int i=0;i<object_cubica_ele_num_;++i)
	// 		for(int k=0;k<9;++k)
	// 			std::cout<<reduced_F[9*object_cubica_ele_num_+i]<<",";
	// std::cout<<"computeReducedF done\n";
	// for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
	// {
	// 	int ele=object_cubica_elements_[cubica_idx];
	// 	//computation of F not done yet
	// 	Mat3d F;
	// 	for(int i=0;i<3;++i)
	// 		for(int j=0;j<3;++j)
	// 			F[i][j]=reduced_F[9*cubica_idx+3*i+j];
	// 	Mat3d temp=trans(F)*F;
	// 	double trace_c=temp[0][0]+temp[1][1]+temp[2][2];
	// 	double lnJ=log(det(F));
	// 	double element_energy=0.5*mu_*(trace_c-3)-mu_*lnJ+0.5*lamda_*lnJ*lnJ;
	// 	energy += object_cubica_weights_[cubica_idx]*element_energy;
	// 	std::cout<<"energy:"<<energy<<",";
	// }
	// delete[] reduced_F;
	// std::cout<<"energy:"<<energy<<"c\n";
	// std::cout<<"computeReducedEnergy:-end\n";

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

// void RealTimeExampleBasedDeformer::computepFpu(const int &ele,Matrix &pFpu) const
// {//9*12 matrix
// 	const Matrix matInv=computeDmInv(ele);
// 	const double m = matInv(1,1);
// 	const double n = matInv(1,2);
// 	const double o = matInv(1,3);
// 	const double p = matInv(2,1);
// 	const double q = matInv(2,2);
// 	const double r = matInv(2,3);
// 	const double s = matInv(3,1);
// 	const double t = matInv(3,2);
// 	const double u = matInv(3,3);
//
// 	const double t1 = -m-p-s;
// 	const double t2 = -n-q-t;
// 	const double t3 = -o-r-u;
// 	for(int i=1;i<=9;++i)
// 		for(int j=1;j<=12;++j)
// 			pFpu(i,j)=0.0;
// 	pFpu(1,1)=pFpu(2,2)=pFpu(3,3)=m;
// 	pFpu(1,4)=pFpu(2,5)=pFpu(3,6)=p;
// 	pFpu(1,7)=pFpu(2,8)=pFpu(3,9)=s;
// 	pFpu(1,10)=pFpu(2,11)=pFpu(3,12)=t1;
// 	pFpu(4,1)=pFpu(5,2)=pFpu(6,3)=n;
// 	pFpu(4,4)=pFpu(5,5)=pFpu(6,6)=q;
// 	pFpu(4,7)=pFpu(5,8)=pFpu(6,9)=t;
// 	pFpu(4,10)=pFpu(5,11)=pFpu(6,12)=t2;
// 	pFpu(7,1)=pFpu(8,2)=pFpu(9,3)=o;
// 	pFpu(7,4)=pFpu(8,5)=pFpu(9,6)=r;
// 	pFpu(7,7)=pFpu(8,8)=pFpu(9,9)=u;
// 	pFpu(7,10)=pFpu(8,11)=pFpu(9,12)=t3;
// }
// void RealTimeExampleBasedDeformer::generateH()
// {//H_:9n*r
// 	std::cout<<"generateH:\n";
// 	H_.ReSize(9*object_cubica_ele_num_,interpolate_eigenfunction_num_);
// 	for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
// 	{
// 		int ele=object_cubica_elements_[cubica_idx];
// 		Matrix pFpu(9,12);
// 		computepFpu(ele,pFpu);//9*12
// 		std::cout<<"pFpu done\n";
// 		Matrix subbasis = tetSubBasis(ele);//12*r
// 		std::cout<<"subbasis done\n";
// 		Matrix temp_matrix(9,interpolate_eigenfunction_num_);
// 		temp_matrix=object_cubica_weights_[ele]*pFpu*subbasis;//9*r
// 		for(int i=1;i<=9;++i)
// 		{
// 			for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 			{
// 				H_(9*cubica_idx+i,j+1)=temp_matrix(i,j+1);
// 			}
// 		}
// 	}
// 	std::cout<<"generateH-end:\n";
// }
// void RealTimeExampleBasedDeformer::computeReducedF(const Vec3d *reduced_dis,double *reduced_F) const
// {//9n*1 vector,F=I+Eq, E:9n*r, q:r*1
// 	memset(reduced_F,0.0,sizeof(double)*9*object_cubica_ele_num_);
// 	for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
// 	{
// 		reduced_F[9*cubica_idx] += 1.0;
// 		reduced_F[9*cubica_idx+4] += 1.0;
// 		reduced_F[9*cubica_idx+8] += 1.0;
// 	}
// 	// std::cout<<"~~~~~~~~E_:\n";
// 	// for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
// 	// {
// 	// 	for(int i=0;i<3;++i)
// 	// 	{
// 	// 		for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 	// 			std::cout<<E_(3*cubica_idx+i+1,j+1)<<"\n";
// 	// 	}
// 	// }
// 	std::cout<<"FFFFFFFFFFFFFFFFFF:reduced_dis:\n";
// 	for(int i=0;i<interpolate_eigenfunction_num_;++i)
// 		std::cout<<reduced_dis[i][0]<<","<<reduced_dis[i][1]<<","<<reduced_dis[i][2]<<",";
// 	for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
// 	{
// 		//int ele=object_cubica_elements_[cubica_idx];
// 		//E_*reduced_dis to make a 9*1 vector for each element
// 		double *temp_value=new double[9];
// 		memset(temp_value,0.0,sizeof(double)*9);
// 		for(int i=0;i<3;++i)
// 		{
// 			for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 			{
// 				std::cout<<"E:"<<E_(3*cubica_idx+i+1,j+1)<<"\n";
// 				std::cout<<"reduced_dis:"<<reduced_dis[j][0]<<"\n";
// 					std::cout<<"reduced_dis:"<<reduced_dis[j][0]<<"\n";
// 						std::cout<<"reduced_dis"<<reduced_dis[j][0]<<"\n";
// 				temp_value[3*i]+=E_(3*cubica_idx+i+1,j+1)*reduced_dis[j][0];
// 				temp_value[3*i+1]+=E_(3*cubica_idx+i+1,j+1)*reduced_dis[j][1];
// 				temp_value[3*i+2]+=E_(3*cubica_idx+i+1,j+1)*reduced_dis[j][2];
// 				std::cout<<temp_value[3*i]<<","<<temp_value[3*i+1]<<","<<temp_value<<",";
// 			}
// 		}
// 		std::cout<<"~~~~~~~~~~c\n";
// 		for(int i=0;i<9;++i)
// 		//	for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 		{
// 			reduced_F[9*cubica_idx+i] += temp_value[i];
// 			std::cout<<"reduced_F:"<<reduced_F[9*cubica_idx+i]<<",";
// 		}
//
// 		delete[] temp_value;
// 	}
// }
// void RealTimeExampleBasedDeformer::flatten(Mat3d &mat,double *flat_mat) const
// {//mat3d to 9x1 vector
// 	for(int i=0;i<3;++i)
// 	{
// 		flat_mat[3*i+0]=mat[i][0];
// 		flat_mat[3*i+1]=mat[i][1];
// 		flat_mat[3*i+2]=mat[i][2];
// 	}
// }
// void RealTimeExampleBasedDeformer::reback(const double *flat_mat, Mat3d &mat)
// {//9x1 vector to matrix(3,3)
// 	for(int i=0;i<3;++i)
// 	{
// 		mat[i][0]=flat_mat[3*i];
// 		mat[i][1]=flat_mat[3*i+1];
// 		mat[i][2]=flat_mat[3*i+2];
// 	}
// }

//compute basis matrix Du
// void RealTimeExampleBasedDeformer::generateE()
// {//9nxr
// 	std::cout<<"generateE:\n";
// 	E_.ReSize(3*object_cubica_ele_num_,interpolate_eigenfunction_num_);
// 	std::cout<<"object_cubica_ele_num_:\n";
// 	for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
// 	{
// 		int ele=object_cubica_elements_[cubica_idx];
// 		std::cout<<"ele:"<<ele<<"\n";
// 		Matrix DmInv=computeDmInv(ele);
// 		std::cout<<"DmInv:"<<DmInv(1,1)<<"\n";
// 		std::cout<<simulation_mesh_->getVertexIndex(ele,0)<<":\n";
// 		Matrix subBasis0=vertexSubBasis(simulation_mesh_->getVertexIndex(ele,0));
// 		std::cout<<simulation_mesh_->getVertexIndex(ele,1)<<":\n";
// 		Matrix subBasis1=vertexSubBasis(simulation_mesh_->getVertexIndex(ele,1));
// 		std::cout<<simulation_mesh_->getVertexIndex(ele,2)<<":\n";
// 		Matrix subBasis2=vertexSubBasis(simulation_mesh_->getVertexIndex(ele,2));
// 		std::cout<<simulation_mesh_->getVertexIndex(ele,3)<<":\n";
// 		Matrix subBasis3=vertexSubBasis(simulation_mesh_->getVertexIndex(ele,3));
//
// 		std::cout<<"subBasis0:"<<subBasis0(1,1)<<"\n";
// 		std::cout<<"subBasis1:"<<subBasis1(1,1)<<"\n";
// 		std::cout<<"subBasis2:"<<subBasis2(1,1)<<"\n";
// 		std::cout<<"subBasis3:"<<subBasis3(1,1)<<"\n";
// 		subBasis0 -= subBasis3;
// 		subBasis1 -= subBasis3;
// 		subBasis2 -= subBasis3;
// 		Matrix result(1,interpolate_eigenfunction_num_);
// 		for (int i = 1; i < 4; i++)
// 	    {
// 			result = DmInv(1,i)*subBasis0+DmInv(2,i)*subBasis1+DmInv(3,i)*subBasis2;
// 			for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 				std::cout<<result(1,j+1)<<"\n";
// 	      // each tet takes up nine rows, with three rows per subbasis
// 			for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 			{
// 				E_(3*cubica_idx+i,j+1)=result(1,j+1);
// 			}
//
// 	    }
// 	}
// 	for(int cubica_idx=0;cubica_idx<object_cubica_ele_num_;++cubica_idx)
// 		{
// 			for (int i = 1; i < 4; i++)
// 			{
// 				for(int j=0;j<interpolate_eigenfunction_num_;++j)
// 					std::cout<<E_(3*cubica_idx+i,j+1)<<",";
// 						std::cout<<"\n";
// 			}
// 				std::cout<<"-----------------\n";
// 		}
// 	std::cout<<"generateE-end:\n";
// }
// double RealTimeExampleBasedDeformer::computeReducedDis(double *reduced_dis) const
// {
//
// 	return dis;
// }
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

// bool RealTimeExampleBasedDeformer::loadExampleCubicaData(const std::string &file_name_prefix)
// {
// 	std::cout<<"Load cubica data for examples...\n";
// 	if(example_num_==0)
// 	{
// 		std::cout<<"examle cubica file unloaded.\n";
// 		return false;
// 	}
// 	example_cubica_ele_num_ = new unsigned int[example_num_];
// 	example_cubica_elements_ = new unsigned int*[example_num_];
// 	example_cubica_weights_ = new double*[example_num_];
// 	for(unsigned int i=0;i<example_num_;++i)
// 	{
//         std::string file_num_str,file_name;
//         std::stringstream adaptor;
// 		adaptor.str("");
// 		adaptor.clear();
// 		adaptor<<i;
// 		adaptor>>file_num_str;
// 		//read file, file format is .eigen,
// 		file_name=file_name_prefix+file_num_str+std::string(".cubature");
//         std::fstream input_file(file_name.c_str());
// 		if(!input_file)
// 		{
// 			std::cout<<"Error: Cannot open file "<<file_name<<std::endl;
// 			return false;
// 		}
// 		string temp_str;
// 		std::getline(input_file,temp_str);
// 		example_cubica_ele_num_[i]=atoi(temp_str.c_str());
// 		example_cubica_elements_[i] = new unsigned int[example_cubica_ele_num_[i]];
// 		example_cubica_weights_[i] = new double [example_cubica_ele_num_[i]];
// 		while(std::getline(input_file,temp_str))
// 		{
// 			if(temp_str.compare(0,6,string("*tetID"))==0)
// 				break;
// 		}
// 		unsigned int str_num=0;
// 		while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
// 		{
// 			unsigned int temp_value;
// 			input_file>>temp_value;
// 			example_cubica_elements_[i][str_num]=temp_value;
// 			str_num++;
// 			if(str_num>=example_cubica_ele_num_[i])
// 				break;
// 		}
// 		while(std::getline(input_file,temp_str))
// 		{
// 			if(temp_str.compare(0,10,string("*tetWeight"))==0)
// 				break;
// 		}
// 		str_num=0;
// 		while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
// 		{
// 			double temp_value;
// 			input_file>>temp_value;
// 			example_cubica_weights_[i][str_num]=temp_value;
// 			++str_num;
// 			if(str_num>=example_cubica_ele_num_[i])
// 				break;
// 		}
// 	}
// 	std::cout<<"Load example cubica data succeed.\n";
// 	//isload_example_cubica_=true;
// 	return true;
// }
}  //namespace RTLB
