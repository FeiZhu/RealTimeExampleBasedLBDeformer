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
#include "sceneObjectDeformable.h"
#include "coupled_quasi_harmonics.h"
#include "planes.h"
#include "real_time_example_based_deformer.h"

namespace RTLB{

RealTimeExampleBasedDeformer::RealTimeExampleBasedDeformer()
{
    //TO DO
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
	for(unsigned int i=0;i<object_eigenfunction_num_;++i)
		if(object_eigenfunctions_[i])
			delete[] object_eigenfunctions_[i];
	for(unsigned int i=0;i<example_num_;++i)
	{
		if(example_eigenfunctions_[i])
		{
			for(unsigned int j=0;j<example_eigenfunction_num_;++j)
				delete[] example_eigenfunctions_[i][j];
			delete[] example_eigenfunctions_[i];
		}
		delete[] example_eigenvalues_[i];
		delete[] example_eigencoefs_[i];
		delete[] example_vertex_volume_[i];
	}
	delete[] example_volume_;
	delete[] target_eigencoefs_;
	//delete registration information
	for(unsigned int i=0;i<corresponding_function_num_;++i)
		if(object_corresponding_functions_[i])
			delete[] object_corresponding_functions_[i];
	for(unsigned int i=0;i<example_num_;++i)
	{
		if(example_corresponding_functions_[i])
		{
			for(unsigned int j=0;j<example_eigenfunction_num_;++j)
				delete[] example_corresponding_functions_[i][j];
			delete[] example_corresponding_functions_[i];
		}
	}
	if(planes_==NULL)
		delete[] planes_;

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
			unsigned int global_idx=simulation_mesh_->getVertexIndex(ele_idx,i);
			object_vertex_volume_[global_idx]+=simulation_mesh_->getElementVolume(ele_idx);
		}
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
    std::string object_eigenfunction_num_str;
	std::getline(input_file,object_eigenfunction_num_str);
	object_eigenfunction_num_=atoi(object_eigenfunction_num_str.c_str());
	object_eigenvalues_=new double[object_eigenfunction_num_];
	unsigned int str_num=0;
	while((!input_file.eof())&&(input_file.peek()!=std::ifstream::traits_type::eof()))
	{
		double temp_value;
		input_file>>temp_value;
		object_eigenvalues_[str_num]=temp_value;
		str_num++;
		if(str_num>=object_eigenfunction_num_)
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
		if(str_num>=object_eigenfunction_num_-1)
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
        std::string example_eigenfunction_num_str;
		std::getline(input_file,example_eigenfunction_num_str);
		example_eigenfunction_num_=atoi(example_eigenfunction_num_str.c_str());
		example_eigenvalues_[ex_num]=new double[example_eigenfunction_num_];
		for(unsigned int j=0;j<example_eigenfunction_num_;++j)
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
			if(str_num>=example_eigenfunction_num_)
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
			if(str_num>=example_eigenfunction_num_-1)
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
	}
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

bool RealTimeExampleBasedDeformer::loadObjectCubicaData(const std::string &file_name)
{
	//TO DO
	return false;
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
	if(!output_file)
	{
		std::cout<<"Error:unable to open file "<<file_name<<".\n";
		return false;
	}
	//write eigenvalues first
	output_file<<"*eigenValues"<<std::endl;
	output_file<<object_eigenfunction_num_<<std::endl;
	for(unsigned int i=0;i<object_eigenfunction_num_;++i)
	{
		output_file<<object_eigenvalues_[i]<<" ";
	}
	output_file<<std::endl;
	//write eigenvectors
	output_file<<"*eigenVectors"<<std::endl;
	unsigned int vert_num=simulation_mesh_->getNumVertices();
	output_file<<vert_num<<std::endl;
	output_file<<object_eigenfunction_num_<<std::endl;
	for(unsigned int i=0;i<vert_num;++i)
	{
		for(unsigned int j=0;j<object_eigenfunction_num_;++j)
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
		output_file<<example_eigenfunction_num_<<std::endl;
		for(unsigned int j=0;j<example_eigenfunction_num_;++j)
		{
			output_file<<example_eigenvalues_[i][j]<<" ";
		}
		output_file<<std::endl;
		//write eigenvectors
		output_file<<"*eigenVectors"<<std::endl;
		unsigned int vert_num=examples_[i]->getNumVertices();
		output_file<<vert_num<<std::endl;
		output_file<<example_eigenfunction_num_<<std::endl;
		for(unsigned int j=0;j<vert_num;++j)
		{
			for(unsigned int k=0;k<example_eigenfunction_num_;++k)
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

const VolumetricMesh* RealTimeExampleBasedDeformer::exampleMesh(unsigned int example_idx) const
{
	assert(example_idx < example_num_);
	return examples_[example_idx];
}

//file format of corresponding function file:
//first line is an integer p indicating the number of corresponding functions (regions on the object)
//the following are (1+example_pos_num)*p  lines, every (1+example_pose_num) consecutive lines are a group
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
			object_corresponding_functions_[idx][i]=1;
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
				example_corresponding_functions_[ex_idx][idx][i]=1;
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
		std::cout<<"Error: simulation tet mesh is not loaded.\n";
		return false;
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
		std::cout<<"Error: object corresponding functions is not loaded.\n";
		return false;
	}
	int row_b =simulation_mesh_->getNumVertices();
	object_vertex_volume_=new double [row_b];
	int col_b = object_eigenfunction_num_;
	//change to col order,each col is an eigenfunction
	double **obj_eigenfunction_col=new double *[row_b];
	double **new_obj_eigenfunction_col=new double *[row_b];
	for(unsigned int i=0;i<row_b;++i)
	{
		obj_eigenfunction_col[i]=new double[col_b];
		new_obj_eigenfunction_col[i]=new double[col_b];
	}
	for(unsigned int i=0;i<col_b;++i)
	{
		for(unsigned int j=0;j<row_b;++j)
		{
			obj_eigenfunction_col[j][i]=object_eigenfunctions_[i][j];
			//std::cout<<obj_eigenfunction_col[j][i]<<",";
		}
		//std::cout<<endl;
	}
	for(unsigned int i=0;i<example_num_;++i)
	{
		double **G=object_corresponding_functions_,**F=example_corresponding_functions_[i];
		int p=corresponding_function_num_;
		int row_a=examples_[i]->getNumVertices();
		int col_a=example_eigenfunction_num_;
		double **ex_eigenfunction_col=new double *[row_a];
		double **new_ex_eigenfunction_col=new double *[row_a];
		for(unsigned int j=0;j<row_a;++j)
		{
			ex_eigenfunction_col[j]=new double[col_a];
			new_ex_eigenfunction_col[j]=new double[col_a];
		}
		for(unsigned int j=0;j<col_a;++j)
		{
			for(unsigned int k=0;k<row_a;++k)
			{
				ex_eigenfunction_col[k][j]=example_eigenfunctions_[i][j][k];
				//std::cout<<ex_eigenfunction_col[k][j]<<",";
			}
			//std::cout<<endl;
		}
		example_vertex_volume_[i]=new double [row_a];

		int col=col_a<col_b?col_a:col_b;
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
		quasi_base_generator.getCoupledBases(corresponding_function_num_,new_ex_eigenfunction_col,new_obj_eigenfunction_col);
		std::cout<<"Done.\n";

		//transform col-order eigen-bases to row-order so that each row is one eigenfunction
		for(unsigned int j=0;j<row_a;++j)
			for(unsigned int k=0;k<corresponding_function_num_;++k)
			{
				example_eigenfunctions_[i][j][k]=new_ex_eigenfunction_col[j][k];
				//std::cout<<new_ex_eigenfunction_col[j][k]<<",";
			}


		for(unsigned int j=0;j<row_a;++j)
		{
			delete[] ex_eigenfunction_col[j];
			delete[] new_ex_eigenfunction_col[j];
		}
		delete[] ex_eigenfunction_col;
		delete[] new_ex_eigenfunction_col;
	}


	for(unsigned int i=0;i<example_num_;++i)
	{
		for(unsigned int j=0;j<example_eigenfunction_num_;++j)
		{
			for(unsigned int k=0;k<10;++k)
			{
				std::cout<<example_eigenfunctions_[i][j][k]<<",";
			}
			std::cout<<endl;
		}
		std::cout<<"-----------------"<<endl;
	}
    return true;
}

void RealTimeExampleBasedDeformer::projectOnEigenFunctions(const VolumetricMesh *mesh, const double *displacement, const double *vertex_volume,
                                                           const double **eigenfunctions, const double *eigenvalues, unsigned int eigenfunction_num,
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
			std::cout<<mesh->getVertex(vert_idx)<<std::endl;
			std::cout<<mesh->getVertex(vert_idx)[vert_idx]<<std::endl;
			for(unsigned int dim=0;dim<3;++dim)
			{
				vert_pos[dim]=(*mesh->getVertex(vert_idx))[dim]+displacement[3*vert_idx+dim];
				eigencoefs[i][dim]+=vert_pos[dim]*eigenfunctions[i][vert_idx]*vertex_volume[vert_idx]*scale_factor;
			}
		}
	}
}

void RealTimeExampleBasedDeformer::reconstructFromEigenCoefs(const double **eigenfunctions, const double *eigenvalues,const Vec3d *eigencoefs,
                                                             int eigenfunction_num, int vert_num, Vec3d *vert_pos)
{
    double scale_factor;
	for(unsigned int vert_idx=0;vert_idx<vert_num;++vert_idx)
	{
		vert_pos[vert_idx][0]=vert_pos[vert_idx][1]=vert_pos[vert_idx][2]=0.0;
		for(unsigned int i=0;i<eigenfunction_num;++i)
		{
			scale_factor=1.0/eigenvalues[i];
			for(unsigned int dim=0;dim<3;++dim)
			{
				vert_pos[vert_idx][dim]+=eigencoefs[vert_idx][dim]*eigenfunctions[i][vert_idx]*scale_factor;
			}
		}
	}
}

}  //namespace RTLB
