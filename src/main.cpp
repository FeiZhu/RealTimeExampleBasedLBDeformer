/*
 * @file: main.cpp
 * @brief: entry of the program
 * @author: Fei Zhu,Mirror
 *
 */

 #include "real_time_example_based_deformer.h"

 int main()
 {
 	const std::string file_name="models/a.veg";
 	std::string write_file_name="models/b.smesh";
 	std::string file_name_prefix="models/armadillo";
 	std::string write_example_file_name_prefix="models/ex_armadillo";
 	std::string surface_mesh_file_name="models/a.obj";
 	std::string write_obj_file_name="models/b.obj";
 	std::string object_eigen_file_name="models/a.eigen";
 	std::string write_object_eigen_file_name="models/b.eigen";
 	std::string example_eigen_file_name_predix="models/ex_a";
 	std::string write_example_eigen_file_name_predix="models/ex_b";
 	std::string reduced_file_name="models/reduced.basis";
    std::string corresponding_file_name="models/cor.ini";
 	RealTimeExampleBasedDeformer real_time_example_based_deformer;
 	//real_time_example_based_deformer.loadSimulationMesh(file_name);
 	//real_time_example_based_deformer.saveSimulationMesh(write_file_name);
 	//real_time_example_based_deformer.loadExamples(file_name_prefix,2);
 	//real_time_example_based_deformer.saveExamples(write_example_file_name_prefix);
 	//real_time_example_based_deformer.loadVisualMesh(surface_mesh_file_name);
 	//real_time_example_based_deformer.saveVisualMesh(write_obj_file_name);
 	//real_time_example_based_deformer.loadObjectEigenfunctions(object_eigen_file_name);
 	//real_time_example_based_deformer.saveObjectEigenfunctions(write_object_eigen_file_name);
 	//real_time_example_based_deformer.loadExampleEigenFunctions(example_eigen_file_name_predix);
 	//real_time_example_based_deformer.saveExampleEigenfunctions(write_example_eigen_file_name_predix);
 	//real_time_example_based_deformer.loadReducedBasis(reduced_file_name);
    //corresponding functions test:
    //first step
    real_time_example_based_deformer.loadCorrespondenceData(corresponding_file_name);
    //second step:generate new eigenfunction for examples
    real_time_example_based_deformer.loadSimulationMesh(file_name);
    real_time_example_based_deformer.loadExamples(file_name_prefix);
    real_time_example_based_deformer.loadObjectEigenfunctions(object_eigen_file_name);
    real_time_example_based_deformer.loadExampleEigenFunctions(example_eigen_file_name_predix);
    real_time_example_based_deformer.registerEigenfunctions();
    //third step:save new eigenfunctions for examples_
    real_time_example_based_deformer.saveExampleEigenfunctions(write_example_eigen_file_name_predix);
 	return 1;
 }
