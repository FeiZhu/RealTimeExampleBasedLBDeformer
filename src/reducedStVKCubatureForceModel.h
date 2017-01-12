/*
 * @file: ReducedStVKCubatureForceModel.h
 * @brief: reduced force model for stvk material using cubature optimization
 * @author: Mirror
 *
 */

#ifndef _REDUCEDSTVKCUBATUREFORCE_MODEL_H_
#define _REDUCEDSTVKCUBATUREFORCE_MODEL_H_

#include <stdlib.h>
#include "reducedForceModel.h"
#include "matrix.h"
#include "volumetricMesh.h"
#include "tetMesh.h"
#include "mat3d.h"
#include "modalMatrix.h"
class VolumetricMesh;

class ReducedStVKCubatureForceModel : public virtual ReducedForceModel
{
public:
    ReducedStVKCubatureForceModel(const int &r,VolumetricMesh *volumetricMesh,double *U,
                                const int &cubica_num, const double *cubica_weights,const unsigned int *cubica_elements,
                                double **restpos,bool addGravity=false, double g=9.81);
    virtual ~ReducedStVKCubatureForceModel();
    virtual void GetInternalForce(double * q, double * internalForces);
    virtual void GetTangentStiffnessMatrix(double * q, double * tangentStiffnessMatrix);

    //virtual void GetReducedInternalForceClass() { return f_; }
//    virtual void GetReducedStiffnessMatrixClass() { return K_; }
    void SetGravity(bool addGravity,double g,double *U) {this->add_gravity_ = addGravity;}
    void testEnergyGradients();
    void testObjectiveGradients();
    void computeReducedEnergy(const double *q,double &energy) const;
    void computeReducedInternalForce(const double *q,double *forces) const;
    //compute elastic pos for given reference_configration, q is reduced displacement of (current-reference)
    void computeReducedElasticInternalForce(const double *dis,double *forces,const double *reference_pos) const;
protected:
    void InitGravity(double *U); //aux function
    // Matrix<double> vertexSubBasis(const int &vert_idx) const;//3*r
    Matrix<double> tetSubBasis(const int &ele) const;//12*r
    Mat3d computeDs(const double *reduced_pos) const;
    Mat3d computeDmInv(const int &ele) const;
    Mat3d computeF(const int &cubica_idx,const double *q) const;
    Mat3d firstPiolaKirchhoff(Mat3d &F) const;
    Mat3d computeF_gradient(const int &ele,const int &vert_idx,const int &vert_idx_dim) const;
    Mat3d computeP_gradient(const int &ele,const Mat3d &F,const int &vert_idx,const int &vert_idx_dim) const;
    void computeReducedDis(const double *x, double *q) const;

    void computeReducedStiffnessMatrix(const double *q,double *reduced_K/*Matrix<double> &reduced_K*/) const;

    //used for example-based elastic energy
    Mat3d computeElasticDs(const double *ele_deformed_pos) const;
    Mat3d computeElasticDmInv(const double *ele_reference_pos) const;
    Mat3d computeReducedElasticF(const double *ele_dis,const double *ele_reference_pos) const;
    void computeReducedElasticEnergy(const double *dis,double &energy,const double *reference_pos) const;
    void FindOrthonormalVector(Vec3d & v, Vec3d & result) const;
    int ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const;


protected:
    int r_;
    VolumetricMesh *volumetric_mesh_=NULL;
    TetMesh *tet_mesh_=NULL;
    double **U_=NULL;
    int cubica_num_;
    double *cubica_elements_=NULL;
    double *cubica_weights_=NULL;
    double ***cubica_subBasis_ = NULL;
    bool add_gravity_;
    double g_;
    double *gravity_force_=NULL;
    // double *f_=NULL;
    // double *K_=NULL;
    double **restpos_=NULL;
    double *deformed_ = NULL;
    //material
    double mu_=0.0;
    double lamda_=0.0;
    bool own_neoHookean_stiffness_matrix_=false;
    double *gf_=NULL;
    double *test_=NULL;
    ModalMatrix *modal_matrix_ = NULL;
};

#endif
