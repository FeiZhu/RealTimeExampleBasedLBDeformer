/*
 * @file: ReducedStVKCubatureForceModel.h
 * @brief: reduced force model for stvk material using cubature optimization
 * @author: Mirror
 *
 */

#include <iomanip>
#include <iostream>
#include "reducedStVKCubatureForceModel.h"
#include "volumetricMeshENuMaterial.h"
#include "matrix.h"
#include "matrixProjection.h"
#include "performanceCounter.h"

ReducedStVKCubatureForceModel::ReducedStVKCubatureForceModel(const int &r,VolumetricMesh *volumetricMesh,double *U,
                            const int &cubica_num, const double *cubica_weights,const unsigned int *cubica_elements,
                            double **restpos,bool addGravity, double g)
{
    r_=r;
    add_gravity_=addGravity;
    g_=g;
    cubica_num_=cubica_num;
    cubica_weights_ = new double[cubica_num_];
    cubica_elements_ = new double[cubica_num_];
    // std::cout<<"b\n";
    for(int i=0;i<cubica_num_;++i)
    {
        cubica_weights_[i]=cubica_weights[i];
        cubica_elements_[i]=cubica_elements[i];
    }
    // std::cout<<"c\n";
    volumetric_mesh_ = volumetricMesh;
    U_ = new double*[3*volumetricMesh->getNumVertices()];
    for(int i=0;i<3*volumetricMesh->getNumVertices();++i)
    {
        U_[i] = new double[r_];
        for(int j=0;j<r_;++j)
            U_[i][j]=0.0;
    }
    // std::cout<<"d\n";
    int num=volumetricMesh->getNumVertices();
    //change *U to **U_,U_[vert_idx][basis_idx]
// Matrix<double> Umatrix((int)(3*num),(int)r_);
// Matrix<double> M((int)(3*num),(int)(3*num));
    for(int j=0;j<r_;++j)
        for(int i=0;i<3*num;++i)
        {
            U_[i][j]=U[3*num*j+i];
            // Umatrix(i,j)=U_[i][j];
        }
    // for(int i=0;i<3*num;++i)
    //     for(int j=0;j<3*num;++j)
    //     {
    //         if(i==j)
    //             M(i,j)=1.0/8144;
    //         else
    //             M(i,j)=0.0;
    //     }
    //
    //         Matrix<double> temp=Umatrix.MultiplyT(M);
    //         Matrix<double> temp1=temp*Umatrix;
    //         for(int i=0;i<r_;++i)
    //             {
    //                 for(int j=0;j<r_;++j)
    //                 {
    //                     std::cout<<temp1(i,j)<<",";
    //                 }
    //                 std::cout<<"\n";
    //             }
    //             getchar();
    //test U is normalized or not
    // double test=0.0;
    // for(int i=0;i<3*num;++i)
    // {
    //     test+=U_[1][i]*U_[2][i]*(1.0/8144);
    // }
    // std::cout<<test<<"\n";
    restpos_ = new double*[cubica_num_];
    for(int i=0;i<cubica_num_;++i)
    {
        restpos_[i]=new double[12];
        for(int j=0;j<12;++j)
            restpos_[i][j]=restpos[i][j];
    }
    //init cubica subbasis
    cubica_subBasis_=new double**[(int)cubica_num_];
    for(int i=0;i<cubica_num_;++i)
    {
        cubica_subBasis_[i]=new double*[12];
        for(int j=0;j<12;++j)
        {
            cubica_subBasis_[i][j]=new double[(int)r_];
        }
    }
    for(int i=0;i<cubica_num_;++i)
    {
        int ele=cubica_elements_[i];
        for(int j=0;j<4;++j)
        {
            int vertID=volumetric_mesh_->getVertexIndex(ele,j);
            for(int k=0;k<r_;++k)
            {
                cubica_subBasis_[i][3*j][k]=U_[3*vertID][k];
                cubica_subBasis_[i][3*j+1][k]=U_[3*vertID+1][k];
                cubica_subBasis_[i][3*j+2][k]=U_[3*vertID+2][k];
            }
        }
    }
    deformed_=new double[12];
	memset(deformed_,0.0,sizeof(double)*12);
    // modal_matrix_ = new ModalMatrix(volumetric_mesh_->getNumVertices(),r_,U);
    VolumetricMesh::Material * material = volumetric_mesh_->getElementMaterial(0);
	VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
	if (eNuMaterial == NULL)
	{
		std::cout<<"Error: mesh does not consist of E, nu materials.\n";
		exit(0);
	}
	lamda_ = eNuMaterial->getLambda();
	mu_ = eNuMaterial->getMu();
    gf_=new double[(int)r];
    // std::cout<<"f\n";
    InitGravity(U);
}

ReducedStVKCubatureForceModel::~ReducedStVKCubatureForceModel()
{
    if(cubica_weights_)
        delete[] cubica_weights_;
    if(cubica_elements_)
        delete[] cubica_elements_;
    for(int i=0;i<cubica_num_;++i)
        delete[] restpos_[i];
    delete[] restpos_;
    for(int i=0;i<3*volumetric_mesh_->getNumVertices();++i)
        delete[] U_[i];
    delete[] U_;
    if(gravity_force_)
        delete[] gravity_force_;
    if(volumetric_mesh_)
        delete[] volumetric_mesh_;
    if(gf_)
        delete[] gf_;
    for(int i=0;i<cubica_num_;++i)
    {
        if(cubica_subBasis_[i])
        {
            for(int j=0;j<12;++j)
            {
                if(cubica_subBasis_[i][j])
                    delete[] cubica_subBasis_[i][j];
            }
            delete[] cubica_subBasis_[i];
        }
    }
}

void ReducedStVKCubatureForceModel::GetInternalForce(double * q, double * internalForces)
{
    computeReducedInternalForce(q,internalForces);
}

void ReducedStVKCubatureForceModel::GetTangentStiffnessMatrix(double * q, double * tangentStiffnessMatrix)
{
    computeReducedStiffnessMatrix(q,tangentStiffnessMatrix);
}
void ReducedStVKCubatureForceModel::InitGravity(double *U)
{
    int n=volumetric_mesh_->getNumVertices();
    double *full_gravity_force=new double[3*n];
    memset(full_gravity_force,0.0,sizeof(double)*3*n);
    double *unit_reduced_gravity_force=new double[r_];
    memset(unit_reduced_gravity_force,0.0,sizeof(double)*r_);
    gravity_force_=new double[r_];
    memset(gravity_force_,0.0,sizeof(double)*r_);
    volumetric_mesh_->computeGravity(full_gravity_force,1.0);
    ProjectVector(3*n,r_,U,unit_reduced_gravity_force,full_gravity_force);
    for(int i=0;i<r_;++i)
        gravity_force_[i] = g_*unit_reduced_gravity_force[i];
    for(int i=0;i<r_;++i)
    std::cout<<gravity_force_[i]<<",";
    delete[] full_gravity_force;
    delete[] unit_reduced_gravity_force;
}

Matrix<double> ReducedStVKCubatureForceModel::tetSubBasis(const int &cubica_idx) const
{//tet_subBasis is 12*r
	Matrix<double> tet_subBasis(12,(int)r_,true);
	// int *global_idx=new int[4];
	// for(int i=0;i<4;++i)
	// 	global_idx[i]=volumetric_mesh_->getVertexIndex(ele,i);
    // std::cout<<"tetBasis:\n";
	for(int j=0;j<4;++j)
	{
        // std::cout<<"vertex:"<<global_idx[j]*3<<",\n";
        // int vertID=volumetric_mesh_->getVertexIndex(ele,j);
		for(int i=0;i<r_;++i)
		{
			tet_subBasis(3*j,i)=cubica_subBasis_[cubica_idx][3*j][i];
			tet_subBasis(3*j+1,i)=cubica_subBasis_[cubica_idx][3*j+1][i];
			tet_subBasis(3*j+2,i)=cubica_subBasis_[cubica_idx][3*j+2][i];

			// tet_subBasis(3*j,i)=U_[3*global_idx[j]][i];
			// tet_subBasis(3*j+1,i)=U_[3*global_idx[j]+1][i];
			// tet_subBasis(3*j+2,i)=U_[3*global_idx[j]+2][i];
            // std::cout<<tet_subBasis(3*j,i)<<","<<tet_subBasis(3*j+1,i)<<","<<tet_subBasis(3*j+2,i)<<",\n";
		}
	}
    // getchar();
	// delete[] global_idx;
	return tet_subBasis;
}
Mat3d ReducedStVKCubatureForceModel::computeDs(const double *reduced_pos) const
{//reduced dis is 12*1 for each element
	Mat3d Ds(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
	{
		Ds[0][i]=reduced_pos[3*i]-reduced_pos[9];
		Ds[1][i]=reduced_pos[3*i+1]-reduced_pos[10];
		Ds[2][i]=reduced_pos[3*i+2]-reduced_pos[11];

	}
	return Ds;
}
Mat3d ReducedStVKCubatureForceModel::computeDmInv(const int &ele) const
{//3x3
	int *global_idx=new int[4];
	Vec3d *vert_pos=new Vec3d[4];
	for(int j=0;j<4;++j)
	{
		global_idx[j]=volumetric_mesh_->getVertexIndex(ele,j);
		vert_pos[j]=*volumetric_mesh_->getVertex(global_idx[j]);
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
Mat3d ReducedStVKCubatureForceModel::computeF(const int &cubica_idx,const double *q) const//q:rx1
{//for all cubica elements
    int ele=cubica_elements_[cubica_idx];
    double *temp=new double[12];
    memset(temp,0.0,sizeof(double)*12);
    for(int i=0;i<12;++i)
        for(int j=0;j<r_;++j)
            temp[i]+=cubica_subBasis_[cubica_idx][i][j]*q[j];

	for(int i=0;i<12;++i)
		deformed_[i]=restpos_[cubica_idx][i]+temp[i];
	Mat3d F=computeDs(deformed_)*computeDmInv(ele);//F for each ele
    //temp handle invert F_
	// Mat3d U,V;
    // Vec3d Fhat;
    // ModifiedSVD(F,U,Fhat,V);
    // //clamphat if below the principle stretch threshold
    // double principle_threshold = 0.6;
    // for(unsigned int i = 0; i < 3 ; ++i)
    //     if(Fhat[i] < principle_threshold)
    //         Fhat[i] = principle_threshold;
    // Mat3d Fhat_mat;
    // for(unsigned int i = 0; i < 3; ++i)
    //     for(unsigned int j = 0; j < 3; ++j)
    //         Fhat_mat[i][j] = (i==j)?Fhat[i]:0;
    // F = U*Fhat_mat*trans(V);
    delete[] temp;
    return F;
}
Mat3d ReducedStVKCubatureForceModel::computeF_gradient(const int &ele,const int &vert_idx,const int &vert_idx_dim) const
{
	//vert_idx denotes the j-th vertex, vert_idx_dim is the k-th coordinate of the j-th vertex
	//we get the result as dF_i/dx_j^k, which is the derivative force of vertex i to the vertex j on the coordinate k
	//identity vector definition
	if((vert_idx>4)||(vert_idx_dim>3))
	{
		std::cout<<"the vert_idx or the vert_idx_dim is out of range, they should be smaller than 3";
	}
	Mat3d result_matrix(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	std::vector<Vec3d> e_vector;
	e_vector.resize(3,Vec3d(0.0,0.0,0.0));
	e_vector[0][0]=1.0;
	e_vector[1][1]=1.0;
	e_vector[2][2]=1.0;
	Mat3d e_matrix(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    Mat3d DmInv=computeDmInv(ele);
    //Mat3d DmInv(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0);
	//for j=1,2,3 we have dF/dx_j^k=e_k*e_j^T*Dm^-1
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
		Mat3d vert_cord_matrix(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
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
		result_matrix=vert_cord_matrix*DmInv;
	}
	else
	{
		//compute dF/dx_j^k,j=1,2,3;k=1,2,3
		result_matrix=e_matrix*DmInv;
	}
	return result_matrix;
}

Mat3d ReducedStVKCubatureForceModel::firstPiolaKirchhoff(Mat3d &F) const
{//3*3
	Mat3d P(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    Mat3d I(1.0);
    Mat3d temp=trans(F)*F;
	Mat3d E=0.5*(temp-I);
    double trace_E=E[0][0]+E[1][1]+E[2][2];
    P=F*(2*mu_*E+lamda_*trace_E*I);
	return P;
}
Mat3d ReducedStVKCubatureForceModel::computeP_gradient(const int &ele,const Mat3d &F,const int &vert_idx,
														const int &vert_idx_dim) const
{
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
void ReducedStVKCubatureForceModel::computeReducedDis(const double *x, double *q) const
{
	memset(q,0.0,sizeof(double)*r_);
    for(unsigned int i=0;i<r_;++i)
    {
        for(unsigned int j=0;j<3*volumetric_mesh_->getNumVertices();++j)
        {
            q[i]+=x[j]*U_[j][i];
        }
    }
}
void ReducedStVKCubatureForceModel::computeReducedEnergy(const double *q,double &energy) const
{
	energy=0.0;
	//computeF(reduced_dis);//compute F for all cubica elements
	for(int cubica_idx=0;cubica_idx<cubica_num_;++cubica_idx)
	{
		Mat3d F=computeF(cubica_idx,q);
        Mat3d I(1.0);
        Mat3d temp=trans(F)*F;
    	Mat3d E=0.5*(temp-I);
        double trace_E=E[0][0]+E[1][1]+E[2][2];
        double doubleE=0.0;
        for(int i=0;i<3;++i)
            for(int j=0;j<3;++j)
                doubleE+=E[i][j]*E[i][j];
        double element_energy=0.5*lamda_*trace_E*trace_E+mu_*doubleE;
		energy += cubica_weights_[cubica_idx]*element_energy;
	}
}
void ReducedStVKCubatureForceModel::computeReducedInternalForce(const double *q,double *forces) const
{//q:r*1
	// PerformanceCounter counter2;
	// counter2.StartCounter();
	memset(forces,0.0,sizeof(double)*r_);
    // double total_time=0.0,other_total_time=0.0,F_time=0.0,assemble_f=0.0;
	for(int cubica_idx=0;cubica_idx<cubica_num_;++cubica_idx)
	{
		int ele=cubica_elements_[cubica_idx];
    	Mat3d F=computeF(cubica_idx,q);
		Mat3d P=firstPiolaKirchhoff(F);
        Mat3d temp1=trans(computeDmInv(ele));
        Mat3d temp=P*temp1;
        double *ele_force=new double[12];
        memset(ele_force,0.0,sizeof(double)*12);
		for(int i=0;i<4;++i)
		{
			for(int j=0;j<3;++j)
			{
				if(i==3)
					ele_force[3*i+j]=(-1.0)*(temp[j][0]+temp[j][1]+temp[j][2]);
				else
					ele_force[3*i+j]=temp[j][i];
			}
		}
        memset(gf_,0.0,sizeof(double)*r_);
        for(int i=0;i<r_;++i)
            for(int j=0;j<12;++j)
                gf_[i]+=cubica_subBasis_[cubica_idx][j][i]*ele_force[j];
        delete[] ele_force;
		for(int i=0;i<r_;++i)
			forces[i] += cubica_weights_[cubica_idx]*gf_[i];
	}
    if(add_gravity_)
        for(int i=0;i<r_;++i)
            forces[i] -= gravity_force_[i];
    // counter2.StopCounter();
    // std::cout<<"integrator compute internal force:"<<counter2.GetElapsedTime()<<"\n";
}

void ReducedStVKCubatureForceModel::computeReducedStiffnessMatrix(const double *q,double *reduced_K/*Matrix<double> &reduced_K*/) const
{
    double total_time=0.0,else_time=0.0,other_count_time=0.0,assemble_k=0.0,F_time=0.0,whole_time=0.0;
    static int count=1;
    // ++count;
    Matrix<double> K((int)r_,(int)r_);
    // PerformanceCounter counter1;
    // counter1.StartCounter();
    // std::cout<<cubica_num_<<"\n";
	for(int cubica_idx=0;cubica_idx<cubica_num_;++cubica_idx)
	{

    	// PerformanceCounter counter21;
    	// counter21.StartCounter();
        //
    	// PerformanceCounter counter33;
    	// counter33.StartCounter();
		int ele=cubica_elements_[cubica_idx];
		Matrix<double> subU=tetSubBasis(cubica_idx);//12*r

        // counter33.StopCounter();
        // else_time+=counter33.GetElapsedTime();
		Matrix<double> ele_K(12,12);
    	// PerformanceCounter counter2;
    	// counter2.StartCounter();
        Mat3d F=computeF(cubica_idx,q);
        // counter2.StopCounter();
        // F_time+=counter2.GetElapsedTime();

    	// PerformanceCounter counter3;
    	// counter3.StartCounter();
		Mat3d trans_DmInv=trans(computeDmInv(ele));
        // counter3.StopCounter();
        // std::cout<<"computeK-counter3:"<<counter3.GetElapsedTime()<<"\n";
        // PerformanceCounter counter4;
    	// counter4.StartCounter();
		for(int i=0;i<4;++i)
		{
			std::vector<Mat3d> g_derivative(3);//computes dg/dx_j^0,dg/dx_j^1,dg/dx_j^2
			g_derivative.clear();
			for(int j=0;j<3;++j)
			{
                g_derivative[j].set(0.0);
				g_derivative[j]=computeP_gradient(ele,F,i,j)*trans_DmInv;
			}
            for(int j=0;j<4;++j)
			{
				// Mat3d f_derivative(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
				for(int row=0;row<3;++row)
				{
					for(int col=0;col<3;++col)
					{
						if(j==3)
                        {
                            ele_K(3*i+row,3*j+col)=(-1.0)*(g_derivative[row][col][0]+g_derivative[row][col][1]+g_derivative[row][col][2]);
                        }
						else
							ele_K(3*i+row,3*j+col)=g_derivative[row][col][j];
					}
				}
			}
		}
        // counter4.StopCounter();
        // assemble_k+=counter4.GetElapsedTime();
        //
        // counter21.StopCounter();
        // other_count_time+=counter21.GetElapsedTime();
        // PerformanceCounter counter5;
    	// counter5.StartCounter();
        // MultiplyMatrix()
        // modal_matrix_ = new ModalMatrix(volumetric_mesh_->getNumVertices(),r_,U);
        // modal_matrix_->ProjectMatrix(r_,ele_K,)
        Matrix<double> temp=subU.MultiplyT(ele_K);

        // for(int i=0;i<r_;++i)
        // {
        //     for(int j=0;j<12;++j)
        //         std::cout<<temp(i,j)<<",";
        //     std::cout<<"\n";
        // }
        // std::cout<<"\n";
        // std::cout<<"\n";
        // if(cubica_idx==0)
        //     getchar();
        Matrix<double> temp1=temp*subU;
        K+=(cubica_weights_[cubica_idx]*temp1);
        // counter5.StopCounter();
        // total_time+=counter5.GetElapsedTime();
	}
    for(int i=0;i<r_;++i)
        for(int j=0;j<r_;++j)
           reduced_K[i*r_+j]=K(i,j);
            // reduced_K(i,j)=K(i,j);
    // counter1.StopCounter();
    // whole_time+=counter1.GetElapsedTime();
    // std::cout<<"computeK--F:"<<F_time/count<<"\n";
    // std::cout<<"computeK--assemble K:"<<assemble_k/count<<"\n";
    // std::cout<<"computeK--else time:"<<else_time/count<<"\n";
    // std::cout<<"computeK--other time:"<<other_count_time/count<<"\n";    //
    // std::cout<<"computeK--multi matrix:"<<total_time/count<<"\n";
    // std::cout<<"computeK--for all cubica elements:"<<whole_time/count<<"\n";
}

Mat3d ReducedStVKCubatureForceModel::computeElasticDs(const double *ele_deformed_pos) const
{//reduced dis is 12*1 for each element
	Mat3d Ds(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
	{
		// Ds[0][i]=ele_deformed_pos[3*i]-ele_deformed_pos[9];
		// Ds[1][i]=ele_deformed_pos[3*i+1]-ele_deformed_pos[10];
		// Ds[2][i]=ele_deformed_pos[3*i+2]-ele_deformed_pos[11];
        for(int j=0;j<3;++j)
            Ds[j][i]=ele_deformed_pos[3*i+j]-ele_deformed_pos[9+j];
	}
    // std::cout<<"ele_deformed_pos:\n";
    // for(int i=0;i<12;++i)
    // {
    //     std::cout<<ele_deformed_pos[i]<<"\n";
    // }
	return Ds;
}
Mat3d ReducedStVKCubatureForceModel::computeElasticDmInv(const double *reference_pos) const
{
    //reference_pos is 12x1 vector
	Mat3d Dm(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
	{
        for(int j=0;j<3;++j)
            Dm[j][i]=reference_pos[3*i+j]-reference_pos[9+j];
    	// Dm[i][1]=reference_pos[3+i]-reference_pos[9+i];
    	// Dm[i][2]=reference_pos[6+i]-reference_pos[9+i];
	}
	Mat3d DmInv=inv(Dm);
	return DmInv;
}
Mat3d ReducedStVKCubatureForceModel::computeReducedElasticF(const double *ele_dis,const double *ele_reference_pos) const
{
	//compute displacement x=X+Uq
	double *deformed=new double[12];
	memset(deformed,0.0,sizeof(double)*12);
    // Matrix<double> subU=tetSubBasis(ele);//12xr
    // double *temp=new double[12];
    // memset(temp,0.0,sizeof(double)*12);
    // for(int i=0;i<12;++i)
    //     for(int j=0;j<r_;++j)
    //         temp[i]+=cubica_subBasis_[cubica_idx][i][j]*q[j];//project reduced displacement q to full space

//     std::cout<<"q:\n";
//     for(int j=0;j<r_;++j)
//         std::cout<<q[j]<<",";
//     std::cout<<"\n";
//     std::cout<<"reference_pos:\n";
//     for(int j=0;j<3;++j)
// 	{
//         std::cout<<reference_pos[j]<<",";
//     }
// std::cout<<"\n";
	for(int j=0;j<12;++j)
	{
		deformed[j]=ele_reference_pos[j]+ele_dis[j];
        // std::cout<<deformed[j]<<",";
	}
    // std::cout<<"deformed:"<<deformed[3*vert+0]<<","<<reference_pos[3*vert+1]<<","<<reference_pos[3*vert+0]<<"\n";
	Mat3d F=computeElasticDs(deformed)*computeElasticDmInv(ele_reference_pos);
    // std::cout<<"ds:"<<computeElasticDs(deformed)<<"\n";
    // std::cout<<"DmInv:"<<computeElasticDmInv(ele_reference_pos)<<"\n";
    //temp handle invert F_
	Mat3d U,V;
    Vec3d Fhat;
    ModifiedSVD(F,U,Fhat,V);
    //clamphat if below the principle stretch threshold
    double principle_threshold = 0.6;
    for(unsigned int i = 0; i < 3 ; ++i)
        if(Fhat[i] < principle_threshold)
            Fhat[i] = principle_threshold;
    Mat3d Fhat_mat;
    for(unsigned int i = 0; i < 3; ++i)
        for(unsigned int j = 0; j < 3; ++j)
            Fhat_mat[i][j] = (i==j)?Fhat[i]:0;
    F = U*Fhat_mat*trans(V);
    // delete[] temp;
    // std::cout<<F<<"\n";
	delete[] deformed;
    return F;
}

void ReducedStVKCubatureForceModel::computeReducedElasticEnergy(const double *dis,double &energy,const double *reference_pos) const
{
	energy=0.0;
	//computeF(reduced_dis);//compute F for all cubica elements
	// for(int cubica_idx=0;cubica_idx<cubica_num_;++cubica_idx)
    for(int ele=0;ele<volumetric_mesh_->getNumElements();++ele)
	{
        // int ele=cubica_elements_[cubica_idx];
        double *ele_pos=new double[12];
        double *ele_dis=new double[12];
        for(int i=0;i<4;++i)
        {
            int vertID=volumetric_mesh_->getVertexIndex(ele,i);
            // std::cout<<"vertID:"<<vertID<<"\n";
            // getchar();
            for(int j=0;j<3;++j)
            {
                ele_pos[3*i+j]=reference_pos[3*vertID+j];
                ele_dis[3*i+j]=dis[3*vertID+j];
            }
        }
		Mat3d F=computeReducedElasticF(ele_dis,ele_pos);
        Mat3d I(1.0);
        Mat3d temp=trans(F)*F;
    	Mat3d E=0.5*(temp-I);
        double trace_E=E[0][0]+E[1][1]+E[2][2];
        double doubleE=0.0;
        for(int i=0;i<3;++i)
            for(int j=0;j<3;++j)
                doubleE+=E[i][j]*E[i][j];
        double element_energy=0.5*lamda_*trace_E*trace_E+mu_*doubleE;
		energy += element_energy;
	}
}
void ReducedStVKCubatureForceModel::computeReducedElasticInternalForce(const double *dis,double *forces,const double *reference_pos) const
{//q:r*1
    std::cout<<"compute reduced elastic internal force with reference pos, dis, to get a reduced force\n";
	// PerformanceCounter counter2;
	// counter2.StartCounter();
    // std::cout<<"begin:\n";
	// memset(forces,0.0,sizeof(double)*3*volumetric_mesh_->getNumVertices());
    memset(forces,0.0,sizeof(double)*r_);
    // double *global_dis=new double[3*volumetric_mesh_->getNumVertices()];
    // memset(global_dis,0.0,sizeof(double)*3*volumetric_mesh_->getNumVertices());
    // modal_matrix_->AssembleVector(dis,global_dis);
    // for(int i=0;i<r_;++i)
    // std::cout<<"initial-forces:"<<forces[i]<<"\n";
    // getchar();
    // std::cout<<3*volumetric_mesh_->getNumVertices()<<"\n";
    // getchar();
    // double *global_forces=new double[3*volumetric_mesh_->getNumVertices()];
    // memset(global_forces,0.0,sizeof(double)*3*volumetric_mesh_->getNumVertices());
    // double total_time=0.0,other_total_time=0.0,F_time=0.0,assemble_f=0.0;
    std::cout<<"cubica_num_:"<<cubica_num_<<"\n";
    for(int cubica_idx=0;cubica_idx<cubica_num_;++cubica_idx)
	{
		int ele=cubica_elements_[cubica_idx];
    // for(int ele=0;ele<volumetric_mesh_->getNumElements();++ele)
    // {
        // std::cout<<"ele:"<<ele<<"\n";
        double *ele_pos=new double[12];
        double *ele_dis=new double[12];
        for(int i=0;i<4;++i)
        {
            int vertID=volumetric_mesh_->getVertexIndex(ele,i);
            // std::cout<<"vertID:"<<vertID<<"\n";
            // getchar();
            for(int j=0;j<3;++j)
            {
                ele_pos[3*i+j]=reference_pos[3*vertID+j];
                ele_dis[3*i+j]=dis[3*vertID+j];
            }
        }
    	Mat3d F=computeReducedElasticF(ele_dis,ele_pos);
        // if(cubica_idx==0)
        // std::cout<<det(F)<<"\n";
		Mat3d P=firstPiolaKirchhoff(F);
        Mat3d temp1=trans(computeElasticDmInv(ele_pos));
        Mat3d temp=P*temp1;
		// Matrix<double> ele_force(12,1);
        double *ele_force=new double[12];
        memset(ele_force,0.0,sizeof(double)*12);
		for(int i=0;i<4;++i)
		{
            // int vertID=volumetric_mesh_->getVertexIndex(ele,i);
			for(int j=0;j<3;++j)
			{
				if(i==3)
					ele_force[3*i+j]=(-1.0)*(temp[j][0]+temp[j][1]+temp[j][2]);
				else
					ele_force[3*i+j]=temp[j][i];
                // std::cout<<"ele_force:"<<ele_force[3*i+j]<<",";
			}
		}
        // std::cout<<"---a---"<<ele<<"\n";
        // for(int i=0;i<4;++i)
        // {
        //     int vertID=volumetric_mesh_->getVertexIndex(ele,i);
        //     // std::cout<<"vertID:"<<vertID<<"\n";
        //     for(int j=0;j<3;++j)
        //         forces[3*vertID+j] += ele_force[3*i+j];
        // }
        // std::cout<<"---b---\n";
        double *g=new double[(int)r_];
        memset(g,0.0,sizeof(double)*r_);
        for(int i=0;i<r_;++i)
            for(int j=0;j<12;++j)
                g[i]+=cubica_subBasis_[cubica_idx][j][i]*ele_force[j];

		for(int i=0;i<r_;++i)
        {
            forces[i] += cubica_weights_[cubica_idx]*g[i];
        }

        delete[] ele_force;
        delete[] g;
        delete[] ele_pos;
        delete[] ele_dis;
	}
    // modal_matrix_->ProjectVector(global_forces,forces);
    // std::cout<<"----end\n";
    // delete[] global_forces;
}
void ReducedStVKCubatureForceModel::testEnergyGradients()
{
    // double *x=new double[r_];
    // double *dx=new double[r_];
    // double *grad = new double[r_];
    // srand((unsigned)time(0));
    // int lowest=1,highest=10;
	// int range=(highest-lowest)+1;
	// for(int i=0;i<r_;++i)
	// {
	// 	x[i]=(lowest+rand()%range)/10.0;
	// 	std::cout<<x[i]<<",";
	// }
    // computeReducedInternalForce(x,grad);
	// for(int i=0;i<r_;++i)
	// {
	// 	dx[i]=(lowest+rand()%range)/1.0e7;
	// 	x[i]+=dx[i];
	// }
    // double f_plus,f_min;
    // computeReducedEnergy(x,f_plus);
    // for(int i=0;i<r_;++i)
    // {
    //     x[i]-=2*dx[i];
    // }
    // computeReducedEnergy(x,f_min);
    // double f_grad_times_dx=0.0;
    // for(int i=0;i<r_;++i)
    //     f_grad_times_dx+=grad[i]*2*dx[i];
    //     std::cout<<f_plus<<"...\n";
    // std::cout<<"Objective, df, analytic: "<<std::setprecision(15)<<f_grad_times_dx<<", numerial: "<<std::setprecision(15)<<f_plus-f_min;
    // std::cout<<"f_plus:"<<f_plus<<",f_min:"<<f_min<<", absolute_error= "<<f_plus-f_min-f_grad_times_dx;
    // std::cout<<", rel_error= "<<(f_plus-f_min-f_grad_times_dx)/(fabs(f_grad_times_dx)>1e-20?fabs(f_grad_times_dx):1e-20)<<"\n";
    // delete[] x;
    // delete[] dx;
    // delete[] grad;
    double *x=new double[3*volumetric_mesh_->getNumVertices()];
    double *dx=new double[3*volumetric_mesh_->getNumVertices()];
    double *y=new double[3*volumetric_mesh_->getNumVertices()];
    double *grad = new double[3*volumetric_mesh_->getNumVertices()];
    srand((unsigned)time(0));
    int lowest=1,highest=10;
	int range=(highest-lowest)+1;
	for(int i=0;i<3*volumetric_mesh_->getNumVertices();++i)
	{
		y[i]=(lowest+rand()%range);
		// std::cout<<x[i]<<",";
	}
    // modal_matrix_->AssembleVector(x,qx);
	for(int i=0;i<3*volumetric_mesh_->getNumVertices();++i)
	{
		x[i]=(lowest+rand()%range)/10.0;
		// std::cout<<x[i]<<",";
	}
    computeReducedElasticInternalForce(x,grad,y);
	for(int i=0;i<3*volumetric_mesh_->getNumVertices();++i)
	{
		dx[i]=(lowest+rand()%range)/1.0e7;
		x[i]+=dx[i];
	}
    double f_plus,f_min;
    computeReducedElasticEnergy(x,f_plus,y);
    for(int i=0;i<3*volumetric_mesh_->getNumVertices();++i)
    {
        x[i]-=2*dx[i];
    }
    computeReducedElasticEnergy(x,f_min,y);
    double f_grad_times_dx=0.0;
    for(int i=0;i<3*volumetric_mesh_->getNumVertices();++i)
        f_grad_times_dx+=grad[i]*2*dx[i];
        std::cout<<f_plus<<"...\n";
    std::cout<<"Objective, df, analytic: "<<std::setprecision(15)<<f_grad_times_dx<<", numerial: "<<std::setprecision(15)<<f_plus-f_min;
    std::cout<<"f_plus:"<<f_plus<<",f_min:"<<f_min<<", absolute_error= "<<f_plus-f_min-f_grad_times_dx;
    std::cout<<", rel_error= "<<(f_plus-f_min-f_grad_times_dx)/(fabs(f_grad_times_dx)>1e-20?fabs(f_grad_times_dx):1e-20)<<"\n";
    delete[] x;
    delete[] dx;
    delete[] grad;
}
void ReducedStVKCubatureForceModel::testObjectiveGradients()
{
        //test internal force and force gradient--stiffness matrix
	// int num=r_;
	// double *q=new double[num];
	// //double *x=new double[num];
	// double *dx=new double[num];
	// double *f=new double[num];
	// double *f_plus=new double[num];
	// double *f_min=new double[num];
	// Matrix<double> K((int)num,(int)num);
	// double *dK=new double[num];
	// srand((unsigned)time(0));
    // double *x=new double[(int)num];
	// //randomly generate displacement in range [0.1,1]
	// // double x[24]={0.3,0.1,0.4,0.3,0.8,1,0.5,1,0.7,0.2,0.2,0.3,0.2,0.4,0.6,0.2,0.7,0.3,0.8,0.6,1,1,0.3,0.5};
	// for(int i=0;i<num;++i)
	// {
	// 	x[i]=(1+rand()%10)/10.0;
	// 	std::cout<<x[i]<<",";
	// }
    // std::cout<<"a\n";
	// computeReducedStiffnessMatrix(x,K);
    // // std::cout<<"KKKKKK:\n";
	// // for(int j=0;j<r_;++j)
	// // {
	// // 	for(int k=0;k<r_;++k)
	// // 	{
	// // 		std::cout<<K(j,k)<<",";
	// // 	}
	// // 	std::cout<<"\n";
	// // }
	// //perturb a little bit
	// for(int i=0;i<num;++i)
	// {
	// 	dx[i]=(1+rand()%10)/1.0e7;
	// 	x[i]+=dx[i];
    //     std::cout<<"x+:"<<x[i]<<",";
	// }
    // std::cout<<"b\n";
	// computeReducedInternalForce(x,f_plus);
    // std::cout<<"b-end\n";
	// for(int i=0;i<num;++i)
	// {
	// 	//dx[i]*=1.0;
	// 	x[i]-=2*dx[i];
	// }
	// computeReducedInternalForce(x,f_min);
	// Matrix<double> dis_matrix((int)num,1);
	// for(int i=0;i<num;++i)
	// 	dis_matrix(i,0)=2.0*dx[i];
	// Matrix<double> temp_matrix=K*dis_matrix;
	// for(int i=0;i<num;++i)
	// {
	// 	dK[i]=temp_matrix(i,0);
	// }
    // std::cout<<"c\n";
	// double max_rel_error=1000;
	// double df,dk;
	// for(int i=0;i<num;++i)
	// {
	// 	double rel_error = fabs((f_plus[i]-f_min[i]-dK[i])/dK[i]);
	// 	//std::cout<<"f_plus:"<<f_plus[i]<<"f_min:"<<f_min[i]<<"dk:"<<dK[i]<<"\n";
	// 	if(rel_error < max_rel_error)
	// 	{
	// 		df=f_plus[i]-f_min[i];
	// 		dk = dK[i];
	// 		max_rel_error = rel_error;
	// 	}
	// }
	// std::cout<<"Max relative error:\n";
	// std::cout<<"df: "<<std::setprecision(15)<<df<<"\n";
	// std::cout<<"k*dk: "<<std::setprecision(15)<<dk<<"\n";
	// std::cout<<"rel error: "<<max_rel_error<<"\n";
	// //delete[] x;
	// delete[] dx;
	// delete[] f;
    // delete[] dK;
	// delete[] f_plus;
	// delete[] f_min;
    // delete[] x;
}

int ReducedStVKCubatureForceModel::ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const
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

void ReducedStVKCubatureForceModel::FindOrthonormalVector(Vec3d & v, Vec3d & result) const
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
