/*
 * @file: neohookeanReducedForceModel.cpp
 * @brief: reduced force model for neohookean material
 * @author: Mirror
 *
 */

#include <iomanip>
#include <iostream>
#include "reducedNeoHookeanForceModel.h"
#include "volumetricMeshENuMaterial.h"
#include "matrix.h"
#include "matrixProjection.h"
#include "generateMassMatrix.h"
#include "performanceCounter.h"

ReducedNeoHookeanForceModel::ReducedNeoHookeanForceModel(const int &r,VolumetricMesh *volumetricMesh,double *U,
                            const int &cubica_num, const double *cubica_weights,const unsigned int *cubica_elements,
                            double **restpos,bool addGravity, double g)
{
    // std::cout<<"a\n";
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
    for(int j=0;j<r_;++j)
        for(int i=0;i<3*num;++i)
        {
            U_[i][j]=U[3*num*j+i];
        }
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

ReducedNeoHookeanForceModel::~ReducedNeoHookeanForceModel()
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
        for(int j=0;j<12;++j)
        {
            if(cubica_subBasis_[i][j])
                delete[] cubica_subBasis_[i][j];
        }
        if(cubica_subBasis_[i])
            delete[] cubica_subBasis_[i];
    }
}

void ReducedNeoHookeanForceModel::GetInternalForce(double * q, double * internalForces)
{
    computeReducedInternalForce(q,internalForces);
}

void ReducedNeoHookeanForceModel::GetTangentStiffnessMatrix(double * q, double * tangentStiffnessMatrix)
{
    computeReducedStiffnessMatrix(q,tangentStiffnessMatrix);
}
void ReducedNeoHookeanForceModel::InitGravity(double *U)
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

Matrix<double> ReducedNeoHookeanForceModel::tetSubBasis(const int &ele) const
{//tet_subBasis is 12*r
	Matrix<double> tet_subBasis(12,(int)r_,true);
	int *global_idx=new int[4];
	for(int i=0;i<4;++i)
		global_idx[i]=volumetric_mesh_->getVertexIndex(ele,i);
    // std::cout<<"tetBasis:\n";
	for(int j=0;j<4;++j)
	{
        // std::cout<<"vertex:"<<global_idx[j]*3<<",\n";
        int vertID=volumetric_mesh_->getVertexIndex(ele,j);
		for(int i=0;i<r_;++i)
		{
			tet_subBasis(3*j,i)=U_[3*global_idx[j]][i];
			tet_subBasis(3*j+1,i)=U_[3*global_idx[j]+1][i];
			tet_subBasis(3*j+2,i)=U_[3*global_idx[j]+2][i];
            // std::cout<<tet_subBasis(3*j,i)<<","<<tet_subBasis(3*j+1,i)<<","<<tet_subBasis(3*j+2,i)<<",\n";
		}
	}
    // getchar();
	delete[] global_idx;
	return tet_subBasis;
}
Mat3d ReducedNeoHookeanForceModel::computeDs(const double *reduced_dis) const
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
Mat3d ReducedNeoHookeanForceModel::computeDmInv(const int &ele) const
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
Mat3d ReducedNeoHookeanForceModel::computeF(const int &cubica_idx,const double *q) const//q:rx1
{//for all cubica elements
    int ele=cubica_elements_[cubica_idx];
    // Matrix<double> reduced_dis_matrix((int)r_,1);//rx1
	// for(int i=0;i<r_;++i)
	// {
	// 	reduced_dis_matrix(i,0)=q[i];
	// }
	//compute displacement x=X+Uq
	double *deformed=new double[12];
	memset(deformed,0.0,sizeof(double)*12);
    // Matrix<double> subU=tetSubBasis(ele);//12xr
    //Matrix<double> temp=subU*reduced_dis_matrix;//12xr,rx1->12x1
    double *temp=new double[12];
    memset(temp,0.0,sizeof(double)*12);
    for(int i=0;i<12;++i)
        for(int j=0;j<r_;++j)
            temp[i]+=cubica_subBasis_[cubica_idx][i][j]*q[j];

	for(int i=0;i<12;++i)
		deformed[i]=restpos_[cubica_idx][i]+temp[i];
	Mat3d F=computeDs(deformed)*computeDmInv(ele);//F for each ele

    delete[] temp;
	delete[] deformed;
    return F;
}
Mat3d ReducedNeoHookeanForceModel::computeF_gradient(const int &ele,const int &vert_idx,const int &vert_idx_dim) const
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

Mat3d ReducedNeoHookeanForceModel::firstPiolaKirchhoff(Mat3d &F) const
{//3*3
	Mat3d P(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	P=mu_*(F-trans(inv(F)))+lamda_*log(det(F))*trans(inv(F));
	return P;
}
Mat3d ReducedNeoHookeanForceModel::computeP_gradient(const int &ele,const Mat3d &F,const int &vert_idx,
														const int &vert_idx_dim) const
{
	Mat3d pFpx=computeF_gradient(ele,vert_idx,vert_idx_dim);
	Mat3d temp=inv(F)*pFpx;
	double trace=temp[0][0]+temp[1][1]+temp[2][2];
	Mat3d tiF=trans(inv(F));
	Mat3d pPpx=mu_*pFpx+(mu_-lamda_*log(det(F)))*tiF*trans(pFpx)*tiF+lamda_*trace*tiF;
	return pPpx;
}

void ReducedNeoHookeanForceModel::computeReducedEnergy(const double *q,double &energy) const
{
	energy=0.0;
	//computeF(reduced_dis);//compute F for all cubica elements
	for(int cubica_idx=0;cubica_idx<cubica_num_;++cubica_idx)
	{
	//	int ele=object_cubica_elements_[cubica_idx];
		Mat3d F=computeF(cubica_idx,q);
		Mat3d temp=trans(F)*F;
		double trace_c=temp[0][0]+temp[1][1]+temp[2][2];
		double lnJ=log(det(F));
		double element_energy=0.5*mu_*(trace_c-3)-mu_*lnJ+0.5*lamda_*lnJ*lnJ;
		energy += cubica_weights_[cubica_idx]*element_energy;
        // energy = element_energy;
	}

}
void ReducedNeoHookeanForceModel::computeReducedInternalForce(const double *q,double *forces) const
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
		//Matrix<double> ele_force(12,1);
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
		// Matrix<double> subU=tetSubBasis(ele);//12xr
        //Matrix<double> g=Transpose(subU)*ele_force;//rx1
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

void ReducedNeoHookeanForceModel::computeReducedStiffnessMatrix(const double *q,double *reduced_K/*Matrix<double> &reduced_K*/) const
{
    double total_time=0.0,other_count_time=0.0,assemble_k=0.0,F_time=0.0;
    Matrix<double> K((int)r_,(int)r_);
    // PerformanceCounter counter1;
    // counter1.StartCounter();
    // std::cout<<cubica_num_<<"\n";
	for(int cubica_idx=0;cubica_idx<cubica_num_;++cubica_idx)
	{

    	// PerformanceCounter counter21;
    	// counter21.StartCounter();
		int ele=cubica_elements_[cubica_idx];
		Matrix<double> subU=tetSubBasis(ele);//12*r
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
							ele_K(3*i+row,3*j+col)=(-1.0)*(g_derivative[row][col][0]+g_derivative[row][col][1]+g_derivative[row][col][2]);
						else
							ele_K(3*i+row,3*j+col)=g_derivative[row][col][j];
						// ele_K(3*i+row,3*j+col)=f_derivative[row][col];
					}
				}
			}
		}
        // counter4.StopCounter();
        // assemble_k+=counter4.GetElapsedTime();


        // counter21.StopCounter();
        // other_count_time+=counter21.GetElapsedTime();
        // PerformanceCounter counter5;
    	// counter5.StartCounter();
        // MultiplyMatrix()
        Matrix<double> temp=subU.MultiplyT(ele_K);
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
    // std::cout<<"computeK--F:"<<F_time<<"\n";
    // std::cout<<"computeK--assemble K:"<<assemble_k<<"\n";
    // std::cout<<"computeK--other time:"<<other_count_time<<"\n";    //
    // std::cout<<"computeK--multi matrix:"<<total_time<<"\n";
    // std::cout<<"computeK--for all cubica elements:"<<counter1.GetElapsedTime()<<"\n";
}
Mat3d ReducedNeoHookeanForceModel::computeElasticDmInv(const int &ele,const double *u) const
{//3x3int *global_idx=new int[4];
	Vec3d *vert_pos=new Vec3d[4];
	for(int j=0;j<4;++j)
	{
		int vertID=volumetric_mesh_->getVertexIndex(ele,j);
        for(int i=0;i<3;++i)
		    vert_pos[j][i]=(*volumetric_mesh_->getVertex(vertID))[i]+u[3*j+i];
	}
	Mat3d Dm(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	for(int i=0;i<3;++i)
	{
		Dm[i][0]=vert_pos[0][i]-vert_pos[3][i];
		Dm[i][1]=vert_pos[1][i]-vert_pos[3][i];
		Dm[i][2]=vert_pos[2][i]-vert_pos[3][i];
	}
	Mat3d DmInv=inv(Dm);
	delete[] vert_pos;
	return DmInv;
}
Mat3d ReducedNeoHookeanForceModel::computeReducedElasticF(const int &cubica_idx,const double *q,const double *u) const
{
    int ele=cubica_elements_[cubica_idx];
	//compute displacement x=X+Uq
	double *deformed=new double[12];
	memset(deformed,0.0,sizeof(double)*12);
    // Matrix<double> subU=tetSubBasis(ele);//12xr
    double *temp=new double[12];
    memset(temp,0.0,sizeof(double)*12);
    for(int i=0;i<12;++i)
        for(int j=0;j<r_;++j)
            temp[i]+=cubica_subBasis_[cubica_idx][i][j]*q[j];
	for(int j=0;j<12;++j)
	{
		deformed[j]=restpos_[cubica_idx][j]+u[j]+temp[j];
	}
	Mat3d F=computeDs(deformed)*computeElasticDmInv(ele,u);//F for each ele
    delete[] temp;
	delete[] deformed;
    return F;
}
void ReducedNeoHookeanForceModel::computeReducedElasticInternalForce(const double *q,double *forces,const double *u) const
{//q:r*1

	// PerformanceCounter counter2;
	// counter2.StartCounter();
	memset(forces,0.0,sizeof(double)*r_);
    double total_time=0.0,other_total_time=0.0,F_time=0.0,assemble_f=0.0;

    for(int cubica_idx=0;cubica_idx<cubica_num_;++cubica_idx)
	{
		int ele=cubica_elements_[cubica_idx];
        double *temp_u=new double[12];
        for(int i=0;i<4;++i)
        {
            int vertID=volumetric_mesh_->getVertexIndex(ele,i);
            for(int j=0;j<3;++j)
            {
                temp_u[3*i+j]=u[3*vertID+j];
            }
        }
    	Mat3d F=computeReducedElasticF(cubica_idx,q,temp_u);
		Mat3d P=firstPiolaKirchhoff(F);
        Mat3d temp1=trans(computeElasticDmInv(ele,temp_u));
        Mat3d temp=P*temp1;
		// Matrix<double> ele_force(12,1);
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
		// Matrix<double> subU=tetSubBasis(ele);//12xr
        // Matrix<double> g=Transpose(subU)*ele_force;//rx1
        double *g=new double[(int)r_];
        memset(g,0.0,sizeof(double)*r_);
        for(int i=0;i<r_;++i)
            for(int j=0;j<12;++j)
                g[i]+=cubica_subBasis_[cubica_idx][j][i]*ele_force[j];
		for(int i=0;i<r_;++i)
			forces[i] += cubica_weights_[cubica_idx]*g[i];
        delete[] ele_force;
        delete[] g;
        delete[] temp_u;
	}
}
void ReducedNeoHookeanForceModel::testEnergyGradients()
{
    double *x=new double[r_];
    double *dx=new double[r_];
    double *grad = new double[r_];
    srand((unsigned)time(0));
    int lowest=1,highest=10;
	int range=(highest-lowest)+1;
	for(int i=0;i<r_;++i)
	{
		x[i]=(lowest+rand()%range)/10.0;
		std::cout<<x[i]<<",";
	}
    computeReducedInternalForce(x,grad);
	for(int i=0;i<r_;++i)
	{
		dx[i]=(lowest+rand()%range)/1.0e7;
		x[i]+=dx[i];
	}
    double f_plus,f_min;
    computeReducedEnergy(x,f_plus);
    for(int i=0;i<r_;++i)
    {
        x[i]-=2*dx[i];
    }
    computeReducedEnergy(x,f_min);
    double f_grad_times_dx=0.0;
    for(int i=0;i<r_;++i)
        f_grad_times_dx+=grad[i]*2*dx[i];
        std::cout<<f_plus<<"...\n";
    std::cout<<"Objective, df, analytic: "<<std::setprecision(15)<<f_grad_times_dx<<", numerial: "<<std::setprecision(15)<<f_plus-f_min;
    std::cout<<"f_plus:"<<f_plus<<",f_min:"<<f_min<<", absolute_error= "<<f_plus-f_min-f_grad_times_dx;
    std::cout<<", rel_error= "<<(f_plus-f_min-f_grad_times_dx)/(fabs(f_grad_times_dx)>1e-20?fabs(f_grad_times_dx):1e-20)<<"\n";
    delete[] x;
    delete[] dx;
    delete[] grad;
}
void ReducedNeoHookeanForceModel::testObjectiveGradients()
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
