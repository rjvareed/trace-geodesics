#include <iostream>
#include <ginac/ginac.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

using namespace GiNaC;

int func(double t, const double y[], double f[], void *params){
	double a = *(double*)params;
	f[0] = y[1];
	f[1] = a*y[0];
	return GSL_SUCCESS;
}

int jac(double t,const double y[],double *dfdy, double dfdt[], void *params){
	double a = *(double*)params;
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,2,2);
	gsl_matrix *m = &dfdy_mat.matrix;
	
	//set jacobian to
	//( 0.0 1.0 )
	//( a   0.0 )
	gsl_matrix_set(m,0,0,0.0);
	gsl_matrix_set(m,0,1,1.0);
	gsl_matrix_set(m,1,0,a);
	gsl_matrix_set(m,1,1,0.0);
	
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	return GSL_SUCCESS;
}

int main(void){
	symbol x("x");
	symbol y("y");
	symbol X[2] = {x,y};
	symtab table;
	table["x"] = X[0];
	table["y"] = X[1];
	parser reader(table);
	ex g00 = reader("y^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/k))+1/(1-k/(x^2+y^2)^(1/2))");
	ex g01 = reader("x*y/((x^2+y^2)*((x^2+y^2)^(1/2)/k-1))");
	ex g10 = reader("x*y/((x^2+y^2)*((x^2+y^2)^(1/2)/k-1))");
	ex g11 = reader("x^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/k))+1/(1-k/(x^2+y^2)^(1/2))");
	matrix g_mat_covariant = {{g00,g01},{g10,g11}};
	matrix g_mat_contravariant = g_mat_covariant.inverse();
	ex g_covariant[2][2] = {{g_mat_covariant(0,0),g_mat_covariant(0,1)},{g_mat_covariant(1,0),g_mat_covariant(1,1)}};
	ex g_contravariant[2][2] = {{g_mat_contravariant(0,0),g_mat_contravariant(0,1)},{g_mat_contravariant(1,0),g_mat_contravariant(1,1)}};
	std::cout << g_covariant[0][0] <<"\n\n"<<g_covariant[0][1]<<"\n\n"<<g_covariant[1][0]<<"\n\n"<<g_covariant[1][1]<<"\n\n\n\n";
	std::cout << g_contravariant[0][0] <<"\n\n"<<g_contravariant[0][1]<<"\n\n"<<g_contravariant[1][0]<<"\n\n"<<g_contravariant[1][1]<<"\n\n";
	
	//christoffel symbols
	//G^a_bc=sum(g^ad(d_b g_dc + d_c g_db - d_d g_bc))
	ex christoffel_symbol[2][2][2];
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				for(int d=0;d<2;d++)
					christoffel_symbol[a][b][c] = 0.5*g_contravariant[a][d]*(g_covariant[d][c].diff(X[b])+g_covariant[d][b].diff(X[c])-g_covariant[b][c].diff(X[d]));
	
	symbol vx("vx");
	symbol vy("vy");
	symbol V[2] = {vx,vy};
	ex F_ex[4];
	
	// dx/ds = v_x
	F_ex[0] = vx;
	// dy/ds = v_y
	F_ex[1] = vy;
	// dv^a/ds = -Gamma^a_bc v^b v^c
	// dv_x/ds = -Gamma^0_ab v^a v^b
	F_ex[2] = 0.0;
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			F_ex[2]-=christoffel_symbol[0][a][b]*V[a]*V[b];
	// dv_y/ds = -Gamma^1_ab v^a v^b
	F_ex[3] = 0.0;
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			F_ex[3]-=christoffel_symbol[1][a][b]*V[a]*V[b];
	
	symbol Y[4] = {x,y,vx,vy};
	//set up jacobian
	//J_ij=df_i/dy_j
	ex jacobian_ex[4][4];
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			jacobian_ex[i][j] = F_ex[i].diff(Y[j]);
	
	std::cout << "\n\n";
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			std::cout << jacobian_ex[i][j] << "\n\n";
		}
	}

	//FUNCP_1P fp;
	//compile_ex(expression,x,fp);
	//std::cout << fp(2.0) <<"\n";
	
	
	
	/*
	double a = -1.0;
	gsl_odeiv2_system sys = {func    , jac,      2,         &a        };
	//                      {function, jacobian, dimension, parameters};
	
	gsl_odeiv2_driver *driver= gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rk8pd,1e-6,1e-6,0.0);
	
	int i=0;
	double t=0.0,t1=100.0;
	//initial conditions
	double y[2] = {1.0,0.0};
	
	for(i=1;i<=100;i++){
		double ti = i * t1 / 100.0;
		int status = gsl_odeiv2_driver_apply(driver,&t,ti,y);
		if(status != GSL_SUCCESS){
			printf("error, return value = %d\n",status);
			break;
		}

		printf("%lf %lf %lf\n",t,y[0],y[1]);
	}
	gsl_odeiv2_driver_free(driver);
	*/
	
	return 0;
}
