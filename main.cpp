#include <iostream>
#include <ginac/ginac.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "DrawingArea.h"

using namespace GiNaC;

FUNCP_2P christoffel_symbol[2][2][2];
FUNCP_2P christoffel_symbol_derivative[2][2][2][2];

int func(double s, const double y[], double f[], void *params){
	// dx/ds=v_x
	f[0] = y[2];
	// dy/ds=v_y
	f[1] = y[3];
	//dv_x/ds=-G^0_ab v^a v^b
	f[2] = 0.0;
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			f[2]-=christoffel_symbol[0][a][b](y[0],y[1])*y[2+a]*y[2+b];
	//dv_y/ds=-G^1_ab v^a v^b
	f[3] = 0.0;
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			f[3]-=christoffel_symbol[1][a][b](y[0],y[1])*y[2+a]*y[2+b];
	return GSL_SUCCESS;
}

int jac(double t,const double y[],double *dfdy, double dfds[], void *params){
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,4,4);
	gsl_matrix *m = &dfdy_mat.matrix;
	
	//set jacobian to
	//( 0 0 1 0 )
	//( 0 0 0 1 )
	//( df2/dy0 df2/dy1 df2/dy2 df2/dy3 )
	//( df3/dy0 df3/dy1 df3/dy2 df3/dy3 )
	gsl_matrix_set(m,0,0,0.0);
	gsl_matrix_set(m,0,1,0.0);
	gsl_matrix_set(m,0,2,1.0);
	gsl_matrix_set(m,0,3,0.0);
	
	gsl_matrix_set(m,1,0,0.0);
	gsl_matrix_set(m,1,1,0.0);
	gsl_matrix_set(m,1,2,0.0);
	gsl_matrix_set(m,1,3,1.0);
	
	double df2dy0 = 0;
	double df2dy1 = 0;
	double df3dy0 = 0;
	double df3dy1 = 0;
	
	//derivatives for d_i (-G^a_bc v^b v^c)
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++){
			df2dy0-=y[2+a]*y[2+b]*christoffel_symbol_derivative[0][0][a][b] (y[0],y[1]);
			df2dy1-=y[2+a]*y[2+b]*christoffel_symbol_derivative[1][0][a][b] (y[0],y[1]);
			df3dy0-=y[2+a]*y[2+b]*christoffel_symbol_derivative[0][1][a][b] (y[0],y[1]);
			df3dy1-=y[2+a]*y[2+b]*christoffel_symbol_derivative[1][1][a][b] (y[0],y[1]);
		}

	//manually calculated derivatives for d/dv^i (G^a_bc v^a v^b)
	double df2dy2 = -2*christoffel_symbol[0][0][0] (y[0],y[1])*y[2]-y[3]*(christoffel_symbol[0][0][1] (y[0],y[1])+ christoffel_symbol[0][1][0] (y[0],y[1]));
	double df2dy3 = -2*christoffel_symbol[0][1][1] (y[0],y[1])*y[3]-y[2]*(christoffel_symbol[0][0][1] (y[0],y[1])+ christoffel_symbol[0][1][0] (y[0],y[1]));
	double df3dy2 = -2*christoffel_symbol[1][0][0] (y[0],y[1])*y[2]-y[3]*(christoffel_symbol[1][0][1] (y[0],y[1])+ christoffel_symbol[1][1][0] (y[0],y[1]));
	double df3dy3 = -2*christoffel_symbol[1][1][1] (y[0],y[1])*y[3]-y[2]*(christoffel_symbol[1][0][1] (y[0],y[1])+ christoffel_symbol[1][1][0] (y[0],y[1]));
	
	gsl_matrix_set(m,2,0,df2dy0);
	gsl_matrix_set(m,2,1,df2dy1);
	gsl_matrix_set(m,2,2,df2dy2);
	gsl_matrix_set(m,2,3,df2dy3);
	
	gsl_matrix_set(m,3,0,df3dy0);
	gsl_matrix_set(m,3,1,df3dy1);
	gsl_matrix_set(m,3,2,df3dy2);
	gsl_matrix_set(m,3,3,df3dy3);
	
	dfds[0] = 0.0;
	dfds[1] = 0.0;
	dfds[2] = 0.0;
	dfds[3] = 0.0;
	return GSL_SUCCESS;
}

int main(int argc, char **argv){
	auto app = Gtk::Application::create();
	std::vector<Point> paths[NUM_PATHS];
	for(int i=0;i<NUM_PATHS;i++){
		for(int j=0;j<POINTS_PER_PATH;j++){
			Point p;
			p.y=i*50.0;
			p.x=j*500.0/POINTS_PER_PATH;
			paths[i].push_back(p);
		}
	}

	return app->make_window_and_run<DrawingArea>(argc, argv, paths);
	
	symbol x("x");
	symbol y("y");
	symbol X[2] = {x,y};
	symtab table;
	table["x"] = X[0];
	table["y"] = X[1];
	parser reader(table);
	//:s/k/1.234/g
	//ex g00 = reader("y^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/k))+1/(1-k/(x^2+y^2)^(1/2))");
	//ex g01 = reader("x*y/((x^2+y^2)*((x^2+y^2)^(1/2)/k-1))");
	//ex g10 = reader("x*y/((x^2+y^2)*((x^2+y^2)^(1/2)/k-1))");
	//ex g11 = reader("x^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/k))+1/(1-k/(x^2+y^2)^(1/2))");
	ex g00 = reader("y^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/0.1))+1/(1-0.1/(x^2+y^2)^(1/2))");
	ex g01 = reader("x*y/((x^2+y^2)*((x^2+y^2)^(1/2)/0.1-1))");
	ex g10 = reader("x*y/((x^2+y^2)*((x^2+y^2)^(1/2)/0.1-1))");
	ex g11 = reader("x^2/((x^2+y^2)*(1-(x^2+y^2)^(1/2)/0.1))+1/(1-0.1/(x^2+y^2)^(1/2))");
	matrix g_mat_covariant = {{g00,g01},{g10,g11}};
	matrix g_mat_contravariant = g_mat_covariant.inverse();
	ex g_covariant[2][2] = {{g_mat_covariant(0,0),g_mat_covariant(0,1)},{g_mat_covariant(1,0),g_mat_covariant(1,1)}};
	ex g_contravariant[2][2] = {{g_mat_contravariant(0,0),g_mat_contravariant(0,1)},{g_mat_contravariant(1,0),g_mat_contravariant(1,1)}};
	
	//christoffel symbols
	//G^a_bc=1/2*sum(g^ad(d_b g_dc + d_c g_db - d_d g_bc))
	ex christoffel_symbol_ex[2][2][2];
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				for(int d=0;d<2;d++)
					christoffel_symbol_ex[a][b][c] = 0.5*g_contravariant[a][d]*(g_covariant[d][c].diff(X[b])+g_covariant[d][b].diff(X[c])-g_covariant[b][c].diff(X[d]));
	//christoffel symbol derivatives
	//d_a Gamma^b_cd
	ex christoffel_symbol_derivative_ex[2][2][2][2];
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				for(int d=0;d<2;d++)
					christoffel_symbol_derivative_ex[a][b][c][d] = christoffel_symbol_ex[b][c][d].diff(X[a]);
	
	//compile christoffel symbols
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				compile_ex(christoffel_symbol_ex[a][b][c],x,y,christoffel_symbol[a][b][c]);
	
	//compile christoffel symbol derivatives
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				for(int d=0;d<2;d++)
					compile_ex(christoffel_symbol_derivative_ex[a][b][c][d],x,y,christoffel_symbol_derivative[a][b][c][d]);
	
	gsl_odeiv2_system sys = {func    , jac,      4,         NULL        };
	//                      {function, jacobian, dimension, parameters};
	
	gsl_odeiv2_driver *driver= gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rk8pd,1e-6,1e-6,0.0);
	
	int i=0;
	double t=0.0,t1=100.0;
	//initial conditions
	double Y[4] = {-1.0,1.0,0.2,0.0};
	
	for(i=1;i<=100;i++){
		double ti = i * t1 / 1000.0;
		int status = gsl_odeiv2_driver_apply(driver,&t,ti,Y);
		if(status != GSL_SUCCESS){
			printf("error, return value = %d\n",status);
			break;
		}

		printf("%lf %lf %lf\n",t,Y[0],Y[1]);
	}
	gsl_odeiv2_driver_free(driver);
	
	return 0;
}
