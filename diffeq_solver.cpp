#include "diffeq_solver.h"
#include "point.h"

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
	
	//df0/ds=d/ds(v_x)=-Gamma^0_ab v^a v^b
	//df1/ds=d/ds(v_y)=-Gamma^1_ab v^a v^b
	double df0ds = 0.0;
	double df1ds = 0.0;
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++){
			df0ds-=christoffel_symbol[0][a][b] (y[0],y[1]) * y[2+a] * y[2+b];
			df1ds-=christoffel_symbol[1][a][b] (y[0],y[1]) * y[2+a] * y[2+b];
		}
	dfds[0] = df0ds;
	dfds[1] = df1ds;
	
	//dv^a/ds = -Gamma^a_bc v^b v^c
	//df2/ds=d/ds(dv_x/ds)=d/ds(-Gamma^0_ab v^a v^b)=-Gamma^0_ab (v^a dv^b/ds + v^b dv^a/ds) =Gamma^0_ab v^a Gamma^b_cd v^c v^d + Gamma^0_ab v^b Gamma^a_cd v^c v^d
	//df3/ds=d/ds(dv_y/ds)=d/ds(-Gamma^1_ab v^a v^b)=-Gamma^1_ab (v^a dv^b/ds + v^b dv^a/ds) =Gamma^1_ab v^a Gamma^b_cd v^c v^d + Gamma^1_ab v^b Gamma^a_cd v^c v^d
	double df2ds = 0.0;
	double df3ds = 0.0;
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				for(int d=0;d<2;d++){
					df2ds += christoffel_symbol[0][a][b] (y[0],y[1])*y[2+a]*christoffel_symbol[b][c][d] (y[0],y[1])*y[2+c]*y[2+d]+christoffel_symbol[0][a][b] (y[0],y[1])*y[2+b]*christoffel_symbol[a][c][d] (y[0],y[1])*y[2+c]*y[2+d];
					df3ds += christoffel_symbol[1][a][b] (y[0],y[1])*y[2+a]*christoffel_symbol[b][c][d] (y[0],y[1])*y[2+c]*y[2+d]+christoffel_symbol[1][a][b] (y[0],y[1])*y[2+b]*christoffel_symbol[a][c][d] (y[0],y[1])*y[2+c]*y[2+d];
				}
	dfds[2] = df2ds;
	dfds[3] = df3ds;
	return GSL_SUCCESS;
}
//map point x [X_MIN,X_MAX] to grid coordinates [0,size_x]
double map_point_x(double x,int size_x){
		return (x-X_MIN)/(X_MAX-X_MIN)*(double)size_x;
}
double map_point_y(double y,int size_y){
		return (double)size_y - (y-Y_MIN)/(Y_MAX-Y_MIN)*(double)size_y;
}
void parse_symbols(std::string s_g00,std::string s_g01,std::string s_g10,std::string s_g11){
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
	ex g00 = reader(s_g00);
	ex g01 = reader(s_g01);
	ex g10 = reader(s_g10);
	ex g11 = reader(s_g11);
	matrix g_mat_covariant = {{g00,g01},{g10,g11}};
	matrix g_mat_contravariant = g_mat_covariant.inverse();
	ex g_covariant[2][2] = {{g_mat_covariant(0,0),g_mat_covariant(0,1)},{g_mat_covariant(1,0),g_mat_covariant(1,1)}};
	ex g_contravariant[2][2] = {{g_mat_contravariant(0,0),g_mat_contravariant(0,1)},{g_mat_contravariant(1,0),g_mat_contravariant(1,1)}};
	
	std::cout << "Calculating Christoffel symbols...\n";
	//christoffel symbols
	//G^a_bc=1/2*sum(g^ad(d_b g_dc + d_c g_db - d_d g_bc))
	ex christoffel_symbol_ex[2][2][2];
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				christoffel_symbol_ex[a][b][c] = 0.0;
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				for(int d=0;d<2;d++)
					christoffel_symbol_ex[a][b][c] += 0.5*g_contravariant[a][d]*(g_covariant[d][c].diff(X[b])+g_covariant[d][b].diff(X[c])-g_covariant[b][c].diff(X[d]));
	
	std::cout << "Calculating Christoffel symbol derivatives...\n";
	//christoffel symbol derivatives
	//d_a Gamma^b_cd
	ex christoffel_symbol_derivative_ex[2][2][2][2];
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				for(int d=0;d<2;d++)
					christoffel_symbol_derivative_ex[a][b][c][d] = christoffel_symbol_ex[b][c][d].diff(X[a]);
	
	std::cout << "Compiling Christoffel symbols...\n";
	//compile christoffel symbols
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				compile_ex(christoffel_symbol_ex[a][b][c],x,y,christoffel_symbol[a][b][c]);
	
	std::cout << "Compiling Christoffel symbol derivatives...\n";
	//compile christoffel symbol derivatives
	for(int a=0;a<2;a++)
		for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				for(int d=0;d<2;d++)
					compile_ex(christoffel_symbol_derivative_ex[a][b][c][d],x,y,christoffel_symbol_derivative[a][b][c][d]);
}
void calculate_paths(std::vector<Point> *paths,int size_x,int size_y){
	gsl_odeiv2_system sys = {func    , jac,      4,         NULL        };
	//                      {function, jacobian, dimension, parameters};
	const int PATHS_PER_RUN = NUM_PATHS / 2;
	const int NUM_STEPS = 1000;
	const double EDGE_TOLERANCE = (Y_MAX-Y_MIN+X_MAX-X_MIN)/4.0 / 20.0;
	const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_rkf45;
	for(int j=0;j<PATHS_PER_RUN;j++){
		std::cout << "Solving differential equation " << j+1 << "/" << NUM_PATHS << "\n";
		gsl_odeiv2_driver *driver= gsl_odeiv2_driver_alloc_y_new(&sys,step_type,1e-6,1e-6,0.0);
		
		int i=0;
		double t=0.0,t1=100.0;
		//initial conditions
		double Y[4] = {X_MIN - EDGE_TOLERANCE,Y_MIN - EDGE_TOLERANCE + (double)j/PATHS_PER_RUN*(Y_MAX-Y_MIN),(Y_MAX-Y_MIN)/10.0,0.0};
		
		for(i=1;i<=NUM_STEPS;i++){
			double ti = i * t1 / (double)NUM_STEPS;
			int status = gsl_odeiv2_driver_apply(driver,&t,ti,Y);
			if(status != GSL_SUCCESS){
				printf("error, return value = %d\n",status);
				break;
			}
			if(Y[0] > X_MAX+EDGE_TOLERANCE || Y[0] < X_MIN-EDGE_TOLERANCE || Y[1] > Y_MAX+EDGE_TOLERANCE || Y[1] < Y_MIN-EDGE_TOLERANCE || std::isnan(Y[0]) || std::isnan(Y[1]))
				break;
			Point p;
			p.x=map_point_x(Y[0],size_x);
			p.y=map_point_y(Y[1],size_y);
			paths[j].push_back(p);
		}
		gsl_odeiv2_driver_free(driver);
	}
	for(int j=0;j<PATHS_PER_RUN;j++){
		std::cout << "Solving differential equation " << j+1+PATHS_PER_RUN << "/" << NUM_PATHS << "\n";
		gsl_odeiv2_driver *driver= gsl_odeiv2_driver_alloc_y_new(&sys,step_type,1e-6,1e-6,0.0);
		
		int i=0;
		double t=0.0,t1=100.0;
		//initial conditions
		double Y[4] = {X_MIN - EDGE_TOLERANCE + (double)j/PATHS_PER_RUN*(X_MAX-X_MIN),Y_MIN - EDGE_TOLERANCE,0.0,(X_MAX-X_MIN)/10.0};
		
		for(i=1;i<=NUM_STEPS;i++){
			double ti = i * t1 / (double)NUM_STEPS;
			int status = gsl_odeiv2_driver_apply(driver,&t,ti,Y);
			if(status != GSL_SUCCESS){
				printf("error, return value = %d\n",status);
				break;
			}
			if(Y[0] > X_MAX+EDGE_TOLERANCE || Y[0] < X_MIN-EDGE_TOLERANCE || Y[1] > Y_MAX+EDGE_TOLERANCE || Y[1] < Y_MIN-EDGE_TOLERANCE || std::isnan(Y[0]) || std::isnan(Y[1]))
				break;
			Point p;
			p.x=map_point_x(Y[0],size_x);
			p.y=map_point_y(Y[1],size_y);
			paths[PATHS_PER_RUN+j].push_back(p);
		}
		gsl_odeiv2_driver_free(driver);
	}
}

