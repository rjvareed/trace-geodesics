#include <ginac/ginac.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <vector>
#include "point.h"

#ifndef __DIFFEQ_SOLVER_H
#define __DIFFEQ_SOLVER_H

int func(double s, const double y[], double f[], void *params);
int jac(double t,const double y[],double *dfdy, double dfds[], void *params);
void calculate_paths(std::vector<Point> *paths,int size_x,int size_y);
double map_point_x(double x,int size_x);
double map_point_y(double y,int size_y);

#endif
