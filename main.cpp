#include <iostream>
#include "diffeq_solver.h"
#include "DrawingArea.h"
#include "point.h"

int main(int argc, char **argv){

	int size_x = 500, size_y = 500;
	auto app = Gtk::Application::create();
	std::vector<Point> *paths = new std::vector<Point>[NUM_PATHS]();
	
	calculate_paths(paths,size_x,size_y);

	return app->make_window_and_run<DrawingArea>(argc, argv, size_x, size_y, paths);
}
