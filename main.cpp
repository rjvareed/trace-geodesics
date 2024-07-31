#include <iostream>
#include "diffeq_solver.h"
#include "DrawingArea.h"
#include "point.h"

int main(int argc, char **argv){
	int size_x = 500, size_y = 500;
	auto app = Gtk::Application::create();
	std::string g00,g01,g10,g11;
	std::cout << "Enter four covariant metric tensor entries ({{g00,g01},{g10,g11}} for Cartesian coordinates, separated by newlines:\n";
	std::getline(std::cin,g00);
	std::getline(std::cin,g01);
	std::getline(std::cin,g10);
	std::getline(std::cin,g11);
	parse_symbols(g00,g01,g10,g11);
	
	std::vector<Point> *paths = new std::vector<Point>[NUM_PATHS]();
	calculate_paths(paths,size_x,size_y);
	
	return app->make_window_and_run<DrawingArea>(argc, argv, size_x, size_y, paths);
}
