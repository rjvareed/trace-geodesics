#include <gtkmm.h>
#include <vector>
#include "point.h"

#ifndef __DRAWINGAREA_H
#define __DRAWINGAREA_H

class DrawingArea : public Gtk::Window{
public:
	DrawingArea(int size_x, int size_y, std::vector<Point> *paths);
	~DrawingArea() override;
	
protected:
	//draw functions:
	void on_drawingarea_scribble_draw(const Cairo::RefPtr<Cairo::Context>& cr, int width, int height);
	//signal handlers:
	void on_drawingarea_scribble_resize(int width, int height);
	void scribble_create_surface();
	void scribble_draw_brush(double x, double y);
	//Member widgets:
	Gtk::Frame m_Frame_Scribble;
	Gtk::DrawingArea m_DrawingArea_Scribble;
	Cairo::RefPtr<Cairo::ImageSurface> m_surface;
	
private:
	std::vector<Point> *paths;
};

#endif
