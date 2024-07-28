#include <gtkmm.h>
#include <vector>

#ifndef __DRAWINGAREA_H
#define __DRAWINGAREA_H

#define NUM_PATHS 10
#define POINTS_PER_PATH 10
struct Point{
	double x;
	double y;
};

class DrawingArea : public Gtk::Window{
public:
	DrawingArea(std::vector<Point> *paths);
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
