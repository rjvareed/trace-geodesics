#include "DrawingArea.h"

DrawingArea::DrawingArea(int size_x, int size_y, std::vector<Point> *paths){
	set_title("Drawing Area");
	set_default_size(size_x, size_y);
	
	m_Frame_Scribble.set_expand(true);
	set_child(m_Frame_Scribble);
	
	m_Frame_Scribble.set_child(m_DrawingArea_Scribble);
	
	m_DrawingArea_Scribble.set_draw_func(
			sigc::mem_fun(*this, &DrawingArea::on_drawingarea_scribble_draw));
	
	/* Connect signal handlers */
	m_DrawingArea_Scribble.signal_resize().connect(
			sigc::mem_fun(*this, &DrawingArea::on_drawingarea_scribble_resize));
	this->paths = paths;
}

DrawingArea::~DrawingArea(){}

void DrawingArea::on_drawingarea_scribble_draw(const Cairo::RefPtr<Cairo::Context>& cr,int, int){
	cr->set_source(m_surface, 0, 0);
	cr->paint();
}

// Create a new surface of the appropriate size to store our scribbles.
void DrawingArea::scribble_create_surface(){
	int width=m_DrawingArea_Scribble.get_width(),height=m_DrawingArea_Scribble.get_height();
	m_surface = Cairo::ImageSurface::create(Cairo::Surface::Format::ARGB32,width,height);
	
	// Initialize the surface to white.
	auto cr = Cairo::Context::create(m_surface);
	cr->set_source_rgb(1,1,1);
	cr->rectangle(0,0,width,height);
	cr->fill();
	cr->set_source_rgb(0, 0, 0);
	cr->set_line_width(1.0);
	for(int i=0;i<NUM_PATHS;i++){
		for(unsigned long j=0;j<paths[i].size()-1;j++){
			Point p1 = paths[i][j];
			cr->move_to(p1.x,p1.y);
			Point p2 = paths[i][j+1];
			cr->line_to(p2.x,p2.y);
		}
	}
	cr->stroke();
}

void DrawingArea::on_drawingarea_scribble_resize(int /* width */, int /* height */){
	scribble_create_surface();
}
