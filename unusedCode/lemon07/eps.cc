/* -*- C++ -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library
 *
 * Copyright (C) 2003-2008
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#include <lemon/eps.h>

namespace lemon {
  
  void EpsDrawer::defMacros()
  {
    out << "/clmode true def\n" <<
      "/cshowmode false def\n" <<
      "/defont (Helvetica) findfont def\n" <<
      "/fosi 12 def\n" <<
      "\n" <<
      "/st { clmode { currentpoint stroke newpath moveto } if } bind def\n" <<
      "/str { currentpoint stroke newpath moveto /clmode true def } bind def\n"
	<<
      "/fl { currentpoint fill newpath moveto /clmode true def } bind def\n" <<
      "/eofl { currentpoint eofill newpath moveto /clmode true def } bind def\n"
	<<
      "/cl { currentpoint clip newpath moveto /clmode true def } bind def\n"
	<<
      "/eocl { currentpoint eoclip newpath moveto /clmode true def } bind def\n"
	<<
      "\n" <<
      "/l { moveto lineto st } bind def\n" <<
      "/lt { lineto st } bind def\n" <<
      "/mt { moveto } bind def\n" <<
      "/c { dup 3 index add 2 index moveto 0 360 arc st } bind def\n" <<
      "/collect { /clmode false def currentpoint newpath moveto } bind def\n" <<
      "\n" <<
      "/fontset { defont fosi scalefont setfont } bind def\n" <<
      "/stfs { /fosi exch def fontset } bind def\n" <<
      "/cshow { dup stringwidth pop\n" <<
      "   neg 2 div 0 rmoveto show } bind def\n" <<
      "/xshow { cshowmode { cshow } { show } ifelse } def\n" <<
      "\n" <<
      "fontset\n" <<
      "newpath\n" <<
      "0 0 moveto\n" <<
      "1 setlinecap\n";
  }

  void EpsDrawer::init(int x1,int y1,int x2,int y2)
  {
    out << "%!PS-Adobe-2.0 EPSF-2.0\n" <<
      "%%BoundingBox: " << x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 <<
      "\n%%EndComments\n";
    defMacros();
  }

  void EpsDrawer::init(double x1,double y1,double x2,double y2)
  {
    out << "%!PS-Adobe-2.0\n" <<
      "%%HiResBoundingBox: " << 
      x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 <<
      "\n%%EndComments\n";
    defMacros();
  }


  EpsDrawer::EpsDrawer(std::ostream &os,int x,int y) : local_stream(false),
						       out(os)
  {
    init(0,0,x,y);
  }

  EpsDrawer::EpsDrawer(std::ostream &os,int x1,int y1,int x2,int y2) : 
    local_stream(false),
    out(os)
  {
    init(x1,y1,x2,y2);
  }

  EpsDrawer::EpsDrawer(std::ostream &os,dim2::Point<double> s) : local_stream(false),
							out(os)
  {
    init(0.0,0.0,s.x,s.y);
  }

  EpsDrawer::EpsDrawer(std::ostream &os,dim2::Point<double> a, dim2::Point<double> b) :
    local_stream(false),
    out(os)
  {
    init(a.x,a.y,b.x,b.y);
  }


  EpsDrawer::EpsDrawer(const std::string &name,int x,int y) :
    local_stream(true),
    out(*new std::ofstream(name.c_str()))
  {
    init(0,0,x,y);
  }

  EpsDrawer::EpsDrawer(const std::string &name,int x1,int y1,int x2,int y2) : 
    local_stream(true),
    out(*new std::ofstream(name.c_str()))
  {
    init(x1,y1,x2,y2);
  }
  
  EpsDrawer::EpsDrawer(const std::string &name,dim2::Point<double> s) :
    local_stream(true),
    out(*new std::ofstream(name.c_str()))
  {
    init(0.0,0.0,s.x,s.y);
  }

  EpsDrawer::EpsDrawer(const std::string &name,dim2::Point<double> a, dim2::Point<double> b) :
    local_stream(true),
    out(*new std::ofstream(name.c_str()))
  {
    init(a.x,a.y,b.x,b.y);
  }


  EpsDrawer::~EpsDrawer()
  {
    out << "showpage\n";
    if(local_stream) delete &out;
  }

  EpsDrawer &EpsDrawer::save()
  {
    out << "gsave\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::restore()
  {
    out << "grestore\n";
    return *this;  
  }
 
  EpsDrawer &EpsDrawer::line(double x1,double y1,double x2,double y2)
  {
    out << x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 << " l\n";
    return *this;
  
  }

  EpsDrawer &EpsDrawer::lineTo(double x1,double y1)
  {
    out << x1 << ' ' << y1 << " lt\n";
    return *this;
  
  }

  EpsDrawer &EpsDrawer::moveTo(double x1,double y1)
  {
    out << x1 << ' ' << y1 << " mt\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::circle(double x,double y, double r)
  {
    out << x << ' ' << y << ' ' << r << " c\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::operator<<(const std::string &s)
  {
    out << "(" << s <<") xshow\n";
    return *this;
  }

  EpsDrawer &EpsDrawer::operator<<(const char *s)
  {
    out << "(" << s <<") xshow\n";
    return *this;
  }

  EpsDrawer &EpsDrawer::operator<<(int i)
  {
    out << "(" << i <<") xshow\n";
    return *this;
  }

  EpsDrawer &EpsDrawer::operator<<(double d)
  {
    out << "(" << d <<") xshow\n";
    return *this;
  }

  EpsDrawer &EpsDrawer::fontSize(double si)
  {
    out << si << " stfs\n";
    return *this;
  }
  EpsDrawer &EpsDrawer::font(std::string s)
  {
    out << "/defont ("<<s<<") findfont def fontset\n";
    return *this;
  }


  EpsDrawer &EpsDrawer::collect()
  {
    out << "collect\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::closePath()
  {
    out << "closepath\n";
    return *this;
  }

  EpsDrawer &EpsDrawer::stroke()
  {
    out << "str\n";
    return *this;  
  }
  EpsDrawer &EpsDrawer::fill()
  {
    out << "fl\n";
    return *this;  
  }
  EpsDrawer &EpsDrawer::eoFill()
  {
    out << "eofl\n";
    return *this;  
  }
  EpsDrawer &EpsDrawer::clip()
  {
    out << "cl\n";
    return *this;  
  }
  EpsDrawer &EpsDrawer::eoClip()
  {
    out << "eocl\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::lineWidth(double w)
  {
    out << w << " setlinewidth\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::lineCap(int i)
  {
    out << i << " setlinecap\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::lineJoin(int i)
  {
    out << i << " setlinejoin\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::miterLimit(double w)
  {
    out << w << " setmiterlimit\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::color(double r, double g, double b)
  {
    out << r << ' ' << g << ' ' << b << " setrgbcolor\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::translate(double x,double y)
  {
    out << x << ' ' << y << " translate\n";
    return *this;  
  }

  EpsDrawer &EpsDrawer::rotate(double r)
  {
    out << r << " rotate\n";
    return *this;  
  }
  EpsDrawer &EpsDrawer::scale(double sx, double sy)
  {
    out << sx << ' ' << sy << " scale\n";
    return *this;  
  }
  
  EpsDrawer &EpsDrawer::clear()
  {
    out << "erasepage\n";
    return *this;  
  }
  
  EpsDrawer &EpsDrawer::centerMode(bool m)
  {
    if(m) out << "/cshowmode true def\n";
    else out << "/cshowmode false def\n";

    return *this;  
  }
  
  EpsDrawer &EpsDrawer::flush()
  {
    out << "flush\n";
    //  fflush(fp);
    return *this;
  }

  EpsDrawer &EpsDrawer::node(NodeShapes t, double x, double y, double r,
			     Color col, Color brd)
  {
    out << "gsave\n"
	<< brd.red() << ' ' << brd.green() << ' ' << brd.blue() 
	<< " setrgbcolor\n";
    switch(t) {
    case CIRCLE:
      out << "newpath " << x << ' ' << y << ' ' << r 
	  << " dup 3 index add 2 index moveto 0 360 arc fill\n";
      break;
    case SQUARE:
      out << "newpath\n"
	  << x-r << ' ' << y-r << " moveto\n"
	  << x-r << ' ' << y+r << " lineto\n"
	  << x+r << ' ' << y+r << " lineto\n"
	  << x+r << ' ' << y-r << " lineto closepath fill\n";
      break;
    case DIAMOND:
      out << "newpath\n"
	  << x-r << ' ' << y   << " moveto\n"
	  << x   << ' ' << y+r << " lineto\n"
	  << x+r << ' ' << y   << " lineto\n"
	  << x   << ' ' << y-r << " lineto closepath fill\n";
      break;
    case MALE:
      break;
    case FEMALE:
      break;
    }
    r/=1.1;
    out << col.red() << ' ' << col.green() << ' ' << col.blue() 
	<< " setrgbcolor\n";
    switch(t) {
    case CIRCLE:
      out << "newpath " << x << ' ' << y << ' ' << r 
	  << " dup 3 index add 2 index moveto 0 360 arc fill\n";
      break;
    case SQUARE:
      out << "newpath\n"
	  << x-r << ' ' << y-r << " moveto\n"
	  << x-r << ' ' << y+r << " lineto\n"
	  << x+r << ' ' << y+r << " lineto\n"
	  << x+r << ' ' << y-r << " lineto closepath fill\n";
      break;
    case DIAMOND:
      out << "newpath\n"
	  << x-r << ' ' << y   << " moveto\n"
	  << x   << ' ' << y+r << " lineto\n"
	  << x+r << ' ' << y   << " lineto\n"
	  << x   << ' ' << y-r << " lineto closepath fill\n";
      break;
    case MALE:
      break;
    case FEMALE:
      break;
    }

    out << "grestore\n";
    return *this;
  }
  
}
