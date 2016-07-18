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

#ifndef LEMON_EPS_H
#define LEMON_EPS_H

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<lemon/color.h>
#include<lemon/dim2.h>

  ///\ingroup eps_io
  ///\file
  ///\brief Simple tool to create \c .eps files
  ///
  ///\author Alpar Juttner

namespace lemon {
  
  ///\ingroup eps_io
  ///\brief A simple tool to create \c .eps files
  ///
  ///A simple tool to create \c .eps files
  ///\author Alpar Juttner
  class EpsDrawer
  {
    void defMacros();
    void init(int x1,int y1,int x2,int y2);
    void init(double x1,double y1,double x2,double y2);
    bool local_stream;
  public:
    
    std::ostream &out;
    
    ///Node shapes
    ///
    enum NodeShapes { 
      /// = 0
      ///\image html nodeshape_0.png
      ///\image latex nodeshape_0.eps "CIRCLE shape (0)" width=2cm
      CIRCLE=0, 
      /// = 1
      ///\image html nodeshape_1.png
      ///\image latex nodeshape_1.eps "SQUARE shape (1)" width=2cm
      ///
      SQUARE=1, 
      /// = 2
      ///\image html nodeshape_2.png
      ///\image latex nodeshape_2.eps "DIAMOND shape (2)" width=2cm
      ///
      DIAMOND=2,
      /// = 3
      ///\image html nodeshape_3.png
      ///\image latex nodeshape_2.eps "MALE shape (4)" width=2cm
      ///
      ///\warning Not implemented
      MALE=3,
      /// = 4
      ///\image html nodeshape_4.png
      ///\image latex nodeshape_2.eps "FEMALE shape (4)" width=2cm
      ///
      ///\warning Not implemented
      FEMALE=4
    };
    ///\e 

    ///The generated file is put to \c os.
    ///
    /// \c x and \c y determine the upper
    ///right corner of the bounding box. The lower left corner is (0,0).
    EpsDrawer(std::ostream &os,int x,int y);
    ///\e

    ///The generated file is put to \c os.
    ///
    ///(x1,y1) and (x2,y2)
    /// determine the lower left and the upper right corners of
    ///the bounding box, respectively.
    EpsDrawer(std::ostream &os,int x1,int y1,int x2,int y2);
    ///\e

    ///The generated file is put to \c os.
    ///
    ///\c s determines the upper
    ///right corner of the bounding box. The lower left corner is (0,0).
    EpsDrawer(std::ostream &os,dim2::Point<double> s);
    ///\e

    ///The generated file is put to \c os.
    ///
    ///\c a and \c b
    /// determine the lower left and the upper right corners of
    ///the bounding box, respectively.
    EpsDrawer(std::ostream &os,dim2::Point<double> a, dim2::Point<double> b);
    ///\e

    ///The generated picture is put to the file \c name.
    ///
    ///\c x and \c y determine the upper
    ///right corner of the bounding box. The lower left corner is (0,0).
    EpsDrawer(const std::string &name,int x,int y);
    ///\e

    ///The generated picture is put to the file \c name.
    ///
    ///(x1,y1) and (x2,y2)
    /// determine the lower left and the upper right corners of
    ///the bounding box, respectively.
    EpsDrawer(const std::string &name,int x1,int y1,int x2,int y2);
    ///\e

    ///The generated picture is put to the file \c name.
    ///
    ///\c s determines the upper
    ///right corner of the bounding box. The lower left corner is (0,0).
    EpsDrawer(const std::string &name,dim2::Point<double> s);
    ///\e

    ///The generated picture is put to the file \c name.
    ///
    ///\c a and \c b
    /// determine the lower left and the upper right corners of
    ///the bounding box, respectively.
    EpsDrawer(const std::string &name,dim2::Point<double> a, dim2::Point<double> b);

//     template<class T> EpsDrawer(std::ostream &os,BoundingBox<T> b) 
//     template<class T> EpsDrawer(std::ostream &os,BoundingBox<T> b);
    
    ~EpsDrawer();
    
    ///Save the current graphic state.

    ///This function saves the current coordinate system, and the parameters
    ///set by \ref color(), \ref lineWidth() etc.
    ///The can be \ref restore "restore()"d later.
    ///
    ///The \ref save() - \ref restore() pairs can be nested.
    ///
    EpsDrawer &save();
    ///Restore the saves graphic state.

    EpsDrawer &restore();
    
    ///Draw a line
    EpsDrawer &line(double x1,double y1,double x2,double y2);
    ///Draw a line from the current point
    EpsDrawer &lineTo(double x1,double y1);
    ///Move the current point
    EpsDrawer &moveTo(double x1,double y1);
    ///Draw a circle
    EpsDrawer &circle(double x,double y, double r);
    
    ///Draw a line
    template<class T> EpsDrawer &line(dim2::Point<T> p1,dim2::Point<T> p2) 
    {
      return line(p1.x,p1.y,p2.x,p2.y);
    }
    ///Draw a line from the current point
    template<class T> EpsDrawer &lineTo(dim2::Point<T> p)
    {
      return lineTo(p.x,p.y);
    }
    ///Move the current point
    template<class T> EpsDrawer &moveTo(dim2::Point<T> p)
    {
      return moveTo(p.x,p.y);
    }
    ///Draw a circle
    template<class T> EpsDrawer &circle(dim2::Point<T> p, double r)
    {
      return circle(p.x,p.y,r);
    }
    
    ///Set the font size
    EpsDrawer &fontSize(double si);
    ///Set the fint type
    EpsDrawer &font(std::string );
    ///Sets whether text output is centerized of not

    ///Sets whether text output is centerized of not.
    ///
    ///\warning \ref save() doesn't save this setting.
    ///
    EpsDrawer &centerMode(bool m);
    ///Turn to collect mode.

    ///If you call this function, then the drawing operations like \ref line(),
    ///\ref lineTo(), \ref moveTo() will not take place immediately, but istead
    ///they
    ///are collected. These operations form a \e path.
    ///Then you can \ref stroke(), \ref fill(), \ref eoFill(), \ref clip() or
    ///\ref eoClip() it.
    ///When drawing, you can also use \ref closePath() to - surprise - close the
    ///current path.
    ///
    ///This example draws a red filled diamond.
    ///\code
    ///  EpsDraw ed("diamond.eps",-1,-1,1,1);
    ///  ed.color(1,0,0,).collect().line(0,-1,1,0).lineTo(0,1)
    ///    .lineTo(-1,0).closePath().fill();
    ///\endcode
    EpsDrawer &collect();
    ///Close the current drawing path
    EpsDrawer &closePath();
    ///Stroke (draw) a path

    ///Stroke (draw) a path.
    ///\sa collect()
    ///
    EpsDrawer &stroke();
    ///Fill a path

    ///Fill a path.
    ///\sa collect()
    ///
    EpsDrawer &fill();
    ///Even-odd fill a path

    ///Even-odd fill a path.
    ///\sa collect()
    ///
    EpsDrawer &eoFill();
    ///Set a clipping area.

    ///This function sets a clipping area. After that, the drawing operations
    ///will affect only this area.
    ///\sa collect()
    ///
    EpsDrawer &clip();
    ///Set a clipping area using even-odd rule.

    ///This function sets a clipping area using even-odd rule.
    ///After that, the drawing operations
    ///will affect only this area.
    ///\sa collect()
    ///
    EpsDrawer &eoClip();
    
    ///Set the line width.
    EpsDrawer &lineWidth(double w);
    ///Set the style of the line ends

    ///\param i It can be 0, 1 or 2
    ///
    EpsDrawer &lineCap(int i);
    ///Set the style of the line joins

    ///\param i It can be 0, 1 or 2
    ///
    EpsDrawer &lineJoin(int i);
    ///Set the cut length of near parallel joining lines.
    EpsDrawer &miterLimit(double w);
    
    ///Set the drawing color
    EpsDrawer &color(double r, double g, double b);
    ///Set the drawing color
    EpsDrawer &color(Color c)
    {
      return color(c.red(),c.green(),c.blue());
    }
    
    ///Draw a node shape

    ///Draw a node shape.
    ///
    ///\param t The shape of the drawn object
    ///\param x The \c x coordinate of the node
    ///\param y The \c y coordinate of the node
    ///\param r The size (radius) of the node
    ///\param col Color of the node. The default color is white
    ///\param brd Color of the node border. The default color is black
    EpsDrawer &node(NodeShapes t, double x, double y, double r,
		    Color col=WHITE, Color brd=BLACK);
    ///Draw a node shape
    
    ///Draw a node shape.
    ///
    ///\param t The shape of the drawn object
    ///\param pos Position of the node
    ///\param r The size (radius) of the node
    ///\param col Color of the node. The default color is white
    ///\param brd Color of the node border. The default color is black
    template<class T>
    EpsDrawer &node(NodeShapes t, dim2::Point<T> pos, double r,
		    Color col=WHITE, Color brd=BLACK)
    {
      return node(t,pos.x,pos.y,r,col,brd);
    }

    ///Translate the coordinate system
    EpsDrawer &translate(double x,double y);
    ///Translate the coordinate system
    template<class T> EpsDrawer &translate(dim2::Point<T> p)
    {
      return translate(p.x,p.y);
    }
    ///Rotate the coordinate system
    EpsDrawer &rotate(double r);
    ///Scale the coordinate system
    EpsDrawer &scale(double sx, double sy);
    ///Scale the coordinate system
    EpsDrawer &scale(double s) { return scale(s,s); }
    ///Scale the coordinate system
    template<class T> EpsDrawer &scale(dim2::Point<T> p)
    {
      return scale(p.x,p.y);
    }
    
    ///\e
    EpsDrawer &flush();
    ///Clear the image
    EpsDrawer &clear();
    
    ///Print a text at the current point
    EpsDrawer &operator<<(const std::string &s);
    ///Print a text at the current point
    EpsDrawer &operator<<(const char *s);
    ///Print a number at the current point
    EpsDrawer &operator<<(int i);
    ///Print a number at the current point
    EpsDrawer &operator<<(double d);
    ///Print a coordinate at the current point
    template<class T>
    EpsDrawer &operator<<(dim2::Point<T> p) 
    {
      out << "((" << p.x << ',' << p.y <<")) show\n";
      return *this;
    }
    
  };
  
}

#endif
