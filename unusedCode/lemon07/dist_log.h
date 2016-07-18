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

#ifndef LEMON_RANDOM_H
#define LEMON_RANDOM_H

#include<iostream>
#include<fstream>
#include<string>

#include <lemon/dim2.h>

///\ingroup misc
///\file
///\brief Measure a Distribution
///
///\todo Needs lot more docs
///


namespace lemon {

  ///Measure a distribution
  class DistLog
  {
    std::vector<int> _dist;
    double _lo,_up;
    int _count;
    bool _cut;
  public:
    ///\e
    Dist(double lo,double up,int gr,bool cut=true) 
      : _dist(gr,0),_lo(lo),_up(up),_count(0),_cut(cut) {}
    ///\e
    void operator()(double v)
    {
      if(_cut) {
	if(_lo<=v && v<_up)
	  _dist[int((v-_lo)/(_up-_lo)*_dist.size())]++;
      }
      else {
	_dist[std::max(0,std::min(int(_dist.size())-1,
				  int((v-_lo)/(_up-_lo)*_dist.size())
				  ))]++;
      }
      _count++;
    }
    ///\e
    void dump(std::ostream& os=std::cout)
    {
      for(int i=0;i<_dist.size();i++)
	os << _lo+(i+0.5)*(_up-_lo)/_dist.size() << ' '
	   << double(_dist[i])/_count << std::endl;
    }
    ///\e
    void dump(const std::string& file_name)
    {
      dump(std::ofstream(file_name.c_str()));
    }
  };
  
  ///Measure a two dimensional distribution
  class DistLog2
  {
  public:
    typedef dim2::Point<double> Point;
  private:
    std::vector<int> _dist;
    int _gr;
    Point _lo,_up;
    int _count;
    bool _cut;
  public:  
    ///\e
    Dist2(Point a,Point b,int gr,bool cut=true) :
      _dist(gr*gr,0),_gr(gr),
      _lo(a),_up(b),_count(0),_cut(cut) {}
    ///\e
    Dist2(double lox,double upx,double loy,double upy,int gr,bool cut=true) :
      _dist(gr*gr,0),_gr(gr),
      _lo(Point(lox,loy)),_up(Point(upx,upy)),_count(0),_cut(cut) {}
    ///\e
    void operator()(Point v)
    {
      if(_cut)
	{
	  if(v.x>=_lo.x && v.x<_up.x && v.y>=_lo.y && v.y<_up.y)
	    _dist[int((v.x-_lo.x)/(_up.x-_lo.x)*_gr)*_gr+
		  int((v.y-_lo.y)/(_up.y-_lo.y)*_gr)]++;
	}
      else {
	_dist[std::max(0,std::min(_gr-1,
				  int((v.x-_lo.x)/(_up.x-_lo.x)*_gr)
				  ))*_gr+
	      std::max(0,std::min(_gr-1,
				  int((v.y-_lo.y)/(_up.y-_lo.y)*_gr)
				  ))
	      ]++;
      }
      _count++;
    }
    ///\e
    void dump(std::ostream& os=std::cout)
    {
      for(int i=0;i<_gr;i++)
	{
	  for(int j=0;j<_gr;j++)
	    os << _lo.x+(i+0.5)*(_up.x-_lo.x)/_gr << ' '
	       << _lo.y+(j+0.5)*(_up.y-_lo.y)/_gr << ' '
	       << double(_dist[i*_gr+j])/_count << std::endl;
	  os << std::endl;
	}
    }
    ///\e
    void dump(const std::string& file_name)
    {
      dump(std::ofstream(file_name.c_str()));
    }
  };
}

#endif
