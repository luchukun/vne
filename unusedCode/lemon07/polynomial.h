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

#ifndef LEMON_BEZIER_H
#define LEMON_BEZIER_H

///\ingroup misc
///\file
///\brief A simple class implementing polynomials.
///
///\author Alpar Juttner

#include<vector>

namespace lemon {

  /// \addtogroup misc
  /// @{

  ///Simple polinomial class

  ///This class implements a polynomial where the coefficients are of
  ///type \c T.
  ///
  ///The coefficients are stored in an std::vector.
  template<class T>
  class Polynomial
  {
    std::vector<T> _coeff;
  public:
    ///Construct a polynomial of degree \c d.
    explicit Polynomial(int d=0) : _coeff(d+1) {}
    ///\e
    template<class U> Polynomial(const U &u) : _coeff(1,u) {}
    ///\e
    template<class U> Polynomial(const Polynomial<U> &u) : _coeff(u.deg()+1)
    {
      for(int i=0;i<int(_coeff.size());i++) _coeff[i]=u[i];
    }
    ///Query the degree of the polynomial.
    
    ///Query the degree of the polynomial.
    ///\warning This number differs from real degree of the polinomial if
    ///the coefficient of highest degree is 0.
    int deg() const { return _coeff.size()-1; }
    ///Set the degree of the polynomial.

    ///Set the degree of the polynomial. In fact it resizes the
    ///coefficient vector.
    void deg(int d) { _coeff.resize(d+1);}

    ///Returns (as a reference) the coefficient of degree \c d.
    typename std::vector<T>::reference operator[](int d) { return _coeff[d]; }
    ///Returns (as a const reference) the coefficient of degree \c d.
    typename std::vector<T>::const_reference
    operator[](int d) const {return _coeff[d];}
    
    ///Substitute the value u into the polinomial.

    ///Substitute the value u into the polinomial.
    ///The calculation will be done using type \c R.
    ///The following examples shows the usage of the template parameter \c R.
    ///\code
    ///  Polynomial<dim2::Point<double> > line(1);
    ///  line[0]=dim2::Point<double>(12,25);
    ///  line[1]=dim2::Point<double>(2,7);
    ///  ...
    ///  dim2::Point<double> d = line.subst<dim2::Point<double> >(23.2);
    ///\endcode
    ///
    ///\code
    ///  Polynomial<double> p;
    ///  Polynomial<double> q;
    ///  ...
    ///  Polynomial<double> s = p.subst<Polynomial<double> >(q);
    ///\endcode
    template<class R,class U>
    R subst(const U &u) const
    {
      typename std::vector<T>::const_reverse_iterator i=_coeff.rbegin();
      R v=*i;
      for(++i;i!=_coeff.rend();++i) v=v*u+*i;
      return v;
    }
    ///Substitute the value u into the polinomial.

    ///Substitute the value u into the polinomial.
    ///The calculation will be done using type \c T
    ///(i.e. using the type of the coefficients.)
    template<class U>
    T operator()(const U &u) const 
    {
      return subst<T>(u);
    }
    
    ///Derivate the polynomial (in place)
    Polynomial &derivateMyself()
    {
      for(int i=1;i<int(_coeff.size());i++) _coeff[i-1]=i*_coeff[i];
      _coeff.pop_back();
      return *this;
    }
    
    ///Return the derivate of the polynomial
    Polynomial derivate() const
    {
      Polynomial tmp(deg()-1);
      for(int i=1;i<int(_coeff.size());i++) tmp[i-1]=i*_coeff[i];
      return tmp;
    }

    ///Integrate the polynomial (in place)
    Polynomial &integrateMyself()
    {
      _coeff.push_back(T());
      for(int i=_coeff.size()-1;i>0;i--) _coeff[i]=_coeff[i-1]/i;
      _coeff[0]=0;
      return *this;
    }
    
    ///Return the integrate of the polynomial
    Polynomial integrate() const
    {
      Polynomial tmp(deg()+1);
      tmp[0]=0;
      for(int i=0;i<int(_coeff.size());i++) tmp[i+1]=_coeff[i]/(i+1);
      return tmp;
    }

    ///\e
    template<class U>
    Polynomial &operator+=(const Polynomial<U> &p)
    {
      if(p.deg()>deg()) _coeff.resize(p.deg()+1);
      for(int i=0;i<=int(std::min(deg(),p.deg()));i++)
	_coeff[i]+=p[i];
      return *this;
    }
    ///\e
    template<class U>
    Polynomial &operator-=(const Polynomial<U> &p)
    {
      if(p.deg()>deg()) _coeff.resize(p.deg()+1);
      for(int i=0;i<=std::min(deg(),p.deg());i++) _coeff[i]-=p[i];
      return *this;
    }
    ///\e
    template<class U>
    Polynomial &operator+=(const U &u)
    {
      _coeff[0]+=u;
      return *this;
    }
    ///\e
    template<class U>
    Polynomial &operator-=(const U &u)
    {
      _coeff[0]+=u;
      return *this;
    }
    ///\e
    template<class U>
    Polynomial &operator*=(const U &u)
    {
      for(typename std::vector<T>::iterator i=_coeff.begin();i!=_coeff.end();++i)
	*i*=u;
      return *this;
    }
    ///\e
    template<class U>
    Polynomial &operator/=(const U &u)
    {
      for(typename std::vector<T>::iterator i=_coeff.begin();i!=_coeff.end();++i)
	*i/=u;
      return *this;
    }

  };
  
  ///Equality comparison

  ///\relates Polynomial
  ///\warning Two polynomials are defined to be unequal if their degrees differ,
  ///even if the non-zero coefficients are the same.
  template<class U,class V>
  bool operator==(const Polynomial<U> &u,const Polynomial<V> &v)
  {
    if(u.deg()!=v.deg()) return false;
    for(int i=0;i<=u.deg();i++) if(u[i]!=v[i]) return false;
    return true;
  }

  ///Non-equality comparison

  ///\relates Polynomial
  ///\warning Two polynomials are defined to be unequal if their degrees differ,
  ///even if the non-zero coefficients are the same.
  template<class U,class V>
  bool operator!=(const Polynomial<U> &u,const Polynomial<V> &v)
  {
    return !(u==v);
  }

  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator+(const Polynomial<U> &u,const Polynomial<V> &v)
  {
    Polynomial<U> tmp=u;
    tmp+=v;
    return tmp;
  }

  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator-(const Polynomial<U> &u,const Polynomial<V> &v)
  {
    Polynomial<U> tmp=u;
    tmp-=v;
    return tmp;
  }

  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator*(const Polynomial<U> &u,const Polynomial<V> &v)
  {
    Polynomial<U> tmp(u.deg()+v.deg());
    for(int i=0;i<=v.deg();i++)
      for(int j=0;j<=u.deg();j++)
	tmp[i+j]+=v[i]*u[j];
    return tmp;
  }
  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator+(const Polynomial<U> &u,const V &v)
  {
    Polynomial<U> tmp=u;
    tmp+=v;
    return tmp;
  }
  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator+(const V &v,const Polynomial<U> &u)
  {
    Polynomial<U> tmp=u;
    tmp+=v;
    return tmp;
  }
  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator-(const Polynomial<U> &u,const V &v)
  {
    Polynomial<U> tmp=u;
    tmp-=v;
    return tmp;
  }
  ///\e

  ///\relates Polynomial
  ///
  template<class U>
  Polynomial<U> operator-(const Polynomial<U> &u)
  {
    Polynomial<U> tmp(u.deg());
    for(int i=0;i<=u.deg();i++) tmp[i]=-u[i];
    return tmp;
  }
  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator-(const V &v,const Polynomial<U> &u)
  {
    Polynomial<U> tmp=-u;
    tmp+=v;
    return tmp;
  }
  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator*(const Polynomial<U> &u,const V &v)
  {
    Polynomial<U> tmp=u;
    tmp*=v;
    return tmp;
  }
  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator*(const V &v,const Polynomial<U> &u)
  {
    Polynomial<U> tmp=u;
    tmp*=v;
    return tmp;
  }
  ///\e

  ///\relates Polynomial
  ///
  template<class U,class V>
  Polynomial<U> operator/(const Polynomial<U> &u,const V &v)
  {
    Polynomial<U> tmp=u;
    tmp/=v;
    return tmp;
  }
    
  /// @}

} //END OF NAMESPACE LEMON

#endif // LEMON_POLYNOMIAL_H
