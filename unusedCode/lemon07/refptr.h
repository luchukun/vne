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

#ifndef LEMON_REFPTR_H
#define LEMON_REFPTR_H

///\ingroup misc
///\file
///\brief A reference counted pointer implementation.
///
///\todo Undocumented


namespace lemon {

  
  ///Reference counted pointer

  ///This is a simple implementation of a reference counted pointer.
  ///
  ///\warning Current implementation is far from being thread-safe.
  template<class T>
  class RefPtr 
  {
    mutable RefPtr *prev, *next;
    T *ref;

    void lock() {}
    void unlock() {}
    
    void attach(RefPtr &r) 
    {
      if(r.ref) {
	prev=&r; next=r.next; ref=r.ref;
	r.next=this;
      }
    }
    void attach(const T *p) 
    {
      prev=0; next=0; ref=p;
    }
    void release() 
    {
      if(ref) {
	bool fr=true;
	if(prev) { fr=false; prev->next=next; }
	if(next) { fr=false; next->prev=prev; }
	if(fr) delete ref;
	ref=0;
      }
    }
  
  public:
    ///\e
    RefPtr() : ref(0) {}

    ///\e
    RefPtr(const RefPtr &r) {
      lock();
      attach(const_cast<RefPtr&>(r));
      unlock();
    }

    ///\e
    RefPtr(T *p) : prev(0), next(0), ref(p) {}

    ///\e
    ~RefPtr() {
      lock();
      release();
      unlock();
    }

    ///\e
    const RefPtr &operator=(const RefPtr &r) { 
      if(ref!=r.ref) {
	lock();
	release(); attach(const_cast<RefPtr&>(r));
	unlock();
      }
      return *this;
    }
  
    ///\e
    const RefPtr &operator=(const T* &p) { 
      if(ref!=p) { lock(); release(); attach(p); unlock(); }
      return *this;
    }
  
    ///\e
    void swap(RefPtr &r) {
      RefPtr *p;
      T *tp;
      lock();
      p=prev; prev=r.prev; r.prev=p;
      p=next; next=r.next; r.next=p;
      tp=ref; ref=r.ref; r.ref=tp;
      unlock();
    }

    ///\e
    void clear() { lock(); release(); unlock(); }

    ///\e
    T * operator->() { return ref; }
    ///\e
    const T * operator->() const { return ref; }
    ///\e
    operator T *() { return ref; }
    ///\e
    operator const T *() const { return ref; }
    
    ///\e
    bool operator<(const RefPtr &r) const { return this->ref < r.ref; }
    ///\e
    bool operator<=(const RefPtr &r) const { return this->ref <= r.ref; }
    ///\e
    bool operator==(const RefPtr &r) const { return this->ref == r.ref; }
    ///\e
    bool operator>=(const RefPtr &r) const { return this->ref >= r.ref; }
    ///\e
    bool operator>(const RefPtr &r) const { return this->ref > r.ref; }
    ///\e
    bool operator!=(const RefPtr &r) const { return this->ref != r.ref; }
    
    ///\e
    operator bool() const { return ref; }

    ///\e
    const RefPtr &borrow(const T* &p) { 
      lock();
      if(ref==p) {
	if(prev) prev->next=next;
	if(next) next->prev=prev;
      }
      else release();
      ref=p;
      next=prev=this;
      unlock();
      return *this;
    }
    
    ///\e
    const RefPtr &borrow() { 
      lock();
      if(prev) prev->next=next;
      if(next) next->prev=prev;
      next=prev=this;
      unlock();
      return *this;
    }
    
  };  //END OF CLASS REFPTR
  
} //END OF NAMESPACE LEMON

#endif
