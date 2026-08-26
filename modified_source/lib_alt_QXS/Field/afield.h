/*!
      @file    afield.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef QXS_AFIELD_INCLUDED
#define QXS_AFIELD_INCLUDED

#include <cstddef>
#include <cstdio>
#include <vector>
#include <string>

#include  "lib_alt/Field/afield_base.h"  // primary template

#include "bridge_defs.h"
#include "bridge_complex.h"
#include "complexTraits.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Field/field.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/alignment_size_qxs.h"
#include "lib_alt_QXS/inline/define_params.h"
#include "aligned_allocator_impl.h"

template<typename REALTYPE>
class AField<REALTYPE, QXS>
{
 public:
  typedef REALTYPE                                      real_t;
  typedef typename ComplexTraits<REALTYPE>::complex_t   complex_t;
  typedef Element_type::type                            element_type; // see bridge_defs.h
  static const Impl IMPL = QXS;
  template<typename T>
  using aligned_allocator = aligned_allocator_offset_impl<T, alignment_size<IMPL>(), alignment_size<IMPL>()>;
  static const std::string class_name;

 protected:
  int m_nin;
  int m_nvol;
  int m_nex;
  element_type m_element_type;
  Bridge::VerboseLevel m_vl;

  std::size_t m_nsize;
  std::size_t m_nsize_off; // size with offset
  std::size_t m_nsizev;    // m_nsize/VLEN
  std::size_t m_offset;
  std::vector<real_t, aligned_allocator<real_t> > m_field;

  int m_size_unit;                                          //!< unit of reduction
  int m_size_rem;                                           //!< nsize % m_size_unit
  std::size_t m_size_quot;                                  //!< nsize/m_size_unit
  std::vector<real_t, aligned_allocator<real_t> > m_reduce; //!< array for reduction

 public:
  //! constructor without argument
  AField() { init(0, 0, 0, Element_type::COMPLEX); }

  //! constructor
  AField(const int nin, const int nvol, const int nex,
         const element_type cmpl = Element_type::COMPLEX)
  { init(nin, nvol, nex, cmpl); }

  //! copy constructor
  AField(const AField<real_t, QXS>& w)
  {
    init(w.nin(), w.nvol(), w.nex(), w.field_element_type());
    copy(w);
  }

  //! copy constructor
  AField(const Field& w)
  {
    init(w.nin(), w.nvol(), w.nex(), w.field_element_type());
    copy(w);
  }

  //! destructor
  ~AField() { tidyup(); }

 private:
  // void init(const size_t nsize);
  void init(const int nin, const int nvol, const int nex,
            const element_type cmpl);

  void tidyup();

 public:

  // resetting object.
  void reset(const int nin, const int nvol, const int nex,
             const element_type cmpl = Element_type::COMPLEX)
  {
    if (check_size(nin, nvol, nex) && (m_element_type == cmpl)) return;

    init(nin, nvol, nex, cmpl);
  }

  //! returning size of inner (on site) d.o.f.
  int nin() const { return m_nin; }

  //! returning size of site d.o.f.
  int nvol() const { return m_nvol; }

  //! returning size of extra d.o.f.
  int nex() const { return m_nex; }

  //! returning element_type (real or complex).
  element_type field_element_type() const { return m_element_type; }

  //! checking size parameters.
  bool check_size(const int nin, const int nvol, const int nex) const
  {
    bool chk = true;
    if ((m_nin != nin) || (m_nvol != nvol) || (m_nex != nex)) chk = false;
    return chk;
  }

  bool check_size(const AField<REALTYPE, QXS>& w) const
  {
    bool chk = true;
    if ((m_nin != w.nin()) || (m_nvol != w.nvol()) || (m_nex != w.nex()))
      chk = false;
    return chk;
  }

  //! reference of data element.  to be discarded.
  inline real_t& e(const int index)
  { return m_field[m_offset + index]; }

  //! return the array size
  inline int size(void) const { return m_nsize; }

  //! return the array size
  inline int ntot(void) const { return m_nsize; }

  real_t *ptr(int i) { return &m_field[m_offset + i]; }

  const real_t *ptr(int i) const { return &m_field[m_offset + i]; }

  complex_t *ptr_complex(int i)
  { return (complex_t *)&m_field[m_offset + 2 * i]; }

  const complex_t *ptr_complex(int i) const
  { return (const complex_t *)&m_field[m_offset + 2 * i]; }

  //! reference of data element
  real_t cmp(const int index) const
  { return m_field[m_offset + index]; }

  void set(const int index, const real_t a)
  { m_field[m_offset + index] = a; }

  void set(const real_t a);

  void copy(const Field& w);

  void copy(const AField<real_t, QXS>& w);

  void copy(const int ex,
            const AField<real_t, QXS>& w, const int ex_w);

  //! this is to be discarded.

  /*
  void copy(const int index,
            const AField<real_t, QXS>& w, const int index_w, const int size);
  */

  void axpy(const real_t, const AField<real_t, QXS>&);


  //! y = y + a*x and return ||y||^2 in one traversal.
  real_t axpy_norm2(const real_t, const AField<real_t, QXS>&);

  //! z = x + a*y in one traversal (destination z is *this).
  void copy_axpy(const AField<real_t, QXS>&, const real_t,
                 const AField<real_t, QXS>&);
  

  void axpy(const int ex, const real_t a,
            const AField<real_t, QXS>& w, const int ex_w);

  void axpy(const int ex, const real_t ar, const real_t ai,
            const AField<real_t, QXS>& w, const int ex_w);

  void axpy(const real_t ar, const real_t ai,
            const AField<real_t, QXS>& w);

  void aypx(const int ex, const real_t ar, const real_t ai,
            const AField<real_t, QXS>& w, const int ex_w);

  void aypx(const real_t ar, const real_t ai,
            const AField<real_t, QXS>& w);


  //! this is to be discarded.
  void axpy(const int index, const real_t a,
            const AField<real_t, QXS>& w, const int index_w,
            const int size);

  void aypx(const real_t, const AField<real_t, QXS>&);

  void scal(const real_t);

  void scal(const real_t, const real_t);

  real_t dot(const AField<real_t, QXS>&) const;

  void dotc(real_t&, real_t&, const AField<real_t, QXS>&) const;

  real_t norm2(void) const;

  complex_t dotc_and_norm2(real_t&, real_t&, const AField<real_t, QXS>&) const;
};

#endif
