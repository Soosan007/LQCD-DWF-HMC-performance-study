/*!
        @file    afield-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2023-04-04 15:28:35 #$
        @version $LastChangedRevision: 2507 $
*/

#ifndef QXS_AFIELD_INC_INCLUDED
#define QXS_AFIELD_INC_INCLUDED

#include <cstdlib>

#include "complexTraits.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"


//====================================================================
template<typename REALTYPE>
void copy(AField<REALTYPE, QXS>& v, const AField<REALTYPE, QXS>& w)
{
  v.copy(w);
}


//====================================================================
template<typename REALTYPE>
void copy(AField<REALTYPE, QXS>& v, const int ex,
          const AField<REALTYPE, QXS>& w, const int ex_w)
{
  v.copy(ex, w, ex_w);
}


//====================================================================
template<typename REALTYPE>
void axpy(AField<REALTYPE, QXS>& v, const int exv,
          const typename AField<REALTYPE, QXS>::real_t a,
          const AField<REALTYPE, QXS>& w, const int exw)
{
  v.axpy(exv, a, w, exw);
}


//====================================================================
template<typename REALTYPE>
void axpy(AField<REALTYPE, QXS>& v,
          const typename AField<REALTYPE, QXS>::real_t a,
          const AField<REALTYPE, QXS>& w)
{
  v.axpy(a, w);
}

// axpy norm2メモリ削減のために融合させたもの
//====================================================================
template<typename REALTYPE>
REALTYPE axpy_norm2(AField<REALTYPE, QXS>& v,
                    const typename AField<REALTYPE, QXS>::real_t a,
                    const AField<REALTYPE, QXS>& w)
{
  return v.axpy_norm2(a, w);
}


// z = x + a*y を1回のメモリ走査で計算する。
//====================================================================
template<typename REALTYPE>
void copy_axpy(AField<REALTYPE, QXS>& z,
               const AField<REALTYPE, QXS>& x,
               const typename AField<REALTYPE, QXS>::real_t a,
               const AField<REALTYPE, QXS>& y)
{
  z.copy_axpy(x, a, y);
}


//====================================================================
template<typename REALTYPE>
void axpy(AField<REALTYPE, QXS>& v, const typename AField<REALTYPE, QXS>::complex_t a, const AField<REALTYPE, QXS>& w)
{
  v.axpy(real(a), imag(a), w);
}


//====================================================================
template<typename REALTYPE>
void axpy(AField<REALTYPE, QXS>& v, const int ex,
          const typename AField<REALTYPE, QXS>::complex_t a,
          const AField<REALTYPE, QXS>& w, const int ex_w)
{
  v.axpy(ex, real(a), imag(a), w, ex_w);
}


//====================================================================
template<typename REALTYPE>
void aypx(const typename AField<REALTYPE, QXS>::real_t a,
          AField<REALTYPE, QXS>& v, const AField<REALTYPE, QXS>& w)
{
  v.aypx(a, w);
}


//====================================================================
template<typename REALTYPE>
void aypx(const typename AField<REALTYPE, QXS>::complex_t a, AField<REALTYPE, QXS>& v, const AField<REALTYPE, QXS>& w)
{
  v.aypx(real(a), imag(a), w);
}


//====================================================================
template<typename REALTYPE>
void scal(AField<REALTYPE, QXS>& v, REALTYPE a)
{
  v.scal(a);
}


//====================================================================
template<typename REALTYPE>
void scal(AField<REALTYPE, QXS>& v, const typename AField<REALTYPE, QXS>::complex_t a)
{
  v.scal(real(a), imag(a));
}


//====================================================================
template<typename REALTYPE>
REALTYPE dot(AField<REALTYPE, QXS>& v, AField<REALTYPE, QXS>& w)
{
  return v.dot(w);
}


//====================================================================
template<typename REALTYPE>
typename AField<REALTYPE, QXS>::complex_t dotc(const AField<REALTYPE, QXS>& v, const AField<REALTYPE, QXS>& w)
{
  REALTYPE vw_r, vw_i;
  v.dotc(vw_r, vw_i, w);
  return typename AField<REALTYPE, QXS>::complex_t(vw_r, vw_i);
}


//====================================================================
template<class INDEX, class AFIELD>
void convert(INDEX& index, AFIELD& v, const Field& w)
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  typename AFIELD::real_t *v2 = const_cast<AFIELD *>(&v)->ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int in = 0; in < Nin; ++in) {
        int iw = in + Nin * (site + Nvol * ex);
        int iv = index.idx(in, Nin, site, ex);
        v2[iv] = w.cmp(iw);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<class INDEX, class AFIELD>
void reverse(INDEX& index, Field& v, const AFIELD& w)
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int in = 0; in < Nin; ++in) {
        int iw = in + Nin * (site + Nvol * ex);
        int iv = index.idx(in, Nin, site, ex);
        v.set(iv, double(w.cmp(iw)));
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<class INDEX, class FIELD>
void convert_spinor(INDEX& index, FIELD& v, const Field& w)
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  assert(Nin == 2 * Nc * Nd);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int id = 0; id < Nd; ++id) {
        for (int ic = 0; ic < Nc; ++ic) {
          int iwr = 2 * (ic + Nc * id) + Nin * (site + Nvol * ex);
          int iwi = 1 + 2 * (ic + Nc * id) + Nin * (site + Nvol * ex);
          int ivr = index.idx_SPr(ic, id, site, ex);
          int ivi = index.idx_SPi(ic, id, site, ex);
          v.e(ivr) = w.cmp(iwr);
          v.e(ivi) = w.cmp(iwi);
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<class INDEX, class FIELD>
void convert_gauge(INDEX& index, FIELD& v, const Field& w)
{
  int Nin  = w.nin();
  int Nvol = w.nvol();
  int Nex  = w.nex();
  assert(v.check_size(Nin, Nvol, Nex));

  int Nc = CommonParameters::Nc();
  assert(Nin == 2 * Nc * Nc);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int ic2 = 0; ic2 < Nc; ++ic2) {
        for (int ic1 = 0; ic1 < Nc; ++ic1) {
          int iwr = 2 * (ic1 + Nc * ic2) + Nin * (site + Nvol * ex);
          int iwi = 1 + 2 * (ic1 + Nc * ic2) + Nin * (site + Nvol * ex);
          int ivr = index.idx_Gr(ic1, ic2, site, ex);
          int ivi = index.idx_Gi(ic1, ic2, site, ex);
          v.e(ivr) = w.cmp(iwr);
          v.e(ivi) = w.cmp(iwi);
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<class INDEX2, class FIELD2, class INDEX1, class FIELD1>
void convert(INDEX2& index2, FIELD2& v2,
             const INDEX1& index1, const FIELD1& v1)
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  int ith, nth, is, ns;

  if ((sizeof(typename FIELD1::real_t) == 4) || (sizeof(typename FIELD2::real_t) == 4)) {
    set_threadtask(ith, nth, is, ns, Nvol / VLENS);
    for (int ex = 0; ex < Nex; ++ex) {
      for (int vsite = is; vsite < ns; ++vsite) {
        for (int in = 0; in < Nin; ++in) {
          for (int vin = 0; vin < VLENS; ++vin) {
            int site = VLENS * vsite + vin;
            int iv1  = index1.idx(in, Nin, site, ex);
            int iv2  = index2.idx(in, Nin, site, ex);
            v2.e(iv2) = v1.cmp(iv1);
          }
        }
      }
    }
  } else {
    set_threadtask(ith, nth, is, ns, Nvol / VLEND);
    for (int ex = 0; ex < Nex; ++ex) {
      for (int vsite = is; vsite < ns; ++vsite) {
        for (int in = 0; in < Nin; ++in) {
          for (int vin = 0; vin < VLEND; ++vin) {
            int site = VLEND * vsite + vin;
            int iv1  = index1.idx(in, Nin, site, ex);
            int iv2  = index2.idx(in, Nin, site, ex);
            v2.e(iv2) = v1.cmp(iv1);
          }
        }
      }
    }
  }


#pragma omp barrier
}


//====================================================================
template<class INDEX2, class FIELD2, class INDEX1, class FIELD1>
void convert_h(INDEX2& index2, FIELD2& v2,
               const INDEX1& index1, const FIELD1& v1)
{
  int Nin  = v1.nin();
  int Nvol = v1.nvol();
  int Nex  = v1.nex();
  assert(v2.check_size(Nin, Nvol, Nex));

  int ith, nth, is, ns;
  if ((sizeof(typename FIELD1::real_t) == 4) || (sizeof(typename FIELD2::real_t) == 4)) {
    set_threadtask(ith, nth, is, ns, Nvol / VLENS);
    for (int ex = 0; ex < Nex; ++ex) {
      for (int vsite = is; vsite < ns; ++vsite) {
        for (int in = 0; in < Nin; ++in) {
          for (int vin = 0; vin < VLENS; ++vin) {
            int site = VLENS * vsite + vin;
            int iv1  = index1.idxh(in, Nin, site, ex);
            int iv2  = index2.idxh(in, Nin, site, ex);
            v2.e(iv2) = v1.cmp(iv1);
          }
        }
      }
    }
  } else {
    set_threadtask(ith, nth, is, ns, Nvol / VLEND);
    for (int ex = 0; ex < Nex; ++ex) {
      for (int vsite = is; vsite < ns; ++vsite) {
        for (int in = 0; in < Nin; ++in) {
          for (int vin = 0; vin < VLEND; ++vin) {
            int site = VLEND * vsite + vin;
            int iv1  = index1.idxh(in, Nin, site, ex);
            int iv2  = index2.idxh(in, Nin, site, ex);
            v2.e(iv2) = v1.cmp(iv1);
          }
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<class INDEX, class FIELD>
void reverse(INDEX& index, Field& v, FIELD& w)
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  w.check_size(Nin, Nvol, Nex);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int in = 0; in < Nin; ++in) {
        int iv = in + Nin * (site + Nvol * ex);
        int iw = index.idx(in, Nin, site, ex);
        v.set(iv, double(w.cmp(iw)));
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<class INDEX, class FIELD>
void reverse_spinor(INDEX& index, Field& v, FIELD& w)
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  w.check_size(Nin, Nvol, Nex);

  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  assert(Nin == 2 * Nc * Nd);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int id = 0; id < Nd; ++id) {
        for (int ic = 0; ic < Nc; ++ic) {
          int ivr = 2 * (ic + Nc * id) + Nin * (site + Nvol * ex);
          int ivi = 1 + 2 * (ic + Nc * id) + Nin * (site + Nvol * ex);
          int iwr = index.idx_SPr(ic, id, site, ex);
          int iwi = index.idx_SPi(ic, id, site, ex);
          v.set(ivr, double(w.cmp(iwr)));
          v.set(ivi, double(w.cmp(iwi)));
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<class INDEX, class FIELD>
void reverse_gauge(INDEX& index, Field& v, FIELD& w)
{
  int Nin  = v.nin();
  int Nvol = v.nvol();
  int Nex  = v.nex();
  w.check_size(Nin, Nvol, Nex);

  int Nc = CommonParameters::Nc();
  assert(Nin == 2 * Nc * Nc);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int ic2 = 0; ic2 < Nc; ++ic2) {
        for (int ic1 = 0; ic1 < Nc; ++ic1) {
          int ivr = 2 * (ic1 + Nc * ic2) + Nin * (site + Nvol * ex);
          int ivi = 1 + 2 * (ic1 + Nc * ic2) + Nin * (site + Nvol * ex);
          int iwr = index.idx_Gr(ic1, ic2, site, ex);
          int iwi = index.idx_Gi(ic1, ic2, site, ex);
          v.set(ivr, double(w.cmp(iwr)));
          v.set(ivi, double(w.cmp(iwi)));
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
REALTYPE norm2(const AField<REALTYPE, QXS>& v)
{
  AField<REALTYPE, QXS> *vp = const_cast<AField<REALTYPE, QXS> *>(&v);
  return vp->norm2();
  //  return v.norm2();
}


//====================================================================
template<typename REALTYPE>
void dotc_and_norm2(typename AField<REALTYPE, QXS>::complex_t& dotc,
                    REALTYPE& v_norm2, REALTYPE& w_norm2,
                    const AField<REALTYPE, QXS>& v,
                    const AField<REALTYPE, QXS>& w)
{
  // returns <v|w>, <v|v>, <w|w>
  assert(v.nex() == w.nex());
  assert(v.nvol() == w.nvol());
  assert(v.nin() == w.nin());

  dotc = v.dotc_and_norm2(v_norm2, w_norm2, w);
}


#define AFIELD_HAS_DOTC_AND_NORM2

//============================================================END=====
#endif
