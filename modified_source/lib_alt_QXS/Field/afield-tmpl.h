/*!
      @file    afield-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#include <cassert>

#include "lib_alt_QXS/inline/afield_th-inc.h"

template<typename REALTYPE>
const std::string AField<REALTYPE, QXS>::class_name = "AField";

//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::init(const int nin, const int nvol, const int nex,
                                 const element_type cmpl)
{
  m_nin          = nin;
  m_nvol         = nvol;
  m_nex          = nex;
  m_element_type = cmpl;

  m_offset = 0;

  if (m_nvol % VLEN != 0) {
    vout.crucial("%s: bad nvol (too small?), must be a multiple of VLEN %d (given nvol=%d)\n",
                 class_name.c_str(), VLEN, m_nvol);
    exit(EXIT_FAILURE);
  }

  m_nsize     = m_nin * m_nvol * m_nex;
  m_nsize_off = m_nsize + m_offset + m_offset;
  // offset is set both before and after the net data
  m_nsizev = m_nsize / VLEN;
  if (m_nsize_off != m_field.size()) {
    m_field.resize(m_nsize_off, 0.0f);
    vout.detailed("%s: data resized: nin = %d, nvol = %d, nex = %d\n",
                  class_name.c_str(), m_nin, m_nvol, m_nex);
  }
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::tidyup()
{
  // currently do nothing.
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::set(const REALTYPE a)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, size());
#pragma omp barrier
  for (int i = is; i < ns; ++i) {
    m_field[m_offset + i] = a;
  }
#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::copy(const Field& w)
{
  assert(check_size(w));
#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, size());

  for (int i = is; i < ns; ++i) {
    m_field[m_offset + i] = REALTYPE(w.cmp(i));
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::copy(const AField<REALTYPE, QXS>& w)
{
  assert(check_size(w));

#pragma omp barrier

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = ptr(0);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  real_t *vp = ptr(0);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nsizev);

  svbool_t pg = set_predicate();
  for (int i = is; i < ns; ++i) {
    svreal_t vt;
    load_vec(pg, vt, &wp[VLEN * i]);
    save_vec(pg, &vp[VLEN * i], vt);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::copy(const int ex,
                                 const AField<REALTYPE, QXS>& w, const int ex_w)
{
  assert(nin() == w.nin());
  assert(nvol() == w.nvol());

#pragma omp barrier

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = ptr(0);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  real_t *vp = ptr(0);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

  int sizev = m_nin * (m_nvol / VLEN);
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

  svbool_t pg = set_predicate();
  for (int i = is; i < ns; ++i) {
    svreal_t vt;
    load_vec(pg, vt, &wp[VLEN * (i + sizev * ex_w)]);
    save_vec(pg, &vp[VLEN * (i + sizev * ex)], vt);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::axpy(const REALTYPE a, const AField<REALTYPE, QXS>& w)
{
  assert(check_size(w));

#pragma omp barrier

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = ptr(0);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  real_t *vp = ptr(0);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nsizev);

  svbool_t pg = set_predicate();
  for (int i = is; i < ns; ++i) {
    svreal_t vt, wt;
    load_vec(pg, wt, &wp[VLEN * i]);
    load_vec(pg, vt, &vp[VLEN * i]);
    axpy_vec(pg, vt, a, wt);
    save_vec(pg, &vp[VLEN * i], vt);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::copy_axpy(
    const AField<REALTYPE, QXS>& x, const REALTYPE a,
    const AField<REALTYPE, QXS>& y)
{
  assert(check_size(x));
  assert(check_size(y));

#pragma omp barrier

#ifdef USE_QXS_ACLE
  real_t *__restrict__ zp = ptr(0);
  real_t *__restrict__ xp = const_cast<AField<REALTYPE, QXS> *>(&x)->ptr(0);
  real_t *__restrict__ yp = const_cast<AField<REALTYPE, QXS> *>(&y)->ptr(0);
#else
  real_t *zp = ptr(0);
  real_t *xp = const_cast<AField<REALTYPE, QXS> *>(&x)->ptr(0);
  real_t *yp = const_cast<AField<REALTYPE, QXS> *>(&y)->ptr(0);
#endif
  assert(zp != xp);
  assert(zp != yp);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nsizev);

  svbool_t pg = set_predicate();
  for (int i = is; i < ns; ++i) {
    svreal_t zt, yt;
    load_vec(pg, zt, &xp[VLEN * i]);
    load_vec(pg, yt, &yp[VLEN * i]);
    axpy_vec(pg, zt, a, yt);
    save_vec(pg, &zp[VLEN * i], zt);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE, QXS>::axpy_norm2(
    const REALTYPE a, const AField<REALTYPE, QXS>& w)
{
  assert(check_size(w));

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = ptr(0);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  real_t *vp = ptr(0);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

#pragma omp barrier

  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

  int Nouter = (Nstv > Nin) ? Nstv : Nin;
  int Ninner = (Nstv * Nin) / Nouter;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nouter);

  svbool_t pg = set_predicate();
  svreal_t yt;
  set_vec(pg, yt, 0.0);

  for (int ex = 0; ex < Nex; ++ex) {
    svreal_t tmp_yt;
    set_vec(pg, tmp_yt, 0.0);
    for (int i = is; i < ns; ++i) {
      real_t *v = &vp[VLEN * Ninner * (i + Nouter * ex)];
      real_t *x = &wp[VLEN * Ninner * (i + Nouter * ex)];
      svreal_t xt;
      set_vec(pg, xt, 0.0);
      for (int in = 0; in < Ninner; ++in) {
        svreal_t vt, wt;
        load_vec(pg, vt, &v[VLEN * in]);
        load_vec(pg, wt, &x[VLEN * in]);
        axpy_vec(pg, vt, a, wt);
        save_vec(pg, &v[VLEN * in], vt);
        add_norm2_vec(pg, xt, vt);
      }
      add_vec(pg, tmp_yt, xt);
    }
    add_vec(pg, yt, tmp_yt);
  }

  real_t result;
  reduce_vec(pg, result, yt);

#pragma omp barrier

  ThreadManager::reduce_sum_global(result, ith, nth);

  return result;
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::axpy(const int ex, const REALTYPE a,
                                 const AField<REALTYPE, QXS>& w,
                                 const int ex_w)
{
  assert(nin() == w.nin());
  assert(nvol() == w.nvol());

#pragma omp barrier

  int sizev = m_nin * (m_nvol / VLEN);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, sizev);

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = ptr(VLEN * sizev * ex);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(VLEN * sizev * ex_w);
#else
  real_t *vp = ptr(VLEN * sizev * ex);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(VLEN * sizev * ex_w);
#endif
  assert(vp != wp);

  svbool_t pg = set_predicate();
  for (int i = is; i < ns; ++i) {
    svreal_t vt, wt;
    load_vec(pg, wt, &wp[VLEN * i]);
    load_vec(pg, vt, &vp[VLEN * i]);
    axpy_vec(pg, vt, a, wt);
    save_vec(pg, &vp[VLEN * i], vt);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::axpy(const REALTYPE ar, const REALTYPE ai,
                                 const AField<REALTYPE, QXS>& w)
{
  assert(check_size(w));

#pragma omp barrier

#ifdef USE_QXS_ACLE
  REALTYPE *__restrict__ vp = this->ptr(0);
  REALTYPE *__restrict__ wp
    = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  REALTYPE *vp = this->ptr(0);
  REALTYPE *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif

  svbool_t pg = set_predicate();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nsizev / 2);

  for (int i = is; i < ns; ++i) {
    int      inr = VLEN * (2 * i);
    int      ini = VLEN * (2 * i + 1);
    svreal_t wtr, wti, vtr, vti;
    load_vec(pg, wtr, &wp[inr]);
    load_vec(pg, wti, &wp[ini]);
    load_vec(pg, vtr, &vp[inr]);
    load_vec(pg, vti, &vp[ini]);
    axpy_vec(pg, vtr, ar, wtr);
    axpy_vec(pg, vtr, -ai, wti);
    axpy_vec(pg, vti, ar, wti);
    axpy_vec(pg, vti, ai, wtr);
    save_vec(pg, &vp[inr], vtr);
    save_vec(pg, &vp[ini], vti);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::axpy(const int ex,
                                 const REALTYPE atr, const REALTYPE ati,
                                 const AField<REALTYPE, QXS>& w, const int ex_w)
{
  assert(nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

  //  set_threadtask(ith, nth, is, ns, Nstv);

#pragma omp barrier

#ifdef USE_QXS_ACLE
  REALTYPE *__restrict__ vp = this->ptr(0) + VLEN * Nin * Nstv * ex;
  REALTYPE *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0) + VLEN * Nin * Nstv * ex_w;
#else
  REALTYPE *vp = this->ptr(0) + VLEN * Nin * Nstv * ex;
  REALTYPE *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0)
                 + VLEN * Nin * Nstv * ex_w;
#endif
  assert(vp != wp);

  svbool_t pg      = set_predicate();
  int      Nstvin2 = Nstv * Nin / 2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstvin2);

  for (int i2 = 2 * is * VLEN; i2 < 2 * ns * VLEN; i2 += 2 * VLEN) {
    int      inr = i2;
    int      ini = i2 + VLEN;
    svreal_t wtr, wti, vtr, vti;
    load_vec(pg, wtr, &wp[inr]);
    load_vec(pg, vtr, &vp[inr]);
    load_vec(pg, wti, &wp[ini]);
    load_vec(pg, vti, &vp[ini]);
    axpy_vec(pg, vtr, atr, wtr);
    axpy_vec(pg, vtr, -ati, wti);
    axpy_vec(pg, vti, atr, wti);
    axpy_vec(pg, vti, ati, wtr);
    save_vec(pg, &vp[inr], vtr);
    save_vec(pg, &vp[ini], vti);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::aypx(const REALTYPE a,
                                 const AField<REALTYPE, QXS>& w)
{
  assert(check_size(w));

#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nsizev);

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = ptr(0);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  real_t *vp = ptr(0);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

  svbool_t pg = set_predicate();
  for (int i = is; i < ns; ++i) {
    svreal_t vt, wt;
    load_vec(pg, wt, &wp[VLEN * i]);
    load_vec(pg, vt, &vp[VLEN * i]);
    aypx_vec(pg, a, vt, wt);
    save_vec(pg, &vp[VLEN * i], vt);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::aypx(const REALTYPE ar, const REALTYPE ai,
                                 const AField<REALTYPE, QXS>& w)
{
  assert(check_size(w));

#pragma omp barrier

#ifdef USE_QXS_ACLE
  REALTYPE *__restrict__ vp = this->ptr(0);
  REALTYPE *__restrict__ wp
    = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  REALTYPE *vp = this->ptr(0);
  REALTYPE *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

  svbool_t pg = set_predicate();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_nsizev / 2);

  for (int i = is; i < ns; ++i) {
    int      inr = VLEN * (2 * i);
    int      ini = VLEN * (2 * i + 1);
    svreal_t wtr, wti, vtr, vti;
    load_vec(pg, wtr, &wp[inr]);
    load_vec(pg, wti, &wp[ini]);
    load_vec(pg, vtr, &vp[inr]);
    load_vec(pg, vti, &vp[ini]);
    axpy_vec(pg, wtr, ar, vtr);
    axpy_vec(pg, wti, ar, vti);
    axpy_vec(pg, wtr, -ai, vti);
    axpy_vec(pg, wti, ai, vtr);
    save_vec(pg, &vp[inr], wtr);
    save_vec(pg, &vp[ini], wti);
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::aypx(const int ex,
                                 const REALTYPE atr, const REALTYPE ati,
                                 const AField<REALTYPE, QXS>& w,
                                 const int ex_w)
{
  assert(nin() == w.nin());
  assert(nvol() == w.nvol());
  assert(ex < nex());
  assert(ex_w < w.nex());

  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

#pragma omp barrier

#ifdef USE_QXS_ACLE
  REALTYPE *__restrict__ vp = this->ptr(0) + VLEN * Nin * Nstv * ex;
  REALTYPE *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0) + VLEN * Nin * Nstv * ex_w;
#else
  REALTYPE *vp = this->ptr(0) + VLEN * Nin * Nstv * ex;
  REALTYPE *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0)
                 + VLEN * Nin * Nstv * ex_w;
#endif
  assert(vp != wp);

  svbool_t pg = set_predicate();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv * Nin / 2);

  for (int i2 = 2 * is * VLEN; i2 < 2 * ns * VLEN; i2 += 2 * VLEN) {
    int      inr = i2;
    int      ini = i2 + VLEN;
    svreal_t wtr, wti, vtr, vti;
    load_vec(pg, wtr, &wp[inr]);
    load_vec(pg, wti, &wp[ini]);
    load_vec(pg, vtr, &vp[inr]);
    load_vec(pg, vti, &vp[ini]);
    axpy_vec(pg, wtr, atr, vtr);
    axpy_vec(pg, wti, atr, vti);
    axpy_vec(pg, wtr, -ati, vti);
    axpy_vec(pg, wti, ati, vtr);
    save_vec(pg, &vp[inr], wtr);
    save_vec(pg, &vp[ini], wti);
  }
#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::scal(const REALTYPE a)
{
  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

  int ith, nth, is, ns;

  real_t *vp = ptr(0);

#pragma omp barrier

  svbool_t pg = set_predicate();
  set_threadtask(ith, nth, is, ns, Nstv * Nin);
  for (int ex = 0; ex < Nex; ++ex) {
    real_t *v = vp + VLEN * Nstv * Nin * ex;
    for (int iv = is * VLEN; iv < ns * VLEN; iv += VLEN) {
      svreal_t vt;
      load_vec(pg, vt, &v[iv]);
      scal_vec(pg, vt, a);
      save_vec(pg, &v[iv], vt);
    }
  }
#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::scal(const REALTYPE ar, const REALTYPE ai)
{
  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

#pragma omp barrier

#ifdef USE_QXS_ACLE
  REALTYPE *__restrict__ vp = this->ptr(0);
#else
  REALTYPE *vp = this->ptr(0);
#endif

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv * Nin / 2);
  svbool_t pg = set_predicate();

  for (int ex = 0; ex < Nex; ++ex) {
    real_t *v = vp + VLEN * Nstv * Nin * ex;
    for (int i = is; i < ns; ++i) {
      int      inr = VLEN * (2 * i);
      int      ini = VLEN * (2 * i + 1);
      svreal_t vtr, vti;
      load_vec(pg, vtr, &vp[inr]);
      load_vec(pg, vti, &vp[ini]);
      axpy_vec(pg, vtr, ar, vtr);
      axpy_vec(pg, vtr, -ai, vti);
      axpy_vec(pg, vti, ar, vti);
      axpy_vec(pg, vti, ai, vtr);
      save_vec(pg, &vp[inr], vtr);
      save_vec(pg, &vp[ini], vti);
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE, QXS>::dot(const AField<REALTYPE, QXS>& w) const
{
  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = const_cast<AField<REALTYPE, QXS> *>(this)->ptr(0);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  real_t *vp = const_cast<AField<REALTYPE, QXS> *>(this)->ptr(0);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

#pragma omp barrier

  svbool_t pg = set_predicate();

  svreal_t yt;
  set_vec(pg, yt, 0.0);

  for (int ex = 0; ex < Nex; ++ex) {
    svreal_t tmp_yt;
    set_vec(pg, tmp_yt, 0.0);
    for (int site = is; site < ns; ++site) {
      real_t   *v = &vp[VLEN * Nin * (site + Nstv * ex)];
      real_t   *w = &wp[VLEN * Nin * (site + Nstv * ex)];
      svreal_t xt, vt, wt;
      set_vec(pg, xt, 0.0);
      for (int in = 0; in < Nin; ++in) {
        load_vec(pg, wt, &w[VLEN * in]);
        load_vec(pg, vt, &v[VLEN * in]);
        add_dot_vec(pg, xt, vt, wt);
      }
      add_vec(pg, tmp_yt, xt);
    }
    add_vec(pg, yt, tmp_yt);
  }

  real_t a;
  reduce_vec(pg, a, yt);

#pragma omp barrier

  ThreadManager::reduce_sum_global(a, ith, nth);

  return a;
}


//====================================================================
template<typename REALTYPE>
void AField<REALTYPE, QXS>::dotc(REALTYPE& atr, REALTYPE& ati,
                                 const AField<REALTYPE, QXS>& w) const
{
  assert(check_size(w));

  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;
  int Nin2 = Nin / 2;

  int ith, nth, is, ns;

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = const_cast<AField<REALTYPE, QXS> *>(this)->ptr(0);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  real_t *vp = const_cast<AField<REALTYPE, QXS> *>(this)->ptr(0);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

#pragma omp barrier


  svbool_t pg = set_predicate();

  int Nouter  = (Nstv > Nin2) ? Nstv : Nin2;
  int Ninner2 = (Nstv * Nin2) / Nouter;
  int Ninner  = 2 * Ninner2;
  set_threadtask(ith, nth, is, ns, Nouter); // for better load balance
  svreal_t ytr, yti;
  set_vec(pg, ytr, 0.0);
  set_vec(pg, yti, 0.0);

  for (int ex = 0; ex < Nex; ++ex) {
    svreal_t tmp_ytr, tmp_yti;
    set_vec(pg, tmp_ytr, 0.0);
    set_vec(pg, tmp_yti, 0.0);
    for (int i = is; i < ns; ++i) {
      real_t   *v = &vp[VLEN * Ninner * (i + Nouter * ex)];
      real_t   *w = &wp[VLEN * Ninner * (i + Nouter * ex)];
      svreal_t xtr, xti, vtr, vti, wtr, wti;
      set_vec(pg, xtr, 0.0);
      set_vec(pg, xti, 0.0);
      for (int in = 0; in < Ninner2; ++in) {
        int inr = 2 * in;
        int ini = 2 * in + 1;
        load_vec(pg, wtr, &w[VLEN * inr]);
        load_vec(pg, wti, &w[VLEN * ini]);
        load_vec(pg, vtr, &v[VLEN * inr]);
        load_vec(pg, vti, &v[VLEN * ini]);
        add_dot_vec(pg, xtr, vtr, wtr);
        add_dot_vec(pg, xtr, vti, wti);
        add_dot_vec(pg, xti, vtr, wti);
        sub_dot_vec(pg, xti, vti, wtr);
      }
      add_vec(pg, tmp_ytr, xtr);
      add_vec(pg, tmp_yti, xti);
    }

    add_vec(pg, ytr, tmp_ytr);
    add_vec(pg, yti, tmp_yti);
  }

  reduce_vec(pg, atr, ytr);
  reduce_vec(pg, ati, yti);

#pragma omp barrier
  real_t sum[2] = { atr, ati };
  ThreadManager::reduce_sum_global(sum, 2, ith, nth);
  atr = sum[0];
  ati = sum[1];
}


//====================================================================
template<typename REALTYPE>
REALTYPE AField<REALTYPE, QXS>::norm2() const
{
  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;

  real_t *vp = const_cast<AField<REALTYPE, QXS> *>(this)->ptr(0);

#pragma omp barrier

  int Nouter = (Nstv > Nin) ? Nstv : Nin;
  int Ninner = (Nstv * Nin) / Nouter;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nouter); // for better load balance

  svbool_t pg = set_predicate();

  svreal_t yt;
  set_vec(pg, yt, 0.0);

  for (int ex = 0; ex < Nex; ++ex) {
    svreal_t tmp_yt;
    set_vec(pg, tmp_yt, 0.0);
    for (int i = is; i < ns; ++i) {
      real_t   *v = &vp[VLEN * Ninner * (i + Nouter * ex)];
      svreal_t xt, vt;
      set_vec(pg, xt, 0.0);
      for (int in = 0; in < Ninner; ++in) {
        load_vec(pg, vt, &v[VLEN * in]);
        add_norm2_vec(pg, xt, vt);
      }
      add_vec(pg, tmp_yt, xt);
    }
    add_vec(pg, yt, tmp_yt);
  }

  real_t a;
  reduce_vec(pg, a, yt);

#pragma omp barrier

  ThreadManager::reduce_sum_global(a, ith, nth);

  return a;
}


//====================================================================
template<typename REALTYPE>
typename AField<REALTYPE, QXS>::complex_t
AField<REALTYPE, QXS>::dotc_and_norm2(REALTYPE& norm2, REALTYPE& w_norm2,
                                      const AField<REALTYPE, QXS>& w) const
{
  if (m_element_type == Element_type::REAL) {
    vout.crucial("%s: dotc_and_norm2 for real field is called\n",
                 class_name.c_str());

    exit(EXIT_FAILURE);
  }
  assert(check_size(w));

  int Nin  = this->nin();
  int Nex  = this->nex();
  int Nstv = this->nvol() / VLEN;
  int Nin2 = Nin / 2;

#ifdef USE_QXS_ACLE
  real_t *__restrict__ vp = const_cast<AField<REALTYPE, QXS> *>(this)->ptr(0);
  real_t *__restrict__ wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#else
  real_t *vp = const_cast<AField<REALTYPE, QXS> *>(this)->ptr(0);
  real_t *wp = const_cast<AField<REALTYPE, QXS> *>(&w)->ptr(0);
#endif
  assert(vp != wp);

  real_t atr = 0.0;
  real_t ati = 0.0;
  norm2   = 0.0;
  w_norm2 = 0.0;
#pragma omp barrier

  int Nouter  = (Nstv > Nin2) ? Nstv : Nin2;
  int Ninner2 = (Nstv * Nin2) / Nouter;
  int Ninner  = 2 * Ninner2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nouter); // for better load balance

  svbool_t pg = set_predicate();

  svreal_t ytr, yti, w2, v2;
  set_vec(pg, ytr, 0.0);
  set_vec(pg, yti, 0.0);
  set_vec(pg, w2, 0.0);
  set_vec(pg, v2, 0.0);

  for (int ex = 0; ex < Nex; ++ex) {
    svreal_t tmp_ytr, tmp_yti, tmp_vt2, tmp_wt2;
    set_vec(pg, tmp_ytr, 0.0);
    set_vec(pg, tmp_yti, 0.0);
    set_vec(pg, tmp_vt2, 0.0);
    set_vec(pg, tmp_wt2, 0.0);
    for (int iout = is; iout < ns; ++iout) {
      real_t   *v = &vp[VLEN * Ninner * (iout + Nouter * ex)];
      real_t   *w = &wp[VLEN * Ninner * (iout + Nouter * ex)];
      svreal_t xtr, xti, vtr, vti, wtr, wti, vt2, wt2;
      set_vec(pg, xtr, 0.0);
      set_vec(pg, xti, 0.0);
      set_vec(pg, vt2, 0.0);
      set_vec(pg, wt2, 0.0);
      for (int in = 0; in < Ninner2; ++in) {
        int inr = 2 * in;
        int ini = 2 * in + 1;
        load_vec(pg, wtr, &w[VLEN * inr]);
        load_vec(pg, wti, &w[VLEN * ini]);
        load_vec(pg, vtr, &v[VLEN * inr]);
        load_vec(pg, vti, &v[VLEN * ini]);
        add_dot_vec(pg, xtr, vtr, wtr);
        add_dot_vec(pg, xtr, vti, wti);
        add_dot_vec(pg, xti, vtr, wti);
        sub_dot_vec(pg, xti, vti, wtr);
        add_norm2_vec(pg, vt2, vtr);
        add_norm2_vec(pg, vt2, vti);
        add_norm2_vec(pg, wt2, wtr);
        add_norm2_vec(pg, wt2, wti);
      }
      add_vec(pg, tmp_ytr, xtr);
      add_vec(pg, tmp_yti, xti);
      add_vec(pg, tmp_vt2, vt2);
      add_vec(pg, tmp_wt2, wt2);
    }
    add_vec(pg, ytr, tmp_ytr);
    add_vec(pg, yti, tmp_yti);
    add_vec(pg, v2, tmp_vt2);
    add_vec(pg, w2, tmp_wt2);
  }

  reduce_vec(pg, atr, ytr);
  reduce_vec(pg, ati, yti);
  reduce_vec(pg, w_norm2, w2);
  reduce_vec(pg, norm2, v2);
  real_t sum[4] = { atr, ati, w_norm2, norm2 };

#pragma omp barrier

  ThreadManager::reduce_sum_global(sum, 4, ith, nth);
  complex_t result = { sum[0], sum[1] };

  w_norm2 = sum[2];
  norm2   = sum[3];

  return result;
}


//============================================================END=====
