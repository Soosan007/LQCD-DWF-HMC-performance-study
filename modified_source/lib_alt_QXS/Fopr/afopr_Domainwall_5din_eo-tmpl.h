/*!
      @file    afopr_Domainwall_5din_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2668 $
*/

#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib/Tools/timer.h"

template<typename AFIELD>
const std::string AFopr_Domainwall_5din_eo<AFIELD>::class_name
  = "AFopr_Domainwall_5din_eo";
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

#ifdef USE_DAGD_DETAIL_PROFILE
  m_dagd_profile_calls = 0;
  for (int i = 0; i < 15; ++i) m_dagd_profile_sec[i] = 0.0;
#endif

  // switch for communication
  int req_comm = 1;  // set 1 if communication forced any time
  //int req_comm = 0;  // set 1 if communication forced any time

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());

  m_repr = "Dirac";  // now only the Dirac repr is available.

  std::string repr;
  if (!params.fetch_string("gamma_matrix_type", repr)) {
    if (repr != "Dirac") {
      vout.crucial("  Error at %s: unsupported gamma-matrix type: %s\n",
                   class_name.c_str(), repr.c_str());
      exit(EXIT_FAILURE);
    }
  }

  int Nc = CommonParameters::Nc();
  if (Nc != 3) {
    vout.crucial("%s: only applicable to Nc = 3\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int Nd = CommonParameters::Nd();
  m_Nvcd = 2 * Nc * Nd;
  m_Ndf  = 2 * Nc * Nc;

  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_Nx2  = m_Nx / 2;
  m_Nst2 = m_Nvol / 2;

  m_Nx2v  = m_Nx2 / VLENX;
  m_Nyv   = m_Ny / VLENY;
  m_Nst2v = m_Nst2 / VLEN;

  // condition check
  if (m_Nx % 2 != 0) {
    vout.crucial(m_vl, "%s: Nx must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (m_Nx % (2 * VLENX) != 0) {
    vout.crucial(m_vl, "%s: Nx must be mulriple of 2*VLENX.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  if (m_Ny % (VLENY) != 0) {
    vout.crucial(m_vl, "%s: Ny must be multiple of VLENY.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "  VLENX = %2d  Nx2v  = %d\n", VLENX, m_Nx2v);
  vout.general(m_vl, "  VLENY = %2d  Nyv   = %d\n", VLENY, m_Nyv);
  vout.general(m_vl, "  VLEN  = %2d  Nst2v = %d\n", VLEN, m_Nst2v);

  m_Leo.resize(m_Ny * m_Nz * m_Nt);

  int ipe3 = Communicator::ipe(3);
  int ipe2 = Communicator::ipe(2);
  int ipe1 = Communicator::ipe(1);
  for (int t = 0; t < m_Nt; ++t) {
    for (int z = 0; z < m_Nz; ++z) {
      for (int y = 0; y < m_Ny; ++y) {
        int t2 = ipe3 * m_Nt + t;
        int z2 = ipe2 * m_Nz + z;
        int y2 = ipe1 * m_Ny + y;
        m_Leo[y + m_Ny * (z + m_Nz * t)] = (y2 + z2 + t2) % 2;
      }
    }
  }

  m_Nsize[0] = m_Nx2v;
  m_Nsize[1] = m_Nyv;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  do_comm_any = 0;
  for (int mu = 0; mu < m_Ndim; ++mu) {
    do_comm[mu] = 1;
    if ((req_comm == 0) && (Communicator::npe(mu) == 1)) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  vout.increase_indent();

  m_Ns = 0; // temporary set
  set_parameters(params);

  vout.decrease_indent();

  m_Nbdsize.resize(m_Ndim);
  int Nbdin = (m_Nvcd / 2) * m_Ns;
  m_Nbdsize[0] = Nbdin * ((m_Ny * m_Nz * m_Nt + 1) / 2);
  m_Nbdsize[1] = Nbdin * m_Nx2 * m_Nz * m_Nt;
  m_Nbdsize[2] = Nbdin * m_Nx2 * m_Ny * m_Nt;
  m_Nbdsize[3] = Nbdin * m_Nx2 * m_Ny * m_Nz;

  setup_channels();

  // gauge configuration.
  m_Ueo.reset(m_Ndf, m_Nvol, m_Ndim);
  m_Ulex.reset(m_Ndf, m_Nvol, m_Ndim);

  vout.detailed(m_vl, "%s: initalization finished.\n",
                class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::tidyup()
{
#ifdef USE_DAGD_DETAIL_PROFILE
#ifdef USE_QXS_ACLE
  const char *profile_backend = "ACLE";
#else
  const char *profile_backend = "GENERAL";
#endif
  const char *profile_precision = (sizeof(real_t) == 4) ? "float" : "double";
  if (m_dagd_profile_calls > 0) {
#ifdef USE_DAGD_FUSED_COPY_AXPY
    const int profile_fused_copy_axpy = 1;
#else
    const int profile_fused_copy_axpy = 0;
#endif
    const double profile_calls = static_cast<double>(m_dagd_profile_calls);
    const double profile_dagd_avg = m_dagd_profile_sec[0] / profile_calls;
    const double profile_d_combine_avg = m_dagd_profile_sec[7] / profile_calls;
    const double profile_ddag_combine_avg = m_dagd_profile_sec[14] / profile_calls;
    const double profile_combine_total =
      m_dagd_profile_sec[7] + m_dagd_profile_sec[14];
    const double profile_combine_avg = profile_combine_total / profile_calls;
    const double profile_combine_fraction =
      (m_dagd_profile_sec[0] > 0.0)
      ? 100.0 * profile_combine_total / m_dagd_profile_sec[0] : 0.0;

    vout.general(m_vl,
      "DAGD_FUSION_PROFILE backend=%s rank=%d fopr=%p precision=%s "
      "fused_copy_axpy=%d timing_scheme=single_outer_block "
      "calls=%llu dagd_total_sec=%.9e "
      "dagd_avg_sec=%.9e d_combine_total_sec=%.9e d_combine_avg_sec=%.9e "
      "ddag_combine_total_sec=%.9e ddag_combine_avg_sec=%.9e "
      "combine_total_sec=%.9e combine_avg_sec=%.9e "
      "combine_fraction_percent=%.9e\n",
      profile_backend, Communicator::nodeid(), static_cast<void *>(this),
      profile_precision, profile_fused_copy_axpy, m_dagd_profile_calls,
      m_dagd_profile_sec[0], profile_dagd_avg,
      m_dagd_profile_sec[7], profile_d_combine_avg,
      m_dagd_profile_sec[14], profile_ddag_combine_avg,
      profile_combine_total, profile_combine_avg,
      profile_combine_fraction);

    vout.general(m_vl,
      "DAGD_DETAIL_PROFILE backend=%s rank=%d fopr=%p precision=%s calls=%llu "
      "total_sec=%.9e d_eo_odd_sec=%.9e l_inv_1_sec=%.9e u_inv_1_sec=%.9e "
      "d_eo_even_sec=%.9e l_inv_2_sec=%.9e u_inv_2_sec=%.9e d_combine_sec=%.9e "
      "udag_inv_1_sec=%.9e ldag_inv_1_sec=%.9e ddag_eo_odd_sec=%.9e "
      "udag_inv_2_sec=%.9e ldag_inv_2_sec=%.9e ddag_eo_even_sec=%.9e "
      "ddag_combine_sec=%.9e\n",
      profile_backend, Communicator::nodeid(), static_cast<void *>(this),
      profile_precision, m_dagd_profile_calls,
      m_dagd_profile_sec[0], m_dagd_profile_sec[1], m_dagd_profile_sec[2],
      m_dagd_profile_sec[3], m_dagd_profile_sec[4], m_dagd_profile_sec[5],
      m_dagd_profile_sec[6], m_dagd_profile_sec[7], m_dagd_profile_sec[8],
      m_dagd_profile_sec[9], m_dagd_profile_sec[10], m_dagd_profile_sec[11],
      m_dagd_profile_sec[12], m_dagd_profile_sec[13], m_dagd_profile_sec[14]);
  }
#endif
  // need to do nothing.
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    size_t Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu, 1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu, 1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void *buf_up = (void *)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu, 1, buf_up);
    void *buf_dn = (void *)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if (do_comm[mu] == 1) {
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_parameters(
  const Parameters& params)
{
  string           gmset_type;
  double           mq, M0;
  int              Ns;
  std::vector<int> bc;
  double           b, c;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", gmset_type);
  if (!err_optional) {
    if (gmset_type != m_repr) {
      vout.crucial(m_vl, "%s: unsupported gamma_matrix_type: %s\n",
                   class_name.c_str(), gmset_type.c_str());
      exit(EXIT_FAILURE);
    }
  }

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double("coefficient_b", b);
  err += params.fetch_double("coefficient_c", c);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(mq, M0, Ns, bc, b, c);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_parameters(
  const double mq,
  const double M0,
  const int Ns,
  const vector<int> bc,
  const double b,
  const double c)
{
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_M0   = real_t(M0);
    m_mq   = real_t(mq);
    m_Ns   = Ns;
    m_NinF = m_Nvcd * m_Ns;

    if (m_boundary.size() != m_Ndim) m_boundary.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
    }
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = real_t(b);
      m_c[is] = real_t(c);
    }
  }

  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n", m_M0);
  vout.general(m_vl, "  Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients b = %16.10f  c = %16.10f\n",
               m_b[0], m_c[0]);

  set_precond_parameters();

  // working 5d vectors.
  if (m_y1.nin() != m_NinF) {
    if (ith == 0) {
      m_y1.reset(m_NinF, m_Nst2, 1);
      m_v1.reset(m_NinF, m_Nst2, 1);
      m_v2.reset(m_NinF, m_Nst2, 1);
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_precond_parameters()
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    if (m_dp.size() != m_Ns) {
      m_dp.resize(m_Ns);
      m_dpinv.resize(m_Ns);
      m_dm.resize(m_Ns);
      m_e.resize(m_Ns - 1);
      m_f.resize(m_Ns - 1);
    }

    for (int is = 0; is < m_Ns; ++is) {
      m_dp[is] = 1.0 + m_b[is] * (4.0 - m_M0);
      m_dm[is] = 1.0 - m_c[is] * (4.0 - m_M0);
    }

    m_e[0] = m_mq * m_dm[m_Ns - 1] / m_dp[0];
    m_f[0] = m_mq * m_dm[0];
    for (int is = 1; is < m_Ns - 1; ++is) {
      m_e[is] = m_e[is - 1] * m_dm[is - 1] / m_dp[is];
      m_f[is] = m_f[is - 1] * m_dm[is] / m_dp[is - 1];
    }

    m_g = m_e[m_Ns - 2] * m_dm[m_Ns - 2];

    for (int is = 0; is < m_Ns - 1; ++is) {
      m_dpinv[is] = 1.0 / m_dp[is];
    }
    m_dpinv[m_Ns - 1] = 1.0 / (m_dp[m_Ns - 1] + m_g);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_coefficients(
  const std::vector<double> vec_b,
  const std::vector<double> vec_c)
{
#pragma omp barrier

  if ((vec_b.size() != m_Ns) || (vec_c.size() != m_Ns)) {
    vout.crucial(m_vl, "%s: size of coefficient vectors incorrect.\n",
                 class_name.c_str());
  }

  vout.general(m_vl, "%s: coefficient vectors are set:\n",
               class_name.c_str());

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = real_t(vec_b[is]);
      m_c[is] = real_t(vec_c[is]);
    }
  }

#pragma omp barrier

  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }

  set_precond_parameters();
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config(Field *u)
{
  int nth = ThreadManager::get_num_threads();

  vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                class_name.c_str(), nth);

  if (nth > 1) {
    set_config_impl(u);
  } else {
    set_config_omp(u);
  }

  vout.detailed(m_vl, "%s: set_config finished\n", class_name.c_str());
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config_omp(Field *u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config_impl(Field *u)
{
#pragma omp barrier

  AIndex_lex<real_t, QXS> index_lex;
  convert_gauge(index_lex, m_Ulex, *u);

  QXS_Gauge::set_boundary(m_Ulex, m_boundary);

  AIndex_eo<real_t, QXS> index_eo;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nvol);

  for (int ex = 0; ex < m_Ndim; ++ex) {
    for (int site = is; site < ns; ++site) {
      for (int in = 0; in < m_Ndf; ++in) {
        int iv1 = index_lex.idx(in, m_Ndf, site, ex);
        int iv2 = index_eo.idx(in, m_Ndf, site, ex);
        m_Ueo.set(iv2, m_Ulex.cmp(iv1));
      }
    }
  }

#pragma omp barrier
}


//====================================================================

/*
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_boundary_config(
                                                 AFIELD& U,
                                                 const int mu)
{ // this function may be called inside parallel region

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nvol);

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  int ipe[4], Nsize[4], Lsize[4];
  for(int i = 0; i < m_Ndim; ++i){
    ipe[i]   = Communicator::ipe(i);
    Nsize[i] = CommonParameters::Nsize(i);
    Lsize[i] = CommonParameters::Lsize(i);
  }

  real_t bc = (real_t)(m_boundary[mu]);

  int j[4];
  for(int site = is; site < ns; ++site){
    int x = site % m_Nx;
    int yzt = site/m_Nx;
    int y = yzt % m_Ny;
    int z = (yzt/m_Ny) % m_Nz;
    int t = yzt/(m_Ny*m_Nz);
    j[0] = ipe[0] * Nsize[0] + x;
    j[1] = ipe[1] * Nsize[1] + y;
    j[2] = ipe[2] * Nsize[2] + z;
    j[3] = ipe[3] * Nsize[3] + t;

   if(j[mu] == Lsize[mu]-1){
      for(int in = 0; in < m_Ndf; ++in){
        int i = index_lex.idx_G(in, site, mu);
        real_t uv = bc * U.cmp(i);
        U.set(i, uv);
      }
    }

  }

#pragma omp barrier

}
*/
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::convert(AFIELD& v,
                                               const Field& w)
{ // the following implementation only applies Dirac repr.
  assert(w.nex() == m_Ns);
  assert(w.nin()*w.nex() == v.nin()*v.nex());

  int Nvol = v.nvol();

  int ith, nth, site0, site1;
  set_threadtask_mult(ith, nth, site0, site1, Nvol);

  int Nx = m_Nx2;
  if(Nvol == m_Nvol){ // converting lexical field
    Nx = m_Nx;
  }
  AIndex_lex<real_t, AFIELD::IMPL> index(Nx, m_Ny, m_Nz, m_Nt);

#pragma omp barrier

  for (int site = site0; site < site1; ++site) {
    for (int is = 0; is < m_Ns; ++is) {
      for (int ic = 0; ic < NC; ++ic) {
        for (int id = 0; id < ND2; ++id) {
          for (int iri = 0; iri < 2; ++iri) {
            int in_org = iri + 2 * (ic + NC * id);
            int in_alt = iri + 2 * (id + ND * ic) + NVCD * is;
            v.set(index.idx(in_alt, m_NinF, site, 0),
                  real_t(w.cmp(in_org, site, is)));
          }
        }

        for (int id = ND2; id < ND; ++id) {
          for (int iri = 0; iri < 2; ++iri) {
            int in_org = iri + 2 * (ic + NC * id);
            int in_alt = iri + 2 * (id + ND * ic) + NVCD * is;
            v.set(index.idx(in_alt, m_NinF, site, 0),
                  real_t(-w.cmp(in_org, site, is)));
          }
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::reverse(Field& v,
                                               const AFIELD& w)
{ // the following implementation only applies Dirac repr.
  int Nvol = v.nvol();

  int ith, nth, site0, site1;
  set_threadtask_mult(ith, nth, site0, site1, Nvol);
  int Nx = m_Nx2;
  if(Nvol == m_Nvol){ // converting lexical field
    Nx = m_Nx;
  }
  AIndex_lex<real_t, AFIELD::IMPL> index(Nx, m_Ny, m_Nz, m_Nt);

#pragma omp barrier

  for (int site = site0; site < site1; ++site) {
    for (int is = 0; is < m_Ns; ++is) {
      for (int ic = 0; ic < NC; ++ic) {
        for (int id = 0; id < ND2; ++id) {
          for (int iri = 0; iri < 2; ++iri) {
            int in_org = iri + 2 * (ic + NC * id);
            int in_alt = iri + 2 * (id + ND * ic) + NVCD * is;
            v.set(in_org, site, is,
                  double(  w.cmp(index.idx(in_alt, m_NinF, site, 0))));
          }
        }

        for (int id = ND2; id < ND; ++id) {
          for (int iri = 0; iri < 2; ++iri) {
            int in_org = iri + 2 * (ic + NC * id);
            int in_alt = iri + 2 * (id + ND * ic) + NVCD * is;
            v.set(in_org, site, is,
                  double( -w.cmp(index.idx(in_alt, m_NinF, site, 0))));
          }
        }
      }
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    D(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    H(v, w);
  } else if (m_mode == "Hdag") {
    Hdag(v, w);
  }else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_dag(AFIELD& v,
                                                const AFIELD& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
  } else if (m_mode == "H") {
    Hdag(v, w);
  } else if (m_mode == "Hdag") {
    H(v, w);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult(AFIELD& v,
                                            const AFIELD& w,
                                            std::string mode)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

  if (mode == "Deo") {
    D_eo(v, w, 0);
  } else if (mode == "Doe") {
    D_eo(v, w, 1);
  } else if (mode == "Dee") {
    D_ee(v, w, 0);
  } else if (mode == "Doo") {
    D_ee(v, w, 1);
  } else if (mode == "Dee_inv") {
    D_ee_inv(v, w, 0);
  } else if (mode == "Doo_inv") {
    D_ee_inv(v, w, 1);
  } else {
    std::cout << "mode undeifined in AFopr_Domainwall_eo.\n";
    abort();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_dag(AFIELD& v,
                                                const AFIELD& w,
                                                std::string mode)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

  if (mode == "Deo") {
    Ddag_eo(v, w, 1);
  } else if (mode == "Doe") {
    Ddag_eo(v, w, 0);
  } else if (mode == "Dee") {
    Ddag_ee(v, w, 0);
  } else if (mode == "Doo") {
    Ddag_ee(v, w, 1);
  } else if (mode == "Dee_inv") {
    Ddag_ee_inv(v, w, 0);
  } else if (mode == "Doo_inv") {
    Ddag_ee_inv(v, w, 1);
  } else {
    std::cout << "mode undeifined in AFopr_Domainwall_eo.\n";
    abort();
  }
}

//CG法で最も処理に時間がかかる部分DdagDを作用させるための計算
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::DdagD(AFIELD& v,
                                             const AFIELD& w)
{
#ifdef USE_DAGD_DETAIL_PROFILE
  const int profile_ith = ThreadManager::get_thread_id();
  Timer profile_total;
  Timer profile_stage;
  profile_total.start();

  profile_stage.start();
  D_eo(m_v1, w, 1);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[1] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  L_inv(m_v2, m_v1);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[2] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  U_inv(m_v1, m_v2);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[3] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  D_eo(m_v2, m_v1, 0);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[4] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  L_inv(m_v1, m_v2);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[5] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  U_inv(m_v2, m_v1);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[6] += profile_stage.elapsed_sec();

  // Fair comparison: a single outer timer encloses the same logical
  // z = x - y block in both the fused and unfused branches.
  profile_stage.reset(); profile_stage.start();
#ifdef USE_DAGD_FUSED_COPY_AXPY
  copy_axpy(m_v1, w, real_t(-1.0), m_v2);
#else
  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);
#endif
#pragma omp barrier
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[7] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  Udag_inv(v, m_v1);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[8] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  Ldag_inv(m_v2, v);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[9] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  Ddag_eo(v, m_v2, 1);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[10] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  Udag_inv(m_v2, v);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[11] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  Ldag_inv(v, m_v2);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[12] += profile_stage.elapsed_sec();

  profile_stage.reset(); profile_stage.start();
  Ddag_eo(m_v2, v, 0);
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[13] += profile_stage.elapsed_sec();

  // Fair comparison: use the same single timer start/stop pair here too.
  profile_stage.reset(); profile_stage.start();
#ifdef USE_DAGD_FUSED_COPY_AXPY
  copy_axpy(v, m_v1, real_t(-1.0), m_v2);
#else
  copy(v, m_v1);
  axpy(v, real_t(-1.0), m_v2);
#endif
#pragma omp barrier
  profile_stage.stop();
  if (profile_ith == 0) m_dagd_profile_sec[14] += profile_stage.elapsed_sec();

  profile_total.stop();
  if (profile_ith == 0) {
    m_dagd_profile_sec[0] += profile_total.elapsed_sec();
    ++m_dagd_profile_calls;
  }
#else
  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

#ifdef USE_DAGD_FUSED_COPY_AXPY
  copy_axpy(m_v1, w, real_t(-1.0), m_v2);
#else
  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);
#endif
#pragma omp barrier

  Udag_inv(v, m_v1);
  Ldag_inv(m_v2, v);
  Ddag_eo(v, m_v2, 1);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
  Ddag_eo(m_v2, v, 0);

#ifdef USE_DAGD_FUSED_COPY_AXPY
  copy_axpy(v, m_v1, real_t(-1.0), m_v2);
#else
  copy(v, m_v1);
  axpy(v, real_t(-1.0), m_v2);
#endif
#pragma omp barrier
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  D_eo(m_v1, w, 1);
  L_inv(m_v2, m_v1);
  U_inv(m_v1, m_v2);
  D_eo(m_v2, m_v1, 0);
  L_inv(m_v1, m_v2);
  U_inv(m_v2, m_v1);

  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);
  mult_gm5R(v, m_v1);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  Udag_inv(m_v1, w);
  Ldag_inv(m_v2, m_v1);
  Ddag_eo(m_v1, m_v2, 1);
  Udag_inv(m_v2, m_v1);
  Ldag_inv(m_v1, m_v2);
  Ddag_eo(m_v2, m_v1, 0);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Hdag(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier
  mult_gm5R(v, w);

  Udag_inv(m_v1, v);
  Ldag_inv(m_v2, m_v1);
  Ddag_eo(m_v1, m_v2, 1);
  Udag_inv(m_v2, m_v1);
  Ldag_inv(m_v1, m_v2);
  Ddag_eo(m_v2, m_v1, 0);

  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_eo(AFIELD& v,
                                            const AFIELD& w,
                                            const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  mult_D_eo(v, w, ieo);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_eo(AFIELD& v,
                                               const AFIELD& w,
                                               const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  mult_Ddag_eo(v, w, ieo);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  // gm5 is the same as the lexical case
  BridgeQXS::mult_domainwall_5din_mult_gm5_dirac(vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_R(AFIELD& v, const AFIELD& w)
{
  vout.crucial("%s: mult_R, not tested\n", class_name.c_str());
  abort();

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  // R is the same as the lexical case
  BridgeQXS::mult_domainwall_5din_mult_R(vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_gm5R(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  // gm5R is the same as the lexical case
  BridgeQXS::mult_domainwall_5din_mult_gm5R_dirac(vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_D_eo(AFIELD& v,
                                                 const AFIELD& w,
                                                 const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_y1.ptr(0);  // working vector
  real_t *up = m_Ueo.ptr(0);

  int ith = ThreadManager::get_thread_id();

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_eo_clear(vp, m_Ns, m_Nsize);

  BridgeQXS::mult_domainwall_5din_eo_5dir_dirac(yp, wp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                m_Nsize, do_comm);

#pragma omp barrier

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.start();

    real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
    real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
    real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
    real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
    real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
    real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
    real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
    real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

    BridgeQXS::mult_domainwall_5din_eo_hop1_dirac(
      buf1_xp, buf1_xm, buf1_yp, buf1_ym,
      buf1_zp, buf1_zm, buf1_tp, buf1_tm,
      up, yp,
      m_mq, m_M0, m_Ns, &m_boundary[0],
      &m_Leo[0], m_Nsize, do_comm, ieo);

#pragma omp barrier

    if (ith == 0) chset_send.start();
  }

  BridgeQXS::mult_domainwall_5din_eo_hopb_dirac(vp, up, yp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                &m_Leo[0], m_Nsize, do_comm, ieo);

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.wait();

#pragma omp barrier

    real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
    real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
    real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
    real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
    real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
    real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
    real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
    real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

    BridgeQXS::mult_domainwall_5din_eo_hop2_dirac(vp, up, yp,
                                                  buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                                  buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                                  m_mq, m_M0, m_Ns, &m_boundary[0],
                                                  &m_Leo[0], m_Nsize, do_comm, ieo);

    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_Ddag_eo(AFIELD& v,
                                                    const AFIELD& w,
                                                    const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_y1.ptr(0);  // working vector
  real_t *up = m_Ueo.ptr(0);

  int ith = ThreadManager::get_thread_id();

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_eo_mult_gm5_dirac(vp, wp, m_Ns, m_Nsize);

  BridgeQXS::mult_domainwall_5din_eo_clear(yp, m_Ns, m_Nsize);

#pragma omp barrier

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.start();

    real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
    real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
    real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
    real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
    real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
    real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
    real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
    real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

    BridgeQXS::mult_domainwall_5din_eo_hop1_dirac(
      buf1_xp, buf1_xm, buf1_yp, buf1_ym,
      buf1_zp, buf1_zm, buf1_tp, buf1_tm,
      up, vp,
      m_mq, m_M0, m_Ns, &m_boundary[0],
      &m_Leo[0], m_Nsize, do_comm, ieo);

#pragma omp barrier

    if (ith == 0) chset_send.start();
  }

  BridgeQXS::mult_domainwall_5din_eo_hopb_dirac(yp, up, vp,
                                                m_mq, m_M0, m_Ns, &m_boundary[0],
                                                &m_b[0], &m_c[0],
                                                &m_Leo[0], m_Nsize, do_comm, ieo);

  if (do_comm_any > 0) {
    if (ith == 0) chset_recv.wait();

#pragma omp barrier

    real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
    real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
    real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
    real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
    real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
    real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
    real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
    real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

    BridgeQXS::mult_domainwall_5din_eo_hop2_dirac(yp, up, vp,
                                                  buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                                                  buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                                                  m_mq, m_M0, m_Ns, &m_boundary[0],
                                                  &m_Leo[0], m_Nsize, do_comm, ieo);

    if (ith == 0) chset_send.wait();
  }

#pragma omp barrier

  BridgeQXS::mult_domainwall_5din_eo_5dirdag_dirac(vp, yp,
                                                   m_mq, m_M0, m_Ns, &m_boundary[0],
                                                   &m_b[0], &m_c[0],
                                                   m_Nsize, do_comm);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_ee_inv(AFIELD& v,
                                                const AFIELD& w,
                                                const int ieo)
{
#pragma omp barrier
  L_inv(m_v1, w);
  U_inv(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_ee_inv(AFIELD& v,
                                                   const AFIELD& w,
                                                   const int ieo)
{
#pragma omp barrier
  Udag_inv(m_v1, w);
  Ldag_inv(v, m_v1);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_ee(AFIELD& v,
                                            const AFIELD& w,
                                            const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even (no difference)
  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_mult(ith, nth, site0, site1, m_Nst2v);

#pragma omp barrier

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t vL[NVCD], wL[NVCD];

    for (int is = 0; is < m_Ns; ++is) {
      real_t FF1 = m_b[is] * (4.0 - m_M0) + 1.0;
      real_t FF2 = m_c[is] * (4.0 - m_M0) - 1.0;

      load_vec(wL, &wp2[Nin4 * is], NVCD);
      set_vec(vL, FF1, wL, NVCD);

      int    is_up = (is + 1) % m_Ns;
      real_t Fup   = 0.5 * FF2;
      if (is == m_Ns - 1) Fup *= -m_mq;
      load_vec(wL, &wp2[Nin4 * is_up], NVCD);
      add_aPm5_dirac_vec(vL, Fup, wL, NC);

      int    is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn   = 0.5 * FF2;
      if (is == 0) Fdn *= -m_mq;
      load_vec(wL, &wp2[Nin4 * is_dn], NVCD);
      add_aPp5_dirac_vec(vL, Fdn, wL, NC);

      save_vec(&vp2[Nin4 * is], vL, NVCD);
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_ee(AFIELD& v,
                                               const AFIELD& w,
                                               const int ieo)
{   // ieo = 0: even < odd, ieo = 1: odd <- even
  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_mult(ith, nth, site0, site1, m_Nst2v);

#pragma omp barrier

  for (int site = site0; site < site1; ++site) {
    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t vL[NVCD], wL[NVCD];

    for (int is = 0; is < m_Ns; ++is) {
      real_t FF1 = m_b[is] * (4.0 - m_M0) + 1.0;
      load_vec(vL, &wp2[Nin4 * is], NVCD);
      scal_vec(vL, FF1, NVCD);

      int    is_up = (is + 1) % m_Ns;
      real_t Fup   = 0.5 * (m_c[is_up] * (4.0 - m_M0) - 1.0);
      if (is == m_Ns - 1) Fup *= -m_mq;
      load_vec(wL, &wp2[Nin4 * is_up], NVCD);
      add_aPp5_dirac_vec(vL, Fup, wL, NC);

      int    is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn   = 0.5 * (m_c[is_dn] * (4.0 - m_M0) - 1.0);
      if (is == 0) Fdn *= -m_mq;
      load_vec(wL, &wp2[Nin4 * is_dn], NVCD);
      add_aPm5_dirac_vec(vL, Fdn, wL, NC);

      save_vec(&vp2[Nin4 * is], vL, NVCD);
    }
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::L_inv(AFIELD& v,
                                             const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  BridgeQXS::mult_domainwall_5din_eo_L_inv_dirac(
    vp, wp, m_Ns, m_Nsize,
    &m_e[0], &m_dpinv[0], &m_dm[0]);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::U_inv(AFIELD& v,
                                             const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  BridgeQXS::mult_domainwall_5din_eo_U_inv_dirac(
    vp, wp, m_Ns, m_Nsize,
    &m_f[0], &m_dpinv[0], &m_dm[0]);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Udag_inv(AFIELD& v,
                                                const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  BridgeQXS::mult_domainwall_5din_eo_Udag_inv_dirac(
    vp, wp, m_Ns, m_Nsize,
    &m_f[0], &m_dpinv[0], &m_dm[0]);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ldag_inv(AFIELD& v,
                                                const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  BridgeQXS::mult_domainwall_5din_eo_Ldag_inv_dirac(
    vp, wp, m_Ns, m_Nsize,
    &m_e[0], &m_dpinv[0], &m_dm[0]);
}


//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_5din_eo<AFIELD>::flop_count(std::string mode)
{
  // flop counting of this class is not confirmed.

  int    Lvol  = CommonParameters::Lvol();
  double vsite = static_cast<double>(Lvol);
  double vNs   = static_cast<double>(m_Ns);

  double flop_Wilson_hop;
  double flop_LU_inv;
  double flop_pre;
  int    Nc = CommonParameters::Nc();
  int    Nd = CommonParameters::Nd();
  if (m_repr == "Dirac") {
    flop_Wilson_hop = static_cast<double>(
      vNs * Nc * Nd * (6 * (4 * Nc + 2) + 2 * (4 * Nc + 1))) * vsite;
    flop_pre    = static_cast<double>(vNs * Nc * Nd * 14) * vsite;
    flop_LU_inv = static_cast<double>(Nc * Nd * (2 + (vNs - 1) * 26)) * vsite;
  } else if (m_repr == "Chiral") {
    flop_Wilson_hop = static_cast<double>(
      vNs * Nc * Nd * (8 * (4 * Nc + 2))) * vsite;
    flop_pre    = static_cast<double>(vNs * Nc * Nd * 6) * vsite;
    flop_LU_inv = static_cast<double>(Nc * Nd * (2 + (vNs - 1) * 10)) * vsite;
  }
  double flop_axpy = static_cast<double>(vNs * 2 * Nc * Nd) * vsite;

  //  double axpy1 = static_cast<double>(2 * m_NinF);
  //  double scal1 = static_cast<double>(1 * m_NinF);
  //  double flop_DW = vNs * (flop_Wilson + vsite*(6*axpy1 + 2*scal1));
  // In Ddag case, flop_Wilson + 7 axpy which equals flop_DW.
  double flop_Deo = flop_Wilson_hop + flop_pre;


  //  double flop_LU_inv = 2.0 * vsite *
  //            ((3.0*axpy1 + scal1)*(vNs-1.0) + axpy1 + 2.0*scal1);

  double flop = 0.0;
  if ((mode == "D") || (mode == "Ddag")) {
    flop = flop_Deo + flop_LU_inv + 0.5 * flop_axpy;  // axpy is for 1 - (Deo etc)
  } else if (mode == "DdagD") {
    flop = 2 * flop_Deo + 2 * flop_LU_inv + flop_axpy;
  } else if ((mode == "Dee_inv") || (mode == "Doo_inv")) {
    flop = 0.5 * flop_LU_inv; // 0.5 is for even-odd
  } else if ((mode == "Dee") || (mode == "Doo")) {
    flop = 0.5 * flop_pre;
  } else if ((mode == "Doe") || (mode == "Doe")) {
    flop = 0.5 * flop_Deo;  // 0.5 is for even-odd
  } else {
    vout.crucial(m_vl, "Error at %s: input mode is undefined: %s.\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  return flop;
}


//============================================================END=====
