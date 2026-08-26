#include "lib/fapp_macros.h"
#include "lib/Tools/timer.h"

template<typename AFIELD>
const std::string ASolver_CG<AFIELD>::class_name = "ASolver_CG";

/*
    diff returns relative residual (was absolute squared)
                                 [24 Jul 2023 I.Knamori]
*/

//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_x.reset(nin, nvol, nex);
  m_r.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);
  m_s.reset(nin, nvol, nex);
  m_previous_solution.reset(nin, nvol, nex);

  //  m_vl = Bridge::DETAILED;

  m_nconv = -1;
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::tidyup(void)
{
  // ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  InitialGuess init_guess_mode = InitialGuess::RHS;
  if (params.find_string("initial_guess_mode")) {
    const string initial_guess_mode = params.get_string("initial_guess_mode");
    vout.detailed(m_vl, "  initila_guess_mode %s\n", initial_guess_mode.c_str());
    if (initial_guess_mode == "RHS") {
      init_guess_mode = InitialGuess::RHS;
    } else if (initial_guess_mode == "GIVEN") {
      init_guess_mode = InitialGuess::GIVEN;
    } else if (initial_guess_mode == "ZERO") {
      init_guess_mode = InitialGuess::ZERO;
    } else {
      vout.crucial(m_vl, "Error at %s: unknown initial guess mode, %s\n", class_name.c_str(), initial_guess_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  bool reuse_previous_solution = false;
  if (params.find_bool("reuse_previous_solution")) {
    reuse_previous_solution = params.get_bool("reuse_previous_solution");
  }

  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond, init_guess_mode);
  m_reuse_previous_solution = reuse_previous_solution;
  m_previous_solution_valid = false;
  vout.general(m_vl, "  reuse_previous_solution: %s\n",
               m_reuse_previous_solution ? "true" : "false");
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const int Niter,
                                        const real_t Stop_cond)
{
  set_parameters(Niter, Stop_cond, InitialGuess::RHS);
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::set_parameters(const int Niter,
                                        const real_t Stop_cond,
                                        const InitialGuess init_guess_mode)
{
  ThreadManager::assert_single_thread(class_name);

  m_Niter        = Niter;
  m_Stop_cond    = Stop_cond;
  m_initial_mode = init_guess_mode;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %16.8e\n", m_Stop_cond);
  vout.general(m_vl, "  init_guess_mode: %d\n", m_initial_mode);
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve(AFIELD& xq, const AFIELD& b,
                               int& Nconv, real_t& diff)
{
#ifdef USE_SOLVER_MULT_TRACE
#pragma omp master
  {
    const std::string solver_mult_mode = m_fopr->get_mode();
    vout.general(m_vl,
      "SOLVER_MULT_TRACE event=begin rank=%d solver=%p mode=%s\n",
      Communicator::nodeid(), static_cast<void *>(this),
      solver_mult_mode.c_str());
  }
#pragma omp barrier
#endif

  Timer profile_total;
  Timer profile_part;
  unsigned long long profile_call = 0;
  int iterations_executed = 0;

  profile_total.start();
#pragma omp master
  {
    profile_call = ++m_profile_solve_calls;
    m_profile_setup_copy_b   = 0.0;
    m_profile_rhs_norm2      = 0.0;
    m_profile_init           = 0.0;
    m_profile_steps          = 0.0;
    m_profile_step_mult      = 0.0;
    m_profile_step_dot       = 0.0;
    m_profile_step_axpy_x    = 0.0;
    m_profile_step_axpy_r    = 0.0;
    m_profile_step_norm2_r   = 0.0;
    m_profile_step_axpy_norm2 = 0.0;
    m_profile_step_residual_block = 0.0;
    m_profile_step_aypx      = 0.0;
    m_profile_final          = 0.0;
    m_profile_final_mult     = 0.0;
    m_profile_final_axpy     = 0.0;
    m_profile_final_norm2    = 0.0;
  }
#pragma omp barrier

  profile_part.start();
  copy(m_s, b);
  profile_part.stop();
#pragma omp master
  { m_profile_setup_copy_b = profile_part.elapsed_sec(); }

  profile_part.reset();
  profile_part.start();
  real_t sr    = norm2(m_s);
  profile_part.stop();
#pragma omp master
  { m_profile_rhs_norm2 = profile_part.elapsed_sec(); }

  real_t snorm = 1.0 / sr;
  vout.detailed(m_vl, "  snorm = %22.15e\n", snorm);

  real_t rr, rrp;
  int    nconv = -1;

  const std::string current_solver_mode = m_fopr->get_mode();
  const bool use_previous_solution =
    m_reuse_previous_solution && m_previous_solution_valid
    && (m_previous_solution_mode == current_solver_mode);
  const InitialGuess effective_initial_mode =
    use_previous_solution ? InitialGuess::GIVEN : m_initial_mode;
  const AFIELD& initial_solution =
    use_previous_solution ? m_previous_solution : xq;

  profile_part.reset();
  profile_part.start();
  solve_CG_init(initial_solution, effective_initial_mode, rrp, rr);
  profile_part.stop();
#pragma omp master
  { m_profile_init = profile_part.elapsed_sec(); }
  vout.detailed(m_vl, "  init: %22.15e\n", rr * snorm);

  // A supplied initial guess may already satisfy the stopping criterion.
  // Avoid an unnecessary CG step (and a possible 0/0 in alpha) in that case.
  if (rr * snorm < m_Stop_cond) {
    nconv = 0;
#pragma omp master
    { m_nconv = 0; }
  }

  for (int iter = 0; (nconv == -1) && (iter < m_Niter); ++iter) {
    profile_part.reset();
    profile_part.start();
    solve_CG_step(rrp, rr);
    profile_part.stop();
#pragma omp master
    { m_profile_steps += profile_part.elapsed_sec(); }
    iterations_executed = iter + 1;
    vout.detailed(m_vl, "%6d  %22.15e\n", iter, rr * snorm);

    if (rr * snorm < m_Stop_cond) {
      nconv = iter;
#pragma omp master
      {
        m_nconv = nconv + 1;
      }
      break;
    }
  }

  if (nconv == -1)  {
    vout.crucial(m_vl, "Error at %s: not converged\n",
                 class_name.c_str());
    vout.crucial(m_vl, "  iter(final): %8d  %22.15e\n",
                 m_Niter, rr * snorm);
#pragma omp barrier
    //exit(EXIT_FAILURE);
  }

  vout.detailed(m_vl, "converged:\n");
  vout.detailed(m_vl, "  nconv = %d\n", nconv);

  profile_part.reset();
  profile_part.start();

  copy(xq, m_x);

  if (m_reuse_previous_solution) {
    copy(m_previous_solution, m_x);
#pragma omp master
    {
      m_previous_solution_mode  = current_solver_mode;
      m_previous_solution_valid = true;
    }
#pragma omp barrier
  }

  START_FAPP("mult", 1, 1);
  Timer profile_final_op;
  profile_final_op.start();
  m_fopr->mult(m_s, xq);
  profile_final_op.stop();
#pragma omp master
  { m_profile_final_mult = profile_final_op.elapsed_sec(); }
  STOP_FAPP("mult", 1, 1);

  profile_final_op.reset();
  profile_final_op.start();
  axpy(m_s, real_t(-1.0), b);
  profile_final_op.stop();
#pragma omp master
  { m_profile_final_axpy = profile_final_op.elapsed_sec(); }

  profile_final_op.reset();
  profile_final_op.start();
  real_t diff2 = norm2(m_s) * snorm;
  profile_final_op.stop();
#pragma omp master
  { m_profile_final_norm2 = profile_final_op.elapsed_sec(); }

  profile_part.stop();
#pragma omp master
  { m_profile_final = profile_part.elapsed_sec(); }

#pragma omp master
  {
    diff  = sqrt(diff2);
    Nconv = m_nconv;
  }
#pragma omp barrier

  profile_total.stop();
#pragma omp master
  {
#ifdef USE_QXS_ACLE
    const char *profile_backend = "ACLE";
#else
    const char *profile_backend = "GENERAL";
#endif
#ifdef USE_CG_FUSED_AXPY_NORM2
    const int profile_fused_axpy_norm2 = 1;
#else
    const int profile_fused_axpy_norm2 = 0;
#endif
    const char *profile_precision = (sizeof(real_t) == 4) ? "float" : "double";
    const double profile_step_accounted =
      m_profile_step_mult + m_profile_step_dot + m_profile_step_axpy_x
      + m_profile_step_residual_block + m_profile_step_aypx;
    const double profile_step_other = m_profile_steps - profile_step_accounted;
    const double profile_avg_iteration =
      iterations_executed > 0 ? m_profile_steps / iterations_executed : 0.0;

    vout.general(m_vl,
      "CG_PROFILE backend=%s rank=%d solver=%p call=%llu precision=%s "
      "fused_axpy_norm2=%d timing_scheme=single_outer_block "
      "previous_guess_used=%d converged=%d iterations=%d final_relative_residual=%.15e "
      "total_sec=%.9e setup_copy_b_sec=%.9e rhs_norm2_sec=%.9e init_sec=%.9e "
      "steps_sec=%.9e avg_iteration_sec=%.9e step_mult_sec=%.9e "
      "step_dot_sec=%.9e step_axpy_x_sec=%.9e step_axpy_r_sec=%.9e "
      "step_norm2_r_sec=%.9e step_axpy_norm2_sec=%.9e "
      "step_residual_block_sec=%.9e step_aypx_sec=%.9e "
      "step_other_sec=%.9e final_sec=%.9e final_mult_sec=%.9e "
      "final_axpy_sec=%.9e final_norm2_sec=%.9e\n",
      profile_backend, Communicator::nodeid(), static_cast<void *>(this),
      profile_call, profile_precision, profile_fused_axpy_norm2,
      use_previous_solution ? 1 : 0, (nconv != -1) ? 1 : 0,
      iterations_executed, static_cast<double>(diff),
      profile_total.elapsed_sec(), m_profile_setup_copy_b, m_profile_rhs_norm2,
      m_profile_init, m_profile_steps, profile_avg_iteration,
      m_profile_step_mult, m_profile_step_dot, m_profile_step_axpy_x,
      m_profile_step_axpy_r, m_profile_step_norm2_r,
      m_profile_step_axpy_norm2, m_profile_step_residual_block,
      m_profile_step_aypx, profile_step_other,
      m_profile_final, m_profile_final_mult, m_profile_final_axpy,
      m_profile_final_norm2);
#ifdef USE_SOLVER_MULT_TRACE
    const std::string solver_mult_mode = m_fopr->get_mode();
    const int solver_mult_calls = iterations_executed + 1
      + ((effective_initial_mode == InitialGuess::ZERO) ? 0 : 1);
    vout.general(m_vl,
      "SOLVER_MULT_TRACE event=end rank=%d solver=%p mode=%s "
      "converged=%d iterations=%d mult_calls=%d final_relative_residual=%.15e\n",
      Communicator::nodeid(), static_cast<void *>(this),
      solver_mult_mode.c_str(), (nconv != -1) ? 1 : 0,
      iterations_executed, solver_mult_calls, static_cast<double>(diff));
#endif
  }
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_init(const AFIELD& xq,
                                      const InitialGuess init_mode,
                                      real_t& rrp, real_t& rr)
{
  if (init_mode == InitialGuess::RHS) {
#ifdef DEBUG
    vout.general(m_vl, "%s: using InitialGuess::RHS\n", class_name.c_str());
#endif
    copy(m_r, m_s);
    copy(m_x, m_s);
    m_fopr->mult(m_s, m_x);
    axpy(m_r, real_t(-1.0), m_s);
    copy(m_p, m_r);
    rr  = norm2(m_r);
    rrp = rr;
  } else if (init_mode == InitialGuess::GIVEN) {
#ifdef DEBUG
    vout.general(m_vl, "%s: using InitialGuess::GIVEN\n", class_name.c_str());
#endif
    // Use the solution supplied by the caller as x_0 and construct
    // r_0 = b - A x_0.  The caller can pass the previous solve result in xq.
    copy(m_r, m_s);
    copy(m_x, xq);
    m_fopr->mult(m_s, m_x);
    axpy(m_r, real_t(-1.0), m_s);
    copy(m_p, m_r);
    rr  = norm2(m_r);
    rrp = rr;
  } else if (init_mode == InitialGuess::ZERO) {
#ifdef DEBUG
    vout.general(m_vl, "%s: using InitialGuess::ZERO\n", class_name.c_str());
#endif
    copy(m_r, m_s);
    m_s.set(0.0);
    m_x.set(0.0);
    copy(m_p, m_r);
    rr  = norm2(m_r);
    rrp = rr;
  } else {
    vout.crucial("%s: unkown init guess mode\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}

// 解析用メモ　インターン生が記入したものなので削除してよいです　2026/08/20更新
// multはm_foprに設定されているフェルミオン演算子をm_pに作用せ、結果をm_sに格納する処理をしている
// multがDdagDを作用させている部分なので最も重そう(予測であり確定情報ではない)
//====================================================================
template<typename AFIELD>
void ASolver_CG<AFIELD>::solve_CG_step(real_t& rrp, real_t& rr)
{
  using complex_t = typename AFIELD::complex_t;

  Timer profile_op;
  profile_op.start();
  m_fopr->mult(m_s, m_p);
  profile_op.stop();
#pragma omp master
  { m_profile_step_mult += profile_op.elapsed_sec(); }

  profile_op.reset();
  profile_op.start();
  real_t pap = dot(m_s, m_p);
  profile_op.stop();
#pragma omp master
  { m_profile_step_dot += profile_op.elapsed_sec(); }
  //    m_fopr->mult_normA_dev(pap, m_s, m_p);
  real_t cr = rrp / pap;

  profile_op.reset();
  profile_op.start();
  axpy(m_x, cr, m_p);
  profile_op.stop();
#pragma omp master
  { m_profile_step_axpy_x += profile_op.elapsed_sec(); }

  // Fair performance comparison: both branches use exactly one timer
  // start/stop pair around the same logical residual-update block.
  // The old per-operation timers are intentionally not used here because
  // two timers in the unfused branch and one in the fused branch bias both
  // the block measurement and the enclosing CG/HMC timings.
  profile_op.reset();
  profile_op.start();
#ifdef USE_CG_FUSED_AXPY_NORM2
  rr = axpy_norm2(m_r, -cr, m_s); //! AXPYによるベクトル更新とノルム二乗計算を1回の走査で行う融合処理。
#else
  axpy(m_r, -cr, m_s);
  rr = norm2(m_r);
#endif
  profile_op.stop();
#pragma omp master
  { m_profile_step_residual_block += profile_op.elapsed_sec(); }

  real_t bk = rr / rrp;

  profile_op.reset();
  profile_op.start();
  aypx(bk, m_p, m_r);
  profile_op.stop();
#pragma omp master
  { m_profile_step_aypx += profile_op.elapsed_sec(); }

  rrp = rr;
}


//====================================================================
template<typename AFIELD>
double ASolver_CG<AFIELD>::flop_count()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  int ninit = 1;

  double flop_field  = static_cast<double>(Nin * Nvol * Nex) * NPE;
  double flop_vector = (6 + ninit * 4 + m_nconv * 11) * flop_field;
  double flop_fopr   = (1 + ninit + m_nconv) * m_fopr->flop_count();

  double flop = flop_vector + flop_fopr;

  return flop;
}


//====================================================================
//============================================================END=====
