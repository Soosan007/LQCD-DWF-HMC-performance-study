#ifndef ASOLVER_CG_H
#define ASOLVER_CG_H

#include <cstdio>
#include <cstdlib>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "lib/Fopr/afopr.h"

#include "lib_alt/Solver/asolver.h"

template<typename AFIELD>
class ASolver_CG : public ASolver<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using ASolver<AFIELD>::m_vl;
  static const std::string class_name;
  using InitialGuess = typename ASolver<AFIELD>::InitialGuess;

 protected:

  AFopr<AFIELD> *m_fopr; //!< fermion operator.

  int m_Niter;           //!< maximum iteration number.
  real_t m_Stop_cond;    //!< stopping criterion (squared).

  //! to remember convergence iteration to provide flop count.
  int m_nconv;

  //! mode switch for initial guess
  InitialGuess m_initial_mode;

  //! optional automatic cache for a previous compatible solution.
  bool m_reuse_previous_solution;
  bool m_previous_solution_valid;
  std::string m_previous_solution_mode;
  AFIELD m_previous_solution;

  //! calling constructor without fermion operator is forbidden.
  ASolver_CG() { }
  
  // CG法で計測したものを保存するための変数
  // Per-solve wall-clock profile accumulated by the OpenMP master thread.
  unsigned long long m_profile_solve_calls;
  double m_profile_setup_copy_b;
  double m_profile_rhs_norm2;
  double m_profile_init;
  double m_profile_steps;
  double m_profile_step_mult;
  double m_profile_step_dot;
  double m_profile_step_axpy_x;
  double m_profile_step_axpy_r;
  double m_profile_step_norm2_r;
  double m_profile_step_axpy_norm2;
  //! Fair comparison timer: one start/stop pair around the complete
  //! residual update block in both fused and unfused builds.
  double m_profile_step_residual_block;
  double m_profile_step_aypx;
  double m_profile_final;
  double m_profile_final_mult;
  double m_profile_final_axpy;
  double m_profile_final_norm2;

  //! working vectors.
  AFIELD m_x, m_r, m_p, m_s;

 public:
  //! constructor.
  ASolver_CG(AFopr<AFIELD> *fopr)
    : m_Niter(0), m_Stop_cond(0.0L),
      m_reuse_previous_solution(false), m_previous_solution_valid(false),
      m_profile_solve_calls(0)
  {
    m_fopr = fopr;
    this->init();
  }

  //! destructor.
  ~ASolver_CG() { this->tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond);

  //! setting parameters.
  void set_parameters(const int Niter, const real_t Stop_cond, const InitialGuess init_guess_mode);

  void set_init_mode(const InitialGuess init_guess)
  {
    m_initial_mode = init_guess;
  }

  //! solver main.
  void solve(AFIELD& xq, const AFIELD& b, int& nconv, real_t& diff);

  //! returns the pointer to the fermion operator.
  AFopr<AFIELD> *get_fopr() { return m_fopr; }

  //! returns the floating point operation count.
  double flop_count();


 protected:

  void init(void);

  void tidyup(void);

  void solve_CG_init(const AFIELD& xq, const InitialGuess init_mode,
                     real_t& rrp, real_t& rr);

  void solve_CG_step(real_t& rrp, real_t& rr);


#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD> *fopr)
  { return new ASolver_CG<AFIELD>(fopr); }

 public:
  static bool register_factory()
  {
    bool init = ASolver<AFIELD>::Factory_fopr::Register("CG",
                                                        create_object_with_fopr);
    return init;
  }
#endif
};

#endif // ASOLVER_CG_H
