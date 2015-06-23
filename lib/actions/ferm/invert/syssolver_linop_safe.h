// -*- C++ -*-
// $Id: syssolver_linop_safe.h, 20012-03-06 15:20:54 bglaessle jnajjar Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_safe_h__
#define __syssolver_linop_safe_h__
#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_safe_params.h"
#include "actions/ferm/invert/invcg2.h"
#include "actions/ferm/invert/invbicgstab.h"


namespace Chroma
{

  //! CG system solver namespace
  namespace LinOpSysSolverSafeEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverSafe : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverSafe(Handle< LinearOperator<T> > A_,
       const SysSolverSafeParams& invParam_) :
      A(A_), invParam(invParam_)
      { numNonConverged = new int;
        *numNonConverged = 0;
      }

    //! Destructor is automatic
    ~LinOpSysSolverSafe() {
      delete numNonConverged;
    }

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const
      {
   START_CODE();
   SystemSolverResults_t res;  // initialized by a constructor
   StopWatch swatch;
   swatch.reset();
   swatch.start();

   bool didnot_compute_res=true;
//    QDPIO::cout << *numNonConverged << std::endl;

   try {
      if( *numNonConverged > 0 and invParam.Strategy==0 )
       throw( std::string("Try slow CG muahaha.") );

            res = InvBiCGStab(*A,
         chi,
         psi,
         invParam.RsdBiCGStab,
         invParam.MaxBiCGStab,
         PLUS);
         {
            T r;
            r[A->subset()]=chi;
            T tmp;
            (*A)(tmp, psi, PLUS);
            r[A->subset()] -= tmp;
            res.resid = sqrt(norm2(r, A->subset()));
         }
           didnot_compute_res=false;
         if(toBool(res.resid/sqrt(norm2(chi,A->subset())) > invParam.max_BiCG_relative_Res)){
            didnot_compute_res=true;
            QDPIO::cout << "SAFE_SOLVER: Relative residual too high -> Restart with CG"<<std::endl;
            throw( std::string("Relative residual too high -> Try slow CG muahaha.") );
         }
   }
   catch ( const std::string& err ) {
      QDPIO::cout <<"Caught "<<err<<std::endl;
      ++(*numNonConverged);
      T chi_tmp;
      (*A)(chi_tmp, chi, MINUS);
      res = InvCG2(*A, chi_tmp, psi, invParam.RsdCG, invParam.MaxCG);

#ifdef CHROMA_DO_ONE_CG_RESTART
        // Save existing n_count
      int n_count = res.n_count;

      // One automatic restart (if enabled)
      res = InvCG2(*A, chi_tmp, psi, invParam.RsdCGRestart, invParam.MaxCGRestart);
      res.n_count += n_count;
#endif
   }

//    QDPIO::cout << *numNonConverged << std::endl;

   swatch.stop();
   double time = swatch.getTimeInSeconds();

   if(didnot_compute_res){
     T r;
     r[A->subset()]=chi;
     T tmp;
     (*A)(tmp, psi, PLUS);
     r[A->subset()] -= tmp;
     res.resid = sqrt(norm2(r, A->subset()));
   }
   QDPIO::cout << "SAFE_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << std::endl;
      QDPIO::cout << "SAFE_SOLVER_TIME: "<<time<< " sec" << std::endl;



   END_CODE();

   return res;
      }


  private:
    // Hide default constructor
    LinOpSysSolverSafe() {}

    Handle< LinearOperator<T> > A;
    SysSolverSafeParams invParam;
   int* numNonConverged;
  };

} // End namespace

#endif
