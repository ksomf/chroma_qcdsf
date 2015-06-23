// -*- C++ -*-
//  Andrea Nobile , Luca Castagnini
/*! \file
 *  \brief Solve a M*psi=chi linear system by DDS
 */

#ifndef __syssolver_linop_dds_h__
#define __syssolver_linop_dds_h__


#include "chroma_config.h"
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_dds_params.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "dds_interface.h"


namespace Chroma
{

  //! system solver namespace
  namespace LinOpSysSolverDDSEnv
  {
    //! Register the syssolver
    bool registerAll();
  }

  inline void latticeColorMatrixToDDS(const multi1d<LatticeColorMatrix> &u){
    
    ColorMatrix matrix;
    multi1d<int> site(Nd);  
    LatticeColorMatrix shiftedU;
    int whoami = Layout::nodeNumber();
    int totSites = Layout::sitesOnNode();
    
    for(int mu = 0; mu < Nd; mu++){      
      shiftedU =  shift(u[mu], BACKWARD, mu);
      for(int i = 0 ; i < totSites; i++){
	site = Layout::siteCoords(whoami,i);
	int mu2 = (mu+1)%4; // 3 = time in dds
	// dds wants local coordinates
	int t2 = site[3]; 
	int x2 = site[0];
	int y2 = site[1]; 
	int z2 = site[2]; 
	
	if((t2+x2+y2+z2)&0x1) {	
	  matrix.elem() = (u[mu]).elem(i);
	  double * ptr = dds_uptr(t2,  x2,  y2,  z2,  mu2);	 
	  
	  for(int row = 0; row < Nc; row++) {
	    for(int col = 0; col < Nc; col++){
	      Complex val = peekColor(matrix,row,col);
	      ptr[6*row + 2*col] = toDouble(real(val));
	      ptr[6*row + 2*col + 1] = toDouble(imag(val));  
	    }
	  }
	  // backward
	  ptr+=18;
	  // takes backward from the shifted lattice
	  matrix.elem() = (shiftedU).elem(i);
	  
	  for(int row = 0; row < Nc; row++) {
	    for(int col = 0; col < Nc; col++){
	      Complex val = peekColor(matrix,row,col);
	      ptr[6*row + 2*col] = toDouble(real(val));
	      ptr[6*row + 2*col + 1] = toDouble(imag(val));
	    } 
	  }
	} //if odd
      } //sites
      
    } //mu
    
  }
  
  inline void latticeFermionToDDS(const LatticeFermion &latt_ferm ){   
    multi1d<int> site(Nd);
    int whoami = Layout::nodeNumber();
    int totSites = Layout::sitesOnNode();
    ColorVector cv;
    Fermion ferm ;
    
    for(int i = 0 ; i < totSites; i++){   
      site = Layout::siteCoords(whoami,i);
      ferm.elem() = (latt_ferm).elem(i);
      // dds uses t,x,y,z order
      double * ptr = dds_sptr_src(site[3], site[0] , site[1], site[2]);
      
      for(int spin = 0; spin < Ns; spin++){  
	cv = peekSpin( ferm, spin );
	for(int col = 0; col < Nc; col++){
	  Complex val = peekColor( cv , col );
	  ptr[6*(spin) + 2*col] = toDouble(real(val));
	  ptr[6*(spin) + 2*col + 1] = toDouble(imag(val));
	} // for col      
      } // for spin
    } // for i
  }
  
  inline void ddsToLatticeFermion( LatticeFermion &latt_ferm ){   
    
    multi1d<int> site(Nd); // node coordinates
    ColorVector cv ;
    Fermion ferm ;
    int totSites = Layout::sitesOnNode();
    int whoami = Layout::nodeNumber();
    
    for(int i = 0 ; i < totSites; i++){ 
      site = Layout::siteCoords(whoami,i);
      if(whoami == Layout::nodeNumber(site) ) {
	// dds uses t,x,y,z order
	double * ptr = dds_sptr_dst(site[3], site[0], site[1], site[2]);
	
        for(int spin = 0; spin < Ns; spin++){
	  for(int col = 0; col < Nc; col++){			
	    double a = ptr[6*spin + 2*col];
	    double b = ptr[6*spin + 2*col + 1];
	    
	    Complex val = cmplx(Real(a),Real(b));
	    pokeColor ( cv, val , col );	
	  } // for col	
	  pokeSpin(ferm, cv, spin );	
        } // for spin
        pokeSite(latt_ferm,ferm,site);
      } // if whoami
    } // for i
    
  }
  
  
  //! Solve a M*psi=chi linear system by CG2
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverDDS : public LinOpSystemSolver<T>
  {
  public:
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverDDS(Handle< LinearOperator<T> > A_,
		      Handle< FermState<
		      LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix>
		      >
		      > state_,
		      const SysSolverDDSParams& invParam_) : 
      A(A_), invParam(invParam_) 
    {
      multi1d<int> nrow(Nd);
      nrow =  Layout::lattSize();
      
      /*
	typedef struct
	{
	int np[4]; //number of processes in each dimension t x y z
	int gl[4]; //global lattice size t x y z
	int npr[8]; //tm, tp, xm, xp, ym, yp, zm, zp MPI ranks of neigh. proc.
	int cpr[4]; //t, x, y, z  cartesian coords of this proc.
	int id; // MPI rank of this proc.
	} dds_lat_grid;
      */
      dds_lat_grid init_params;
      multi1d<int> lognrow(Nd);
      multi1d<int> cpr(Nd);
      multi1d<int> tmpcoord(Nd);
      nrow =  Layout::lattSize();
      lognrow = Layout::logicalSize () ;
      cpr = Layout::nodeCoord();
      
      init_params.gl[0] = nrow[3]; // time
      init_params.gl[1] = nrow[0]; // x
      init_params.gl[2] = nrow[1]; // y
      init_params.gl[3] = nrow[2]; // z
      
      init_params.np[0] = lognrow[3];
      init_params.np[1] = lognrow[0];
      init_params.np[2] = lognrow[1];
      init_params.np[3] = lognrow[2];
      
      init_params.cpr[0] = cpr[3]; // time
      init_params.cpr[1] = cpr[0]; // x
      init_params.cpr[2] = cpr[1]; // y
      init_params.cpr[3] = cpr[2]; // z
      
      init_params.id = Layout::nodeNumber ();
      
      // time
      tmpcoord = Layout::nodeCoord();
      tmpcoord[3] = (tmpcoord[3] + lognrow[3] -1)%lognrow[3]  ;
      init_params.npr[0] = Layout::getNodeNumberFrom(tmpcoord);
      tmpcoord = Layout::nodeCoord();
      tmpcoord[3] = (tmpcoord[3] +1)%lognrow[3];  ;
      init_params.npr[1] = Layout::getNodeNumberFrom(tmpcoord);
      // x
      tmpcoord = Layout::nodeCoord();
      tmpcoord[0] = (tmpcoord[0] + lognrow[0] -1)%lognrow[0]  ;
      init_params.npr[2] = Layout::getNodeNumberFrom(tmpcoord);
      tmpcoord = Layout::nodeCoord();
      tmpcoord[0] = (tmpcoord[0] +1)%lognrow[0];  ;
      init_params.npr[3] = Layout::getNodeNumberFrom(tmpcoord);
      // y
      tmpcoord = Layout::nodeCoord();
      tmpcoord[1] = (tmpcoord[1] + lognrow[1] -1)%lognrow[1]  ;
      init_params.npr[4] = Layout::getNodeNumberFrom(tmpcoord);
      tmpcoord = Layout::nodeCoord();
      tmpcoord[1] = (tmpcoord[1] +1)%lognrow[1];  ;
      init_params.npr[5] = Layout::getNodeNumberFrom(tmpcoord);
      // z
      tmpcoord = Layout::nodeCoord();
      tmpcoord[2] = (tmpcoord[2] + lognrow[2] -1)%lognrow[2]  ;
      init_params.npr[6] = Layout::getNodeNumberFrom(tmpcoord);
      tmpcoord = Layout::nodeCoord();
      tmpcoord[2] = (tmpcoord[2] +1)%lognrow[2];  ;
      init_params.npr[7] = Layout::getNodeNumberFrom(tmpcoord);
      
      int ndfl = invParam.DeflatedNV;
      int nkv = invParam.Nkv;
      int ncy = invParam.Ncy;
      int nmr  = invParam.Nmr;
      int blkrel = invParam.BlkRel;
      int bs[4];
      for(int i=0;i<4;i++) {
	bs[i]=invParam.BS[(i+3)%Nd];
      }


      QDPIO::cout << "DD Solver Init." << endl;
      
      dds_init(&init_params);
      dds_tune(nmr, ncy, nkv, ndfl, blkrel, bs); 

      double kappa = toDouble(invParam.Kappa); 
      double csw = toDouble(invParam.Csw);
 
      if(invParam.slrc) {
	thin_fs  = new PeriodicFermState<T,P,Q>( state_.cast< SLICFermState<T, P, Q> >()->getThinLinks() );
      	latticeColorMatrixToDDS((thin_fs->getLinks()));
	QMP_barrier();
	dds_clover(kappa, csw);
	dds_gauge_updated();
	QMP_barrier();
	
	latticeColorMatrixToDDS((state_->getLinks()));
	QMP_barrier();
      } else {
	
	latticeColorMatrixToDDS((state_->getLinks()));
	QMP_barrier();
	dds_clover(kappa, csw);
	dds_gauge_updated();
	QMP_barrier();
      }
      

      
      QDPIO::cout << " Loading gauge field in DD Solver done." << endl;
      /*
      QDPIO::cout << " Dirac operator test..." << endl;
      double kappa = toDouble(invParam.Kappa); 
      double csw = toDouble(invParam.Csw); 


      QDPIO::cout << " Dirac test done: OK! " << endl;
      */
    }
    
    //! Destructor is automatic
    ~LinOpSysSolverDDS() { 
      //dds_terminate();

    }
    
    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}
    
    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const {
      
	START_CODE();	
	SystemSolverResults_t res;  // initialized by a constructor
	StopWatch swatch;
	swatch.reset();
	swatch.start();
	
	
	double kappa = toDouble(invParam.Kappa); 
	double csw = toDouble(invParam.Csw); 
	double tol = toDouble(invParam.RsdDDS); 
	int nmx = invParam.MaxDDS; 
	
	double residual = 0.;
	
	LatticeFermion tmp ;
	//change basis (different gamma rep)
	tmp = Gamma(2)*chi;
	
	latticeFermionToDDS(tmp);
	
	//solve
	int iter_count;
	iter_count  = dds_solve( kappa, csw, tol, nmx, &residual); 
	res.n_count = iter_count;
	res.resid = residual;	
	
	ddsToLatticeFermion(tmp);
	//change basis back
	psi = Gamma(2)*tmp;
	
	swatch.stop();
	double time = swatch.getTimeInSeconds();
	QDPIO::cout << "DD Solver total time: "<<time<< " sec, res: " << residual << " " << endl;
	
	  
	{ // test residual
	  T r;
	  r[A->subset()]=chi;
	  T tmp;
	  (*A)(tmp, psi, PLUS);
	  r[A->subset()] -= tmp;
	  res.resid = sqrt(norm2(r, A->subset()));
	}
	if((residual/toDouble(res.resid)) < 0.2) {
	  QDPIO::cout << "DDS RESIDUAL CHECK FAILED :( DD: " << residual << "CHK: " << res.resid << endl; 	  
	  exit(1);
	} else {
	  QDPIO::cout << "DDS_SOLVER_CHECK: "<< res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
	}
	END_CODE();
	return res;
    }
    
    
  private:
    // Hide default constructor
    LinOpSysSolverDDS() {}
    Handle< FermState<T,P,Q> > thin_fs;
    Handle< LinearOperator<T> > A;
    SysSolverDDSParams invParam;
  };
  
} // End namespace

#endif 

