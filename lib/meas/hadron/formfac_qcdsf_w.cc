// $Id: formfac_w.cc,v 3.1 2006/05/05 03:07:20 edwards Exp $
/*! \file
 *  \brief Form-factors
 *
 *  Form factors constructed from a quark and a sequential quark propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/formfac_qcdsf_w.h"
#include "meas/smear/displace.h"


namespace Chroma
{

  /*!
   * Structures for hadron parts
   *
   * \ingroup hadron
   *
   * @{
   */

  // Read a momenta struct
  // void read(BinaryReader& bin, FormFac_momentaQCDSF_t& mom)
  // {
  //   read(bin, mom.magic);
  //   if (mom.magic != 20301)
  //   {
  //     QDPIO::cerr << "read(FormFac_momenta_t): magic number invalid" << std::endl;
  //     QDP_abort(1);
  //   }
  //   read(bin, mom.inser_mom);
  //   read(bin, mom.local_current);
  //   read(bin, mom.nonlocal_current);
  // }
  // void read(BinaryReader& bin, FormFac_insertionQCDSF_t& mes)
  // {
  //   read(bin, mes.gamma_value);
  //   read(bin, mes.momenta);
  // }
  // void read(BinaryReader& bin, FormFac_insertionsQCDSF_t& form)
  // {
  //   read(bin, form.output_version);
  //   read(bin, form.formFac);
  // }




  void write(BinaryWriter& bin, const FormFac_momentaQCDSF_t& mom)
  {
    //int magic = 20301;
    //write(bin, magic);
    //write(bin, mom.inser_mom);
    write( bin , mom.local_current , mom.local_current.size() );      // at the moment only local current support
    //write(bin, mom.nonlocal_current);
  }
  void write(BinaryWriter& bin, const FormFac_insertionQCDSF_t& mes)
  {
    // write(bin, mes.gamma_value); // already written in XML
    write( bin , mes.momenta , mes.momenta.size() );
  }
  void write(BinaryWriter& bin, const FormFac_insertionsQCDSF_t& form)
  {
    //write(bin, form.output_version); // already written in XML
    write( bin , form.formFac , form.formFac.size() );
  }




  //! Compute contractions for current insertion 3-point functions.
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * \param form               structures holding formfactors ( Write )
   * \param u                  gauge fields (used for non-local currents) ( Read )
   * \param quark_propagator   quark propagator ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   * \param gamma_insertion    extra gamma insertion at source ( Read )
   * \param phases             fourier transform phase factors ( Read )
   * \param t0                 cartesian coordinates of the source ( Read )
   */

  void FormFacQCDSF(FormFac_insertionsQCDSF_t& form,
	       const multi1d<LatticeColorMatrix>& u,
	       const LatticePropagator& quark_propagator,
	       const LatticePropagator& seq_quark_prop,
	       int gamma_insertion,
	       const SftMom& phases,
	       int t0)
  {
    START_CODE();

    // Length of lattice in j_decay direction and 3pt correlations fcns
    int length = phases.numSubsets();

    int G5 = Ns*Ns-1;

    // Construct the anti-quark propagator from the seq. quark prop.
    LatticePropagator anti_quark_prop = Gamma(G5) * seq_quark_prop * Gamma(G5);

    // Rough timings (arbitrary units):
    //   Variant 1: 120
    //   Variant 2: 140
    // See previous cvs versions (before 1.10) for Variant 2 - only keeping Variant 1

    form.formFac.resize(Nd*Nd);

    // Loop over gamma matrices of the insertion current of insertion current
    for(int gamma_value = 0; gamma_value < Nd*Nd; ++gamma_value)
    {
      //  For the case where the gamma value indicates we are evaluating either
      //  the vector or axial vector currents, we will also evaluate
      //  the non-local currents.  The non-local vector current is the conserved
      //  current.  The non-local axial vector current would be partially
      //  conserved but for the Wilson term.  In these cases we will set
      //  mu = corresponding direction.  In all other cases, we will set mu = -1.

      bool compute_nonlocal;
      int mu;

      switch(gamma_value){
      case  1:
      case 14:
	mu = 0;
	compute_nonlocal = true;
	break;
      case  2:
      case 13:
	mu = 1;
	compute_nonlocal = true;
	break;
      case  4:
      case 11:
	mu = 2;
	compute_nonlocal = true;
	break;
      case  8:
      case  7:
	mu = 3;
	compute_nonlocal = true;
	break;
      default:
	mu = -1;
	compute_nonlocal = false;
      }

      // The local non-conserved vector-current matrix element
      LatticeComplex corr_local_fn =
	trace(adj(anti_quark_prop) * Gamma(gamma_value) * quark_propagator * Gamma(gamma_insertion));

      multi2d<DComplex> hsum, hsum_nonlocal;
      hsum = phases.sft(corr_local_fn);

      // Construct the non-local current matrix element
      //
      // The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
      //                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
      // NOTE: the 1/2  is included down below in the sumMulti stuff
      LatticeComplex corr_nonlocal_fn;
      if(compute_nonlocal){
	corr_nonlocal_fn =
	  trace(adj(u[mu] * shift(anti_quark_prop, FORWARD, mu)) *
		(quark_propagator + Gamma(gamma_value) * quark_propagator) * Gamma(gamma_insertion));
	LatticePropagator tmp_prop1 = u[mu] *
	  shift(quark_propagator, FORWARD, mu);
	corr_nonlocal_fn -= trace(adj(anti_quark_prop) *
				  (tmp_prop1 - Gamma(gamma_value) * tmp_prop1) * Gamma(gamma_insertion));

	hsum_nonlocal = phases.sft(corr_nonlocal_fn);
      }


      form.formFac[gamma_value].gamma_value = gamma_value;
      form.formFac[gamma_value].momenta.resize(phases.numMom());  // hold momenta output

      // Loop over insertion momenta and print out results
      for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num)
      {
	form.formFac[gamma_value].momenta[inser_mom_num].inser_mom = phases.numToMom(inser_mom_num);

	multi1d<Complex> local_cur3ptfn(length); // always compute
	multi1d<Complex> nonlocal_cur3ptfn;
	if (compute_nonlocal)
	  nonlocal_cur3ptfn.resize(length);      // possibly compute

	for (int t=0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;

	  local_cur3ptfn[t_eff] = Complex(hsum[inser_mom_num][t]);
	  if (compute_nonlocal)
	    nonlocal_cur3ptfn[t_eff] = 0.5 * Complex(hsum_nonlocal[inser_mom_num][t]);

	} // end for(t)

	form.formFac[gamma_value].momenta[inser_mom_num].local_current    = local_cur3ptfn;
	form.formFac[gamma_value].momenta[inser_mom_num].nonlocal_current = nonlocal_cur3ptfn;

      } // end for(inser_mom_num)
    } // end for(gamma_value)

    END_CODE();
  }








  void FormFac1DQCDSF(FormFac_insertionsQCDSF_t& form,
		 const multi1d<LatticeColorMatrix>& u,
		 const LatticePropagator& quark_propagator,
		 const LatticePropagator& seq_quark_prop,
		 int gamma_insertion,
		 const SftMom& phases,
		 int t0)
  {
    QDPIO::cout << "FormFac1DQCDSF" << std::endl;

    START_CODE();

    // Length of lattice in j_decay direction and 3pt correlations fcns
    int length = phases.numSubsets();

    int G5 = Ns*Ns-1;

    // Construct the anti-quark propagator from the seq. quark prop.
    LatticePropagator seq_src = Gamma(G5) * adj(seq_quark_prop) * Gamma(G5);
    LatticePropagator seq_src_dagger = Gamma(G5) * seq_quark_prop * Gamma(G5);

    // Rough timings (arbitrary units):
    //   Variant 1: 120
    //   Variant 2: 140
    // See previous cvs versions (before 1.10) for Variant 2 - only keeping Variant 1

    form.formFac.resize(Nd*Nd*4);

    // Loop over gamma matrices of the insertion current of insertion current
    for(int nu = 0 ; nu < 4 ; ++nu ) {

      QDPIO::cout << "nu = " << nu << std::endl;

      LatticePropagator quark_prop_der = rightNabla( quark_propagator , u , nu , 1 );
      LatticePropagator seq_src_dagger_der = adj( rightNabla( seq_src_dagger , u , nu , 1 ) );

      for(int gamma_value = 0; gamma_value < Nd*Nd; ++gamma_value)
	{
	  int store = nu*Nd*Nd + gamma_value;
	  QDPIO::cout << "gamma = " << gamma_value << " store = " << store << std::endl;

	  //  For the case where the gamma value indicates we are evaluating either
	  //  the vector or axial vector currents, we will also evaluate
	  //  the non-local currents.  The non-local vector current is the conserved
	  //  current.  The non-local axial vector current would be partially
	  //  conserved but for the Wilson term.  In these cases we will set
	  //  mu = corresponding direction.  In all other cases, we will set mu = -1.


	  // The local non-conserved vector-current matrix element
	  //       LatticeComplex corr_local_fn =
	  // 	trace(adj(anti_quark_prop) * Gamma(gamma_value) * quark_propagator * Gamma(gamma_insertion));

	  LatticeComplex corr_local_fn =
	    0.25 * trace( seq_src * Gamma(gamma_value) * quark_prop_der * Gamma(gamma_insertion) -
			  seq_src_dagger_der * Gamma(gamma_value) * quark_propagator * Gamma(gamma_insertion) );

	  multi2d<DComplex> hsum, hsum_nonlocal;
	  hsum = phases.sft(corr_local_fn);

	  // Construct the non-local current matrix element
	  //
	  // The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
	  //                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
	  // NOTE: the 1/2  is included down below in the sumMulti stuff


	  form.formFac[ store ].gamma_value = store;
	  form.formFac[ store ].momenta.resize(phases.numMom());  // hold momenta output

	  // Loop over insertion momenta and print out results
	  for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num)
	    {
	      form.formFac[ store ].momenta[inser_mom_num].inser_mom = phases.numToMom(inser_mom_num);

	      multi1d<Complex> local_cur3ptfn(length); // always compute

	      for (int t=0; t < length; ++t)
		{
		  int t_eff = (t - t0 + length) % length;

		  local_cur3ptfn[t_eff] = Complex(hsum[inser_mom_num][t]);

		} // end for(t)

	      form.formFac[ store ].momenta[inser_mom_num].local_current    = local_cur3ptfn;
	      //form.formFac[gamma_value].momenta[inser_mom_num].nonlocal_current = nonlocal_cur3ptfn;

	    } // end for(inser_mom_num)
	} // end for(gamma_value)
    }

    END_CODE();
  }

  template <typename T,int N>
  void printprop(PSpinMatrix<T,N> t1,PSpinMatrix<T,N> t2) {
    for (int i=0;i<4;i++)
      for (int q=0;q<4;q++)
	for (int i2=0;i2<3;i2++)
	  for (int q2=0;q2<3;q2++)
	    QDPIO::cout << t1.elem(i,q).elem(i2,q2) << "  " << t2.elem(i,q).elem(i2,q2) << std::endl;
  }

  template<typename T>
  T DmuDnu(const multi1d<LatticeColorMatrix>& u, const T& psi0,const T& psi1,int g,int mu,int nu)
  {
    T term1 = -psi0 * Gamma(g) *
      adj(shift(u[mu],BACKWARD,mu)) * shift(u[nu],BACKWARD,mu) *
      shift(shift(psi1,BACKWARD,mu),FORWARD,nu);

    T term2 = +psi0 * Gamma(g) *
      adj(shift(u[mu],BACKWARD,mu)) *
      adj(shift(shift(u[nu],BACKWARD,mu),BACKWARD,nu)) *
      shift(shift(psi1,BACKWARD,mu),BACKWARD,nu);

    T term3 = +psi0 * Gamma(g) *
      u[mu] * shift(u[nu],FORWARD,mu) *
      shift(shift(psi1,FORWARD,mu),FORWARD,nu);

    T term4 = -psi0 * Gamma(g) *
      u[mu] * adj(shift(shift(u[nu],FORWARD,mu),BACKWARD,nu)) *
      shift(shift(psi1,FORWARD,mu),BACKWARD,nu);

    T tmp;
    tmp = 0.25 * (
		  term1 + term2 + term3 + term4 +

		  shift(term1,FORWARD,mu) + shift(term2,FORWARD,mu) +
		  shift(term3,BACKWARD,mu) + shift(term4,BACKWARD,mu) +

		  shift(term1,BACKWARD,nu) + shift(term2,FORWARD,nu) +
		  shift(term3,BACKWARD,nu) + shift(term4,FORWARD,nu) +

		  shift(shift(term1,FORWARD,mu),BACKWARD,nu) +
		  shift(shift(term1,FORWARD,mu),FORWARD,nu) +
		  shift(shift(term1,BACKWARD,mu),BACKWARD,nu) +
		  shift(shift(term1,BACKWARD,mu),FORWARD,nu)     );


    // {
    //   T t1 = term1 + term2 + term3 + term4;
    //   T t2 = psi0 * Gamma(g) * rightNabla( rightNabla( psi1 , u , nu , 1 ) , u , mu , 1 );
    //   int mis=0;
    //   for (int s = 0 ; s < Layout::sitesOnNode() ; s++)
    // 	for (int i=0;i<4;i++)
    // 	  for (int q=0;q<4;q++)
    // 	    for (int i2=0;i2<3;i2++)
    // 	      for (int q2=0;q2<3;q2++)
    // 		if ( ( t1.elem(s).elem(i,q).elem(i2,q2).real() - t2.elem(s).elem(i,q).elem(i2,q2).real() > 0.001 ) ||
    // 		     ( t1.elem(s).elem(i,q).elem(i2,q2).imag() - t2.elem(s).elem(i,q).elem(i2,q2).imag() > 0.001 ) ) {
    // 		  mis++;
    // 		  QDPIO::cout << t1.elem(s).elem(i,q).elem(i2,q2) << " != " <<  t2.elem(s).elem(i,q).elem(i2,q2) << std::endl;
    // 		}
    //   QDPIO::cout << "term1 mismatch = " << mis << std::endl;
    // }

    // {
    //   T t1 =
    // 	shift(term1,FORWARD,mu) + shift(term2,FORWARD,mu) +
    // 	shift(term3,BACKWARD,mu) + shift(term4,BACKWARD,mu);
    //   T t2 =  adj(rightNabla( psi0 , u , mu , 1 )) * Gamma(g) * rightNabla( psi1 , u , nu , 1 );
    //   int mis=0;
    //   for (int s = 0 ; s < Layout::sitesOnNode() ; s++)
    // 	for (int i=0;i<4;i++)
    // 	  for (int q=0;q<4;q++)
    // 	    for (int i2=0;i2<3;i2++)
    // 	      for (int q2=0;q2<3;q2++)
    // 		if ( ( t1.elem(s).elem(i,q).elem(i2,q2).real() - t2.elem(s).elem(i,q).elem(i2,q2).real() > 0.001 ) ||
    // 		     ( t1.elem(s).elem(i,q).elem(i2,q2).imag() - t2.elem(s).elem(i,q).elem(i2,q2).imag() > 0.001 ) ) {
    // 		  mis++;
    // 		  QDPIO::cout << t1.elem(s).elem(i,q).elem(i2,q2) << " != " <<  t2.elem(s).elem(i,q).elem(i2,q2) << std::endl;
    // 		}
    //   QDPIO::cout << "term2 mismatch = " << mis << std::endl;
    // }

    // {
    //   T t1 =
    // 	shift(shift(term1,FORWARD,mu),BACKWARD,nu) +
    // 	shift(shift(term1,FORWARD,mu),FORWARD,nu) +
    // 	shift(shift(term1,BACKWARD,mu),BACKWARD,nu) +
    // 	shift(shift(term1,BACKWARD,mu),FORWARD,nu);

    //   T t20 = adj(rightNabla( psi0 , u , nu , 1 ));
    //   T t2 =  adj(rightNabla( t20 , u , mu , 1 )) * Gamma(g) * psi1;
    //   int mis=0;
    //   for (int s = 0 ; s < Layout::sitesOnNode() ; s++)
    // 	for (int i=0;i<4;i++)
    // 	  for (int q=0;q<4;q++)
    // 	    for (int i2=0;i2<3;i2++)
    // 	      for (int q2=0;q2<3;q2++)
    // 		if ( ( t1.elem(s).elem(i,q).elem(i2,q2).real() - t2.elem(s).elem(i,q).elem(i2,q2).real() > 0.001 ) ||
    // 		     ( t1.elem(s).elem(i,q).elem(i2,q2).imag() - t2.elem(s).elem(i,q).elem(i2,q2).imag() > 0.001 ) ) {
    // 		  mis++;
    // 		  QDPIO::cout << t1.elem(s).elem(i,q).elem(i2,q2) << " != " <<  t2.elem(s).elem(i,q).elem(i2,q2) << std::endl;
    // 		}
    //   QDPIO::cout << "term4 mismatch = " << mis << std::endl;
    // }


    // for (int s = 0 ; s < 1 ; s++)
    // 	printprop(t1.elem(s),t2.elem(s));

    // 1st term
    // psi0 * Gamma(g) * rightNabla( rightNabla( psi1 , u , nu , 1 ) , u , mu , 1 );

    // // 2nd term
    // leftNabla( quark_propagator , u , mu , 1 ) * rightNabla( quark_propagator , u , nu , 1 );
    // leftNabla( quark_propagator , u , mu , 1 ) = adj(rightNabla( quark_propagator , u , mu , 1 ));

    // // 3rd term more tricky

    // // 4th term
    // leftNabla( leftNabla( quark_propagator , u , mu , 1 ) , u , nu , 1 );


    return tmp;
  }





  void FormFac2DQCDSF(FormFac_insertionsQCDSF_t& form,
		      const multi1d<LatticeColorMatrix>& u,
		      const LatticePropagator& quark_propagator,
		      const LatticePropagator& seq_quark_prop,
		      int gamma_insertion,
		      const SftMom& phases,
		      int t0)
  {
    QDPIO::cout << "FormFac2DQCDSF" << std::endl;

    START_CODE();

    // Length of lattice in j_decay direction and 3pt correlations fcns
    int length = phases.numSubsets();

    int G5 = Ns*Ns-1;

    // Construct the anti-quark propagator from the seq. quark prop.
    // LatticePropagator seq_src = Gamma(G5) * adj(seq_quark_prop) * Gamma(G5);
    // LatticePropagator seq_src_dagger = Gamma(G5) * seq_quark_prop * Gamma(G5);

    LatticePropagator anti_quark_prop = Gamma(G5) * adj(seq_quark_prop) * Gamma(G5);

    // Rough timings (arbitrary units):
    //   Variant 1: 120
    //   Variant 2: 140
    // See previous cvs versions (before 1.10) for Variant 2 - only keeping Variant 1

    form.formFac.resize(Nd*Nd*4*4);

    // Loop over gamma matrices of the insertion current of insertion current
    for(int nu = 0 ; nu < 4 ; ++nu )
      {
	for(int mu = 0 ; mu < 4 ; ++mu )
	  {

	    QDPIO::cout << "nu = " << nu << ",  mu = " << mu << std::endl;

	    for(int gamma_value = 0; gamma_value < Nd*Nd; ++gamma_value)
	      {
		LatticePropagator deriv2D = DmuDnu( u , anti_quark_prop , quark_propagator , gamma_value , mu , nu ) *
		  Gamma(gamma_insertion);

		int store = nu*Nd*Nd*4 + mu*Nd*Nd + gamma_value;
		QDPIO::cout << "gamma = " << gamma_value << " store = " << store << std::endl;

		LatticeComplex corr_local_fn = trace( deriv2D );

		multi2d<DComplex> hsum, hsum_nonlocal;
		hsum = phases.sft(corr_local_fn);

		form.formFac[ store ].gamma_value = store;
		form.formFac[ store ].momenta.resize(phases.numMom());

		for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num)
		  {
		    form.formFac[ store ].momenta[inser_mom_num].inser_mom = phases.numToMom(inser_mom_num);
		    multi1d<Complex> local_cur3ptfn(length);
		    for (int t=0; t < length; ++t)
		      {
			int t_eff = (t - t0 + length) % length;
			local_cur3ptfn[t_eff] = Complex(hsum[inser_mom_num][t]);
		      }
		    form.formFac[ store ].momenta[inser_mom_num].local_current    = local_cur3ptfn;
		  }
	      }
	  }
      }
    END_CODE();
  }







}  // end namespace Chroma
