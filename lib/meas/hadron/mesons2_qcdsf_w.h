// -*- C++ -*-
// $Id: mesons2_w.h,v 1.1 2006-07-10 19:53:36 edwards Exp $
/*! \file
 *  \brief Meson 2-pt functions
 */

#ifndef __mesons2_qcdsf_h__
#define __mesons2_qcdsf_h__

namespace Chroma {

  //! Meson 2-pt functions
  /* This routine is specific to Wilson fermions!
   *
   * Construct meson propagators and writes in COMPLEX 
   *
   * The two propagators can be identical or different.
   *
   * \param quark_prop_1  first quark propagator ( Read )
   * \param quark_prop_2  second (anti-) quark propagator ( Read )
   * \param t0            timeslice coordinate of the source ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   * \param xml           xml file object ( Write )
   * \param xml_group     string used for writing xml data ( Read )
   *
   *        ____
   *        \
   * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
   *        /
   *        ----
   *          x
   */


  struct Mesons_corr_QCDSF_t
  {
    //multi1d<int>     insert_mom;
    multi1d<Complex>       correlator;
  };
  

  struct Mesons_mom_QCDSF_t
  {
    multi1d<Mesons_corr_QCDSF_t >     momentum;
  };


  struct Mesons_gamma2_QCDSF_t
  {
    multi1d<Mesons_mom_QCDSF_t > gamma_value;
  };

  struct MesonsQCDSF_t
  {
    multi1d<Mesons_gamma2_QCDSF_t > gamma_value;
  };



  void write(BinaryWriter& bin, const Mesons_corr_QCDSF_t& mes);
  void write(BinaryWriter& bin, const Mesons_mom_QCDSF_t& mes);
  void write(BinaryWriter& bin, const Mesons_gamma2_QCDSF_t& mes);
  void write(BinaryWriter& bin, const MesonsQCDSF_t& mesons);

  // void read(BinaryReader& bin,  Mesons_corr_QCDSF_t<T>& mes);
  // void read(BinaryReader& bin,  Mesons_mom_QCDSF_t<T>& mes);
  // void read(BinaryReader& bin,  Mesons_gamma2_QCDSF_t<T>& mes);
  // void read(BinaryReader& bin,  MesonsQCDSF_t<T>& mesons);



  void mesons2qcdsf(const LatticePropagator& quark_prop_1,
		    const LatticePropagator& quark_prop_2,
		    const SftMom& phases,
		    int t0,
		    MesonsQCDSF_t& mesons);


  void concur2qcdsf(const multi1d<LatticeColorMatrix>& u,
		    const LatticePropagator& quark_prop_1,
		    const LatticePropagator& quark_prop_2,
		    const SftMom& phases,
		    int t0,
		    MesonsQCDSF_t& mesons);

  void mesons2qcdsfsmall(const LatticePropagator& quark_prop_1,
			 const LatticePropagator& quark_prop_2,
			 const SftMom& phases,
			 int t0,
			 Mesons_gamma2_QCDSF_t& mesons);





  // template<typename T>
  // void read(BinaryReader& bin,  Mesons_corr_QCDSF_t<T>& mes)
  // {
  //   read( bin , mes.insert_mom );
  //   read( bin , mes.correlator );
  // }
  // template<typename T>
  // void read(BinaryReader& bin,  Mesons_mom_QCDSF_t<T>& mes)
  // {
  //   read( bin , mes.momentum );
  // }
  // template<typename T>
  // void read(BinaryReader& bin,  Mesons_gamma2_QCDSF_t<T>& mes)
  // {
  //   read( bin , mes.gamma_value );
  // }
  // template<typename T>
  // void read(BinaryReader& bin,  MesonsQCDSF_t<T>& mesons)
  // {
  //   read( bin , mesons.gamma_value );
  // }







  
}  // end namespace Chroma

#endif
