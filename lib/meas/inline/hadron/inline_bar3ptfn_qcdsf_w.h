// -*- C++ -*-
// $Id: inline_bar3ptfn_w.h,v 3.3 2007/04/18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline measurement of bar3ptfn
 *
 * Form-factors
 */

#ifndef __inline_bar3ptfn_qcdsf_h__
#define __inline_bar3ptfn_qcdsf_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineBar3ptfnEnvQCDSF
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineBar3ptfnParamsQCDSF 
  {
    InlineBar3ptfnParamsQCDSF();
    InlineBar3ptfnParamsQCDSF(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct Param_t
    {
      int              mom2_max;           /*!< (mom)^2 <= mom2_max */
      int              j_decay;
      int              deriv;
    } param;

    struct SeqProp_t
    {
      std::string      seqprop_id;
      int              gamma_insertion;    /*!< second gamma insertion */
    };

    struct NamedObject_t
    {
      std::string           gauge_id;
      std::string           prop_id;
      multi1d<SeqProp_t>    seqprops;
      std::string           bar3ptfn_file;
      std::string           bar3ptfn_1D_file;
      std::string           bar3ptfn_2D_file;
    } named_obj;
  };

  void read(XMLReader& xml, const string& path, InlineBar3ptfnParamsQCDSF::Param_t& param);


  //! Inline measurement of 3pt functions
  /*! \ingroup inlinehadron */
  class InlineBar3ptfnQCDSF : public AbsInlineMeasurement 
  {
  public:
    ~InlineBar3ptfnQCDSF() {}
    InlineBar3ptfnQCDSF(const InlineBar3ptfnParamsQCDSF& p) : params(p) {}
    InlineBar3ptfnQCDSF(const InlineBar3ptfnQCDSF& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineBar3ptfnParamsQCDSF params;
  };

};

#endif
