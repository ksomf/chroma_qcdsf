// -*- C++ -*-
// $Id: inline_barspec_w.h,v 3.5 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations
 */

#ifndef __inline_barspec_qcdsf_h__
#define __inline_barspec_qcdsf_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma
{
  /*! \ingroup inlinehadron */
  namespace InlineBarSpecEnvQCDSF
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineBarspecParamsQCDSF
  {
    InlineBarspecParamsQCDSF();
    InlineBarspecParamsQCDSF(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      bool time_rev;           // Use time reversal in baryon spectroscopy
      bool fwdbwd_average;           // Whether to average the forward and backward baryon 2pt functions
      bool fwdbwd_used;           //the forward backward flag is present in the input xml
      int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
      bool avg_equiv_mom;      // average over equivalent momenta
      bool lime;
      bool xml;
    } param;

    struct NamedObject_t
    {
      std::string  gauge_id;           /*!< Input gauge field */

      struct Props_t
      {
	std::string  first_id;
	std::string  second_id;
	std::string  third_id;
	bool haveThird;
      };

      multi1d<Props_t> sink_pairs;
    } named_obj;

    std::string lime_file;  // QCDSF format
    std::string lime_trev_file;  // QCDSF format
    std::string xml_file;  // Alternate XML file pattern
  };

  void read(XMLReader& xml, const std::string& path, InlineBarspecParamsQCDSF::Param_t& param);

  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlineBarspecQCDSF : public AbsInlineMeasurement
  {
  public:
    ~InlineBarspecQCDSF() {}
    InlineBarspecQCDSF(const InlineBarspecParamsQCDSF& p) : params(p) {}
    InlineBarspecQCDSF(const InlineBarspecQCDSF& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out);

  protected:
    //! Do the measurement
    void func_xml(const unsigned long update_no,
		  XMLWriter& xml_out);

    void func_lime(const unsigned long update_no,
		   std::string& lime_file);

  private:
    InlineBarspecParamsQCDSF params;
  };

};

#endif
