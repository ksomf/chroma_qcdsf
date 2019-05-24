// -*- C++ -*-
// $Id: inline_hadspec_w.h,v 3.5 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations
 */

#ifndef __inline_messpec_qcdsf_h__
#define __inline_messpec_qcdsf_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMesSpecEnvQCDSF
  {
    extern const std::string name;
    bool registerAll();
  }

  namespace InlineMesSpecEnvQCDSFsmall
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMesSpecParamsQCDSF
  {
    InlineMesSpecParamsQCDSF();
    InlineMesSpecParamsQCDSF(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
      bool avg_equiv_mom;      // average over equivalent momenta
      bool xml;
      bool lime;
      bool linkops;
      bool traceflag;
    } param;

    struct NamedObject_t
    {
      std::string  gauge_id;           /*!< Input gauge field */

      struct Props_t
      {
	std::string  first_id;
	std::string  second_id;
      };

      multi1d<Props_t> sink_pairs;
    } named_obj;

    std::string lime_file;  // Alternate XML file pattern
    std::string xml_file;  // Alternate XML file pattern
  };

  void read(XMLReader& xml, const std::string& path, InlineMesSpecParamsQCDSF::Param_t& param);
  void read(XMLReader& xml, const std::string& path, InlineMesSpecParamsQCDSF::NamedObject_t& input);

  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlineMesSpecQCDSF : public AbsInlineMeasurement 
  {
  public:
    ~InlineMesSpecQCDSF() {}
    InlineMesSpecQCDSF(const InlineMesSpecParamsQCDSF& p) : params(p) {}
    InlineMesSpecQCDSF(const InlineMesSpecQCDSF& p) : params(p.params) {}

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
    InlineMesSpecParamsQCDSF params;
  };

  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlineMesSpecQCDSFsmall : public AbsInlineMeasurement 
  {
  public:
    ~InlineMesSpecQCDSFsmall() {}
    InlineMesSpecQCDSFsmall(const InlineMesSpecParamsQCDSF& p) : params(p) {}
    InlineMesSpecQCDSFsmall(const InlineMesSpecQCDSFsmall& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func_lime(const unsigned long update_no,
		   std::string& lime_file); 

  private:
    InlineMesSpecParamsQCDSF params;
  };

};

#endif
