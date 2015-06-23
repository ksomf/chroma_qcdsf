// -*- C++ -*-
/*! \file
 * \brief Inline calculation of the RMS of a wavefunction
 *
 * Calculates the root mean square of the wavefunction (from a propagator source)
 */

#ifndef __inline_rms_h__
#define __inline_rms_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineRMSEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlineRMSParams 
  {
    InlineRMSParams();
    InlineRMSParams(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    struct Param_t
    {
      int version;
      bool psi_wfn;
      bool psi_dag_psi_wfn;
    } param;

    struct NamedObject_t
    {
      std::string     gauge_id;
      std::string     source_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline propagator calculation
  /*! \ingroup inlinehadron */
  class InlineRMS : public AbsInlineMeasurement 
  {
  public:
    ~InlineRMS() {}
    InlineRMS(const InlineRMSParams& p) : params(p) {}
    InlineRMS(const InlineRMS& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineRMSParams params;
  };

}

#endif
