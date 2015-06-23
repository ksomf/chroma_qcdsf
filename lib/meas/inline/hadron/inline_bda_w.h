// $Id: inline_bda_w.h, v 1.0 2011-04-11 chagen, rwschiel $

#ifndef __inline_bda_w_h__
#define __inline_bda_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineBDAEnv 
  {
    extern const std::string name;
    bool registerAll ();
  }


  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineBDAParams 
  {
    InlineBDAParams ();
    InlineBDAParams (XMLReader& xml_in, const std::string& path);
    void write (XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int version;
      int gam_diquark;
      int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
      multi2d <int> moms;      // the momenta to be calculated and written
      bool avg_equiv_mom;      // average over equivalent momenta
    } param;

    struct NamedObject_t
    {
      std::string          gauge_id;  /*!< Input gauge field */
      multi1d<std::string> prop_ids;  /*!< Input forward propagators */
    } named_obj;

    std::string bda_file ;  // binary file to write the qqq object
    std::string xml_file ;  // Alternate XML file pattern
  };


  //! Inline measurement of BDAs
  /*! \ingroup inlinehadron */
  class InlineBDA : public AbsInlineMeasurement 
  {
  public:
    ~InlineBDA () {}
    InlineBDA (const InlineBDAParams& p) : params (p) {}
    InlineBDA (const InlineBDA& p) : params (p.params) {}

    unsigned long getFrequency (void) const {return params.frequency;}

    //! Do the measurement
    void operator() (const unsigned long update_no,
		     XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func (const unsigned long update_no,
	       XMLWriter& xml_out); 

  private:
    InlineBDAParams params;
  };

};

#endif
