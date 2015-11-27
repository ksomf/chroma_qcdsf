/*! \file inline_npr_vertex_nonlocal_w.h
 *  \author Stanislav Kazmin (University Leipzig)
 *  \date 2015-11-24
 *  \brief NPR vertex calculations for the non local axial current.
 */

#ifndef __inline_npr_nonlocal_vertex_h__
#define __inline_npr_nonlocal_vertex_h__

#include "meas/inline/hadron/inline_npr_vertex_w.h"
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
/*! \ingroup inlinehadron */
namespace InlineNprVertexNonlocalEnv
{
extern const std::string name;
bool registerAll();
}

//! Parameter structure
/*! \ingroup inlinehadron */
struct InlineNprVertexNonlocalParams
{
	InlineNprVertexNonlocalParams();
	InlineNprVertexNonlocalParams(XMLReader& xml_in, const std::string& path);
	void write(XMLWriter& xml_out, const std::string& path);
	unsigned long frequency;
	//! Parameters
	struct Param_t
	{
		int          links_max;          /*!< maximum number of links */
		std::string  file_name;          /*!< bb output file name pattern */
		GroupXML_t   cfs;                /*!< Fermion state */
	} param;

	//! Propagators
	struct NamedObject_t
	{
		std::string       gauge_id;        /*!< Input Gauge id */
		std::string       prop_id;         /*!< Input forward prop */
	} named_obj;

	std::string xml_file;  // Alternate XML file pattern
};


//! Inline measurement of NPR vertices
/*! \ingroup inlinehadron */
class InlineNprVertexNonlocal : public AbsInlineMeasurement
{
	public:
		~InlineNprVertexNonlocal() {}
		InlineNprVertexNonlocal(const InlineNprVertexNonlocalParams& p) : params(p) {}
		InlineNprVertexNonlocal(const InlineNprVertexNonlocal& p) : params(p.params) {}

		unsigned long getFrequency(void) const
		{
			return params.frequency;
		}

		//! Do the measurement
		void operator()(const unsigned long update_no, XMLWriter& xml_out);

	protected:
		//! Do the measurement
		void func(const unsigned long update_no, XMLWriter& xml_out);

	private:
		InlineNprVertexNonlocalParams params;
};

};

#endif // __inline_npr_nonlocal_vertex_h__
