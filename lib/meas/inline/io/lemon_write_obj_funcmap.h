// -*- C++ -*-
// $Id: lemon_write_obj_funcmap.h,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#ifndef __lemon_write_obj_funcmap_h__
#define __lemon_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

#ifdef QDP_USE_LEMON

namespace Chroma
{

  //! Write object function map
  /*! \ingroup inlineio */
  namespace LemonWriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function map
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_4(const string&,
			     const string&, 
			     QDP_volfmt_t, QDP_serialparallel_t),
		  void (*)(const string& buffer_id,
			   const string& filename, 
			   QDP_volfmt_t volfmt, QDP_serialparallel_t serpar),
		  StringFunctionMapError> >
    TheLemonWriteObjFuncMap;

    bool registerAll();
  }

} // end namespace Chroma

#endif // QDP_USE_LEMON

#endif // __inline_lemon_write_obj_h__ 

