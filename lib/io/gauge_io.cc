/*! \file
 * \brief Routines associated with Chroma propagator gauge IO
 */

#include "chromabase.h"
#include "io/gauge_io.h"

namespace Chroma {


// Read a Chroma propagator
/*
 * \param file_xml     xml reader holding config info ( Modify )
 * \param record_xml   xml reader holding config info ( Modify )
 * \param u            gauge configuration ( Modify )
 * \param file         path ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */
void readGauge(XMLReader& file_xml,
	       XMLReader& record_xml,
	       multi1d<LatticeColorMatrix>& u,
	       const std::string& file,
	       QDP_serialparallel_t serpar)
{
  QDPFileReader to(file_xml,file,serpar);

  /*
   * This is problematic - the size should
   * come from the read - a resize. Currently, QDPIO does not
   * support this
   */
#if 0
  multi1d<LatticeColorMatrixF> u_f(u.size());
  read(to, record_xml, u_f); // Always read in single precision
#endif
  read(to,record_xml,u);
  if (to.bad())
  {
    QDPIO::cerr << __func__ << ": error reading file " << file << std::endl;
    QDP_abort(1);
  }

#if 0
  for(int mu=0; mu < u.size(); ++mu)
    u[mu] = u_f[mu];
#endif
  close(to);
}




void readGaugeILDG(XMLReader& file_xml,
		   XMLReader& record_xml,
		   multi1d<LatticeColorMatrix>& u,
		   const std::string& file,
		   QDP_serialparallel_t serpar)
{
	QLimeReader w(file.c_str());

  bool found=false;
  while ( !found ) {
	  std::string recname (w.recordName());
    QDP::QDPIO::cout << "lime record name=" << recname << ", size=" << w.recordSize() << "\n";
    found = recname.find("ildg-data-lfn") != std::string::npos ;
    if (!found)
      if ( w.nextRecord() != QLIME_SUCCESS ) {
		  QDP::QDPIO::cout << "end of lime file" << std::endl;
	break;
      }
  }

  if (!found) {
    QDPIO::cerr << "No ILDG lfn string in gauge file found!" << std::endl;
    QDP_abort(1);
  }

  uint64_t recsize = w.recordSize();
  char * tmp = new char[ recsize ];
  w.read( (void*)(tmp) , recsize );
  std::string strLfn(tmp,recsize);
  delete[] tmp;

  XMLBufferWriter f_xml, r_xml;
  push(f_xml, "gauge");
  write(f_xml, "id", int(0));
  pop(f_xml);
  push(r_xml, "ILDG");
  write(r_xml,"LFN",strLfn);
  pop(r_xml);



  QDPFileReader to(file_xml,file,serpar);

  /*
   * This is problematic - the size should
   * come from the read - a resize. Currently, QDPIO does not
   * support this
   */
  multi1d<LatticeColorMatrixF> u_f(u.size());
  read(to,record_xml,u_f);      // Always read in single precision!

  file_xml.open(f_xml);
  record_xml.open(r_xml);

  if (to.bad())
    {
      QDPIO::cerr << __func__ << ": error reading file " << file << std::endl;
      QDP_abort(1);
    }

  for(int mu=0; mu < u.size(); ++mu)
    u[mu] = u_f[mu];

  close(to);
}


// Write a Gauge field in QIO format
/*
 * \param file_xml    xml reader holding config info ( Modify )
 * \param record_xml  xml reader holding config info ( Modify )
 * \param u           gauge configuration ( Modify )
 * \param file        path ( Read )
 * \param volfmt      either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
 * \param serpar      either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */
void writeGauge(XMLBufferWriter& file_xml,
		XMLBufferWriter& record_xml,
		const multi1d<LatticeColorMatrix>& u,
		const std::string& file,
		QDP_volfmt_t volfmt,
		QDP_serialparallel_t serpar)
{
  QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
  if (to.bad())
  {
    QDPIO::cerr << __func__ << ": error writing file " << file << std::endl;
    QDP_abort(1);
  }

#if 0
  multi1d<LatticeColorMatrixF> u_f(u.size());
  for(int mu=0; mu < u.size(); ++mu)
    u_f[mu] = u[mu];

  write(to,record_xml,u_f);         // Always save in single precision!
#endif
  write(to, record_xml, u);         // Write in native precision
  close(to);
}


}  // end namespace Chroma
