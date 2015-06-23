#include "meas/hadron/rms_w.h"
#include "util/ft/sftmom.h"
#include "util/ferm/transf.h"
#include <set>

namespace Chroma
{

  //so far distance defined from the origin
  struct DistanceFunc : public SetFunc
  {
    DistanceFunc(multi1d<int> &dirs, multi1d<int> &sloc){
      lat_dir.resize(Nd);
      lat_dir = dirs;
      srcloc.resize(Nd);
      srcloc = sloc;
      gen_dists();
    }

    int operator()(const multi1d<int>& coord) const
    {
      const multi1d<int>& latt_size = Layout::lattSize();

      int xdist = abs(coord[lat_dir[0]]-srcloc[lat_dir[0]]);
      int ydist = abs(coord[lat_dir[1]]-srcloc[lat_dir[1]]);
      int zdist = abs(coord[lat_dir[2]]-srcloc[lat_dir[2]]);
      if(xdist >latt_size[lat_dir[0]]/2) xdist = latt_size[lat_dir[0]]-xdist;
      if(ydist >latt_size[lat_dir[1]]/2) ydist = latt_size[lat_dir[1]]-ydist;
      if(zdist >latt_size[lat_dir[2]]/2) zdist = latt_size[lat_dir[2]]-zdist;
      int dist = xdist*xdist+ydist*ydist+zdist*zdist;
      if(dist > max) dist = max;

      return dist;
    }

    // The number of subsets is the length of the lattice
    // in direction mu
    int numSubsets()const {return max+1;}

    int numDists()const {return dist_num.size();}

    int max;
    multi1d<int> lat_dir;
    multi1d<int> srcloc;
	std::map<int,int> dist_num;

    void gen_dists(){

      const multi1d<int>& latt_size = Layout::lattSize();

      max= 3*(latt_size[lat_dir[0]]/2)*(latt_size[lat_dir[0]]/2);

      for (int d1 = -latt_size[lat_dir[0]]/2 ; d1 < latt_size[lat_dir[0]]/2 ; d1++){
	for (int d2 = -latt_size[lat_dir[1]]/2 ; d2 < latt_size[lat_dir[1]]/2 ; d2++){
	  for (int d3 = -latt_size[lat_dir[2]]/2 ; d3 < latt_size[lat_dir[2]]/2 ; d3++){
	    int d = d1*d1+d2*d2+d3*d3;
	    if(dist_num.count(d)==0){
	      dist_num.insert(std::pair<int,int>(d,1));
	    }else{
	      int num = dist_num[d];
	      num++;
	      dist_num.erase(d);
	      dist_num.insert(std::pair<int,int>(d,num));
	    }
	  }
	}
      }
    }


  };

  void write(XMLWriter &xml_out, const std::string xml_group, DistanceFunc &dist_func,
	     multi1d<Real> & wave){

    push(xml_out,xml_group);

    int cnt = 0;
	std::map<int,int>::iterator it;
    for(it=dist_func.dist_num.begin(); it != dist_func.dist_num.end(); it++){
      int dist = (*it).first;
      Double r = sqrt(Double(dist));
      write(xml_out,"dist",r);
      Double wv = wave[cnt];
      write(xml_out,"wave",wv);
      cnt++;
    }
    pop(xml_out);
  }

  void norm_dist(DistanceFunc &dist_func, multi1d<Real> & wave){
    int cnt = 0;
	std::map<int,int>::iterator it;
    for(it=dist_func.dist_num.begin(); it != dist_func.dist_num.end(); it++){
      int dist = (*it).first;
      wave[cnt] = wave[cnt]/Double((*it).second);
      cnt++;
    }
  }


  void rms(const LatticePropagator& propagator, const int j_decay,
	   multi1d<int>& srcloc, const bool psi, const bool psi_dagger_psi,
	   XMLWriter& xml_out, const std::string& xml_group){

    START_CODE();


    const multi1d<int>& latt_size = Layout::lattSize();

    int cnt = 0;
    multi1d<int> lat_dir(Nd-1);
    for(int mu=0; mu < Nd; mu++){
      if(mu != j_decay){
	lat_dir[cnt]=mu;
	cnt++;
      }
    }

    DistanceFunc dist_func(lat_dir, srcloc);

    Set distance;
    distance.make(dist_func);

    LatticeComplex psi_wave = zero;
    LatticeComplex qdens = zero;

    for (int color_source = 0; color_source < Nc; color_source++){
      int spin_source = 0;

      LatticeFermion chi;

      PropToFerm(propagator, chi, color_source, spin_source);

      int spin_sink = 0;

      LatticeColorVector cv = peekSpin(chi,spin_sink);

      Real norm = sqrt(norm2(cv));
      cv = cv / norm;

      qdens += localInnerProduct(cv, cv);

      for(int color_sink = 0; color_sink < Nc; color_sink++){
	psi_wave += peekColor(cv,color_sink);
      }
    }

    qdens /= Double(Nc);

    push(xml_out,xml_group);


    multi1d<Real> wave_pdp(dist_func.numDists());
    wave_pdp = 0.;

    Real ms = zero;

    cnt = 0;
	std::map<int,int>::iterator it;

    for(it=dist_func.dist_num.begin(); it != dist_func.dist_num.end(); it++){
      Subset s = distance[(*it).first];
      wave_pdp[cnt]= real(sum(qdens,s));
      ms += ((*it).first)* wave_pdp[cnt];
      cnt++;
    }

    QDPIO::cout << "RMS: " << sqrt(ms) << std::endl;
    write(xml_out, "rm_x2", sqrt(ms));

    if(psi_dagger_psi){

      norm_dist(dist_func, wave_pdp);

      Double norm = 0.;
      for(int c = 0; c < wave_pdp.size(); c++){
	norm += wave_pdp[c];
      }

      wave_pdp /= norm;

      write(xml_out,"psi_dagger_psi",dist_func,wave_pdp);
    }

    if(psi){

      multi1d<Real> wave_p(dist_func.numDists());
      wave_p = 0.;

      cnt = 0;
	  std::map<int,int>::iterator it;

      for(it=dist_func.dist_num.begin(); it != dist_func.dist_num.end(); it++){
	Subset s = distance[(*it).first];
	wave_p[cnt]= real(sum(psi_wave,s));
	cnt++;
      }

      norm_dist(dist_func, wave_p);

      Double norm = 0.;
      for(int c = 0; c < wave_p.size(); c++){
	norm += wave_p[c];
      }

      wave_p /= norm;

      write(xml_out,"psi",dist_func,wave_p);
    }

    pop(xml_out);

    END_CODE();

  }

}
