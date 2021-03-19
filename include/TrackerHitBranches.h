#ifndef TrackerHitBranches_h
#define TrackerHitBranches_h 1

#include "LCTupleConf.h" 

#include "CollectionBranches.h"


class TTree ;

namespace EVENT{
  class LCCollection ;
  class LCCEvent ;
}

/** TrackerHitBranches holds branches created from LCRelations.
 * 
 * @author F. Gaede, DESY
 * @version $Id:$
 */

class TrackerHitBranches : public CollectionBranches {
  
public:
  
  TrackerHitBranches() {} ;
  
  virtual void initBranches( TTree* tree, const std::string& prefix="" ) ;
  
  virtual void fill(const EVENT::LCCollection* col, EVENT::LCEvent* evt ) ;

  //Add fill method to include individual constituents (hits)
  virtual void fill(const EVENT::LCCollection* col, const EVENT::LCCollection* hitsCol, EVENT::LCEvent* evt ) ;
  
  virtual ~TrackerHitBranches() {} ;
  

private:
  
  int    _ntrh    {} ;
  int    _thori[ LCT_TRACKERHIT_MAX ]  {} ;
  int    _thci0[ LCT_TRACKERHIT_MAX ]  {} ;
  int    _thci1[ LCT_TRACKERHIT_MAX ]  {} ;
  double _thpox[ LCT_TRACKERHIT_MAX ]  {} ;
  double _thpoy[ LCT_TRACKERHIT_MAX ]  {} ;
  double _thpoz[ LCT_TRACKERHIT_MAX ]  {} ;
  float  _thedp[ LCT_TRACKERHIT_MAX ]  {} ;
  float  _thtim[ LCT_TRACKERHIT_MAX ]  {} ;

  float  _thcov[ LCT_TRACKERHIT_MAX ][6]  {} ;

  float  _thtyp[ LCT_TRACKERHIT_MAX ]  {} ;
  float  _thqua[ LCT_TRACKERHIT_MAX ]  {} ;
  float  _thede[ LCT_TRACKERHIT_MAX ]  {} ;
  float  _thplen[ LCT_TRACKERHIT_MAX ]  {} ;
  int  _thsrc[ LCT_TRACKERHIT_MAX ]  {} ;


  int _thcidx[ LCT_SIMTRACKERHIT_MAX] {};
  int _thclen[ LCT_SIMTRACKERHIT_MAX] {};

  int _ntrc {};
  float _tcedp[ LCT_TRACKERRAWHIT_MAX ] {};
  float _tctim[ LCT_TRACKERRAWHIT_MAX ] {};
  int _tcrp0[ LCT_TRACKERRAWHIT_MAX ] {};
  int _tcrp1[ LCT_TRACKERRAWHIT_MAX ] {};

} ;

#endif



