#ifndef SimCalorimeterHitBranches_h
#define SimCalorimeterHitBranches_h 1

#include "LCTupleConf.h" 

#include "CollectionBranches.h"


class TTree ;

namespace EVENT{
  class LCCollection ;
  class LCCEvent ;
}

/** SimCalorimeterHitBranches holds branches created from LCRelations.
 * 
 * @author F. Gaede, DESY
 * @version $Id: SimCalorimeterHitBranches.h 2945 2012-01-16 15:05:02Z gaede $
 */

class SimCalorimeterHitBranches : public CollectionBranches {
  
public:
  
  SimCalorimeterHitBranches() {} ;
  
  virtual void initBranches( TTree* tree, const std::string& prefix="" ) ;
  
  virtual void fill(const EVENT::LCCollection* col, EVENT::LCEvent* evt ) ;
  
  virtual ~SimCalorimeterHitBranches() {} ;
  

private:
  
  int    _nsch    {} ;    
  
  int    _scori[ LCT_SIMCALORIMETERHIT_MAX ]  {} ;

  int    _scci0[ LCT_SIMCALORIMETERHIT_MAX ]  {} ;
  int    _scci1[ LCT_SIMCALORIMETERHIT_MAX ]  {} ;
  float  _scpox[ LCT_SIMCALORIMETERHIT_MAX ]  {} ;
  float  _scpoy[ LCT_SIMCALORIMETERHIT_MAX ]  {} ;
  float  _scpoz[ LCT_SIMCALORIMETERHIT_MAX ]  {} ;
  float  _scene[ LCT_SIMCALORIMETERHIT_MAX ]  {} ;
  int    _scmcc[ LCT_SIMCALORIMETERHIT_MAX ]  {} ;
  //float  _sctim[ LCT_SIMCALORIMETERHIT_MAX ][ 100 ]  {} ;
  //int  _scpdg[ LCT_SIMCALORIMETERHIT_MAX ][ 1000 ]  {} ;
  //int _scci0cont[ LCT_SIMCALORIMETERHIT_MAX ][ 1000 ]  {} ;
} ;

#endif



