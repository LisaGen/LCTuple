#include "TrackerHitBranches.h"
#include "Exceptions.h"

#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/MCParticle.h"
#include "EVENT/SimTrackerHit.h"
#include "streamlog/streamlog.h"

#include "TTree.h"


void TrackerHitBranches::initBranches( TTree* tree, const std::string& pre ){

  if( tree == 0 ){

    throw lcio::Exception("  TrackerHitBranches::initBranches - invalid tree pointer !!! " ) ;
  }

  if (_writeparameters) CollectionBranches::initBranches(tree, (pre+"st").c_str());

  // -- Tracker hits (clusters of)
  tree->Branch( (pre+"ntrh").c_str() , &_ntrh , (pre+"ntrh/I").c_str() ) ;

  tree->Branch( (pre+"thori").c_str() , _thori , (pre+"thori["+pre+"ntrh]/I").c_str() ) ;
  tree->Branch( (pre+"thci0").c_str() , _thci0 , (pre+"thci0["+pre+"ntrh]/I").c_str() ) ;
  tree->Branch( (pre+"thci1").c_str() , _thci1 , (pre+"thci1["+pre+"ntrh]/I").c_str() ) ;
  tree->Branch( (pre+"thpox").c_str() , _thpox , (pre+"thpox["+pre+"ntrh]/D").c_str() ) ;
  tree->Branch( (pre+"thpoy").c_str() , _thpoy , (pre+"thpoy["+pre+"ntrh]/D").c_str() ) ;
  tree->Branch( (pre+"thpoz").c_str() , _thpoz , (pre+"thpoz["+pre+"ntrh]/D").c_str() ) ;
  tree->Branch( (pre+"thedp").c_str() , _thedp , (pre+"thedp["+pre+"ntrh]/F").c_str() ) ;
  tree->Branch( (pre+"thtim").c_str() , _thtim , (pre+"thtim["+pre+"ntrh]/F").c_str() ) ;
  
  tree->Branch( (pre+"thcov").c_str() , _thcov , (pre+"thcov["+pre+"ntrh][6]/F").c_str() ) ;

  tree->Branch( (pre+"thtyp").c_str() , _thtyp , (pre+"thtyp["+pre+"ntrh]/F").c_str() ) ;
  tree->Branch( (pre+"thqua").c_str() , _thqua , (pre+"thqua["+pre+"ntrh]/F").c_str() ) ;
  tree->Branch( (pre+"thede").c_str() , _thede , (pre+"thede["+pre+"ntrh]/F").c_str() ) ;
  tree->Branch( (pre+"thplen").c_str(), _thplen, (pre+"thplen["+pre+"ntrh]/F").c_str()) ;
  tree->Branch( (pre+"thsrc").c_str(), _thsrc, (pre+"thsrc["+pre+"ntrh]/I").c_str()) ;

  // index of rawHits() constituents (if stored); loop over [_thidx0;thclen-1]
  tree->Branch( (pre+"thcidx").c_str() , _thcidx , (pre+"thcidx["+pre+"ntrh]/I").c_str() ) ;
  tree->Branch( (pre+"thclen").c_str() , _thclen , (pre+"thclen["+pre+"ntrh]/I").c_str() ) ;

  // -- Tracker rawHits constutuents (individual pixels/strips)
  tree->Branch( (pre+"ntrc").c_str() , &_ntrc , (pre+"ntrc/I").c_str() ) ;

  tree->Branch( (pre+"tcedp").c_str() , _tcedp , (pre+"tcedp["+pre+"ntrc]/F").c_str() ) ;
  tree->Branch( (pre+"tctim").c_str() , _tctim , (pre+"tctim["+pre+"ntrc]/F").c_str() ) ;
  //relative position within the sensitive element in units of segmentation, up to two local coordinates
  tree->Branch( (pre+"tcrp0").c_str() , _tcrp0 , (pre+"tcrp0["+pre+"ntrc]/I").c_str() ) ;
  tree->Branch( (pre+"tcrp1").c_str() , _tcrp1 , (pre+"tcrp1["+pre+"ntrc]/I").c_str() ) ;
  
}
  

void TrackerHitBranches::fill(const EVENT::LCCollection* col, EVENT::LCEvent* evt ){
  //call the other method without a hits collection
  fill(col, 0, evt);
}

void TrackerHitBranches::fill(const EVENT::LCCollection* col, const EVENT::LCCollection* hitsCol, EVENT::LCEvent* evt ){
  

  streamlog_out( DEBUG ) << " TrackerHitBranches::fill called ... (col: " << col << ", hits col: " << hitsCol << ")" << std::endl ;

  if( !col ) return ;


  if( col->getTypeName() != lcio::LCIO::TRACKERHIT && col->getTypeName() != lcio::LCIO::TRACKERHITPLANE && col->getTypeName() != lcio::LCIO::TRACKERHITZCYLINDER ){

    std::string exStr("TrackerHitBranches::fill: invalid collection type : " ) ;

    throw EVENT::Exception( exStr + col->getTypeName() ) ; 
  }

  if (_writeparameters) CollectionBranches::fill(col, evt);

  _ntrh  = col->getNumberOfElements() ;
  _ntrc = 0; //updated below
  
  for(int i=0 ; i < _ntrh ; ++i){
    
    lcio::TrackerHit* hit = static_cast<lcio::TrackerHit*>( col->getElementAt(i) ) ;

    _thori[i] = hit->ext<CollID>();

    _thci0[i] = hit->getCellID0() ;
    _thci1[i] = hit->getCellID1() ;
    _thpox[i] = hit->getPosition()[0] ;
    _thpoy[i] = hit->getPosition()[1] ;
    _thpoz[i] = hit->getPosition()[2] ;
    _thedp[i] = hit->getEDep() ;
    _thtim[i] = hit->getTime() ;

    for(int j=0;j<6;++j){
      _thcov[ i ][ j ] = hit->getCovMatrix()[j] ;
    }

    _thtyp[ i ] = hit->getType()      ;
    _thqua[ i ] = hit->getQuality()   ;
    _thede[ i ] = hit->getEDepError() ;
    
    const lcio::LCObjectVec &rawHits = hit->getRawHits();
    
    _thsrc [ i ] = -1;
    _thplen[ i ] = -1;
    _thcidx[ i ] = _ntrc;
    float maxEDep = -1;
    float sourcePoll[3] = {0.0, 0.0, 0.0}; //see which energy deposit is dominant
    for (size_t j=0; j<rawHits.size() and (j+_ntrc) < LCT_TRACKERRAWHIT_MAX; ++j) {
      lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
      if (hitConstituent) {
        _tcedp[j+_ntrc] = hitConstituent->getEDep();
        _tctim[j+_ntrc] = hitConstituent->getTime();
        if (hitConstituent->getEDep() > maxEDep) {
          maxEDep = hitConstituent->getEDep();
          _thplen[i] = hitConstituent->getPathLength(); //truth-based
        }
        if (hitConstituent->isOverlay()) {
          sourcePoll[2] += hitConstituent->getEDep();
        } else if (hitConstituent->isProducedBySecondary()) {
          sourcePoll[1] += hitConstituent->getEDep();
        } else {
          sourcePoll[0] += hitConstituent->getEDep();
        }
        // missing: incidence angle, momentum

        //compute relative position
        const double *localPos = hitConstituent->getPosition();
        _tcrp0[j+_ntrc] = localPos[0];
        _tcrp1[j+_ntrc] = localPos[1];
      }
      //check which energy deposit dominates
      if ((sourcePoll[0] > sourcePoll[1]) && (sourcePoll[0] > sourcePoll[2])) {
        //dominated by prompt particle
        _thsrc [ i ] = 0;
      } else if (sourcePoll[1] > sourcePoll[2]) {
        //dominated by secondary particle
        _thsrc [ i ] = 1; 
      } else if (sourcePoll[2] > 0) {
        //dominated by Overlay
        _thsrc [ i ] = 2;
      }
    }
    //update total number of constituents
    _ntrc += rawHits.size();
    _thclen [ i ] = rawHits.size();
  }
}






















