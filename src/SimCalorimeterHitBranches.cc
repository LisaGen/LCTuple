#include "SimCalorimeterHitBranches.h"
#include "Exceptions.h"

#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include "streamlog/streamlog.h"

#include "TTree.h"


void SimCalorimeterHitBranches::initBranches( TTree* tree, const std::string& pre){

  if( tree == 0 ){

    throw lcio::Exception("  SimCalorimeterHitBranches::initBranches - invalid tree pointer !!! " ) ;
  }

  if (_writeparameters) CollectionBranches::initBranches(tree, (pre+"sc").c_str());


  tree->Branch( (pre+"nsch").c_str() , &_nsch , (pre+"nsch/I").c_str() ) ;
  
  tree->Branch( (pre+"scori").c_str() , _scori , (pre+"scori["+pre+"nsch]/I").c_str() ) ;

  tree->Branch( (pre+"scci0").c_str() , _scci0 , (pre+"scci0["+pre+"nsch]/I").c_str() ) ;
  tree->Branch( (pre+"scci1").c_str() , _scci1 , (pre+"scci1["+pre+"nsch]/I").c_str() ) ;
  tree->Branch( (pre+"scpox").c_str() , _scpox , (pre+"scpox["+pre+"nsch]/F").c_str() ) ;
  tree->Branch( (pre+"scpoy").c_str() , _scpoy , (pre+"scpoy["+pre+"nsch]/F").c_str() ) ;
  tree->Branch( (pre+"scpoz").c_str() , _scpoz , (pre+"scpoz["+pre+"nsch]/F").c_str() ) ;
  tree->Branch( (pre+"scene").c_str() , _scene , (pre+"scene["+pre+"nsch]/F").c_str() ) ;
  tree->Branch( (pre+"scmcc").c_str() , _scmcc , (pre+"scmcc["+pre+"nsch]/I").c_str() ) ;
  //tree->Branch( (pre+"sctim").c_str() , _sctim , (pre+"sctim["+pre+"nsch][100]/F").c_str() ) ; 
 
  //tree->Branch( (pre+"scpdg").c_str() , _scpdg , (pre+"scpdg["+pre+"nsch][1000]/I").c_str() ) ;  
  //tree->Branch( (pre+"scci0cont").c_str() , _scci0cont , (pre+"scci0cont["+pre+"nsch][1000]/I").c_str() ) ;
  
}


void SimCalorimeterHitBranches::fill(const EVENT::LCCollection* col, EVENT::LCEvent* evt ){
  
  if( !col ) return ;

  streamlog_out( DEBUG2 ) << " SimCalorimeterHitBranches::fill called ... " << std::endl ;

  if( col->getTypeName() != lcio::LCIO::SIMCALORIMETERHIT ){

    std::string exStr("SimCalorimeterHitBranches::fill: invalid collection type : " ) ;

    throw EVENT::Exception( exStr + col->getTypeName() ) ; 
  }

  if (_writeparameters) CollectionBranches::fill(col, evt);

  _nsch  = col->getNumberOfElements() ;

  for(int i=0 ; i < _nsch ; ++i){
    
    //std::cout<<"SimHit: " << i <<std::endl;
    
    lcio::SimCalorimeterHit* hit = static_cast<lcio::SimCalorimeterHit*>( col->getElementAt(i) ) ;

    _scori[i] = hit->ext<CollID>();
    _scci0[i] = hit->getCellID0() ;
    _scci1[i] = hit->getCellID1() ;
    _scpox[i] = hit->getPosition()[0] ;
    _scpoy[i] = hit->getPosition()[1] ;
    _scpoz[i] = hit->getPosition()[2] ;
    _scene[i] = hit->getEnergy() ;
    _scmcc[i] = hit->getNMCContributions() ;
  
    //std::cout<<"Hit PDG_ID" << hit->getPDGCont(j)<<std::endl;

    //for (int j=0; j < hit->getNMCParticles(); ++j) {
     // _sctim[i][j] = hit->getTimeCont(j);
      
      //_scMCpdg[i][j] = hit->getParticleCont(j)->getPDG();
      
      
      //std::cout<<"Hit PDG_ID" << hit->getPDGCont(j)<<std::endl;
    //} 
    //std::cout<<"***********************"<<std::endl; 

    
    /*for (int j=0; j < hit->getNMCContributions(); ++j) {
     // _sctim[i][j] = hit->getTimeCont(j);
      
      _scpdg[i][j] = hit->getPDGCont(j);
      _scci0cont[i][j] = hit->getCellID0() ;
      //std::cout<<"Hit sctim" << hit->getTimeCont(j)<<std::endl;
      std::cout<<"Hit PDG_ID" << _scpdg[i][j]<<std::endl;
    }+*/
    //std::cout<<"***********************"<<std::endl; 
  }
  	
}






















