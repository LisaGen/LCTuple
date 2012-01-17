#ifndef LCTupleConf_h
#define LCTupleConf_h 1


#include "LCRTRelations.h"
//==== define user extension for index in collection

struct CollIndex : public lcrtrel::LCIntExtension<CollIndex> {} ;

// ===== define constants needed for LCTuple =================


// =================================================================
//
//     maximum number of elements in colum wise ntuple
//
//         => ADJUST AS NEEDED !!!!!!
// =================================================================


#define LCT_MCPARTICLE_MAX        1000000
#define LCT_RECOPARTICLE_MAX       100000
#define LCT_TRACK_MAX              100000
#define LCT_TRACKSTATE_MAX        1000000
#define LCT_CLUSTER_MAX            100000
#define LCT_RELATION_MAX          1000000
#define LCT_SIMTRACKERHIT_MAX     1000000
#define LCT_SIMCALORIMETERHIT_MAX 1000000
#define LCT_PARTICLEID_MAX        1000000

#define LCT_STRING_MAX       1024 



//-------------------------------------------------------------------


#endif
