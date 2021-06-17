#include "RecoTracker/SpecialSeedGenerators/interface/SeedGeneratorForPixelTelescope.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"

SeedGeneratorForPixelTelescope::SeedGeneratorForPixelTelescope(const MagneticField* mf,
							   const TrackerGeometry* geom,
			       const TransientTrackingRecHitBuilder* builder,
			       const Propagator* propagatorAlong,
			       const Propagator* propagatorOpposite,
			       const std::vector<int>& charges,
			       bool momFromPSet,
			       double errorRescaling):theMagfield(mf), theTracker(geom), theBuilder(builder), thePropagatorAlong(propagatorAlong), thePropagatorOpposite(propagatorOpposite), theSetMomentum(momFromPSet), theCharges(charges), theErrorRescaling(errorRescaling){}


std::vector<TrajectorySeed*> SeedGeneratorForPixelTelescope::seed(const SeedingHitSet& hits,
						    const PropagationDirection& dir,
							const NavigationDirection&  seedDir,
							const edm::EventSetup& iSetup){
    std::vector<TrajectorySeed*> seeds;
    if (hits.size() == 2){
	if(!theSetMomentum){
	    TrajectorySeed* seed = seedFromPair(hits, dir, seedDir);
			if (seed) seeds.push_back(seed);
		} else {
			for (std::vector<int>::const_iterator iCh = theCharges.begin(); iCh != theCharges.end(); iCh++){
		TrajectorySeed* seed = seedFromPair(hits, dir, seedDir, *iCh);
				if (seed) seeds.push_back(seed);
			}
		}
    } else {
	throw cms::Exception("CombinatorialSeedGeneratorForCosmics") << " Wrong number of hits in Set: "
							<< hits.size() << ", should be 2  ";    
    }
    return seeds;
}



TrajectorySeed* SeedGeneratorForPixelTelescope::seedFromPair(const SeedingHitSet& hits,
							   const PropagationDirection& dir,
			       const NavigationDirection&  seedDir, int charge) const{
    if (hits.size() != 2) {
		throw cms::Exception("CombinatorialSeedGeneratorForCosmics") <<
			"call to SeedGeneratorForPixelTelescope::seedFromPair with " << hits.size() << " hits ";
	}
    auto innerHit = hits[0];
	auto outerHit = hits[1];
	GlobalPoint inner  = theTracker->idToDet(innerHit->geographicalId() )->surface().toGlobal(innerHit->localPosition() );
	GlobalPoint outer  = theTracker->idToDet(outerHit->geographicalId() )->surface().toGlobal(outerHit->localPosition() );
    LogDebug("SeedGeneratorForPixelTelescope") <<
			"Using the following hits: outer(r, phi, theta) (" << outer.perp()
			<< ", " << outer.phi()
			<< "," << outer.theta()
			<<")	inner ("
			<< inner.perp()
			<< ", " << inner.phi()
			<< "," << inner.theta() <<")";
    GlobalPoint* firstPoint  = nullptr;
	GlobalPoint* secondPoint = nullptr;
	int momentumSign	 = 1;
	//const TrackingRecHit* firstHit  = 0;
	//const TrackingRecHit* secondHit = 0;
    std::vector<const BaseTrackerRecHit*> trHits;
	//choose the prop dir and hit order accordingly to where the seed is made
	if (seedDir == outsideIn){
	LogDebug("SeedGeneratorForPixelTelescope")
			<< "Seed from outsideIn alongMomentum OR insideOut oppositeToMomentum";
		firstPoint = &outer;
		secondPoint = &inner;
		//firstHit  = outerHit;
	//secondHit = innerHit;
	trHits.push_back(outerHit);
	trHits.push_back(innerHit);
	} else {
	LogDebug("SeedGeneratorForPixelTelescope")
			<< "Seed from outsideIn oppositeToMomentum OR insideOut alongMomentum";
		firstPoint = &inner;
		secondPoint = &outer;
		momentumSign = -1;
		//firstHit  = innerHit;
	//secondHit = outerHit;
	trHits.push_back(innerHit);
	trHits.push_back(outerHit);
	}
    if (dir == oppositeToMomentum) momentumSign = -1;
    GlobalVector momentum = momentumSign*theP*(*secondPoint-*firstPoint).unit();
	/*GlobalTrajectoryParameters gtp(*firstPoint, 
				       momentum,
				       charge,
				       &(*theMagfield));*/
	//FreeTrajectoryState* startingState = new FreeTrajectoryState(gtp,initialError(trHits[0]));//trHits[1]));
    //if (!qualityFilter(startingState, hits)) return 0;
	//TrajectorySeed* seed = buildSeed(startingState, firstHit, dir);
	//TrajectorySeed* seed = buildSeed(startingState, trHits, dir);
	TrajectorySeed* seed = buildSeed(momentum, charge, trHits, dir);
	return seed;
}


TrajectorySeed* SeedGeneratorForPixelTelescope::buildSeed(const GlobalVector& momentum,
			    int charge, 
			    //const TrackingRecHit* firsthit,
			    std::vector<const BaseTrackerRecHit*>& trHits,
			    const PropagationDirection& dir) const {
    const Propagator* propagator = thePropagatorAlong;
    if (dir == oppositeToMomentum) propagator = thePropagatorOpposite;

    //debug
    GlobalPoint first = theTracker->idToDet(trHits[0]->geographicalId() )->surface().toGlobal(trHits[0]->localPosition() );
	GlobalPoint second  = theTracker->idToDet(trHits[1]->geographicalId() )->surface().toGlobal(trHits[1]->localPosition() );
	LogDebug("SeedGeneratorForPixelTelescope") <<
			"Using the following hits: first(r, phi, theta) (" << first.perp()
			<< ", " << first.phi()
			<< "," << first.theta()
			<<")	second ("
			<< second.perp()
			<< ", " << second.phi()
			<< "," << second.theta() <<")";
    //build an initial trajectory state with big errors, 
    //then update it with the first hit
    auto transHit = trHits[0];
    LocalVector lmom = transHit->surface()->toLocal(momentum);
    LocalTrajectoryParameters ltp(charge/momentum.mag(),
			      lmom.x()/lmom.z(),lmom.y()/lmom.z(),
			      transHit->localPosition().x(),
			      transHit->localPosition().y(),
			      lmom.z()>0.?1.:-1.); 

    LocalTrajectoryError lte(100.,100.,
			(1.+(lmom.x()/lmom.z())*(lmom.x()/lmom.z())),
			(1.+(lmom.y()/lmom.z())*(lmom.y()/lmom.z())),
			1./momentum.mag());
    
    TrajectoryStateOnSurface initialState(ltp,lte,*transHit->surface(),&(*theMagfield));
    KFUpdator updator;
    TrajectoryStateOnSurface startingState = updator.update(initialState,*transHit);

    auto transHit2 = trHits[1];

    //TrajectoryStateOnSurface seedTSOS(*startingState, *plane);
    const TrajectoryStateOnSurface propTSOS = propagator->propagate(startingState,*(transHit2->surface()));
    if (!propTSOS.isValid()){
	LogDebug("SeedGeneratorForPixelTelescope") << "first propagation failed";
	return nullptr;
    }
    TrajectoryStateOnSurface seedTSOS = updator.update(propTSOS, *transHit2);
    if (!seedTSOS.isValid()){
	LogDebug("SeedGeneratorForPixelTelescope") << "first update failed";
		return nullptr;
    }	
    LogDebug("SeedGeneratorForPixelTelescope") << "starting TSOS " << seedTSOS ;
    seedTSOS.rescaleError(theErrorRescaling);
    PTrajectoryStateOnDet PTraj=
	trajectoryStateTransform::persistentState(seedTSOS, trHits[1]->geographicalId().rawId());
    edm::OwnVector<TrackingRecHit> seed_hits;
    //build the transientTrackingRecHit for the starting hit, then call the method clone to rematch if needed
    for (auto ihits = trHits.begin(); ihits != trHits.end(); ihits++){
	seed_hits.push_back((*ihits)->clone());
    }
    TrajectorySeed* seed = new TrajectorySeed(PTraj,seed_hits,dir);
    return seed;
}

bool SeedGeneratorForPixelTelescope::qualityFilter(const SeedingHitSet& hits) const{
    //if we are setting the momentum with the PSet we look for aligned hits
	if (theSetMomentum){
		if (hits.size()==3){
			std::vector<GlobalPoint> gPoints;
			unsigned int nHits = hits.size();
			for (unsigned int iHit = 0; iHit < nHits; ++iHit) gPoints.push_back(hits[iHit]->globalPosition());
			unsigned int subid=(*hits[0]).geographicalId().subdetId();
	    if(subid == StripSubdetector::TEC || subid == StripSubdetector::TID){
				LogDebug("SeedGeneratorForPixelTelescope") 
		    << "In the endcaps we cannot decide if hits are aligned with only phi and z";
		return true;
			}
			FastCircle circle(gPoints[0],
					  gPoints[1],
					  gPoints[2]);
			if (circle.rho() < 500 && circle.rho() != 0) {
				LogDebug("SeedGeneratorForPixelTelescope") <<
					"Seed qualityFilter rejected because rho = " << circle.rho();
				return false;
			}

		}
	}
	return true;
}

bool SeedGeneratorForPixelTelescope::qualityFilter(const GlobalVector& momentum) const{
    if (!theSetMomentum){
	if (momentum.perp() < theP) return false;
    }
    return true; 
}
