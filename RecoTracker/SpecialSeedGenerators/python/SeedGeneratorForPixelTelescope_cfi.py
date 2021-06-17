# The following comments couldn't be translated into the new config version:

#"StandardHitPairGenerator"
import FWCore.ParameterSet.Config as cms

from RecoTracker.TkTrackingRegions.GlobalTrackingRegion_cfi import *
from RecoLocalTracker.SiStripClusterizer.SiStripClusterChargeCut_cfi import *

layerInfo = cms.PSet(
        FPix = cms.PSet(
	 TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelPairs'),
	 HitProducer = cms.string('siPixelRecHits'),
     )
)

layerList = cms.vstring(
 'FPix4_neg+FPix3_pos', 
)



pixelTelescopeSpecialSeedGenerator = cms.EDProducer("PixelTelescopeSpecialSeedGenerator",
    SeedMomentum = cms.double(0.0), ##initial momentum in GeV !!!set to a lower value for slice test data

    ErrorRescaling = cms.double(50.0),
    RegionFactoryPSet = cms.PSet(
        RegionPSetBlock,
        ComponentName = cms.string('GlobalRegionProducer')
    ),
     Charges = cms.vint32(-1, 1),
    OrderedHitsFactoryPSets = cms.VPSet(
	cms.PSet(
	    ComponentName = cms.string('PixelTelescopeSpecialSeedGenerator'),
	    maxTheta = cms.double(0.1),
	    PropagationDirection = cms.string('alongMomentum'),
	    NavigationDirection = cms.string('outsideIn'),
	    LayerSrc = cms.InputTag("PixelTelescopeLayerPairs") # replace with something with pairs of pixel layers
	), 
	cms.PSet(
	    ComponentName = cms.string('PixelTelescopeSpecialSeedGenerator'),
	    maxTheta = cms.double(0.1),
	    PropagationDirection = cms.string('oppositeToMomentum'),
	    NavigationDirection = cms.string('outsideIn'),
	    LayerSrc = cms.InputTag("PixelTelescopeLayerPairs")
	)),
    UseScintillatorsConstraint = cms.bool(False),
    TTRHBuilder = cms.string('WithTrackAngle'),
    SeedsFromPositiveY = cms.bool(False),
    SeedsFromNegativeY = cms.bool(False),
    doClusterCheck = cms.bool(True),
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(10000),
    MaxNumberOfPixelClusters = cms.uint32(10000),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    CheckHitsAreOnDifferentLayers = cms.bool(False),
    SetMomentum = cms.bool(True),
    requireBOFF = cms.bool(False),
    maxSeeds = cms.int32(10),

    
    
)


