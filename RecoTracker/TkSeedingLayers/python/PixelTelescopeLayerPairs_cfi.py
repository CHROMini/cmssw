import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *

pixelTelescopeLayerPairs = seedingLayersEDProducer.clone()

pixelTelescopeLayerPairs.layerList = cms.vstring('FPix4+BPix3'
)

pixelTelescopeLayerPairs.BPix = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    HitProducer = cms.string('siPixelRecHits'),
)

pixelTelescopeLayerPairs.FPix = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    HitProducer = cms.string('siPixelRecHits'),
)



