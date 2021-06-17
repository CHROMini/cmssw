import FWCore.ParameterSet.Config as cms

from Geometry.PixelTelescope.PixelTelescopeGeometryXML_cfi import *

trackerNumberingGeometry = cms.ESProducer('TrackerGeometricDetESModule', 
  fromDDD = cms.bool(True),
  appendToDataLabel = cms.string('')
)

from Geometry.CommonDetUnit.globalTrackingGeometry_cfi import *

trackerGeometry = cms.ESProducer('TrackerDigiGeometryESModule',
  appendToDataLabel = cms.string(''),
  fromDDD = cms.bool(True),
  applyAlignment = cms.bool(False),
  alignmentsLabel = cms.string('')
)

trackerTopology = cms.ESProducer('TrackerTopologyEP',
  appendToDataLabel = cms.string('')
)

trackerRecoGeometry = cms.ESProducer("TrackerRecoGeometryESProducer",
    appendToDataLabel = cms.string(''),
    trackerGeometryLabel = cms.untracked.string('')
)
#hltESPTrackerRecoGeometryESProducer
#from RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi import *

