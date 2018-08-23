import FWCore.ParameterSet.Config as cms

from Geometry.PixelTelescope.PixelTelescopeGeometryXML_cfi import *

trackerNumberingGeometry = cms.ESProducer('TrackerGeometricDetESModule', 
  fromDDD = cms.bool(True),
  appendToDataLabel = cms.string('')
)

trackerGeometry = cms.ESProducer('TrackerDigiGeometryESModule',
  appendToDataLabel = cms.string(''),
  fromDDD = cms.bool(True),
  applyAlignment = cms.bool(False),
  alignmentsLabel = cms.string('')
)

trackerTopology = cms.ESProducer('TrackerTopologyEP',
  appendToDataLabel = cms.string('')
)
