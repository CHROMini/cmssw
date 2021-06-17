import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
#process.load('Configuration.Geometry.GeometryExtended2019_cff')
process.load('Geometry.PixelTelescope.PixelTelescopeRecoGeometry_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('FWCore.MessageLogger.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource")

process.MessageLogger = cms.Service("MessageLogger")

process.test = cms.EDAnalyzer("GlobalTrackingGeometryTest")

process.p = cms.Path(process.test)

