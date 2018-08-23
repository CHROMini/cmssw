import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('PixelTelescope'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
