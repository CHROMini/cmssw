import FWCore.ParameterSet.Config as cms


from RecoTracker.TkTrackingRegions.GlobalTrackingRegion_cfi import *
from RecoLocalTracker.SiStripClusterizer.SiStripClusterChargeCut_cfi import *

from Configuration.StandardSequences.Eras import eras
process = cms.Process('RECO',eras.Run2_25ns)

import FWCore.ParameterSet.VarParsing as opts
opt = opts.VarParsing ('analysis')

opt.register('inputDir',  '/afs/cern.ch/work/m/mersi/public/Chromie/RAW/run100208',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'Directory of input raw files')

opt.register('outputFileName', 'PixelTelescope_BeamData_OUTPUT.root',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'Name of output reco file')

opt.register('outputDQMFileName', 'PixelTelescope_BeamData_DQM.root',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'Name of output dqm file')

opt.parseArguments()

 
#process.MessageLogger = cms.Service("MessageLogger",
#        destinations = cms.untracked.vstring(                           #1
#                'myOutputFile'                                          #2
#        ),
#        myOutputFile = cms.untracked.PSet(                              #3
#                threshold = cms.untracked.string( 'INFO' )          #4
#        ),
#)        


import os
##my_path = "/data/veszpv/project/TelescopeData/beam/run001055/"
my_path = opt.inputDir
my_extensions = ['root']
file_names = ["file:"+os.path.join(my_path, fn) for fn in os.listdir(my_path)
              if any(fn.endswith(ext) for ext in my_extensions)]

process.load('Geometry.PixelTelescope.PixelTelescopeRecoGeometry_cfi')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi')
process.siPixelDigis.UsePhase1 = cms.bool(True)
process.siPixelDigis.InputLabel = cms.InputTag("rawDataCollector")

process.siPixelClusters = cms.EDProducer("SiPixelClusterProducer",
    SiPixelGainCalibrationServiceParameters = cms.PSet(),
    src = cms.InputTag("siPixelDigis"),
    ChannelThreshold = cms.int32(10),
    MissCalibrate = cms.untracked.bool(False),#True
    SplitClusters = cms.bool(False),
    VCaltoElectronGain    = cms.int32(47),
    VCaltoElectronGain_L1 = cms.int32(50),
    VCaltoElectronOffset    = cms.int32(-60),  
    VCaltoElectronOffset_L1 = cms.int32(-670),  
    # **************************************
    # ****  payLoadType Options         ****
    # ****  HLT - column granularity    ****
    # ****  Offline - gain:col/ped:pix  ****
    # **************************************
    payloadType = cms.string('Offline'),
    #payloadType = cms.string('Full'),
    SeedThreshold = cms.int32(1000),
    ClusterThreshold    = cms.int32(4000),
    ClusterThreshold_L1 = cms.int32(2000),
    # **************************************
    maxNumberOfClusters = cms.int32(-1), # -1 means no limit.
)

process.bysipixelclustmulteventfilter = cms.EDFilter('BySiPixelClusterMultiplicityEventFilter',
                                                     multiplicityConfig = cms.PSet(
                                                         collectionName = cms.InputTag("siPixelClusters"),
                                                         moduleThreshold = cms.untracked.int32(-1),
                                                         useQuality = cms.untracked.bool(False),
                                                         qualityLabel = cms.untracked.string("")
                                                         ),
                                                     cut = cms.string("mult > 0")
                                                     )

process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
process.load('RecoLocalTracker.SiStripRecHitConverter.StripCPEESProducer_cfi') # if needed...
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = cms.bool(False)
process.PixelCPEGenericESProducer.TruncatePixelCharge = cms.bool(False)
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = cms.bool(False)
process.PixelCPEGenericESProducer.DoCosmics = cms.bool(False)
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = cms.bool(False)
process.PixelCPEGenericESProducer.useLAWidthFromDB = cms.bool(False)
process.PixelCPEGenericESProducer.lAOffset = cms.double(0.1)
process.PixelCPEGenericESProducer.lAWidthFPix = cms.double(0.1)

process.siPixelRecHits = cms.EDProducer("SiPixelRecHitConverter",
    src = cms.InputTag("siPixelClusters"),
    CPE = cms.string('PixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
)

process.pixeltrackerlocalreco = cms.Sequence(process.siPixelDigis*process.siPixelClusters*process.siPixelRecHits)

process.reconstruction_pixelOnly = cms.Sequence(
    process.pixeltrackerlocalreco
)

process.reconstruction = cms.Path(process.reconstruction_pixelOnly)


process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Geometry.PixelTelescope.PixelTelescopeDBConditions_cfi')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(file_names),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('-s nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string(opt.outputFileName),
    outputCommands = cms.untracked.vstring(
     'keep *',
     'drop *_simCastorDigis_*_*',
     'drop *_simEcalUnsuppressedDigis_*_*',
     'drop *_simHcalUnsuppressedDigis_*_*',
     'drop *_simSiPixelDigis_*_*',
     'drop *_simSiStripDigis_*_*',
    ),                                         
    splitLevel = cms.untracked.int32(0)
)

## DQM

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string(opt.outputDQMFileName),
      closeFileFast = cms.untracked.bool(True)
)

#process.DQMData = cms.EDAnalyzer('PixelTelescope', 
process.DQMData = cms.EDAnalyzer('Ana3D', 
 	tracks = cms.untracked.InputTag('ctfWithMaterialTracks'),
	PixelDigisLabel = cms.InputTag("siPixelDigis"),
	PixelClustersLabel = cms.InputTag("siPixelClusters"),
	PixelHitsLabel = cms.InputTag("siPixelRecHits"),
)

process.DQM = cms.Path(process.DQMData)



############## Track reconstruction ##############
############## Track reconstruction ##############
############## Track reconstruction ##############
############## Track reconstruction ##############
############## Track reconstruction ##############


# Tracking configuration file fragment for P5 cosmic running
#


#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

#process.load('Geometry.PixelTelescope.PixelTelescopeDBConditions_cfi')
#process.load('Geometry.PixelTelescope.PixelTelescopeRecoGeometry_cfi')

process.hltESPTrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer",
    appendToDataLabel = cms.string(''),
    trackerGeometryLabel = cms.untracked.string('')
)


from RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cff import *
# TTRHBuilders
from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *
process.load('RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff')

#from RecoTracker.SpecialSeedGenerators.CombinatorialSeedGeneratorForCosmicsP5_cff import *
process.load('RecoTracker.SpecialSeedGenerators.CombinatorialSeedGeneratorForCosmicsP5_cff')

from RecoTracker.SpecialSeedGenerators.SimpleCosmicBONSeeder_cff import *
from RecoTracker.TkSeedGenerator.GlobalCombinedSeeds_cff import *
combinedP5SeedsForCTF = RecoTracker.TkSeedGenerator.GlobalCombinedSeeds_cfi.globalCombinedSeeds.clone()
combinedP5SeedsForCTF.seedCollections = cms.VInputTag(
    cms.InputTag('combinatorialcosmicseedfinderP5'),
    cms.InputTag('simpleCosmicBONSeeds'),
)
#backward compatibility 2.2/3.1
combinedP5SeedsForCTF.PairCollection = cms.InputTag('combinatorialcosmicseedfinderP5')
combinedP5SeedsForCTF.TripletCollection = cms.InputTag('simpleCosmicBONSeeds')

from RecoTracker.CkfPattern.CkfTrackCandidatesP5_cff import *
# HEREEEEEEEEEE
process.ckfTrackCandidatesP5 = ckfTrackCandidatesP5.clone()
# HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	    

ckfTrackCandidatesP5.src = cms.InputTag('combinedP5SeedsForCTF')
#backward compatibility 2.2/3.1
#ckfTrackCandidatesP5.SeedProducer = 'combinedP5SeedsForCTF'
#import RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cfi

from RecoTracker.TrackProducer.CTFFinalFitWithMaterialP5_cff import *
# Final Track Selector for CTF
from RecoTracker.FinalTrackSelectors.CTFFinalTrackSelectorP5_cff import *

process.load("TrackingTools.TrackRefitter.ctfWithMaterialTrajectoriesP5_cff")
process.load("RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderP5_cff")
process.load("TrackingTools.TrajectoryCleaning.TrajectoryCleanerBySharedHits_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi")


#process.ctftracksP5 = cms.Sequence(combinatorialcosmicseedinglayersP5+combinatorialcosmicseedfinderP5*
#			    ckfTrackCandidatesP5*ctfWithMaterialTracksCosmics*ctfWithMaterialTracksP5)
			    
# HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	    
process.ctftracksP5 = cms.Sequence(process.combinatorialcosmicseedingtripletsP5+process.combinatorialcosmicseedfinderP5*
			    process.ckfTrackCandidatesP5)	
# HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	    

#process.ctftracksP5 = cms.Sequence(process.combinatorialcosmicseedinglayersP5)





from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)



process.reconstruction = cms.Path(process.reconstruction_pixelOnly*process.ctftracksP5)


process.schedule = cms.Schedule(process.reconstruction,process.DQM,process.endjob_step,process.RECOSIMoutput_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

