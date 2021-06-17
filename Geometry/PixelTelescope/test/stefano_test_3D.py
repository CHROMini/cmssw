import FWCore.ParameterSet.Config as cms


#from RecoTracker.SpecialSeedGenerators.SeedGeneratorForPixelTelescope_cfi import *
from RecoTracker.TkSeedingLayers.PixelTelescopeLayerPairs_cfi import *
from RecoTracker.TkTrackingRegions.GlobalTrackingRegion_cfi import *
from RecoLocalTracker.SiStripClusterizer.SiStripClusterChargeCut_cfi import *

from Configuration.StandardSequences.Eras import eras
process = cms.Process('RECO',eras.Run2_25ns)

import FWCore.ParameterSet.VarParsing as opts
opt = opts.VarParsing ('analysis')

opt.register('inputDir',  '/afs/cern.ch/work/m/mersi/public/Chromie/RAW/run100208',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'Directory of input raw files')

opt.register('outputFileName', 'PixelTelescope_BeamData_RAW.root',
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

############## Tracking reco
##### For the seeding
process.load('RecoTracker.SpecialSeedGenerators.SeedGeneratorForPixelTelescope_cfi')
process.load('RecoTracker.TkSeedingLayers.PixelTelescopeLayerPairs_cfi')


process.pixelTelescopeLayerPairs = seedingLayersEDProducer.clone()

process.pixelTelescopeLayerPairs.layerList = cms.vstring('FPix4+BPix3'
)

process.pixelTelescopeLayerPairs.BPix = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    HitProducer = cms.string('siPixelRecHits'),
)

process.pixelTelescopeLayerPairs.FPix = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    HitProducer = cms.string('siPixelRecHits'),
)



layerInfo = cms.PSet(
        FPix = cms.PSet(
	 TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelPairs'),
	 HitProducer = cms.string('siPixelRecHits'),
     )
)

layerList = cms.vstring(
 'FPix4_neg+FPix3_pos', 
)



process.pixelTelescopeSpecialSeedGenerator = cms.EDProducer("PixelTelescopeSpecialSeedGenerator",
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






process.ctftracksPixelTelescope = cms.Sequence( process.pixelTelescopeLayerPairs*process.pixelTelescopeSpecialSeedGenerator)


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)



process.reconstruction = cms.Path(process.reconstruction_pixelOnly*process.ctftracksPixelTelescope)


process.schedule = cms.Schedule(process.reconstruction,process.DQM,process.endjob_step,process.RECOSIMoutput_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

