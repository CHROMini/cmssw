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

opt.register('outputFileName', 'PixelTelescope_BeamData_OUTPUT_RECO.root',
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
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')    ##### Why would that global tag be ok??????

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
	PixelDigisLabel = cms.InputTag("siPixelDigis"),
	PixelClustersLabel = cms.InputTag("siPixelClusters"),
	PixelHitsLabel = cms.InputTag("siPixelRecHits"),
	SeedFinderLabel = cms.InputTag("combinatorialcosmicseedfinderP5"),
	TrackCandidateLabel = cms.InputTag("ckfTrackCandidatesP5"),
	tracks = cms.untracked.InputTag('ctfWithMaterialTracksCosmics'),
)

process.DQM = cms.Path(process.DQMData)










###########################################################################################
############################### Track reconstruction ######################################
###########################################################################################
####### Tracking configuration file fragment for P5 cosmic running ########################



############################# MEASUREMENT AND PROPAGATION ###################################
from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *
process.load('RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff')
process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi")


from RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cff import *
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi")

process.MeasurementTrackerEvent.stripClusterProducer = ''



import TrackingTools.MaterialEffects.MaterialPropagator_cfi
MaterialPropagatorPtMin035 = TrackingTools.MaterialEffects.MaterialPropagator_cfi.MaterialPropagator.clone(
    ComponentName = 'PropagatorWithMaterialPtMin035',
    ptMin = 0.35
)

import RecoTracker.TkSeedGenerator.SeedFromConsecutiveHitsStraightLineCreator_cfi
process.SeedCreatorPSet = RecoTracker.TkSeedGenerator.SeedFromConsecutiveHitsStraightLineCreator_cfi.SeedFromConsecutiveHitsStraightLineCreator.clone(
    propagator = cms.string('PropagatorWithMaterialPtMin035')
)



############################# SEEDING ###################################


# obsolete stuff
#from RecoTracker.SpecialSeedGenerators.CombinatorialSeedGeneratorForCosmicsP5_cff import *
#process.load('RecoTracker.SpecialSeedGenerators.CombinatorialSeedGeneratorForCosmicsP5_cff')
#from RecoTracker.SpecialSeedGenerators.SimpleCosmicBONSeeder_cff import *
#from RecoTracker.TkSeedGenerator.GlobalCombinedSeeds_cff import *
#combinedP5SeedsForCTF = RecoTracker.TkSeedGenerator.GlobalCombinedSeeds_cfi.globalCombinedSeeds.clone()
#combinedP5SeedsForCTF.seedCollections = cms.VInputTag(
#    cms.InputTag('combinatorialcosmicseedfinderP5'),
#    cms.InputTag('simpleCosmicBONSeeds'),
#)
#backward compatibility 2.2/3.1
#combinedP5SeedsForCTF.PairCollection = cms.InputTag('combinatorialcosmicseedfinderP5')
#combinedP5SeedsForCTF.TripletCollection = cms.InputTag('simpleCosmicBONSeeds')


process.load("RecoTracker.SpecialSeedGenerators.CombinatorialSeedGeneratorForCosmicsP5_cff")
process.combinatorialcosmicseedingtripletsP5 = cms.EDProducer("SeedingLayersEDProducer",
    MTIB = cms.PSet(
    ),
    MTOB = cms.PSet(
    ),
    TEC = cms.PSet(
    ),
    TIB = cms.PSet(
    ),
    TOB = cms.PSet(
    ),
    FPix = cms.PSet( 
      hitErrorRPhi = cms.double( 0.0051 ),
      #TTRHBuilder = cms.string( "WithoutRefit" ),
      TTRHBuilder = cms.string( "WithTrackAngle" ),
      #TTRHBuilder = cms.string( "myTTRHBuilderWithoutAngle4MixedTriplets" ),
      useErrorsFromParam = cms.bool( True ),
      hitErrorRZ = cms.double( 0.0036 ),
      HitProducer = cms.string( "siPixelRecHits" )
    ),
    BPix = cms.PSet(
    ),

    layerList = cms.vstring(
        'FPix3_neg+FPix2_neg+FPix1_neg', 
        'FPix4_neg+FPix2_neg+FPix1_neg', 
        'FPix4_neg+FPix3_neg+FPix1_neg', 
        'FPix4_neg+FPix3_neg+FPix2_neg', 
    )
)
process.combinatorialcosmicseedingpairsTECnegP5.layerList = cms.vstring()
process.combinatorialcosmicseedingpairsTECposP5.layerList = cms.vstring()
process.combinatorialcosmicseedingpairsTOBP5.layerList = cms.vstring()







# TUNEEEEEEEEEEE
process.combinatorialcosmicseedfinderP5 = cms.EDProducer("CtfSpecialSeedGenerator",
    Charges = cms.vint32(-1),
    CheckHitsAreOnDifferentLayers = cms.bool(False),
    ClusterCollectionLabel = cms.InputTag(""),
    DontCountDetsAboveNClusters = cms.uint32(20),
    ErrorRescaling = cms.double(50.0),

    MaxNumberOfCosmicClusters = cms.uint32(300),
    MaxNumberOfPixelClusters = cms.uint32(300000),
    OrderedHitsFactoryPSets = cms.VPSet(cms.PSet(
        ComponentName = cms.string('GenericTripletGenerator'),
        LayerSrc = cms.InputTag("combinatorialcosmicseedingtripletsP5"),
        maxTheta = cms.double(0.1),                             #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        NavigationDirection = cms.string('outsideIn'),       #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        PropagationDirection = cms.string('alongMomentum')    #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    )),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            originHalfLength = cms.double(50.0),   #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            originRadius = cms.double(50.0),
            originXPos = cms.double(50.0),
            originYPos = cms.double(0.0),
            originZPos = cms.double(-100.0),
            precise = cms.bool(False),            #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            ptMin = cms.double(0.9),              #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedMomentum = cms.double(5.0),
    SeedsFromNegativeY = cms.bool(False),          #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    SeedsFromPositiveY = cms.bool(False),          #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    SetMomentum = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    UseScintillatorsConstraint = cms.bool(False),
    LowerScintillatorParameters = cms.PSet(
        GlobalX = cms.double(50.0),
        GlobalY = cms.double(0.0),
        GlobalZ = cms.double(-100.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    UpperScintillatorParameters = cms.PSet(
        GlobalX = cms.double(50.0),
        GlobalY = cms.double(0.0),
        GlobalZ = cms.double(300.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    doClusterCheck = cms.bool(True),
    maxSeeds = cms.int32(10000),
    requireBOFF = cms.bool(False)          #TUNEEEEEEEEEEEEEEEEEEEEEEEEEEEE
)





############################# TRACKS CANDIDATES ###################################


#process.load('RecoTracker.TkNavigation.CosmicsNavigationSchoolESProducer_cfi')
process.load('RecoTracker.TkNavigation.TelescopeNavigationSchoolESProducer_cfi')


from RecoTracker.CkfPattern.CkfTrackCandidatesP5_cff import *
process.ckfTrackCandidatesP5 = ckfTrackCandidatesP5.clone()
process.ckfTrackCandidatesP5.src = cms.InputTag('combinatorialcosmicseedfinderP5')
process.ckfTrackCandidatesP5.NavigationSchool = cms.string("TelescopeNavigationSchool")
#backward compatibility 2.2/3.1
#ckfTrackCandidatesP5.SeedProducer = 'combinedP5SeedsForCTF'
#import RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cfi




############################# TRACKS PRODUCER ###################################

process.load("RecoTracker.TrackProducer.CTFFinalFitWithMaterialP5_cff")

process.ctfWithMaterialTracksCosmics = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('ctf'),
    Fitter = cms.string('FittingSmootherRKP5'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('TelescopeNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("pixel_telescope_beamspot_tag"),  # NB: Beam spot is not valid, NEED TO FIX THIS. See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFindingBeamSpot
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("ckfTrackCandidatesP5"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfWithMaterialTracksP5 = process.ctfWithMaterialTracksCosmics.clone(
    src = "ctfWithMaterialTracksCosmics"
)

#process.ctfWithMaterialTracksP5 = cms.EDProducer("TrackProducer",
#    AlgorithmName = cms.string('ctf'),
#    Fitter = cms.string('FittingSmootherRKP5'),
#    GeometricInnerState = cms.bool(True),
#    MeasurementTracker = cms.string(''),
#    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
#    NavigationSchool = cms.string('TelescopeNavigationSchool'),
#    Propagator = cms.string('RungeKuttaTrackerPropagator'),
#    SimpleMagneticField = cms.string(''),
 #   TTRHBuilder = cms.string('WithAngleAndTemplate'),
 #   TrajectoryInEvent = cms.bool(False),
#    alias = cms.untracked.string('ctfWithMaterialTracks'),
#    beamSpot = cms.InputTag("offlineBeamSpot"),
#    clusterRemovalInfo = cms.InputTag(""),
 #   src = cms.InputTag("ctfWithMaterialTracksCosmics"),
#    useHitsSplitting = cms.bool(False),
 #   useSimpleMF = cms.bool(False)
#)





############################# TRACKS FITTING ###################################


#process.load("TrackingTools.TrackRefitter.ctfWithMaterialTrajectoriesP5_cff")

process.RKTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('RKFitter'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)

process.RKTrajectorySmoother.RecoGeometry = cms.string("GlobalDetLayerGeometry")

process.load("RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderP5_cff")

process.ckfBaseTrajectoryFilterP5.minimumNumberOfHits = 0
process.ckfBaseTrajectoryFilterP5.minPt = 0.0001
process.ckfBaseTrajectoryFilterP5.maxLostHits = 9999
process.ckfBaseTrajectoryFilterP5.maxConsecLostHits = 8

process.CkfBaseTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(8),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.001),
    minimumNumberOfHits = cms.int32(0),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.load("TrackingTools.TrajectoryCleaning.TrajectoryCleanerBySharedHits_cfi")



############################# TRACKS SELECTOR ###################################

#import RecoTracker.FinalTrackSelectors.cosmictrackSelector_cfi
#process.cosmictrackSelector = cms.EDProducer("CosmicTrackSelector",
#                                     src = cms.InputTag("ctfWithMaterialTracksCosmics"),
#                                     keepAllTracks = cms.bool(False), ## if set to true tracks failing this filter are kept in the output
#                                     beamspot = cms.InputTag(""),
#                                     #untracked bool copyTrajectories = true // when doing retracking before
#                                     copyTrajectories = cms.untracked.bool(False),
#                                     copyExtras = cms.untracked.bool(True), ## set to false on AOD
#                                     qualityBit = cms.string(''), # set to '' or comment out if you don't want to set the bit                                     
#                                     # parameters for adapted optimal cuts on chi2
#                                     chi2n_par = cms.double(10.0),
#                                     # Impact parameter absolute cuts.
#                                     max_d0 = cms.double(110.),
#                                     max_z0 = cms.double(300.),
#                                     # track parameter cuts 
#                                     max_eta = cms.double(2.0),
#                                     min_pt = cms.double(1.0),
                                     # Cut on numbers of valid hits
#                                     min_nHit = cms.uint32(5),
                                     # Cut on number of Pixel Hit 
#                                     min_nPixelHit = cms.uint32(0),
                                     # Cuts on numbers of layers with hits/3D hits/lost hits. 
#                                     minNumberLayers = cms.uint32(0),
#                                     minNumber3DLayers = cms.uint32(0),
#                                    maxNumberLostLayers = cms.uint32(999)
#                                     )

# Final Track Selector for CTF
#from RecoTracker.FinalTrackSelectors.CTFFinalTrackSelectorP5_cff import *



############################# DEBUG ANALYZER ###################################
### and an analyzer
#process.trajCout = cms.EDAnalyzer('TrajectoryAnalyzer',
#   trajectoryInput=cms.InputTag('ctfWithMaterialTracksCosmics')
#)
#process.trajCout = cms.EDAnalyzer('TrackValidator')
#process.Ana = cms.Path(process.trajCout)




########################################## DEFINE SEQUENCE ##########################################

#process.ctftracksP5 = cms.Sequence(process.combinatorialcosmicseedinglayersP5)

process.ctftracksP5 = cms.Sequence(process.combinatorialcosmicseedingtripletsP5+process.combinatorialcosmicseedfinderP5*process.MeasurementTrackerEvent*
			    process.ckfTrackCandidatesP5*process.ctfWithMaterialTracksCosmics)
			    
#process.ctftracksP5 = cms.Sequence(process.combinatorialcosmicseedingtripletsP5+process.combinatorialcosmicseedfinderP5*process.MeasurementTrackerEvent*
#			    process.ckfTrackCandidatesP5*process.ctfWithMaterialTracksCosmics*process.ctfWithMaterialTracksP5)			    




process.reconstruction = cms.Path(process.reconstruction_pixelOnly*process.ctftracksP5)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)


process.schedule = cms.Schedule(process.reconstruction,process.DQM,process.endjob_step,process.RECOSIMoutput_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
print process.dumpPython()
