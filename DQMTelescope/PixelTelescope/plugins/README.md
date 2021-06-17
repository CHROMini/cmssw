
List of files important to home-made tracking in Geometry/PixelTelescope/plugins:

- SimpleTracking.cc   : main code to do the home-made tracking, 
     after the definition of the SimplePlane (one per module -> SimplePlaneR and SimplePlaneL 
     for the module ont the right side and on the left side respectively), 
     doSeeding searches for the pair of clusters which can be taken as seeds for the tracks,
     doPatternReco aims to add clusters from other layers to build the tracks,
     some distances from track extrapolation to clusters are plotted in DQM_TrackPull*_X,
     ComputeResiduals allows to properly compute the residuals without using the considered
     cluster in the track fit. 
- SimplePlane.h       : definition of the SimplePlane class
- SimplePlane.cc   
- TelescopeTracks.h   : definition of the class for the TelescopeTracks class
- TelescopeTracks.cc  : where is the fitTrack function which minimizes the SumDistance2 function


To run the code, go in Geometry/PixelTelescope/test:

cmsRun PixelTelescope_BeamData_RECO_tracking.py


Change the input directory where to find the raw file on which to run:
opt.register('inputDir',  '/opt/sbg/data/safe1/cms/ccollard/TrackerTelescope/ui6/PixelGeomV2/CMSSW_10_1_11/src/Geometry/PixelTelescope/test/RAW',
             opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
             'Directory of input raw files')

Some parameters can be configured in the python, like 
xmove_param if you want to apply manually some mis-alignment to one module:
process.DQMData = cms.EDAnalyzer('SimpleTracking',
        tracks = cms.untracked.InputTag('ctfWithMaterialTracks'),
        PixelDigisLabel = cms.InputTag("siPixelDigis"),
        PixelClustersLabel = cms.InputTag("siPixelClusters"),
        PixelHitsLabel = cms.InputTag("siPixelRecHits"),
        xmove_param = cms.double(0.),
)

Note that in order to run without trouble, change in SimpleTracking.cc the link to the 
noisy_list.txt present in this Geometry/PixelTelescope/test

