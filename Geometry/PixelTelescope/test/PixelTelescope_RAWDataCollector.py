import FWCore.ParameterSet.Config as cms

process = cms.Process("RAW")

import FWCore.ParameterSet.VarParsing as opts
opt = opts.VarParsing ('analysis')

opt.register('inputDir',  './',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'Directory of input raw files')

opt.register('outputFileName', 'PixelTelescope_BeamData_RAW.root',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'Name of merged output raw file')

opt.parseArguments()

import os
##my_path = "/data/veszpv/project/TelescopeData/beam/run001055/"
my_path = opt.inputDir
my_extensions = ['raw']
file_names = ["file:"+os.path.join(my_path, fn) for fn in os.listdir(my_path)
              if any(fn.endswith(ext) for ext in my_extensions)]


process.source = cms.Source("FRDStreamSource",
                            verbosity = cms.untracked.uint32(2),
                            fileNames = cms.untracked.vstring(file_names),
                            useL1EventID = cms.untracked.bool(True)
)

#print process.source.fileNames

process.rawDataCollector = cms.EDProducer("RemapRawDataCollection", # "RawDataCollectorByLabel",
    verbose = cms.untracked.int32(0),     # 0 = quiet, 1 = collection list, 2 = FED list
    RawCollectionList = cms.VInputTag(cms.InputTag('source')),
    FEDMapping = cms.untracked.VPSet(
        cms.PSet(
            src = cms.int32(1),
            dest = cms.int32(1200),
        ),
    )
)

process.collect_raw = cms.Sequence(process.rawDataCollector)
process.produce_raw = cms.Path(process.collect_raw)

process.RAWOutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string(opt.outputFileName),
    outputCommands = cms.untracked.vstring(
     'drop *',
     'keep *_rawDataCollector_*_*'
    ),
    splitLevel = cms.untracked.int32(0)
)

process.output_step = cms.EndPath(process.RAWOutput)
process.schedule = cms.Schedule(process.produce_raw,
                                process.output_step)
 
