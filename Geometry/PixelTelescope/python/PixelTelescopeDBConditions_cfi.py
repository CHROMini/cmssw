import FWCore.ParameterSet.Config as cms
import CalibTracker.Configuration.Common.PoolDBESSource_cfi
# ----------------------- CONDITIONS -------------------

# --------------------- SiPixelQuality -----------------
SiPixelQualityDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelQualityFromDbRcd'),
      tag = cms.string('SiPixelQuality_phase1_2017_v5'),  
    )
  )
)
es_prefer_Quality = cms.ESPrefer("PoolDBESSource","SiPixelQualityDBReader")
#-------------------------------------------------------

# --------------------- CablingMap ---------------------
CablingMapDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  connect = cms.string('sqlite_file:./pixeltelescope_cabling.db'), #local DB
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelFedCablingMapRcd'),
      tag = cms.string('SiPixelFedCablingMap_pixeltelescope_v1'),    
    )
  )
)
es_prefer_CablingReader = cms.ESPrefer("PoolDBESSource","CablingMapDBReader")
#-------------------------------------------------------

# --------------------- Gain Calib DB ------------------
GainDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),#from Prep
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelGainCalibrationOfflineRcd'),
      tag = cms.string('SiPixelGainCalibration_phase1_mc_v2')
    )
  )
)
es_prefer_Gain = cms.ESPrefer("PoolDBESSource","GainDBReader")
#-------------------------------------------------------

# --------------------- Gain Calib HLT DB ------------------
GainDBReaderHLT = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),#from Prep
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelGainCalibrationForHLTRcd'),
      tag = cms.string('SiPixelGainCalibration_hlt_phase1_mc_v2')
    )
  )
)
es_prefer_GainHLT = cms.ESPrefer("PoolDBESSource","GainDBReaderHLT")
#-------------------------------------------------------

# ----------------------- GenError ---------------------
GenErrReader = cms.ESSource("PoolDBESSource",
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),#from Prep
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelGenErrorDBObjectRcd'),
      tag = cms.string('SiPixelGenErrorDBObject_phase1_38T_mc_v2'), 
    )
  )
)
es_prefer_GenErr = cms.ESPrefer("PoolDBESSource","GenErrReader")
#-------------------------------------------------------

# --------------------- LorentzAngle -------------------

LorentzAngleDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  connect = cms.string ('sqlite_file:./SiPixelLorentzAngle.db'), #local DB
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelLorentzAngleRcd'),
      label = cms.untracked.string(''),
      tag = cms.string('SiPixelLorentzAngle_pixeltelescope_v1'),
    )
  )
)
es_prefer_LA = cms.ESPrefer("PoolDBESSource","LorentzAngleDBReader")
#-------------------------------------------------------

LorentzAngleSimDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  connect = cms.string ('sqlite_file:./SiPixelLorentzAngleSim.db'), #local DB
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelLorentzAngleSimRcd'),
      label = cms.untracked.string(''),
      tag = cms.string('SiPixelLorentzAngle_pixeltelescope_mc_v1'),     )
  )
)
es_prefer_LASim = cms.ESPrefer("PoolDBESSource","LorentzAngleSimDBReader")
#-------------------------------------------------------
