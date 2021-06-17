import FWCore.ParameterSet.Config as cms

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    rootNodeName = cms.string('telescope:BeamArea'),
    geomXMLFiles = cms.vstring(
        # Default CMS materials
        'Geometry/CMSCommonData/data/materials.xml',

        # Define standalone Phase 1 BPIX module and associated materials
        'Geometry/TrackerCommonData/data/trackermaterial.xml',
        'Geometry/TrackerCommonData/data/PhaseI/pixbarmaterial.xml',
        'Geometry/TrackerCommonData/data/PhaseI/pixfwdMaterials.xml',
        'Geometry/PixelTelescope/data/Phase1BPIXLayer4Module.xml',

        # Define telescope	
        'Geometry/PixelTelescope/data/telescope.xml',

        # Define DUT
        'Geometry/PixelTelescope/data/DUT.xml',
        
        # Tune geometry parameters
        'Geometry/PixelTelescope/data/tuneParameters.xml',
        
        # Sim files
        'Geometry/PixelTelescope/data/Sim/pixelTelescopeTopology.xml',
        'Geometry/PixelTelescope/data/Sim/pixelTelescopeSensitive.xml',
        'Geometry/PixelTelescope/data/Sim/pixelTelescopeProdCuts.xml',
        
        # Reco file
        'Geometry/PixelTelescope/data/Reco/pixelTelescopeRecoMaterial.xml',
        
        # DetId Scheme file
        'Geometry/TrackerCommonData/data/PhaseI/trackerParameters.xml'    
        #'Geometry/TrackerCommonData/data/PhaseII/trackerParameters.xml'  # Load this instead to load Phase 1 AND Phase 2 schemes.
                                                                          # Warning: take care with Phase 1 BIG_PIX values, look at 
                                                                          # Geometry/TrackerGeometryBuilder/src/TrackerGeomBuilderFromGeometricDet.cc


        # Local copies made in case release changes, but loading them from the release seems to be better:
        # 'Geometry/PixelTelescope/data/LocalCopy/materials.xml',
        # 'Geometry/PixelTelescope/data/LocalCopy/trackermaterial.xml',
        # 'Geometry/PixelTelescope/data/LocalCopy/pixbarmaterial.xml',
        # 'Geometry/PixelTelescope/data/LocalCopy/pixfwdMaterials.xml'
    )
)
