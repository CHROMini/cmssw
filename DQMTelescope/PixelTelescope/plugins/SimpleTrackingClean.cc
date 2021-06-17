// Package:    DQMTelescope/SimpleTrackingClean
// Class:      SimpleTrackingClean
//
// class SimpleTrackingClean SimpleTrackingClean.cc DQMTelescope/SimpleTrackingClean/plugins/SimpleTrackingClean.cc

// Original Author:  Jeremy Andrea
//      Created:     03.08.2018
//      Updated by:  Caroline Collard
//      Updated by:  Nikkie Deelen
//      Updated by:  Patrick Asenov

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"

#include "DQMTelescope/PixelTelescope/plugins/TelescopeTracks.h"
#include "DQMTelescope/PixelTelescope/plugins/SimplePlane.h"

#include <cstring>
#include <string>
#include <TH2F.h>
#include <TTree.h>
#include <TString.h>
#include <TVector3.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TVirtualFitter.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <iostream>
#include <fstream>
//
// class declaration
//

using reco::TrackCollection ;
using namespace ROOT::Math;

class SimpleTrackingClean : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:

    explicit SimpleTrackingClean ( const edm::ParameterSet& ) ;
    ~SimpleTrackingClean ( ) ;

    static void fillDescriptions ( edm::ConfigurationDescriptions& descriptions ) ;

private:

    virtual void beginJob ( ) override ;
    virtual void analyze ( const edm::Event&, const edm::EventSetup& ) override ;
    virtual void endJob ( ) override ;

    // ----------member data ---------------------------
    edm::EDGetTokenT< TrackCollection >                        tracksToken_ ;
    edm::EDGetTokenT< edm::DetSetVector< PixelDigi > >         pixeldigiToken_ ;
    edm::EDGetTokenT< edmNew::DetSetVector< SiPixelCluster > > pixelclusterToken_ ;
    edm::EDGetTokenT< edmNew::DetSetVector< SiPixelRecHit > >  pixelhitToken_ ;

    edm::Service<TFileService> fs ;

    // information to stored in the output file
    std::map< uint32_t, TH1F* > DQM_ClusterCharge ;
    std::map< uint32_t, TH1F* > DQM_ClusterSize_X ;
    std::map< uint32_t, TH1F* > DQM_ClusterSize_Y ;
    std::map< uint32_t, TH1F* > DQM_ClusterSize_XY ;
    std::map< uint32_t, TH1F* > DQM_NumbOfClusters_per_Event;
    std::map< uint32_t, TH2F* > DQM_ClusterPosition ;
    std::map< uint32_t, TH1F* > DQM_Chi2 ;
    std::map< uint32_t, TH1F* > DQM2_Chi2 ;
    std::map< uint32_t, TH1F* > DQM2_Sum2 ;

    // Correlation plots for the telescope
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation_X ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation_Y ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation2_X ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation2_Y ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation3_X ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation3_Y ;

    std::map< uint32_t, TH2F*> DQM_Corr_dX_X ;
    std::map< uint32_t, TH2F*> DQM_Corr_dY_Y ;

    // PixelNoisePlot
    std::map< uint32_t, TH1F* > DQM_PixelNoise ;

    // Residuals
    std::map< uint32_t, TH1F* > DQM_TrackPull_X ;
    std::map< uint32_t, TH1F* > DQM_TrackPull_Y ;
    std::map< uint32_t, TH1F* > DQM_TrackPull2_X ;
    std::map< uint32_t, TH1F* > DQM_TrackPull2_Y ;
    std::map< uint32_t, TH1F* > DQM_TrackPull3_X ;
    std::map< uint32_t, TH1F* > DQM_TrackPull3_Y ;
    std::map< uint32_t, TH1F* > DQM2_TrackPull_X ;
    std::map< uint32_t, TH1F* > DQM2_TrackPull_Y ;

    TH1F* DQM_NumbOfSeeds_per_Event;
    TH1F* DQM_NumbOfTracks_per_Event;
    TH1F* DQM_NumbOfClus_per_Event;

    TH1F* DQM_NumbOfSeeds_4cl_per_Event;
    TH1F* DQM_NumbOfTracks_4cl_per_Event;
    TH1F* distanceR_per_Track;
    TH1F* deltaX_per_Track;
    TH1F* deltaY_per_Track;

    //TH2F* testTrack;
    //TH2F* testModule;
    //TH2F* testModule2;
    //TH2F* testModuleLayer;

    TH1F* DQM_TrackPull4_X_M3124;
    TH1F* DQM_TrackPull4_Y_M3124;
    TH1F* DQM_TrackPull5_X_M3124;
    TH1F* DQM_TrackPull5_Y_M3124;
    TH1F* DQM_TrackPull4_X_M3074;
    TH1F* DQM_TrackPull4_Y_M3074;
    TH1F* DQM_TrackPull5_X_M3074;
    TH1F* DQM_TrackPull5_Y_M3074;

    TH2F* DQM_dead_M3124;
    TH2F* DQM_dead_M3074;
    TH2F* DQM_2npeak_M3124;
    TH2F* DQM_2npeak_M3074;
    TH2F* DQM_1speak_M3124;
    TH2F* DQM_1speak_M3074;

    // 3D Tree
    TTree* cluster3DTree ;
    Int_t      tree_runNumber;
    Int_t      tree_lumiSection;
    Int_t      tree_event;
    Int_t      tree_detId;
    Int_t      tree_nclusters;
    Double_t   tree_clusterSize;
    Double_t   tree_clusterSizeX;
    Double_t   tree_clusterSizeY;
    Double_t   tree_x;
    Double_t   tree_y;
    Double_t   tree_z;
    Double_t   tree_x2;
    Double_t   tree_y2;
    Double_t   tree_z2;
    Double_t   tree_xloc;
    Double_t   tree_yloc;
    Double_t   tree_zloc;
    Double_t   tree_xbary;
    Double_t   tree_ybary;
    Double_t   tree_zbary;
    Double_t   tree_charge;
    TString    tree_modName ;
    Long64_t   tree_maxEntries = 100000000 ;

    // 3D Tree
    TTree* TrackTree ;
    Int_t      tree_trackevent;
    Double_t   tree_trackParam0;
    Double_t   tree_trackParam1;
    Double_t   tree_trackParam2;
    Double_t   tree_trackParam3;
    Double_t   tree_trackParam4;
    Double_t   tree_trackParam5;
    Double_t   tree_kx;
    Double_t   tree_ky;
    Double_t   tree_chi2;
    Int_t      tree_npoints;
    Int_t      tree_npointsL;
    Int_t      tree_npointsR;

    // Residuals tree
    TTree* residualsTree;
    Double_t   tree_eventN;
    Double_t   tree_detId2;
    Double_t   tree_biasedResidualX;
    Double_t   tree_biasedResidualY;
    Double_t   tree_unbiasedResidualX;
    Double_t   tree_unbiasedResidualY;


    // detId versus moduleName and other definitions
    std::map<int , TString> detId_to_moduleName ;
    std::map<TString , int> moduleName_to_detID ;
    std::map<TString , std::vector<double> > moduleName_to_position ;
    std::vector<std::pair<int, int> > LayersDefinition;

    // simple implementation of the geometry for home-made tracking
    std::vector<SimplePlane > theSimpleLayersR;
    std::vector<SimplePlane > theSimpleLayersL;
    // test if we want to artificially apply some mis-alignment
    double xmove_par;

    // functions used for home-made tracking
    std::vector<TelescopeTracks> theTeleTrackCollection;

    void doSeeding(edm::Handle<edmNew::DetSetVector<SiPixelCluster> > , const TrackerGeometry *,
                   const PixelClusterParameterEstimator &, edm::ESHandle<TrackerGeometry> );

    void doPatternReco(edm::Handle<edmNew::DetSetVector<SiPixelCluster> > , const TrackerGeometry *,
                       const PixelClusterParameterEstimator &, edm::ESHandle<TrackerGeometry> );

    void ApplyAlignment(int detId, double &xaligned, double &yaligned);
    void ApplyXmove(int detId, double &xaligned, double &yaligned);

    void ComputeResiduals(int detId_to_study, int eventNumber);

    // information about noisy clusters used in home-made tracking
    int noisy_det[1000];
    int noisy_barx[1000];
    int noisy_bary[1000];
    int inoisy;

} ;

/////////////////////
// Class functions //
/////////////////////

SimpleTrackingClean::SimpleTrackingClean( const edm::ParameterSet& iConfig ) : tracksToken_( consumes<TrackCollection>( iConfig.getUntrackedParameter<edm::InputTag>( "tracks" ) ) ) {

    // detId vs Module name
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344200196, "M3090") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344201220, "M3124") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344462340, "M3082") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344463364, "M3175") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344724484, "M3009") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344725508, "M3057") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344986628, "M3027") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(344987652, "M3074") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(352588804, "M3192") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(352589828, "M3204") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(352850948, "M3226") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(352851972, "M3265") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(353113092, "M3023") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(353114116, "M3239") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(353375236, "M3164") ) ;
    detId_to_moduleName.insert( std::pair<uint32_t, TString>(353376260, "M3173") ) ;

    // Module name to detID
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3090", 344200196) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3124", 344201220) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3082", 344462340) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3175", 344463364) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3009", 344724484) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3057", 344725508) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3027", 344986628) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3074", 344987652) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3192", 352588804) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3204", 352589828) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3226", 352850948) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3265", 352851972) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3023", 353113092) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3239", 353114116) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3164", 353375236) ) ;
    moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3173", 353376260) ) ;

    LayersDefinition.push_back(std::pair<int, int> (moduleName_to_detID["M3173"], moduleName_to_detID["M3164"]) );
    LayersDefinition.push_back(std::pair<int, int> (moduleName_to_detID["M3239"], moduleName_to_detID["M3023"]) );
    LayersDefinition.push_back(std::pair<int, int> (moduleName_to_detID["M3265"], moduleName_to_detID["M3226"]) );
    LayersDefinition.push_back(std::pair<int, int> (moduleName_to_detID["M3204"], moduleName_to_detID["M3192"]) );
    LayersDefinition.push_back(std::pair<int, int> (moduleName_to_detID["M3090"], moduleName_to_detID["M3124"]) );
    LayersDefinition.push_back(std::pair<int, int> (moduleName_to_detID["M3082"], moduleName_to_detID["M3175"]) );
    LayersDefinition.push_back(std::pair<int, int> (moduleName_to_detID["M3009"], moduleName_to_detID["M3057"]) );
    LayersDefinition.push_back(std::pair<int, int> (moduleName_to_detID["M3027"], moduleName_to_detID["M3074"]) );

    //layer 3 : left 353376260, right  353375236
    //layer 2 : left  353114116 , right  353113092
    //layer 1 : left 352851972 , right  352850948
    //layer 0 : left   352589828 , right 352588804
    //layer -3 : left 344200196, right  344201220
    //layer -2 : left 344462340 , right 344463364
    //layer -1 : left  344724484 , right 344725508
    //Layer -0 : left 344986628 , right  344987652

    // geometry: definition of the telescope plane for home-made tracking
    // for run 100208 -- old version of SimplePlane initialization
    /*
    theSimpleLayers.push_back( SimplePlane(14.0 + 10.9 + 12.3 + 5.0,  1, 353376260,  353375236));
    theSimpleLayers.push_back( SimplePlane(14.0 + 10.9 + 12.3,        2, 353114116,  353113092));
    theSimpleLayers.push_back( SimplePlane(14.0 + 10.9,               3, 352851972,  352850948));
    theSimpleLayers.push_back( SimplePlane(14.0,                      4, 352589828,  352588804));
    theSimpleLayers.push_back( SimplePlane(-20.5,                     5, 344200196,  344201220));
    theSimpleLayers.push_back( SimplePlane(-20.5 - 5.0,               6, 344462340,  344463364));
    theSimpleLayers.push_back( SimplePlane(-20.5 - 5.0 - 13.5,        7, 344724484,  344725508));
    theSimpleLayers.push_back( SimplePlane(-20.5 - 5.0 - 13.5 - 10.2, 8, 344986628,  344987652));
    */
    // for run 100328
    // separation of the planes: now only 1 module per plane
    theSimpleLayersR.push_back( SimplePlane(50., 0., 17.3 + 5.0 + 5.0 + 5.0,  -30., 20., 1,  353376260) );
    theSimpleLayersR.push_back( SimplePlane(50., 0., 17.3 + 5.0 + 5.0,        -30., 20., 2,  353114116) );
    theSimpleLayersR.push_back( SimplePlane(50., 0., 17.3 + 5.0,              -30., 20., 3,  352851972) );
    theSimpleLayersR.push_back( SimplePlane(50., 0., 17.3,                    -30., 20., 4,  352589828) );
    theSimpleLayersR.push_back( SimplePlane(50., 0., -17.3,                   -30., 20., 5,  344200196) );
    theSimpleLayersR.push_back( SimplePlane(50., 0., -17.3 - 5.0,             -30., 20., 6,  344462340) );
    theSimpleLayersR.push_back( SimplePlane(50., 0., -17.3 - 5.0 - 5.0,       -30., 20., 7,  344724484) );
    theSimpleLayersR.push_back( SimplePlane(50., 0., -17.3 - 5.0 - 5.0 - 5.0, -30., 20., 8,  344986628) );

    // possibility to add some mis-alignment of one plane, configuration in the python
    xmove_par = iConfig.getParameter<double>("xmove_param");
    std::cout << " xmove_par " << xmove_par << std::endl;
    //theSimpleLayersL.push_back( SimplePlane(50., 0., 17.3 + 5.0 + 5.0 + 5.0,  -30., 20., 1,  353375236) );
    theSimpleLayersL.push_back( SimplePlane(50.+xmove_par, 0., 17.3 + 5.0 + 5.0 + 5.0,  -30., 20., 1,  353375236) );
    theSimpleLayersL.push_back( SimplePlane(50., 0., 17.3 + 5.0 + 5.0,        -30., 20., 2,  353113092) );
    theSimpleLayersL.push_back( SimplePlane(50., 0., 17.3 + 5.0,              -30., 20., 3,  352850948) );
    theSimpleLayersL.push_back( SimplePlane(50., 0., 17.3,                    -30., 20., 4,  352588804) );
    theSimpleLayersL.push_back( SimplePlane(50., 0., -17.3,                   -30., 20., 5,  344201220) );
    theSimpleLayersL.push_back( SimplePlane(50., 0., -17.3 - 5.0,             -30., 20., 6,  344463364) );
    theSimpleLayersL.push_back( SimplePlane(50., 0., -17.3 - 5.0 - 5.0,       -30., 20., 7,  344725508) );
    theSimpleLayersL.push_back( SimplePlane(50., 0., -17.3 - 5.0 - 5.0 - 5.0, -30., 20., 8,  344987652) );

    //apply alignment
    for(unsigned int i=0; i<theSimpleLayersR.size(); i++) theSimpleLayersR[i].applyAlignment();
    for(unsigned int i=0; i<theSimpleLayersL.size(); i++) theSimpleLayersL[i].applyAlignment();

    //apply rotation
    for(unsigned int i=0; i<theSimpleLayersR.size(); i++) theSimpleLayersR[i].doRotation();
    for(unsigned int i=0; i<theSimpleLayersL.size(); i++) theSimpleLayersL[i].doRotation();

    TFileDirectory sub1 = fs->mkdir(  "run100000" ); // This does not make sense yet

    cluster3DTree = sub1.make<TTree>("cluster3DTree", "3D Cluster Tree");
    TrackTree = sub1.make<TTree>("trackTree", "parameters of reco tracks");
    residualsTree = sub1.make<TTree>("residualsTree","Tree with residuals after tracking");

    TFileDirectory sub2 = sub1.mkdir( "dqmPlots" ) ;
    TFileDirectory sub3 = sub1.mkdir( "correlationPlots" ) ;

    // Set branch addresses.
    cluster3DTree -> Branch ( "runNumber", &tree_runNumber ) ;
    cluster3DTree -> Branch ( "lumiSection", &tree_lumiSection ) ;
    cluster3DTree -> Branch ( "event", &tree_event ) ;
    cluster3DTree -> Branch ( "detId", &tree_detId ) ;
    cluster3DTree -> Branch ( "modName", &tree_modName ) ;
    cluster3DTree -> Branch ( "nclusters", &tree_nclusters ) ;
    cluster3DTree -> Branch ( "clusterSize", &tree_clusterSize ) ;
    cluster3DTree -> Branch ( "clusterSizeX", &tree_clusterSizeX ) ;
    cluster3DTree -> Branch ( "clusterSizeY", &tree_clusterSizeY ) ;
    cluster3DTree -> Branch ( "x", &tree_x ) ;
    cluster3DTree -> Branch ( "y", &tree_y ) ;
    cluster3DTree -> Branch ( "z", &tree_z ) ;
    cluster3DTree -> Branch ( "x2", &tree_x2 ) ;
    cluster3DTree -> Branch ( "y2", &tree_y2 ) ;
    cluster3DTree -> Branch ( "z2", &tree_z2 ) ;
    cluster3DTree -> Branch ( "xloc", &tree_xloc ) ;
    cluster3DTree -> Branch ( "yloc", &tree_yloc ) ;
    cluster3DTree -> Branch ( "zloc", &tree_zloc ) ;
    cluster3DTree -> Branch ( "xbary", &tree_xbary ) ;
    cluster3DTree -> Branch ( "ybary", &tree_ybary ) ;
    cluster3DTree -> Branch ( "zbary", &tree_zbary ) ;
    cluster3DTree -> Branch ( "clusterCharge", &tree_charge ) ;
    cluster3DTree -> SetCircular ( tree_maxEntries ) ;

    //std::cout << "ST3: Branch addresses set" << std::endl;
    TrackTree -> Branch ( "tree_trackevent",  &tree_trackevent ) ;
    TrackTree -> Branch ( "tree_trackParam0", &tree_trackParam0 ) ;
    TrackTree -> Branch ( "tree_trackParam1", &tree_trackParam1 ) ;
    TrackTree -> Branch ( "tree_trackParam2", &tree_trackParam2 ) ;
    TrackTree -> Branch ( "tree_trackParam3", &tree_trackParam3 ) ;
    TrackTree -> Branch ( "tree_trackParam4", &tree_trackParam4 ) ;
    TrackTree -> Branch ( "tree_trackParam5", &tree_trackParam5 ) ;
    TrackTree -> Branch ( "tree_kX", &tree_kx ) ;
    TrackTree -> Branch ( "tree_kY", &tree_ky ) ;
    TrackTree -> Branch ( "tree_chi2", &tree_chi2 ) ;
    TrackTree -> Branch ( "tree_npoints", &tree_npoints ) ;
    TrackTree -> Branch ( "tree_npointsL", &tree_npointsL ) ;
    TrackTree -> Branch ( "tree_npointsR", &tree_npointsR ) ;

    residualsTree -> Branch ( "eventNumber", &tree_eventN );
    residualsTree -> Branch ( "detectorID", &tree_detId2 );
    residualsTree -> Branch ( "biasedResidualX", &tree_biasedResidualX );
    residualsTree -> Branch ( "biasedResidualY", &tree_biasedResidualY );
    residualsTree -> Branch ( "unbiasedResidualX", &tree_unbiasedResidualX );
    residualsTree -> Branch ( "unbiasedResidualY", &tree_unbiasedResidualY );

    DQM_NumbOfSeeds_per_Event = sub2.make<TH1F>( "DQM_NumbOfSeeds_per_Event", "Number of seeds/event",  30, 0., 30.0  ) ;
    DQM_NumbOfSeeds_per_Event -> GetXaxis ( ) -> SetTitle ( "Number of seeds / event " ) ;
    DQM_NumbOfSeeds_per_Event -> GetYaxis ( ) -> SetTitle ( "Count" ) ;

    DQM_NumbOfTracks_per_Event = sub2.make<TH1F>( "DQM_NumbOfTracks_per_Event", "Number of tracks/event",30, 0., 30.0  ) ;
    DQM_NumbOfTracks_per_Event -> GetXaxis ( ) -> SetTitle ( "Number of tracks / event " ) ;
    DQM_NumbOfTracks_per_Event -> GetYaxis ( ) -> SetTitle ( "Count" ) ;

    DQM_NumbOfClus_per_Event = sub2.make<TH1F>( "DQM_NumbOfClus_per_Event", "Number of clustes/event",30, 0., 30.0  ) ;
    DQM_NumbOfClus_per_Event -> GetXaxis ( ) -> SetTitle ( "Number of clusters / event " ) ;
    DQM_NumbOfClus_per_Event -> GetYaxis ( ) -> SetTitle ( "Count" ) ;

    DQM_NumbOfSeeds_4cl_per_Event = sub2.make<TH1F>( "DQM_NumbOfSeeds_4cl_per_Event", "Number of seeds/event",  30, 0., 30.0  ) ;
    DQM_NumbOfSeeds_4cl_per_Event -> GetXaxis ( ) -> SetTitle ( "Number of seeds / event " ) ;
    DQM_NumbOfSeeds_4cl_per_Event -> GetYaxis ( ) -> SetTitle ( "Count" ) ;

    DQM_NumbOfTracks_4cl_per_Event = sub2.make<TH1F>( "DQM_NumbOfTracks_4cl_per_Event", "Number of tracks/event",30, 0., 30.0  ) ;
    DQM_NumbOfTracks_4cl_per_Event -> GetXaxis ( ) -> SetTitle ( "Number of tracks / event " ) ;
    DQM_NumbOfTracks_4cl_per_Event -> GetYaxis ( ) -> SetTitle ( "Count" ) ;

    distanceR_per_Track = sub2.make<TH1F>(  "distance-R", "distance-R", 100, 0., 20.0  ) ;
    distanceR_per_Track -> GetXaxis ( ) -> SetTitle ( "distance-R (cm)" ) ;
    distanceR_per_Track -> GetYaxis ( ) -> SetTitle ( "Count" ) ;

    deltaX_per_Track = sub2.make<TH1F>(  "deltaX", "deltaX", 200, -20., 20.0  ) ;
    deltaX_per_Track -> GetXaxis ( ) -> SetTitle ( "deltaX (cm)" ) ;
    deltaX_per_Track -> GetYaxis ( ) -> SetTitle ( "Count" ) ;
    deltaY_per_Track = sub2.make<TH1F>(  "deltaY", "deltaY", 200, -20., 20.0  ) ;
    deltaY_per_Track -> GetXaxis ( ) -> SetTitle ( "deltaY (cm)" ) ;
    deltaY_per_Track -> GetYaxis ( ) -> SetTitle ( "Count" ) ;

    DQM_TrackPull4_X_M3124 = sub2.make<TH1F>(  "DQM_TrackPull4_X_M3124",  "pull track X for M3124" , 500, -0.05, 0.05 );
    DQM_TrackPull4_Y_M3124 = sub2.make<TH1F>(  "DQM_TrackPull4_Y_M3124",  "pull track Y for M3124" , 500, -0.05, 0.05 );
    DQM_TrackPull5_X_M3124 = sub2.make<TH1F>(  "DQM_TrackPull5_X_M3124",  "pull track X for M3124" , 500, -0.05, 0.05 );
    DQM_TrackPull5_Y_M3124 = sub2.make<TH1F>(  "DQM_TrackPull5_Y_M3124",  "pull track Y for M3124" , 500, -0.05, 0.05 );
    DQM_TrackPull4_X_M3074 = sub2.make<TH1F>(  "DQM_TrackPull4_X_M3074",  "pull track X for M3074" , 500, -0.05, 0.05 );
    DQM_TrackPull4_Y_M3074 = sub2.make<TH1F>(  "DQM_TrackPull4_Y_M3074",  "pull track Y for M3074" , 500, -0.05, 0.05 );
    DQM_TrackPull5_X_M3074 = sub2.make<TH1F>(  "DQM_TrackPull5_X_M3074",  "pull track X for M3074" , 500, -0.05, 0.05 );
    DQM_TrackPull5_Y_M3074 = sub2.make<TH1F>(  "DQM_TrackPull5_Y_M3074",  "pull track Y for M3074" , 500, -0.05, 0.05 );

    DQM_dead_M3124 = sub2.make<TH2F>(  "Track_posX_dead_M3124",  "Track_posY_dead_M3124" , 100, -2., 0. , 100, -3.5,3.5);
    DQM_dead_M3074 = sub2.make<TH2F>(  "Track_posX_dead_M3074",  "Track_posY_dead_M3074" , 100, -2., 0. , 100, -3.5,3.5);
    DQM_2npeak_M3124 = sub2.make<TH2F>(  "Track_posX_2npeak_M3124",  "Track_posY_2npeak_M3124" , 100, -2., 0. , 100, -3.5,3.5);
    DQM_2npeak_M3074 = sub2.make<TH2F>(  "Track_posX_2npeak_M3074",  "Track_posY_2npeak_M3074" , 100, -2., 0. , 100, -3.5,3.5);
    DQM_1speak_M3124 = sub2.make<TH2F>(  "Track_posX_1speak_M3124",  "Track_posY_1speak_M3124" , 100, -2., 0. , 100, -3.5,3.5);
    DQM_1speak_M3074 = sub2.make<TH2F>(  "Track_posX_1speak_M3074",  "Track_posY_1speak_M3074" , 100, -2., 0. , 100, -3.5,3.5);

    for ( std::map<int, TString>::iterator it = detId_to_moduleName.begin(); it != detId_to_moduleName.end(); it++ ) {

        TString modulename = it -> second ;

        std::vector<TH2F *> tmp_vec_x ; // for the corr plots, hence the extra for loop
        std::vector<TH2F *> tmp_vec_y ; // for the corr plots, hence the extra for loop

        // Make the correlation plots
        for ( std::map<int, TString>::iterator jt = it; jt != detId_to_moduleName.end(); jt++ ) { // jt=it to make sure we do not have double plots.

            TString modulename0 = jt -> second ;

            TH2F* DQM_Correlation_X_tmp = sub3.make<TH2F>( ( "DQM_Correlation_X_" + modulename + "_" + modulename0).Data(), ( "X-Correlation between " + modulename + " and " + modulename0 ).Data(), 160., 0., 160., 160., 0., 160. ) ;
            TH2F* DQM_Correlation_Y_tmp = sub3.make<TH2F>( ( "DQM_Correlation_Y_" + modulename + "_" + modulename0).Data(), ( "Y-Correlation between " + modulename + " and " + modulename0 ).Data(), 416., 0., 416., 416., 0., 416. ) ;
            TH2F* DQM_Correlation2_X_tmp = sub3.make<TH2F>( ( "DQM_Correlation2_X_" + modulename + "_" + modulename0).Data(), ( "X-Correlation between " + modulename + " and " + modulename0 ).Data(), 88, 47.8, 52.2, 88, 47.8, 52.2 ) ;
            TH2F* DQM_Correlation2_Y_tmp = sub3.make<TH2F>( ( "DQM_Correlation2_Y_" + modulename + "_" + modulename0).Data(), ( "Y-Correlation between " + modulename + " and " + modulename0 ).Data(), 60, -3., 3., 60, -3., 3. ) ;
            TH2F* DQM_Correlation3_X_tmp = sub3.make<TH2F>( ( "DQM_Correlation3_X_" + modulename + "_" + modulename0).Data(), ( "X-Correlation between " + modulename + " and " + modulename0 ).Data(), 160., 0., 160., 160., 0., 160. ) ;
            TH2F* DQM_Correlation3_Y_tmp = sub3.make<TH2F>( ( "DQM_Correlation3_Y_" + modulename + "_" + modulename0).Data(), ( "Y-Correlation between " + modulename + " and " + modulename0 ).Data(), 416., 0., 416., 416., 0., 416. ) ;

            DQM_Correlation_X_tmp->GetXaxis()->SetTitle("x_" + modulename) ;
            DQM_Correlation_X_tmp->GetYaxis()->SetTitle("x_" + modulename0) ;
            DQM_Correlation_Y_tmp->GetXaxis()->SetTitle("y_" + modulename) ;
            DQM_Correlation_Y_tmp->GetYaxis()->SetTitle("y_" + modulename0) ;

            DQM_Correlation2_X_tmp->GetXaxis()->SetTitle("x_" + modulename) ;
            DQM_Correlation2_X_tmp->GetYaxis()->SetTitle("x_" + modulename0) ;
            DQM_Correlation2_Y_tmp->GetXaxis()->SetTitle("y_" + modulename) ;
            DQM_Correlation2_Y_tmp->GetYaxis()->SetTitle("y_" + modulename0) ;

            DQM_Correlation3_X_tmp->GetXaxis()->SetTitle("x_" + modulename) ;
            DQM_Correlation3_X_tmp->GetYaxis()->SetTitle("x_" + modulename0) ;
            DQM_Correlation3_Y_tmp->GetXaxis()->SetTitle("y_" + modulename) ;
            DQM_Correlation3_Y_tmp->GetYaxis()->SetTitle("y_" + modulename0) ;

            std::pair<uint32_t, uint32_t> modulePair = std::make_pair ( it->first, jt->first ) ;

            DQM_Correlation_X.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation_X_tmp ) ) ;
            DQM_Correlation_Y.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation_Y_tmp ) ) ;

            DQM_Correlation2_X.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation2_X_tmp ) ) ;
            DQM_Correlation2_Y.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation2_Y_tmp ) ) ;

            DQM_Correlation3_X.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation3_X_tmp ) ) ;
            DQM_Correlation3_Y.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation3_Y_tmp ) ) ;

        }//end for j

        // Make the DQM plots
        TH1F* DQM_PixelNoise_tmp = sub2.make<TH1F>( ( "DQM_PixelNoise_" + modulename ).Data(), ( "Pixel hit distribution " + modulename).Data(), 66560, 0., 66560. ) ;

        TH1F* DQM_ClusterCharge_tmp = sub2.make<TH1F>( ( "DQM_ClusterCharge_" + modulename ).Data(), ( "Cluster charge for " + modulename).Data(), 100, 0., 100000. );
        TH1F* DQM_ClusterSize_X_tmp = sub2.make<TH1F>( ( "DQM_ClusterSize_X_" + modulename ).Data(), ( "X cluster size for " + modulename).Data(), 30, 0., 30. );
        TH1F* DQM_ClusterSize_Y_tmp = sub2.make<TH1F>( ( "DQM_ClusterSize_Y_" + modulename ).Data(), ( "Y cluster size for " + modulename ).Data(), 30, 0., 30. ) ;
        TH1F* DQM_ClusterSize_XY_tmp = sub2.make<TH1F>( ( "DQM_ClusterSize_XY_" + modulename ).Data(), ( "Cluster Size for "  + modulename ).Data(), 30, 0., 30. ) ;
        TH1F* DQM_NumbOfClusters_per_Event_tmp = sub2.make<TH1F>( ("DQM_NumbOfClusters_per_Event_" + modulename).Data(), ("number of clusters for "  + modulename).Data(), 30, 0., 30. );
        TH2F* DQM_ClusterPosition_tmp = sub2.make<TH2F>( ( "DQM_ClusterPosition_" + modulename ).Data(), ( "Cluster occupancy per col per row for " + modulename ).Data(), 416, 0. - 0.5, 416. - 0.5, 160, 0. - 0.5, 160. - 0.5 );
        TH1F* DQM_Chi2_tmp = sub2.make<TH1F>( ( "DQM_Chi2_" + modulename ).Data(), ( "Chi2 for tracks passing by " + modulename).Data(), 100, 0., 20. );
        TH1F* DQM2_Chi2_tmp = sub2.make<TH1F>( ( "DQM2_Chi2_" + modulename ).Data(), ( "Chi2 for tracks passing by " + modulename).Data(), 100, 0., 20. );
        TH1F* DQM2_Sum2_tmp = sub2.make<TH1F>( ( "DQM2_Sum2_" + modulename ).Data(), ( "Sum2 for tracks passing by " + modulename).Data(), 100, 0., 0.05 );

        DQM_PixelNoise_tmp -> GetXaxis() -> SetTitle("PixelNumber");
        DQM_PixelNoise_tmp -> GetYaxis() -> SetTitle("Number of hits per pixel");

        DQM_ClusterCharge_tmp -> GetXaxis ( ) ->SetTitle ( "Charge (electrons)" ) ;
        DQM_ClusterSize_X_tmp -> GetXaxis ( ) ->SetTitle ( "size (pixels)" ) ;
        DQM_ClusterSize_Y_tmp -> GetXaxis ( ) ->SetTitle ( "size (pixels)" ) ;
        DQM_ClusterSize_XY_tmp -> GetXaxis ( ) ->SetTitle ( "size (pixels)" ) ;
        DQM_ClusterPosition_tmp -> GetXaxis ( ) ->SetTitle ( "col" ) ;
        DQM_ClusterPosition_tmp -> GetYaxis ( ) ->SetTitle ( "row" ) ;
        DQM_NumbOfClusters_per_Event_tmp -> GetXaxis ( ) -> SetTitle ( "Number of clusters / event " ) ;
        DQM_NumbOfClusters_per_Event_tmp -> GetYaxis ( ) -> SetTitle ( "Count" ) ;
        DQM_Chi2_tmp -> GetXaxis ( ) ->SetTitle ( "Chi2" ) ;
        DQM2_Chi2_tmp -> GetXaxis ( ) ->SetTitle ( "Chi2" ) ;
        DQM2_Sum2_tmp -> GetXaxis ( ) ->SetTitle ( "Sum2" ) ;

        DQM_PixelNoise.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_PixelNoise_tmp ) ) ;
        DQM_ClusterCharge.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterCharge_tmp ) ) ;
        DQM_ClusterSize_X.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_X_tmp ) ) ;
        DQM_ClusterSize_Y.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_Y_tmp ) ) ;
        DQM_ClusterSize_XY.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_XY_tmp ) ) ;
        DQM_NumbOfClusters_per_Event.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_NumbOfClusters_per_Event_tmp ) );
        DQM_ClusterPosition.insert ( std::pair< uint32_t, TH2F* >( it->first, DQM_ClusterPosition_tmp ) ) ;
        DQM_Chi2.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_Chi2_tmp ) ) ;
        DQM2_Chi2.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM2_Chi2_tmp ) ) ;
        DQM2_Sum2.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM2_Sum2_tmp ) ) ;

        TH1F* DQM_TrackPull_X_tmp = sub2.make<TH1F>( ( "DQM_TrackPull_X_" + modulename ).Data(), ( "pull track X for " + modulename).Data(), 500, -0.5, 0.5 );
        TH1F* DQM_TrackPull_Y_tmp = sub2.make<TH1F>( ( "DQM_TrackPull_Y_" + modulename ).Data(), ( "pull track Y for " + modulename).Data(), 500, -0.5, 0.5 );
        TH1F* DQM_TrackPull_X_tmp2 = sub2.make<TH1F>( ( "DQM_TrackPull2_X_" + modulename ).Data(), ( "pull track X for " + modulename).Data(), 500, -0.05, 0.05 );
        TH1F* DQM_TrackPull_Y_tmp2 = sub2.make<TH1F>( ( "DQM_TrackPull2_Y_" + modulename ).Data(), ( "pull track Y for " + modulename).Data(), 500, -0.05, 0.05 );

        TH1F* DQM_TrackPull_X_tmp3 = sub2.make<TH1F>( ( "DQM_TrackPull3_X_" + modulename ).Data(), ( "pull track X for " + modulename).Data(), 500, -0.05, 0.05 );
        TH1F* DQM_TrackPull_Y_tmp3 = sub2.make<TH1F>( ( "DQM_TrackPull3_Y_" + modulename ).Data(), ( "pull track Y for " + modulename).Data(), 500, -0.05, 0.05 );

        TH2F* DQM_Corr_dX_X_tmp = sub2.make<TH2F>( ( "DQM_Corr_dX_X_" + modulename ).Data(), ( "dXloc vs Xloc for " + modulename).Data(), 32, -1.6, 1.6, 500, -0.1, 0.1 ) ;
        TH2F* DQM_Corr_dY_Y_tmp = sub2.make<TH2F>( ( "DQM_Corr_dY_Y_" + modulename ).Data(), ( "dYloc vs Yloc for " + modulename).Data(), 64, -3.2, 3.2, 500, -0.1, 0.1 ) ;

        TH1F* DQM2_TrackPull_X_tmp = sub2.make<TH1F>( ( "DQM2_TrackPull_X_" + modulename ).Data(), ( "pull track X for " + modulename).Data(), 500, -0.5, 0.5 );
        TH1F* DQM2_TrackPull_Y_tmp = sub2.make<TH1F>( ( "DQM2_TrackPull_Y_" + modulename ).Data(), ( "pull track Y for " + modulename).Data(), 500, -0.5, 0.5 );

        DQM_TrackPull_X[it->first] = DQM_TrackPull_X_tmp;
        DQM_TrackPull_Y[it->first] = DQM_TrackPull_Y_tmp;
        DQM_TrackPull2_X[it->first] = DQM_TrackPull_X_tmp2;
        DQM_TrackPull2_Y[it->first] = DQM_TrackPull_Y_tmp2;
        DQM_TrackPull3_X[it->first] = DQM_TrackPull_X_tmp3;
        DQM_TrackPull3_Y[it->first] = DQM_TrackPull_Y_tmp3;
        DQM_Corr_dX_X[it->first] = DQM_Corr_dX_X_tmp;
        DQM_Corr_dY_Y[it->first] = DQM_Corr_dY_Y_tmp;
        DQM2_TrackPull_X[it->first] = DQM2_TrackPull_X_tmp;
        DQM2_TrackPull_Y[it->first] = DQM2_TrackPull_Y_tmp;

    }//end for it

    pixeldigiToken_    = consumes<edm::DetSetVector<PixelDigi> >        (iConfig.getParameter<edm::InputTag>("PixelDigisLabel"))   ;
    pixelclusterToken_ = consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("PixelClustersLabel"));
    pixelhitToken_     = consumes<edmNew::DetSetVector<SiPixelRecHit> > (iConfig.getParameter<edm::InputTag>("PixelHitsLabel"))    ;

}//end SimpleTrackingClean()

SimpleTrackingClean::~SimpleTrackingClean()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void SimpleTrackingClean::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

    using namespace edm;

    EventID myEvId = iEvent.id();

    //get collection of digi
    edm::Handle<edm::DetSetVector<PixelDigi> > pixeldigis;
    iEvent.getByToken(pixeldigiToken_,pixeldigis  );

    //get collection of cluster
    edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelclusters;
    iEvent.getByToken(pixelclusterToken_,pixelclusters  );

    //get collection or RecHits
    edm::Handle< edmNew::DetSetVector<SiPixelRecHit> > pixelhits;
    iEvent.getByToken(pixelhitToken_,pixelhits  );

    // Get the geometry of the tracker for converting the LocalPoint to a GlobalPoint
    edm::ESHandle<TrackerGeometry> tracker;
    iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
    const TrackerGeometry *tkgeom = &(*tracker);

    // // Choose the CPE Estimator that will be used to estimate the LocalPoint of the cluster
    edm::ESHandle<PixelClusterParameterEstimator> cpEstimator;
    iSetup.get<TkPixelCPERecord>().get("PixelCPEGeneric", cpEstimator);
    const PixelClusterParameterEstimator &cpe(*cpEstimator);

    //---------------------------------
    //loop on digis
    //---------------------------------

    for( edm::DetSetVector<PixelDigi>::const_iterator DSViter=pixeldigis->begin(); DSViter!=pixeldigis->end(); DSViter++   ) {

        edm::DetSet<PixelDigi>::const_iterator begin=DSViter->begin();
        edm::DetSet<PixelDigi>::const_iterator end  =DSViter->end();

        //auto id = DetId(DSViter->detId());

        std::map< int, int > numHits_per_Channel ;

        for(edm::DetSet<PixelDigi>::const_iterator iter=begin; iter!=end; ++iter) {

            // First get the channel number of the hits
            //int channel = iter -> channel ( ) ;
            int rowN = iter -> row () ;
            int colN = iter -> column () ;
            auto id = DetId(DSViter->detId());
            DQM_PixelNoise[ id.rawId() ] -> Fill ( rowN*colN ) ;

            int channel = 0 ;
            if ( rowN == 0 ) channel = colN ;
            else channel = 160 + rowN * colN ;

            // Add one to the number of hits of the right channel
            if ( numHits_per_Channel.count ( channel ) > 0 ) {
                numHits_per_Channel[ channel ] ++ ;
            } else {
                numHits_per_Channel.insert ( std::pair< int, int >( channel, 1 ) ) ;
            }//end if channel

        }//end for Digis

    }//end for Detectors

    //---------------------------------
    // Loop over the clusters
    //---------------------------------

    for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end(); DSViter++   ) {

        edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
        edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();

        auto id = DetId(DSViter->detId());

        for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter2=pixelclusters->begin(); DSViter2!=pixelclusters->end(); DSViter2++   ) {

            edmNew::DetSet<SiPixelCluster>::const_iterator begin2=DSViter2->begin();
            edmNew::DetSet<SiPixelCluster>::const_iterator end2  =DSViter2->end();

            auto id2 = DetId(DSViter2->detId());

            // test
            int nclus_check1=0;
            int nclus_check2=0;
            for(edmNew::DetSet<SiPixelCluster>::const_iterator iter=begin; iter!=end; ++iter) {
                nclus_check1++;
            }
            for(edmNew::DetSet<SiPixelCluster>::const_iterator iter2=begin2; iter2!=end2; ++iter2) {
                nclus_check2++;
            }

            for(edmNew::DetSet<SiPixelCluster>::const_iterator iter=begin; iter!=end; ++iter) {

                float x = iter->x();                   // barycenter x position
                float y = iter->y();                   // barycenter y position
                int row = x , col = y ;

                const PixelGeomDetUnit *pixdet1 = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(id);
                LocalPoint lp1(-9999., -9999., -9999.);
                PixelClusterParameterEstimator::ReturnType params1 = cpe.getParameters(*iter,*pixdet1);
                lp1 = std::get<0>(params1);
                const Surface& surface1 = tracker->idToDet(id)->surface();
                GlobalPoint gp1 = surface1.toGlobal(lp1);
                double xgp1=0, ygp1=0;
                xgp1 = gp1.x();
                ygp1 = gp1.y();

                for(edmNew::DetSet<SiPixelCluster>::const_iterator iter2=begin2; iter2!=end2; ++iter2) {

                    float x2 = iter2->x();                   // barycenter x position
                    float y2 = iter2->y();                   // barycenter y position
                    int row2 = x2, col2 = y2 ;

                    std::pair<uint32_t, uint32_t> modulePair = std::make_pair ( id.rawId(), id2.rawId() ) ;

                    auto itHistMap = DQM_Correlation_X.find(modulePair);

                    if ( itHistMap == DQM_Correlation_X.end() ) continue;
                    itHistMap->second->Fill ( row, row2 ) ;

                    auto itHistMap2 = DQM_Correlation_Y.find(modulePair);
                    if ( itHistMap2 == DQM_Correlation_Y.end() ) continue;
                    itHistMap2->second->Fill ( col, col2 ) ;

                    const PixelGeomDetUnit *pixdet2 = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(id2);
                    LocalPoint lp2(-9999., -9999., -9999.);
                    PixelClusterParameterEstimator::ReturnType params2 = cpe.getParameters(*iter2,*pixdet2);
                    lp2 = std::get<0>(params2);
                    const Surface& surface2 = tracker->idToDet(id2)->surface();
                    GlobalPoint gp2 = surface2.toGlobal(lp2);

                    double xgp2=0, ygp2=0;
                    xgp2 = gp2.x();
                    ygp2 = gp2.y();
                    auto itHistMap3 = DQM_Correlation2_X.find(modulePair);
                    if ( itHistMap3 == DQM_Correlation2_X.end() ) continue;
                    itHistMap3->second->Fill ( xgp1, xgp2 ) ;

                    auto itHistMap4 = DQM_Correlation2_Y.find(modulePair);
                    if ( itHistMap4 == DQM_Correlation2_Y.end() ) continue;
                    itHistMap4->second->Fill ( ygp1, ygp2 ) ;

                    if (nclus_check1==1 && nclus_check2==1) {
                        auto itHistMap5 = DQM_Correlation3_X.find(modulePair);
                        if ( itHistMap5 == DQM_Correlation3_X.end() ) continue;
                        itHistMap5->second->Fill ( row, row2 ) ;

                        auto itHistMap6 = DQM_Correlation3_Y.find(modulePair);
                        if ( itHistMap6 == DQM_Correlation3_Y.end() ) continue;
                        itHistMap6->second->Fill ( col, col2 ) ;
                    }//end if cluster check
                }//end for clusters2
            }//end for cluster
        }//end for first detectors2
    }//end for first detectors

    // Counting the number of clusters for this detector and this event
    std::map< uint32_t, int > clustersPerModule ;
    for ( std::map<int, TString>::iterator it = detId_to_moduleName.begin(); it != detId_to_moduleName.end(); it++ ) clustersPerModule.insert( std::pair< uint32_t, int >( it->first, 0 ) ) ;

    for ( edmNew::DetSetVector< SiPixelCluster >::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end(); DSViter++ ) {

        edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
        edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();

        auto id = DetId(DSViter->detId());

        for ( edmNew::DetSet< SiPixelCluster >::const_iterator iter=begin; iter!=end; ++iter ) {

            float x = iter->x();                   // barycenter x position
            float y = iter->y();                   // barycenter y position
            int size = iter->size();               // total size of cluster (in pixels)
            int sizeX = iter->sizeX();             // size of cluster in x-iterrection
            int sizeY = iter->sizeY();             // size of cluster in y-iterrection

            int row = x-0.5, col = y -0.5;
            DQM_ClusterCharge[ id.rawId() ] -> Fill ( iter->charge() ) ;
            DQM_ClusterSize_X[ id.rawId() ] -> Fill ( sizeX ) ;
            DQM_ClusterSize_Y[ id.rawId() ] -> Fill ( sizeY ) ;
            DQM_ClusterSize_XY[ id.rawId() ] -> Fill ( size ) ;
            DQM_ClusterPosition[ id.rawId() ] -> Fill ( col, row ) ;

            clustersPerModule[ id.rawId() ] ++ ;

        }//end for clusters in detector
    }//end for detectors

    for ( std::map<uint32_t, int>::iterator it = clustersPerModule.begin(); it != clustersPerModule.end(); it++ ) DQM_NumbOfClusters_per_Event[ it->first ] -> Fill ( it->second ) ;

    //---------------------------------
    // Loop over the hits
    //---------------------------------

    /*  for ( edmNew::DetSetVector< SiPixelRecHit >::const_iterator DSViter=pixelhits->begin(); DSViter!=pixelhits->end(); DSViter++ ) {

        edmNew::DetSet<SiPixelRecHit>::const_iterator begin=DSViter->begin();
        edmNew::DetSet<SiPixelRecHit>::const_iterator end  =DSViter->end();

        for ( edmNew::DetSet< SiPixelRecHit >::const_iterator iter=begin; iter!=end; ++iter ) {

          // Here we do something with the hits.

        }//end for DetSet
      }//end for DetSetVector
    */

    //---------------------------------
    // Fill the 3D tree
    //---------------------------------

    std::vector<TVector3> vectClust;
    std::vector<TVector3> vectClust_err;

    int iclu_count=0;
    for ( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end(); DSViter++ ) {

        edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
        edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();

        // Surface of (detId) and use the surface to convert 2D to 3D
        auto detId = DetId(DSViter->detId());
        int iCluster=0;
        const PixelGeomDetUnit *pixdet = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(detId);
        LocalPoint lp(-9999., -9999., -9999.);

        // Then loop on the clusters of the module
        for ( edmNew::DetSet<SiPixelCluster>::const_iterator itCluster=begin; itCluster!=end; ++itCluster ) {

            PixelClusterParameterEstimator::ReturnType params = cpe.getParameters(*itCluster,*pixdet);
            lp = std::get<0>(params);

            const Surface& surface = tracker->idToDet(detId)->surface();
            GlobalPoint gp = surface.toGlobal(lp);

            double x=0, y=0, z=0;
            x = gp.x();
            y = gp.y();
            z = gp.z();


            // check the noisy clusters
            bool noisyflag=false;
            for (int in=0; in<inoisy; in++) {
                if (int(detId)==noisy_det[in]) {
                    double xbary=itCluster->x();
                    double ybary=itCluster->y();
                    if (int(xbary)==noisy_barx[in] && int(ybary)==noisy_bary[in]) {
                        noisyflag=true;
                    }
                }
            }
            if (!noisyflag) iclu_count++;
            // end check

            double ex=0.05, ey=0.05, ez=5;

            TVector3 tmpGP(x, y, z);
            TVector3 tmpGP_err(ex, ey, ez);
            vectClust.push_back(tmpGP);
            vectClust_err.push_back(tmpGP_err);

            // Let's fill in the tree
            tree_runNumber = myEvId.run();
            tree_lumiSection = myEvId.luminosityBlock();
            tree_event = myEvId.event();
            tree_detId = detId;
            tree_nclusters = iCluster++;
            tree_clusterSize = itCluster->size();
            tree_clusterSizeX = itCluster->sizeX();
            tree_clusterSizeY = itCluster->sizeY();
            tree_charge = itCluster->charge();
            tree_x = x;
            tree_y = y;
            tree_z = z;
            tree_modName = detId_to_moduleName[detId];

            tree_xbary= itCluster->x();                   // barycenter x position
            tree_ybary= itCluster->y();                   // barycenter y position
            tree_zbary= 0;

            // compute the local positions according to theSimpleLayersR or theSimpleLayersL definitions
            int ilayernumber=-1;
            bool isleft=false;
            bool isright=false;
            double x2=0,y2=0,z2=0;
            for (int il1=0; il1<8; il1++) {
                //if (int(detId) == LayersDefinition[il1].first || int(detId) == LayersDefinition[il1].second) ilayernumber=il1;
                if (int(detId) == LayersDefinition[il1].first)  {
                    ilayernumber=il1;
                    isright=true;
                }
                else if (int(detId) == LayersDefinition[il1].second) {
                    ilayernumber=il1;
                    isleft=true;
                }

            }
            if (ilayernumber>-1) {
                TVector3 gptemp(x,y,z);
                TVector3 lptemp;
                //lptemp= theSimpleLayers[ilayernumber].getLocalPointPosition(gptemp);
                if (isright) lptemp= theSimpleLayersR[ilayernumber].getLocalPointPosition(gptemp);
                else if (isleft) lptemp= theSimpleLayersL[ilayernumber].getLocalPointPosition(gptemp);
                x2= lptemp.x();
                y2= lptemp.y();
                z2= lptemp.z();
            }
            tree_x2 = x2;
            tree_y2 = y2;
            tree_z2 = z2;

            // store the local position from CMSSW
            tree_xloc=lp.x();
            tree_yloc=lp.y();
            tree_zloc=lp.z();

            /*
                  // check the transformatio of the CMSSW local positions to the local positions according to theSimpleLayersR and theSimpleLayersL definitions
                  if (x-0.03>50-1.1*(y+0.03)/5.6 && x+0.03<51.6-1.1*(y-0.03)/5.6) { // some margin in case of misalignment
                    tree_xloc+=0.81;
                  }
                  else {
                    tree_xloc+=0.81;
                    tree_xloc*=-1.;
                  }
            */
            cluster3DTree->Fill();

            //if(int(detId) == 344200196 || int(detId) == 344201220) testModuleLayer->Fill( gp.x(),  gp.y()  );
        } //end for clusters of the first detector
    } //end for first detectors
    DQM_NumbOfClus_per_Event-> Fill(iclu_count) ;

    //**********************************************************************
    //**********************************************************************
    //**************** This is for Home-Made Tracking***********************
    //**********************************************************************
    //**********************************************************************
    //**********************************************************************


    //**************** Get collection of "Seeds" **************************
    //********  From pair of clusters in the 2 first layers **************

    doSeeding( pixelclusters, tkgeom, cpe,  tracker);
//  std::cout << "number of seeds " << theTeleTrackCollection.size() << std::endl;

    /*
      for(unsigned int itrack = 0; itrack< theTeleTrackCollection.size(); itrack++){
        std::vector<int > modulesID = theTeleTrackCollection[itrack].getModDetId();
        std::vector<TVector3 > theGP = theTeleTrackCollection[itrack].getGlobalPoints();
        //std::cout << "truc " << std::endl;
        for(unsigned int iHit=0; iHit < theGP.size(); iHit++){
          //std::cout << "det ID " << modulesID[iHit] << "  x " <<  theGP[iHit].X() << "  y " << theGP[iHit].Y()<< "  z " << theGP[iHit].Z()<< std::endl;
        }
      }
    */

    //********  Add clusters on the next layer to build the track ***********
    doPatternReco( pixelclusters, tkgeom, cpe,  tracker);
    //std::cout << "number of tracks " << theTeleTrackCollection.size() << std::endl;

    DQM_NumbOfSeeds_per_Event-> Fill(theTeleTrackCollection.size()) ;
    if (iclu_count>=4) DQM_NumbOfSeeds_4cl_per_Event-> Fill(theTeleTrackCollection.size()) ;
    int icount_trk=0;

    // loop over the tracks
    for(unsigned int itrack = 0; itrack< theTeleTrackCollection.size(); itrack++) {

        std::vector<int > modulesID = theTeleTrackCollection[itrack].getModDetId();
        std::vector<TVector3 > theGP = theTeleTrackCollection[itrack].getGlobalPoints();
        std::vector<TVector3 > theLP = theTeleTrackCollection[itrack].getPseudoLocalPoints();
        std::vector<double > plane_paramA = theTeleTrackCollection[itrack].getAssPlaneA();
        std::vector<double > plane_paramB = theTeleTrackCollection[itrack].getAssPlaneB();
        std::vector<double > plane_paramC = theTeleTrackCollection[itrack].getAssPlaneC();
        std::vector<double > plane_paramD = theTeleTrackCollection[itrack].getAssPlaneD();

        //if(theGP.size() < 4) continue;   // remove too short tracks
        icount_trk++;
        //if(theTeleTrackCollection[itrack].getChi2() > 3) continue;   // temporary check by Caroline
        //std::cout << "********** new track ******* " << std::endl;
        /*
            std::cout << "parameters " << theTeleTrackCollection[itrack].getParameter(0)<< " " <<
            theTeleTrackCollection[itrack].getParameter(1) << " " <<
            theTeleTrackCollection[itrack].getParameter(2)<< " " <<
            theTeleTrackCollection[itrack].getParameter(3)<< " " <<
            theTeleTrackCollection[itrack].getParameter(4)<< " " <<
            theTeleTrackCollection[itrack].getParameter(5)<<  std::endl;
            std::cout << " chi2 " << theTeleTrackCollection[itrack].getChi2() <<  " norm chi2 " << theTeleTrackCollection[itrack].getNormChi2() <<  std::endl;
        */

        tree_trackevent = myEvId.event();
        tree_trackParam0=theTeleTrackCollection[itrack].getParameter(0);
        tree_trackParam1=theTeleTrackCollection[itrack].getParameter(1);
        tree_trackParam2=theTeleTrackCollection[itrack].getParameter(2);
        tree_trackParam3=theTeleTrackCollection[itrack].getParameter(3);
        tree_trackParam4=theTeleTrackCollection[itrack].getParameter(4);
        tree_trackParam5=theTeleTrackCollection[itrack].getParameter(5);
        tree_kx=tree_trackParam1/tree_trackParam5;
        tree_ky=tree_trackParam3/tree_trackParam5;
        tree_chi2=theTeleTrackCollection[itrack].getChi2();
        tree_npoints=theGP.size();
        tree_npointsL=0;
        tree_npointsR=0;
        for(unsigned int iHit=0; iHit < theGP.size(); iHit++) {
            for (int il1=0; il1<8; il1++) {
                if (modulesID[iHit] == LayersDefinition[il1].first)  {
                    tree_npointsR++;
                } else if (modulesID[iHit] == LayersDefinition[il1].second) {
                    tree_npointsL++;
                }
            }
        }

        TrackTree->Fill();

        for(unsigned int iHit=0; iHit < theGP.size(); iHit++) {
            //std::cout << "modulesID " << detId_to_moduleName[modulesID[iHit]] << std::endl;
            double xtemp, ytemp, ztemp;
            double parFit[6] = {theTeleTrackCollection[itrack].getParameter(0), theTeleTrackCollection[itrack].getParameter(1), theTeleTrackCollection[itrack].getParameter(2), theTeleTrackCollection[itrack].getParameter(3), theTeleTrackCollection[itrack].getParameter(4), theTeleTrackCollection[itrack].getParameter(5)};

            //theTeleTrackCollection[itrack].line(theGP[iHit].Z(), parFit, xtemp, ytemp, ztemp);
            double planeqz[4];


            bool isleft=false;
            bool isright=false;
            int ilayernumber=-1;
            for (int il1=0; il1<8; il1++) {
                //if (modulesID[iHit] == LayersDefinition[il1].first || modulesID[iHit] == LayersDefinition[il1].second) ilayernumber=il1;
                if (modulesID[iHit] == LayersDefinition[il1].first)  {
                    ilayernumber=il1;
                    isright=true;
                }
                else if (modulesID[iHit] == LayersDefinition[il1].second) {
                    ilayernumber=il1;
                    isleft=true;
                }
            }

            if (ilayernumber>-1) {
                planeqz[0]= plane_paramA[iHit];
                planeqz[1]= plane_paramB[iHit];
                planeqz[2]= plane_paramC[iHit];
                planeqz[3]= plane_paramD[iHit];
                theTeleTrackCollection[itrack].intersection(planeqz, parFit, xtemp, ytemp, ztemp);
                TVector3 gptemp(xtemp, ytemp, ztemp);
                TVector3 lptemp;
                if (isright) lptemp= theSimpleLayersR[ilayernumber].getLocalPointPosition(gptemp);
                else if (isleft) lptemp= theSimpleLayersL[ilayernumber].getLocalPointPosition(gptemp);
                DQM_TrackPull_X[modulesID[iHit]]->Fill(lptemp.X() - theLP[iHit].X());
                DQM_TrackPull_Y[modulesID[iHit]]->Fill(lptemp.Y() - theLP[iHit].Y());
                DQM_TrackPull2_X[modulesID[iHit]]->Fill(lptemp.X() - theLP[iHit].X());
                DQM_TrackPull2_Y[modulesID[iHit]]->Fill(lptemp.Y() - theLP[iHit].Y());
                DQM_Chi2[modulesID[iHit]]->Fill(theTeleTrackCollection[itrack].getChi2());

                if (-0.6 < lptemp.X() && lptemp.X() < 0.2 && -1.1 < lptemp.Y() && lptemp.Y() < -0.3) {
                    // test central position of the beam
                    DQM_TrackPull3_X[modulesID[iHit]]->Fill(lptemp.X() - theLP[iHit].X());
                    DQM_TrackPull3_Y[modulesID[iHit]]->Fill(lptemp.Y() - theLP[iHit].Y());
                }
                // some test to understand the impact of the inactive area on the fit
                if (modulesID[iHit] == 344201220) {  // M3124 on the left
                    // if inside dead zone :
                    float xlim1=0.   +2.24462e-02 +0.1;  // x2, y2 val + misaligned + margin
                    float xlim2=-0.8 +2.24462e-02 +0.1;
                    float ylim1=0.   +8.72844e-02 -0.1;
                    float ylim2=1.5  +8.72844e-02 +0.1;
                    if ( (lptemp.X() < xlim1 && ylim1 < lptemp.Y() && lptemp.Y() <ylim2)  ||
                            (lptemp.X() < xlim2 && lptemp.Y() < ylim1) ) {
                        DQM_TrackPull5_X_M3124->Fill(lptemp.X() - theLP[iHit].X());
                        DQM_TrackPull5_Y_M3124->Fill(lptemp.Y() - theLP[iHit].Y());
                        DQM_dead_M3124->Fill(lptemp.X(), lptemp.Y());
                    } else {
                        // if in active zone:
                        DQM_TrackPull4_X_M3124->Fill(lptemp.X() - theLP[iHit].X());
                        DQM_TrackPull4_Y_M3124->Fill(lptemp.Y() - theLP[iHit].Y());
                    }
                    if (lptemp.Y() - theLP[iHit].Y() <-0.002) {
                        DQM_2npeak_M3124->Fill(lptemp.X(), lptemp.Y());
                    } else {
                        DQM_1speak_M3124->Fill(lptemp.X(), lptemp.Y());
                    }
                } else if (modulesID[iHit] == 344987652) { // M3074 on the left
                    // if inside dead zone :
                    float xlim1=0.   -9.50454e-02 +0.1;
                    float xlim2=-0.8 -9.50454e-02 -0.1;
                    float ylim1=-2.5 -9.50454e-02 +0.1;
                    float ylim2=0.   +1.33751e-01 -0.1;
                    float ylim3=0.8  +1.33751e-01 +0.1;
                    float ylim4=1.6  +1.33751e-01 -0.1;
                    float ylim5=2.4  +1.33751e-01 +0.1;
                    if ( (lptemp.X() < xlim1 && lptemp.Y() <ylim1) ||
                            (xlim2 < lptemp.X() && lptemp.X() < xlim1 && ylim2< lptemp.Y() &&  lptemp.Y() <ylim3) ||
                            (xlim2 < lptemp.X() && lptemp.X() < xlim1 && ylim4< lptemp.Y() &&  lptemp.Y() <ylim5) ) {
                        DQM_TrackPull5_X_M3074->Fill(lptemp.X() - theLP[iHit].X());
                        DQM_TrackPull5_Y_M3074->Fill(lptemp.Y() - theLP[iHit].Y());
                        DQM_dead_M3074->Fill(lptemp.X(), lptemp.Y());
                    }
                    // if in active zone:
                    else {
                        DQM_TrackPull4_X_M3074->Fill(lptemp.X() - theLP[iHit].X());
                        DQM_TrackPull4_Y_M3074->Fill(lptemp.Y() - theLP[iHit].Y());
                    }
                    if (lptemp.Y() - theLP[iHit].Y() >0.002) {
                        DQM_2npeak_M3074->Fill(lptemp.X(), lptemp.Y());
                    }
                    else {
                        DQM_1speak_M3074->Fill(lptemp.X(), lptemp.Y());
                    }

                }
                // end of test

                DQM_Corr_dX_X[modulesID[iHit]]->Fill(theLP[iHit].X(),lptemp.X() - theLP[iHit].X());
                DQM_Corr_dY_Y[modulesID[iHit]]->Fill(theLP[iHit].Y(),lptemp.Y() - theLP[iHit].Y());
            }
            else {
                //std::cout << " not found out on which layer is GP " << ilayernumber << std::endl;
                DQM_TrackPull_X[modulesID[iHit]]->Fill(-10 );
                DQM_TrackPull_Y[modulesID[iHit]]->Fill(-10 );
                DQM_TrackPull2_X[modulesID[iHit]]->Fill(-10 );
                DQM_TrackPull2_Y[modulesID[iHit]]->Fill(-10 );
            }

            //must add number of tracks per event
        }
    }

    // compute the residuals without including the cluster in the track GP list.
    for (unsigned int il1=0; il1<8; il1++) {
        int idetR = LayersDefinition[il1].first;
        ComputeResiduals(idetR, myEvId.event());
        int idetL = LayersDefinition[il1].second;
        ComputeResiduals(idetL, myEvId.event());
    }

    DQM_NumbOfTracks_per_Event->Fill(icount_trk);
    if (iclu_count>=4) DQM_NumbOfTracks_4cl_per_Event-> Fill(icount_trk) ;

    theTeleTrackCollection.clear();

}//end analyze()

// ------------ method called once each job just before starting event loop  ------------
void SimpleTrackingClean::beginJob ( ) {
    //std::cout << "ST14: beginJob" << std::endl;


    // load the noisy channel list
    std::ifstream file_noisy;
    //file_noisy.open("/opt/sbg/data/safe1/cms/ccollard/TrackerTelescope/ui6/PixelGeomV2/CMSSW_10_1_11/src/Geometry/PixelTelescope/test/noisy_list.txt");
    file_noisy.open("/eos/cms/store/group/dpg_tracker_upgrade/BeamTestTelescope/2018-TELESCOPE-COMMISSIONING/noisy_list.txt");
    inoisy=0;
    if (!file_noisy) {
        std::cerr << "cannot open file  noisy_list.txt " << std::endl;
    }
    else
    {
        while (!file_noisy.eof () && inoisy<1000) {
            file_noisy >> noisy_det[inoisy] >> noisy_barx[inoisy] >> noisy_bary[inoisy];
            inoisy++;
        }
    }
    std::cout << " file_noisy : " << inoisy << " entries " << std::endl;
    if (inoisy>0) std::cout << " example detId "<< noisy_det[0] << " x "<< noisy_barx[0] << " y " << noisy_bary[0] <<  std::endl;
    file_noisy.close ();
    //
    //
}//end void beginJob ( )

// ------------ method called once each job just after ending the event loop  ------------
void SimpleTrackingClean::endJob ( ) {
    //std::cout << "ST15: endJob" << std::endl;
}//end void endJobd ( )

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SimpleTrackingClean::fillDescriptions ( edm::ConfigurationDescriptions& descriptions ) {

    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);

}//end void fillDescriptions ( )

void SimpleTrackingClean::ApplyXmove(int detId,  double &xaligned, double &yaligned) {
    // example of manual misalignment shift with parameter in the configuration file
    if (detId==353375236) {
        xaligned+=xmove_par;
    }
}

void SimpleTrackingClean::ApplyAlignment(int detId, double &xaligned, double &yaligned) {
    // apply coarse alignement
    if (detId==353114116) {
        xaligned-=-2.56917e-02;
        yaligned-=-2.38370e-02;
    }
    else if (detId==353113092) {
        xaligned-=-2.02377e-02;
        yaligned-=-2.70933e-02;
    }
    if (detId==352851972) {
        xaligned-=2.47805e-02;
        yaligned-=-6.97738e-03;
    }
    else if (detId==352850948) {
        xaligned-=1.85113e-02;
        yaligned-=-9.17572e-03;
    }
    else if (detId==352589828) {
        xaligned-=5.50173e-02;
        yaligned-=-3.05621e-02;
    }
    else if (detId==352588804) {
        xaligned-=6.73395e-02;
        yaligned-=-3.06608e-02;
    }
    else if (detId==344200196) {
        xaligned-=-2.77220e-02;
        yaligned-=-9.31173e-02;
    }
    else if (detId==344201220) {
        xaligned-=-2.24462e-02;
        yaligned-=-8.72844e-02;
    }
    else if (detId==344462340) {
        xaligned-=7.14737e-03;
        yaligned-=-1.09290e-01;
    }
    else if (detId==344724484) {
        xaligned-=-7.60673e-02;
        yaligned-=-1.20672e-01;
    }
    else if (detId==344986628) {
        xaligned-=7.71767e-02;
        yaligned-=-1.20759e-01;
    }
    else if (detId==344987652) {
        xaligned-=9.50454e-02;
        yaligned-=-1.33751e-01;
    }

}

//******************************************************
// SimpleTrackingClean::doSeeding
// pair of clusters
// from the 2 first layers,
// with delta(X) and delta(Y) compatible within 0.1 cm
// computed in the local positions
//******************************************************

void SimpleTrackingClean::doSeeding(edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelclusters,
                               const TrackerGeometry *tkgeom,  const PixelClusterParameterEstimator &cpe, edm::ESHandle<TrackerGeometry> tracker) {

    //std::cout << "-------------" << std::endl;
    for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end(); DSViter++   ) {

        edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
        edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();

        auto detid = DetId(DSViter->detId());
        //std::cout << "detID "  << int(detid) << std::endl;

        /*
            std::cout << " caro debug position " << std::endl;

               for ( edmNew::DetSet<SiPixelCluster>::const_iterator itCluster=begin; itCluster!=end; ++itCluster ) {
                 const PixelGeomDetUnit *pixdet = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(detid);
                LocalPoint lp(-9999., -9999., -9999.);

                PixelClusterParameterEstimator::ReturnType params = cpe.getParameters(*itCluster,*pixdet);
                lp = std::get<0>(params);

                double xl=0, yl=0, zl=0;
                xl = lp.x();
                yl = lp.y();
                zl = lp.z();

                const Surface& surface = tracker->idToDet(detid)->surface();
                GlobalPoint gp = surface.toGlobal(lp);

                double x=0, y=0, z=0;
                x = gp.x();
                y = gp.y();
                z = gp.z();
                std::cout << " caro local position    x   "  << xl  << "   y "  << yl  <<  "   z   "  << zl  << std::endl;
                std::cout << " caro global position   x   "  << x  << "   y "  << y  <<  "   z   "  << z  << std::endl;
                std::cout << " caro barycenter position   x   "  << itCluster->x()  << "   y "  << itCluster->y() << std::endl;

                TVector3 gp2(x, y, z);
                TVector3 lp2 = theSimpleLayers[0].getLocalPointPosition(gp2);
                double xl2=0, yl2=0, zl2=0;
                xl2 = lp2.x();
                yl2 = lp2.y();
                zl2 = lp2.z();
                std::cout << " caro global position*rot    x   "  << xl2  << "   y "  << yl2  <<  "   z   "  << zl2  << std::endl;

                TVector3 lp3(xl, yl, zl);
                TVector3 gp3 = theSimpleLayers[0].getLocalPointPosition(lp3);
                double x3=0, y3=0, z3=0;
                x3 = gp3.x();
                y3 = gp3.y();
                z3 = gp3.z();
                std::cout << " caro local position*rot    x   "  << x3  << "   y "  << y3  <<  "   z   "  << z3  << std::endl;


               }
               std::cout << " end caro debug position " << std::endl;
        */

        // first layer called layer0
        if(int(detid) != 353376260 && int(detid) != 353375236) continue;

        // Then loop on the clusters of the module
        for ( edmNew::DetSet<SiPixelCluster>::const_iterator itCluster=begin; itCluster!=end; ++itCluster ) {

            const PixelGeomDetUnit *pixdet = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(detid);
            LocalPoint lp(-9999., -9999., -9999.);

            PixelClusterParameterEstimator::ReturnType params = cpe.getParameters(*itCluster,*pixdet);
            lp = std::get<0>(params);

            const Surface& surface = tracker->idToDet(detid)->surface();
            GlobalPoint gp = surface.toGlobal(lp);

            double x=0, y=0, z=0;
            x = gp.x();
            y = gp.y();
            z = gp.z();

            // check the noisy clusters
            bool noisyflag=false;
            for (int in=0; in<inoisy; in++) {
                if (int(detid)==noisy_det[in]) {
                    double xbary=itCluster->x();
                    double ybary=itCluster->y();
                    if (int(xbary)==noisy_barx[in] && int(ybary)==noisy_bary[in]) {
                        noisyflag=true;
                    }
                }
            }
            if (noisyflag) continue;
            // end check

            /*
            	if ((x < -0.3) && (x > 0.1)) continue;
            	if ((y < -1.1) && (y > 0.0)) continue;
            */

            for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter2=pixelclusters->begin(); DSViter2!=pixelclusters->end(); DSViter2++   ) {

                auto detid2 = DetId(DSViter2->detId());

                if(int(detid2) == 353376260 || int(detid2) == 353375236) continue;   // reject layer0


                if(int(detid2) !=353114116  && int(detid2) != 353113092) continue;  // second layer, called layer1

                // idea for next step: try to consider also the combination layer 0 + layer 2 or layer 1 + layer 2
                // but in that case, we probably need a protection for the seeding not to get duplication of seed tracks.

                edmNew::DetSet<SiPixelCluster>::const_iterator begin2=DSViter2->begin();
                edmNew::DetSet<SiPixelCluster>::const_iterator end2  =DSViter2->end();

                for ( edmNew::DetSet<SiPixelCluster>::const_iterator itCluster2=begin2; itCluster2!=end2; ++itCluster2 ) {

                    const PixelGeomDetUnit *pixdet2 = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(detid2);
                    LocalPoint lp2(-9999., -9999., -9999.);

                    PixelClusterParameterEstimator::ReturnType params2 = cpe.getParameters(*itCluster2,*pixdet2);
                    lp2 = std::get<0>(params2);

                    const Surface& surface2 = tracker->idToDet(detid2)->surface();
                    GlobalPoint gp2 = surface2.toGlobal(lp2);

                    double x2=0, y2=0, z2=0;
                    x2 = gp2.x();
                    y2 = gp2.y();
                    z2 = gp2.z();

                    // check the noisy clusters
                    bool noisyflag2=false;
                    for (int in=0; in<inoisy; in++) {
                        if (int(detid2)==noisy_det[in]) {
                            double xbary2=itCluster2->x();
                            double ybary2=itCluster2->y();
                            if (int(xbary2)==noisy_barx[in] && int(ybary2)==noisy_bary[in]) {
                                noisyflag2=true;
                            }
                        }
                    }
                    if (noisyflag2) continue;
                    // end check


                    ApplyAlignment(int(detid2),x2,y2);
                    //ApplyXmove(int(detid2),x2,y2);
                    
                    //std::cout << "detid  " << int(detid)  << "  x  " << x << " y " << y << " z " << z << std::endl;
                    //std::cout << "detid2 " << int(detid2) << " x2 " << x2 << " y2 " << y2 << " z2 "<< z2 << std::endl;

                    if(fabs(x-x2) < 0.1 && fabs(y-y2) < 0.1 && fabs(z-z2) < 1000 ) {

                        TelescopeTracks theseed;
                        TVector3 clust1(x,   y,  z);
                        TVector3 clust2(x2, y2, z2);
                        //TVector3 clust_err(0.05, 0.05, 0.05);
                        TVector3 clust_err(0.0029, 0.0043, 0.0029);  // to be checked 150microns/sqrt(12), 100microns/sqrt(12);

                        theseed.addGlobalPoint(clust1);
                        theseed.addGlobalPointErr(clust_err);
                        theseed.addCluster(*itCluster);
                        theseed.addModDetId(int(detid));

                        // for local position, transform the cmssw point to local axes obtained with the SimplePlane definition
                        // as using x2/y2 after misalignement is not correct because local position should remain invariant
                        // and using cmssw coordinate is to complex to deal with (test on 24-25/01/2019 -> total mess at the end)
                        double xlocclust1=lp.x();
                        if (x-0.03>50-1.1*(y+0.03)/5.6 && x+0.03<51.6-1.1*(y-0.03)/5.6) { // some margin in case of misalignment
                            xlocclust1+=0.81;
                        }
                        else {
                            xlocclust1+=0.81;
                            xlocclust1*=-1.;
                        }

                        TVector3 locclust1(xlocclust1,lp.y(),lp.z());
                        if (int(detid) == 353376260) {
                            //locclust1 = theSimpleLayersR[0].getLocalPointPosition(clust1);
                            theseed.addPseudoLocalPoint(locclust1);
                            theseed.addPseudoLocalPointErr(clust_err);
                            theseed.addAssPlaneA(theSimpleLayersR[0].getParamA());
                            theseed.addAssPlaneB(theSimpleLayersR[0].getParamB());
                            theseed.addAssPlaneC(theSimpleLayersR[0].getParamC());
                            theseed.addAssPlaneD(theSimpleLayersR[0].getParamD());
                            theseed.addAssPlaneX0(theSimpleLayersR[0].getOriginX());
                            theseed.addAssPlaneY0(theSimpleLayersR[0].getOriginY());
                            theseed.addAssPlaneZ0(theSimpleLayersR[0].getOriginZ());
                            theseed.addAssPlaneTheta(theSimpleLayersR[0].getTheta());
                            theseed.addAssPlanePhi(theSimpleLayersR[0].getPhi());
                        }
                        else if (int(detid) == 353375236) {
                            //locclust1 = theSimpleLayersL[0].getLocalPointPosition(clust1);
                            theseed.addPseudoLocalPoint(locclust1);
                            theseed.addPseudoLocalPointErr(clust_err);
                            theseed.addAssPlaneA(theSimpleLayersL[0].getParamA());
                            theseed.addAssPlaneB(theSimpleLayersL[0].getParamB());
                            theseed.addAssPlaneC(theSimpleLayersL[0].getParamC());
                            theseed.addAssPlaneD(theSimpleLayersL[0].getParamD());
                            theseed.addAssPlaneX0(theSimpleLayersL[0].getOriginX());
                            theseed.addAssPlaneY0(theSimpleLayersL[0].getOriginY());
                            theseed.addAssPlaneZ0(theSimpleLayersL[0].getOriginZ());
                            theseed.addAssPlaneTheta(theSimpleLayersL[0].getTheta());
                            theseed.addAssPlanePhi(theSimpleLayersL[0].getPhi());
                        }


                        theseed.addGlobalPoint(clust2);
                        theseed.addGlobalPointErr(clust_err);
                        theseed.addCluster(*itCluster2);
                        theseed.addModDetId(int(detid2));

                        double xlocclust2=lp2.x();
                        if (x2-0.03>50-1.1*(y2+0.03)/5.6 && x2+0.03<51.6-1.1*(y2-0.03)/5.6) { // some margin in case of misalignment
                            xlocclust2+=0.81;
                        }
                        else {
                            xlocclust2+=0.81;
                            xlocclust2*=-1.;
                        }
                        TVector3 locclust2(xlocclust2,lp2.y(),lp2.z());
                        if(int(detid2) ==353114116) {
                            //locclust2 = theSimpleLayersR[1].getLocalPointPosition(clust2);
                            theseed.addPseudoLocalPoint(locclust2);
                            theseed.addPseudoLocalPointErr(clust_err);
                            theseed.addAssPlaneA(theSimpleLayersR[1].getParamA());
                            theseed.addAssPlaneB(theSimpleLayersR[1].getParamB());
                            theseed.addAssPlaneC(theSimpleLayersR[1].getParamC());
                            theseed.addAssPlaneD(theSimpleLayersR[1].getParamD());
                            theseed.addAssPlaneX0(theSimpleLayersR[1].getOriginX());
                            theseed.addAssPlaneY0(theSimpleLayersR[1].getOriginY());
                            theseed.addAssPlaneZ0(theSimpleLayersR[1].getOriginZ());
                            theseed.addAssPlaneTheta(theSimpleLayersR[1].getTheta());
                            theseed.addAssPlanePhi(theSimpleLayersR[1].getPhi());
                        }
                        if(int(detid2) == 353113092) {
                            //locclust2 = theSimpleLayersL[1].getLocalPointPosition(clust2);
                            theseed.addPseudoLocalPoint(locclust2);
                            theseed.addPseudoLocalPointErr(clust_err);
                            theseed.addAssPlaneA(theSimpleLayersL[1].getParamA());
                            theseed.addAssPlaneB(theSimpleLayersL[1].getParamB());
                            theseed.addAssPlaneC(theSimpleLayersL[1].getParamC());
                            theseed.addAssPlaneD(theSimpleLayersL[1].getParamD());
                            theseed.addAssPlaneX0(theSimpleLayersL[1].getOriginX());
                            theseed.addAssPlaneY0(theSimpleLayersL[1].getOriginY());
                            theseed.addAssPlaneZ0(theSimpleLayersL[1].getOriginZ());
                            theseed.addAssPlaneTheta(theSimpleLayersL[1].getTheta());
                            theseed.addAssPlanePhi(theSimpleLayersL[1].getPhi());
                        }

                        //std::cout << "in fit tracks" << std::endl;
                        theseed.fitTrack();

                        theTeleTrackCollection.push_back(theseed);

                        double xtemp, ytemp, ztemp;
                        double parFit[6] = {theseed.getParameter(0), theseed.getParameter(1), theseed.getParameter(2), theseed.getParameter(3), theseed.getParameter(4), theseed.getParameter(5)};
                        //theseed.line(z, parFit, xtemp, ytemp, ztemp);
                        double planeq0[4];
                        if (int(detid) == 353376260) {
                            planeq0[0]= theSimpleLayersR[0].getParamA();
                            planeq0[1]= theSimpleLayersR[0].getParamB();
                            planeq0[2]= theSimpleLayersR[0].getParamC();
                            planeq0[3]= theSimpleLayersR[0].getParamD();
                        }
                        else if (int(detid) == 353375236) {
                            planeq0[0]= theSimpleLayersL[0].getParamA();
                            planeq0[1]= theSimpleLayersL[0].getParamB();
                            planeq0[2]= theSimpleLayersL[0].getParamC();
                            planeq0[3]= theSimpleLayersL[0].getParamD();
                        }
                        theseed.intersection(planeq0, parFit, xtemp, ytemp, ztemp);
                        //std::cout << "  check the fit " << xtemp << " " << ytemp << "  " << ztemp << std::endl;
                        //std::cout << "  to be compared with " << x  << " " << y  << "  " << z  << std::endl;
                        //std::cout << "-------------" << std::endl;
                        //theseed.line(z2, parFit, xtemp, ytemp, ztemp);
                        double planeq1[4];
                        if(int(detid2) ==353114116) {
                            planeq1[0]= theSimpleLayersR[1].getParamA();
                            planeq1[1]= theSimpleLayersR[1].getParamB();
                            planeq1[2]= theSimpleLayersR[1].getParamC();
                            planeq1[3]= theSimpleLayersR[1].getParamD();
                        }
                        if(int(detid2) == 353113092) {
                            planeq1[0]= theSimpleLayersL[1].getParamA();
                            planeq1[1]= theSimpleLayersL[1].getParamB();
                            planeq1[2]= theSimpleLayersL[1].getParamC();
                            planeq1[3]= theSimpleLayersL[1].getParamD();
                        }
                        theseed.intersection(planeq1, parFit, xtemp, ytemp, ztemp);
                        //std::cout << "  check the fit 2 " << xtemp << " " << ytemp << "  " << ztemp << std::endl;
                        //std::cout << "  to be compared with " << x2  << " " << y2  << "  " << z2  << std::endl; // ??
                    }
                }
            }
        }
    }
}

//*************************************************************
// SimpleTrackingClean::doPatternReco
// loop on each layers to add clusters to the tracks :
// extrapolate the tracks in the considered layer,
// compare the distance in the 2D plane of the telescope
// between the extrapolation of the tracks and the cluster,
// add the closest cluster satisfying DeltaX<0.3 and DeltaY<0.3
//**************************************************************

void SimpleTrackingClean::doPatternReco(edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelclusters,
                                   const TrackerGeometry *tkgeom,  const PixelClusterParameterEstimator &cpe, edm::ESHandle<TrackerGeometry> tracker) {

    for(unsigned itrack=0; itrack<theTeleTrackCollection.size(); itrack++) {
        double distanceR = 0.;
        double deltaX = 0.;
        double deltaY = 0.;

        for(int ilayer = 2; ilayer <8; ilayer++) {

            //if (ilayer==6) continue;     // just noise in M3009 on the right side, and M3057 is dead module on the left side

            for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end(); DSViter++   ) {

                edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
                edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();

                auto detid = DetId(DSViter->detId());


                if( int(detid) != LayersDefinition[ilayer].first && int(detid) != LayersDefinition[ilayer].second) continue;

                bool isleft=false;
                bool isright=false;
                if (int(detid) == LayersDefinition[ilayer].first)  {
                    isright=true;
                }
                else if (int(detid) == LayersDefinition[ilayer].second) {
                    isleft=true;
                }

                TVector3 closestClust;
                TVector3 closestClust2;
                SiPixelCluster closustSiPixelCulster;
                int closestR = 1000;
                TVector3 clust_err(0.0029, 0.0043, 0.0029);  // to be checked 150microns/sqrt(12), 100microns/sqrt(12);

                // Then loop on the clusters of the module
                for ( edmNew::DetSet<SiPixelCluster>::const_iterator itCluster=begin; itCluster!=end; ++itCluster ) {

                    const PixelGeomDetUnit *pixdet = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(detid);
                    LocalPoint lp(-9999., -9999., -9999.);

                    PixelClusterParameterEstimator::ReturnType params = cpe.getParameters(*itCluster,*pixdet);
                    lp = std::get<0>(params);

                    const Surface& surface = tracker->idToDet(detid)->surface();
                    GlobalPoint gp = surface.toGlobal(lp);

                    double x3=0, y3=0, z3=0;
                    x3 = gp.x();
                    y3 = gp.y();
                    z3 = gp.z();

                    // check the noisy clusters
                    bool noisyflag=false;
                    for (int in=0; in<inoisy; in++) {
                        if (int(detid)==noisy_det[in]) {
                            double xbary=itCluster->x();
                            double ybary=itCluster->y();
                            if (int(xbary)==noisy_barx[in] && int(ybary)==noisy_bary[in]) {
                                noisyflag=true;
                            }
                        }
                    }
                    if (noisyflag) continue;
                    // end check


                    ApplyAlignment(int(detid),x3,y3);
                    //ApplyXmove(int(detid),x3,y3);

                    TVector3 clust(x3,   y3,  z3);
                    double xlocclust3=lp.x();
                    if (x3-0.03>50-1.1*(y3+0.03)/5.6 && x3+0.03<51.6-1.1*(y3-0.03)/5.6) { // some margin in case of misalignment
                        xlocclust3+=0.81;
                    }
                    else {
                        xlocclust3+=0.81;
                        xlocclust3*=-1.;
                    }

                    TVector3 locclust(xlocclust3,lp.y(),lp.z());
                    /*
                              TVector3 locclust;
                              if (isright) {
                               locclust= theSimpleLayersR[ilayer].getLocalPointPosition(clust);
                              }
                              else if (isleft)  {
                               locclust= theSimpleLayersL[ilayer].getLocalPointPosition(clust);
                              }
                    */

                    // compute the extrapolation of the track in the considered plane
                    double xtemp, ytemp, ztemp;
                    double parFit[6] = {theTeleTrackCollection[itrack].getParameter(0), theTeleTrackCollection[itrack].getParameter(1), theTeleTrackCollection[itrack].getParameter(2), theTeleTrackCollection[itrack].getParameter(3), theTeleTrackCollection[itrack].getParameter(4), theTeleTrackCollection[itrack].getParameter(5)};
                    //theTeleTrackCollection[itrack].line(z3, parFit, xtemp, ytemp, ztemp);
                    double planeq[4];
                    if (isright) {
                        planeq[0]= theSimpleLayersR[ilayer].getParamA();
                        planeq[1]= theSimpleLayersR[ilayer].getParamB();
                        planeq[2]= theSimpleLayersR[ilayer].getParamC();
                        planeq[3]= theSimpleLayersR[ilayer].getParamD();
                    }
                    else if (isleft)  {
                        planeq[0]= theSimpleLayersL[ilayer].getParamA();
                        planeq[1]= theSimpleLayersL[ilayer].getParamB();
                        planeq[2]= theSimpleLayersL[ilayer].getParamC();
                        planeq[3]= theSimpleLayersL[ilayer].getParamD();
                    }
                    theTeleTrackCollection[itrack].intersection(planeq, parFit, xtemp, ytemp, ztemp);
                    // transform the global position into local position
                    TVector3 gptemp(xtemp, ytemp, ztemp);
                    TVector3 pseudoloc;
                    if (isright) {
                        pseudoloc = theSimpleLayersR[ilayer].getLocalPointPosition(gptemp);
                    }
                    else if (isleft)  {
                        pseudoloc = theSimpleLayersL[ilayer].getLocalPointPosition(gptemp);
                    }


                    // compute the distance in the 2D plane of the telescope
                    distanceR = pow( pow(pseudoloc.X() - locclust.X() , 2)+ pow(pseudoloc.Y() - locclust.Y(), 2),  0.5);
                    deltaX = pseudoloc.X() - locclust.X();
                    deltaY = pseudoloc.Y() - locclust.Y();

                    if(fabs(deltaX) < 0.3 &&  fabs(deltaY) < 0.3 && closestR > distanceR) {
                        // limit the region where to find a close cluster
                        closestR=distanceR;
                        closestClust.SetX(x3);
                        closestClust.SetY(y3);
                        closestClust.SetZ(z3);
                        closestClust2 = locclust;
                        closustSiPixelCulster=*itCluster;
                    }
                }  // end loop on cluster

                if (closestR<1000) {
                    // don't add the point (0,0,0) if no close cluster has been found
                    theTeleTrackCollection[itrack].addGlobalPoint(closestClust);
                    theTeleTrackCollection[itrack].addGlobalPointErr(clust_err);
                    theTeleTrackCollection[itrack].addCluster(closustSiPixelCulster);
                    theTeleTrackCollection[itrack].addModDetId(int(detid));
                    theTeleTrackCollection[itrack].addPseudoLocalPoint(closestClust2);
                    theTeleTrackCollection[itrack].addPseudoLocalPointErr(clust_err);
                    if (isright) {
                        theTeleTrackCollection[itrack].addAssPlaneA(theSimpleLayersR[ilayer].getParamA());
                        theTeleTrackCollection[itrack].addAssPlaneB(theSimpleLayersR[ilayer].getParamB());
                        theTeleTrackCollection[itrack].addAssPlaneC(theSimpleLayersR[ilayer].getParamC());
                        theTeleTrackCollection[itrack].addAssPlaneD(theSimpleLayersR[ilayer].getParamD());
                        theTeleTrackCollection[itrack].addAssPlaneX0(theSimpleLayersR[ilayer].getOriginX());
                        theTeleTrackCollection[itrack].addAssPlaneY0(theSimpleLayersR[ilayer].getOriginY());
                        theTeleTrackCollection[itrack].addAssPlaneZ0(theSimpleLayersR[ilayer].getOriginZ());
                        theTeleTrackCollection[itrack].addAssPlaneTheta(theSimpleLayersR[ilayer].getTheta());
                        theTeleTrackCollection[itrack].addAssPlanePhi(theSimpleLayersR[ilayer].getPhi());
                    }
                    else if (isleft) {
                        theTeleTrackCollection[itrack].addAssPlaneA(theSimpleLayersL[ilayer].getParamA());
                        theTeleTrackCollection[itrack].addAssPlaneB(theSimpleLayersL[ilayer].getParamB());
                        theTeleTrackCollection[itrack].addAssPlaneC(theSimpleLayersL[ilayer].getParamC());
                        theTeleTrackCollection[itrack].addAssPlaneD(theSimpleLayersL[ilayer].getParamD());
                        theTeleTrackCollection[itrack].addAssPlaneX0(theSimpleLayersL[ilayer].getOriginX());
                        theTeleTrackCollection[itrack].addAssPlaneY0(theSimpleLayersL[ilayer].getOriginY());
                        theTeleTrackCollection[itrack].addAssPlaneZ0(theSimpleLayersL[ilayer].getOriginZ());
                        theTeleTrackCollection[itrack].addAssPlaneTheta(theSimpleLayersL[ilayer].getTheta());
                        theTeleTrackCollection[itrack].addAssPlanePhi(theSimpleLayersL[ilayer].getPhi());
                    }
                    theTeleTrackCollection[itrack].fitTrack();
                }
            }// end loop on DSV
        }// end loop on Layer
        distanceR_per_Track-> Fill (distanceR) ;
        deltaX_per_Track-> Fill (deltaX) ;
        deltaY_per_Track-> Fill (deltaY) ;
    }// end loop on TrackCollection
}

//*************************************************************
// SimpleTrackingClean::ComputeResiduals
// compute the residual for a point after having removed it
// from the track extrapolation
//*************************************************************

void SimpleTrackingClean::ComputeResiduals(int detId_to_study, int eventNumber) {
// goal of this function: compute the residual for a point after having removed it from the track extrapolation
// for the purpose of the alignment, also store the biased residual in the root file

    // create a new track without the point from detId_to_study included
    for(unsigned int itrack = 0; itrack< theTeleTrackCollection.size(); itrack++) {
        std::vector<int > modulesID = theTeleTrackCollection[itrack].getModDetId();
        std::vector<TVector3 > theGP = theTeleTrackCollection[itrack].getGlobalPoints();
        std::vector<TVector3 > theLP = theTeleTrackCollection[itrack].getPseudoLocalPoints();
        std::vector<double > plane_paramA = theTeleTrackCollection[itrack].getAssPlaneA();
        std::vector<double > plane_paramB = theTeleTrackCollection[itrack].getAssPlaneB();
        std::vector<double > plane_paramC = theTeleTrackCollection[itrack].getAssPlaneC();
        std::vector<double > plane_paramD = theTeleTrackCollection[itrack].getAssPlaneD();
        std::vector<double > plane_x0 = theTeleTrackCollection[itrack].getAssPlaneX0();
        std::vector<double > plane_y0 = theTeleTrackCollection[itrack].getAssPlaneY0();
        std::vector<double > plane_z0 = theTeleTrackCollection[itrack].getAssPlaneZ0();
        std::vector<double > plane_th = theTeleTrackCollection[itrack].getAssPlaneTheta();
        std::vector<double > plane_ph = theTeleTrackCollection[itrack].getAssPlanePhi();
        std::vector<SiPixelCluster> clu_list = theTeleTrackCollection[itrack].getclusterList();
        std::vector<TVector3> glob_pos = theTeleTrackCollection[itrack].getGlobalPoints();
        std::vector<TVector3> glob_er = theTeleTrackCollection[itrack].getGlobalPointsErr();
        std::vector<TVector3> pseudo_pos = theTeleTrackCollection[itrack].getPseudoLocalPoints();
        std::vector<TVector3> pseudo_er = theTeleTrackCollection[itrack].getPseudoLocalPointsErr();
        TVector3 theGP_to_study;
        TVector3 theLP_to_study;
        TelescopeTracks theModTeleTrack;
        bool found_cluster_to_study=false;
        int counter_theModTeleTrack=0;
        double planeq0[4];
        double parFitbiased[6] = {theTeleTrackCollection[itrack].getParameter(0), theTeleTrackCollection[itrack].getParameter(1), theTeleTrackCollection[itrack].getParameter(2), theTeleTrackCollection[itrack].getParameter(3), theTeleTrackCollection[itrack].getParameter(4), theTeleTrackCollection[itrack].getParameter(5)};
        for(unsigned int iHit=0; iHit < modulesID.size(); iHit++) {
            if (detId_to_study==modulesID[iHit]) {
                // the module to study
                theGP_to_study=theGP[iHit];
                theLP_to_study=theLP[iHit];
                planeq0[0]= plane_paramA[iHit];
                planeq0[1]= plane_paramB[iHit];
                planeq0[2]= plane_paramC[iHit];
                planeq0[3]= plane_paramD[iHit];
                found_cluster_to_study=true;

            }
            else {
                // to be included in the theModTeleTrackCollection
                theModTeleTrack.addGlobalPoint(glob_pos[iHit]);
                theModTeleTrack.addGlobalPointErr(glob_er[iHit]);
                theModTeleTrack.addCluster(clu_list[iHit]);
                theModTeleTrack.addModDetId(modulesID[iHit]);
                theModTeleTrack.addPseudoLocalPoint(pseudo_pos[iHit]);
                theModTeleTrack.addPseudoLocalPointErr(pseudo_er[iHit]);
                theModTeleTrack.addAssPlaneA(plane_paramA[iHit]);
                theModTeleTrack.addAssPlaneB(plane_paramB[iHit]);
                theModTeleTrack.addAssPlaneC(plane_paramC[iHit]);
                theModTeleTrack.addAssPlaneD(plane_paramD[iHit]);
                theModTeleTrack.addAssPlaneX0(plane_x0[iHit]);
                theModTeleTrack.addAssPlaneY0(plane_y0[iHit]);
                theModTeleTrack.addAssPlaneZ0(plane_z0[iHit]);
                theModTeleTrack.addAssPlaneTheta(plane_th[iHit]);
                theModTeleTrack.addAssPlanePhi(plane_ph[iHit]);
                counter_theModTeleTrack++;
            }
        }
        if (!found_cluster_to_study) continue;
        if (counter_theModTeleTrack<3) continue;

        int ilayernumber=-1;
        bool isleft=false;
        bool isright=false;
        for (int il1=0; il1<8; il1++) {
            if (detId_to_study == LayersDefinition[il1].first)  {
                ilayernumber=il1;
                isright=true;
            }
            else if (detId_to_study == LayersDefinition[il1].second) {
                ilayernumber=il1;
                isleft=true;
            }
        }

        // biased residual
        double xbiased, ybiased, zbiased;
        theTeleTrackCollection[itrack].intersection(planeq0, parFitbiased, xbiased, ybiased, zbiased);
        TVector3 gpbiased(xbiased, ybiased, zbiased);
        TVector3 lpbiased;
        if (isright) lpbiased= theSimpleLayersR[ilayernumber].getLocalPointPosition(gpbiased);
        else if (isleft) lpbiased= theSimpleLayersL[ilayernumber].getLocalPointPosition(gpbiased);

        // fit of the modified track
        theModTeleTrack.fitTrack();
        // compute of the unbiased residual for theGP_to_study
        double xtemp, ytemp, ztemp;
        double parFit[6] = {theModTeleTrack.getParameter(0), theModTeleTrack.getParameter(1), theModTeleTrack.getParameter(2), theModTeleTrack.getParameter(3), theModTeleTrack.getParameter(4), theModTeleTrack.getParameter(5)};
        theModTeleTrack.intersection(planeq0, parFit, xtemp, ytemp, ztemp);
        TVector3 gptemp(xtemp, ytemp, ztemp);

        DQM2_Chi2[detId_to_study]->Fill(theModTeleTrack.getChi2());
        // xtemp = global coordinates
        // --> transformation into LP
        TVector3 lptemp;
        if (isright) lptemp= theSimpleLayersR[ilayernumber].getLocalPointPosition(gptemp);
        else if (isleft) lptemp= theSimpleLayersL[ilayernumber].getLocalPointPosition(gptemp);

        // store the info
        DQM2_TrackPull_X[detId_to_study]->Fill(lptemp.X() - theLP_to_study.X());
        DQM2_TrackPull_Y[detId_to_study]->Fill(lptemp.Y() - theLP_to_study.Y());

        // Fill the tree
        tree_eventN = eventNumber;
        tree_detId2 = detId_to_study;
        tree_biasedResidualX = lpbiased.X() - theLP_to_study.X() ;
        tree_biasedResidualY = lpbiased.Y() - theLP_to_study.Y() ;
        tree_unbiasedResidualX = lptemp.X() - theLP_to_study.X() ;
        tree_unbiasedResidualY = lptemp.Y() - theLP_to_study.Y() ;
        residualsTree->Fill();

        double Sum2study= (lptemp.X() - theLP_to_study.X())*(lptemp.X() - theLP_to_study.X())
                          + (lptemp.Y() - theLP_to_study.Y())*(lptemp.Y() - theLP_to_study.Y());
        DQM2_Sum2[detId_to_study]->Fill(Sum2study);
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE ( SimpleTrackingClean ) ;

/*
void SimpleTrackingClean::ComputeResiduals(int detId_to_study, int eventNumber, double brX, double brY) {
// goal of this function: compute the residual for a point after having removed it from the track extrapolation

    // create a new track without the point from detId_to_study included
    for(unsigned int itrack = 0; itrack< theTeleTrackCollection.size(); itrack++) {
        std::vector<int > modulesID = theTeleTrackCollection[itrack].getModDetId();
        std::vector<TVector3 > theGP = theTeleTrackCollection[itrack].getGlobalPoints();
        std::vector<TVector3 > theLP = theTeleTrackCollection[itrack].getPseudoLocalPoints();
        std::vector<double > plane_paramA = theTeleTrackCollection[itrack].getAssPlaneA();
        std::vector<double > plane_paramB = theTeleTrackCollection[itrack].getAssPlaneB();
        std::vector<double > plane_paramC = theTeleTrackCollection[itrack].getAssPlaneC();
        std::vector<double > plane_paramD = theTeleTrackCollection[itrack].getAssPlaneD();
        std::vector<double > plane_x0 = theTeleTrackCollection[itrack].getAssPlaneX0();
        std::vector<double > plane_y0 = theTeleTrackCollection[itrack].getAssPlaneY0();
        std::vector<double > plane_z0 = theTeleTrackCollection[itrack].getAssPlaneZ0();
        std::vector<double > plane_th = theTeleTrackCollection[itrack].getAssPlaneTheta();
        std::vector<double > plane_ph = theTeleTrackCollection[itrack].getAssPlanePhi();
        std::vector<SiPixelCluster> clu_list = theTeleTrackCollection[itrack].getclusterList();
        std::vector<TVector3> glob_pos = theTeleTrackCollection[itrack].getGlobalPoints();
        std::vector<TVector3> glob_er = theTeleTrackCollection[itrack].getGlobalPointsErr();
        std::vector<TVector3> pseudo_pos = theTeleTrackCollection[itrack].getPseudoLocalPoints();
        std::vector<TVector3> pseudo_er = theTeleTrackCollection[itrack].getPseudoLocalPointsErr();
        TVector3 theGP_to_study;
        TVector3 theLP_to_study;
        TelescopeTracks theModTeleTrack;
        bool found_cluster_to_study=false;
        int counter_theModTeleTrack=0;
        double planeq0[4];
        for(unsigned int iHit=0; iHit < modulesID.size(); iHit++) {
            if (detId_to_study==modulesID[iHit]) {
                // the module to study
                theGP_to_study=theGP[iHit];
                theLP_to_study=theLP[iHit];
                planeq0[0]= plane_paramA[iHit];
                planeq0[1]= plane_paramB[iHit];
                planeq0[2]= plane_paramC[iHit];
                planeq0[3]= plane_paramD[iHit];
                found_cluster_to_study=true;
            } else {
                // to be included in the theModTeleTrackCollection
                theModTeleTrack.addGlobalPoint(glob_pos[iHit]);
                theModTeleTrack.addGlobalPointErr(glob_er[iHit]);
                theModTeleTrack.addCluster(clu_list[iHit]);
                theModTeleTrack.addModDetId(modulesID[iHit]);
                theModTeleTrack.addPseudoLocalPoint(pseudo_pos[iHit]);
                theModTeleTrack.addPseudoLocalPointErr(pseudo_er[iHit]);
                theModTeleTrack.addAssPlaneA(plane_paramA[iHit]);
                theModTeleTrack.addAssPlaneB(plane_paramB[iHit]);
                theModTeleTrack.addAssPlaneC(plane_paramC[iHit]);
                theModTeleTrack.addAssPlaneD(plane_paramD[iHit]);
                theModTeleTrack.addAssPlaneX0(plane_x0[iHit]);
                theModTeleTrack.addAssPlaneY0(plane_y0[iHit]);
                theModTeleTrack.addAssPlaneZ0(plane_z0[iHit]);
                theModTeleTrack.addAssPlaneTheta(plane_th[iHit]);
                theModTeleTrack.addAssPlanePhi(plane_ph[iHit]);
                counter_theModTeleTrack++;
            }
        }
        if (!found_cluster_to_study) continue;
        if (counter_theModTeleTrack<3) continue;
        theModTeleTrack.fitTrack();

        // compute of the residual for theGP_to_study
        double xtemp, ytemp, ztemp;
        double parFit[6] = {theModTeleTrack.getParameter(0), theModTeleTrack.getParameter(1), theModTeleTrack.getParameter(2), theModTeleTrack.getParameter(3), theModTeleTrack.getParameter(4), theModTeleTrack.getParameter(5)};
        theModTeleTrack.intersection(planeq0, parFit, xtemp, ytemp, ztemp);
        TVector3 gptemp(xtemp, ytemp, ztemp);
        int ilayernumber=-1;
        bool isleft=false;
        bool isright=false;
        for (int il1=0; il1<8; il1++) {
            if (detId_to_study == LayersDefinition[il1].first)  {
                ilayernumber=il1;
                isright=true;
            } else if (detId_to_study == LayersDefinition[il1].second) {
                ilayernumber=il1;
                isleft=true;
            }

        }
        DQM2_Chi2[detId_to_study]->Fill(theModTeleTrack.getChi2());
        // xtemp = global coordinates
        // --> transformation into LP
        TVector3 lptemp;
        if (isright) lptemp= theSimpleLayersR[ilayernumber].getLocalPointPosition(gptemp);
        else if (isleft) lptemp= theSimpleLayersL[ilayernumber].getLocalPointPosition(gptemp);

        // store the info
        DQM2_TrackPull_X[detId_to_study]->Fill(lptemp.X() - theLP_to_study.X());
        DQM2_TrackPull_Y[detId_to_study]->Fill(lptemp.Y() - theLP_to_study.Y());

        // Fill the tree
        tree_eventN = eventNumber;
        tree_detId2 = detId_to_study;
        tree_biasedResidualX = brX;
        tree_biasedResidualY = brY;
        tree_unbiasedResidualX = lptemp.X() - theLP_to_study.X() ;
        tree_unbiasedResidualY = lptemp.Y() - theLP_to_study.Y() ;

        residualsTree->Fill();

        double Sum2study= (lptemp.X() - theLP_to_study.X())*(lptemp.X() - theLP_to_study.X())
                          + (lptemp.Y() - theLP_to_study.Y())*(lptemp.Y() - theLP_to_study.Y());
        DQM2_Sum2[detId_to_study]->Fill(Sum2study);
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE ( SimpleTrackingClean ) ;
*/
