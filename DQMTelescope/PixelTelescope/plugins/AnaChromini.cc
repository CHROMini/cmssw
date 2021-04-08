// -*- C++ -*-
//
// Package:    DQMTelescope/AnaChromini
// Class:      AnaChromini
//
// class AnaChromini AnaChromini.cc DQMTelescope/PixelTelescope/plugins/AnaChromini.cc

// Description: [one line class summary]

// Implementation:
//     [Notes on implementation]

// Original Author:  Caroline Collard
//      Created:  18.12.2019
//      Starting from a copy of AnaNikkie.cc
//      Updated by:  

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
// test
#include "DQM/SiPixelPhase1Common/interface/SiPixelCoordinates.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"

#include <cstring>
#include <string> 
#include <TH2F.h>
#include <TTree.h>
#include <TString.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection ;

class AnaChromini : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:

    explicit AnaChromini ( const edm::ParameterSet& ) ;
    ~AnaChromini ( ) ;

    static void fillDescriptions ( edm::ConfigurationDescriptions& descriptions ) ;
   
  private:

    virtual void beginJob ( ) override ;
    virtual void analyze ( const edm::Event&, const edm::EventSetup& ) override ;
    virtual void endJob ( ) override ;

    // ----------member data ---------------------------

    edm::EDGetTokenT<TrackCollection> tracksToken_ ;  //used to select what tracks to read from configuration file
      
    edm::Service<TFileService> fs ;
     
    // information to stored in the output file
    std::map< uint32_t, TH1F* > DQM_ClusterCharge ;
    std::map< uint32_t, TH1F* > DQM_ClusterSize_X ;
    std::map< uint32_t, TH1F* > DQM_ClusterSize_Y ;
    std::map< uint32_t, TH1F* > DQM_ClusterSize_XY ;
    std::map< uint32_t, TH1F* > DQM_NumbOfClusters_per_Event;
    std::map< uint32_t, TH2F* > DQM_ClusterPosition ;

    std::map< uint32_t, TH2F* >  DQM_DigiPosition ;
    std::map< uint32_t, TH1F* > DQM_NumbOfDigi ;

    TH1F * DQM_NumbOfCluster_Tot ;
    TH1F * DQM_ClusterCharge_Tot ;
    TH1F * DQM_NumbOfDigi_Tot ;
      
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi> >         pixeldigiToken_ ;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelclusterToken_ ;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit> >  pixelhitToken_ ;
      

    // detId versus moduleName and other definitions
    std::vector<uint32_t > list_of_modules ;
    std::map<int , TString> detId_to_moduleName ;
    std::map<TString , int> moduleName_to_detID ;
    std::map<TString , std::vector<double> > moduleName_to_position ;
    std::vector<TString> names_of_modules ;

    // Correlation plots for the telescope
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation_X ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation_Y ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Distance_XY ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Distance_XYcm ;

    TTree* cluster3DTree ;

    // Declaration of leaves types
    Int_t      tree_runNumber ;
    Int_t      tree_lumiSection ;
    Int_t      tree_event ;
    Int_t      tree_detId ;
    Int_t      tree_cluster ;
    Double_t   tree_charge ;
    Int_t      tree_size ;
    Double_t   tree_x ;
    Double_t   tree_y ;
    Double_t   tree_z ;
    TString    tree_modName;

};

/////////////////////
// Class functions //
/////////////////////

AnaChromini::AnaChromini(const edm::ParameterSet& iConfig)
 :
  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
 //now do what ever initialization is needed
 


/* Version with Chromini4Modul
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344463364, "M3704") ) ;
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344724484, "M3672") ) ;
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344725508, "M3558") ) ;
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344462340, "M3028") ) ;
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3704", 344463364) ) ;
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3672", 344724484) ) ;
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3558", 344725508) ) ;
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3028", 344462340) ) ;
  list_of_modules.push_back(344463364) ;
  list_of_modules.push_back(344724484) ; 
  list_of_modules.push_back(344725508) ; 
  list_of_modules.push_back(344462340) ; 
  names_of_modules.push_back("M3704") ;
  names_of_modules.push_back("M3672") ;
  names_of_modules.push_back("M3558") ;
  names_of_modules.push_back("M3028") ;
*/

/* old version === last one before Cyrce
  // map for the detId vs Module name:
  //  fiber 6 (white) = 344463364
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344463364, "M3558") ) ;
  //   fiber 3 (green) = 344462340
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344462340, "M3211") ) ;
  //  fiber 1 (blue) = 344725508
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344725508, "M3028") ) ;
  //  fiber 12 (aqua)  = 344724484
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344724484, "M3704") ) ;
  //  fiber 11 (pink)  = 344987652
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344987652, "M4643") ) ;
  //  fiber 10 (purple)  = 344986628
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344986628, "M3672") ) ;
*/

 
  // Cyrce Design
  // map for the detId vs Module name:
  //  fiber 6 (white) = 344463364
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344463364, "M3558") ) ;
  //   fiber 3 (green) = 344462340
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344462340, "M4643") ) ;
  //  fiber 1 (blue) = 344725508
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344725508, "M3028") ) ;
  //  fiber 12 (aqua)  = 344724484
  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344724484, "M3672") ) ;
  //  fiber 11 (pink)  = 344987652
//  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344987652, "M3211") ) ;
  //  fiber 10 (purple)  = 344986628
//  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344986628, "M3704") ) ;

/*
//  
//  fiber 8 = 344201220
//  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344201220, "M3704") ) ;
//  fiber 5 = 344200196
//  detId_to_moduleName.insert( std::pair<uint32_t, TString>(344200196, "M3704") ) ;
//  */

/*
  // Module name to detID
  // fiber 6 
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3558", 344463364) ) ;
  //  fiber 3 
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3211", 344462340) ) ;
  //  fiber 1
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3028", 344725508) ) ;
  //  fiber 12 
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3704", 344724484) ) ;
  // fiber 11
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M4643", 344987652) ) ;
  // fiber 10
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3672", 344986628) ) ;
*/

// Cyrce Design
  // Module name to detID
  // fiber 6 
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3558", 344463364) ) ;
  //  fiber 3 
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M4643", 344462340) ) ;
  //  fiber 1
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3028", 344725508) ) ;
  //  fiber 12 
  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3672", 344724484) ) ;
  // fiber 11
//  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3211", 344987652) ) ;
  // fiber 10
//  moduleName_to_detID.insert( std::pair<TString, uint32_t>("M3704", 344986628) ) ;


  list_of_modules.push_back(344463364) ;
  list_of_modules.push_back(344462340) ;
  list_of_modules.push_back(344725508) ; 
  list_of_modules.push_back(344724484) ; 
//  list_of_modules.push_back(344987652) ; 
//  list_of_modules.push_back(344986628) ; 

  names_of_modules.push_back("M3558") ;
  names_of_modules.push_back("M4643") ;
  names_of_modules.push_back("M3028") ;
  names_of_modules.push_back("M3672") ;
//  names_of_modules.push_back("M3211") ;
//  names_of_modules.push_back("M3704") ;



  TFileDirectory sub1 = fs->mkdir(  "run100000" ); // This does not make sense yet

  cluster3DTree = sub1.make<TTree>("cluster3DTree", "3D Cluster Tree");
  TFileDirectory sub2 = sub1.mkdir( "dqmPlots" ) ;
  TFileDirectory sub3 = sub2.mkdir( "runSummary" ) ;
  TFileDirectory sub4 = sub1.mkdir( "correlationPlots" ) ;  

  // Set branch addresses.
  cluster3DTree->Branch("runNumber",&tree_runNumber);
  cluster3DTree->Branch("lumiSection",&tree_lumiSection);
  cluster3DTree->Branch("event",&tree_event);
  cluster3DTree->Branch("detId",&tree_detId);
  cluster3DTree->Branch("modName",&tree_modName);
  cluster3DTree->Branch("cluster",&tree_cluster);
  cluster3DTree->Branch("charge",&tree_charge);
  cluster3DTree->Branch("size",&tree_size);
  cluster3DTree->Branch("x",&tree_x);
  cluster3DTree->Branch("y",&tree_y);
  cluster3DTree->Branch("z",&tree_z);

  //for(unsigned int i=0; i<list_of_modules.size(); i++) modulesNbr_to_idx[list_of_modules[i]] = i;

  for ( std::map<int, TString>::iterator it = detId_to_moduleName.begin(); it != detId_to_moduleName.end(); it++ ) {
    
    TString modulename = it -> second ;

    std::vector<TH2F *> tmp_vec_x ; // for the corr plots, hence the extra for loop
    std::vector<TH2F *> tmp_vec_y ; // for the corr plots, hence the extra for loop

    // Make the correlation plots
    for ( std::map<int, TString>::iterator jt = it; jt != detId_to_moduleName.end(); jt++ ) { // jt=it to make sure we do not have double plots.
   
      TString modulename0 = jt -> second ;

      // inversion with respect to Chromie
      // in Chromini : Y from 0 to 160 and X from 0 to 416 !
      TH2F* DQM_Correlation_Y_tmp = sub4.make<TH2F>( ( "DQM_Correlation_Y_" + modulename + "_" + modulename0).Data(), ( "Y-Correlation between " + modulename + " and " + modulename0 ).Data(), 160., 0., 160., 160., 0., 160. ) ;
      TH2F* DQM_Correlation_X_tmp = sub4.make<TH2F>( ( "DQM_Correlation_X_" + modulename + "_" + modulename0).Data(), ( "X-Correlation between " + modulename + " and " + modulename0 ).Data(), 416., 0., 416., 416., 0., 416. ) ;
      TH2F* DQM_Distance_XY_tmp = sub4.make<TH2F>( ( "DQM_Distance_XY_" + modulename + "_" + modulename0).Data(), ( "Y distance vs X distance between " + modulename + " and " + modulename0 ).Data(), 100., -50., 50., 100., -50., 50. ) ;
      TH2F* DQM_Distance_XYcm_tmp = sub4.make<TH2F>( ( "DQM_Distance_XYcm_" + modulename + "_" + modulename0).Data(), ( "Y distance vs X distance between " + modulename + " and " + modulename0 ).Data(), 100., -0.5, 0.5, 100., -0.5, 0.5 ) ;

      DQM_Correlation_X_tmp->GetXaxis()->SetTitle("x_" + modulename) ;
      DQM_Correlation_X_tmp->GetYaxis()->SetTitle("x_" + modulename0) ;
      DQM_Correlation_Y_tmp->GetXaxis()->SetTitle("y_" + modulename) ;
      DQM_Correlation_Y_tmp->GetYaxis()->SetTitle("y_" + modulename0) ;
      DQM_Distance_XY_tmp->GetXaxis()->SetTitle("x_" + modulename + " - x_" + modulename0) ;
      DQM_Distance_XY_tmp->GetYaxis()->SetTitle("y_" + modulename + " - y_" + modulename0) ;
      DQM_Distance_XYcm_tmp->GetXaxis()->SetTitle("x_" + modulename + " - x_" + modulename0) ;
      DQM_Distance_XYcm_tmp->GetYaxis()->SetTitle("y_" + modulename + " - y_" + modulename0) ;

      std::pair<uint32_t, uint32_t> modulePair = std::make_pair ( it->first, jt->first ) ;

      DQM_Correlation_X.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation_X_tmp ) ) ;
      DQM_Correlation_Y.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation_Y_tmp ) ) ;
      DQM_Distance_XY.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Distance_XY_tmp ) ) ;
      DQM_Distance_XYcm.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Distance_XYcm_tmp ) ) ;

    }//end for j 

    // Make the DQM plots
    TH1F* DQM_ClusterCharge_tmp = sub3.make<TH1F>( ("DQM_ClusterCharge_"+ modulename).Data(), ("Cluster charge for "+ modulename).Data(), 200, 0., 200000. );
    TH1F* DQM_ClusterSize_Y_tmp = sub3.make<TH1F>( ("DQM_ClusterSize_Y_"+ modulename).Data(), ("X cluster size for "+ modulename).Data(), 30, 0., 30. );
    TH1F* DQM_ClusterSize_X_tmp = sub3.make<TH1F>( ("DQM_ClusterSize_X_"+ modulename).Data(), ("Y cluster size for "+ modulename).Data(), 30, 0., 30. );
    TH1F* DQM_ClusterSize_XY_tmp = sub3.make<TH1F>( ("DQM_ClusterSize_XY_"+ modulename).Data(), ("Cluster Size for "  + modulename).Data(), 30, 0., 30. );
    TH1F* DQM_NumbOfClusters_per_Event_tmp = sub3.make<TH1F>( ("DQM_NumbOfClusters_per_Event_" + modulename).Data(), ("number of clusters for "  + modulename).Data(), 30, 0., 30. );
    TH2F* DQM_ClusterPosition_tmp = sub3.make<TH2F>( ("DQM_ClusterPosition_"+ modulename).Data(), ("Cluster occupancy per col per row for "+ modulename).Data(), 416, 0., 416., 160, 0, 160 );
      
    DQM_ClusterCharge_tmp->GetXaxis()->SetTitle("Charge (electrons)");
    DQM_ClusterSize_X_tmp->GetXaxis()->SetTitle("size (pixels)");
    DQM_ClusterSize_Y_tmp->GetXaxis()->SetTitle("size (pixels)");	
    DQM_ClusterSize_XY_tmp->GetXaxis()->SetTitle("size (pixels)"); 
    DQM_ClusterPosition_tmp->GetXaxis()->SetTitle("col"); 
    DQM_ClusterPosition_tmp->GetYaxis()->SetTitle("row"); 
    DQM_NumbOfClusters_per_Event_tmp -> GetXaxis ( ) -> SetTitle ( "Number of clusters / event " ) ;
    DQM_NumbOfClusters_per_Event_tmp -> GetYaxis ( ) -> SetTitle ( "Count" ) ;

    DQM_ClusterCharge.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterCharge_tmp ) ) ;
    DQM_ClusterSize_X.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_X_tmp ) ) ;
    DQM_ClusterSize_Y.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_Y_tmp ) ) ;
    DQM_ClusterSize_XY.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_XY_tmp ) ) ;
    DQM_NumbOfClusters_per_Event.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_NumbOfClusters_per_Event_tmp ) );
    DQM_ClusterPosition.insert ( std::pair< uint32_t, TH2F* >( it->first, DQM_ClusterPosition_tmp ) ) ;

      
    TH2F* DQM_DigiPosition_tmp = sub3.make<TH2F>( ("DQM_DigiPosition_"+ modulename).Data(), ("Digi occupancy per col per row for "+ modulename).Data(),  416,  0., 416., 160, 0, 160	);
    TH1F* DQM_NumbOfDigi_tmp = sub3.make<TH1F>( ("DQM_NumbOfDigi"+ modulename).Data(),    ("Number of cluster for "  + modulename).Data(), 30,  0., 30. );
 
    DQM_DigiPosition.insert ( std::pair< uint32_t, TH2F* >( it->first, DQM_DigiPosition_tmp)); 
    DQM_NumbOfDigi.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_NumbOfDigi_tmp));

    DQM_DigiPosition_tmp->GetXaxis()->SetTitle("col"); 
    DQM_DigiPosition_tmp->GetYaxis()->SetTitle("row"); 

  }//end for it
    
  DQM_NumbOfDigi_Tot    = sub3.make<TH1F>( "DQM_NumbOfDigi_Tot",    "total number of digi"   , 30,  0., 30.);
  DQM_NumbOfCluster_Tot = sub3.make<TH1F>( "DQM_NumbOfCluster_Tot", "total number of cluster", 30,  0., 30.);
  
  pixeldigiToken_    = consumes<edm::DetSetVector<PixelDigi> >        (iConfig.getParameter<edm::InputTag>("PixelDigisLabel"))   ;
  pixelclusterToken_ = consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("PixelClustersLabel"));
  pixelhitToken_     = consumes<edmNew::DetSetVector<SiPixelRecHit> > (iConfig.getParameter<edm::InputTag>("PixelHitsLabel"))    ;
   
}//end AnaChromini()

AnaChromini::~AnaChromini()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
AnaChromini::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
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
  
  //define iterations (in a map) to count the number of cluster per module in the event
  std::map<int, int> numberofDigi_per_module;  
  int numberofDigi_total = 0;
  for ( std::map<int, TString>::iterator it = detId_to_moduleName.begin(); it != detId_to_moduleName.end(); it++ ) numberofDigi_per_module.insert( std::pair< int, int >( it->first, 0 ) ) ;;
  
  for( edm::DetSetVector<PixelDigi>::const_iterator DSViter=pixeldigis->begin(); DSViter!=pixeldigis->end(); DSViter++   ) {
      
    edm::DetSet<PixelDigi>::const_iterator begin=DSViter->begin();
    edm::DetSet<PixelDigi>::const_iterator end  =DSViter->end();
      
    auto id = DetId(DSViter->detId());
//    if (myEvId.event()==1 || myEvId.event()==1000) {
        bool trouve=0;
        for (unsigned int jj=0; jj<list_of_modules.size(); jj++) {
              if ((uint32_t) id==list_of_modules[jj]) trouve=1;
        }
      if (trouve==0)
        std::cout << " DetId " << (int) id << std::endl;
//    }
      
    for(edm::DetSet<PixelDigi>::const_iterator iter=begin;iter!=end;++iter) {
         
	
      float x = iter->column(); // barycenter x position
      float y = iter->row();    // barycenter y position

      // Apply Mask
      if ((uint32_t) id==344724484 && y>=80 && 312<= x && x<364) continue;
      
      DQM_DigiPosition[ id.rawId() ] -> Fill ( x, y ) ;
      numberofDigi_per_module[ id.rawId() ]++ ;

//      const DetId& detaidi=DSViter->detId();
//      const PixelDigi* digiprov = iter;
//        int       pseudo_roc_num = uint64_t(1<<16) * id.rawId() + (1<<8) * (y/80) + x/52;
//      std::cout << "  Test Digi Position " << id.rawId()  << " x  "   << x  <<   "  y "  << y << "    ROC "  << pseudo_roc_num << std::endl;
//SiPixelCoordinates::roc(detaidi, digiprov) << std::endl;

      numberofDigi_total++;

//      if (numberofDigi_per_module[ id.rawId()]>0) std::cout << "Digi in " << int(id.rawId()) << " =  " <<  detId_to_moduleName[int(id.rawId())] << std::endl;

    }//end for DetSet
  }//end for DetSetVector

  DQM_NumbOfDigi_Tot->Fill(numberofDigi_total);
  
  //fill the number of cluster in the event per module
  for ( std::map<int, int>::iterator it = numberofDigi_per_module.begin(); it != numberofDigi_per_module.end(); it++ ) DQM_NumbOfDigi[ it->first ] -> Fill ( it->second ) ;
  
  //---------------------------------
  //loop on clusters
  //---------------------------------
  unsigned int numberOfClusters = 0;
  
  //define iterations (in a map) to count the number of cluster per module in the event
  std::map<int, int> numberofCluster_per_module; 
  int numberofCulster_total = 0; 
  for(unsigned int i=0; i<list_of_modules.size(); i++) numberofCluster_per_module[ list_of_modules[i] ] = 0;


//
//
 for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end();DSViter++   ) {

    edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();

    auto id = DetId(DSViter->detId());
    const PixelGeomDetUnit *pixdet1 = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(id);
    LocalPoint lp1(-9999., -9999., -9999.);

    for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter2=pixelclusters->begin(); DSViter2!=pixelclusters->end();DSViter2++   ) {

      edmNew::DetSet<SiPixelCluster>::const_iterator begin2=DSViter2->begin();
      edmNew::DetSet<SiPixelCluster>::const_iterator end2  =DSViter2->end();

      auto id2 = DetId(DSViter2->detId());
      const PixelGeomDetUnit *pixdet2 = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(id2);
      LocalPoint lp2(-9999., -9999., -9999.);

       for(edmNew::DetSet<SiPixelCluster>::const_iterator iter=begin;iter!=end;++iter) {
        float x = iter->x();                   // barycenter x position for Chromie, y for Chromini
        float y = iter->y();                   // barycenter y position for Chromie, x for Chromini
        int row = x , col = y ;

        // Apply MASK
        if ((uint32_t) id==344724484 && x>=80 && 312<= y && y<364) continue;


        PixelClusterParameterEstimator::ReturnType params1 = cpe.getParameters(*iter,*pixdet1);
        lp1 = std::get<0>(params1);
        const Surface& surface1 = tracker->idToDet(id)->surface();
        GlobalPoint gp1 = surface1.toGlobal(lp1);
        float xgb1 = gp1.x();
        float ygb1 = gp1.y();

        for(edmNew::DetSet<SiPixelCluster>::const_iterator iter2=begin2;iter2!=end2;++iter2) {

          float x2 = iter2->x();                   // barycenter x position for Chromie, y for Chromini
          float y2 = iter2->y();                   // barycenter y position for Chromie, x for Chromini
          int row2 = x2, col2 = y2 ;

          // Apply MASK
          if ((uint32_t) id2==344724484 && x2>=80 && 312<= y2 && y2<364) continue;
	  
          std::pair<uint32_t, uint32_t> modulePair = std::make_pair ( id.rawId(), id2.rawId() ) ;

          bool test_roc_malade=false;
/*
          if (((uint32_t) id==344724484 && (uint32_t) id2==344462340)||((uint32_t) id2==344724484 && (uint32_t) id==344462340) ) { //3672 - 4643
               if ((uint32_t) id==344462340 && x<80 && 208 <= y && y< 260) test_roc_malade=true;
               if ((uint32_t) id2==344462340 && x2<80 && 208 <= y2 && y2< 260) test_roc_malade=true;
          }
*/

          PixelClusterParameterEstimator::ReturnType params2 = cpe.getParameters(*iter2,*pixdet2);
          lp2 = std::get<0>(params2);
          const Surface& surface2 = tracker->idToDet(id2)->surface();
          GlobalPoint gp2 = surface2.toGlobal(lp2);
          float xgb2 = gp2.x();
          float ygb2 = gp2.y();

          auto itHistMap = DQM_Correlation_Y.find(modulePair);
	        
          if ( itHistMap == DQM_Correlation_Y.end() ) continue;
          if(!test_roc_malade) itHistMap->second->Fill ( row, row2 ) ;

          auto itHistMap2 = DQM_Correlation_X.find(modulePair);
          if ( itHistMap2 == DQM_Correlation_X.end() ) continue;
          if(!test_roc_malade) itHistMap2->second->Fill ( col, col2 ) ;

          auto itHistMap3 = DQM_Distance_XY.find(modulePair);
          if ( itHistMap3 == DQM_Distance_XY.end() ) continue;
          if(!test_roc_malade) itHistMap3->second->Fill (col- col2, row-row2) ;

          auto itHistMap4 = DQM_Distance_XYcm.find(modulePair);
          if ( itHistMap4 == DQM_Distance_XYcm.end() ) continue;
          if(!test_roc_malade) itHistMap4->second->Fill (xgb1- xgb2, ygb1-ygb2) ;

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

        // Apply MASK
        if ((uint32_t) id==344724484 && x>=80 && 312<= y && y<364) continue;

      DQM_ClusterCharge[ id.rawId() ] -> Fill ( iter->charge() ) ;
      DQM_ClusterSize_X[ id.rawId() ] -> Fill ( sizeY ) ;
      DQM_ClusterSize_Y[ id.rawId() ] -> Fill ( sizeX ) ;
      DQM_ClusterSize_XY[ id.rawId() ] -> Fill ( size ) ;
      DQM_ClusterPosition[ id.rawId() ] -> Fill ( col, row ) ;

      numberofCulster_total++;

      clustersPerModule[ id.rawId() ] ++ ;

    }//end for clusters in detector  
  }//end for detectors

  for ( std::map<uint32_t, int>::iterator it = clustersPerModule.begin(); it != clustersPerModule.end(); it++ ) DQM_NumbOfClusters_per_Event[ it->first ] -> Fill ( it->second ) ;


  // Now we fill the 3D cluster tree
  for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end();DSViter++   ) {

    edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();

    // Surface of (detId) and use the surface to convert 2D to 3D      
    auto detId = DetId(DSViter->detId());
    int iCluster=0;
    const PixelGeomDetUnit *pixdet = (const PixelGeomDetUnit*) tkgeom->idToDetUnit(detId);
    LocalPoint lp(-9999., -9999., -9999.);

    // Then loop on the clusters of the module
    for(edmNew::DetSet<SiPixelCluster>::const_iterator itCluster=begin;itCluster!=end;++itCluster) {
      PixelClusterParameterEstimator::ReturnType params = cpe.getParameters(*itCluster,*pixdet);
      lp = std::get<0>(params);

      const Surface& surface = tracker->idToDet(detId)->surface();
      GlobalPoint gp = surface.toGlobal(lp);

        // Apply MASK
        double xbary=itCluster->x();
        double ybary=itCluster->y();
        if ( (uint32_t) detId==344724484 && xbary>=80 && 312<= ybary && ybary<364) continue;

      double x=0, y=0, z=0;
      x = gp.x();
      y = gp.y();
      z = gp.z();
	
      // Let's fill in the tree
      tree_runNumber = myEvId.run();
      tree_lumiSection = myEvId.luminosityBlock();
      tree_event = myEvId.event();
      tree_detId = detId;
      tree_cluster = iCluster++;	
      tree_charge = itCluster->charge();	
      tree_size = itCluster->size();	
      tree_x = x;
      tree_y = y;
      tree_z = z;
      tree_modName = detId_to_moduleName[detId];
      cluster3DTree->Fill();

      } //end for clusters of the first detector
    } //end for first detectors
      
  DQM_NumbOfCluster_Tot->Fill(numberofCulster_total);
//  if(numberofCulster_total != 0 ) std::cout << "number of cluster " << numberOfClusters << std::endl;
  
  
  //---------------------------------
  //loop on hits
  //---------------------------------
  
  for( edmNew::DetSetVector<SiPixelRecHit>::const_iterator DSViter=pixelhits->begin(); DSViter!=pixelhits->end(); DSViter++   ) {
      
    edmNew::DetSet<SiPixelRecHit>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiPixelRecHit>::const_iterator end  =DSViter->end();
    
    for(edmNew::DetSet<SiPixelRecHit>::const_iterator iter=begin;iter!=end;++iter) {

      // Here we do something with the hits.
         
    }//end for DetSet
  }//end for DetSetVector
}//end analyze()


// ------------ method called once each job just before starting event loop  ------------
void
AnaChromini::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
AnaChromini::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnaChromini::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnaChromini);
