// -*- C++ -*-
//
// Package:    DQMTelescope/AnaStefano
// Class:      AnaStefano
//
// class AnaStefano AnaStefano.cc DQMTelescope/AnaStefano/plugins/AnaStefano.cc

// Description: [one line class summary]

// Implementation:
//     [Notes on implementation]

// Original Author:  Jeremy Andrea
//      Updated by:  Nikkie Deelen
//         Created:  03.08.2018

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

#include <string> 
#include "TH2F.h"
#include "TTree.h"
#include <TCanvas.h>
#include <TStyle.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection ;

class AnaStefano : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:

    explicit AnaStefano ( const edm::ParameterSet& ) ;
    ~AnaStefano ( ) ;

    static void fillDescriptions ( edm::ConfigurationDescriptions& descriptions ) ;
   
  private:

    virtual void beginJob ( ) override ;
    virtual void analyze ( const edm::Event&, const edm::EventSetup& ) override ;
    virtual void endJob ( ) override ;

    // ----------member data ---------------------------

    edm::EDGetTokenT<TrackCollection> tracksToken_ ;  //used to select what tracks to read from configuration file
      
    edm::Service<TFileService> fs ;
     
    std::vector<TH1F *> DQM_ClusterCharge;
    std::vector<TH1F *> DQM_ClusterSize_X   ;  
    std::vector<TH1F *> DQM_ClusterSize_Y   ; 
    std::vector<TH1F *> DQM_ClusterSize_XY ;
    std::vector<TH1F *> DQM_NumbOfCluster;
    std::vector<TH2F *> DQM_ClusterPosition ;
      
    std::vector<TH2F *> DQM_DigiPosition ;
    std::vector<TH1F *> DQM_NumbOfDigi ;

    TH1F * DQM_NumbOfCluster_Tot ;
    TH1F * DQM_ClusterCharge_Tot ;
      
    TH1F * DQM_NumbOfDigi_Tot ;
      
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi> >         pixeldigiToken_ ;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelclusterToken_ ;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit> >  pixelhitToken_ ;
      
    std::vector<uint32_t > list_of_modules ;
    std::map<int, int> modulesNbr_to_idx ;
    std::map<int , TString> detId_to_moduleName ;
    std::vector<TString> names_of_modules ;

    TTree* clusterTree ;

    // Declaration of leaves types
    Int_t      tree_runNumber ;
    Int_t      tree_lumiSection ;
    Int_t      tree_event ;
    Int_t      tree_detId ;
    Int_t      tree_cluster ;
    Double_t   tree_x ;
    Double_t   tree_y ;

    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation_X ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation_Y ;

   //std::vector<std::vector<TH2F *>> DQM_Correlation_X ;
   //std::vector<std::vector<TH2F *>> DQM_Correlation_Y ;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AnaStefano::AnaStefano(const edm::ParameterSet& iConfig)
 :
  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
   //now do what ever initialization is needed

    list_of_modules.push_back(344200196) ;
    list_of_modules.push_back(344201220) ;
    list_of_modules.push_back(344462340) ;
    list_of_modules.push_back(344463364) ;
    list_of_modules.push_back(344724484) ;
    list_of_modules.push_back(344725508) ;
    list_of_modules.push_back(344986628) ;
    list_of_modules.push_back(344987652) ;
    list_of_modules.push_back(352588804) ;
    list_of_modules.push_back(352589828) ;
    list_of_modules.push_back(352850948) ;
    list_of_modules.push_back(352851972) ;
    list_of_modules.push_back(353113092) ;
    list_of_modules.push_back(353114116) ;
    list_of_modules.push_back(353375236) ;
    list_of_modules.push_back(353376260) ;

    names_of_modules.push_back("M3090") ;
    names_of_modules.push_back("M3124") ;
    names_of_modules.push_back("M3082") ;
    names_of_modules.push_back("M3175") ;
    names_of_modules.push_back("M3009") ;
    names_of_modules.push_back("M3057") ;
    names_of_modules.push_back("M3027") ;
    names_of_modules.push_back("M3074") ;
    names_of_modules.push_back("M3192") ;
    names_of_modules.push_back("M3204") ;
    names_of_modules.push_back("M3226") ;
    names_of_modules.push_back("M3265") ;
    names_of_modules.push_back("M3023") ;
    names_of_modules.push_back("M3239") ;
    names_of_modules.push_back("M3164") ;
    names_of_modules.push_back("M3173") ;

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

    for(unsigned int i=0; i<list_of_modules.size(); i++) modulesNbr_to_idx[list_of_modules[i]] = i;
    
    TFileDirectory sub1 = fs->mkdir(  "run100000" ); // This does not make sense yet

    clusterTree = sub1.make<TTree>("clusterTree", "Cluster Tree");

    // Set branch addresses.
    clusterTree->Branch("runNumber",&tree_runNumber);
    clusterTree->Branch("lumiSection",&tree_lumiSection);
    clusterTree->Branch("event",&tree_event);
    clusterTree->Branch("detId",&tree_detId);
    clusterTree->Branch("cluster",&tree_cluster);
    clusterTree->Branch("x",&tree_x);
    clusterTree->Branch("y",&tree_y);

    TFileDirectory sub2 = sub1.mkdir( "dqmPlots" ) ;
    TFileDirectory sub3 = sub2.mkdir( "runSummary" ) ;
    TFileDirectory sub4 = sub1.mkdir( "correlationPlots" ) ;    

    for(unsigned int i=0; i<list_of_modules.size(); i++){
      
      //TString modulename = std::to_string(names_of_modules[i]) ;
      TString modulename = names_of_modules[i] ;
      TString modulenameNext = "";
      //if ( i < names_of_modules.size() ) TString modulenameNext = std::to_string(names_of_modules[i+1]) ;
      if ( i < names_of_modules.size() ) TString modulenameNext = names_of_modules[i+1] ;

        //TString modulename = std::to_string(list_of_modules[i]) ;
        //auto itModName = detId_to_moduleName.find( list_of_modules[i] );
        //if ( itModName == detId_to_moduleName.end() )
          //continue;

        //TString modulename = itModName->second ;
        //TString modulenameNext = "";
        //if ( i < list_of_modules.size() ) { 
        //  TString modulenameNext = std::to_string(list_of_modules[i+1]) ;
        //  auto itModNameNext = detId_to_moduleName.find( list_of_modules[i+1] );
	//  if ( itModNameNext == detId_to_moduleName.end() )
        //    continue;
        //  modulenameNext = itModNameNext->second ;
        //}//end if list size
      
      std::vector<TH2F *> tmp_vec_x ; // for the corr plots, hence the extra for loop
      std::vector<TH2F *> tmp_vec_y ; // for the corr plots, hence the extra for loop

      for ( unsigned int j=i; j<list_of_modules.size(); j++ ) { // To make sure we do not have double plots.

        //TString modulename0 = std::to_string(names_of_modules[j]) ;
        TString modulename0 = names_of_modules[j] ;
        TString modulenameNext0 = "";
        //if ( j < names_of_modules.size() ) TString modulenameNext0 = std::to_string(names_of_modules[j+1]) ;
        if ( j < names_of_modules.size() ) TString modulenameNext0 = names_of_modules[j+1] ;

        //TString modulename0 = std::to_string(names_of_modules[j]) ;
        //auto itModName0 = detId_to_moduleName.find( list_of_modules[j] );
        //if ( itModName0 == detId_to_moduleName.end() )
        //  continue;

        //TString modulename0 = itModName0->second ;
        //TString modulenameNext0 = "";
        //if ( j < list_of_modules.size() ) { 
        //  //TString modulenameNext0 = std::to_string(names_of_modules[j+1]) ;
        //  auto itModNameNext0 = detId_to_moduleName.find( list_of_modules[j+1] );
	//  if ( itModNameNext0 == detId_to_moduleName.end() )
        //    continue;
        //  modulenameNext0 = itModNameNext0->second ;
        //}//end if list size

        //TH2F* DQM_Correlation_X_tmp = sub4.make<TH2F>( ( "DQM_Correlation_X_" + modulename + "_" + modulename0).Data(), ( "X-Correlation between " + modulename + " " + modulename0 ).Data(), 416., 0., 416., 416., 0., 416. ) ;
        //TH2F* DQM_Correlation_Y_tmp = sub4.make<TH2F>( ( "DQM_Correlation_Y_" + modulename + "_" + modulename0).Data(), ( "Y-Correlation between " + modulename + " " + modulename0 ).Data(), 160., 0., 160., 160., 0., 160. ) ;

        TH2F* DQM_Correlation_X_tmp = sub4.make<TH2F>( ( "DQM_Correlation_X_" + modulename + "_" + modulename0).Data(), ( "X-Correlation between " + modulename + " and " + modulename0 ).Data(), 160., 0., 160., 160., 0., 160. ) ;
        TH2F* DQM_Correlation_Y_tmp = sub4.make<TH2F>( ( "DQM_Correlation_Y_" + modulename + "_" + modulename0).Data(), ( "Y-Correlation between " + modulename + " and " + modulename0 ).Data(), 416., 0., 416., 416., 0., 416. ) ;

        DQM_Correlation_X_tmp->GetXaxis()->SetTitle("x_" + modulename) ;
        DQM_Correlation_X_tmp->GetYaxis()->SetTitle("x_" + modulename0) ;
        DQM_Correlation_Y_tmp->GetXaxis()->SetTitle("y_" + modulename) ;
        DQM_Correlation_Y_tmp->GetYaxis()->SetTitle("y_" + modulename0) ;

        std::pair<uint32_t, uint32_t> modulePair = std::make_pair ( list_of_modules[i], list_of_modules[j] ) ;
        //std::pair<std::pair, TH2F*> mapPair = std::make_pair ( modulePair, DQM_Correlation_X_tmp ) ;

	DQM_Correlation_X.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation_X_tmp ) ) ;
	DQM_Correlation_Y.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation_Y_tmp ) ) ;

        //DQM_Correlation_X.insert( modulePair, DQM_Correlation_X_tmp ) ;
        //DQM_Correlation_Y.insert( modulePair, DQM_Correlation_Y_tmp ) ;

        //tmp_vec_x.push_back( DQM_Correlation_X_tmp ) ;
        //tmp_vec_y.push_back( DQM_Correlation_Y_tmp ) ;

      }//end for j 
      //DQM_Correlation_X.push_back ( tmp_vec_x ) ;
      //DQM_Correlation_Y.push_back ( tmp_vec_y ) ;

      TH1F* DQM_ClusterCharge_tmp       = sub3.make<TH1F>( ("DQM_ClusterCharge_"+ modulename).Data()  , ("Cluster charge for "+ modulename).Data(),      100,  0., 100000. );
      TH1F* DQM_ClusterSize_X_tmp       = sub3.make<TH1F>( ("DQM_ClusterSize_X_"+ modulename).Data()  , ("X cluster size for "+ modulename).Data(),      30,  0., 30.     );
      TH1F* DQM_ClusterSize_Y_tmp       = sub3.make<TH1F>( ("DQM_ClusterSize_Y_"+ modulename).Data()  , ("Y cluster size for "+ modulename).Data(),      30,  0., 30.     );
      TH1F* DQM_ClusterSize_XY_tmp      = sub3.make<TH1F>( ("DQM_ClusterSize_XY_"+ modulename).Data(),  ("Cluster Size for "  + modulename).Data(),      30,  0., 30.     );
      TH1F* DQM_NumbOfCluster_tmp       = sub3.make<TH1F>( ("DQM_NumbOfCluster_"+ modulename).Data(),   ("number of cluster for "  + modulename).Data(), 30,  0., 30.     );
      TH2F* DQM_ClusterPosition_tmp     = sub3.make<TH2F>( ("DQM_ClusterPosition_"+ modulename).Data(), ("Cluster occupancy per col per row for "+ modulename).Data(),   416,  0., 416., 160, 0, 160	);
      
      DQM_ClusterCharge_tmp->GetXaxis()->SetTitle("Charge (electrons)");
      DQM_ClusterSize_X_tmp->GetXaxis()->SetTitle("size (pixels)");
      DQM_ClusterSize_Y_tmp->GetXaxis()->SetTitle("size (pixels)");	
      DQM_ClusterSize_XY_tmp->GetXaxis()->SetTitle("size (pixels)"); 
      DQM_ClusterPosition_tmp->GetXaxis()->SetTitle("col"); 
      DQM_ClusterPosition_tmp->GetYaxis()->SetTitle("row"); 
      
      DQM_ClusterCharge.push_back(DQM_ClusterCharge_tmp);  
      DQM_ClusterSize_X.push_back(DQM_ClusterSize_X_tmp);   
      DQM_ClusterSize_Y.push_back(DQM_ClusterSize_Y_tmp);    
      DQM_ClusterSize_XY.push_back(DQM_ClusterSize_XY_tmp); 
      DQM_NumbOfCluster.push_back(DQM_NumbOfCluster_tmp);
      DQM_ClusterPosition.push_back(DQM_ClusterPosition_tmp);      
      
      TH2F* DQM_DigiPosition_tmp   = sub3.make<TH2F>( ("DQM_DigiPosition_"+ modulename).Data(), ("Digi occupancy per col per row for "+ modulename).Data(),  416,  0., 416., 160, 0, 160	);
      TH1F* DQM_NumbOfDigi_tmp     = sub3.make<TH1F>( ("DQM_NumbOfDigi"+ modulename).Data(),    ("Number of Digis for "  + modulename).Data(),     30,  0., 30.     );
 

      DQM_DigiPosition.push_back(DQM_DigiPosition_tmp); 
      DQM_NumbOfDigi.push_back(DQM_NumbOfDigi_tmp);

      DQM_DigiPosition_tmp->GetXaxis()->SetTitle("col"); 
      DQM_DigiPosition_tmp->GetYaxis()->SetTitle("row"); 

    }
    
    DQM_NumbOfDigi_Tot    = sub3.make<TH1F>( "DQM_NumbOfDigi_Tot",    "total number of digi"   , 30,  0., 30.);
    DQM_NumbOfCluster_Tot = sub3.make<TH1F>( "DQM_NumbOfCluster_Tot", "total number of cluster", 30,  0., 30.);
    
    pixeldigiToken_    = consumes<edm::DetSetVector<PixelDigi> >        (iConfig.getParameter<edm::InputTag>("PixelDigisLabel"))   ;
    pixelclusterToken_ = consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("PixelClustersLabel"));
    pixelhitToken_     = consumes<edmNew::DetSetVector<SiPixelRecHit> > (iConfig.getParameter<edm::InputTag>("PixelHitsLabel"))    ;
   
}

AnaStefano::~AnaStefano()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
AnaStefano::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  
  //---------------------------------
  //loop on digis
  //---------------------------------
  
  // define iterations (in a map) to count the number of cluster per module in the event
  std::map<int, int> numberofDigi_per_module;  
  int numberofDigi_total = 0;
  for(unsigned int i=0; i<list_of_modules.size(); i++) numberofDigi_per_module[ list_of_modules[i]  ] = 0;
 
  //TCanvas* aCanvas = new TCanvas("eventDisplay", "Event Display", 2000, 500); 
  for( edm::DetSetVector<PixelDigi>::const_iterator DSViter=pixeldigis->begin(); DSViter!=pixeldigis->end(); DSViter++   ) {
      numberofDigi_total++;
      //auto id = DetId(DSViter->detId());
      auto nDigisHere = DSViter->size();
      if (nDigisHere>15) {
        //int tree_runNumber = myEvId.run();
        int tree_lumiSection = myEvId.luminosityBlock();
        int tree_event = myEvId.event();
        //int tree_id = id;
        if ((tree_lumiSection==9)&&(tree_event==14965)) {
	//DQM_DigiPosition[ modulesNbr_to_idx[int(id.rawId())]]->Reset();
          edm::DetSet<PixelDigi>::const_iterator begin=DSViter->begin();
          edm::DetSet<PixelDigi>::const_iterator end  =DSViter->end();
          for(edm::DetSet<PixelDigi>::const_iterator iter=begin;iter!=end;++iter) {
            float x = iter->column();                // x position
            float y = iter->row();                   // y position
            float adc = iter->adc();                 // ~ energy release
	    //DQM_DigiPosition[ modulesNbr_to_idx[int(id.rawId())]]->Fill(x, y, adc);
	    std::cout << "Nice event hit x=" << x << " y=" <<y<<" adc=" <<adc <<std::endl;
          }
        //aCanvas->cd();
        //DQM_DigiPosition[ modulesNbr_to_idx[int(id.rawId())]]->Draw("colz");
        }
        //aCanvas->SaveAs(Form("event_bigCluster_Run%d_LS%d_Ev%d_detId%d.png", tree_runNumber, tree_lumiSection, tree_event, tree_id));
      }
   }
  //delete aCanvas;
  
}


// ------------ method called once each job just before starting event loop  ------------
void
AnaStefano::beginJob()
{
  gStyle->SetPalette(56);
  gStyle->SetOptStat(0);
}

// ------------ method called once each job just after ending the event loop  ------------
void
AnaStefano::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnaStefano::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(AnaStefano);
