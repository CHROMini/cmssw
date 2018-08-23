// -*- C++ -*-
//
// Package:    DQMTelescope/PixelTelescope
// Class:      PixelTelescope
//
/**\class PixelTelescope PixelTelescope.cc DQMTelescope/PixelTelescope/plugins/PixelTelescope.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jeremy Andrea
//         Created:  Thu, 14 Jun 2018 14:43:12 GMT
//
//


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

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class PixelTelescope : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PixelTelescope(const edm::ParameterSet&);
      ~PixelTelescope();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      
      edm::Service<TFileService> fs;
      
      std::vector<TH1F *> DQM_ClusterCharge;
      std::vector<TH1F *> DQM_ClusterSize_X   ;  
      std::vector<TH1F *> DQM_ClusterSize_Y   ; 
      std::vector<TH1F *> DQM_ClusterSize_XY ;
      std::vector<TH1F *> DQM_NumbOfCluster;
      std::vector<TH2F *> DQM_ClusterPosition ;
//      TTree* clusterTree; 
      
      std::vector<TH2F *> DQM_DigiPosition ;
      std::vector<TH1F *> DQM_NumbOfDigi;
      
      
      TH1F * DQM_NumbOfCluster_Tot;
      TH1F * DQM_ClusterCharge_Tot;
      
      
      TH1F * DQM_NumbOfDigi_Tot;
      
      edm::EDGetTokenT<edm::DetSetVector<PixelDigi> >         pixeldigiToken_;
      edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelclusterToken_;
      edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit> >  pixelhitToken_;
      
      std::vector<uint32_t > list_of_modules;
      std::map<int, int> modulesNbr_to_idx;

      // Declaration of leaves types
//      Int_t           tree_runNumber;
//      Int_t           tree_lumiSection;
//      Int_t           tree_event;
//      Int_t           tree_detId;
//      Int_t           tree_cluster;
//      Double_t        tree_x;
//      Double_t        tree_y;
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
PixelTelescope::PixelTelescope(const edm::ParameterSet& iConfig)
 :
  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
   //now do what ever initialization is needed
   
   
    /*TFileDirectory subDQMData   = fs->mkdir( "DQMData" );
    TFileDirectory subRunNum    = subDQMData.mkdir( "Run 1" );
    TFileDirectory subWorkSpace = subRunNum.mkdir( "PixelTelescope" );
    TFileDirectory subRunSumm   = subWorkSpace.mkdir( "Run summary" );*/
    
    
    list_of_modules.push_back(344200196);
    list_of_modules.push_back(344201220); 
    list_of_modules.push_back(344462340);
    list_of_modules.push_back(344463364); 
    list_of_modules.push_back(344724484);
    list_of_modules.push_back(344725508); 
    list_of_modules.push_back(344986628);
    list_of_modules.push_back(344987652); 
    list_of_modules.push_back(352588804);
    list_of_modules.push_back(352589828);
    list_of_modules.push_back(352850948); 
    list_of_modules.push_back(352851972);
    list_of_modules.push_back(353113092); 
    list_of_modules.push_back(353114116); 
    list_of_modules.push_back(353375236); 
    list_of_modules.push_back(353376260); 

 
    for(unsigned int i=0; i<list_of_modules.size(); i++) modulesNbr_to_idx[list_of_modules[i]] = i;
    
    TFileDirectory sub1 = fs->mkdir(  "Run 100000" );

//    clusterTree = sub1.make<TTree>("clusterTree", "Cluster Tree");

    // Set branch addresses.
//    clusterTree->Branch("runNumber",&tree_runNumber);
//    clusterTree->Branch("lumiSection",&tree_lumiSection);
//    clusterTree->Branch("event",&tree_event);
//    clusterTree->Branch("detId",&tree_detId);
//    clusterTree->Branch("cluster",&tree_cluster);
//    clusterTree->Branch("x",&tree_x);
//    clusterTree->Branch("y",&tree_y);

    TFileDirectory sub2 = sub1.mkdir( "PixelTelescope" );
    TFileDirectory sub3 = sub2.mkdir( "Run summary" );
    
    for(unsigned int i=0; i<list_of_modules.size(); i++){
      
      TString modulename = std::to_string(list_of_modules[i]) ;
      
      TH1F* DQM_ClusterCharge_tmp       = sub3.make<TH1F>( ("DQM_ClusterCharge_"+ modulename).Data()  , ("Cluster charge for detID "+ modulename).Data(),      100,  0., 100000. );
      TH1F* DQM_ClusterSize_X_tmp       = sub3.make<TH1F>( ("DQM_ClusterSize_X_"+ modulename).Data()  , ("X cluster size for detId "+ modulename).Data(),      30,  0., 30.     );
      TH1F* DQM_ClusterSize_Y_tmp       = sub3.make<TH1F>( ("DQM_ClusterSize_Y_"+ modulename).Data()  , ("Y cluster size for detId "+ modulename).Data(),      30,  0., 30.     );
      TH1F* DQM_ClusterSize_XY_tmp      = sub3.make<TH1F>( ("DQM_ClusterSize_XY_"+ modulename).Data(),  ("Cluster Size for detId "  + modulename).Data(),      30,  0., 30.     );
      TH1F* DQM_NumbOfCluster_tmp       = sub3.make<TH1F>( ("DQM_NumbOfCluster_"+ modulename).Data(),   ("number of cluster for detId "  + modulename).Data(), 30,  0., 30.     );
      TH2F* DQM_ClusterPosition_tmp     = sub3.make<TH2F>( ("DQM_ClusterPosition_"+ modulename).Data(), ("Cluster occupancy per col per row for detId "+ modulename).Data(),   416,  0., 416., 160, 0, 160	);
      
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
      
      TH2F* DQM_DigiPosition_tmp   = sub3.make<TH2F>( ("DQM_DigiPosition_"+ modulename).Data(), ("Digi occupancy per col per row for detId "+ modulename).Data(),  416,  0., 416., 160, 0, 160	);
      TH1F* DQM_NumbOfDigi_tmp     = sub3.make<TH1F>( ("DQM_NumbOfDigi"+ modulename).Data(),    ("Number of cluster for detId "  + modulename).Data(),     30,  0., 30.     );
 
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


PixelTelescope::~PixelTelescope()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PixelTelescope::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
//   EventID myEvId = iEvent.id();
   
   /* Handle<TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    for(TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    }
   
  */
  
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
  
    
  
  //define iterations (in a map) to count the number of cluster per module in the event
  std::map<int, int> numberofDigi_per_module;  
  int numberofDigi_total = 0;
  for(unsigned int i=0; i<list_of_modules.size(); i++) numberofDigi_per_module[ list_of_modules[i]  ] = 0;
  
  
  for( edm::DetSetVector<PixelDigi>::const_iterator DSViter=pixeldigis->begin(); DSViter!=pixeldigis->end(); DSViter++   ) {
      
      edm::DetSet<PixelDigi>::const_iterator begin=DSViter->begin();
      edm::DetSet<PixelDigi>::const_iterator end  =DSViter->end();
      
      auto id = DetId(DSViter->detId());
      //std::cout << "id " << id.rawId ()  << std::endl;
      
      for(edm::DetSet<PixelDigi>::const_iterator iter=begin;iter!=end;++iter) {
         
	
         float x = iter->column();                   // barycenter x position
         float y = iter->row();                   // barycenter y position
	 DQM_DigiPosition[ modulesNbr_to_idx[int(id.rawId())]]->Fill(x, y);  
	 
	 std::cout << " modulesNbr_to_idx[int(id.rawId())] " <<  modulesNbr_to_idx[int(id.rawId())] << "  " << int(id.rawId()) <<  std::endl;
	 numberofDigi_per_module[ modulesNbr_to_idx[int(id.rawId())]]++;
	 
	 
	 //std::cout << "   " << numberofDigi_per_module[ modulesNbr_to_idx[int(id.rawId())]] << "  " <<  modulesNbr_to_idx[int(id.rawId())] << "  " << id.rawId() << std::endl;
	 numberofDigi_total++;

      }
  
  }
  DQM_NumbOfDigi_Tot->Fill(numberofDigi_total);
  
  
  //fill the number of cluster in the event per module
  for(unsigned int i=0; i<list_of_modules.size(); i++){
    DQM_NumbOfDigi[i]->Fill( numberofDigi_per_module[i]) ;
    //std::cout <<  numberofDigi_per_module[i] << std::endl;
  }
  
  
  
  
  
  
  //---------------------------------
  //loop on clusters
  //---------------------------------
  unsigned int numberOfClusters = 0;
  //unsigned int numberOfFpixClusters = 0;
  
  
  //define iterations (in a map) to count the number of cluster per module in the event
  std::map<int, int> numberofCluster_per_module; 
  int numberofCulster_total = 0; 
  for(unsigned int i=0; i<list_of_modules.size(); i++) numberofCluster_per_module[ list_of_modules[i]  ] = 0;
  
  for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end();DSViter++   ) {
      numberOfClusters++;
      
      edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
      edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();
      


      
      auto id = DetId(DSViter->detId());
      //std::cout << "id " << id.rawId ()  << std::endl;
      //int iCluster = 0;
      for(edmNew::DetSet<SiPixelCluster>::const_iterator iter=begin;iter!=end;++iter) {
        
	
        float x = iter->x();                   // barycenter x position
        float y = iter->y();                   // barycenter y position
        int size = iter->size();               // total size of cluster (in pixels)
        int sizeX = iter->sizeX();             // size of cluster in x-iterrection
        int sizeY = iter->sizeY();             // size of cluster in y-iterrection
        /*int minPixelRow = iter->minPixelRow(); // min x index
        int maxPixelRow = iter->maxPixelRow(); // max x index
        int minPixelCol = iter->minPixelCol(); // min y index
        int maxPixelCol = iter->maxPixelCol(); // max y index
	*/
	
        // Let's fill in the tree
//        tree_runNumber = myEvId.run();
//        tree_lumiSection = myEvId.luminosityBlock();
//        tree_event = myEvId.event();
//        tree_detId = id;
//        tree_cluster = iCluster++;	
//        tree_x = x;
//        tree_y = y;
//        clusterTree->Fill();
	
        int row = x-0.5, col = y -0.5;
	
        DQM_ClusterCharge[ modulesNbr_to_idx[id.rawId()]]->Fill(iter->charge());
	DQM_ClusterSize_X[ modulesNbr_to_idx[id.rawId()]]->Fill(sizeX);   
	DQM_ClusterSize_Y[ modulesNbr_to_idx[id.rawId()]]->Fill(sizeY);     
	DQM_ClusterSize_XY[ modulesNbr_to_idx[id.rawId()]]->Fill(size);   
	DQM_ClusterPosition[ modulesNbr_to_idx[id.rawId()]]->Fill(col, row);  
	
	numberofCluster_per_module[ modulesNbr_to_idx[id.rawId()]]++;
	
	numberofCulster_total++;
      }
    }
      //if(endcap) numberOfFpixClusters++;
      //float charge = 0.001*(di->charge()); // total charge of cluster
      
  DQM_NumbOfCluster_Tot->Fill(numberofCulster_total);
  if(numberofCulster_total != 0 ) std::cout << "number of cluster " << numberOfClusters << std::endl;
  
  
  //fill the number of cluster in the event per module
  for(unsigned int i=0; i<list_of_modules.size(); i++) DQM_NumbOfCluster[i]->Fill( numberofCluster_per_module[i]) ;
  
  

  //---------------------------------
  //loop on hits
  //---------------------------------
  
    
  
  for( edmNew::DetSetVector<SiPixelRecHit>::const_iterator DSViter=pixelhits->begin(); DSViter!=pixelhits->end(); DSViter++   ) {
      
      edmNew::DetSet<SiPixelRecHit>::const_iterator begin=DSViter->begin();
      edmNew::DetSet<SiPixelRecHit>::const_iterator end  =DSViter->end();
    
      for(edmNew::DetSet<SiPixelRecHit>::const_iterator iter=begin;iter!=end;++iter) {
         
      }
  
  }

}


// ------------ method called once each job just before starting event loop  ------------
void
PixelTelescope::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PixelTelescope::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PixelTelescope::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(PixelTelescope);
