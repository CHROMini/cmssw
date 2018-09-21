// Package:    DQMTelescope/AnaNikkie
// Class:      AnaNikkie
//
// class AnaNikkie AnaNikkie.cc DQMTelescope/AnaNikkie/plugins/AnaNikkie.cc

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

using reco::TrackCollection ;

class AnaNikkie : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:

    explicit AnaNikkie ( const edm::ParameterSet& ) ;
    ~AnaNikkie ( ) ;

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
     
    std::map< uint32_t, TH1F* > DQM_ClusterCharge ;
    std::map< uint32_t, TH1F* > DQM_ClusterSize_X ;    
    std::map< uint32_t, TH1F* > DQM_ClusterSize_Y ;    
    std::map< uint32_t, TH1F* > DQM_ClusterSize_XY ;
    std::map< uint32_t, TH1F* > DQM_NumbOfClusters_per_Event;
    std::map< uint32_t, TH2F* > DQM_ClusterPosition ;
    std::map< uint32_t, TH2F* > DQM_Hits_per_Pixel_per_Event ;

    // Correlation plots for the telescope
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation_X ;
    std::map< std::pair<uint32_t, uint32_t>, TH2F*> DQM_Correlation_Y ;

    // 3D Tree 
    TTree* cluster3DTree ;

    Int_t      tree_runNumber ;
    Int_t      tree_lumiSection ;
    Int_t      tree_event ;
    Int_t      tree_detId ;
    Int_t      tree_cluster ;
    Double_t   tree_x ;
    Double_t   tree_y ;
    Double_t   tree_z ;
    TString    tree_modName ;
    Long64_t   tree_maxEntries = 10000 ;

    // detId versus moduleName
    std::map<int , TString> detId_to_moduleName ;

} ;

/////////////////////
// Class functions //
/////////////////////

AnaNikkie::AnaNikkie( const edm::ParameterSet& iConfig ) : tracksToken_( consumes<TrackCollection>( iConfig.getUntrackedParameter<edm::InputTag>( "tracks" ) ) ) {

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

  TFileDirectory sub1 = fs -> mkdir ( "clusterTree3D" ) ; 
  
  cluster3DTree = sub1.make<TTree> ("clusterTree3D", "3D Cluster Tree") ;
  TFileDirectory sub2 = fs -> mkdir ( "dqmPlots" ) ;
  TFileDirectory sub3 = fs -> mkdir ( "correlationPlots" ) ;  

  // Set branch addresses.
  cluster3DTree -> Branch ( "runNumber", &tree_runNumber ) ;
  cluster3DTree -> Branch ( "lumiSection", &tree_lumiSection ) ;
  cluster3DTree -> Branch ( "event", &tree_event ) ;
  cluster3DTree -> Branch ( "detId", &tree_detId ) ;
  cluster3DTree -> Branch ( "modName", &tree_modName ) ;
  cluster3DTree -> Branch ( "cluster", &tree_cluster ) ;
  cluster3DTree -> Branch ( "x", &tree_x ) ;
  cluster3DTree -> Branch ( "y", &tree_y ) ;
  cluster3DTree -> Branch ( "z", &tree_z ) ;
  cluster3DTree -> SetCircular ( tree_maxEntries ) ;

  for ( std::map<int, TString>::iterator it = detId_to_moduleName.begin(); it != detId_to_moduleName.end(); it++ ) {
    
    TString modulename = it -> second ;

    std::vector<TH2F *> tmp_vec_x ; // for the corr plots, hence the extra for loop
    std::vector<TH2F *> tmp_vec_y ; // for the corr plots, hence the extra for loop

    // Make the correlation plots
    for ( std::map<int, TString>::iterator jt = it; jt != detId_to_moduleName.end(); jt++ ) { // jt=it to make sure we do not have double plots.
   
      TString modulename0 = jt -> second ;

      TH2F* DQM_Correlation_X_tmp = sub3.make<TH2F>( ( "DQM_Correlation_X_" + modulename + "_" + modulename0).Data(), ( "X-Correlation between " + modulename + " and " + modulename0 ).Data(), 160., 0., 160., 160., 0., 160. ) ;
      TH2F* DQM_Correlation_Y_tmp = sub3.make<TH2F>( ( "DQM_Correlation_Y_" + modulename + "_" + modulename0).Data(), ( "Y-Correlation between " + modulename + " and " + modulename0 ).Data(), 416., 0., 416., 416., 0., 416. ) ;

      DQM_Correlation_X_tmp->GetXaxis()->SetTitle("x_" + modulename) ;
      DQM_Correlation_X_tmp->GetYaxis()->SetTitle("x_" + modulename0) ;
      DQM_Correlation_Y_tmp->GetXaxis()->SetTitle("y_" + modulename) ;
      DQM_Correlation_Y_tmp->GetYaxis()->SetTitle("y_" + modulename0) ;

      std::pair<uint32_t, uint32_t> modulePair = std::make_pair ( it->first, jt->first ) ;

      DQM_Correlation_X.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation_X_tmp ) ) ;
      DQM_Correlation_Y.insert ( std::pair < std::pair<uint32_t, uint32_t>, TH2F*>( modulePair, DQM_Correlation_Y_tmp ) ) ;

    }//end for j 

    // Make the DQM plots
    TH1F* DQM_ClusterCharge_tmp = sub2.make<TH1F>( ( "DQM_ClusterCharge_" + modulename ).Data(), ( "Cluster charge for " + modulename).Data(), 100, 0., 100000. );
    TH1F* DQM_ClusterSize_X_tmp = sub2.make<TH1F>( ( "DQM_ClusterSize_X_" + modulename ).Data(), ( "X cluster size for " + modulename).Data(), 30, 0., 30. );
    TH1F* DQM_ClusterSize_Y_tmp = sub2.make<TH1F>( ( "DQM_ClusterSize_Y_" + modulename ).Data(), ( "Y cluster size for " + modulename ).Data(), 30, 0., 30. ) ;
    TH1F* DQM_ClusterSize_XY_tmp = sub2.make<TH1F>( ( "DQM_ClusterSize_XY_" + modulename ).Data(), ( "Cluster Size for "  + modulename ).Data(), 30, 0., 30. ) ;
    TH1F* DQM_NumbOfClusters_per_Event_tmp = sub2.make<TH1F>( ("DQM_NumbOfClusters_per_Event_" + modulename).Data(), ("number of clusters for "  + modulename).Data(), 30, 0., 30. );
    TH2F* DQM_ClusterPosition_tmp = sub2.make<TH2F>( ( "DQM_ClusterPosition_" + modulename ).Data(), ( "Cluster occupancy per col per row for " + modulename ).Data(), 416, 0., 416., 160, 0., 160. ) ;
    TH2F* DQM_Hits_per_Pixel_per_Event_tmp = sub2.make<TH2F>( ( "DQM_HitsPerEventPerPixel_" + modulename ).Data(), ("Hits per event per pixel for " + modulename ).Data(), 4160, 0., 4160., 100000, 0., 2. ) ;  
    DQM_ClusterCharge_tmp -> GetXaxis ( ) ->SetTitle ( "Charge (electrons)" ) ;
    DQM_ClusterSize_X_tmp -> GetXaxis ( ) ->SetTitle ( "size (pixels)" ) ;
    DQM_ClusterSize_Y_tmp -> GetXaxis ( ) ->SetTitle ( "size (pixels)" ) ;	
    DQM_ClusterSize_XY_tmp -> GetXaxis ( ) ->SetTitle ( "size (pixels)" ) ; 
    DQM_ClusterPosition_tmp -> GetXaxis ( ) ->SetTitle ( "col" ) ; 
    DQM_ClusterPosition_tmp -> GetYaxis ( ) ->SetTitle ( "row" ) ; 
    DQM_NumbOfClusters_per_Event_tmp -> GetXaxis ( ) -> SetTitle ( "Number of clusters / event " ) ;
    DQM_NumbOfClusters_per_Event_tmp -> GetYaxis ( ) -> SetTitle ( "Count" ) ;
    DQM_Hits_per_Pixel_per_Event_tmp -> GetXaxis ( ) -> SetTitle ( "Pixel" ) ;
    DQM_Hits_per_Pixel_per_Event_tmp -> GetYaxis ( ) -> SetTitle ( "Number of hits / event" ) ;

    DQM_ClusterCharge.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterCharge_tmp ) ) ;  
    DQM_ClusterSize_X.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_X_tmp ) ) ;
    DQM_ClusterSize_Y.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_Y_tmp ) ) ;
    DQM_ClusterSize_XY.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_ClusterSize_XY_tmp ) ) ;
    DQM_NumbOfClusters_per_Event.insert ( std::pair< uint32_t, TH1F* >( it->first, DQM_NumbOfClusters_per_Event_tmp ) );
    DQM_ClusterPosition.insert ( std::pair< uint32_t, TH2F* >( it->first, DQM_ClusterPosition_tmp ) ) ;      
    DQM_Hits_per_Pixel_per_Event.insert ( std::pair< uint32_t, TH2F* >( it->first, DQM_Hits_per_Pixel_per_Event_tmp ) ) ;

  }//end for it
  
  pixeldigiToken_    = consumes<edm::DetSetVector<PixelDigi> >        (iConfig.getParameter<edm::InputTag>("PixelDigisLabel"))   ;
  pixelclusterToken_ = consumes<edmNew::DetSetVector<SiPixelCluster> >(iConfig.getParameter<edm::InputTag>("PixelClustersLabel"));
  pixelhitToken_     = consumes<edmNew::DetSetVector<SiPixelRecHit> > (iConfig.getParameter<edm::InputTag>("PixelHitsLabel"))    ;
   
}//end AnaNikkie()

AnaNikkie::~AnaNikkie()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void AnaNikkie::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  using namespace edm;
   
  EventID myEvId = iEvent.id();
   
  /* Handle<TrackCollection> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  for(TrackCollection::const_iterator itTrack = tracks->begin();
    itTrack != tracks->end();
    ++itTrack) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = itTrack->charge();
  }*/
  
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

  /*for( edm::DetSetVector<PixelDigi>::const_iterator DSViter=pixeldigis->begin(); DSViter!=pixeldigis->end(); DSViter++   ) { 
    
    edm::DetSet<PixelDigi>::const_iterator begin=DSViter->begin();
    edm::DetSet<PixelDigi>::const_iterator end  =DSViter->end();
    
    auto id = DetId(DSViter->detId());
    
    std::map< int, int > numHits_per_Channel ;
    
    for(edm::DetSet<PixelDigi>::const_iterator iter=begin;iter!=end;++iter) {
  
      // First get the channel number of the hits
      int channel = iter -> channel ( ) ;
      
      // Add one to the number of hits of the right channel
      if ( numHits_per_Channel.count ( channel ) > 0 ) {
        numHits_per_Channel[ channel ] ++ ;
      } else {
        numHits_per_Channel.insert ( std::pair< int, int >( channel, 1 ) ) ;
      }//end if channel

    }//end for Digis   

    for ( std::map< int, int >::iterator it = numHits_per_Channel.begin(); it != numHits_per_Channel.end(); it++  )  DQM_Hits_per_Pixel_per_Event[ id.rawId ( ) ] -> Fill ( it -> first, it -> second ) ;

  }//end for Detectors
*/
  //---------------------------------
  // Loop over the clusters
  //---------------------------------

  for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end();DSViter++   ) {

    edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();
      
    auto id = DetId(DSViter->detId());

    for( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter2=pixelclusters->begin(); DSViter2!=pixelclusters->end();DSViter2++   ) {

      edmNew::DetSet<SiPixelCluster>::const_iterator begin2=DSViter2->begin();
      edmNew::DetSet<SiPixelCluster>::const_iterator end2  =DSViter2->end();

      auto id2 = DetId(DSViter2->detId());

      for(edmNew::DetSet<SiPixelCluster>::const_iterator iter=begin;iter!=end;++iter) {

        float x = iter->x();                   // barycenter x position
        float y = iter->y();                   // barycenter y position
        int row = x - 0.5, col = y - 0.5;

        for(edmNew::DetSet<SiPixelCluster>::const_iterator iter2=begin2;iter2!=end2;++iter2) {

          float x2 = iter2->x();                   // barycenter x position
          float y2 = iter2->y();                   // barycenter y position
          int row2 = x2 - 0.5, col2 = y2 - 0.5;
	  
          std::pair<uint32_t, uint32_t> modulePair = std::make_pair ( id.rawId(), id2.rawId() ) ;

          auto itHistMap = DQM_Correlation_X.find(modulePair);
	        
          if ( itHistMap == DQM_Correlation_X.end() ) continue;
          itHistMap->second->Fill ( row, row2 ) ;

          auto itHistMap2 = DQM_Correlation_Y.find(modulePair);
          if ( itHistMap2 == DQM_Correlation_Y.end() ) continue;
          itHistMap2->second->Fill ( col, col2 ) ;

        }//end for clusters2
      }//end for cluster
    }//end for first detectors2
  }//end for first detectors

  for ( edmNew::DetSetVector< SiPixelCluster >::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end(); DSViter++ ) {

    edmNew::DetSet<SiPixelCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiPixelCluster>::const_iterator end  =DSViter->end();
      
    auto id = DetId(DSViter->detId());

    // Counting the number of clusters for this detector and this event
    float numberOfClustersForThisModule = 0 ;

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
	
      numberOfClustersForThisModule++ ;

    }//end for clusters in detector

    // So now we fille the number of clusters per module
    DQM_NumbOfClusters_per_Event[ id.rawId() ] -> Fill ( numberOfClustersForThisModule ) ;
    
  }//end for detectors

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

  for ( edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=pixelclusters->begin(); DSViter!=pixelclusters->end();DSViter++ ) {

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
	
      // Let's fill in the tree
      tree_runNumber = myEvId.run();
      tree_lumiSection = myEvId.luminosityBlock();
      tree_event = myEvId.event();
      tree_detId = detId;
      tree_cluster = iCluster++;	
      tree_x = x;
      tree_y = y;
      tree_z = z;
      tree_modName = detId_to_moduleName[detId];
      cluster3DTree->Fill();

    } //end for clusters of the first detector
  } //end for first detectors

}//end analyze()

// ------------ method called once each job just before starting event loop  ------------
void AnaNikkie::beginJob ( ) {
}//end void beginJob ( ) 

// ------------ method called once each job just after ending the event loop  ------------
void AnaNikkie::endJob ( ) {
}//end void endJobd ( ) 

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AnaNikkie::fillDescriptions ( edm::ConfigurationDescriptions& descriptions ) {

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

//define this as a plug-in
DEFINE_FWK_MODULE ( AnaNikkie ) ;
