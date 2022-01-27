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

    edm::Service<TFileService> fs ;
    
    // Get specific data
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> pixelDigiToken_ ;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit>> pixelHitToken_ ;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterToken_ ;
    edm::EDGetTokenT<TrackCollection> trackToken_ ;
     
    // Module IDs and names
    std::vector<TString> module_names;
    std::vector<uint32_t> module_det_ids;
    std::map<int, TString> det_id_to_module_name;
    std::map<TString, int> module_name_to_det_id;
    std::map<TString, std::vector<double>> module_name_to_position;
     
    // Content of the output ROOT file
    
    // Histograms per module
    std::map<uint32_t, TH1F*> hist_dqm_digi_column_position;
    std::map<uint32_t, TH1F*> hist_dqm_digi_row_position;
    std::map<uint32_t, TH2F*> hist_dqm_digi_2d_position;
    std::map<uint32_t, TH1F*> hist_dqm_n_digis;
    
    std::map<uint32_t, TH1F*> hist_dqm_cluster_local_x_position;
    std::map<uint32_t, TH1F*> hist_dqm_cluster_local_y_position;
    std::map<uint32_t, TH2F*> hist_dqm_cluster_local_2d_position;
    
    std::map<uint32_t, TH1F*> hist_dqm_cluster_global_x_position;
    std::map<uint32_t, TH1F*> hist_dqm_cluster_global_y_position;
    std::map<uint32_t, TH2F*> hist_dqm_cluster_global_2d_position;
    
    std::map<uint32_t, TH1F*> hist_dqm_n_clusters;
    std::map<uint32_t, TH1F*> hist_dqm_cluster_charge;
    std::map<uint32_t, TH1F*> hist_dqm_cluster_size;
    
    // Histograms showing correlations per module pair
    std::map<std::pair<uint32_t, uint32_t>, TH2F*> hist_local_x_correlation;
    std::map<std::pair<uint32_t, uint32_t>, TH2F*> hist_global_x_correlation;
    std::map<std::pair<uint32_t, uint32_t>, TH2F*> hist_local_y_correlation;
    std::map<std::pair<uint32_t, uint32_t>, TH2F*> hist_global_y_correlation;
    std::map<std::pair<uint32_t, uint32_t>, TH1F*> hist_local_delta_x;
    std::map<std::pair<uint32_t, uint32_t>, TH1F*> hist_global_delta_x;
    std::map<std::pair<uint32_t, uint32_t>, TH1F*> hist_local_delta_y;
    std::map<std::pair<uint32_t, uint32_t>, TH1F*> hist_global_delta_y;
    std::map<std::pair<uint32_t, uint32_t>, TH2F*> hist_local_delta_x_vs_delta_y;
    std::map<std::pair<uint32_t, uint32_t>, TH2F*> hist_global_delta_x_vs_delta_y;
    
    // Histograms for all modules
    TH1F* hist_n_digis_tot;
    TH1F* hist_n_clusters_tot;
    TH1F* hist_cluster_charge_tot;
    
    // Cluster tree
    TTree* cluster_tree;
    
    // Declaration of leaves types
    Int_t tree_run_number;
    Int_t tree_lumi_section;
    Int_t tree_event;
    Int_t tree_det_id;
    TString tree_module_name;
    Int_t tree_n_clusters;
    Double_t tree_cluster_charge;
    Int_t tree_cluster_size;
    Double_t tree_cluster_local_x_position;
    Double_t tree_cluster_local_y_position;
    Double_t tree_cluster_global_x_position;
    Double_t tree_cluster_global_y_position;
    Double_t tree_cluster_global_z_position;
     
    // information to stored in the output file
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::map<uint32_t, TH1F*> DQM_ClusterSize_X ;
    std::map<uint32_t, TH1F*> DQM_ClusterSize_Y ;
    std::map<uint32_t, TH1F*> DQM_ClusterSize_XY ;
   
      
     
    

};

/////////////////////
// Class functions //
/////////////////////

AnaChromini::AnaChromini(const edm::ParameterSet& iConfig)
 :
  trackToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
 //now do what ever initialization is needed
 
 
  pixelDigiToken_ = consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("PixelDigisLabel"));
  pixelHitToken_ = consumes<edmNew::DetSetVector<SiPixelRecHit>>(iConfig.getParameter<edm::InputTag>("PixelHitsLabel"));
  pixelClusterToken_ = consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("PixelClustersLabel"));
 
  // CHROMini setup
    
  // Det ID to module name
  // fiber 6
  det_id_to_module_name.insert(std::pair<uint32_t, TString>(344463364, "M3558"));
  // fiber 3
  det_id_to_module_name.insert(std::pair<uint32_t, TString>(344462340, "M4643"));
  // fiber 1
  det_id_to_module_name.insert(std::pair<uint32_t, TString>(344725508, "M3028"));
  // fiber 12
  det_id_to_module_name.insert(std::pair<uint32_t, TString>(344724484, "M3672"));
  
  TFileDirectory cluster_dir = fs->mkdir("cluster");
  cluster_tree = cluster_dir.make<TTree>("cluster_tree", "Cluster tree");
  
  // Set branch addresses.
  cluster_tree->Branch("run_number", &tree_run_number);
  cluster_tree->Branch("lumi_section", &tree_lumi_section);
  cluster_tree->Branch("event", &tree_event);
  cluster_tree->Branch("det_id", &tree_det_id);
  cluster_tree->Branch("module_name", &tree_module_name);
  cluster_tree->Branch("n_clusters", &tree_n_clusters);
  cluster_tree->Branch("cluster_charge", &tree_cluster_charge);
  cluster_tree->Branch("cluster_size", &tree_cluster_size);
  cluster_tree->Branch("cluster_local_x", &tree_cluster_local_x_position);
  cluster_tree->Branch("cluster_local_y", &tree_cluster_local_y_position);
  cluster_tree->Branch("cluster_global_x", &tree_cluster_global_x_position);
  cluster_tree->Branch("cluster_global_y", &tree_cluster_global_y_position);
  cluster_tree->Branch("cluster_global_z", &tree_cluster_global_z_position);
  
  TFileDirectory dqm_dir = fs->mkdir("dqm");
  TFileDirectory correlation_dir = fs->mkdir("correlation");
  TFileDirectory run_dir = fs->mkdir("run_summary");
  
  // inversion with respect to Chromie
  // in Chromini : Y from 0 to 160 and X from 0 to 416 !
  int n_columns = 416;
  int n_rows = 160;
  
  int n_x_bins = 200;
  float min_x = -1.5;
  float max_x = 1.5;
  
  int n_y_bins = 100;
  float min_y = -0.5;
  float max_y = 1.0;
  
  int n_delta_x_bins = 100;
  int min_delta_x = -50;
  int max_delta_x = 50;
  //float min_delta_x_in_cm = -0.5;
  //float max_delta_x_in_cm = 0.5;
  
  int n_delta_y_bins = 100;
  int min_delta_y = -50;
  int max_delta_y = 50;
  //float min_delta_y_in_cm = -0.5;
  //float max_delta_y_in_cm = 0.5;
  
  int n_digis_bins = 20;
  int n_clusters_bins = 5;
  
  for (auto iter=det_id_to_module_name.begin(); iter!=det_id_to_module_name.end(); iter++)
  {
  	// Access key
    	auto imodule_det_id = iter->first;
	// Access value
	auto imodule_name = iter->second;
	
	module_det_ids.push_back(imodule_det_id);
	module_names.push_back(imodule_name);
	module_name_to_det_id.insert(std::pair<TString, uint32_t>(imodule_name, imodule_det_id));
	
	// DQM plots
	
	TH1F* tmp_hist_dqm_digi_column_position = dqm_dir.make<TH1F>(("digi_column_position_"+ imodule_name).Data(), ("Digi Column for "+imodule_name).Data(),n_columns , -0.5, n_columns-0.5);
	TH1F* tmp_hist_dqm_digi_row_position = dqm_dir.make<TH1F>(("digi_row_position_"+ imodule_name).Data(), ("Digi Row for "+imodule_name).Data(),n_rows , -0.5, n_rows-0.5);
	TH2F* tmp_hist_dqm_digi_2d_position = dqm_dir.make<TH2F>(("digi_2d_position_"+ imodule_name).Data(), ("Digi 2D position for "+imodule_name).Data(), n_columns, -0.5, n_columns-0.5, n_rows, -0.5, n_rows-0.5);
	TH1F* tmp_hist_dqm_n_digis = dqm_dir.make<TH1F>(("n_digis_"+ imodule_name).Data(), ("Number of digis for "+imodule_name).Data(),n_digis_bins , -0.5, n_digis_bins-0.5);
	
	TH1F* tmp_hist_dqm_cluster_local_x_position = dqm_dir.make<TH1F>(("cluster_local_x_position_"+imodule_name).Data(), ("Cluster X_{"+imodule_name+"} [local unity]").Data(),n_columns , -0.5, n_columns-0.5);
	TH1F* tmp_hist_dqm_cluster_local_y_position = dqm_dir.make<TH1F>(("cluster_local_y_position_"+imodule_name).Data(), ("Cluster Y_{"+imodule_name+"} [local unity]").Data(),n_rows , -0.5, n_rows-0.5);
	TH2F* tmp_hist_dqm_cluster_local_2d_position = dqm_dir.make<TH2F>(("cluster_local_2d_position_"+imodule_name).Data(), ("Cluster 2D position for "+imodule_name+" [local unity]").Data(), n_columns, -0.5, n_columns-0.5, n_rows, -0.5, n_rows-0.5);
	
	TH1F* tmp_hist_dqm_cluster_global_x_position = dqm_dir.make<TH1F>(("cluster_global_x_position_"+imodule_name).Data(), ("Cluster X_{"+imodule_name+"} [cm]").Data(),n_x_bins , min_x, max_x);
	TH1F* tmp_hist_dqm_cluster_global_y_position = dqm_dir.make<TH1F>(("cluster_global_y_position_"+imodule_name).Data(), ("Cluster Y_{"+imodule_name+"} [cm]").Data(),n_y_bins , min_y, max_y);
	TH2F* tmp_hist_dqm_cluster_global_2d_position = dqm_dir.make<TH2F>(("cluster_global_2d_position_"+imodule_name).Data(), ("Cluster 2D position for "+imodule_name+" [cm]").Data(), n_x_bins, min_x, max_x, n_y_bins, min_y, max_y);
	
	
	TH1F* tmp_hist_dqm_n_clusters = dqm_dir.make<TH1F>(("n_clusters_"+ imodule_name).Data(), ("Number of clusters for "+imodule_name).Data(),n_clusters_bins , -0.5, n_clusters_bins-0.5);
	TH1F* tmp_hist_dqm_cluster_charge = dqm_dir.make<TH1F>(("cluster_charge_"+ imodule_name).Data(), ("Cluster charge for "+imodule_name).Data(),100 , 0, 1000);
	TH1F* tmp_hist_dqm_cluster_size = dqm_dir.make<TH1F>(("cluster_size_"+ imodule_name).Data(), ("Size of clusters for "+imodule_name).Data(),n_clusters_bins , -0.5, n_clusters_bins-0.5);
	
	
	tmp_hist_dqm_digi_column_position->GetXaxis()->SetTitle("Column Number for "+imodule_name);
	tmp_hist_dqm_digi_row_position->GetXaxis()->SetTitle("Row Number for "+imodule_name);
	
	tmp_hist_dqm_digi_2d_position->GetXaxis()->SetTitle("Column Number for "+imodule_name);
      	tmp_hist_dqm_digi_2d_position->GetYaxis()->SetTitle("Row Number for "+imodule_name);
	
	tmp_hist_dqm_cluster_local_x_position->GetXaxis()->SetTitle("X_{"+imodule_name+"} [local unity]");
	tmp_hist_dqm_cluster_local_y_position->GetXaxis()->SetTitle("Y_{"+imodule_name+"} [local unity]");
	
	tmp_hist_dqm_cluster_local_2d_position->GetXaxis()->SetTitle("X_{"+imodule_name+"} [local unity]");
      	tmp_hist_dqm_cluster_local_2d_position->GetYaxis()->SetTitle("Y_{"+imodule_name+"} [local unity]");
	
	tmp_hist_dqm_cluster_global_x_position->GetXaxis()->SetTitle("X_{"+imodule_name+"} [cm]");
	tmp_hist_dqm_cluster_global_y_position->GetXaxis()->SetTitle("Y_{"+imodule_name+"} [cm]");
	
	tmp_hist_dqm_cluster_global_2d_position->GetXaxis()->SetTitle("X_{"+imodule_name+"} [cm]");
      	tmp_hist_dqm_cluster_global_2d_position->GetYaxis()->SetTitle("Y_{"+imodule_name+"} [cm]");
	
	tmp_hist_dqm_cluster_charge->GetXaxis()->SetTitle("Cluster charge [electrons]");
	tmp_hist_dqm_cluster_size->GetXaxis()->SetTitle("Cluster size [pixels]");
      
    	
	hist_dqm_digi_column_position.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_digi_column_position));
    	hist_dqm_digi_row_position.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_digi_row_position));
    	hist_dqm_digi_2d_position.insert(std::pair<uint32_t, TH2F*>(imodule_det_id, tmp_hist_dqm_digi_2d_position));
    	hist_dqm_n_digis.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_n_digis));
    
    	hist_dqm_cluster_local_x_position.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_cluster_local_x_position));
    	hist_dqm_cluster_local_y_position.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_cluster_local_y_position));
    	hist_dqm_cluster_local_2d_position.insert(std::pair<uint32_t, TH2F*>(imodule_det_id, tmp_hist_dqm_cluster_local_2d_position));
	
	hist_dqm_cluster_global_x_position.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_cluster_global_x_position));
    	hist_dqm_cluster_global_y_position.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_cluster_global_y_position));
    	hist_dqm_cluster_global_2d_position.insert(std::pair<uint32_t, TH2F*>(imodule_det_id, tmp_hist_dqm_cluster_global_2d_position));
	
	
    	hist_dqm_n_clusters.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_n_clusters));
    	hist_dqm_cluster_charge.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_cluster_charge));
    	hist_dqm_cluster_size.insert(std::pair<uint32_t, TH1F*>(imodule_det_id, tmp_hist_dqm_cluster_size));
	
	
	for (auto jter=det_id_to_module_name.begin(); jter!=det_id_to_module_name.end(); jter++)
	{
	    // Access key
    	    auto jmodule_det_id = jter->first;
	    // Access value
	    auto jmodule_name = jter->second;
		
	    if (imodule_det_id != jmodule_det_id)
	    {
    
    		TH2F* tmp_hist_local_x_correlation = correlation_dir.make<TH2F>(("local_x_correlation_"+imodule_name+"_"+jmodule_name).Data(), ("X correlation between "+imodule_name+" and "+ jmodule_name+" [local unity]").Data(), n_columns, -0.5, n_columns-0.5, n_columns, -0.5, n_columns-0.5);
    
		TH2F* tmp_hist_local_y_correlation = correlation_dir.make<TH2F>(("local_y_correlation_"+imodule_name+"_"+jmodule_name).Data(), ("Y correlation between "+imodule_name+" and "+ jmodule_name+" [local unity]").Data(), n_rows, -0.5, n_rows-0.5, n_rows, -0.5, n_rows-0.5);
		
		TH1F* tmp_hist_local_delta_x = correlation_dir.make<TH1F>(("local_delta_x_"+imodule_name+"_"+jmodule_name).Data(), ("#Delta X between "+imodule_name+" and "+ jmodule_name+" [local unity]").Data(), n_delta_x_bins, min_delta_x, max_delta_x);
    
		TH1F* tmp_hist_local_delta_y = correlation_dir.make<TH1F>(("local_delta_y_"+imodule_name+"_"+jmodule_name).Data(), ("#Delta Y between "+imodule_name+" and "+ jmodule_name+" [local unity]").Data(), n_delta_y_bins, min_delta_y, max_delta_y);
		
		TH2F* tmp_hist_local_delta_x_vs_delta_y = correlation_dir.make<TH2F>(("local_delta_x_vs_delta_y_"+imodule_name+"_"+jmodule_name).Data(), ("#Delta X vs #Delta Y between "+imodule_name+" and "+ jmodule_name+" [local unity]").Data(), n_delta_x_bins, min_delta_x, max_delta_x, n_delta_y_bins, min_delta_y, max_delta_y);
		
		
		TH2F* tmp_hist_global_x_correlation = correlation_dir.make<TH2F>(("global_x_correlation_"+imodule_name+"_"+jmodule_name).Data(), ("X correlation between "+imodule_name+" and "+ jmodule_name+" [cm]").Data(), n_x_bins, min_x, max_x, n_x_bins, min_x, max_x);
    
		TH2F* tmp_hist_global_y_correlation = correlation_dir.make<TH2F>(("global_y_correlation_"+imodule_name+"_"+jmodule_name).Data(), ("Y correlation between "+imodule_name+" and "+ jmodule_name+" [cm]").Data(), n_y_bins, min_y, max_y, n_y_bins, min_y, max_y);
		
		TH1F* tmp_hist_global_delta_x = correlation_dir.make<TH1F>(("global_delta_x_"+imodule_name+"_"+jmodule_name).Data(), ("#Delta X between "+imodule_name+" and "+ jmodule_name+" [cm]").Data(), 100, -1, 1);
    
		TH1F* tmp_hist_global_delta_y = correlation_dir.make<TH1F>(("global_delta_y_"+imodule_name+"_"+jmodule_name).Data(), ("#Delta Y between "+imodule_name+" and "+ jmodule_name+" [cm]").Data(), 100, -1, 1);
		
		TH2F* tmp_hist_global_delta_x_vs_delta_y = correlation_dir.make<TH2F>(("global_delta_x_vs_delta_y_"+imodule_name+"_"+jmodule_name).Data(), ("#Delta X vs #Delta Y between "+imodule_name+" and "+ jmodule_name+" [cm]").Data(), 100, -1, 1, 100, -1, 1);
		
		
		tmp_hist_local_x_correlation->GetXaxis()->SetTitle("X_{"+imodule_name+"} [local unity]");
		tmp_hist_local_x_correlation->GetYaxis()->SetTitle("X_{"+jmodule_name+"} [local unity]");
		
		tmp_hist_global_x_correlation->GetXaxis()->SetTitle("X_{"+imodule_name+"} [cm]");
		tmp_hist_global_x_correlation->GetYaxis()->SetTitle("X_{"+jmodule_name+"} [cm]");
		
		tmp_hist_local_y_correlation->GetXaxis()->SetTitle("Y_{"+imodule_name+"} [local unity]");
		tmp_hist_local_y_correlation->GetYaxis()->SetTitle("Y_{"+jmodule_name+"} [local unity]");
		
		tmp_hist_global_y_correlation->GetXaxis()->SetTitle("Y_{"+imodule_name+"} [cm]");
		tmp_hist_global_y_correlation->GetYaxis()->SetTitle("Y_{"+jmodule_name+"} [cm]");
		
		tmp_hist_local_delta_x->GetXaxis()->SetTitle("X_{"+imodule_name+"} - X_{"+jmodule_name+"} [local unity]");
		tmp_hist_global_delta_x->GetXaxis()->SetTitle("X_{"+imodule_name+"} - X_{"+jmodule_name+"} [cm]");
		
		tmp_hist_local_delta_y->GetXaxis()->SetTitle("Y_{"+imodule_name+"} - Y_{"+jmodule_name+"} [local unity]");
		tmp_hist_global_delta_y->GetXaxis()->SetTitle("Y_{"+imodule_name+"} - Y_{"+jmodule_name+"} [cm]");
		
		tmp_hist_local_delta_x_vs_delta_y->GetXaxis()->SetTitle("X_{"+imodule_name+"} - X_{"+jmodule_name+"} [local unity]");
		tmp_hist_local_delta_x_vs_delta_y->GetYaxis()->SetTitle("Y_{"+imodule_name+"} - Y_{"+jmodule_name+"} [local unity]" );
		
		tmp_hist_global_delta_x_vs_delta_y->GetXaxis()->SetTitle("X_{"+imodule_name+"} - X_{"+jmodule_name+"} [cm]");
		tmp_hist_global_delta_x_vs_delta_y->GetYaxis()->SetTitle("Y_{"+imodule_name+"} - Y_{"+jmodule_name+"} [cm]");
		
      
      		std::pair<uint32_t, uint32_t> module_pair = std::make_pair(imodule_det_id, jmodule_det_id);

      		hist_local_x_correlation.insert(std::pair<std::pair<uint32_t, uint32_t>, TH2F*>(module_pair, tmp_hist_local_x_correlation));
		hist_global_x_correlation.insert(std::pair<std::pair<uint32_t, uint32_t>, TH2F*>(module_pair, tmp_hist_global_x_correlation));
		
      		hist_local_y_correlation.insert(std::pair<std::pair<uint32_t, uint32_t>, TH2F*>(module_pair, tmp_hist_local_y_correlation));
		hist_global_y_correlation.insert(std::pair<std::pair<uint32_t, uint32_t>, TH2F*>(module_pair, tmp_hist_global_y_correlation));
      		
		hist_local_delta_x.insert(std::pair<std::pair<uint32_t, uint32_t>, TH1F*>(module_pair, tmp_hist_local_delta_x));
		hist_global_delta_x.insert(std::pair<std::pair<uint32_t, uint32_t>, TH1F*>(module_pair, tmp_hist_global_delta_x));
		
      		hist_local_delta_y.insert(std::pair<std::pair<uint32_t, uint32_t>, TH1F*>(module_pair,  tmp_hist_local_delta_y));
		hist_global_delta_y.insert(std::pair<std::pair<uint32_t, uint32_t>, TH1F*>(module_pair,  tmp_hist_global_delta_y));
		
      		hist_local_delta_x_vs_delta_y.insert(std::pair<std::pair<uint32_t, uint32_t>, TH2F*>(module_pair, tmp_hist_local_delta_x_vs_delta_y));
      		hist_global_delta_x_vs_delta_y.insert(std::pair<std::pair<uint32_t, uint32_t>, TH2F*>(module_pair, tmp_hist_global_delta_x_vs_delta_y));
	
	
	    }
	}
    
  }
  
  hist_n_digis_tot = run_dir.make<TH1F>("hist_n_digis_tot", "Total number of digis", 30,  0., 30.);
  hist_n_clusters_tot = run_dir.make<TH1F>("hist_n_clusters_tot", "Total number of clusters", 10,  -0.5, 9.5);
  hist_cluster_charge_tot = run_dir.make<TH1F>("hist_cluster_charge_tot", "Total cluster charge", 30,  0., 30.);
    
   
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
   
  EventID event_id = iEvent.id();
   
  //get collection of digis
  edm::Handle<edm::DetSetVector<PixelDigi>> pixelDigis;
  iEvent.getByToken(pixelDigiToken_, pixelDigis);
  
  //get collection of clusters
  edm::Handle<edmNew::DetSetVector<SiPixelCluster>> pixelClusters;
  iEvent.getByToken(pixelClusterToken_, pixelClusters);
 
  //get collection or RecHits
  edm::Handle< edmNew::DetSetVector<SiPixelRecHit>> pixelHits;
  iEvent.getByToken(pixelHitToken_, pixelHits);

  // Get the geometry of the tracker for converting the LocalPoint to a GlobalPoint
  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  const TrackerGeometry *tk_geom = &(*tracker); 
  
  // // Choose the CP Estimator that will be used to estimate the LocalPoint of the cluster
  edm::ESHandle<PixelClusterParameterEstimator> cpEstimator;
  iSetup.get<TkPixelCPERecord>().get("PixelCPEGeneric", cpEstimator);
  const PixelClusterParameterEstimator &cpe(*cpEstimator); 
  
  //---------------------------------
  //loop on digis
  //---------------------------------
  
  //define iterations (in a map) to count the number of digi per module in the event
  std::map<int, int> n_digis_per_module;  
  int n_digis_total = 0;
  
  for (std::map<int, TString>::iterator iter=det_id_to_module_name.begin(); iter!=det_id_to_module_name.end(); iter++) n_digis_per_module.insert(std::pair<int, int>(iter->first, 0));
  
  for(edm::DetSetVector<PixelDigi>::const_iterator iDSV=pixelDigis->begin(); iDSV!=pixelDigis->end(); iDSV++) 
  {
      
    edm::DetSet<PixelDigi>::const_iterator DS_begin = iDSV->begin();
    edm::DetSet<PixelDigi>::const_iterator DS_end = iDSV->end();
      
    auto id = DetId(iDSV->detId());
    
    for(edm::DetSet<PixelDigi>::const_iterator idigi=DS_begin; idigi!=DS_end; ++idigi) {
         
	
      float column = idigi->column(); // column
      float row = idigi->row(); // row

      // Apply Mask
      uint32_t det_id = id.rawId();
      if (det_id==344724484 && row>=80 && 312<=column && column<364) continue;
      
      hist_dqm_digi_column_position[det_id]->Fill(column);
      hist_dqm_digi_row_position[det_id]->Fill(row);
      hist_dqm_digi_2d_position[det_id]->Fill(column, row);
      
      n_digis_per_module[det_id]++;
      n_digis_total++;
      
    }//end for DetSet
  }//end for DetSetVector

  hist_n_digis_tot->Fill(n_digis_total);
  
  //fill the number of cluster in the event per module
  for (std::map<int, int>::iterator iter=n_digis_per_module.begin(); iter!=n_digis_per_module.end(); iter++) hist_dqm_n_digis[iter->first]->Fill(iter->second);
  
  //---------------------------------
  //loop on clusters
  //---------------------------------
  
  for(edmNew::DetSetVector<SiPixelCluster>::const_iterator iDSV=pixelClusters->begin(); iDSV!=pixelClusters->end(); iDSV++) 
  {

    edmNew::DetSet<SiPixelCluster>::const_iterator iDS_begin = iDSV->begin();
    edmNew::DetSet<SiPixelCluster>::const_iterator iDS_end = iDSV->end();

    auto id_1 = DetId(iDSV->detId());
    
    const PixelGeomDetUnit *pixel_geom_1 = (const PixelGeomDetUnit*) tk_geom->idToDetUnit(id_1);
    LocalPoint local_point_1(-9999., -9999., -9999.);

    for(edmNew::DetSetVector<SiPixelCluster>::const_iterator jDSV=pixelClusters->begin(); jDSV!=pixelClusters->end(); jDSV++) 
    {

      edmNew::DetSet<SiPixelCluster>::const_iterator jDS_begin = jDSV->begin();
      edmNew::DetSet<SiPixelCluster>::const_iterator jDS_end = jDSV->end();

      auto id_2 = DetId(jDSV->detId());
      const PixelGeomDetUnit *pixel_geom_2 = (const PixelGeomDetUnit*) tk_geom->idToDetUnit(id_2);
      LocalPoint local_point_2(-9999., -9999., -9999.);

       for(edmNew::DetSet<SiPixelCluster>::const_iterator icluster=iDS_begin; icluster!=iDS_end; ++icluster) 
       {
        // Inverse x-y coordinates for CHROMini
       	float local_y_1 = icluster->x(); // barycenter x position for Chromie, y for Chromini
        float local_x_1 = icluster->y(); // barycenter y position for Chromie, x for Chromini
        //float row_1 = x_1 , col_1 = y_1;

        // Apply MASK
        //if ((uint32_t) id_1==344724484 && local_x_1>=80 && 312<=local_y_1 && local_y_1<364) continue;
	if ((uint32_t) id_1==344724484 && local_y_1>=80 && 312<=local_x_1 && local_x_1<364) continue; // inverse x-y


        PixelClusterParameterEstimator::ReturnType parameters_1 = cpe.getParameters(*icluster, *pixel_geom_1);
        local_point_1 = std::get<0>(parameters_1);
        const Surface& surface_1 = tracker->idToDet(id_1)->surface();
        GlobalPoint global_point_1 = surface_1.toGlobal(local_point_1);
        float global_x_1 = global_point_1.x();
        float global_y_1 = global_point_1.y();

        for(edmNew::DetSet<SiPixelCluster>::const_iterator jcluster=jDS_begin; jcluster!=jDS_end; ++jcluster) 
	{
	  // Inverse x-y coordinates for CHROMini
          float local_y_2 = jcluster->x(); // barycenter x position for Chromie, y for Chromini
          float local_x_2 = jcluster->y(); // barycenter y position for Chromie, x for Chromini
          //float row_2 = x_2, col_2 = y_2;

          // Apply MASK
          //if ((uint32_t) id_2==344724484 && local_x_2>=80 && 312<=local_y_2 && local_y_2<364) continue;
	  if ((uint32_t) id_2==344724484 && local_y_2>=80 && 312<=local_x_2 && local_x_2<364) continue; // inverse x-y
	  
          std::pair<uint32_t, uint32_t> module_pair = std::make_pair(id_1.rawId(), id_2.rawId());

  
          PixelClusterParameterEstimator::ReturnType parameters_2 = cpe.getParameters(*jcluster,*pixel_geom_2);
          local_point_2 = std::get<0>(parameters_2);
          const Surface& surface_2 = tracker->idToDet(id_2)->surface();
          GlobalPoint global_point_2 = surface_2.toGlobal(local_point_2);
          float global_x_2 = global_point_2.x();
          float global_y_2 = global_point_2.y();

	  auto hist_map_local_x_correlation = hist_local_x_correlation.find(module_pair);
	  if (hist_map_local_x_correlation == hist_local_x_correlation.end()) continue;
	  hist_map_local_x_correlation->second->Fill(local_x_1,local_x_2);
	  
	  auto hist_map_local_y_correlation = hist_local_y_correlation.find(module_pair);
	  if (hist_map_local_y_correlation == hist_local_y_correlation.end()) continue;
	  hist_map_local_y_correlation->second->Fill(local_y_1,local_y_2);
	  
	  auto hist_map_local_delta_x = hist_local_delta_x.find(module_pair);
	  if (hist_map_local_delta_x == hist_local_delta_x.end()) continue;
	  hist_map_local_delta_x->second->Fill(local_x_1 - local_x_2);
	  
	  auto hist_map_local_delta_y = hist_local_delta_y.find(module_pair);
	  if (hist_map_local_delta_y == hist_local_delta_y.end()) continue;
	  hist_map_local_delta_y->second->Fill(local_y_1 - local_y_2);
	  
	  auto hist_map_local_delta_x_vs_delta_y = hist_local_delta_x_vs_delta_y.find(module_pair);
	  if (hist_map_local_delta_x_vs_delta_y == hist_local_delta_x_vs_delta_y.end()) continue;
	  hist_map_local_delta_x_vs_delta_y->second->Fill(local_x_1 - local_x_2, local_y_1 - local_y_2);
	  
	  auto hist_map_global_x_correlation = hist_global_x_correlation.find(module_pair);
	  if (hist_map_global_x_correlation == hist_global_x_correlation.end()) continue;
	  hist_map_global_x_correlation->second->Fill(global_x_1,global_x_2);
	  
	  auto hist_map_global_y_correlation = hist_global_y_correlation.find(module_pair);
	  if (hist_map_global_y_correlation == hist_global_y_correlation.end()) continue;
	  hist_map_global_y_correlation->second->Fill(global_y_1,global_y_2);
	  
	  auto hist_map_global_delta_x = hist_global_delta_x.find(module_pair);
	  if (hist_map_global_delta_x == hist_global_delta_x.end()) continue;
	  hist_map_global_delta_x->second->Fill(global_x_1 - global_x_2);
	  
	  auto hist_map_global_delta_y = hist_global_delta_y.find(module_pair);
	  if (hist_map_global_delta_y == hist_global_delta_y.end()) continue;
	  hist_map_global_delta_y->second->Fill(global_y_1 - global_y_2);
	  
	  auto hist_map_global_delta_x_vs_delta_y = hist_global_delta_x_vs_delta_y.find(module_pair);
	  if (hist_map_global_delta_x_vs_delta_y == hist_global_delta_x_vs_delta_y.end()) continue;
	  hist_map_global_delta_x_vs_delta_y->second->Fill(global_x_1 - global_x_2, global_y_1 - global_y_2);

        }//end for jcluster
      }//end for icluster
    }//end for first jdet
  }//end for first idet

  
  //define iterations (in a map) to count the number of cluster per module in the event
  std::map<int, int> n_clusters_per_module;  
  int n_clusters_total = 0;

  for (std::map<int, TString>::iterator iter=det_id_to_module_name.begin(); iter!=det_id_to_module_name.end(); iter++) n_clusters_per_module.insert(std::pair<int, int>(iter->first, 0));
  
  for( edmNew::DetSetVector<SiPixelCluster>::const_iterator iDSV=pixelClusters->begin(); iDSV!=pixelClusters->end(); iDSV++) 
  {

    edmNew::DetSet<SiPixelCluster>::const_iterator iDS_begin = iDSV->begin();
    edmNew::DetSet<SiPixelCluster>::const_iterator iDS_end = iDSV->end();

    // Surface of (detId) and use the surface to convert 2D to 3D      
    auto det_id = DetId(iDSV->detId());
    int n_cluster = 0;
    
    const PixelGeomDetUnit *pixel_geom = (const PixelGeomDetUnit*) tk_geom->idToDetUnit(det_id);
    LocalPoint local_point(-9999., -9999., -9999.);

    // Then loop on the clusters of the module
    for(edmNew::DetSet<SiPixelCluster>::const_iterator icluster=iDS_begin; icluster!=iDS_end; ++icluster) 
    {
      PixelClusterParameterEstimator::ReturnType parameters = cpe.getParameters(*icluster, *pixel_geom);
      local_point = std::get<0>(parameters);

      const Surface& surface = tracker->idToDet(det_id)->surface();
      GlobalPoint global_point = surface.toGlobal(local_point);

      // Apply MASK
      // Inverse x-y coordinates for CHROMini
      double local_y = icluster->x(); // barycenter x position
      double local_x = icluster->y();// barycenter y position
      
      //if ( (uint32_t) det_id==344724484 && local_x>=80 && 312<=local_y && local_y<364) continue;
      if ( (uint32_t) det_id==344724484 && local_y>=80 && 312<=local_x && local_x<364) continue; // inverse x-y

      double global_x=0, global_y=0, global_z=0;
      global_x = global_point.x();
      global_y = global_point.y();
      global_z = global_point.z();
      
      hist_dqm_cluster_local_x_position[det_id]->Fill(local_x);
      hist_dqm_cluster_local_y_position[det_id]->Fill(local_y);
      hist_dqm_cluster_local_2d_position[det_id]->Fill(local_x, local_y);
      
      hist_dqm_cluster_global_x_position[det_id]->Fill(global_x);
      hist_dqm_cluster_global_y_position[det_id]->Fill(global_y);
      hist_dqm_cluster_global_2d_position[det_id]->Fill(global_x, global_y);
      
      hist_dqm_cluster_charge[det_id]->Fill(icluster->charge());
      hist_dqm_cluster_size[det_id]->Fill(icluster->size());
      //sizeX()
      //sizeY()

      n_clusters_per_module[det_id]++;
      n_clusters_total++;

	
      // Let's fill in the tree
      tree_run_number = event_id.run();
      tree_lumi_section = event_id.luminosityBlock();
      tree_event = event_id.event();
      tree_det_id = det_id;
      tree_module_name = det_id_to_module_name[det_id];
      tree_n_clusters = n_cluster++;	
      tree_cluster_charge = icluster->charge();	
      tree_cluster_size = icluster->size();
      tree_cluster_local_x_position = local_x;
      tree_cluster_local_y_position = local_y;
      tree_cluster_global_x_position = global_x;
      tree_cluster_global_y_position = global_y;
      tree_cluster_global_z_position = global_z;
      cluster_tree->Fill();
  

     } //end for cluster of the first detector
    } //end for first detector
    
  hist_n_clusters_tot->Fill(n_clusters_total);
  
  //fill the number of cluster in the event per module
  for (std::map<int, int>::iterator iter=n_clusters_per_module.begin(); iter!=n_clusters_per_module.end(); iter++) hist_dqm_n_clusters[iter->first]->Fill(iter->second);
  
  //---------------------------------
  //loop on hits
  //---------------------------------
  /*
  for( edmNew::DetSetVector<SiPixelRecHit>::const_iterator DSViter=pixelhits->begin(); DSViter!=pixelhits->end(); DSViter++   ) {
      
    edmNew::DetSet<SiPixelRecHit>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiPixelRecHit>::const_iterator end  =DSViter->end();
    
    for(edmNew::DetSet<SiPixelRecHit>::const_iterator iter=begin;iter!=end;++iter) {

      // Here we do something with the hits.
         
    }//end for DetSet
  }//end for DetSetVector
  
  */
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
