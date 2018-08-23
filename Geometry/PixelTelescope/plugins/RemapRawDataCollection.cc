/** \file
 * Implementation of class RemapRawDataCollection
 *
 */

#include "RemapRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h" 
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Provenance/interface/ProcessHistory.h" 
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

using namespace edm;

RemapRawDataCollection::RemapRawDataCollection(const edm::ParameterSet& pset) {

  inputTags_ = pset.getParameter<std::vector<InputTag> >("RawCollectionList");
  verbose_ = pset.getUntrackedParameter<int>("verbose",0);

  inputTokens_.reserve(inputTags_.size());
  for(tag_iterator_t inputTag = inputTags_.begin(); inputTag != inputTags_.end(); ++inputTag ) {
    inputTokens_.push_back(consumes<FEDRawDataCollection>(*inputTag));
  }

  mappingList_ = pset.getUntrackedParameter<Parameters>("FEDMapping");
  for(auto it = mappingList_.begin(); it != mappingList_.end(); ++it) {
    int src=it->getParameter<int>("src");
    int dest=it->getParameter<int>("dest");
    auto check = dests_.emplace(dest);
    if (!check.second) {
      std::cout << "*** FED ID "<< dest <<" is already defined as destination!\n";
      continue;
    }
    if (mapping_.find(src)==mapping_.end()) {
      std::vector<int> v;
      mapping_[src]=v;
    }
    mapping_[src].push_back(dest);
//         FEDMapping = cms.untracked.VPSet(
//           cms.PSet(
//             src = cms.uint32(1),
//             dest = cms.uint32(1200),
//           ),
//         )
  }

  for (auto it = mapping_.begin(); it != mapping_.end(); ++it) {
    std::cout << it->first;
    for (size_t i=0; i<it->second.size(); ++i) {
      std::cout << "\t-> " << it->second[i] << std::endl;
    }
  }

  produces<FEDRawDataCollection>();
}

RemapRawDataCollection::~RemapRawDataCollection(){

}


void RemapRawDataCollection::produce(Event & e, const EventSetup& c){

 /// Get Data from all FEDs
 std::vector< Handle<FEDRawDataCollection> > rawData;
 rawData.reserve(inputTokens_.size());
 for(tok_iterator_t inputTok = inputTokens_.begin(); inputTok != inputTokens_.end(); ++inputTok ) {
   Handle<FEDRawDataCollection> input;
   if (e.getByToken(*inputTok,input)){
     rawData.push_back(input);
   }
   //else{     //skipping the inputtag requested. but this is a normal operation to bare data & MC. silent warning   }
 }

 auto producedData = std::make_unique<FEDRawDataCollection>();

 for (unsigned int i=0; i< rawData.size(); ++i ) { 

   const FEDRawDataCollection *rdc=rawData[i].product();

   if ( verbose_ > 0 ) {
     std::cout << "\nRAW collection #" << i+1 << std::endl;
     std::cout << "branch name = " << rawData[i].provenance()->branchName() << std::endl;
     std::cout << "process index = " << rawData[i].provenance()->productID().processIndex() << std::endl;
   }

   for ( int j=0; j< FEDNumbering::MAXFEDID; ++j ) {
     const FEDRawData & fedData = rdc->FEDData(j);
     if (fedData.size() == 0) continue;
     auto s=mapping_.find(j);
     if (dests_.find(j)==dests_.end() && s==mapping_.end()) {
       CopyFEDRawData(j, fedData, j, producedData->FEDData(j));
     }
     if (s==mapping_.end()) continue;
     for (size_t k=0; k<s->second.size(); k++) {
       CopyFEDRawData(j, fedData, s->second[k], producedData->FEDData(s->second[k]));
     }
   }
 }
 
 // Insert the new product in the event  
 e.put(std::move(producedData));

}


void RemapRawDataCollection::CopyFEDRawData(int src, const FEDRawData& fedData, int dest, FEDRawData& fedDataProd) {
  if(verbose_ > 1) std::cout << "Copying data from FED #" << src << std::endl;
  if ( fedDataProd.size() != 0 ) {
    if(verbose_ > 1) {
      std::cout << " More than one FEDRawDataCollection with data in FED ";
      std::cout << dest << " Skipping the 2nd\n";
    }
    return;
  } 
  fedDataProd.resize(fedData.size());
  unsigned char *dataProd=fedDataProd.data();
  const unsigned char *data=fedData.data();
  for ( unsigned int k=0; k<fedData.size(); ++k ) {
    dataProd[k]=data[k];
  }
}
