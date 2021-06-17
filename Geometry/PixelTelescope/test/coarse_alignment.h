//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 10 12:07:32 2018 by ROOT version 6.10/09
// from TTree cluster3DTree/3D Cluster Tree
// found on file: PixelTelescope_BeamData_DQM_caro.root
//////////////////////////////////////////////////////////

#ifndef coarse_alignment_h
#define coarse_alignment_h
#define noisy_max 10000

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TString.h"

class coarse_alignment {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   int noisy_det[noisy_max];
   int noisy_barx[noisy_max];
   int noisy_bary[noisy_max];
   int inoisy;

   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           lumiSection;
   Int_t           event;
   Int_t           detId;
   TString         *modName;
   Int_t           cluster;
   Double_t        x;
   Double_t        y;
   Double_t        z;
   Double_t        x2;
   Double_t        y2;
   Double_t        z2;
   Double_t        xbary;
   Double_t        ybary;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_event;   //!
   TBranch        *b_detId;   //!
   TBranch        *b_modName;   //!
   TBranch        *b_cluster;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_x2;   //!
   TBranch        *b_y2;   //!
   TBranch        *b_z2;   //!
   TBranch        *b_xbary;   //!
   TBranch        *b_ybary;   //!

   coarse_alignment(TTree *tree=0);
   virtual ~coarse_alignment();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     LoadNoisy();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef coarse_alignment_cxx
coarse_alignment::coarse_alignment(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PixelTelescope_BeamData_DQM_caro.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("PixelTelescope_BeamData_DQM_caro.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("PixelTelescope_BeamData_DQM_caro.root:/DQMData/run100000");
      dir->GetObject("cluster3DTree",tree);

   }
   Init(tree);
}

coarse_alignment::~coarse_alignment()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t coarse_alignment::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t coarse_alignment::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void coarse_alignment::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   modName = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("detId", &detId, &b_detId);
   fChain->SetBranchAddress("modName", &modName, &b_modName);
   fChain->SetBranchAddress("cluster", &cluster, &b_cluster);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("x2", &x2, &b_x2);
   fChain->SetBranchAddress("y2", &y2, &b_y2);
   fChain->SetBranchAddress("z2", &z2, &b_z2);
   fChain->SetBranchAddress("xbary", &xbary, &b_xbary);
   fChain->SetBranchAddress("ybary", &ybary, &b_ybary);
   Notify();
}

Bool_t coarse_alignment::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void coarse_alignment::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t coarse_alignment::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef coarse_alignment_cxx
