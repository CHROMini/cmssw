#define coarse_alignment_cxx
#include "coarse_alignment.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

using namespace std;

void coarse_alignment::LoadNoisy()
{
   ifstream file_noisy;
   file_noisy.open ("noisy_list.txt");
   inoisy=0;
   if (!file_noisy) { cerr << "cannot open file  noisy_list.txt " << endl; }
   else
   {
      while (!file_noisy.eof () && inoisy<noisy_max) {
       file_noisy >> noisy_det[inoisy] >> noisy_barx[inoisy] >> noisy_bary[inoisy];
       inoisy++;
      }
   }
   //file_noisy.close ();
   cout << " file_noisy : " << inoisy << " entries " << endl;
   if (inoisy>0) cout << " example detId "<< noisy_det[0] << " x "<< noisy_barx[0] << " y " << noisy_bary[0] <<  endl;
   file_noisy.close ();
}
void coarse_alignment::Loop()
{
//   In a ROOT session, you can do:
//      root> .L coarse_alignment.C
//      root> coarse_alignment t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   cout << "in Loop" << endl;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   TFile *myfile=new TFile("test_coarse.root",      "recreate");
   TH1F* h_dx_x2pos_12 = new TH1F("h_dx_x2pos_12","h_dx_x2pos_12",100,-0.2,0.2);
   TH1F* h_dy_x2pos_12 = new TH1F("h_dy_x2pos_12","h_dy_x2pos_12",100,-0.2,0.2);
   TH1F* h_dx_x2neg_12 = new TH1F("h_dx_x2neg_12","h_dx_x2neg_12",100,-0.2,0.2);
   TH1F* h_dy_x2neg_12 = new TH1F("h_dy_x2neg_12","h_dy_x2neg_12",100,-0.2,0.2);

   TH1F* h_dx_11 = new TH1F("h_dx_11","h_dx_11",250,-0.5,0.5);
   TH1F* h_dy_11 = new TH1F("h_dy_11","h_dy_11",250,-0.5,0.5);
   TH1F* h_dx_22 = new TH1F("h_dx_22","h_dx_22",250,-0.5,0.5);
   TH1F* h_dy_22 = new TH1F("h_dy_22","h_dy_22",250,-0.5,0.5);

   TH1F* h_dx_x2pos_23 = new TH1F("h_dx_x2pos_23","h_dx_x2pos_23",100,-0.2,0.2);
   TH1F* h_dy_x2pos_23 = new TH1F("h_dy_x2pos_23","h_dy_x2pos_23",100,-0.2,0.2);
   TH1F* h_dx_x2neg_23 = new TH1F("h_dx_x2neg_23","h_dx_x2neg_23",100,-0.2,0.2);
   TH1F* h_dy_x2neg_23 = new TH1F("h_dy_x2neg_23","h_dy_x2neg_23",100,-0.2,0.2);

   TH1F* h_dx_x2pos_34 = new TH1F("h_dx_x2pos_34","h_dx_x2pos_34",100,-0.2,0.2);
   TH1F* h_dy_x2pos_34 = new TH1F("h_dy_x2pos_34","h_dy_x2pos_34",100,-0.2,0.2);
   TH1F* h_dx_x2neg_34 = new TH1F("h_dx_x2neg_34","h_dx_x2neg_34",100,-0.2,0.2);
   TH1F* h_dy_x2neg_34 = new TH1F("h_dy_x2neg_34","h_dy_x2neg_34",100,-0.2,0.2);


   TH1F* h_dx_x2pos_45 = new TH1F("h_dx_x2pos_45","h_dx_x2pos_45",100,-0.2,0.2);
   TH1F* h_dy_x2pos_45 = new TH1F("h_dy_x2pos_45","h_dy_x2pos_45",100,-0.2,0.2);
   TH1F* h_dx_x2pos_45prim = new TH1F("h_dx_x2pos_45prim","h_dx_x2pos_45prim",500,-1.,1.);
   TH1F* h_dy_x2pos_45prim = new TH1F("h_dy_x2pos_45prim","h_dy_x2pos_45prim",500,-1.,1.);
   TH1F* h_dx_x2neg_45 = new TH1F("h_dx_x2neg_45","h_dx_x2neg_45",100,-0.2,0.2);
   TH1F* h_dy_x2neg_45 = new TH1F("h_dy_x2neg_45","h_dy_x2neg_45",100,-0.2,0.2);
//   TH1F* h_dx_x2neg_45 = new TH1F("h_dx_x2neg_45","h_dx_x2neg_45",500,-1.0,1.0);
//   TH1F* h_dy_x2neg_45 = new TH1F("h_dy_x2neg_45","h_dy_x2neg_45",500,-1.0,1.0);

   TH1F* h_dx_x2neg_46 = new TH1F("h_dx_x2neg_46","h_dx_x2neg_46",100,-0.2,0.2);
   TH1F* h_dy_x2neg_46 = new TH1F("h_dy_x2neg_46","h_dy_x2neg_46",100,-0.2,0.2);
   TH1F* h_dx_x2pos_46 = new TH1F("h_dx_x2pos_46","h_dx_x2pos_46",100,-0.2,0.2);
   TH1F* h_dy_x2pos_46 = new TH1F("h_dy_x2pos_46","h_dy_x2pos_46",100,-0.2,0.2);
//   TH1F* h_dx_x2neg_46 = new TH1F("h_dx_x2neg_46","h_dx_x2neg_46",500,-1.0,1.0);
//   TH1F* h_dy_x2neg_46 = new TH1F("h_dy_x2neg_46","h_dy_x2neg_46",500,-1.0,1.0);
   TH1F* h_dx_x2pos_56 = new TH1F("h_dx_x2pos_56","h_dx_x2pos_56",100,-0.2,0.2);
   TH1F* h_dy_x2pos_56 = new TH1F("h_dy_x2pos_56","h_dy_x2pos_56",100,-0.2,0.2);
   TH1F* h_dx_x2neg_56 = new TH1F("h_dx_x2neg_56","h_dx_x2neg_56",100,-0.2,0.2);
   TH1F* h_dy_x2neg_56 = new TH1F("h_dy_x2neg_56","h_dy_x2neg_56",100,-0.2,0.2);
   TH1F* h_dx_x2pos_67 = new TH1F("h_dx_x2pos_67","h_dx_x2pos_67",100,-0.2,0.2);
   TH1F* h_dy_x2pos_67 = new TH1F("h_dy_x2pos_67","h_dy_x2pos_67",100,-0.2,0.2);
   TH1F* h_dx_x2neg_67 = new TH1F("h_dx_x2neg_67","h_dx_x2neg_67",100,-0.2,0.2);
   TH1F* h_dy_x2neg_67 = new TH1F("h_dy_x2neg_67","h_dy_x2neg_67",100,-0.2,0.2);

   TH1F* h_dx_x2pos_47 = new TH1F("h_dx_x2pos_47","h_dx_x2pos_47",100,-0.2,0.2);
   TH1F* h_dy_x2pos_47 = new TH1F("h_dy_x2pos_47","h_dy_x2pos_47",100,-0.2,0.2);
   TH1F* h_dx_x2neg_47 = new TH1F("h_dx_x2neg_47","h_dx_x2neg_47",100,-0.2,0.2);
   TH1F* h_dy_x2neg_47 = new TH1F("h_dy_x2neg_47","h_dy_x2neg_47",100,-0.2,0.2);

   TH1F* h_dx_x2pos_78 = new TH1F("h_dx_x2pos_78","h_dx_x2pos_78",100,-0.2,0.2);
   TH1F* h_dy_x2pos_78 = new TH1F("h_dy_x2pos_78","h_dy_x2pos_78",100,-0.2,0.2);
   TH1F* h_dx_x2neg_78 = new TH1F("h_dx_x2neg_78","h_dx_x2neg_78",100,-0.2,0.2);
   TH1F* h_dy_x2neg_78 = new TH1F("h_dy_x2neg_78","h_dy_x2neg_78",100,-0.2,0.2);

   TH1F* h_dx_x2pos_48 = new TH1F("h_dx_x2pos_48","h_dx_x2pos_48",100,-0.2,0.2);
   TH1F* h_dy_x2pos_48 = new TH1F("h_dy_x2pos_48","h_dy_x2pos_48",100,-0.2,0.2);
   TH1F* h_dx_x2neg_48 = new TH1F("h_dx_x2neg_48","h_dx_x2neg_48",100,-0.2,0.2);
   TH1F* h_dy_x2neg_48 = new TH1F("h_dy_x2neg_48","h_dy_x2neg_48",100,-0.2,0.2);
   TH1F* h_dx_x2neg_58 = new TH1F("h_dx_x2neg_58","h_dx_x2neg_58",100,-0.2,0.2);
   TH1F* h_dy_x2neg_58 = new TH1F("h_dy_x2neg_58","h_dy_x2neg_58",100,-0.2,0.2);
   TH1F* h_dx_x2pos_68 = new TH1F("h_dx_x2pos_68","h_dx_x2pos_68",100,-0.2,0.2);
   TH1F* h_dy_x2pos_68 = new TH1F("h_dy_x2pos_68","h_dy_x2pos_68",100,-0.2,0.2);

   int nmax = 300000;
   double xclus[nmax];  
   double yclus[nmax];  
   int detclus[nmax];  
   int evlist[nmax];
   int index=0;  


   cout << "nentries " << nentries << endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      bool noisyflag=false;

      for (int in=0; in<inoisy; in++) {
         if (detId==noisy_det[in]) {
            if (int(xbary)==noisy_barx[in] && int(ybary)==noisy_bary[in]) {
                   noisyflag=true;
//                   cout << " find noisy channel "  << detId << " " << xbary << " "  << ybary << endl;
            }
         }
      }

      if (index<nmax && !noisyflag) {
        evlist[index]=event;
        detclus[index]=detId;
        xclus[index]=x;
        yclus[index]=y;
        if (detId==353114116) {
          xclus[index]-=-2.56917e-02;
          yclus[index]-=-2.38370e-02;
        }
        else if (detId==353113092) {
          xclus[index]-=-2.02377e-02;
          yclus[index]-=-2.70933e-02;
        }
        else if (detId==352851972) {
          xclus[index]-=2.47805e-02;
          yclus[index]-=-6.97738e-03;
        }
        else if (detId==352850948) {
          xclus[index]-=1.85113e-02;
          yclus[index]-=-9.17572e-03;
        }
        else if (detId==352589828) {
          xclus[index]-=5.50173e-02;
          yclus[index]-=-3.05621e-02;
        }
        else if (detId==352588804) {
          xclus[index]-=6.73395e-02;
          yclus[index]-=-3.06608e-02;
        }
        else if (detId==344200196) {
           xclus[index]-=-2.77220e-02;
           yclus[index]-=-9.31173e-02;
        }
        else if (detId==344201220) {
           xclus[index]-=-2.24462e-02;
           yclus[index]-=-8.72844e-02;
        }
        else if (detId==344462340) {
          xclus[index]-=7.14737e-03;
          yclus[index]-=-1.09290e-01;
        }
        else if (detId==344724484) {
          xclus[index]-=-7.60673e-02;
          yclus[index]-=-1.20672e-01;
        }
        else if (detId==344986628) {
          xclus[index]-=7.71767e-02;
          yclus[index]-=-1.20759e-01;
        }
        else if (detId==344987652) {
          xclus[index]-=9.50454e-02;
          yclus[index]-=-1.33751e-01;
        }

        index++;
      }
   }
   cout << "number of clusters not noisy " << index << endl;


   for (int j=0; j<index; j++) {
      if (detclus[j]==353376260 || detclus[j]==353375236) {
       // layer 1
       for (int k=0; k<index; k++) {
         if (evlist[k]==evlist[j]) {
          if (detclus[j]==353376260 && detclus[k]==353114116) {
             h_dx_x2pos_12->Fill(xclus[k]-xclus[j]);
             h_dy_x2pos_12->Fill(yclus[k]-yclus[j]);
          }

          if (detclus[j]==353375236 && detclus[k]==353113092) {
             h_dx_x2neg_12->Fill(xclus[k]-xclus[j]);
             h_dy_x2neg_12->Fill(yclus[k]-yclus[j]);
          }  
          if (detclus[j]==353376260 && detclus[k]==353375236) {
             h_dx_11->Fill(xclus[k]-xclus[j]);
             h_dy_11->Fill(yclus[k]-yclus[j]);
          }
            
         }
       }
      }
      if (detclus[j]==353114116 || detclus[j]==353113092) {
       // layer 2
       for (int k=0; k<index; k++) {
	 if (evlist[k]==evlist[j]) {
	  if (detclus[j]==353114116 && detclus[k]==352851972) {
	     h_dx_x2pos_23->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_23->Fill(yclus[k]-yclus[j]);
	  }

	  if (detclus[j]==353113092 && detclus[k]==352850948) {
	     h_dx_x2neg_23->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_23->Fill(yclus[k]-yclus[j]);
	  }  
          if (detclus[j]==353114116 && detclus[k]==353113092) {
             h_dx_22->Fill(xclus[k]-xclus[j]);
             h_dy_22->Fill(yclus[k]-yclus[j]);
          }
	    
	 }
       }
      }     
      if (detclus[j]==352851972 || detclus[j]==352850948) {
       // layer 3
       for (int k=0; k<index; k++) {
	 if (evlist[k]==evlist[j]) {
	  if (detclus[j]==352851972 && detclus[k]==352589828) {
	     h_dx_x2pos_34->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_34->Fill(yclus[k]-yclus[j]);
	  }
	  if (detclus[j]==352850948 && detclus[k]==352588804) {
	     h_dx_x2neg_34->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_34->Fill(yclus[k]-yclus[j]);
	  }  
	    
	 }
       }
      }     
      if (detclus[j]==352589828 || detclus[j]==352588804) {
       // layer 4
       for (int k=0; k<index; k++) {
	 if (evlist[k]==evlist[j]) {
//	  if (detclus[j]==352589828 && detclus[k]==344725508) {
	  if (detclus[j]==352589828 && detclus[k]==344200196) {
	     h_dx_x2pos_45->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_45->Fill(yclus[k]-yclus[j]);
	  }
//	  if (detclus[j]==352588804 && detclus[k]==344724484) {
	  if (detclus[j]==352588804 && detclus[k]==344201220) {
	     h_dx_x2neg_45->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_45->Fill(yclus[k]-yclus[j]);
	  }  
/*
	  if (detclus[j]==352589828 && detclus[k]==344724484) {
	     h_dx_x2pos_45prim->Fill(-1.*xclus[k]-xclus[j]);
	     h_dy_x2pos_45prim->Fill(yclus[k]-yclus[j]);
	  }
*/
//	  if (detclus[j]==352588804 && detclus[k]==344462340) {
	  if (detclus[j]==352589828 && detclus[k]==344462340) {
	     h_dx_x2pos_46->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_46->Fill(yclus[k]-yclus[j]);
	  }  
	  if (detclus[j]==352588804 && detclus[k]==344463364) {
	     h_dx_x2neg_46->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_46->Fill(yclus[k]-yclus[j]);
	  }  
/*
	  if (detclus[j]==352589828 && detclus[k]==344462340) {
	     h_dx_x2pos_46prim->Fill(-1.*xclus[k]-xclus[j]);
	     h_dy_x2pos_46prim->Fill(yclus[k]-yclus[j]);
	  }
*/
//	  if (detclus[j]==352589828 && detclus[k]==344200196) {
	  if (detclus[j]==352589828 && detclus[k]==344724484) {
	     h_dx_x2pos_47->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_47->Fill(yclus[k]-yclus[j]);
	  }
//	  if (detclus[j]==352588804 && detclus[k]==344201220) {
	  if (detclus[j]==352588804 && detclus[k]==344725508) {
	     h_dx_x2neg_47->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_47->Fill(yclus[k]-yclus[j]);
	  } 
	  if (detclus[j]==352589828 && detclus[k]==344986628) {
	     h_dx_x2pos_48->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_48->Fill(yclus[k]-yclus[j]);
	  }
	  if (detclus[j]==352588804 && detclus[k]==344987652) {
	     h_dx_x2neg_48->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_48->Fill(yclus[k]-yclus[j]);
	  }  
	    
	 }
       }
      }     
//      if (detclus[j]==344725508 || detclus[j]==344724484) {
      if (detclus[j]==344200196 || detclus[j]==344201220) {
       // layer 5
       for (int k=0; k<index; k++) {
	 if (evlist[k]==evlist[j]) {
//	  if (detclus[j]==344724484 && detclus[k]==344462340) {
	  if (detclus[j]==344201220 && detclus[k]==344463364) {
	     h_dx_x2neg_56->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_56->Fill(yclus[k]-yclus[j]);
	  }
	  if (detclus[j]==344200196 && detclus[k]==344462340) {
	     h_dx_x2pos_56->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_56->Fill(yclus[k]-yclus[j]);
	  }
	  if (detclus[j]==344201220 && detclus[k]==344987652) {
	     h_dx_x2neg_58->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_58->Fill(yclus[k]-yclus[j]);
	  }  
	 }
       }
      }     
      if (detclus[j]==344463364 || detclus[j]==344462340) {
       // layer 6
       for (int k=0; k<index; k++) {
	 if (evlist[k]==evlist[j]) {
//	  if (detclus[j]==344462340 && detclus[k]==344201220) {
	  if (detclus[j]==344463364 && detclus[k]==344725508) {
	     h_dx_x2neg_67->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_67->Fill(yclus[k]-yclus[j]);
	  }
	  if (detclus[j]==344462340 && detclus[k]==344724484) {
	     h_dx_x2pos_67->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_67->Fill(yclus[k]-yclus[j]);
	  }
	  if (detclus[j]==344462340 && detclus[k]==344986628) {
	     h_dx_x2pos_68->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_68->Fill(yclus[k]-yclus[j]);
	  }
	 }
       }
      }     
//      if (detclus[j]==344201220 || detclus[j]==344200196) {
      if (detclus[j]==344724484 || detclus[j]==344725508) {
       // layer 7
       for (int k=0; k<index; k++) {
	 if (evlist[k]==evlist[j]) {
//	  if (detclus[j]==344200196 && detclus[k]==344986628) {
	  if (detclus[j]==344724484 && detclus[k]==344986628) {
	     h_dx_x2pos_78->Fill(xclus[k]-xclus[j]);
	     h_dy_x2pos_78->Fill(yclus[k]-yclus[j]);
	  }
//	  if (detclus[j]==344201220 && detclus[k]==344987652) {
	  if (detclus[j]==344725508 && detclus[k]==344987652) {
	     h_dx_x2neg_78->Fill(xclus[k]-xclus[j]);
	     h_dy_x2neg_78->Fill(yclus[k]-yclus[j]);
	  }  
	 }
       }
      }     
        
   }



   myfile->cd();
   h_dx_x2pos_12->Write();
   h_dy_x2pos_12->Write();
   h_dx_x2neg_12->Write();
   h_dy_x2neg_12->Write();

   h_dx_11->Write();
   h_dy_11->Write();
   h_dx_22->Write();
   h_dy_22->Write();

   h_dx_x2pos_23->Write();
   h_dy_x2pos_23->Write();
   h_dx_x2neg_23->Write();
   h_dy_x2neg_23->Write();

   h_dx_x2pos_34->Write();
   h_dy_x2pos_34->Write();
   h_dx_x2neg_34->Write();
   h_dy_x2neg_34->Write();
   h_dx_x2pos_45->Write();
   h_dy_x2pos_45->Write();
   h_dx_x2pos_45prim->Write();
   h_dy_x2pos_45prim->Write();
   h_dx_x2neg_45->Write();
   h_dy_x2neg_45->Write();
   h_dx_x2pos_46->Write();
   h_dy_x2pos_46->Write();
   h_dx_x2neg_46->Write();
   h_dy_x2neg_46->Write();

   h_dx_x2pos_56->Write();
   h_dy_x2pos_56->Write();
   h_dx_x2neg_56->Write();
   h_dy_x2neg_56->Write();
   h_dx_x2pos_47->Write();
   h_dy_x2pos_47->Write();
   h_dx_x2neg_47->Write();
   h_dy_x2neg_47->Write();
   h_dx_x2pos_67->Write();
   h_dy_x2pos_67->Write();
   h_dx_x2neg_67->Write();
   h_dy_x2neg_67->Write();
   h_dx_x2pos_78->Write();
   h_dy_x2pos_78->Write();
   h_dx_x2neg_78->Write();
   h_dy_x2neg_78->Write();
   h_dx_x2pos_48->Write();
   h_dy_x2pos_48->Write();
   h_dx_x2neg_48->Write();
   h_dy_x2neg_48->Write();
   h_dx_x2neg_58->Write();
   h_dy_x2neg_58->Write();
   h_dx_x2pos_68->Write();
   h_dy_x2pos_68->Write();
   myfile->Close();
}
