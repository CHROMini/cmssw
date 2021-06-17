{

  FILE * fout1 = fopen("noisy_list.txt","w");

  TChain *super = new TChain("DQMData/run100000/cluster3DTree");
  super->Add("PixelTelescope_BeamData_DQM_caro.root");

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

/*
  float xmin=-1.72;
  float xmax=1.65;
  float ymin=-3.28;
  float ymax=3.18;
*/
  float xmin=0;
  float xmax=160;
  float ymin=0;
  float ymax=420;
  int nbinx=160;
  int nbiny=420;
  std::vector<TH2F*> hdet;

  for (int i=0; i<16; i++) {
    TString name = Form("hdet%d",i);
    TH2F *hdettemp= new TH2F(name,name,nbinx, xmin,xmax,nbiny, ymin,ymax);
    hdet.push_back(hdettemp);
  }

  super->Draw("ybary:xbary >>hdet0 ","detId==344987652");
  super->Draw("ybary:xbary >>hdet1 ","detId==344986628");
  super->Draw("ybary:xbary >>hdet2 ","detId==344725508");
  super->Draw("ybary:xbary >>hdet3 ","detId==344724484");

  super->Draw("ybary:xbary >>hdet4 ","detId==344463364");
  super->Draw("ybary:xbary >>hdet5 ","detId==344462340");
  super->Draw("ybary:xbary >>hdet6 ","detId==344201220");
  super->Draw("ybary:xbary >>hdet7 ","detId==344200196");

  super->Draw("ybary:xbary >>hdet8 ","detId==352588804");
  super->Draw("ybary:xbary >>hdet9 ","detId==352589828");
  super->Draw("ybary:xbary >>hdet10 ","detId==352850948");
  super->Draw("ybary:xbary >>hdet11 ","detId==352851972");
  super->Draw("ybary:xbary >>hdet12 ","detId==353113092");
  super->Draw("ybary:xbary >>hdet13 ","detId==353114116");
  super->Draw("ybary:xbary >>hdet14 ","detId==353375236");
  super->Draw("ybary:xbary >>hdet15 ","detId==353376260");

/*
  std::vector<TCanvas*> can;
  for (int i=0; i<16; i++) {
    TString name = Form("c_%d",i);
    TCanvas *c1=new TCanvas(name,name,10,32,782,552);
  }  
*/

//  hdet[1]->Draw("colz");

  for (int i=0; i<16; i++) {
//   can[i]->cd();
//   hdet[i]->Draw("colz");
    float sum=0;
    float nbinpos=0;
    for (int ix=0;ix<nbinx+1;ix++) {
      for (int iy=0;iy<nbiny+1;iy++) {
        sum+=hdet[i]->GetBinContent(ix+1,iy+1);
        if (hdet[i]->GetBinContent(ix+1,iy+1)>0) nbinpos++;
      }
    }
    sum/=nbinpos;
    cout << " hdet " << i << " mean entry per bin " << sum << endl;
    for (int ix=0;ix<nbinx+1;ix++) {
      for (int iy=0;iy<nbiny+1;iy++) {
        float bincont=hdet[i]->GetBinContent(ix+1,iy+1);
        if (bincont>10*sum) {
        cout << " hdet " << i << " x " << ix << " y "<< iy << " bincont " << bincont << endl; 
        int detval=0;
        if (i==0) detval=344987652;
        if (i==1) detval=344986628;
        if (i==2) detval=344725508;
        if (i==3) detval=344724484;
        if (i==4) detval=344463364;
        if (i==5) detval=344462340;
        if (i==6) detval=344201220;
        if (i==7) detval=344200196;
        if (i==8) detval=352588804;
        if (i==9) detval=352589828;
        if (i==10) detval=352850948;
        if (i==11) detval=352851972;
        if (i==12) detval=353113092;
        if (i==13) detval=353114116;
        if (i==14) detval=353375236;
        if (i==15) detval=353376260;
        fprintf(fout1,"%d       %d       %d\n",detval,ix,iy);
        }
      } 
    }
  }
  fclose(fout1);

//   hdet[8]->Draw("colz");

/*
  c1->SaveAs("c1.pdf");
  c2->SaveAs("c2.pdf");
  c3->SaveAs("c3.pdf");
  c4->SaveAs("c4.pdf");
  c5->SaveAs("c5.pdf");
  c6->SaveAs("c6.pdf");
  c7->SaveAs("c7.pdf");
  c8->SaveAs("c8.pdf");
*/

}
