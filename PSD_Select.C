#define PSD_Select_cxx
#include "PSD_Select.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>

using namespace std;

void PSD_Select::Loop(Long64_t entries)
{
//   In a ROOT session, you can do:
//      root> .L PSD_Select.C
//      root> PSD_Select t
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


  Int_t gamma_points=0;
  Int_t neutron_points=0;
  vector <unsigned int>* trace;

  if(fChannel==0)
    trace=trace_0;
  else
    trace=trace_3;

  g_gammas->Set(0);
  g_neutrons->Set(0);

  ofstream *out_n;
  ofstream *out_g;
  if(fWriteFiles){
    out_n=new ofstream(Form("Channel_%d_Neutron_trace_%d-%d.dat",fChannel,Amp_L,Amp_H));
    out_g=new ofstream(Form("Channel_%d_Gamma_trace_%d-%d.dat",fChannel,Amp_L,Amp_H));
  }
  Double_t voltage;
  fsampling_time=4e-9;
  Double_t trace_time_n=0;
  Double_t trace_time_g=0;
  //TCanvas *c2;
  //if(fVerbosity>0)
  //c2=new TCanvas();

   if (fChain == 0) return;

   Long64_t nentries=0;
   if(entries==-1)
     nentries = fChain->GetEntriesFast();
   else
     nentries=entries;
   Long64_t counter_g=0;
   Long64_t counter_n=0;
   theTrace->Set(0);
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     neutron=kFALSE;
     
     if(!fVerbosity){
     if(jentry%10000==0)
       std::cout<<jentry<<" traces analyzed: "<<counter_n<<"/"<<counter_g<<" neutron/gamma traces found\r"<<std::flush;
     }
     else
       std::cout<<jentry<<" traces analyzed so far"<<std::endl; 
     
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     //if(jentry>60000)
     //fVerbosity=2;
     //if (Cut(ientry) < 0 || trace->size()==0) continue;
     if ( trace->size()==0 ||trace->size()<400){ 
       if(fVerbosity>0)
	 std::cout<<jentry<<" event Skipped!"<<std::endl;
       continue;
     }
     if(fVerbosity>0)
       std::cout<<"PSD_Select::QDCcalc("<<jentry<<")"<<std::endl;
     QDCcalc(trace,0);
     if(false){
      for(UInt_t j=0;j<trace->size();j++){
	theTrace->SetPoint(j,j,trace->at(j));
	if(Cut(ientry) == 0){
	  theTrace->SetLineColor(1);
	  theTrace->SetMarkerColor(1);
	  g_neutrons->SetPoint(neutron_points,j*fsampling_time*1e9,trace->at(j));
	  if(true){
	    trace_time_n+=fsampling_time;
	    voltage=int(trace->at(j)-6600)/65355.;
	    if(fWriteFiles)
	      *out_n<<trace_time_n<<","<<voltage<<",";
	    //cout<<trace->at(j)<<","<<voltage<<",";
	    }
	  neutron_points++;
	  
	}
	else if(false){
	  theTrace->SetLineColor(2);
	  theTrace->SetMarkerColor(2);
	  g_gammas->SetPoint(gamma_points,j*fsampling_time*1e9,trace->at(j));
	  if(true){
	    trace_time_g+=fsampling_time;
	    voltage=int(trace->at(j)-6600)/65355.;
	    if(fWriteFiles)
	    *out_g<<trace_time_g<<","<<voltage<<",";
	  }
	  gamma_points++;
	}
      }
      if(Cut(ientry) == 0){
	counter_n++;
	neutron=kTRUE;
      }
      else	  counter_g++;	  
      theTrace->SetMarkerStyle(20);
      if(fVerbosity>0){
	//c2->WaitPrimitive();
	//c2->Update();
	theTrace->Draw("ALP");
	//std::cin.get();     
      }
     }//if false
      outTree->Fill();
      
      //theTrace.Fit(&fun,"RQN");
      //  if (Cut(ientry) == 0 ){
       // 	 hslope_1n->Fill(fun.GetParameter(1));
       // 	 hslope_2n->Fill(fun.GetParameter(3));
       // }
       // else{
       // 	 hslope_1g->Fill(fun.GetParameter(1));
       // 	 hslope_2g->Fill(fun.GetParameter(3));
       // }
   }
   TCanvas *c=new TCanvas();
   //TH1F* h=c->DrawFrame(0,Amp_H-5000,450,Amp_H+1000);
   if(g_gammas->GetN()>0)
     g_gammas->Draw("AP");
   g_gammas->GetXaxis()->SetTitle("time [ns]");
   g_gammas->GetYaxis()->SetTitle("Amplitude [ADC units]");
   if(g_neutrons->GetN()>0)
     g_neutrons->Draw("P");
   if(fWriteFiles){
     out_g->close();
     out_n->close();
   }
   outFile->Write();
   std::cout<<std::endl;

}
