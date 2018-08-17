//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug  1 15:14:07 2018 by ROOT version 6.08/02
// from TTree T/Output Tree
// found on file: Output.root
//////////////////////////////////////////////////////////

#ifndef PSD_Select_h
#define PSD_Select_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCutG.h>
#include <TH1F.h>
// Header file for the classes stored in the TTree if any.
#include "vector"

#include <algorithm>
#include <utility>
#include <iterator> 
#include <iostream> 

class PSD_Select {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TFile *outFile;
   TTree *outTree;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       event;
   Double_t        phase[4];
   Double_t        Pmax[4];
   Double_t        Fmax[4];
   Double_t        time[4];
   Double_t        Pixietime[4];
   Double_t        ToF;
   Double_t        qdc[4];
   Double_t        leadqdc[4];
   Double_t        sbase[4];
   Double_t        abase[4];
   Double_t        thresh[4];
   Double_t        uPoint[4];
   Double_t        lPoint[4];
   Double_t        uThresh[4];
   Double_t        lThresh[4];
   Double_t        tailqdc[4];
   Double_t        ratio[4];
   Double_t        slope[4];
   Double_t        dpoint[4];
   Bool_t          k4fold;
   TCutG*          neutron_cut;
   TCutG*          gamma_cut;
   Bool_t          fWriteFiles;
   Int_t           fChannel;
   vector<unsigned int> *trace_0;
   vector<unsigned int> *trace_3;



   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_phase;   //!
   TBranch        *b_Pmax;   //!
   TBranch        *b_Fmax;   //!
   TBranch        *b_time;   //!
   TBranch        *b_Pixietime;   //!
   TBranch        *b_ToF;   //!
   TBranch        *b_qdc;   //!
   TBranch        *b_leadqdc;   //!
   TBranch        *b_sbase;   //!
   TBranch        *b_abase;   //!
   TBranch        *b_thresh;   //!
   TBranch        *b_uPoint;   //!
   TBranch        *b_lPoint;   //!
   TBranch        *b_uThresh;   //!
   TBranch        *b_lThresh;   //!
   TBranch        *b_tailqdc;   //!
   TBranch        *b_ratio;   //!
   TBranch        *b_slope;   //!
   TBranch        *b_dpoint;   //!
   TBranch        *b_k4fold;   //!
   TBranch        *b_trace_0;   //!
   TBranch        *b_trace_3;   //!

   TH1F *hslope_1n;
   TH1F *hslope_2n;
  
   TH1F *hslope_1g;
   TH1F *hslope_2g;

   //PSD variables
   Double_t qdc2[2];
   Double_t tailqdc2[2];
   Double_t ratio2[2];
   Double_t leadqdc2[2];
   Bool_t neutron;
   
   Int_t fDelay;
   Int_t fNbins;


   TGraph *g_gammas;
   TGraph *g_neutrons;
   TGraph *theTrace;

   Int_t Amp_L;
   Int_t Amp_H;
   Double_t fsampling_time=4e-9;

   void SetAmpLimits(Int_t low, Int_t high){Amp_L=low; Amp_H=high; }
   void EnableWrite(){fWriteFiles=kTRUE; }
   void DisableWrite(){fWriteFiles=kFALSE; }
   void SetCCDelay(Int_t CCdelay){fDelay = CCdelay;}   
   void SetNbins(Int_t value){fNbins = value;}   
   void SetVerboseLevel(Int_t value){fVerbosity = value;}   
   void SetChannelNumber(Int_t value){fChannel = value;}   


   PSD_Select(TTree *tree=0);
   virtual ~PSD_Select();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Long64_t entries = -1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   Int_t fVerbosity;
   //PSD functions
   virtual Double_t   CalcBaseline(vector <UInt_t> *dTrace);
   virtual void   QDCcalc(vector <UInt_t> *dTrace, int chan);
  


};

#endif

#ifdef PSD_Select_cxx
PSD_Select::PSD_Select(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Output.root");
      }
      f->GetObject("T",tree);

   }
   TFile *cutFile=new TFile("cuts.root");
   
   neutron_cut=(TCutG*)cutFile->Get("neutrons");
   //neutron_cut->SetVarY(ratio[1]);
   //neutron_cut->SetVarX(qdc[1]);
   gamma_cut=(TCutG*)cutFile->Get("gammas");
   //gamma_cut->SetVarY(ratio[1]);
   //gamma_cut->SetVarX(qdc[1]);

   theTrace=new TGraph();
   g_gammas=new TGraph();
   g_gammas->SetName("g_gammas");
   g_gammas->SetMarkerColor(2);
   g_neutrons=new TGraph();
   g_neutrons->SetName("g_neutrons");
   SetAmpLimits(10000,10100);
   SetChannelNumber(0);
   DisableWrite();
   outFile=new TFile("PSD_Out.root","RECREATE");
   outTree=new TTree("T2","PSD Output Tree");
   outTree->Branch("qdc2[2]",&qdc2,"qdc[2]/D");
   outTree->Branch("leadqdc2[2]",&leadqdc2,"leadqdc[2]/D");
   outTree->Branch("tailqdc2[2]",&tailqdc2,"tailqdc[2]/D");
   outTree->Branch("ratio2[2]",&ratio2,"ratio[2]/D");
   outTree->Branch("neutron",&neutron,"neutron/O");


   hslope_1n=new TH1F("h1","neutron slope1 distribution",200,0,.05);
   hslope_2n=new TH1F("h2","neutron slope2 distribution",200,0,.05);
   
   hslope_1g=new TH1F("h3","gamma slope1 distribution",200,0,.05);
   hslope_2g=new TH1F("h4","gamma slope2 distribution",200,0,.05);


   Init(tree);
}

PSD_Select::~PSD_Select()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PSD_Select::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PSD_Select::LoadTree(Long64_t entry)
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

void PSD_Select::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trace_0 = 0;
   trace_3 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("phase[4]", phase, &b_phase);
   fChain->SetBranchAddress("Pmax[4]", Pmax, &b_Pmax);
   fChain->SetBranchAddress("Fmax[4]", Fmax, &b_Fmax);
   fChain->SetBranchAddress("time[4]", time, &b_time);
   fChain->SetBranchAddress("Pixietime[4]", Pixietime, &b_Pixietime);
   fChain->SetBranchAddress("ToF", &ToF, &b_ToF);
   fChain->SetBranchAddress("qdc[4]", qdc, &b_qdc);
   fChain->SetBranchAddress("leadqdc[4]", leadqdc, &b_leadqdc);
   fChain->SetBranchAddress("sbase[4]", sbase, &b_sbase);
   fChain->SetBranchAddress("abase[4]", abase, &b_abase);
   fChain->SetBranchAddress("thresh[4]", thresh, &b_thresh);
   fChain->SetBranchAddress("uPoint[4]", uPoint, &b_uPoint);
   fChain->SetBranchAddress("lPoint[4]", lPoint, &b_lPoint);
   fChain->SetBranchAddress("uThresh[4]", uThresh, &b_uThresh);
   fChain->SetBranchAddress("lThresh[4]", lThresh, &b_lThresh);
   fChain->SetBranchAddress("tailqdc[4]", tailqdc, &b_tailqdc);
   fChain->SetBranchAddress("ratio[4]", ratio, &b_ratio);
   fChain->SetBranchAddress("slope[4]", slope, &b_slope);
   fChain->SetBranchAddress("dpoint[4]", dpoint, &b_dpoint);
   fChain->SetBranchAddress("k4fold", &k4fold, &b_k4fold);
   fChain->SetBranchAddress("trace_0", &trace_0, &b_trace_0);
   fChain->SetBranchAddress("trace_3", &trace_3, &b_trace_3);
   Notify();
   SetCCDelay(20);
   fNbins=150;
   fVerbosity=0;

}

Bool_t PSD_Select::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PSD_Select::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PSD_Select::Cut(Long64_t entry)
{
  Int_t val=-1;
  if(Pmax[fChannel]>Amp_L && Pmax[fChannel]<Amp_H ){
    if(neutron_cut->IsInside(qdc[1],ratio[1]))
      val=0;
    else if(gamma_cut->IsInside(qdc[1],ratio[1]))
      val=1;
  }
  return val;
}

Double_t PSD_Select::CalcBaseline(vector <UInt_t> *dTrace)
{

  vector <UInt_t>::iterator it;



  UInt_t max_position=0;
  UInt_t initialpos=0,finalpos;
  it=max_element(dTrace->begin(),dTrace->end());
  max_position=distance(dTrace->begin(),it);
  if(fVerbosity>1){
 std::cout<<"PSD_Select::CalcBaseline()->max_position->"<< max_position<<" "<<*it<<std::endl;
  }
  
  finalpos=max_position-50;

  Double_t Baseline=0;
  Int_t count=0;
  //if(int(finalpos)<0) return -1;
  if(fVerbosity>1)
    std::cout<<"PSD_Select::CalcBaseline"<< initialpos<<" "<<int(finalpos)<<std::endl;
  for (int i =initialpos; i < (int)finalpos; i++){
    Baseline+=dTrace->at(i);
    count++;
  } 

  Baseline/=count;

  return Baseline;
}


void PSD_Select::QDCcalc(vector <UInt_t> *dTrace, int chan){
 qdc2[chan]=0;
 tailqdc2[chan]=0;
 leadqdc2[chan]=0;
 vector <UInt_t>::iterator it;
 UInt_t max_position=0;
 it=max_element(dTrace->begin(),dTrace->end());
 max_position=distance(dTrace->begin(),it);
 Double_t cfdpos=max_position;
 if(cfdpos==0) return;
 if(fVerbosity>0)
   std::cout<<"PSD_Select::QDCcalc()->CalcBaseline()"<<std::endl;
 Double_t baseline=CalcBaseline(dTrace);
 int cfdposL = floor(cfdpos);
 int cfdposR = ceil(cfdpos);

 double tL = cfdpos-(double)cfdposL;
 double tR = cfdpos-(double)cfdposR;
 double partialT = ((double)dTrace->at(cfdposR)*tR+(double)dTrace->at(cfdposR-1))*tR;
 tailqdc2[chan] += (partialT+(double)dTrace->at(cfdposR))*tR/2;
 if(fVerbosity>0){
 std::cout<<"***************************************"<<std::endl;
 std::cout<<"baseline->"<< baseline<<std::endl;
 std::cout<<"max_position->"<< max_position<<std::endl;
 std::cout<<"cfdpos->"<<cfdposL<<" "<<cfdposR<<std::endl;
 std::cout<<"cfdpos+fDelay->"<<cfdpos+fDelay<<" +fNbins "<<cfdpos+fNbins<<std::endl;
 }
 if (cfdposR+fNbins>=dTrace->size()-1 || cfdposR<5) return;
 for (int i = cfdposR-5; i < cfdposR+fNbins; i++){
  if(((double)dTrace->at(i)+(double)dTrace->at(i+1))/2.0>baseline){
   qdc2[chan] += ((double)(dTrace->at(i)+dTrace->at(i+1))-baseline*2.0)*0.5; 
   if (i>cfdposR+fDelay) tailqdc2[chan] += ((double)(dTrace->at(i)+dTrace->at(i+1))-baseline*2.0)*0.5;
   if (i<cfdposR+4) leadqdc2[chan] += ((double)(dTrace->at(i)+dTrace->at(i+1))-baseline*2.0)*0.5;
   }
  else continue;
  }
 ratio2[chan] =  tailqdc2[chan]/qdc2[chan];

 if(fVerbosity>0){
 std::cout<<"qdc[]->"<<qdc2[chan]<<std::endl;
 std::cout<<"leadqdc[]->"<<leadqdc2[chan]<<std::endl;
 std::cout<<"tailqdc[]->"<<tailqdc2[chan]<<std::endl;

 std::cout<<"Ratio[]->"<<ratio2[chan]<<std::endl;
 }
// }

}



#endif // #ifdef PSD_Select_cxx
