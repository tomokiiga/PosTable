#include "TROOT.h"
#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
//#include <initializer_list>
using namespace std;


int table(){
    //gROOT->SetStyle("ATLAS");

    // Read tree and set branchs
    	 TFile* file = TFile::Open("./merge.root");
    TTree* reduced = dynamic_cast<TTree*>(file->Get("reduced"));
   
	
    // #include "./SetTree.icc"
   Float_t         Header[5];
   Float_t         NPETotal4us;
   Float_t         Ratio4us;
   Float_t         Position_wm[3];
    

  TBranch        *b_Header;   //!
  TBranch        *b_NPETotal4us;   //!
  TBranch        *b_Ratio4us;   //!
  TBranch        *b_Position_wm;   //!


   reduced->SetBranchAddress("Header", Header, &b_Header);
   reduced->SetBranchAddress("NPETotal4us", &NPETotal4us, &b_NPETotal4us);
   reduced->SetBranchAddress("Position_wm", Position_wm, &b_Position_wm);
   reduced->SetBranchAddress("Ratio4us", &Ratio4us, &b_Ratio4us);

    reduced->SetBranchStatus("*",0);
    reduced->SetBranchStatus("Header");
    reduced->SetBranchStatus("NPETotal4us");
    reduced->SetBranchStatus("Position_wm");
    reduced->SetBranchStatus("Ratio4us");

    Int_t Entries = reduced->GetEntries();
    cout<<"ok"<<endl;

     const double InitFactor = 0.0;
     const double TermFactor = 0.2;

     const double nSigma = 1.5; // Fitting range for position table

    // Set NPETotal4us range
    TH1D* HistNPE = new TH1D("HistNPE","HistNPE",1000,0,5000);
    for(int iEntry=Entries*InitFactor; iEntry<Entries*TermFactor; ++iEntry){
        reduced->GetEntry(iEntry);
        if(iEntry % 1000000 == 0) cout << iEntry << "/" << Entries << "(" << (double)iEntry/Entries*100. << "%)" << endl;
        HistNPE->Fill(NPETotal4us);
    }
    TF1* NPEGaus = new TF1("NPEGaus","gaus");
    const double DummyCenter = HistNPE->GetBinCenter(HistNPE->GetMaximumBin());
    HistNPE->Draw();
    HistNPE->Fit("NPEGaus","S","",DummyCenter*0.9,DummyCenter*1.1);
    const double RangeCenterNPE = NPEGaus->GetParameter(1);
    const double RangeSigmaNPE = 1.;
    const double RangeMinNPE = NPEGaus->GetParameter(1) - RangeSigmaNPE*NPEGaus->GetParameter(2);
    const double RangeMaxNPE = NPEGaus->GetParameter(1) + RangeSigmaNPE*NPEGaus->GetParameter(2);

    // Range of Ratio4us
    const double MinRatio = 0.14;
    const double MaxRatio = 0.18;

    // 1st Fit -- x,y,z-projection fitting
    TH1D* PosPro[3];
    for(int i=0;i<3;++i){PosPro[i] = new TH1D(Form("PosPro%01d",i),Form("PosPro%01d",i),1200,-600,600);}

    for(int iEntry=Entries*InitFactor; iEntry<Entries*TermFactor; ++iEntry){
        reduced->GetEntry(iEntry);
        if(iEntry % 1000000 == 0) cout << iEntry << "/" << Entries << "(" << (double)iEntry/Entries*100. << "%)" << endl;
        if(RangeMinNPE > NPETotal4us) continue;
        if(RangeMaxNPE < NPETotal4us) continue;
        if(MinRatio > Ratio4us) continue;
        if(MaxRatio < Ratio4us) continue;
        PosPro[0]->Fill(Position_wm[0]);
        PosPro[1]->Fill(Position_wm[1]);
        PosPro[2]->Fill(Position_wm[2]);
    }
    cout<<"ok2"<<endl;
       
    TF1* PosTestGaus = new TF1("PosTestGaus","gaus");
    vector<vector<double>> PosCenter = {{-360, -180, 0, 180, 330}, {-300, -100, 100, 300}, {-400, -250, -85, 110, 290, 430}};
    vector<vector<double>> PosSigma  = {{ 40, 40, 40, 40, 40},{ 40, 40, 40, 40},{ 40, 40, 40, 40, 40, 40}};

    for(int i=0;i<PosCenter.size();++i){
        for(int j=0;j<PosCenter.at(i).size();++j){
            PosPro[i]->Fit("PosTestGaus","SQ","",PosCenter.at(i).at(j)-nSigma*PosSigma.at(i).at(j),PosCenter.at(i).at(j)+nSigma*PosSigma.at(i).at(j));
            cout << i << " " << j << " " << PosTestGaus->GetChisquare()/PosTestGaus->GetNDF() << " " 
                << PosTestGaus->GetParameter(1) << " " << PosTestGaus->GetParameter(2) << endl;
            PosCenter.at(i).at(j) = PosTestGaus->GetParameter(1);
            PosSigma.at(i).at(j)  = PosTestGaus->GetParameter(2);
        }
    }
    TCanvas* C1 = new TCanvas();
    PosPro[0]->Draw();

    int CrystalTmp, OrderXTmp, OrderYTmp, OrderZTmp;
    vector<int> OrderX, OrderY, OrderZ;
    OrderX.push_back(0);
    OrderY.push_back(0);
    OrderZ.push_back(0);
    ifstream ifs("./CrystalList.txt");
    while(ifs >> CrystalTmp >> OrderXTmp >> OrderYTmp >> OrderZTmp){
        OrderX.push_back(OrderXTmp);
        OrderY.push_back(OrderYTmp);
        OrderZ.push_back(OrderZTmp);
    }
    // Output file
    ofstream ofs("./out.txt",ios::app);
    TCanvas* C2 = new TCanvas();

    // Final Fitting
    //int FittingCrystal = 59;
    TH1D* Pos[3];
    for(int i=0;i<3;++i){Pos[i] = new TH1D(Form("Pos%01d",i),Form("Pos%01d",i),1200,-600,600);}
    TF1* PosGaus = new TF1("PosGaus","gaus");
    PosGaus->SetLineColor(kRed);

    ofs << FittingCrystal << "  ";
    // X-axis (Dim 0)
    for(int iEntry=Entries*InitFactor; iEntry<Entries*TermFactor; ++iEntry){
        tree->GetEntry(iEntry);
        if(iEntry % 1000000 == 0) cout << iEntry << "/" << Entries << "(" << (double)iEntry/Entries*100. << "%)" << endl;
        if(RangeMinNPE > NPETotal4us) continue;
        if(RangeMaxNPE < NPETotal4us) continue;
        if(MinRatio > Ratio4us) continue;
        if(MaxRatio < Ratio4us) continue;

        //Y
        if(PosCenter.at(1).at(OrderY.at(FittingCrystal)) - nSigma*PosSigma.at(1).at(OrderY.at(FittingCrystal)) > Position_wm[1]) continue;
        if(PosCenter.at(1).at(OrderY.at(FittingCrystal)) + nSigma*PosSigma.at(1).at(OrderY.at(FittingCrystal)) < Position_wm[1]) continue;

        //Z
        if(PosCenter.at(2).at(OrderZ.at(FittingCrystal)) - nSigma*PosSigma.at(2).at(OrderZ.at(FittingCrystal)) > Position_wm[2]) continue;
        if(PosCenter.at(2).at(OrderZ.at(FittingCrystal)) + nSigma*PosSigma.at(2).at(OrderZ.at(FittingCrystal)) < Position_wm[2]) continue;

        //X
        Pos[0]->Fill(Position_wm[0]);
    }
    Pos[0]->Fit("PosGaus","S","",PosCenter.at(0).at(OrderX.at(FittingCrystal)) - nSigma*PosSigma.at(0).at(OrderX.at(FittingCrystal)),
                PosCenter.at(0).at(OrderX.at(FittingCrystal)) + nSigma*PosSigma.at(0).at(OrderX.at(FittingCrystal)));
    Pos[0]->Fit("PosGaus","S","",PosGaus->GetParameter(1)-nSigma*PosGaus->GetParameter(2),PosGaus->GetParameter(1)+nSigma*PosGaus->GetParameter(2));

    ofs << PosGaus->GetParameter(1) << "  ";
    ofs << PosGaus->GetParameter(2) << "  ";
    C2->SaveAs(Form("./figs/Cry%02d-X.pdf",FittingCrystal));


    // Y-axis (Dim 1)
    for(int iEntry=Entries*InitFactor; iEntry<Entries*TermFactor; ++iEntry){
        tree->GetEntry(iEntry);
        if(iEntry % 1000000 == 0) cout << iEntry << "/" << Entries << "(" << (double)iEntry/Entries*100. << "%)" << endl;
        if(RangeMinNPE > NPETotal4us) continue;
        if(RangeMaxNPE < NPETotal4us) continue;
        if(MinRatio > Ratio4us) continue;
        if(MaxRatio < Ratio4us) continue;

        //Z
        if(PosCenter.at(2).at(OrderZ.at(FittingCrystal)) - nSigma*PosSigma.at(2).at(OrderZ.at(FittingCrystal)) > Position_wm[2]) continue;
        if(PosCenter.at(2).at(OrderZ.at(FittingCrystal)) + nSigma*PosSigma.at(2).at(OrderZ.at(FittingCrystal)) < Position_wm[2]) continue;

        //X
        if(PosCenter.at(0).at(OrderX.at(FittingCrystal)) - nSigma*PosSigma.at(0).at(OrderX.at(FittingCrystal)) > Position_wm[0]) continue;
        if(PosCenter.at(0).at(OrderX.at(FittingCrystal)) + nSigma*PosSigma.at(0).at(OrderX.at(FittingCrystal)) < Position_wm[0]) continue;

        //Y
        Pos[1]->Fill(Position_wm[1]);
        
    }
    Pos[1]->Fit("PosGaus","S","",PosCenter.at(1).at(OrderY.at(FittingCrystal)) - nSigma*PosSigma.at(1).at(OrderY.at(FittingCrystal)),
                PosCenter.at(1).at(OrderY.at(FittingCrystal)) + nSigma*PosSigma.at(1).at(OrderY.at(FittingCrystal)));
    Pos[1]->Fit("PosGaus","S","",PosGaus->GetParameter(1)-nSigma*PosGaus->GetParameter(2),PosGaus->GetParameter(1)+nSigma*PosGaus->GetParameter(2));
    ofs << PosGaus->GetParameter(1) << "  ";
    ofs << PosGaus->GetParameter(2) << "  ";
    C2->SaveAs(Form("./figs/Cry%02d-Y.pdf",FittingCrystal));

    // Z-axis (Dim 2)
    for(int iEntry=Entries*InitFactor; iEntry<Entries*TermFactor; ++iEntry){
        tree->GetEntry(iEntry);
        if(iEntry % 1000000 == 0) cout << iEntry << "/" << Entries << "(" << (double)iEntry/Entries*100. << "%)" << endl;
        if(RangeMinNPE > NPETotal4us) continue;
        if(RangeMaxNPE < NPETotal4us) continue;
        if(MinRatio > Ratio4us) continue;
        if(MaxRatio < Ratio4us) continue;

        //X
        if(PosCenter.at(0).at(OrderX.at(FittingCrystal)) - nSigma*PosSigma.at(0).at(OrderX.at(FittingCrystal)) > Position_wm[0]) continue;
        if(PosCenter.at(0).at(OrderX.at(FittingCrystal)) + nSigma*PosSigma.at(0).at(OrderX.at(FittingCrystal)) < Position_wm[0]) continue;

        //Y
        if(PosCenter.at(1).at(OrderY.at(FittingCrystal)) - nSigma*PosSigma.at(1).at(OrderY.at(FittingCrystal)) > Position_wm[1]) continue;
        if(PosCenter.at(1).at(OrderY.at(FittingCrystal)) + nSigma*PosSigma.at(1).at(OrderY.at(FittingCrystal)) < Position_wm[1]) continue;

        //Z
        Pos[2]->Fill(Position_wm[2]);
        
    }
    Pos[2]->Fit("PosGaus","S","",PosCenter.at(2).at(OrderZ.at(FittingCrystal)) - nSigma*PosSigma.at(2).at(OrderZ.at(FittingCrystal)),
                PosCenter.at(2).at(OrderZ.at(FittingCrystal)) + nSigma*PosSigma.at(2).at(OrderZ.at(FittingCrystal)));
    Pos[2]->Fit("PosGaus","S","",PosGaus->GetParameter(1)-nSigma*PosGaus->GetParameter(2),PosGaus->GetParameter(1)+nSigma*PosGaus->GetParameter(2));
    ofs << PosGaus->GetParameter(1) << "  ";
    ofs << PosGaus->GetParameter(2) << "  ";
    ofs << endl;
    C2->SaveAs(Form("./figs/Cry%02d-Z.pdf",FittingCrystal));

    
return 0;
}
