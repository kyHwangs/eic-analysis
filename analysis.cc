#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"
#include "functions.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
  TString filename = argv[1];
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  std::ifstream in;
  int i_tmp;
  double ceren, scint;
  std::vector<std::pair<double,double>> fCalibs;
  in.open("<Ca & Eq Constants filename>.csv",std::ios::in);
  while (true) {
    in >> i_tmp >> ceren >> scint;
    if (!in.good()) break;
    fCalibs.push_back(std::make_pair(ceren, scint));
  }
  in.close();

  gStyle->SetOptFit(1);

  TH1F* tEdep = new TH1F("totEdep","Total Energy deposit;MeV;Evt",100,low*1000.,high*1000.);
  tEdep->Sumw2(); tEdep->SetLineColor(kRed); tEdep->SetLineWidth(2);
  TH1F* tE_C = new TH1F("E_C","Energy of Cerenkov ch.;GeV;Evt",100,low,high);
  tE_C->Sumw2(); tE_C->SetLineColor(kBlue); tE_C->SetLineWidth(2);
  TH1F* tE_S = new TH1F("E_S","Energy of Scintillation ch.;GeV;Evt",100,low,high);
  tE_S->Sumw2(); tE_S->SetLineColor(kRed); tE_S->SetLineWidth(2);
  TH1F* tE_SC = new TH1F("E_SC","E_{S}+E_{C};GeV;Evt",100,2.*low,2.*high);
  tE_SC->Sumw2(); tE_SC->SetLineColor(kBlack); tE_SC->SetLineWidth(2);
  TH1F* tE_DR = new TH1F("E_DR","Dual-readout corrected Energy;GeV;Evt",100,low,high);
  tE_DR->Sumw2(); tE_DR->SetLineColor(kBlack); tE_DR->SetLineWidth(2);
  TH1F* tP_leak = new TH1F("Pleak","Momentum leak;MeV;Evt",100,0.,1000.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2);
  TH1F* tP_leak_nu = new TH1F("Pleak_nu","Neutrino energy leak;MeV;Evt",100,0.,1000.*high);
  tP_leak_nu->Sumw2(); tP_leak_nu->SetLineWidth(2);

  TH1F* tT_C = new TH1F("time_C","Cerenkov time;ns;p.e.",150,10.,70.);
  tT_C->Sumw2(); tT_C->SetLineColor(kBlue); tT_C->SetLineWidth(2);
  TH1F* tT_S = new TH1F("time_S","Scint time;ns;p.e.",150,10.,70.);
  tT_S->Sumw2(); tT_S->SetLineColor(kRed); tT_S->SetLineWidth(2);
  TH1F* tWav_S = new TH1F("wavlen_S","Scint wavelength;nm;p.e.",120,300.,900.);
  tWav_S->Sumw2(); tWav_S->SetLineColor(kRed); tWav_S->SetLineWidth(2);
  TH1F* tWav_C = new TH1F("wavlen_C","Cerenkov wavelength;nm;p.e.",120,300.,900.);
  tWav_C->Sumw2(); tWav_C->SetLineColor(kBlue); tWav_C->SetLineWidth(2);
  TH1F* tNhit_S = new TH1F("nHits_S","Number of Scint p.e./SiPM;p.e.;n",200,0.,200.);
  tNhit_S->Sumw2(); tNhit_S->SetLineColor(kRed); tNhit_S->SetLineWidth(2);
  TH1F* tNhit_C = new TH1F("nHits_C","Number of Cerenkov p.e./SiPM;p.e.;n",50,0.,50.);
  tNhit_C->Sumw2(); tNhit_C->SetLineColor(kBlue); tNhit_C->SetLineWidth(2);

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename)+".root");
  drInterface->set("DRsim","DRsimEventData");
  unsigned int entries = drInterface->entries();

  std::vector<float> E_Ss,E_Cs;

  while (drInterface->numEvt() < entries) {
    if (drInterface->numEvt() % 100 == 0) printf("Analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);

    float Edep = 0.;
    for (auto edepItr = drEvt.Edeps.begin(); edepItr != drEvt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;
      Edep += edep.Edep;
    }
    tEdep->Fill(Edep);

    float Pleak = 0.;
    float Eleak_nu = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        Eleak_nu += leak4vec.P();
      } else {
        Pleak += leak4vec.P();
      }
    }
    tP_leak->Fill(Pleak);
    tP_leak_nu->Fill(Eleak_nu);

    float sEtmp = 0; float cEtmp = 0;
    for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {

      int fEta = 0;
      if ( tower->towerTheta.first >= 0 ) {
        fEta = tower->towerTheta.first;
      } else {
        fEta = std::abs( tower->towerTheta.first + 1 );
      }

      for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {
        if ( DRsimInterface::IsCerenkov(sipm->x,sipm->y) ) {
          tNhit_C->Fill(sipm->count);
          cEtmp += sipm->count * fCalibs.at(fEta).first;

          for (const auto timepair : sipm->timeStruct) {
            tT_C->Fill(timepair.first.first+0.05,timepair.second);
          }
          for (const auto wavpair : sipm->wavlenSpectrum) {
            tWav_C->Fill(wavpair.first.first,wavpair.second);
          }
        } else {
          tNhit_S->Fill(sipm->count);
          sEtmp += sipm->count * fCalibs.at(fEta).second;

          for (const auto timepair : sipm->timeStruct) {
            tT_S->Fill(timepair.first.first+0.05,timepair.second);
          }
          for (const auto wavpair : sipm->wavlenSpectrum) {
            tWav_S->Fill(wavpair.first.first,wavpair.second);
          }
        }
      }
    }

    E_Cs.push_back(cEtmp);
    E_Ss.push_back(sEtmp);

    tE_C->Fill(cEtmp);
    tE_S->Fill(sEtmp);
    tE_SC->Fill(cEtmp+sEtmp);
  } // event loop

  TCanvas* c = new TCanvas("c","");

  tEdep->Draw("Hist"); c->SaveAs(filename+"_Edep.png");

  c->cd();
  tE_S->SetTitle("");
  tE_S->Draw("Hist"); c->Update();
  TPaveStats* statsE_S = (TPaveStats*)c->GetPrimitive("stats");
  statsE_S->SetName("Scint");
  statsE_S->SetTextColor(kRed);
  statsE_S->SetY1NDC(.6); statsE_S->SetY2NDC(.8);

  tE_C->Draw("Hist&sames"); c->Update();
  TPaveStats* statsE_C = (TPaveStats*)c->GetPrimitive("stats");
  statsE_C->SetName("Cerenkov");
  statsE_C->SetTextColor(kBlue);
  statsE_C->SetY1NDC(.8); statsE_C->SetY2NDC(1.);
  c->SaveAs(filename+"_EcsHist.png");

  TF1* grE_C = new TF1("Cfit","gaus",low,high); grE_C->SetLineColor(kBlue);
  TF1* grE_S = new TF1("Sfit","gaus",low,high); grE_S->SetLineColor(kRed);
  tE_C->SetOption("p"); tE_C->Fit(grE_C,"R+&same");
  tE_S->SetOption("p"); tE_S->Fit(grE_S,"R+&same");

  c->cd();
  tE_S->SetTitle("");
  tE_S->Draw(""); c->Update();
  statsE_S->SetName("Scint");
  statsE_S->SetTextColor(kRed);
  statsE_S->SetX1NDC(.7);
  statsE_S->SetY1NDC(.4); statsE_S->SetY2NDC(.7);

  tE_C->Draw("sames"); c->Update();
  statsE_C->SetName("Cerenkov");
  statsE_C->SetTextColor(kBlue);
  statsE_C->SetX1NDC(.7);
  statsE_C->SetY1NDC(.7); statsE_C->SetY2NDC(1.);
  c->SaveAs(filename+"_Ecs.png");

  TF1* grE_SC = new TF1("S+Cfit","gaus",2.*low,2.*high); grE_SC->SetLineColor(kBlack);
  tE_SC->SetOption("p"); tE_SC->Fit(grE_SC,"R+&same");
  tE_SC->Draw(""); c->SaveAs(filename+"_Esum.png");

  c->SetLogy(1);
  tP_leak->Draw("Hist"); c->SaveAs(filename+"_Pleak.png");
  tP_leak_nu->Draw("Hist"); c->SaveAs(filename+"_Pleak_nu.png");
  c->SetLogy(0);

  TGraph* grSvsC = new TGraph(entries,&(E_Ss[0]),&(E_Cs[0]));
  grSvsC->SetTitle("SvsC;E_S;E_C");
  grSvsC->SetMarkerSize(0.5); grSvsC->SetMarkerStyle(20);
  grSvsC->GetXaxis()->SetLimits(0.,high);
  grSvsC->GetYaxis()->SetRangeUser(0.,high);
  grSvsC->SetMaximum(high);
  grSvsC->SetMinimum(0.);
  grSvsC->Draw("ap");
  c->SaveAs(filename+"_SvsC.png");

  tT_C->Draw("Hist"); c->SaveAs(filename+"_tC.png");
  tT_S->Draw("Hist"); c->SaveAs(filename+"_tS.png");
  tWav_C->Draw("Hist"); c->SaveAs(filename+"_wavC.png");
  tWav_S->Draw("Hist"); c->SaveAs(filename+"_wavS.png");
  tNhit_C->Draw("Hist"); c->SaveAs(filename+"_nhitC.png");
  tNhit_S->Draw("Hist"); c->SaveAs(filename+"_nhitS.png");
}
