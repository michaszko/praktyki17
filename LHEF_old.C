#define LHEF_cxx
#include "LHEF.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void LHEF::Loop()
{
  TFile *f = new TFile("Histogramy.root", "RECREATE");

  std::vector<Double_t> kwark1, kwark2;
  std::vector<Double_t> elektorn, mion, eneutrino, mneutrino;
  std::vector<Double_t> aelektorn, amion, aeneutrino, amneutrino;
  std::vector<Int_t> jet1, jet2;
  std::vector<Int_t> lepton;
  std::vector<TLorentzVector> elektorn_v, mion_v, eneutrino_v, mneutrino_v;
  std::vector<TLorentzVector> aelektorn_v, amion_v, aeneutrino_v, amneutrino_v;

  kwark1.clear(); kwark2.clear(); jet1.clear(); jet2.clear(); lepton.clear();
  elektorn.clear(); mion.clear(); eneutrino.clear(); mneutrino.clear();
  aelektorn.clear(); amion.clear(); aeneutrino.clear(); amneutrino.clear();


  TH1F * mjj = new TH1F("mjj", "mjj", 200, 0., 1000.);
  TH1F * mll = new TH1F("mll", "mll", 100, 0., 1000.);
  TH1F * eta_j1  = new TH1F("eta_j1", "eta_j1", 200, -5., 5.);
  TH1F * eta_j2 = new TH1F("eta_j2", "eta_j2", 200, -5., 5.);
  TH1F * Deltaeta_j1j2 = new TH1F("Deltaeta_j1j2", "Deltaeta_j1j2", 100, 0., 10.);
  TH1F * Deltaeta_ll = new TH1F("Deltaeta_ll", "Deltaeta_ll", 100, 0., 5.);
  TH1F * Deltaeta_ll_hf = new TH1F("Deltaeta_ll_hf", "Deltaeta_ll_hf", 100, 0., 50.);
  TH1F * Deltaphi_ww = new TH1F("Deltaphi_ww", "Deltaphi_ww", 100, -1000., 1000);
  TH1F * Deltaphi_ww_hf = new TH1F("Deltaphi_ww_hf", "Deltaphi_ww_hf", 100, 0., 50.);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      for (int i = 0; i < Particle_size; i++)
      {
        //if (Particle_Status[i] == -1) continue;

        if (Particle_PID[i] >= 0 && Particle_PID[i] <= 8 && Particle_Status[i] == 1 )
        {
          jet1.push_back(i);
        }

        if (Particle_PID[i] >= -8 && Particle_PID[i] <= 0 && Particle_Status[i] ==1 )
        {
          jet2.push_back(i);
        }

        if( abs(Particle_PID[i]) == 11 || abs(Particle_PID[i]) == 13 || abs(Particle_PID[i]) == 15 || abs(Particle_PID[i]) == 17)
          lepton.push_back(i);

        // else if ( abs(Particle_PID[i]) == 11 || abs(Particle_PID[i]) == 13 || abs(Particle_PID[i]) == 15 || abs(Particle_PID[i]) == 17 )
        //   mll->Fill(Particle_M[i]);

        if (Particle_PID[i] >= 1 && Particle_PID[i] <= 8 && Particle_Status[i] == 1)
        {
          eta_j1->Fill(Particle_Eta[i]);
          kwark1.push_back(Particle_Eta[i]);
        }

        else if (Particle_PID[i] >= -8 && Particle_PID[i] <= 1 && Particle_Status[i] == 1 )
        {
          eta_j2->Fill(Particle_Eta[i]);
          kwark2.push_back(Particle_Eta[i]);
        }

        TLorentzVector x;
        x.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

        if (Particle_PID[i] == 11)
        {
          elektorn.push_back(Particle_Eta[i]);
          elektorn_v.push_back(x);
        }
        else if (Particle_PID[i] == 12)
        {
          eneutrino.push_back(Particle_Eta[i]);
          eneutrino_v.push_back(x);
        }
        else if (Particle_PID[i] == 13)
        {
          mion.push_back(Particle_Eta[i]);
          mion_v.push_back(x);
        }
        else if (Particle_PID[i] == 14)
        {
          mneutrino.push_back(Particle_Eta[i]);
          mneutrino_v.push_back(x);
        }
        else if (Particle_PID[i] == -11)
        {
          aelektorn.push_back(Particle_Eta[i]);
          aelektorn_v.push_back(x);
        }
        else if (Particle_PID[i] == -12)
        {
          aeneutrino.push_back(Particle_Eta[i]);
          aeneutrino_v.push_back(x);
        }
        else if (Particle_PID[i] == -13)
        {
          amion.push_back(Particle_Eta[i]);
          amion_v.push_back(x);
        }
        else if (Particle_PID[i] == -14)
        {
          amneutrino.push_back(Particle_Eta[i]);
          amneutrino_v.push_back(x);
        }

      }

      for (int i = 0; i < kwark1.size(); i++)
      {
        for (int j = 0; j < kwark2.size(); j++)
        {
          Deltaeta_j1j2->Fill(abs(kwark1[i] - kwark2[j]));
        }
      }

      for(int i = 0; i < elektorn.size(); i++)
      {
        for (int j = 0; j < amion.size(); j++)
        {
          Deltaeta_ll->Fill(abs(elektorn[i] - amion[j]));
        }
      }

      for (int i = 0; i < jet1.size(); i++)
      {
          for (int j = 0; j < jet2.size(); j++)
          {
            TLorentzVector temp1, temp2, sum;
            temp1.SetPxPyPzE(Particle_Px[jet1[i]], Particle_Py[jet1[i]], Particle_Pz[jet1[i]], Particle_E[jet1[i]] );
            temp2.SetPxPyPzE(Particle_Px[jet2[j]], Particle_Py[jet2[j]], Particle_Pz[jet2[j]], Particle_E[jet2[j]] );

            sum = temp1 + temp2;

            mjj->Fill(sum.M());

          }
      }

      for (int i = 0; i < lepton.size(); i+=2)
      {
        TLorentzVector temp1, temp2, sum;
        temp1.SetPxPyPzE(Particle_Px[lepton[i]], Particle_Py[lepton[i]], Particle_Pz[lepton[i]], Particle_E[lepton[i]] );
        temp2.SetPxPyPzE(Particle_Px[lepton[i+1]], Particle_Py[lepton[i+1]], Particle_Pz[lepton[i+1]], Particle_E[lepton[i+1]] );

        sum = temp1 + temp2;

        mjj->Fill(sum.M());
     }

    // TLorentzVector temp1, temp2, sum;
    // temp1.SetPxPyPzE(Particle_Px[lepton[0]], Particle_Py[lepton[0]], Particle_Pz[lepton[0]], Particle_E[lepton[0]] );
    // temp2.SetPxPyPzE(Particle_Px[lepton[1]], Particle_Py[lepton[1]], Particle_Pz[lepton[1]], Particle_E[lepton[1]] );
    //
    // sum = temp1 + temp2;
    //
    // mjj->Fill(sum.M());

      while (!aeneutrino_v.empty())
      {
        TLorentzVector t1 = elektorn_v.pop_back();
        TLorentzVector t2 = aeneutrino_v.pop_back();

        Deltaphi_ww->Fill((t1+t2).M());
      }

      kwark1.clear(); kwark2.clear();
      jet1.clear(); jet2.clear(); lepton.clear();
      elektorn.clear(); mion.clear(); eneutrino.clear(); mneutrino.clear();
      aelektorn.clear(); amion.clear(); aeneutrino.clear(); amneutrino.clear();
   }

   mjj->Write();
   mll->Write();
   eta_j1->Write();
   eta_j2->Write();
   Deltaeta_j1j2->Write();
   Deltaeta_ll->Write();
   Deltaphi_ww->Write();

   f->Close();
}

void Go()
{
  LHEF a;
  a.Loop();
  new TBrowser;
}
