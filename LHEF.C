#define LHEF_cxx
#include "LHEF.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void LHEF::Loop()
{

  // Deklaracja wszystki zmiennych
  TStopwatch * watch = new TStopwatch();

  Double_t phi;

  std::vector<Int_t> elektorn, mion, eneutrino, mneutrino;
  std::vector<Int_t> aelektorn, amion, aeneutrino, amneutrino;
  std::vector<Int_t> jet1, jet2;
  std::vector<Int_t> lepton1, lepton2;
  std::vector<TLorentzVector> W_p, W_n, W_p_hf, W_n_hf;
  TLorentzVector temp1, temp2, sum;

  TLorentzVector *Higgs = new TLorentzVector();
  TVector3 * boost_vector = new TVector3();

  TH1F * mjj = new TH1F("mjj", "mjj", 50, 0., 1000.);
  TH1F * mll = new TH1F("mll", "mll", 50, 0., 1000);

  TH1F * eta_j1  = new TH1F("eta_j1", "eta_j1", 50, -5., 5.);
  TH1F * eta_j2 = new TH1F("eta_j2", "eta_j2", 50, -5., 5.);

  TH1F * eta_l1  = new TH1F("eta_l1", "eta_l1", 50, -5., 5.);
  TH1F * eta_l2 = new TH1F("eta_l2", "eta_l2", 50, -5., 5.);

  TH1F * Deltaeta_j1j2 = new TH1F("Deltaeta_j1j2", "Deltaeta_j1j2", 50, 0., 10.);
  TH1F * Deltaeta_j1j2_hf = new TH1F("Deltaeta_j1j2_hf", "Deltaeta_j1j2_hf", 50, 0., 10.);

  TH1F * Deltaeta_ll = new TH1F("Deltaeta_ll", "Deltaeta_ll", 50, 0., 10.);
  TH1F * Deltaeta_ll_hf = new TH1F("Deltaeta_ll_hf", "Deltaeta_ll_hf", 50, 0., 10.);

  TH1F * Deltaphi_ww = new TH1F("Deltaphi_ww", "Deltaphi_ww", 50, 0., TMath::Pi());
  TH1F * Deltaphi_ww_hf = new TH1F("Deltaphi_ww_hf", "Deltaphi_ww_hf", 50, 0., TMath::Pi());

  TH1F * Deltaphi_j1j2 = new TH1F("Deltaphi_j1j2", "Deltaphi_j1j2", 50, 0., TMath::Pi());
  TH1F * Deltaphi_j1j2_hf = new TH1F("Deltaphi_j1j2_hf", "Deltaphi_j1j2_hf", 50, 0., TMath::Pi());

  TH1F * Deltaphi_ll = new TH1F("Deltaphi_ll", "Deltaphi_ll", 50, 0., TMath::Pi());
  TH1F * Deltaphi_ll_hf = new TH1F("Deltaphi_ll_hf", "Deltaphi_ll_hf", 50, 0., TMath::Pi());

  TH1F * mH = new TH1F("mH", "mH", 50, 0., 1000);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //Pętla po wszystkich przypadkach
   for (Long64_t jentry=0; jentry<nentries*1;jentry++)
   {
     //Czyszczenie wszystkich wektorów
     Higgs->Clear();
     jet1.clear(); jet2.clear(); lepton1.clear(); lepton2.clear();
     elektorn.clear(); mion.clear(); eneutrino.clear(); mneutrino.clear();
     aelektorn.clear(); amion.clear(); aeneutrino.clear(); amneutrino.clear();
     W_p.clear(); W_n.clear(); W_p_hf.clear(); W_n_hf.clear();

      //Zegarek (jeśli działanie programu się przedłuża)
      // if ( jentry%1000 == 0 )
      // {
      //   watch->Stop();
      //   watch->Print();
      //   std::cout << jentry << endl;
      //   watch->Continue();
      // }

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      /***************************************/
      //Pętla po wszystkich particlach w evencie
      for (int i = 0; i < Particle_size; i++)
      {
        //Sprawdzenie czy cząstka jest wychodząca
        if (Particle_Status[i] != 1) continue;

        //W jet1 trzymamy kwarki...
        if (Particle_PID[i] >= 1 && Particle_PID[i] <= 8)
          jet1.push_back(i);

        //... a w jet2 antykwarki
        else if (Particle_PID[i] >= -8 && Particle_PID[i] <= -1)
          jet2.push_back(i);

        //W lepton1 trzymamy leptony...
        else if(  Particle_PID[i] == 11 ||
                  Particle_PID[i] == 13 ||
                  Particle_PID[i] == 15 ||
                  Particle_PID[i] == 17)
          lepton1.push_back(i);

        //... a w lepton2 antyleptony
        else if(  Particle_PID[i] == -11 ||
                  Particle_PID[i] == -13 ||
                  Particle_PID[i] == -15 ||
                  Particle_PID[i] == -17)
          lepton2.push_back(i);

        //Dodatkowo jeszcze segregujemy leptony
        switch (Particle_PID[i]) {
          case 11:
            elektorn.push_back(i);
            break;
          case 12:
            eneutrino.push_back(i);
            break;
          case 13:
            mion.push_back(i);
            break;
          case 14:
            mneutrino.push_back(i);
            break;
          case -11:
            aelektorn.push_back(i);
          break;
          case -12:
            aeneutrino.push_back(i);
          break;
          case -13:
            amion.push_back(i);
          break;
          case -14:
            amneutrino.push_back(i);
          break;
        }

      }

      /**************************************/
      //Sprawdzamy czy jest to nasz proces -
      //czyli na wyjściu jest e- i mu+ ALBO e+ i mu-
      if (elektorn.size() != amion.size() ) continue;
      else if (aelektorn.size() != mion.size() ) continue;

      //std::cout << Particle_PID[jet1[0]] << "\t" << Particle_PID[jet2[0]] << endl;
      //std::cout << "\t" << lepton1.size() << "\t" << lepton2.size() << endl;

      /**************************************/
      //Rekonstruujemy Higgsa w zależności od procesu
      if ( aelektorn.size() == 1 && eneutrino.size() == 1 && mion.size() == 1 && amneutrino.size() == 1 )
        Higgs->SetPxPyPzE(  Particle_Px[aelektorn[0]] +  Particle_Px[eneutrino[0]] + Particle_Px[mion[0]] + Particle_Px[amneutrino[0]],
                            Particle_Py[aelektorn[0]] +  Particle_Py[eneutrino[0]] + Particle_Py[mion[0]] + Particle_Py[amneutrino[0]],
                            Particle_Pz[aelektorn[0]] +  Particle_Pz[eneutrino[0]] + Particle_Pz[mion[0]] + Particle_Pz[amneutrino[0]],
                            Particle_E[aelektorn[0]]  +  Particle_E[eneutrino[0]]  + Particle_E[mion[0]]  + Particle_E[amneutrino[0]]);

      else if ( elektorn.size() == 1 && aeneutrino.size() == 1 && amion.size() == 1 && mneutrino.size() == 1 )
        Higgs->SetPxPyPzE(  Particle_Px[elektorn[0]] +  Particle_Px[aeneutrino[0]] + Particle_Px[amion[0]] + Particle_Px[mneutrino[0]],
                            Particle_Py[elektorn[0]] +  Particle_Py[aeneutrino[0]] + Particle_Py[amion[0]] + Particle_Py[mneutrino[0]],
                            Particle_Pz[elektorn[0]] +  Particle_Pz[aeneutrino[0]] + Particle_Pz[amion[0]] + Particle_Pz[mneutrino[0]],
                            Particle_E[elektorn[0]]  +  Particle_E[aeneutrino[0]]  + Particle_E[amion[0]]  + Particle_E[mneutrino[0]]);

      mH->Fill(Higgs->M());

      //Zapamiętujemy BoostVector naszego Higgsa, żeby nie musieć go za
      //każdym razem liczyć od nowa
      *boost_vector = -(Higgs->BoostVector());

      /***************************************/

      //Pętla po wszystkich jetach i antyjetach
      for (int i = 0; i < jet1.size(); i++)
      {
          for (int j = 0; j < jet2.size(); j++)
          {
            temp1.SetPxPyPzE( Particle_Px[  jet1[i] ],
                              Particle_Py[  jet1[i] ],
                              Particle_Pz[  jet1[i] ],
                              Particle_E[   jet1[i] ] );

            temp2.SetPxPyPzE( Particle_Px[  jet2[j] ],
                              Particle_Py[  jet2[j] ],
                              Particle_Pz[  jet2[j] ],
                              Particle_E[   jet2[j] ] );

            sum = temp1 + temp2;

            //Histogram masy inwariantej dwóch jetów
            mjj->Fill(sum.M());

            ///////////////////////////////
            //Delta eta dla jetów
            Deltaeta_j1j2->Fill( abs(Particle_Eta[jet1[i]] - Particle_Eta[jet2[j]]) );

            //////////////////////////////
            //Delta phi dla jetów
            //z cieciami do Pi
            phi = abs(Particle_Phi[jet1[i]] - Particle_Phi[jet2[j]]);

            if (phi < TMath::Pi()) Deltaphi_j1j2->Fill( phi );
            else Deltaphi_j1j2->Fill( 2 * TMath::Pi() - phi );

            //////////////////////////////
            //Eta dla antykwarkow
            if (i == 0) eta_j2->Fill(Particle_Eta[jet2[j]]);

            /////////////////////////////
            //Boostowanie do układu spoczynkowego Higgsa
            if (Higgs == 0) continue;

            temp1.Boost(*boost_vector);
            temp2.Boost(*boost_vector);

            Deltaeta_j1j2_hf->Fill( abs(temp1.Eta() - temp2.Eta()) );
            Deltaphi_j1j2_hf->Fill( abs(temp1.Phi() - temp2.Phi()) );
          }
          //Eta dla kwarków
          eta_j1->Fill(Particle_Eta[jet1[i]]);
      }

      /***************************************/
      //Pętla po wszystkich leptonach i antyleptonach
      for (int i = 0; i < lepton1.size(); i++)
      {
          for (int j = 0; j < lepton2.size(); j++)
          {
              //if( Particle_PID[lepton1[i]] == -Particle_PID[lepton2[j]] ) continue;

              temp1.SetPxPyPzE( Particle_Px[  lepton1[i]  ],
                                Particle_Py[  lepton1[i]  ],
                                Particle_Pz[  lepton1[i]  ],
                                Particle_E[   lepton1[i] ] );

              temp2.SetPxPyPzE( Particle_Px[  lepton2[j]  ],
                                Particle_Py[  lepton2[j]  ],
                                Particle_Pz[  lepton2[j]  ],
                                Particle_E[   lepton2[j]  ] );

              sum = temp1 + temp2;

              //Histogram masy inwariantnej dwóch leptonów -
              // e+ i mu- ALBO e- i mu+ - bo tylko te wybraliśmy
              mll->Fill(sum.M());

              ////////////////////////////////
              //Delta eta dla leptonów
              Deltaeta_ll->Fill( abs(Particle_Eta[lepton1[i]] - Particle_Eta[lepton2[j]]) );

              ///////////////////////////////
              //Delta phi dla leptonów
              //z cieciami do Pi
              phi = abs(Particle_Phi[lepton1[i]] - Particle_Phi[lepton2[j]]);

              if (phi < TMath::Pi()) Deltaphi_ll->Fill( phi );
              else Deltaphi_ll->Fill( (2 * TMath::Pi() - phi) );

              ///////////////////////////////
              //Eta dla antyleptonów
              if (i == 0) eta_l2->Fill(Particle_Eta[lepton2[j]]);

              ///////////////////////////////

              if (Higgs == 0) continue;
              //Boostowanie do układu spoczynkowego Higgsa
              temp1.Boost(*boost_vector);
              temp2.Boost(*boost_vector);

              Deltaeta_ll_hf->Fill( abs(temp1.PseudoRapidity() - temp2.PseudoRapidity()) );
              Deltaphi_ll_hf->Fill( abs(temp1.Phi() - temp2.Phi()) );
          }
          //Eta dla leptonów
          eta_l1->Fill(Particle_Eta[lepton1[i]]);
      }

      /***************************************/
      //Pętla po wszystkich e+ i ve
      //będziemy z tego odtwarzać W+
      for (int i = 0; i < aelektorn.size(); i++)
      {
          for (int j = 0; j < eneutrino.size(); j++)
          {
              temp1.SetPxPyPzE( Particle_Px[  aelektorn[i]  ],
                                Particle_Py[  aelektorn[i]  ],
                                Particle_Pz[  aelektorn[i]  ],
                                Particle_E[   aelektorn[i]  ] );
              temp2.SetPxPyPzE( Particle_Px[  eneutrino[j]  ],
                                Particle_Py[  eneutrino[j]  ],
                                Particle_Pz[  eneutrino[j]  ],
                                Particle_E[   eneutrino[j]  ] );

              sum = temp1 + temp2;

              //W wektorze zapamiętujemy sumę czterowektorów e+ i ve czyli
              //tak naprawdę czteropęd W+
              W_p.push_back(sum);

              ///////////////////////////////////
              //Boostujemy do ukladu spoczynkowego Higgsa
              if (Higgs == 0) continue;

              temp1.Boost(*boost_vector);
              temp2.Boost(*boost_vector);

              sum = temp1 + temp2;

              //Zapamiętujemy boostowane czterowektory w wektorze
              W_p_hf.push_back( sum );
          }
      }

      /***************************************/
      //Pętla po wszystkich mu- i vm~
      //będziemy z tego odtwarzać W+
      for (int i = 0; i < mion.size(); i++)
      {
          for (int j = 0; j < amneutrino.size(); j++)
          {
              temp1.SetPxPyPzE( Particle_Px[  mion[i]       ],
                                Particle_Py[  mion[i]       ],
                                Particle_Pz[  mion[i]       ],
                                Particle_E[   mion[i]       ] );
              temp2.SetPxPyPzE( Particle_Px[  amneutrino[j] ],
                                Particle_Py[  amneutrino[j] ],
                                Particle_Pz[  amneutrino[j] ],
                                Particle_E[   amneutrino[j] ] );

              sum = temp1 + temp2;
              //W wektorze zapamiętujemy sumę czterowektorów mu- i vm~ czyli
              //tak naprawdę czteropęd W-
              W_n.push_back(sum);

              ///////////////////////////////////
              //Boostujemy do ukladu spoczynkowego Higgsa
              if (Higgs == 0) continue;

              temp1.Boost(*boost_vector);
              temp2.Boost(*boost_vector);

              sum = temp1 + temp2;

              //Zapamiętujemy boostowane czterowektory w wektorze
              W_n_hf.push_back( sum );

          }
      }

      /***************************************/
      //Pętla po bozonach W+ i W-
      for (int i = 0; i < W_n.size(); i++)
      {
          for (int j = 0; j < W_p.size(); j++)
          {
              //Różnica phi między nimi
              phi = abs(W_n[i].Phi() - W_p[j].Phi());

              //Cięcie (normowanie) do Pi
              if (phi < TMath::Pi()) Deltaphi_ww->Fill( phi );
              else Deltaphi_ww->Fill( 2 * TMath::Pi() - phi );

              ////////////////////////////////////////////////
              //Różnica phi między nimi (po booscie)
              phi = abs(W_n_hf[i].Phi() - W_p_hf[j].Phi());

              //Cięcie (normowanie) do Pi
              if (phi < TMath::Pi()) Deltaphi_ww_hf->Fill( phi );
              else Deltaphi_ww_hf->Fill( 2 * TMath::Pi() - phi );
          }
      }
   }

   //Zapisywanie do pliki "Histogramy.root"
   TFile *f = new TFile("Histogramy.root", "RECREATE");

   mjj->Write();
   mll->Write();

   eta_j1->Write();
   eta_j2->Write();
   eta_l1->Write();
   eta_l2->Write();

   Deltaeta_j1j2->Write();
   Deltaeta_j1j2_hf->Write();
   Deltaeta_ll->Write();
   Deltaeta_ll_hf->Write();

   Deltaphi_ll->Write();
   Deltaphi_ll_hf->Write();
   Deltaphi_j1j2->Write();
   Deltaphi_j1j2_hf->Write();
   Deltaphi_ww->Write();
   Deltaphi_ww_hf->Write();

   mH->Write();

   f->Close();
}

void LHEF()
{
  //Dla ROOT 6.x
  class LHEF a;

  //Dla ROOT 5.x
  //LHEF a;

  a.Loop();
  new TBrowser;
}
