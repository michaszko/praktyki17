#define LHEF_cxx
#include "LHEF.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void LHEF::Loop(char* output)
{
  // Deklaracja wszystki zmiennych
  TStopwatch * watch = new TStopwatch();

  Double_t deltaphi, deltaeta_hf, deltaeta, deltaphi_hf, etamax, etamin ;
  Double_t M_h, M_vv, M_fix_2, M_fix_4_, p_vv_z1, p_vv_z2, Delta, p_vv_z;

  Int_t licznik1 = 0, licznik2 = 0, counter = 0;

  std::vector<Int_t> jet, lepton;
  std::vector<Int_t> elektorn, mion, eneutrino, mneutrino;
  std::vector<Int_t> aelektorn, amion, aeneutrino, amneutrino;

  std::vector<TLorentzVector> W_pre, W, W_hf;

  TLorentzVector *temp1 =                new TLorentzVector();
  TLorentzVector *temp2 =                new TLorentzVector();
  TLorentzVector *sum =                 new TLorentzVector();
  TLorentzVector *neutrino_momentum =   new TLorentzVector();
  TLorentzVector *final_momentum =      new TLorentzVector();
  TLorentzVector *ll =                  new TLorentzVector();
  TLorentzVector *jj =                  new TLorentzVector();
  TLorentzVector *Higgs =               new TLorentzVector();

  TVector3 * boost_vector = new TVector3();
  TVector3 * ll1 = new TVector3();
  TVector3 * ll2 = new TVector3();
  TVectorD weight(1);

  TH1F * mjj = new TH1F("mjj", "mjj", 50, 0., 1000.);
  TH1F * mll = new TH1F("mll", "mll", 50, 0., 300);

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

  TH1F * nn_p = new TH1F("nn_p", "neutrinos_momentum", 50, 0., 400);
  TH1F * non_nn_p = new TH1F("non_nn_p", "non-neutrinos_momentum", 50, 0., 400);
  TH1F * all_momentum = new TH1F("all_momentum", "all_momentum", 50, 0., 400);

  TH1F * Pdg_id = new TH1F("Pdg_id", "Pdg_id", 50, -25., 25);
  TH1F * delta = new TH1F("delta", "delta", 100, -1e12, 1e12);

  TH1F * M_vv_h = new TH1F("M_vv_h", "M_vv", 50, 0., 600);
  TH1F * p_vv_z_h_c_1 = new TH1F("p_vv_z_c_1", "p_vv_z_calculated_|psi_vv_z|", 100, -2000., 2000);
  TH1F * p_vv_z_h_c_2 = new TH1F("p_vv_z_c_2", "p_vv_z_calculated_|p_vv_z|", 100, -2000., 2000);
  TH1F * p_vv_z_h = new TH1F("p_vv_z", "p_vv_z", 100, -2000., 2000);

  TH1F * M_h_T = new TH1F("M_h_T", "M_h_T", 100, 0., 400);

  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //Pętla po wszystkich przypadkach
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
     //Czyszczenie wszystkich wektorów
     Higgs->SetPxPyPzE(0,0,0,0); ll->SetPxPyPzE(0,0,0,0); jj->SetPxPyPzE(0,0,0,0);
     neutrino_momentum->SetPxPyPzE(0,0,0,0); final_momentum->SetPxPyPzE(0,0,0,0);
     jet.clear(); lepton.clear();
     elektorn.clear(); mion.clear(); eneutrino.clear(); mneutrino.clear();
     aelektorn.clear(); amion.clear(); aeneutrino.clear(); amneutrino.clear();
     W.clear(); W_hf.clear(); W_pre.clear();

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      /************************************************************************/
      //Pętla po wszystkich particlach w evencie
      for (int i = 0; i < Particle_size; i++)
      {
        //Sprawdzenie czy cząstka jest wychodząca
        if (Particle_Status[i] != 1) continue;

        //W jet trzymamy kwarki
        if (abs(Particle_PID[i]) < 9 || Particle_PID[i] == 21 )
          jet.push_back(i);

        //W lepton trzymamy leptony
        else if(  abs(Particle_PID[i]) == 11 ||
                  abs(Particle_PID[i]) == 13 ||
                  abs(Particle_PID[i]) == 15 ||
                  abs(Particle_PID[i]) == 17)
          lepton.push_back(i);

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
          default:
            break;
        }

      }

      /************************************************************************/
      //Sprawdzamy czy jest to nasz proces -
      //czyli na wyjściu jest e- i mu+ ALBO e+ i mu-
      //chyba w końcu działa
      if (aelektorn.size() == elektorn.size()) continue;
      if (jet.size() % 2 != 0) continue;

      /************************************************************************/
      //Rzeczy potrzebne do cięć
      deltaeta = 0; deltaphi = 0;
      for (int i = 0; i < jet.size(); i++)
      {
        temp1->SetPxPyPzE( Particle_Px[  jet[i] ],
                          Particle_Py[  jet[i] ],
                          Particle_Pz[  jet[i] ],
                          Particle_E[   jet[i] ] );

        *jj += *temp1;

        //Liczę deltaeta dla jetów
        if (deltaeta == 0)  deltaeta += Particle_Eta[ jet[i]  ];
        else                deltaeta -= Particle_Eta[ jet[i]  ];
      }

      //if (temp1->M() > 125 || abs(deltaeta) < 6 || jj->M() < 350 ) continue;
      if (abs(deltaeta) < 4 || jj->M() < 600) continue;

      deltaeta = 0; deltaphi = 0;
      for (int i = 0; i < lepton.size(); i++)
      {
        temp1->SetPxPyPzE( Particle_Px[  lepton[i] ],
                          Particle_Py[  lepton[i] ],
                          Particle_Pz[  lepton[i] ],
                          Particle_E[   lepton[i] ] );

        *ll += *temp1;

        //Liczę deltaphi dla leptonów
        if (deltaphi == 0)  deltaphi += Particle_Phi[ lepton[i]  ];
        else                deltaphi -= Particle_Phi[ lepton[i]  ];

        //Liczę deltaeta dla jetów
        if (deltaeta == 0)  deltaeta += Particle_Eta[ lepton[i]  ];
        else                deltaeta -= Particle_Eta[ lepton[i]  ];
      }

      // if (temp1->M() > 125 || abs(deltaeta) < 4 || ll->M() > 80 || jj->M() < 300) continue;
      //Masa poprzeczna Higgsa
      temp1->SetPxPyPzE(-jj->Px(),-jj->Py(), 0, (*jj + *ll).Pt() + ll->Pt());
      if (temp1->M() > 130 ) continue;

      // if (abs(deltaphi) > 1.5 || ll->M() > 90 ) continue;
      if (ll->M() > 90 || abs(deltaphi) > 1.5) continue;

      counter++;

      /************************************************************************/
      //Obliczenie zetowej skladowej pędu neutrin
      M_h = 126.0; M_vv = 30.0;
      M_fix_2 = (M_h * M_h) - (ll->M() * ll->M()) -
                (M_vv * M_vv) + (2 * ll->Px() * (-(*ll + *jj)).Px()) +
                (2 * ll->Py() * (-(*ll + *jj)).Py());
      M_fix_4_ =  (M_fix_2 * M_fix_2) - (4 * ll->E() * ll->E() * M_vv * M_vv) -
                  (4 * ll->E() * ll->E() * (-(*ll + *jj)).Pt() * (-(*ll + *jj)).Pt());
      Delta = (M_fix_2 * M_fix_2 * ll->Pz() * ll->Pz()) -
              (M_fix_4_ * (ll->Pz() * ll->Pz() - ll->E() * ll->E()));

      delta->Fill(Delta);


      if (Delta < 0 )
      {
        M_vv = 0;
        M_fix_2 = (M_h * M_h) - (ll->M2() * ll->M2()) -
                  (M_vv * M_vv) + (2 * ll->Px() * (-(*ll + *jj)).Px()) +
                  (2 * ll->Py() * (-(*ll + *jj)).Py());
        M_fix_4_ =  (M_fix_2 * M_fix_2) - (4 * ll->E() * ll->E() * M_vv * M_vv) -
                    (4 * ll->E() * ll->E() * (-(*ll + *jj)).Pt() * (-(*ll + *jj)).Pt());
        Delta = (M_fix_2 * M_fix_2 * ll->Pz() * ll->Pz()) -
                (M_fix_4_ * (ll->Pz() * ll->Pz() - ll->E() * ll->E()));

        licznik1++;

        if (Delta < 0 )
        {
          licznik2++;
          continue;
        }
      }

      //////////////////////////////////////////////////////////////////////////

      p_vv_z1 = ((-(ll->Pz()) * M_fix_2) + sqrt(Delta)) /
                ((2 * (ll->Pz() * ll->Pz()) - (ll->E() * ll->E())));
      p_vv_z2 = ((-(ll->Pz()) * M_fix_2) - sqrt(Delta)) /
                ((2 * (ll->Pz() * ll->Pz()) - (ll->E() * ll->E())));

      Double_t pll1, pll2, cosll1, cosll2;

      pll1 =  ((-M_fix_2 * p_vv_z1) +
              sqrt((M_fix_2 * M_fix_2 * p_vv_z1 * p_vv_z1) -
              (p_vv_z1 * p_vv_z1) * ((M_fix_2 * M_fix_2) -
              (4 * ll->E() * ll->E() * p_vv_z1 * p_vv_z1)))) /
              (2 * p_vv_z1 * p_vv_z1);

      pll2 =  ((-M_fix_2 * p_vv_z1) -
              sqrt((M_fix_2 * M_fix_2 * p_vv_z1 * p_vv_z1) -
              (p_vv_z1 * p_vv_z1) * ((M_fix_2 * M_fix_2) -
              (4 * ll->E() * ll->E() * p_vv_z1 * p_vv_z1)))) /
              (2 * p_vv_z1 * p_vv_z1);

      //cout << pll1 - ll->Pz() << "\t" << pll2 - ll->Pz() << endl;

      ll1->SetXYZ(ll->Px(), ll->Py(), pll1);
      ll2->SetXYZ(ll->Px(), ll->Py(), pll2);

      cosll1 = abs(TMath::Sin((ll1->Cross(*ll2)).Angle(TVector3(0,0,1))));

      pll1 =  ((-M_fix_2 * p_vv_z2) +
              sqrt((M_fix_2 * M_fix_2 * p_vv_z2 * p_vv_z2) -
              (p_vv_z2 * p_vv_z2) * ((M_fix_2 * M_fix_2) -
              (4 * ll->E() * ll->E() * p_vv_z2 * p_vv_z2)))) /
              (2 * p_vv_z2 * p_vv_z2);
      pll2 =  ((-M_fix_2 * p_vv_z2) -
              sqrt((M_fix_2 * M_fix_2 * p_vv_z2 * p_vv_z2) -
              (p_vv_z2 * p_vv_z2) * ((M_fix_2 * M_fix_2) -
              (4 * ll->E() * ll->E() * p_vv_z2 * p_vv_z2)))) /
              (2 * p_vv_z2 * p_vv_z2);

      ll1->SetXYZ(ll->Px(), ll->Py(), pll1);
      ll2->SetXYZ(ll->Px(), ll->Py(), pll2);

      cosll2 = abs(TMath::Sin((ll1->Cross(*ll2)).Angle(TVector3(0,0,1))));

      if (cosll1 < cosll2)
        p_vv_z_h_c_1->Fill(p_vv_z1);

      else
        p_vv_z_h_c_1->Fill(p_vv_z2);

        ///////////////////////////////////////////////////////////////////////

      if (abs(p_vv_z1) < abs(p_vv_z2))
      {
        p_vv_z_h_c_2->Fill(p_vv_z1);
        p_vv_z = p_vv_z1;
      }

      else
      {
        p_vv_z_h_c_2->Fill(p_vv_z2);
        p_vv_z = p_vv_z1;
      }

      /************************************************************************/
      //Rekonstruujemy Higgsa
      Higgs->SetPxPyPzE(  -jj->Px(),
                          -jj->Py(),
                          ll->Pz() + p_vv_z,
                          ll->E() + sqrt(30 * 30 + (*jj + *ll).Px() * (*jj + *ll).Px() + (*jj + *ll).Py() * (*jj + *ll).Py()));

      mH->Fill(Higgs->M());

      //Zapamiętujemy BoostVector naszego Higgsa, żeby nie musieć go za
      //każdym razem liczyć od nowa
      *boost_vector = -(Higgs->BoostVector());

      /************************************************************************/
      //Pętla po wszystkich jetach
      //std::cout << jet1.size() << " " << jet2.size() << endl;
      sum->SetPxPyPzE(0,0,0,0); jj->SetPxPyPzE(0,0,0,0);
      etamax = -1e10; etamin = 1e10;
      deltaeta = 0; deltaeta_hf = 0; deltaphi = 0; deltaphi_hf = 0;

      for (int i = 0; i < jet.size(); i++)
      {
        temp1->SetPxPyPzE( Particle_Px[  jet[i] ],
                          Particle_Py[  jet[i] ],
                          Particle_Pz[  jet[i] ],
                          Particle_E[   jet[i] ] );

        *sum += *temp1;

        //Zapisuję czteopęd wszystkich jetów
        *jj += *temp1;

        //Liczę deltaeta dla jetów
        if (deltaeta == 0)  deltaeta += Particle_Eta[ jet[i]  ];
        else                deltaeta -= Particle_Eta[ jet[i]  ];

        //Liczę deltaphi dla jetów
        if (deltaphi == 0)  deltaphi += Particle_Phi[ jet[i]  ];
        else                deltaphi -= Particle_Phi[ jet[i]  ];

        //Liczę eta_j_leading i eta_j_notleading
        if (Particle_Eta[ jet[i]  ] > etamax)
          etamax = Particle_Eta[  jet[i]  ];

        else if (Particle_Eta[ jet[i]  ] < etamin)
          etamin = Particle_Eta[  jet[i]  ];

        //Boostowanie do układu spoczynkowego Higgsa
        temp1->Boost(*boost_vector);

        if (deltaeta_hf == 0) deltaeta_hf += temp1->Eta();
        else                  deltaeta_hf -= temp1->Eta();

        //Liczę deltaphi dla jetów
        if (deltaphi_hf == 0) deltaphi_hf += temp1->Phi();
        else                  deltaphi_hf -= temp1->Phi();

        //Wypełnianie histogramów PDG i final_momentum
        *final_momentum += *temp1;
        Pdg_id->Fill( Particle_PID[ jet[i]  ] );
      }

      //Histogram masy inwariantej dwóch jetów
      mjj->Fill(sum->M());

      //Delta eta dla jetów
      Deltaeta_j1j2->Fill( abs(deltaeta) );

      //Delta phi dla jetów; z cieciami do Pi
      if (abs(deltaphi) < TMath::Pi()) Deltaphi_j1j2->Fill( abs(deltaphi) );
      else Deltaphi_j1j2->Fill( 2 * TMath::Pi() - abs(deltaphi) );

      //Eta dla kwarków leading i sub_leading
      eta_j1->Fill(etamax);
      eta_j2->Fill(etamin);

      //Eta i phi w ukladzie Higgsa
      Deltaeta_j1j2_hf->Fill( abs(deltaeta_hf) );
      Deltaphi_j1j2_hf->Fill( abs(deltaphi_hf) );

      /************************************************************************/
      //Pętla po wszystkich leptonach
      sum->SetPxPyPzE(0,0,0,0); ll->SetPxPyPzE(0,0,0,0);
      etamax = -1e10; etamin = 1e10;
      deltaeta = 0; deltaeta_hf = 0; deltaphi = 0; deltaphi_hf = 0;

      for (int i = 0; i < lepton.size(); i++)
      {
        temp1->SetPxPyPzE( Particle_Px[  lepton[i] ],
                          Particle_Py[  lepton[i] ],
                          Particle_Pz[  lepton[i] ],
                          Particle_E[   lepton[i] ] );

        *sum += *temp1;

        //Zapisuję czteopęd wszystkich jetów
        *ll += *temp1;

        //Liczę deltaeta dla jetów
        if (deltaeta == 0)  deltaeta += Particle_Eta[ lepton[i]  ];
        else                deltaeta -= Particle_Eta[ lepton[i]  ];

        //Liczę deltaphi dla jetów
        if (deltaphi == 0)  deltaphi += Particle_Phi[ lepton[i]  ];
        else                deltaphi -= Particle_Phi[ lepton[i]  ];

        //Liczę eta_j_leading i eta_j_notleading
        if (Particle_Eta[ lepton[i]  ] > etamax)
          etamax = Particle_Eta[  lepton[i]  ];

        else if (Particle_Eta[ lepton[i]  ] < etamin)
          etamin = Particle_Eta[  lepton[i]  ];

        //Boostowanie do układu spoczynkowego Higgsa
        temp1->Boost(*boost_vector);

        if (deltaeta_hf == 0)  deltaeta_hf += temp1->Eta();
        else                   deltaeta_hf -= temp1->Eta();

        //Liczę deltaphi dla jetów
        if (deltaphi_hf == 0)  deltaphi_hf += temp1->Phi();
        else                   deltaphi_hf -= temp1->Phi();

        //Wypełnianie histogramów PDG i final_momentum
        *final_momentum += *temp1;
        Pdg_id->Fill( Particle_PID[ lepton[i]  ] );
      }

      //Histogram masy inwariantej dwóch jetów
      mll->Fill(sum->M());

      //Delta eta dla jetów
      Deltaeta_ll->Fill( abs(deltaeta) );

      //Delta phi dla jetów; z cieciami do Pi
      if (abs(deltaphi) < TMath::Pi()) Deltaphi_ll->Fill( abs(deltaphi) );
      else Deltaphi_ll->Fill( 2 * TMath::Pi() - abs(deltaphi) );

      //Eta dla kwarków leading i sub_leading
      eta_l1->Fill(etamax);
      eta_l2->Fill(etamin);

      //Eta i phi w ukladzie Higgsa
      Deltaeta_ll_hf->Fill( abs(deltaeta_hf) );
      Deltaphi_ll_hf->Fill( abs(deltaphi_hf) );

      /************************************************************************/

      /************************************************************************/
      //Pętla po wszystkich lepotnach w poszukiwaniu W
      //będziemy z tego odtwarzać W+ lub W
      if (!aelektorn.empty() && !eneutrino.empty())
      {
        temp1->SetPxPyPzE( Particle_Px[ aelektorn[0]  ],
                          Particle_Py[  aelektorn[0]  ],
                          Particle_Pz[  aelektorn[0]  ],
                          Particle_E[   aelektorn[0]  ] );
        temp2->SetPxPyPzE( Particle_Px[ eneutrino[0]  ],
                          Particle_Py[  eneutrino[0]  ],
                          Particle_Pz[  eneutrino[0]  ],
                          Particle_E[   eneutrino[0]  ] );

        W_pre.push_back(*temp1);  W_pre.push_back(*temp2);

        //Do pędu porzecznego neutrin i PDGID
        *neutrino_momentum += *temp2;
        Pdg_id->Fill( Particle_PID[eneutrino[0]] );
      }

      if (!amion.empty() && !mneutrino.empty())
      {
        temp1->SetPxPyPzE( Particle_Px[ amion[0] ],
                          Particle_Py[  amion[0]  ],
                          Particle_Pz[  amion[0]  ],
                          Particle_E[   amion[0]  ] );
        temp2->SetPxPyPzE( Particle_Px[ mneutrino[0]  ],
                          Particle_Py[  mneutrino[0]  ],
                          Particle_Pz[  mneutrino[0]  ],
                          Particle_E[   mneutrino[0]  ] );

        W_pre.push_back(*temp1);  W_pre.push_back(*temp2);

        *neutrino_momentum += *temp2;
        Pdg_id->Fill( Particle_PID[mneutrino[0]] );
      }

      if (!mion.empty() && !amneutrino.empty())
      {
        temp1->SetPxPyPzE( Particle_Px[ mion[0]       ],
                          Particle_Py[  mion[0]       ],
                          Particle_Pz[  mion[0]       ],
                          Particle_E[   mion[0]       ] );
        temp2->SetPxPyPzE( Particle_Px[ amneutrino[0] ],
                          Particle_Py[  amneutrino[0] ],
                          Particle_Pz[  amneutrino[0] ],
                          Particle_E[   amneutrino[0] ] );

        W_pre.push_back(*temp1);  W_pre.push_back(*temp2);

        *neutrino_momentum += *temp2;
        Pdg_id->Fill( Particle_PID[amneutrino[0]] );
      }

      if (!elektorn.empty() && !aeneutrino.empty())
      {
        temp1->SetPxPyPzE( Particle_Px[ elektorn[0]   ],
                          Particle_Py[  elektorn[0]   ],
                          Particle_Pz[  elektorn[0]   ],
                          Particle_E[   elektorn[0]   ] );
        temp2->SetPxPyPzE( Particle_Px[ aeneutrino[0] ],
                          Particle_Py[  aeneutrino[0] ],
                          Particle_Pz[  aeneutrino[0] ],
                          Particle_E[   aeneutrino[0] ] );

        W_pre.push_back(*temp1);  W_pre.push_back(*temp2);

        *neutrino_momentum += *temp2;
        Pdg_id->Fill( Particle_PID[aeneutrino[0]] );
      }

      for (int i = 0; i < W_pre.size(); i += 2)
      {
        *sum = W_pre[i] + W_pre[i+1];
        //W wektorze zapamiętujemy sumę czterowektorów mu- i vm~ czyli
        //tak naprawdę czteropęd W-
        W.push_back(*sum);

        //Boostujemy do ukladu spoczynkowego Higgsa
        temp1->Boost(*boost_vector);
        temp2->Boost(*boost_vector);

        *sum = W_pre[i] + W_pre[i+1];

        //Zapamiętujemy boostowane czterowektory w wektorze
        W_hf.push_back(*sum);
      }


      /************************************************************************/
      //Pętla po bozonach W+ i W-
      for (int i = 0; i < W.size(); i++)
      {
        //Różnica phi
        if (deltaphi == 0)  deltaphi += W[i].Phi();
        else                deltaphi -= W[i].Phi();

        //Różnica phi po boostcie
        if (deltaphi_hf == 0) deltaphi_hf += W[i].Phi();
        else                  deltaphi_hf -= W[i].Phi();
      }

      //Cięcie (normowanie) do Pi
      if (abs(deltaphi) < TMath::Pi()) Deltaphi_ww->Fill( abs(deltaphi) );
      else Deltaphi_ww->Fill( 2 * TMath::Pi() - abs(deltaphi) );

      //Cięcie (normowanie) do Pi
      if (abs(deltaphi_hf) < TMath::Pi()) Deltaphi_ww_hf->Fill( abs(deltaphi_hf) );
      else Deltaphi_ww_hf->Fill( 2 * TMath::Pi() - abs(deltaphi_hf) );

      /***********************************************************************/
      //Pęd poprzeczny neutrin
      nn_p->Fill(neutrino_momentum->Pt());
      p_vv_z_h->Fill(neutrino_momentum->Pz());
      M_vv_h->Fill(neutrino_momentum->M());

      /************************************************************************/
      //Pęd poprzeczny wszystkich mozliwych do detekcji czastek
      non_nn_p->Fill(final_momentum->Pt());

      /************************************************************************/
      //Pęd wszystkich cząstek (neutrina + rest)
      all_momentum->Fill((*final_momentum + *neutrino_momentum).Pt());

      /************************************************************************/


      /************************************************************************/
      //Masa poprzeczna Higgsa
      temp1->SetPxPyPzE(-jj->Px(),-jj->Py(), 0, (*jj + *ll).Pt() + ll->Pt());
      M_h_T->Fill(temp1->M());

      //std::cout << counter/nentries << endl;
      //Zapamiętanie przekroju czynnego w wektorze
      //if (jentry == nentries)
      {
        Double_t lol = (Double_t)counter/(Double_t)nentries;
        //cout << *Event_Weight<< " " << ((*Event_Weight) * lol) << endl;
        weight[0] = (*Event_Weight) * lol;
        //cout << weight[0] << endl;
      }

   }

   cout << "Delta < 0 - 1 raz  " << licznik1 << endl << "Delta < 0 - 2 razy  " << licznik2 << endl;

   //Zapisywanie do pliki "output"
   TFile *f = new TFile(output, "RECREATE");

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

   nn_p->Write();
   //non_nn_p->Write();
   //all_momentum->Write();

   //Pdg_id->Write();
   //delta->Write();

   M_vv_h->Write();
   p_vv_z_h_c_1->Write();
   p_vv_z_h_c_2->Write();
   p_vv_z_h->Write();

   M_h_T->Write();

   weight.Write("weight");

   f->Close();
}

void ReadLine(char *namefile, int setlog = 0)
{
	ifstream input;
	input.open(namefile);
	string line, name;
	double norm;
	vector<string> alllines;
	vector<string> legendtitles;
	vector<TFile *> files;
	TObject *obj;
	TH1F *h;
	TKey *key;
	vector<TCanvas *> cvec;
	vector<TIter> nexts;
	TVectorT<double> weight;
	vector<double> weights;
	weight.ResizeTo(1);
	int num1 = 0;
	double scale = 1, max = 0;

	while (getline(input, line))
	{
		if (line.front() == '#')
		{
			cout << line << endl;
			continue;
		}

		alllines.push_back(line);
		// line.resize(line.size() - 5);
		cout << line << endl;
		legendtitles.push_back(line);
		num1++;
	}

	cout << alllines.size() << endl;

	if (num1 % 3 != 0)
	{
		cout << "Fatal error - invalid parameter file!" << endl;
		cout << "Please read first lines for proper form" << endl;
	}

	cout << num1 << endl;

	vector<TH1F *> objs[200];

	for (int i = 0; i < alllines.size(); i++)
	{
		// cout << alllines.at(i) << endl;
		name = alllines.at(i);
		files.push_back(new TFile(name.c_str(), "READ"));
		nexts.push_back(TIter(files.at(i)->GetListOfKeys()));
		while ((key = (TKey *)nexts.at(i)()))
		{
			obj = files.at(i)->Get(key->GetName());
			if (obj->InheritsFrom("TH1"))
				objs[i].push_back((TH1F *)obj);
			if (obj->InheritsFrom("TVectorD"))
			{
				cout << "Doby found a weight!" << endl;
				cout << obj << endl;
				obj->Print();
				weight = *(TVectorD *)obj;
				//cout << weight[0] << endl;
				weights.push_back(weight[0]);
			}
			// cout << " found object: " << obj->GetName() << endl;
		}
	}
	for (int i = 0; i < objs[0].size(); i++)
	{
		max = 0;
		if (i % 4 == 0)
		{
			cvec.push_back(new TCanvas(Form("c%d", i + 2), "just a canvas", 1000, 750));
			cvec.back()->Divide(2, 2);
			gStyle->SetOptStat(0);
		}
		cvec.back()->cd(i % 4 + 1);
		for (int j = 0; j < num1; j++)
		{
			// cout << objs[j].at(i)->Integral() << endl;
			objs[j].at(i)->ComputeIntegral();
			scale = weights.at(j) / (objs[j].at(i)->Integral());
			objs[j].at(i)->Scale((Double_t)scale, "width");
			// cout << scale << endl;
			// if (objs[j].at(i)->Integral() != 0)
			// objs[j].at(i)->Scale(1. / objs[j].at(i)->Integral());
			if (objs[j].at(i)->GetMaximum() > max)
				max = objs[j].at(i)->GetMaximum();
			// if (objs[j].at(i)->)
		}

		TLegend *legend = new TLegend(0.9, 0.7, 0.52, 0.9);
		legend->SetHeader("Procesy");

    THStack *hs = new THStack("hs","");

    for (int j = 0; j < num1; j++)
		{
			if (setlog == 1)
				gPad->SetLogy();
			h = objs[j].at(i);
			legend->AddEntry(h, legendtitles.at(j).c_str(), "l");
			// cout << objs[j].at(i) << endl;
			objs[j].at(i)->SetMaximum(1.05 * max);
			// cout << max << endl;
			objs[j].at(i)->SetLineColor(j + 2);
			objs[j].at(i)->SetFillColor(j + 2);
			// objs[j].at(i)->SetLineStyle(j + 1);
			hs->Add(objs[j].at(i));
      hs->SetTitle(objs[j].at(i)->GetTitle());
			//legend->Draw();
			// if (i % 4 == 0 && j == num1 - 1)
			// cvec.back()->SaveAs(Form("picture%d.jpg", i / 4));
		}

    hs->Draw("hist nostackb");
    legend->Draw();

    // TLegend *legend = new TLegend(0.9, 0.7, 0.52, 0.9);
    // legend->SetHeader("Procesy");
    //
    // for (int j = 0; j < num1; j++)
    // {
    //   if (setlog == 1)
    //     gPad->SetLogy();
    //   h = objs[j].at(i);
    //   legend->AddEntry(h, legendtitles.at(j).c_str(), "l");
    //   // cout << objs[j].at(i) << endl;
    //   objs[j].at(i)->SetMaximum(1.05 * max);
    //   // cout << max << endl;
    //   objs[j].at(i)->SetLineColor(j + 2);
    //   // objs[j].at(i)->SetLineStyle(j + 1);
    //   if (j == 0)
    //     objs[j].at(i)->Draw("hist");
    //   else
    //     objs[j].at(i)->Draw("same hist");
    //   legend->Draw();
    //   // if (i % 4 == 0 && j == num1 - 1)
    //   // cvec.back()->SaveAs(Form("picture%d.jpg", i / 4));
    // }
	}
}


void LHEF(char* input, char* output)
{
  //Dla ROOT 6.x
  class LHEF a(input);

  //Dla ROOT 5.x
  //LHEF a;

  a.Loop(output);
  //new TBrowser;
}

void Go()
{
  LHEF("unweighted_events.root", "Histogramy_1.root");
  LHEF("VBF_AL1.0AT1.0.root", "Histogramy_1010.root");
  LHEF("QCD.root", "Histogramy_QCD.root");
  LHEF("EW.root", "Histogramy_ew.root");

  ReadLine("draw.config");
}
