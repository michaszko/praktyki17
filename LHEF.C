#define LHEF_cxx
#include "LHEF.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void LHEF::Loop(char* output)
{
  // Deklaracja wszystki zmiennych
  TStopwatch * watch = new TStopwatch();

  Double_t phi;
  Double_t M_h, M_vv, M_fix_2, M_fix_4_, p_vv_z;


  std::vector<Int_t> elektorn, mion, eneutrino, mneutrino;
  std::vector<Int_t> aelektorn, amion, aeneutrino, amneutrino;
  std::vector<Int_t> jet1, jet2;
  std::vector<Int_t> lepton1, lepton2;

  std::vector<TLorentzVector> W_p, W_n, W_p_hf, W_n_hf;
  TLorentzVector temp1, temp2, sum, final_momentum, neutrino_momentum, end_momentum;

  TLorentzVector *ll = new TLorentzVector();
  TLorentzVector *jj = new TLorentzVector();
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

  TH1F * nn_p = new TH1F("nn_p", "neutrinos_momentum", 50, 0., 400);
  TH1F * non_nn_p = new TH1F("non_nn_p", "non-neutrinos_momentum", 50, 0., 400);
  TH1F * all_momentum = new TH1F("all_momentum", "all_momentum", 50, 0., 400);

  TH1F * Pdg_id = new TH1F("Pdg_id", "Pdg_id", 30, -15., 15);

  TH1F * p_vv_z_h = new TH1F("p_vv_z", "p_vv_z", 50, -1000., 1000);

  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //Pętla po wszystkich przypadkach
   for (Long64_t jentry=0; jentry<nentries*1;jentry++)
   {
     //Czyszczenie wszystkich wektorów
     Higgs->SetPxPyPzE(0,0,0,0); ll->SetPxPyPzE(0,0,0,0); jj->SetPxPyPzE(0,0,0,0);
     jet1.clear(); jet2.clear(); lepton1.clear(); lepton2.clear();
     elektorn.clear(); mion.clear(); eneutrino.clear(); mneutrino.clear();
     aelektorn.clear(); amion.clear(); aeneutrino.clear(); amneutrino.clear();
     W_p.clear(); W_n.clear(); W_p_hf.clear(); W_n_hf.clear();

      //Zegarek (jeśli działanie programu się przedłuża)
      if ( jentry%1000 == 0 )
      {
        watch->Stop();
        watch->Print();
        std::cout << jentry << endl;
        watch->Continue();
      }

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

      /************************************************************************/
      //Sprawdzamy czy jest to nasz proces -
      //czyli na wyjściu jest e- i mu+ ALBO e+ i mu-
      //chyba w końcu działa
      if (aelektorn.size() == elektorn.size()) continue;
      if (jet1.empty() || jet2.empty()) continue;

      /************************************************************************/
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

      /************************************************************************/
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

            //Zapisuję czteopęd wszystkich jetów
            *jj += sum;

            //Histogram masy inwariantej dwóch jetów
            mjj->Fill(sum.M());

            ////////////////////////////////////////////////////////////////////
            //Delta eta dla jetów
            Deltaeta_j1j2->Fill( abs(Particle_Eta[jet1[i]] - Particle_Eta[jet2[j]]) );

            ////////////////////////////////////////////////////////////////////
            //Delta phi dla jetów
            //z cieciami do Pi
            phi = abs(Particle_Phi[jet1[i]] - Particle_Phi[jet2[j]]);

            if (phi < TMath::Pi()) Deltaphi_j1j2->Fill( phi );
            else Deltaphi_j1j2->Fill( 2 * TMath::Pi() - phi );

            ////////////////////////////////////////////////////////////////////
            //Eta dla antykwarkow
            if (i == 0) eta_j2->Fill(Particle_Eta[jet2[j]]);

            ///////////////////////////////////////////////////////////////////
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

      /************************************************************************/
      //Pętla po wszystkich leptonach i antyleptonach
      for (int i = 0; i < lepton1.size(); i++)
      {
          for (int j = 0; j < lepton2.size(); j++)
          {
              temp1.SetPxPyPzE( Particle_Px[  lepton1[i]  ],
                                Particle_Py[  lepton1[i]  ],
                                Particle_Pz[  lepton1[i]  ],
                                Particle_E[   lepton1[i] ] );

              temp2.SetPxPyPzE( Particle_Px[  lepton2[j]  ],
                                Particle_Py[  lepton2[j]  ],
                                Particle_Pz[  lepton2[j]  ],
                                Particle_E[   lepton2[j]  ] );

              sum = temp1 + temp2;

              //Zapamiętujemy czteropęd leptonów
              *ll += sum;

              //Histogram masy inwariantnej dwóch leptonów -
              // e+ i mu- ALBO e- i mu+ - bo tylko te wybraliśmy
              mll->Fill(sum.M());

              //////////////////////////////////////////////////////////////////
              //Delta eta dla leptonów
              Deltaeta_ll->Fill( abs(Particle_Eta[lepton1[i]] - Particle_Eta[lepton2[j]]) );

              //////////////////////////////////////////////////////////////////
              //Delta phi dla leptonów
              //z cieciami do Pi
              phi = abs(Particle_Phi[lepton1[i]] - Particle_Phi[lepton2[j]]);

              if (phi < TMath::Pi()) Deltaphi_ll->Fill( phi );
              else Deltaphi_ll->Fill( (2 * TMath::Pi() - phi) );

              //////////////////////////////////////////////////////////////////
              //Eta dla antyleptonów
              if (i == 0) eta_l2->Fill(Particle_Eta[lepton2[j]]);

              //////////////////////////////////////////////////////////////////

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

      /************************************************************************/
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

              /////////////////////////////////////////////////////////////////
              //Boostujemy do ukladu spoczynkowego Higgsa
              if (Higgs == 0) continue;

              temp1.Boost(*boost_vector);
              temp2.Boost(*boost_vector);

              sum = temp1 + temp2;

              //Zapamiętujemy boostowane czterowektory w wektorze
              W_p_hf.push_back( sum );
          }
      }

      /************************************************************************/
      //Pętla po wszystkich mu+ i vm
      //będziemy z tego odtwarzać W+
      for (int i = 0; i < amion.size(); i++)
      {
          for (int j = 0; j < mneutrino.size(); j++)
          {
              temp1.SetPxPyPzE( Particle_Px[  amion[i]  ],
                                Particle_Py[  amion[i]  ],
                                Particle_Pz[  amion[i]  ],
                                Particle_E[   amion[i]  ] );
              temp2.SetPxPyPzE( Particle_Px[  mneutrino[j]  ],
                                Particle_Py[  mneutrino[j]  ],
                                Particle_Pz[  mneutrino[j]  ],
                                Particle_E[   mneutrino[j]  ] );

              sum = temp1 + temp2;

              //W wektorze zapamiętujemy sumę czterowektorów e+ i ve czyli
              //tak naprawdę czteropęd W+
              W_p.push_back(sum);

              /////////////////////////////////////////////////////////////////
              //Boostujemy do ukladu spoczynkowego Higgsa
              if (Higgs == 0) continue;

              temp1.Boost(*boost_vector);
              temp2.Boost(*boost_vector);

              sum = temp1 + temp2;

              //Zapamiętujemy boostowane czterowektory w wektorze
              W_p_hf.push_back( sum );
          }
      }

      /************************************************************************/
      //Pętla po wszystkich mu- i vm~
      //będziemy z tego odtwarzać W-
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

              /////////////////////////////////////////////////////////////////
              //Boostujemy do ukladu spoczynkowego Higgsa
              if (Higgs == 0) continue;

              temp1.Boost(*boost_vector);
              temp2.Boost(*boost_vector);

              sum = temp1 + temp2;

              //Zapamiętujemy boostowane czterowektory w wektorze
              W_n_hf.push_back( sum );

          }
      }

      /************************************************************************/
      //Pętla po wszystkich e- i ve~
      //będziemy z tego odtwarzać W-
      for (int i = 0; i < elektorn.size(); i++)
      {
          for (int j = 0; j < aeneutrino.size(); j++)
          {
              temp1.SetPxPyPzE( Particle_Px[  elektorn[i]       ],
                                Particle_Py[  elektorn[i]       ],
                                Particle_Pz[  elektorn[i]       ],
                                Particle_E[   elektorn[i]       ] );
              temp2.SetPxPyPzE( Particle_Px[  aeneutrino[j] ],
                                Particle_Py[  aeneutrino[j] ],
                                Particle_Pz[  aeneutrino[j] ],
                                Particle_E[   aeneutrino[j] ] );

              sum = temp1 + temp2;
              //W wektorze zapamiętujemy sumę czterowektorów mu- i vm~ czyli
              //tak naprawdę czteropęd W-
              W_n.push_back(sum);

              /////////////////////////////////////////////////////////////////
              //Boostujemy do ukladu spoczynkowego Higgsa
              if (Higgs == 0) continue;

              temp1.Boost(*boost_vector);
              temp2.Boost(*boost_vector);

              sum = temp1 + temp2;

              //Zapamiętujemy boostowane czterowektory w wektorze
              W_n_hf.push_back( sum );

          }
      }

      /************************************************************************/
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

              /////////////////////////////////////////////////////////////////
              //Różnica phi między nimi (po booscie)
              phi = abs(W_n_hf[i].Phi() - W_p_hf[j].Phi());

              //Cięcie (normowanie) do Pi
              if (phi < TMath::Pi()) Deltaphi_ww_hf->Fill( phi );
              else Deltaphi_ww_hf->Fill( 2 * TMath::Pi() - phi );
          }
      }

      /***********************************************************************/
      //Pęd poprzeczny neutrin
        neutrino_momentum.SetPxPyPzE(0,0,0,0);

        for (int i = 0; i < aeneutrino.size(); i++)
        {
              temp1.SetPxPyPzE( Particle_Px[  aeneutrino[i] ],
                                Particle_Py[  aeneutrino[i] ],
                                Particle_Pz[  aeneutrino[i] ],
                                Particle_E[   aeneutrino[i] ] );

              neutrino_momentum += temp1;
              Pdg_id->Fill( Particle_PID[aeneutrino[i]] );
        }

        for (int i = 0; i < mneutrino.size(); i++)
        {
              temp1.SetPxPyPzE( Particle_Px[  mneutrino[i] ],
                                Particle_Py[  mneutrino[i] ],
                                Particle_Pz[  mneutrino[i] ],
                                Particle_E[   mneutrino[i] ] );

              neutrino_momentum += temp1;
              Pdg_id->Fill( Particle_PID[mneutrino[i]] );
        }

        for (int i = 0; i < eneutrino.size(); i++)
        {
              temp1.SetPxPyPzE( Particle_Px[  eneutrino[i] ],
                                Particle_Py[  eneutrino[i] ],
                                Particle_Pz[  eneutrino[i] ],
                                Particle_E[   eneutrino[i] ] );

              neutrino_momentum += temp1;
              Pdg_id->Fill( Particle_PID[eneutrino[i]] );
        }

        for (int i = 0; i < amneutrino.size(); i++)
        {
              temp1.SetPxPyPzE( Particle_Px[  amneutrino[i] ],
                                Particle_Py[  amneutrino[i] ],
                                Particle_Pz[  amneutrino[i] ],
                                Particle_E[   amneutrino[i] ] );

              neutrino_momentum += temp1;
              Pdg_id->Fill( Particle_PID[amneutrino[i]] );
        }

        nn_p->Fill(neutrino_momentum.Pt());

      /************************************************************************/
      //Pęd poprzeczny wszystkich mozliwych do detekcji czastek
      final_momentum.SetPxPyPzE(0,0,0,0);

      for (int i = 0; i < aelektorn.size(); i++)
      {
            temp1.SetPxPyPzE( Particle_Px[  aelektorn[i] ],
                              Particle_Py[  aelektorn[i] ],
                              Particle_Pz[  aelektorn[i] ],
                              Particle_E[   aelektorn[i] ] );

            final_momentum += temp1;
            Pdg_id->Fill( Particle_PID[aelektorn[i]] );
      }

      for (int i = 0; i < elektorn.size(); i++)
      {
            temp1.SetPxPyPzE( Particle_Px[  elektorn[i] ],
                              Particle_Py[  elektorn[i] ],
                              Particle_Pz[  elektorn[i] ],
                              Particle_E[   elektorn[i] ] );

            final_momentum += temp1;
            Pdg_id->Fill( Particle_PID[elektorn[i]] );
      }

      for (int i = 0; i < amion.size(); i++)
      {
            temp1.SetPxPyPzE( Particle_Px[  amion[i] ],
                              Particle_Py[  amion[i] ],
                              Particle_Pz[  amion[i] ],
                              Particle_E[   amion[i] ] );

            final_momentum += temp1;
            Pdg_id->Fill( Particle_PID[amion[i]] );
      }

      for (int i = 0; i < mion.size(); i++)
      {
            temp1.SetPxPyPzE( Particle_Px[  mion[i] ],
                              Particle_Py[  mion[i] ],
                              Particle_Pz[  mion[i] ],
                              Particle_E[   mion[i] ] );

            final_momentum += temp1;
            Pdg_id->Fill( Particle_PID[mion[i]] );
      }

      for (int i = 0; i < jet1.size(); i++)
      {
            temp1.SetPxPyPzE( Particle_Px[  jet1[i] ],
                              Particle_Py[  jet1[i] ],
                              Particle_Pz[  jet1[i] ],
                              Particle_E[   jet1[i] ] );

            final_momentum += temp1;
            Pdg_id->Fill( Particle_PID[jet1[i]] );
      }

      for (int i = 0; i < jet2.size(); i++)
      {
            temp1.SetPxPyPzE( Particle_Px[  jet2[i] ],
                              Particle_Py[  jet2[i] ],
                              Particle_Pz[  jet2[i] ],
                              Particle_E[   jet2[i] ] );

            final_momentum += temp1;
            Pdg_id->Fill( Particle_PID[jet2[i]] );
      }

      non_nn_p->Fill(final_momentum.Pt());

      /************************************************************************/
      //Pęd wszystkich cząstek (neutrina + rest)
      all_momentum->Fill((final_momentum + neutrino_momentum).Pt());

      /************************************************************************/
      //Obliczenie zetowej skladowej pędu neutrin
      M_h = 126.0; M_vv = 30.0;
      M_fix_2 = M_h * M_h - ll->E() * ll->E() - M_vv * M_vv + 2 * ll->Px() * (-(*ll + *jj)).Px() + 2 * ll->Py() * (-(*ll + *jj)).Py();
      M_fix_4_ = M_fix_2 * M_fix_2 - 4 * ll->E() * ll->E() * M_vv * M_vv - 4 * ll->E() * ll->E() * (-(*ll + *jj)).Pt() * (-(*ll + *jj)).Pt();

      if (M_fix_2 * M_fix_2 *ll->Pz() * ll->Pz() - M_fix_4_ * (ll->Pz() * ll->Pz() - ll->E() * ll->E()) < 0) continue;

      p_vv_z = (-(ll->Pz() * ll->Pz()) - sqrt(M_fix_2 * M_fix_2 *ll->Pz() * ll->Pz() - M_fix_4_ * (ll->Pz() * ll->Pz() - ll->E() * ll->E())))
      / (2 * (ll->Pz() * ll->Pz() - ll->E() * ll->E()));

      p_vv_z_h->Fill(p_vv_z);

      p_vv_z = (-(ll->Pz() * ll->Pz()) + sqrt(M_fix_2 * M_fix_2 *ll->Pz() * ll->Pz() - M_fix_4_ * (ll->Pz() * ll->Pz() - ll->E() * ll->E())))
      / (2 * (ll->Pz() * ll->Pz() - ll->E() * ll->E()));

      p_vv_z_h->Fill(p_vv_z);
   }

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
   non_nn_p->Write();
   all_momentum->Write();

   Pdg_id->Write();

   p_vv_z_h->Write();

   f->Close();
}

void LHEF(char* input, char* output)
{
  //Dla ROOT 6.x
  class LHEF a(input);

  //Dla ROOT 5.x
  //LHEF a;

  a.Loop(output);
  new TBrowser;
}
