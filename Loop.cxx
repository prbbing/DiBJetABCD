/*
 *      Loop.cxx
 *
 *      Copyright 2010 Sergei Chekanov <chakanau@hep.anl.gov> ANL
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */



#include "Ana.h"
#include "Global.h"
#include "SystemOfUnits.h"
#include "Histo.h"
#include "TSpline.h"
#include<iostream>
#include<fstream>
#include<stdlib.h>
using namespace std;

extern Global glob;
extern Histo  h;

// Event loop.
void Ana::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   cout << " -> number of entries =" << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t n=0; n<nentries; n++) {
      
      if (glob.nev>=glob.MaxEvents) break;

      Long64_t ientry = LoadTree(n);
      if (ientry < 0) break;
      nb = fChain->GetEntry(n);   nbytes += nb;

       if (glob.nev<=10 &&
       glob.nev<=100 && (glob.nev%10) == 0 ||
       glob.nev<=10000 && (glob.nev%10000) == 0  ||
       glob.nev>=1000000 && (glob.nev%1000000) == 0 ||
       glob.nev>=100000000 && (glob.nev%10000000) == 0 )  {
       cout << "Events= " << glob.nev << endl; 
        }

       // counter
       glob.nev++;
//      if(NPV<15)continue;
      // skip bad events
      h.debug->Fill("good events",1.);
      if (CutEvent(ientry) <0) continue;
      h.debug->Fill("tileOnly",1.);

      // weight in MC (data=1)
      double weight=1.0;
//      cout<<glob.m_run<<endl;

//=====================================================================
//btagging scale factor and uncertainty
//=====================================================================

//if(jet_pt->at(0)<430)continue;
//if(jet_pt->at(1)<430)continue;

double btagSF=weight, btagErr=weight;
double bJES = 0.026;//FIXME
// 85% btagCut 77% 0.4803 70% 0.7110 60% 0.8867
double btagCut=0.175848;//0.175848;0.645925,77;0.824427,70;0.934906,60//repro mv2c10 85%
#if montecarlo
      weight=mcEventWeight*weight_pileup;
//if((glob.nev%10000) == 0)cout<<"mc weight: "<<weight<<endl;
/*
      if(weight_btag_fix_85->size()==0)continue;
      double weight_btag = weight_btag_fix_85->at(0);
      double weight_btag_err = 0;
      double bJES = 0.026;//FIXME
      for (int i = 1; i< weight_btag_fix_85->size();i++)
         weight_btag_err += (weight_btag_fix_85->at(i)/weight_btag-1)*(weight_btag_fix_85->at(i)/weight_btag-1);
      weight_btag_err = sqrt(weight_btag_err/(weight_btag_fix_85->size()-1));
      if(weight_btag==-999)weight_btag=1;
      btagSF=weight*weight_btag;
      btagErr=weight_btag_err;
//if((glob.nev%10000) == 0)cout<<"weight: "<<weight<<endl;
      if(glob.m_run<0&&glob.dobtagSF){
        weight=weight*weight_btag;
//if((glob.nev%10000) == 0)cout<<"btag SF weight: "<<weight<<endl;
      }
      if(glob.m_run<0&&glob.dobtagSF&&glob.dobtagSF1up){
        weight=weight*weight_btag*(1+weight_btag_err);
//if((glob.nev%10000) == 0)cout<<"btag SF 1up weight: "<<weight<<endl;
      }
      if(glob.m_run<0&&glob.dobtagSF&&glob.dobtagSF1down){
        weight=weight*weight_btag*(1-weight_btag_err);
//if((glob.nev%10000) == 0)cout<<"btag SF 1down weight: "<<weight<<endl;
      }
*/
#endif
//=====================================================================
//get all the jets in the event
//=====================================================================
      vector<LParticle> selected;
      vector<LParticle> alljets;
      for(unsigned int i = 0; i < njets; i++){
        if (jet_pt->at(i)<60)continue;
        double pt=jet_pt->at(i);
        #if montecarlo
          if (glob.dobJES1up)pt=pt*(1+bJES);
          if (glob.dobJES1down)pt=pt*(1-bJES);
        #endif
        double phi=jet_phi->at(i);
        double eta=jet_eta->at(i);

        TLorentzVector l;
        l.SetPtEtaPhiE(pt,eta,phi,jet_E->at(i));
        LParticle p;
        p.SetP(l);
        p.SetType( 0  );
        p.SetStatus( jet_clean_passLooseBad->at(i) );
        p.SetParent(1 );
        p.SetCharge( 0 );
       // p.SetParameter( jet_chfrac->at(i)   );
       // p.SetParameter( jet_tilecal->at(i)   );
        p.SetParameter( jet_Timing->at(i)   );
        p.SetParameter( jet_MV2c10->at(i)   );
        p.SetParameter( jet_SV1->at(i)   );
        p.SetParameter( jet_IP3D->at(i)   );
        p.SetParameter( jet_LArQuality->at(i)   );
        p.SetParameter( jet_HECQuality->at(i)   );
        p.SetParameter( jet_NegativeE->at(i)   );
        #if montecarlo
          if(glob.dofavorplot)p.SetParameter( jet_HadronConeExclTruthLabelID->at(i)  );
        #else
          p.SetParameter( 0  );
        #endif
        p.SetParameter( 0  );
        p.SetParameter( 0  );

        alljets.push_back(p);
        selected.push_back(p);
      } // end loop over jets 
      
      //=====================================================================
      //calculate jets parameter
      //=====================================================================
      std::sort(selected.begin(), selected.end(), greater<LParticle>() ) ;

      // good jets //
      if (selected.size()<=1) continue;

      std::sort(selected.begin(), selected.end(), greater<LParticle>() ) ;

      LParticle p1=selected.at(0);
      LParticle p2=selected.at(1);
      
      double pt1=p1.GetP().Perp();
      double pt2=p2.GetP().Perp();
      double eta1=p1.GetP().Eta();
      double eta2=p2.GetP().Eta();
      TLorentzVector LP1=p1.GetP();
      TLorentzVector LP2=p2.GetP();
      TLorentzVector PP=LP1+LP2;
      //cout << PP.M() << endl;
      double mass_jj=PP.M();
      double yStar=(p1.GetP().Rapidity()-p2.GetP().Rapidity())/2;
     //===================================================================================
     //cuts
     //===================================================================================
      if(p1.GetStatus()==0||p2.GetStatus()==0)continue;
      if (find(passedTriggers->begin(), passedTriggers->end(), "HLT_j380") == passedTriggers->end())continue;
      if (pt1<430) continue;
      if (pt2<80) continue;
      if( abs(eta1) > glob.ETA_CUT|| abs(eta2) > glob.ETA_CUT) continue;
      if (abs(yStar)>0.8) continue;
      if (mass_jj<1200)continue;

      //=====================================================================
      //fill histograms
      //=====================================================================
      h.NPV->Fill(NPV,weight);
      vector<float> effs = SF(pt1,eta1,pt2,eta2);
      float weightb = weight*effs[2];
      float weightb_err =  weight*effs[3];
      h.jetjetmass_b[0]->Fill(mass_jj,weightb); 
      h.jetjetmass_b[0]->SetBinError(h.jetjetmass_b[0]->FindBin(mass_jj),sqrt(pow(h.jetjetmass_b[0]->GetBinError(h.jetjetmass_b[0]->FindBin(mass_jj)),2)+pow(weightb_err,2)));
      float weightbb = weight*effs[0];
      float weightbb_err =  weight*effs[1];
      h.jetjetmass_bb[0]->Fill(mass_jj,weightbb); 
      h.jetjetmass_bb[0]->SetBinError(h.jetjetmass_bb[0]->FindBin(mass_jj),sqrt(pow(h.jetjetmass_bb[0]->GetBinError(h.jetjetmass_bb[0]->FindBin(mass_jj)),2)+pow(weightbb_err,2)));
    glob.TotalEvents++;

  }//end loop over all events
}//end loop method

vector<float> Ana::SF(float pt1,float eta1,float pt2,float eta2) {
  vector<float> effs;
  TFile f("/users/bingxuan.liu/DiBjet/2016Frozen/analysis_sys_dijet_Scale/efficiency_mc_highmass.root");
  TH2F* leadingEff = (TH2F*)f.Get("leadingJetEff");
  TH2F* secondEff = (TH2F*)f.Get("secondJetEff");
  TH2F* secondEff_con = (TH2F*)f.Get("secondJetEff_con");
  effs.push_back(leadingEff->GetBinContent(leadingEff->GetXaxis()->FindBin(pt1), leadingEff->GetYaxis()->FindBin(abs(eta1)))*secondEff_con->GetBinContent(secondEff_con->GetXaxis()->FindBin(pt2), secondEff_con->GetYaxis()->FindBin(abs(eta2))));
  effs.push_back(sqrt(pow(leadingEff->GetBinContent(leadingEff->GetXaxis()->FindBin(pt1), leadingEff->GetYaxis()->FindBin(abs(eta1))),2)*pow(secondEff_con->GetBinError(secondEff_con->GetXaxis()->FindBin(pt2), secondEff_con->GetYaxis()->FindBin(abs(eta2))),2)+pow(secondEff_con->GetBinContent(secondEff_con->GetXaxis()->FindBin(pt2), secondEff_con->GetYaxis()->FindBin(abs(eta2))),2)*pow(leadingEff->GetBinError(leadingEff->GetXaxis()->FindBin(pt1), leadingEff->GetYaxis()->FindBin(abs(eta1))),2)));
  effs.push_back(leadingEff->GetBinContent(leadingEff->GetXaxis()->FindBin(pt1), leadingEff->GetYaxis()->FindBin(abs(eta1))) + secondEff->GetBinContent(secondEff->GetXaxis()->FindBin(pt2), secondEff->GetYaxis()->FindBin(abs(eta2))) - leadingEff->GetBinContent(leadingEff->GetXaxis()->FindBin(pt1), leadingEff->GetYaxis()->FindBin(abs(eta1)))*secondEff_con->GetBinContent(secondEff_con->GetXaxis()->FindBin(pt2), secondEff_con->GetYaxis()->FindBin(abs(eta2))));
  effs.push_back(sqrt(pow(leadingEff->GetBinError(leadingEff->GetXaxis()->FindBin(pt1), leadingEff->GetYaxis()->FindBin(abs(eta1))),2) + pow(secondEff->GetBinError(secondEff->GetXaxis()->FindBin(pt2), secondEff->GetYaxis()->FindBin(abs(eta2))),2) + pow(leadingEff->GetBinContent(leadingEff->GetXaxis()->FindBin(pt1), leadingEff->GetYaxis()->FindBin(abs(eta1))),2)*pow(secondEff_con->GetBinError(secondEff_con->GetXaxis()->FindBin(pt2), secondEff_con->GetYaxis()->FindBin(abs(eta2))),2)+pow(secondEff_con->GetBinContent(secondEff_con->GetXaxis()->FindBin(pt2), secondEff_con->GetYaxis()->FindBin(abs(eta2))),2)*pow(leadingEff->GetBinError(leadingEff->GetXaxis()->FindBin(pt1), leadingEff->GetYaxis()->FindBin(abs(eta1))),2)));
  return effs;
}
