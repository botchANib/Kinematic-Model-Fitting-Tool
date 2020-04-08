#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMath.h"

#include "include/particlefromtree.hpp"
#include "include/KinematicTools.h"
#include "include/SmearingTool.h"

#include <iostream>

/*
   To compile outside ROOT env do:

   g++ `root-config --cflags --libs` testRead.cpp -o testRead.exe

   Inside ROOT env

   .L testRead.cpp++
   testRead()
*/


void testRead(){

  TFile* tfile = TFile::Open("../BKstTauMu.root");
  TTree* ttree = (TTree*) tfile->Get("BKstTauMuTuple/MCDecayTree");
  
  std::cout << "Entries = " << ttree->GetEntries() << std::endl;
  

  Particle< Double_t > muplus ( "muplus" , ttree );
  Particle< Double_t > kplus  ( "Kplus"  , ttree );
  Particle< Double_t > piminus( "piminus", ttree );

  Particle< Double_t > tau( "tauminus", ttree );

  Particle< Double_t > pi1("piplus"  , ttree );
  Particle< Double_t > pi2("piminus0", ttree );
  Particle< Double_t > pi3("piminus1", ttree );

  Vertex< Double_t >   pv( "B0_TRUEORIGINVERTEX"   , ttree );
  Vertex< Double_t >   sv( "B0_TRUEENDVERTEX"      , ttree );
  Vertex< Double_t >   tv( "tauminus_TRUEENDVERTEX", ttree );

  TH1D* hist_mB      = new TH1D("hist_mB"     ,"hist_mB"     ,100,0,20);
  TH1D* hist_mCorr   = new TH1D("hist_mCorr"  ,"hist_mCor"   ,100,0,20);
  TH1D* hist_mTau    = new TH1D("hist_mTau"   ,"hist_mTau"   ,100,0,4);
  TH1D* hist_mTauRec = new TH1D("hist_mTauRec","hist_mTauRec",100,0,4);
  TH1D* hist_mBRec   = new TH1D("hist_mBRec"  ,"hist_mBRec"  ,100,0,20);
  TH1D* hist_mKstar  = new TH1D("hist_mKstar" ,"hist_mKstar" ,100,0,20);
  TH1D* hist_mTauTau = new TH1D("hist_mTauTau","hist_mTauTau",100,0,4);
  TH1D* hist_mBTau   = new TH1D("hist_mBTau"  ,"hist_mBTau"  ,100,0,20);
  
  const bool applyPositionSmearing = true;
  const bool applyMomentumSmearing = false;

  SmearingTool smearing;
  //smearing.setPositionSmear( 0.01, 0.01, 0.00 );

  for ( Long64_t i = 0; i < ttree->GetEntries(); i++ ){

    ttree->GetEntry(i);

    TLorentzVector vpi1 = pi1.getVec();
    TLorentzVector vpi2 = pi2.getVec();
    TLorentzVector vpi3 = pi3.getVec();

    if ( applyMomentumSmearing ){
      smearing.smearMomentum( vpi1 );
      smearing.smearMomentum( vpi2 );
      smearing.smearMomentum( vpi3 );
    }

    TLorentzVector vtau  = vpi1 + vpi2 + vpi3;

    TLorentzVector vkaon = kplus.getVec();
    TLorentzVector vpion = piminus.getVec();
    TLorentzVector vmuon = muplus.getVec();

    if ( applyMomentumSmearing ){
      smearing.smearMomentum( vkaon );
      smearing.smearMomentum( vpion );
      smearing.smearMomentum( vmuon );
    }

    TLorentzVector vkstar = vkaon  + vpion;
    TLorentzVector vtag   = vkstar + vmuon;
    TLorentzVector vrec   = vtau   + vtag;

    hist_mKstar->Fill( vkstar.M() );
    hist_mB  ->Fill( vrec.M() );
    hist_mTau->Fill( vtau.M() );

    TVector3 vsv = sv.getPos();
    TVector3 vpv = pv.getPos();
    TVector3 vtv = tv.getPos();

    if ( applyPositionSmearing ){
      smearing.smearPosition( vsv );
      smearing.smearPosition( vtv );
    }

    TVector3 dirB   = (vsv - vpv).Unit();
    TVector3 dirTau = (vtv - vsv).Unit();

    TVector3 ptau   = dirTau * KinematicTools::scale( vtag, dirB, dirTau );
    TVector3 pnu    = ptau - vtau.Vect();

    TLorentzVector vnu( pnu, pnu.Mag() );

    TLorentzVector vtaurec = vtau + vnu;

    
    hist_mTauRec->Fill( vtaurec.M() );

    const double mcorr = KinematicTools::mcor( vrec, dirB );
    hist_mCorr  ->Fill( mcorr );

    TLorentzVector vBrec = vtaurec + vtag;

    hist_mBRec->Fill( vBrec.M() );
    

    /*
    //TLorentzVector vNuTau  = KinematicTools::nuFromTauDir( vtau, dirTau, +1);
    TLorentzVector vNuTauP = KinematicTools::nuFromBDir( vtau, vtag, dirB, +1);
    TLorentzVector vNuTauM = KinematicTools::nuFromBDir( vtau, vtag, dirB, -1);

    TLorentzVector vTauTau = vtau    + vNuTauP;
    TLorentzVector vBTauP  = vTauTau + vtag;
    TLorentzVector vBTauM  = vNuTauM + vtau + vtag;

    if ( vNuTauP.P() > 0 ){
      hist_mTauTau->Fill( vTauTau.M() );
      hist_mBTau  ->Fill( vBTauP.M() );
    }

    hist_mBscatter->Fill( vNuTauP.P() > 0 ? vBTauP.M() : 0,
                          vNuTauM.P() > 0 ? vBTauP.M() : 0, )
    */
  }

  TCanvas* can = new TCanvas("can","can");
  can->Divide(2,2);

  can->cd(1);
  hist_mKstar->SetXTitle("m(K#pi)");
  hist_mKstar->Draw();

  can->cd(2);
  hist_mTauRec->Draw();
  hist_mTauRec->SetXTitle("m(#tau)");

  hist_mTauRec->SetMaximum(1000.0);
  hist_mTauRec->SetLineColor( kRed );
  hist_mTau->Draw("SAME");

//  hist_mTauTau->Draw("SAME");
//  hist_mTauTau->SetLineColor( kRed );

  can->cd(3);
  hist_mBRec->Draw();
  hist_mBRec->SetXTitle("m(B)");

  hist_mBRec->SetMaximum(1000.0);
  hist_mBRec->SetLineColor( kRed );
  hist_mB->Draw("SAME");
  hist_mCorr->Draw("SAME");
  hist_mCorr->SetLineColor( kBlue );
//  hist_mBTau->Draw("SAME");

//  hist_mBTau->SetLineColor( kRed );

  can->SaveAs("testRead.pdf");

  return;
}

int main( void ){

  testRead();

  return 0;
}
