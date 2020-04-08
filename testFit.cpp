#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMath.h"


#include "include/particlefromtree.hpp"
#include "include/SmearingTool.h"
#include "include/KinematicTools.h"
#include "include/Fitter.h"
//#include "include/particle.hpp"

#include <iostream>

/*
   To compile outside ROOT env do:

   g++ `root-config --cflags --libs` testRead.cpp -o testRead.exe

   Inside ROOT env

   .L testRead.cpp++
   testRead()
*/


void testFit(){

  TFile* tfile = TFile::Open("../BKstTauMu.root");
  TTree* ttree = (TTree*) tfile->Get("BKstTauMuTuple/MCDecayTree");

  Particle< Double_t > muplus ( "muplus" , ttree );
  Particle< Double_t > kplus  ( "Kplus"  , ttree );
  Particle< Double_t > piminus( "piminus", ttree );

  Particle< Double_t > tau( "tauminus", ttree );

  Particle< Double_t > pi1( "piplus"  , ttree );
  Particle< Double_t > pi2( "piminus0", ttree );
  Particle< Double_t > pi3( "piminus1", ttree );

  Particle< Double_t > nu( "nu_tau" , ttree );

  Vertex< Double_t >   pv( "B0_TRUEORIGINVERTEX"   , ttree );
  Vertex< Double_t >   sv( "B0_TRUEENDVERTEX"      , ttree );
  Vertex< Double_t >   tv( "tauminus_TRUEENDVERTEX", ttree );


  TH1D* hist_tau_mass = new TH1D("hist_tau_mass","",200,0.,20000.);
  TH1D* hist_tau_mass_post = new TH1D("hist_tau_mass_post","",200,0.,20000.);
  TH1D* hist_B_mass  = new TH1D("hist_B_mass" ,"",200,0.,20000.);
  TH1D* hist_B_angle = new TH1D("hist_B_angle" ,"",100,-0.1,0.3);
  TH1D* hist_B_mass_post  = new TH1D("hist_B_mass_post" ,"",200,0.,20000.);
  TH1D* hist_B_angle_post  = new TH1D("hist_B_angle_post" ,"",100,-0.1,0.300);
  TH1D* hist_tau_angle_post  = new TH1D("hist_tau_angle_post" ,"",100,-0.1,0.300);
  TH1D* hist_tau_angle  = new TH1D("hist_tau_angle" ,"",100,-0.1,0.300);

  TH1D* hist_B_Z_post = new TH1D("hist_B_Z","",100,-1.,1.);

  const bool applyPositionSmearing = true;
  const bool applyMomentumSmearing = false;

  SmearingTool smearing;

  for ( Long64_t i = 0; i < ttree->GetEntries(); i++ ){

    Fitter fitter;

    ttree->GetEntry(i);

    TLorentzVector vpi1  = pi1.getVec();
    TLorentzVector vpi2  = pi2.getVec();
    TLorentzVector vpi3  = pi3.getVec();
    TLorentzVector vkaon = kplus.getVec();
    TLorentzVector vpion = piminus.getVec();
    TLorentzVector vmuon = muplus.getVec();
    TLorentzVector vnu   = nu.getVec();

    if ( applyMomentumSmearing ){
      smearing.smearMomentum( vpi1 );
      smearing.smearMomentum( vpi2 );
      smearing.smearMomentum( vpi3 );
      smearing.smearMomentum( vkaon );
      smearing.smearMomentum( vpion );
      smearing.smearMomentum( vmuon );
    }

    TVector3 vsv = sv.getPos();
    TVector3 vpv = pv.getPos();
    TVector3 vtv = tv.getPos();

    if ( applyPositionSmearing ){
      smearing.smearPosition( vsv );
      //smearing.smearPosition( vtv );
    }

    fitter.setOrigin( vpv );
    fitter.setBVertex( vsv );
    fitter.setTauVertex( vtv );
    fitter.setPion1( vpi1 );
    fitter.setPion2( vpi2 );
    fitter.setPion3( vpi3 );
    fitter.setPion( vpion );
    fitter.setKaon( vkaon );
    fitter.setMuon( vmuon );
    fitter.setNeutrino( 0., 0., 0. );
    fitter.setChain();

  //  fitter.m_parent_vtx.Print();

    hist_tau_mass->Fill( fitter.getTauMass() );
    hist_B_mass->Fill( fitter.getBMass() );
    hist_B_angle->Fill( fitter.getBPointingAngle() );
    hist_tau_angle->Fill( fitter.getTauPointingAngle() );

    fitter.setNeutrino( 0., 0., 10. );
    fitter.setChain();

    fitter.progressive_kalman( 50 );

    hist_tau_mass_post->Fill( fitter.getTauMass() );
    hist_B_angle_post->Fill( fitter.getBPointingAngle() );
    hist_B_mass_post->Fill( fitter.getBMass() );
    hist_B_Z_post->Fill( fitter.getBDeltaZ() );
    hist_tau_angle_post->Fill( fitter.getTauPointingAngle() );

  //  fitter.m_parent_vtx.Print();
  }
  //TFile* tfile = TFile::Open("BKstTauTau.root");
  //TTree* ttree = (TTree*) tfile->Get("BKstTauTauTuple/MCDecayTree");

  // Defining 4-momenta of particles using hpp file.
  Particle< Double_t > Kplus_M(    "Kplus",    ttree );
  Particle< Double_t > piminus_M(  "piminus",  ttree );
  Particle< Double_t > muplus_M(   "muplus",   ttree );
  Particle< Double_t > piplus_M(   "piplus",   ttree );
  Particle< Double_t > piminus0_M( "piminus0", ttree );
  Particle< Double_t > piminus1_M( "piminus1", ttree );

  // hpp file does not currently work for end vertex positions, so old method being used.
  // Defining B origin vertex components
  double B0_TRUEORIGINVERTEX_X, B0_TRUEORIGINVERTEX_Y, B0_TRUEORIGINVERTEX_Z;
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_X",&B0_TRUEORIGINVERTEX_X);
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_Y",&B0_TRUEORIGINVERTEX_Y);
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_Z",&B0_TRUEORIGINVERTEX_Z);
  // Defining tauminus origin vertex components
  double tauminus_TRUEORIGINVERTEX_X, tauminus_TRUEORIGINVERTEX_Y, tauminus_TRUEORIGINVERTEX_Z;
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_X", &tauminus_TRUEORIGINVERTEX_X);
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_Y", &tauminus_TRUEORIGINVERTEX_Y);
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_Z", &tauminus_TRUEORIGINVERTEX_Z);
  // Defining tauminus end vertex components
  double tauminus_TRUEENDVERTEX_X, tauminus_TRUEENDVERTEX_Y, tauminus_TRUEENDVERTEX_Z;
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_X", &tauminus_TRUEENDVERTEX_X);
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_Y", &tauminus_TRUEENDVERTEX_Y);
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_Z", &tauminus_TRUEENDVERTEX_Z);


  TH1D* hist_Bmass   = new TH1D("B Meson Mass","Reconstructed B Meson Mass",200,0.,20000.);
  hist_Bmass->SetTitle(";Mass, MeV c^{-2};Number of Entries");
  //hist_Bmass->SetFillColor(kRed);
  TH1D* hist_taumass = new TH1D("Tauon Mass","Reconstructed Tau Lepton Meson Mass",200,0.,20000.);
  //hist_taumass->SetTitle(";Mass, MeV c^{-2};Number of Entries");
  //hist_taumass->SetFillColor(kRed);
  TH1D* hist_p9_X    = new TH1D("","",100,-0.1,0.300);
    hist_p9_X->SetTitle(";#theta(B);Number of Entries");
  
  //hist_p9_X->SetFillColor(kRed);
  TH1D* hist_p9_Y    = new TH1D("","",100,-0.1,0.300);
      hist_p9_Y->SetTitle(";#theta(#tau);Number of Entries");
  //hist_p9_Y->SetFillColor(kRed);
  TH1D* hist_p9_Z    = new TH1D("","",100,-0.1,0.300);
  
  //hist_p9_Z->SetFillColor(kRed);
  TH1D* hist_p9_E    = new TH1D("","",100,-0.1,0.300);
  //hist_p9_E->SetFillColor(kRed);

  TRandom3* rand = new TRandom3();
  
  for (Long64_t i = 0; i < ttree->GetEntries(); i++ ) {
    ttree->GetEntry(i);
    TLorentzVector Kplus_Mvec    = Kplus_M.getVec();
    TLorentzVector piminus_Mvec  = piminus_M.getVec();
    TLorentzVector muplus_Mvec   = muplus_M.getVec();
    TLorentzVector piplus_Mvec   = piplus_M.getVec();
    TLorentzVector piminus0_Mvec = piminus0_M.getVec();
    TLorentzVector piminus1_Mvec = piminus1_M.getVec();
    
    // Particle 2-4 and 6-8 mass in tau or mode change mass!
    double s_Kplus    = rand->Gaus(493.677, 0.013);      // from B J et al (2012) particle listings
    double s_piminus  = rand->Gaus(139.570018, 0.00035); // from C Amsler et al (2008) particle listings
    double s_muplus   =  105.65837;                       // from Beringer J et al (particle data group) 2012 particle summary
    double s_piplus   = rand->Gaus(139.570018, 0.00035);
    double s_piminus0 = rand->Gaus(139.570018, 0.00035);
    double s_piminus1 = rand->Gaus(139.570018, 0.00035);
      
    double a = 0;
    double b = 0.04;
    double c = 0.2;
      
    // Adding errors to particles 2-4
    double Kplus_abs             = sqrt(pow(Kplus_Mvec.X(), 2) + pow(Kplus_Mvec.Y(), 2) + pow(Kplus_Mvec.Z(), 2));
    double Kplus_abs_SMEARED     = rand->Gaus(Kplus_abs, a*Kplus_abs);
    double K_smear_factor        = Kplus_abs_SMEARED/Kplus_abs;
    //double Kplus_X               = rand->Gaus(Kplus_Mvec.X(), K_smear_factor);
    //double Kplus_Y               = rand->Gaus(Kplus_Mvec.Y(), K_smear_factor);
    //double Kplus_Z               = rand->Gaus(Kplus_Mvec.Z(), K_smear_factor);
    double Kplus_X               = (Kplus_Mvec.X()*K_smear_factor);
    double Kplus_Y               = (Kplus_Mvec.Y()*K_smear_factor);
    double Kplus_Z               = (Kplus_Mvec.Z()*K_smear_factor);
    double Kplus_Mvec_E          = sqrt((s_Kplus)*(s_Kplus) + (Kplus_X)*(Kplus_X) + (Kplus_Y)*(Kplus_Y) + (Kplus_Z)*(Kplus_Z));
      
    
    double piminus_abs           = sqrt(pow(piminus_Mvec.X(), 2) + pow(piminus_Mvec.Y(), 2) + pow(piminus_Mvec.Z(), 2));
    double piminus_abs_SMEARED   = rand->Gaus(piminus_abs, a*piminus_abs);
    double piminus_smear_factor  = piminus_abs_SMEARED/piminus_abs;
    double piminus_X             = (piminus_Mvec.X()*piminus_smear_factor);
    double piminus_Y             = (piminus_Mvec.Y()*piminus_smear_factor);
    double piminus_Z             = (piminus_Mvec.Z()*piminus_smear_factor);
    double piminus_Mvec_E        = sqrt((s_piminus)*(s_piminus) + (piminus_X)*(piminus_X) + (piminus_Y)*(piminus_Y) + (piminus_Z)*(piminus_Z));
    
    double muplus_abs            = sqrt(pow(muplus_Mvec.X(), 2) + pow(muplus_Mvec.Y(), 2) + pow(muplus_Mvec.Z(), 2));
    double muplus_abs_SMEARED    = rand->Gaus(muplus_abs, a*muplus_abs);
    double muplus_smear_factor   = muplus_abs_SMEARED/muplus_abs;
    double muplus_X              = (muplus_Mvec.X()*muplus_smear_factor);
    double muplus_Y              = (muplus_Mvec.Y()*muplus_smear_factor);
    double muplus_Z              = (muplus_Mvec.Z()*muplus_smear_factor);
    double muplus_Mvec_E         = sqrt((s_muplus)*(s_muplus) + (muplus_X)*(muplus_X) + (muplus_Y)*(muplus_Y) + (muplus_Z)*(muplus_Z));
    
    // Adding errors to particles 6-8
    double piplus_abs            = sqrt(piplus_Mvec.X()*piplus_Mvec.X() + piplus_Mvec.Y()*piplus_Mvec.Y() + piplus_Mvec.Z()*piplus_Mvec.Z());
    double piplus_abs_smeared    = rand->Gaus(piplus_abs, a*piplus_abs);
    double piplus_smear_factor   = piplus_abs_smeared/piplus_abs;
    double piplus_X              = (piplus_Mvec.X()*piplus_smear_factor);
    double piplus_Y              = (piplus_Mvec.Y()*piplus_smear_factor);
    double piplus_Z              = (piplus_Mvec.Z()*piplus_smear_factor);
    double piplus_Mvec_E         = sqrt((s_piplus)*(s_piplus) + (piplus_X)*(piplus_X) + (piplus_Y)*(piplus_Y) + (piplus_Z)*(piplus_Z));
      
    double piminus0_abs          = sqrt(piminus0_Mvec.X()*piminus0_Mvec.X() + piminus0_Mvec.Y()*piminus0_Mvec.Y() + piminus0_Mvec.Z()*piminus0_Mvec.Z());
    double piminus0_abs_smeared  = rand->Gaus(piminus0_abs, a*piminus0_abs);
    double piminus0_smear_factor = piminus0_abs_smeared/piminus0_abs;
    double piminus0_X            = (piminus0_Mvec.X()*piminus0_smear_factor);
    double piminus0_Y            = (piminus0_Mvec.Y()*piminus0_smear_factor);
    double piminus0_Z            = (piminus0_Mvec.Z()*piminus0_smear_factor);
    double piminus0_Mvec_E       = sqrt((s_piplus)*(s_piplus) + (piminus0_X)*(piminus0_X) + (piminus0_Y)*(piminus0_Y) + (piminus0_Z)*(piminus0_Z));
      
    double piminus1_abs          = sqrt(piminus1_Mvec.X()*piminus1_Mvec.X() + piminus1_Mvec.Y()*piminus1_Mvec.Y() + piminus1_Mvec.Z()*piminus1_Mvec.Z());
    double piminus1_abs_smeared  = rand->Gaus(piminus1_abs, a*piminus1_abs);
    double piminus1_smear_factor = piminus1_abs_smeared/piminus1_abs;
    double piminus1_X            = (piminus1_Mvec.X()*piminus1_smear_factor);
    double piminus1_Y            = (piminus1_Mvec.Y()*piminus1_smear_factor);
    double piminus1_Z            = (piminus1_Mvec.Z()*piminus1_smear_factor);
    double piminus1_Mvec_E       = sqrt((s_piplus)*(s_piplus) + (piminus1_X)*(piminus1_X) + (piminus1_Y)*(piminus1_Y) + (piminus1_Z)*(piminus1_Z));
      
    // Errors on verticies
    double B0_VERTEX_X = rand->Gaus(B0_TRUEORIGINVERTEX_X, 0);
    double B0_VERTEX_Y = rand->Gaus(B0_TRUEORIGINVERTEX_Y, 0);
    double B0_VERTEX_Z = B0_TRUEORIGINVERTEX_Z;
    
    double tauminus_VERTEX_X = rand->Gaus(tauminus_TRUEORIGINVERTEX_X, b);
    double tauminus_VERTEX_Y = rand->Gaus(tauminus_TRUEORIGINVERTEX_Y, b);
    double tauminus_VERTEX_Z = rand->Gaus(tauminus_TRUEORIGINVERTEX_Z, c);

    double tauminus_ENDVERTEX_X = rand->Gaus(tauminus_TRUEENDVERTEX_X, b);
    double tauminus_ENDVERTEX_Y = rand->Gaus(tauminus_TRUEENDVERTEX_Y, b);
    double tauminus_ENDVERTEX_Z = rand->Gaus(tauminus_TRUEENDVERTEX_Z, c);

    // Main Program.
    // Reconstruction calculations.
    double s5unit_X = tauminus_ENDVERTEX_X - tauminus_VERTEX_X;
    double s5unit_Y = tauminus_ENDVERTEX_Y - tauminus_VERTEX_Y;
    double s5unit_Z = tauminus_ENDVERTEX_Z - tauminus_VERTEX_Z;

    double s1unit_X = tauminus_VERTEX_X - B0_VERTEX_X;
    double s1unit_Y = tauminus_VERTEX_Y - B0_VERTEX_Y;
   double s1unit_Z = tauminus_VERTEX_Z - B0_VERTEX_Z;
    
    double xi = ((Kplus_X + piminus_X + muplus_X)*s1unit_Y - (Kplus_Y + piminus_Y + muplus_Y)*s1unit_X)/((s5unit_Y*s1unit_X) - (s5unit_X*s1unit_Y));
                    
    double p9_X = xi*s5unit_X - (piplus_X + piminus0_X + piminus1_X);
    double p9_Y = xi*s5unit_Y - (piplus_Y + piminus0_Y + piminus1_Y);
    double p9_Z = xi*s5unit_Z - (piplus_Z + piminus0_Z + piminus1_Z);
    double p9_E = sqrt(pow(p9_X, 2) + pow(p9_Y, 2) + pow(p9_Z, 2));
        
    double s5_0 = pow(p9_E + piplus_Mvec_E + piminus0_Mvec_E + piminus1_Mvec_E, 2);
    double s5_1 = pow(p9_X + piplus_X + piminus0_X + piminus1_X, 2) + pow(p9_Y + piplus_Y + piminus0_Y + piminus1_Y, 2) + pow(p9_Z + piplus_Z + piminus0_Z + piminus1_Z, 2);
                    
    double s5 = sqrt(s5_0-s5_1);

    double s4_0 = pow((sqrt(pow(s5, 2) + (pow(s5unit_X, 2) + pow(s5unit_Y, 2) + pow(s5unit_Z, 2))*pow(xi, 2)) + Kplus_Mvec_E + piminus_Mvec_E + muplus_Mvec_E), 2);

    double s4_1 = pow((s5unit_X*xi + Kplus_X + piminus_X + muplus_X), 2) + pow((s5unit_Y*xi + Kplus_Y + piminus_Y + muplus_Y), 2) + pow((s5unit_Z*xi + Kplus_Z + piminus_Z + muplus_Z), 2);

    double s1 = sqrt((s4_0-s4_1));
      
      

    double xi_N = ((Kplus_Mvec.X()+piminus_Mvec.X()+muplus_Mvec.X())*s1unit_Y - (Kplus_Mvec.Y()+piminus_Mvec.Y()+muplus_Mvec.Y())*s1unit_X) / ((s5unit_Y*s1unit_X)-(s5unit_X*s1unit_Y));
       double gamma_N = ((Kplus_Mvec.Y()+piminus_Mvec.Y()+muplus_Mvec.Y()) / s1unit_Y) + (s5unit_Y/s1unit_Y) * xi;
      TVector3 part1(gamma_N*s1unit_X, gamma_N*s1unit_Y, gamma_N*s1unit_Z);
      TVector3 part5(xi_N*s5unit_X,    xi_N*s5unit_Y,    xi_N*s5unit_Z);
       TVector3 parpar(tauminus_VERTEX_X, tauminus_VERTEX_Y, tauminus_VERTEX_Z);
      TVector3 parpar2(tauminus_ENDVERTEX_X, tauminus_ENDVERTEX_Y, tauminus_ENDVERTEX_Z);
      
      double x_arr_1_abs = sqrt(s1unit_X *s1unit_X  + s1unit_Y*s1unit_Y +s1unit_Z*s1unit_Z);
      double part1_abs = sqrt(gamma_N*s1unit_X*gamma_N*s1unit_X + gamma_N*s1unit_Y*gamma_N*s1unit_Y + gamma_N*s1unit_Z*gamma_N*s1unit_Z);
      
      double x_arr_2_abs = sqrt(s5unit_X*s5unit_X + s5unit_Y*s5unit_Y + s5unit_Z*s5unit_Z);
      double part5_abs   = sqrt(xi_N*s5unit_X*xi_N*s5unit_X + xi_N*s5unit_Y*xi_N*s5unit_Y + xi_N*s5unit_Z*xi_N*s5unit_Z);
      
      TVector3 g0((tauminus_VERTEX_X - B0_VERTEX_X)*part1_abs - part1.X()*x_arr_1_abs, (tauminus_VERTEX_Y - B0_VERTEX_Y)*part1_abs - part1.X()*x_arr_1_abs, (tauminus_VERTEX_Z - B0_VERTEX_Z)*part1_abs - part1.X()*x_arr_1_abs);
     // cout<<s1<<"\n";
      
        double angle = g0.Angle(parpar);
      
      TVector3 sgo (s5unit_X*part5_abs - part5.X()*x_arr_2_abs, s5unit_Y*part5_abs - part5.Y()*x_arr_2_abs, s5unit_Z*part5_abs - part5.Z()*x_arr_2_abs );
      
      double angle2 = sgo.Angle(parpar2);
    
    // Histogram being filled with b meson mass values.
    hist_Bmass->Fill(s1);
    hist_taumass->Fill(s5);
    
      
    hist_p9_X->Fill(angle);
    hist_p9_Y->Fill(angle2);
    
      
      hist_p9_Z->Fill(p9_Z);
    hist_p9_E->Fill(p9_E);
  }
  TCanvas* can = new TCanvas("can","can");
  can->Divide(2,2);
  can->cd(1);
  
  //hist_B_mass_post->Draw();
  double m_nE = hist_B_mass_post->GetMean();
  double stdDev_B = hist_B_mass_post->GetStdDev();
  std::cout<< "mean: " << m_nE << " "<<"std"<<std::endl;
  hist_B_mass_post->Draw();
    hist_B_mass_post->SetLineStyle(2);
  hist_Bmass->Draw("SAME");
  hist_Bmass->SetLineColor(kBlack);
  hist_B_mass_post->SetLineColor( kRed );
  hist_B_mass_post->SetXTitle("Mass, MeV c^{-2}");
  hist_B_mass_post->SetYTitle("Number of Entries");

  can->cd(2);
  hist_tau_mass_post->Draw();
    hist_tau_mass_post->SetLineStyle(2);
  hist_taumass->Draw("Same");
  hist_tau_mass_post->SetLineColor( kRed );
hist_taumass->SetLineColor( kBlack);
  hist_tau_mass_post->SetXTitle("Mass, MeV c^{-2}");
  hist_tau_mass_post->SetYTitle("Number of Entries");

  can->cd(3);
  hist_p9_X->Draw();
  hist_B_angle_post->Draw("SAME");
    hist_B_angle_post->SetLineStyle(2);
  //hist_p9_X->Draw("SAME");
  hist_B_angle_post->SetLineColor( kRed );
    hist_p9_X->SetLineColor(kBlack);
  hist_B_angle_post->SetXTitle("#theta(B)");
  hist_B_angle_post->SetYTitle("Number of Entries");

  can->cd(4);
  hist_p9_Y->Draw();
  hist_tau_angle_post->Draw("SAME");
  hist_tau_angle_post->SetLineStyle(2);
  //  hist_p9_Y->Draw("SAME");
  hist_tau_angle_post->SetLineColor( kRed );
hist_p9_Y->SetLineColor(kBlack);
  hist_tau_angle_post->SetXTitle("#theta(#tau)");
  hist_tau_angle_post->SetYTitle("Number of Entries");
  
  return;
}


