/*
 * utility.h
 *
 *  Created on: Dec 4, 2018
 *      Author: suyong
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "correction.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include <string>
using namespace ROOT;
using namespace ROOT::VecOps;
using namespace std;

using floats =  ROOT::VecOps::RVec<float>;
using doubles =  ROOT::VecOps::RVec<double>;
using ints =  ROOT::VecOps::RVec<int>;
using bools = ROOT::VecOps::RVec<bool>;
using uchars = ROOT::VecOps::RVec<unsigned char>;

using FourVector = ROOT::Math::PtEtaPhiMVector;
using FourVectorVec = std::vector<FourVector>;
using FourVectorRVec = ROOT::VecOps::RVec<FourVector>;


struct hist1dinfo
{
  ROOT::RDF::TH1DModel hmodel;
  std::string varname;
  std::string weightname;
  std::string mincutstep;
} ;

//for 2D histograms
struct hist2dinfo
{
    ROOT::RDF::TH2DModel hmodel;
    std::string varname1;
    std::string varname2;
    std::string weightname;
    std::string mincutstep;
} ;

struct varinfo
{
  std::string varname;
  std::string vardefinition;
  std::string mincutstep;
};

struct cutinfo
{
  std::string cutdefinition;
  std::string idx;
};

// generates vectors of 4 vectors given vectors of pt, eta, phi, mass
FourVectorVec generate_4vec(floats &pt, floats &eta, floats &phi, floats &mass);
FourVectorVec generate_4vec_2(float pt, float eta, float phi, float mass);
FourVectorVec generate_4vec_3(floats &pt, floats &eta, floats &phi);
FourVectorVec generate_4vec_4(float pt, float eta, float phi);
FourVectorVec generate_4vec_5(floats pt, double eta, floats phi);
FourVectorRVec generate_4Rvec(floats &pt, floats &eta, floats &phi, floats &mass);
FourVectorRVec generate_4Rvec_2(double pt, double eta, double phi, double mass);
FourVectorRVec generate_4Rvec_4(float pt, float eta, float phi);
float FourVecPt(FourVector fv);
float FourVecEta(FourVector fv);
float FourVecPhi(FourVector fv);
float FourVecMass(FourVector fv);
float FourVecEnergy(FourVector fv);
float FourVecPx(FourVector fv);
float FourVecPy(FourVector fv);
float FourVecPz(FourVector fv);
// floats FourVecEnergy(floats &pt, floats &eta, floats &phi, floats &mass);
FourVectorVec genmet4vec(float met_pt, float met_phi);

// return a vector size equal to length of x all filled with evWeight value
floats weightv(floats &x, float evWeight);
floats sphericity(FourVectorVec &p);
double foxwolframmoment(int l, FourVectorVec &p, int minj=0, int maxj=-1);
floats btvcorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string name, std::string syst, ints &hadflav, floats &etas,floats &pts,  floats &btags);
floats btv_case1(std::unique_ptr<correction::CorrectionSet> &cset, std::string type, std::string sys,std::string wp, ints &hadflav, floats &etas, floats &pts );
floats btv_case2(std::unique_ptr<correction::CorrectionSet> &cset, std::string type, std::string sys,std::string wp, ints &hadflav, floats &etas, floats &pts );
floats btv_casetest(std::unique_ptr<correction::CorrectionSet>& cset, std::string type1, std::string sys, std::string wp, ints& hadflav, floats& etas, floats& pts);
floats muoncorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string type, std::string year,floats &etas,floats &pts, std::string sys);
float pucorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string name, std::string syst, float ntruepileup);
floats chi2(float smtop_mass, float smw_mass, float lfvtop_mass);
floats top_reconstruction_whad(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons);
floats top_reco_products(FourVectorVec &jets, FourVectorVec &muons, floats topreco);
ints good_idx(ints good); 
float calculate_deltaEta( FourVector &p1, FourVector &p2);
float calculate_deltaPhi( FourVector &p1, FourVector &p2);
float calculate_deltaR( FourVector &p1, FourVector &p2);
float calculate_invMass( FourVector &p1, FourVector &p2);
FourVector sum_4vec( FourVector &p1, FourVector &p2);
FourVector sum_4Rvec( FourVector p1, FourVector p2);
floats concatenate( double first, double second);
floats sort_discriminant( floats discr, floats obj );
FourVector select_leadingvec( FourVectorVec &v );
floats PrintVector(floats myvector);
floats w_reconstruction (FourVectorVec &jets);
floats compute_DR (FourVectorVec &muons, ints goodMuons_charge);
floats H_reconstruction (FourVectorVec &jets, FourVectorVec &bjets);
//floats dimuon (FourVectorVec &muons, ints Muon_charge);
bools isLfromtopHiggs (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx);//, floats &leptonpt);
FourVectorVec isLfromtoporHiggs (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx, FourVectorVec &genleptons4vecs);
bools isLfromW (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx);
ints MotherfromL (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx);
bools isWfromHiggs (ints &genpdgId, ints &wgenidx, ints &midx, floats &genmass);
bools isWfromtop (ints &genpdgId, ints &wgenidx, ints &midx, floats &genmass);
bools isLfromtopantitop (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx);//, floats &leptonpt);
FourVectorVec isLfromtoporantitop (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx, FourVectorVec &genleptons4vecs);
bools isWfromtoponly (ints &genpdgId, ints &wgenidx, ints &midx, floats &genmass);
bools isWfromantitop (ints &genpdgId, ints &wgenidx, ints &midx, floats &genmass);
bools CheckOverlaps(FourVectorVec &jets, FourVectorVec &leps);
RVec< RVec<int> > distinct_comb_3(int n);
FourVectorRVec Threejets (FourVectorRVec jet1, FourVectorRVec jet2, FourVectorRVec jet3, int jet_number, int bjet_number);
FourVectorRVec Twojets (FourVectorRVec jet1, FourVectorRVec jet2, FourVectorRVec jet3, floats &jet1_btag, floats &jet2_btag, floats &jet3_btag, int jet_number, int bjet_number);
bools isPfromtopHiggs (ints &genpdgId, ints &midx, ints &lgenidx);
bools isPfromHiggs (ints &genpdgId, ints &midx, ints &lgenidx);
bools isPfromtop (ints &genpdgId, ints &midx, ints &lgenidx);
bools isPfromtoponly (ints &genpdgId, ints &midx, ints &lgenidx);
bools isPfromantitop (ints &genpdgId, ints &midx, ints &lgenidx);
bools isdR04 (floats &eta1, floats &eta2, floats &phi1, floats &phi2);
bools isJetFromgJet (floats &jet_pt, floats &gjet_pt, floats &gjetTp_pt, ints &jet_gidx);//, floats &jet_eta, floats &jet_phi, floats &jet_mass, ints &jet_jetid, floats &gjet_eta, floats &gjet_phi, floats &gjet_mass)
bools isJetFromgJet_2 (double &jet_pt, floats &gjet_pt, floats &gjetTp_pt, int &jet_gidx);
floats DR_1 (floats &eta1, floats &eta2, floats &phi1, floats &phi2);
floats DR_2 (double &eta1, floats &eta2, double &phi1, floats &phi2);
int maxDR (double &eta1, floats &eta2, double &phi1, floats &phi2);
ints minDR_1 (floats &eta1, floats &eta2, floats &phi1, floats &phi2);
ints minDR_2 (floats &eta1, floats &eta2, floats &phi1, floats &phi2);
int minDR_3 (double &eta1, floats &eta2, double &phi1, floats &phi2);
ints minDR_4 (floats &eta1, floats &eta2, floats &phi1, floats &phi2, floats &jet_pt, double &bjet_pt);
float minDR_5 (floats &eta1, floats &eta2, floats &phi1, floats &phi2);
float minDR_6 (double &eta1, floats &eta2, double &phi1, floats &phi2);
int Max (floats &variable);
ints mass_W (floats jet_pt, floats jet_eta, floats jet_phi, floats jet_mass);
ints mass_W_offshell (floats jet_pt, floats jet_eta, floats jet_phi, floats jet_mass);
bool false_mass_dilepton (FourVectorVec &leptons);
floats ratio_1 (floats &lepton_pt, floats &jet_pt);
floats ratio_2 (double &lepton_pt, floats &jet_pt);
float minratio_1 (floats &lepton_pt, floats &jet_pt);
float minratio_2 (double &lepton_pt, floats &jet_pt);
FourVectorRVec removebjet(FourVectorRVec jet4vecs, floats jetpt, double bjetpt);
RVec< RVec<int> > distinct_comb_2(int nl, int nj);
FourVectorRVec sum_two_4Rvecs (FourVectorRVec lepton, FourVectorRVec bjet);
RVec< RVec<int> > distinct_comb_2_bis(int nj);
FourVectorRVec reconstruction_neutrino(double leptontoppt, double leptontopeta, double leptontopphi, double leptontopmass, double bjetpt, double bjeteta, double bjetphi, double bjetmass, double leptonHiggspt, double leptonHiggseta, double leptonHiggsphi, double leptonHiggsmass, float MET_px, float MET_py);
double eta_neutrino_Higgs(double leptontopeta, double leptontopphi, double leptonHiggseta, double leptonHiggsphi, floats &nutopeta, floats &nutopphi, floats &nuHiggsphi);
double eta_neutrino_Higgs_2(double leptontopeta, double leptontopphi, double leptonHiggseta, double leptonHiggsphi, double nutopeta, double nutopphi, double nuHiggsphi);
// double Basic_Mt2_332_Calculator::mT2Fcn::operator()(const std::vector<double>& par);
// double mt2_332_Sq(FourVectorRVec visA, FourVectorRVec visB, float ptmiss, float phimiss, const double mEachInvisible);
// double mt2_332(FourVectorRVec visA, FourVectorRVec visB, float ptmiss, float phimiss, const double mEachInvisible);
// floats solve(float metpx, float metpy, floats &be, floats &bpx, floats &bpy, floats &bpz, floats &je, floats &jpx, floats &jpy, floats &jpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &mj, floatks &ml, floats &pnue, floats &pnux, floats &pnuy, floats &pnuz, floats &nupt, floats &nueta, floats &nuphi, floats &num, floats &wmasses);
floats solve(float metpx, float metpy, floats &be, floats &bpx, floats &bpy, floats &bpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &ml, floats &pnue, floats &pnux, floats &pnuy, floats &pnuz, floats &nupt, floats &nueta, floats &nuphi, floats &num, floats &wtmass, floats &watmass);
void quartic(vector<double> polx, vector<double> *pnux, int& cubic_single_root_cmplx);
void cubic(vector<double> polx, vector<double> *pnux);
void quadratic(vector<double> polx, vector<double> *pnux);
// int algebraic_pz_1(floats &be, floats &bpx, floats &bpy, floats &bpz, floats &je, floats &jpx, floats &jpy, floats &jpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &mj, floats &ml, double pnux, double pnuy, double* pnuz, double epsilon);
// int algebraic_pz_2(floats &be, floats &bpx, floats &bpy, floats &bpz, floats &je, floats &jpx, floats &jpy, floats &jpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &mj, floats &ml, double pnux, double pnuy, double* pnuz, double epsilon, float wmass);
int algebraic_pz_1(floats &be, floats &bpx, floats &bpy, floats &bpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &ml, double pnux, double pnuy, double* pnuz, double epsilon, float wmass);
int algebraic_pz_2(floats &be, floats &bpx, floats &bpy, floats &bpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &ml, double pnux, double pnuy, double* pnuz, double epsilon, float wmass);
double evalterm1(vector<double> *a1, double pnux, double pnuy);
double evalterm2(vector<double> *a2, double pnux, double pnuy);
RVec< RVec<int> > distinct_comb_1_even(int nn);
RVec< RVec<int> > distinct_comb_1_odd(int nn);
RVec< RVec<int> > distinct_comb_2_even_odd(int nn);
FourVectorRVec Leptonneutrino (FourVectorRVec lepton, FourVectorRVec neutrino, int nn);
floats deltaeta(double &leta, floats &neta);
floats deltaphi(double &lphi, floats &nphi);
floats deltaR(double &leta, floats &neta, double &lphi, floats &nphi);
floats transverse_mass(double &lmass, floats &te, double &lpt, floats &npt, floats &dphi);
FourVectorRVec Leptonneutrinobjet (FourVectorRVec bjet, FourVectorRVec leptonneutrino, int nn);
floats transverse_mass_bjet(double &bmass, floats &lnmass, floats &te, double &jpt, floats &lnpt, floats &dphi);

#endif /* UTILITY_H_ */
