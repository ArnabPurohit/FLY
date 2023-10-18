/*
 * utility.cpp
 *
 *  Created on: Dec 4, 2018
 *      Author: suyong
 */
#include "utility.h"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/Math.h"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "correction.h"
// #include "Mt2/Basic_Mt2_332_Calculator.h"
// #include "Mt2/Mt2Units.h"
// #include "Minuit2/MnUserParameterState.h"
// #include "Minuit2/MnMigrad.h"
// #include "Minuit2/MnSimplex.h"
// #include "Minuit2/MnScan.h"
// #include "Minuit2/FunctionMinimum.h"
#include <iostream>
#include <algorithm>
using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;

using floats =  ROOT::VecOps::RVec<float>;
using ints =  ROOT::VecOps::RVec<int>;
using FourVector = ROOT::Math::PtEtaPhiMVector;
using FourVectorRVec = ROOT::VecOps::RVec<FourVector>;
using FourVectorVec = vector<FourVector>;
// using ROOT::Minuit2::MnUserParameters;
// using ROOT::Minuit2::MnMigrad;
// using ROOT::Minuit2::MnSimplex;
// using ROOT::Minuit2::MnScan;
// using ROOT::Minuit2::FunctionMinimum;

// Utility function to generate fourvector objects for thigs that pass selections

using namespace std;

FourVectorVec generate_4vec(floats &pt, floats &eta, floats &phi, floats &mass)
{
  const int nsize = pt.size();
  FourVectorVec fourvecs;
  fourvecs.reserve(nsize);
  for (auto i=0; i<nsize; i++)
  {
    // cout<<"pt = "<<pt<<", eta = "<<eta<<", phi = "<<phi<<", mass = "<<mass<<endl;
    fourvecs.emplace_back(pt[i], eta[i], phi[i], fabs(mass[i]));
  }
  
  return fourvecs;
}

FourVectorVec generate_4vec_2(float pt, float eta, float phi, float mass)
{
  FourVectorVec fourvecs;
  fourvecs.reserve(1);
  // cout<<"pt = "<<pt<<", eta = "<<eta<<", phi = "<<phi<<", mass = "<<mass<<endl;
  fourvecs.emplace_back(pt, eta, phi, fabs(mass));
  return fourvecs;
}

FourVectorVec generate_4vec_3(floats &pt, floats &eta, floats &phi)
{
  const int nsize = pt.size();
  FourVectorVec fourvecs;
  fourvecs.reserve(nsize);
  for (auto i=0; i<nsize; i++)
  {
    // cout<<"pt = "<<pt<<", eta = "<<eta<<", phi = "<<phi<<endl;
    fourvecs.emplace_back(pt[i], eta[i], phi[i], phi[i]-phi[i]);
  }
  return fourvecs;
}

FourVectorVec generate_4vec_4(float pt, float eta, float phi)
{
  FourVectorVec fourvecs;
  fourvecs.reserve(1);
  // cout<<"pt = "<<pt<<", eta = "<<eta<<", phi = "<<phi<<", mass = "<<mass<<endl;
  fourvecs.emplace_back(pt, eta, phi, phi-phi);
  return fourvecs;
}

FourVectorVec generate_4vec_5(floats pt, double eta, floats phi)
{
  FourVectorVec fourvecs;
  fourvecs.reserve(1);
  // cout<<"pt = "<<pt<<", eta = "<<eta<<", phi = "<<phi<<", mass = "<<mass<<endl;
  fourvecs.emplace_back(pt[0], (float) eta, phi[0], phi[0]-phi[0]);
  return fourvecs;
}

FourVectorRVec generate_4Rvec(floats &pt, floats &eta, floats &phi, floats &mass)
{
    const int nsize = pt.size();
    FourVectorRVec fourvecs;
    fourvecs.reserve(nsize);
    for (auto i=0; i<nsize; i++)
    {
        // cout<<"pt = "<<pt<<", eta = "<<eta<<", phi = "<<phi<<", mass = "<<mass<<endl;
        fourvecs.emplace_back(pt[i], eta[i], phi[i], fabs(mass[i]));
    }

    return fourvecs;
}

FourVectorRVec generate_4Rvec_2(double pt, double eta, double phi, double mass)
{
    FourVectorRVec fourvecs;
    fourvecs.reserve(1);
    // cout<<"pt = "<<pt<<", eta = "<<eta<<", phi = "<<phi<<", mass = "<<mass<<endl;
    fourvecs.emplace_back(pt, eta, phi, fabs(mass));
    return fourvecs;
}

FourVectorRVec generate_4Rvec_4(double pt, double eta, double phi)
{
  FourVectorRVec fourvecs;
  fourvecs.reserve(1);
  // cout<<"pt = "<<pt<<", eta = "<<eta<<", phi = "<<phi<<", mass = "<<mass<<endl;
  fourvecs.emplace_back(pt, eta, phi, phi-phi);
  return fourvecs;
}

float FourVecPt(FourVector fv)
{
  return fv.Pt();
}

float FourVecEta(FourVector fv)
{
  return fv.Eta();
}

float FourVecPhi(FourVector fv)
{
  return fv.Phi();
}

float FourVecMass(FourVector fv)
{
  return fv.M();
}

float FourVecEnergy(FourVector fv)
{
  return fv.E();
}

float FourVecPx(FourVector fv)
{
  return fv.Px();
}

float FourVecPy(FourVector fv)
{
  return fv.Py();
}

float FourVecPz(FourVector fv)
{
  return fv.Pz();
}

// floats FourVecEnergy(floats &pt, floats &eta, floats &phi, floats &mass)
// {
//   const int nsize = pt.size();
//   FourVectorRVec fourvecs;
//   fourvecs.reserve(nsize);
//   for (auto i=0; i<nsize; i++)
//   {
//     fourvecs.emplace_back(pt[i], eta[i], phi[i], fabs(mass[i]));
//   }

//   return fourvecs.E();
// }

FourVectorVec genmet4vec(float met_pt, float met_phi)
{
        FourVectorVec vecs;
        float met_px = met_pt*cos(met_phi);
        float met_py = met_pt*sin(met_phi);
        for(int i = -500; i <= 500; i+=50){
                FourVector metfourvec;
                metfourvec.SetPxPyPzE(met_px, met_py, i, met_pt);
                vecs.emplace_back(metfourvec);
        }
        return vecs;
}

floats weightv(floats &x, float evWeight)
{
  const int nsize = x.size();
  floats weightvector(nsize, evWeight);
  return weightvector;
}

floats sphericity(FourVectorVec &p)
{
  TMatrixDSym NormMomTensor(3);
  
  NormMomTensor = 0.0;
  double p2sum = 0.0;
  for (auto x: p)
  {
    p2sum += x.P2();
    double mom[3] = {x.Px(), x.Py(), x.Pz()};
    for (int irow=0; irow<3; irow++)
    {
      for (int icol=irow; icol<3; icol++)
      {
	NormMomTensor(irow, icol) += mom[irow] * mom[icol];
      }
    }
  }
  NormMomTensor *= (1.0/p2sum);
  TVectorT<double> Qrev;
  NormMomTensor.EigenVectors(Qrev);
  floats Q(3);
  for (auto i=0; i<3; i++) Q[i] = Qrev[2-i];

  return Q;
}

double foxwolframmoment(int l, FourVectorVec &p, int minj, int maxj)
{   // PRD 87, 073014 (2013)
  double answer = 0.0;
  
  double ptsum=0.0;
  
  if (maxj==-1) // process everything
  {
    maxj = p.size();
  }
  //for (auto x: p)
  for (auto i=minj; i<maxj; i++)
  {
    auto x = p[i];
    ptsum += x.Pt();
    //for (auto y: p)
    for (auto j=minj; j<maxj; j++)
    {
      auto y = p[j];
      double wij = x.Pt() * y.Pt();
      double cosdOmega = x.Vect().Dot(y.Vect()) / (x.P() * y.P());
      if (cosdOmega>1.0) cosdOmega=1.0;
      if (cosdOmega<-1.0) cosdOmega=-1.0;
      answer += wij * ROOT::Math::legendre(l, cosdOmega);
    }
  }
  answer /= ptsum*ptsum;
  if (fabs(answer)>1.0) cout << "FW>1 " << answer << endl;
  return answer;
}

floats chi2(float smtop_mass, float smw_mass, float lfvtop_mass)
{
  floats out;
  // Theory values
  //        const float MT_LFV = 172.5;
  //        const float MT_SM = 172.5;
  //        const float MW = 80.4;
  //        const float WT_LFV = 1.41;
  //        const float WT_SM = 1.41;
  //        const float WW = 2.085;
  
  // Resolution applied values
  const float MT_LFV = 150.5;
  const float MT_SM = 165.2;
  const float MW = 80.8;
  const float WT_LFV = 17.8;
  const float WT_SM = 21.3;
  const float WW = 11.71;    
  
  float chi2_SMTop = pow((MT_SM-smtop_mass)/WT_SM, 2);
  float chi2_SMW = pow((MW-smw_mass)/WW, 2);
  float chi2_LFVTop = pow((MT_LFV-lfvtop_mass)/WT_LFV, 2);
  float chi2 = chi2_SMTop + chi2_SMW + chi2_LFVTop;
  
  out.emplace_back(chi2);
  out.emplace_back(chi2_SMTop);
  out.emplace_back(chi2_SMW);
  out.emplace_back(chi2_LFVTop);
  
  return out;
}

//floats top_reconstruction_whad(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons, FourVectorVec &taus){
floats top_reconstruction_whad(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons)
{
  floats out;
  
  float LFVtop_mass, SMW_mass, SMtop_mass;
  float X_LFVtop, X_SMW, X_SMtop;
  float X_min=9999999999, X_min_LFVtop_mass=-1, X_min_SMW_mass=-1, X_min_SMtop_mass=-1;
  float X_min_LFVtop=999999999, X_min_SMW=999999999, X_min_SMtop=999999999;
  float c_idx=-1, wj1_idx=-1, wj2_idx=-1;
  
  // Mass and Width of top and w boson
  const float MT_LFV = 150.5;
  const float MT_SM = 165.2;
  const float MW = 80.8;
  const float WT_LFV = 17.8;
  const float WT_SM = 21.3;
  const float WW = 11.71;    
  
  // U or C jet
  for(unsigned int j1 = 0; j1<jets.size(); j1++)
  {
    //LFVtop_mass = (jets[j1]+taus[0]+muons[0]).M();
    LFVtop_mass = (jets[j1]+muons[0]).M();
    X_LFVtop = pow((MT_LFV-LFVtop_mass)/WT_LFV,2);
    // Jet1 from W
    for(unsigned int j2 = 0; j2<jets.size(); j2++)
    {
      if(jets[j2].Pt() == bjets[0].Pt() || jets[j2].Pt() == jets[j1].Pt()) continue;
      // Jet2 from W
      for(unsigned int j3 = 0; j3<jets.size(); j3++)
      {
	if(jets[j3].Pt() == jets[j2].Pt() || jets[j3].Pt() == bjets[0].Pt() || jets[j3].Pt() == jets[j1].Pt()) continue;
	SMW_mass = (jets[j2]+jets[j3]).M();
	X_SMW = pow((MW-SMW_mass)/WW,2);
	SMtop_mass = (bjets[0]+jets[j2]+jets[j3]).M();
	X_SMtop = pow((MT_SM-SMtop_mass)/WT_SM,2);
	//min error
	if (X_LFVtop + X_SMW + X_SMtop < X_min)
	{
	  X_min = X_LFVtop + X_SMW + X_SMtop;
	  X_min_LFVtop = X_LFVtop;
	  X_min_SMW = X_SMW;
	  X_min_SMtop = X_SMtop;
	  X_min_LFVtop_mass = LFVtop_mass;
	  X_min_SMW_mass = SMW_mass;
	  X_min_SMtop_mass = SMtop_mass;
	  c_idx = float(j1);
	  wj1_idx = float(j2);
	  wj2_idx = float(j3);
	}
      }
    }
  }
  out.push_back(X_min);               // 0
  out.push_back(X_min_LFVtop_mass);   // 1
  out.push_back(X_min_SMW_mass);      // 2
  out.push_back(X_min_SMtop_mass);    // 3
  out.push_back(c_idx);               // 4
  out.push_back(wj1_idx);             // 5
  out.push_back(wj2_idx);             // 6
  out.push_back(X_min_LFVtop);        // 7
  out.push_back(X_min_SMW);           // 8
  out.push_back(X_min_SMtop);         // 9

  return out;
}

//floats top_reco_products(FourVectorVec &jets, FourVectorVec &muons, FourVectorVec &taus, floats topreco){
floats top_reco_products(FourVectorVec &jets, FourVectorVec &muons, floats topreco)
{
  floats out;
  int j_idx = topreco[4];
  int wjet1_idx = topreco[5];
  int wjet2_idx = topreco[6];
  
  FourVector lfvjet = jets[j_idx];
  FourVector wjet1 = jets[wjet1_idx];
  FourVector wjet2 = jets[wjet2_idx];
  //FourVector tau = taus[0];
  FourVector muon = muons[0];
  
  float wqq_dEta = wjet1.Eta() - wjet2.Eta();
  float wqq_dPhi = ROOT::Math::VectorUtil::DeltaPhi(wjet1, wjet2);
  float wqq_dR = ROOT::Math::VectorUtil::DeltaR(wjet1, wjet2);
  
  float lfvjmu_dEta = lfvjet.Eta() - muon.Eta();
  float lfvjmu_dPhi = ROOT::Math::VectorUtil::DeltaPhi(lfvjet, muon);
  float lfvjmu_dR = ROOT::Math::VectorUtil::DeltaR(lfvjet, muon);
  float lfvjmu_mass = ROOT::Math::VectorUtil::InvariantMass(lfvjet, muon);
  
  /*float lfvjtau_dEta = lfvjet.Eta() - tau.Eta();
    float lfvjtau_dPhi = ROOT::Math::VectorUtil::DeltaPhi(lfvjet, tau);
    float lfvjtau_dR = ROOT::Math::VectorUtil::DeltaR(lfvjet, tau);
    float lfvjtau_mass = ROOT::Math::VectorUtil::InvariantMass(lfvjet, tau);
    
    FourVector mutau = muon + tau;
    float lfvjmutau_dEta = lfvjet.Eta() - mutau.Eta();
    float lfvjmutau_dPhi = ROOT::Math::VectorUtil::DeltaPhi(lfvjet, mutau);
    float lfvjmutau_dR = ROOT::Math::VectorUtil::DeltaR(lfvjet, mutau);
    float lfvjmutau_mass = ROOT::Math::VectorUtil::InvariantMass(lfvjet, mutau);
    
  */
  out.emplace_back(wqq_dEta);         //0
  out.emplace_back(wqq_dPhi);         //1
  out.emplace_back(wqq_dR);           //2
  out.emplace_back(lfvjmu_dEta);      //3
  out.emplace_back(lfvjmu_dPhi);      //4
  out.emplace_back(lfvjmu_dR);        //5
  out.emplace_back(lfvjmu_mass);      //6
  //out.emplace_back(lfvjtau_dEta);     //7
  //out.emplace_back(lfvjtau_dPhi);     //8
  //out.emplace_back(lfvjtau_dR);       //9
  //out.emplace_back(lfvjtau_mass);     //10
  //out.emplace_back(lfvjmutau_dEta);   //11
  //out.emplace_back(lfvjmutau_dPhi);   //12
  //out.emplace_back(lfvjmutau_dR);     //13
  //out.emplace_back(lfvjmutau_mass);   //14
  
  return out;
}

floats btvcorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string type, std::string sys, ints &hadflav,floats &etas,floats &pts,   floats &btags)
{
	floats scalefactors;
	auto nvecs = pts.size();
	scalefactors.reserve(nvecs);
	for (auto i=0; i<nvecs; i++)
	{
		//std::cout << "sys: " << sys << ", hadflav: " << hadflav[i] << ", etas: " << fabs(float(etas[i])) << ", pts: " << float(pts[i]) <<" btag discriminator : " << float(btags[i]) << '\n';
		//for 2018UL working points l: 0.0494, m: 0.2770, t: 0.7264
		//if (btags[i]>0.7264){
		float sfi = cset->at(type)->evaluate({sys, int(hadflav[i]), fabs(float(etas[i])), float(pts[i]), float(btags[i])});
		scalefactors.emplace_back(sfi);
		cout<<" jets central scale factors == "<< sfi << endl;
		//}
	}
	return scalefactors;
}
//muonID_SF
floats muoncorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string type, std::string year, floats &etas, floats &pts, std::string sys)
{
    floats sf_muon;
    auto nvecs = pts.size();
    //cout<<"NVECS====== == "<< nvecs << endl;
    sf_muon.reserve(nvecs);
    
    if (etas.size() != nvecs) {
        throw std::invalid_argument("etas and pts vectors must have the same size!");
    }
    
    for (auto i=0; i<nvecs; i++)
    {
        //std::cout << "year: " << year  << ", etas: " << fabs(float(etas[i])) << ", pts: " << float(pts[i]) << " sys : " << sys << '\n';
        
        if (pts[i] < 0 ) {
            throw std::invalid_argument("Invalid value of pT detected!");
        }
        
        if (fabs(float(etas[i])) > 2.5) {
            throw std::invalid_argument("Invalid value of eta detected!");
        }
		//if (pts[i] < 15 || fabs(float(etas[i])) > 2.4) {
        //sf_muon.emplace_back(-1);
    	//	}
        float sfm = cset->at(type)->evaluate({year, fabs(float(etas[i])), float(pts[i]), sys});
        sf_muon.emplace_back(sfm);

    }
	//cout<<"SF MUON === " << sf_muon << endl;
    return sf_muon;
}

// case 1: Evaluate central SFs 
//fixedWP correction with mujets (here medium WP)
// evaluate('systematic', 'working_point', 'flavor', 'abseta', 'pt')
floats btv_case1(std::unique_ptr<correction::CorrectionSet>& cset, std::string type, std::string sys, std::string wp, ints& hadflav, floats& etas, floats& pts)
{
    floats scalefactors_case1;
    const auto nvecs = pts.size();
    //cout << "NVECS======== " << nvecs << endl;
    scalefactors_case1.reserve(nvecs);
    const auto abs_etas = [etas]() {
        floats res;
        res.reserve(etas.size());
        std::transform(etas.begin(), etas.end(), std::back_inserter(res), [](const auto& e) { return std::fabs(e); });
        return res;
    }();
    const auto cast_pts = [pts]() {
        floats res;
        res.reserve(pts.size());
        std::transform(pts.begin(), pts.end(), std::back_inserter(res), [](const auto& p) { return static_cast<float>(p); });
        return res;
    }();
  	for (auto i = 0; i < nvecs; i++) {
        //std::cout << "sys: " << sys << ", wp: " << wp << ", hadflav: " << hadflav[i] << ", etas: " << abs_etas[i] << ", pts: " << cast_pts[i] << '\n';
        if (hadflav[i] != 0) {
            // const auto bc_jets = cset->at("deepJet_mujets")->evaluate({sys, wp, hadflav[i], abs_etas[i], cast_pts[i]});
            const auto bc_jets = cset->at("deepJet_comb")->evaluate({sys, wp, hadflav[i], abs_etas[i], cast_pts[i]});
            scalefactors_case1.emplace_back(bc_jets);
            //std::cout << "\njet SFs from deepJe_mujets at medium WP\n";
            //std::cout << "SF b/c jets : " << bc_jets << '\n';
        } else{ 
           const auto bc_jets = cset->at("deepJet_incl")->evaluate({sys, wp, hadflav[i], abs_etas[i], cast_pts[i]});
           scalefactors_case1.emplace_back(bc_jets);
            //std::cout << "\njet SFs from deepJet_incl at medium WP\n";
            //std::cout << "SF light jets : " << bc_jets << '\n';
		}
        //
		
    }

    return scalefactors_case1;
}

// case 2: Evaluate varied SFs 
//fixedWP correction uncertainty (here tight WP and comb SF)
// evaluate('systematic', 'working_point', 'flavor', 'abseta', 'pt')
floats btv_case2(std::unique_ptr<correction::CorrectionSet>& cset, std::string type, std::string sys, std::string wp, ints& hadflav, floats& etas, floats& pts)
{
    floats scalefactors_case2;
    const auto nvecs = pts.size();
    //cout << "NVECS======== " << nvecs << endl;
    scalefactors_case2.reserve(nvecs);
    const auto abs_etas = [etas]() {
        floats res;
        res.reserve(etas.size());
        std::transform(etas.begin(), etas.end(), std::back_inserter(res), [](const auto& e) { return std::fabs(e); });
        return res;
    }();
    const auto cast_pts = [pts]() {
        floats res;
        res.reserve(pts.size());
        std::transform(pts.begin(), pts.end(), std::back_inserter(res), [](const auto& p) { return static_cast<float>(p); });
        return res;
    }();
  for (auto i = 0; i < nvecs; i++) {

		float bweight = 1.0;
		 
        if (hadflav[i] != 0) {
			//std::string type = "deepJet_comb" ;
            const auto bc_jets = cset->at("deepJet_comb")->evaluate({sys, wp, hadflav[i], abs_etas[i], cast_pts[i]});
			bweight = bc_jets;
            //std::cout << "\njet SFs up_correlated for comb at tight WP\n";
            //std::cout << "SF b/c : " << bc_jets << '\n';
        } else{ 
		//std::string type = "depJet_incl" ;
           const auto bc_jets = cset->at("deepJet_incl")->evaluate({sys, wp, hadflav[i], abs_etas[i], cast_pts[i]});
			bweight = bc_jets;
			
            //std::cout << "\njet up_correlated for comb at tight  WP\n";
            //std::cout << "SF light jets : " << bc_jets << '\n';
		}
		scalefactors_case2.emplace_back(bweight);
    }

    return scalefactors_case2;

}


float pucorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string name, std::string syst, float ntruepileup)
{
	return float(cset->at(name)->evaluate({ntruepileup, syst.c_str()}));
}

ints good_idx(ints good)
{
  ints out;
  for(unsigned int i = 0; i < good.size(); i++){
    if( good[i] )
    {
      out.emplace_back(i);
    }
  }
  return out;

}

float calculate_deltaEta( FourVector &p1, FourVector &p2)
{
  return p1.Eta() - p2.Eta();
}

float calculate_deltaPhi( FourVector &p1, FourVector &p2)
{
  return ROOT::Math::VectorUtil::DeltaPhi(p1, p2);
}

float calculate_deltaR( FourVector &p1, FourVector &p2)
{
  return ROOT::Math::VectorUtil::DeltaR(p1, p2);
}

float calculate_invMass( FourVector &p1, FourVector &p2)
{
  return ROOT::Math::VectorUtil::InvariantMass(p1, p2);
}

FourVector sum_4vec( FourVector &p1, FourVector &p2)
{
  return p1+p2;
}

FourVector sum_4Rvec( FourVector p1, FourVector p2)
{
  return p1+p2;
}

floats concatenate( double first, double second)
{
  // cout<<"A1A1A1"<<endl;
  floats out;
  out.emplace_back(first);
  out.emplace_back(second);
  // cout<<"out = "<<out<<endl;
  return out;
}

//Get indices that sort the object vectors in descending order
floats sort_discriminant( floats discr, floats obj )
{
  auto sorted_discr = Reverse(Argsort(discr));
  floats out;
  for (auto idx : sorted_discr)
  {
    out.emplace_back(obj[idx]);
  }
  return out;
}

FourVector select_leadingvec( FourVectorVec &v )
{
  FourVector vout;
  if(v.size() > 0) return v[0];
  else return vout;
}

floats PrintVector(floats myvector)
{
  for (size_t i = 0; i < myvector.size(); i++)
  {
    cout<<myvector[i]<<"\n";
  }
}

//=====================build pair example=============================================//
//W reconstruction : 2 jets
//====================================================================================//
floats w_reconstruction (FourVectorVec &jets)
{
  floats out;
  float dijetMass;
  float dijetMass_out;
  float dijetDR;
  const float Mass_W = 80.9; //W mass 
  const float Width_W = 10.8; //W width
  float X_Min = 99999;
  float X_Recwidth_W; 
  
  //Loop on all selected jets
  for(unsigned int j1 = 0; j1<jets.size()-1; j1++)
  {
    for(unsigned int j2 = j1+1; j2<jets.size(); j2++)
    {
      //select 2 jets
      dijetMass = (jets[j1]+jets[j2]).M();
      //select a best W candidate in min width in a event
      X_Recwidth_W = pow((Mass_W-dijetMass)/Width_W,2);
      if (X_Recwidth_W<X_Min)
      {
	dijetDR = ROOT::Math::VectorUtil::DeltaR(jets[j1],jets[j2]);
	X_Min = X_Recwidth_W;
	dijetMass_out = dijetMass;
      }    
    }
  }
  out.push_back(dijetMass_out);   //0:w_mass
  out.push_back(dijetDR);  
  return out;
}

floats compute_DR (FourVectorVec &muons, ints goodMuons_charge)
{
  floats out;
  float mu_ss_DR;
  float mu_os_DR;
  //cout<<"Muonsize: " << muons.size()<<endl;
  if(muons.size()>0)
    //Loop on all selected muons
    for(unsigned int mu1 = 0; mu1<muons.size()-1; mu1++)
    {
      for(unsigned int mu2 = mu1+1; mu2<muons.size(); mu2++)
      {
	//select 2 muons with same sign
	if (goodMuons_charge[mu1]!=goodMuons_charge[mu2]) continue; //check charge of muons
	mu_ss_DR = ROOT::Math::VectorUtil::DeltaR(muons[mu1],muons[mu2]);
	//select 2 muons with same sign
	if (goodMuons_charge[mu1]==goodMuons_charge[mu2]) continue;
	mu_os_DR = ROOT::Math::VectorUtil::DeltaR(muons[mu1],muons[mu2]);
	
      }
    }
  out.push_back(mu_ss_DR);		//0: same sign dimuon DR
  out.push_back(mu_os_DR);        //1: opposite sign dimuon dR
  return out;
}

floats H_reconstruction (FourVectorVec &jets, FourVectorVec &bjets)
{
  floats out;
  unsigned int hj1_idx=-1, hj2_idx=-1, wj1_idx=-1, wj2_idx=-1, tjb_idx=-1, oj_idx=-1;
  
  float H_mass, W_mass, Top_mass;
  float Chi2_H, Chi2_W, Chi2_Top;
  float Chi2_min=FLT_MAX;
  float Chi2_min_H=FLT_MAX, Chi2_min_W=FLT_MAX, Chi2_min_Top=FLT_MAX;
  // Mass and Width - 2018UL
  //const float trueSenario_mass;
  const float trueH_mass = 120.2;
  const float trueZ_mass = 90.9;
  const float trueTop_mass = 175.9;
  const float trueW_mass = 83.9;
  const float widthZ = 11.3;
  const float widthTop = 17.2;
  const float widthW = 10.8;
  
  const float widthH = 14.3;
  float H_mass_out;
        
  //Loop on all selected jets
  for(unsigned int bj1 = 0; bj1<bjets.size()-1; bj1++)
  {
    for(unsigned int bj2 = bj1+1; bj2<bjets.size(); bj2++)
    {
      //select 2 jets
      H_mass = (bjets[bj1]+bjets[bj2]).M();
      //X_Recwidth_W = pow((Mass_W-dijetMass)/Width_W,2);
      Chi2_H = pow((trueH_mass-H_mass)/widthH,2);
      if (Chi2_H < Chi2_min_H)
      {
	Chi2_min_H = Chi2_H;
	H_mass_out = H_mass;
	hj1_idx = bj1;
	hj2_idx = bj2;
      }
    }
  }
  out.push_back(H_mass_out);   //0:H_mass
  out.push_back(hj1_idx);
  out.push_back(hj2_idx);
  
  return out;
}

floats dimuon (FourVectorVec &muons, ints Muon_charge)
{
  //cout<<"AAAA"<<endl;
  floats out;
  //         float dimuonMass;
  //         float dimuonMass_out;
  float dimuonDphi;
  float dimuonDR;
  
  //Loop on all selected muons
  //	if(muons.size()>0){
  for(unsigned int j1 = 0; j1<muons.size(); j1++)
  {
    for(unsigned int j2 = j1+1; j2<muons.size(); j2++)
    {
      //dimuonMass = (muon[j1]+muon[j2]).M();
      //X_Min = X_Recwidth_W;
      //dimuonMass_out = dimuonMass;
      //select 2 muons of same charge
      //if(Muon_charge[j1]==Muon_charge[j2]){
      dimuonDphi = ROOT::Math::VectorUtil::DeltaPhi(muons[j1],muons[j2]);
      dimuonDR = ROOT::Math::VectorUtil::DeltaR(muons[j1],muons[j2]);
      //}
    }
  }
  //	}
  //         out.push_back(dijetMass_out);   //0:dilepton_mass
  out.push_back(dimuonDphi);         //1:dimuon_dphi (0 if the mass is commented)
  out.push_back(dimuonDR);         //2:dimuon_dr (1 if the mass is commented)
  
  return out;
}

bools isLfromtopHiggs (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx)//, floats &leptonpt)
{
  //cout<<"B1B1B1"<<endl;
  bools out;
  unsigned int nl = pdgId.size();
  //cout<<"nl = "<<nl<<endl;
  bool flag = false;
  bool wrong_decay = false;
  bool decay_W = false;
  bool decay_top = false;
  bool decay_Higgs = false;

  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    wrong_decay = false;
    decay_W = false;
    //cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = lgenidx[i];
    //cout<<"Nature of the particle = "<<pdgId[i]<<" "<<endl;
    //cout<<"Lepton_pt = "<<leptonpt[i]<<" "<<endl;
    //cout<<"genlidx = "<<genlidx<<" "<<endl;
    //cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<" "<<endl;
    //if(pdgId[i] == -genpdgId[genlidx])cout<<"The particle and the Gen particle don't have the same charge."<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    //cout<<"moid = "<<moid<<" nature = "<<genpdgId[moid]<<endl;
    if ( moid > genpdgId.size() || moid < -1 ) break;
    //cout<<"moidx = "<<moid<<" "<<endl;
    //cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) != 24 && (abs( genpdgId[moid] ) == 13 || abs( genpdgId[moid] ) == 11))
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      //cout<<"moidx = "<<moid<<" "<<endl;
      //cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
      if( abs( genpdgId[moid] ) == 5 )
      {
	      //cout<<"The lepton comes from a quark b."<<endl;
	      wrong_decay = true;
	      break;
      }
    }
    
    if( abs( genpdgId[moid] ) == 24 && wrong_decay == false )
    {
      decay_W = true;
      //cout<<"The lepton comes from a W."<<endl;
      int gmoid = midx[moid];
      //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
      if ( gmoid > genpdgId.size() || gmoid < -1 ) break;
      //cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;
      
      // make sure the lepton isn't coming from a quark b
      if( abs( genpdgId[gmoid] ) == 5 )
      {
	      //cout<<"The lepton comes from a quark b."<<endl;
	      wrong_decay = true;
	      break;
      }
      if( abs( genpdgId[gmoid] ) == 6 )
      {
	      //cout<<"The lepton comes from a quark t."<<endl;
	      decay_top = true;
      }
      if( abs( genpdgId[gmoid] ) == 25 )
      {
	      //cout<<"The lepton comes from a Higgs boson."<<endl;
	      decay_Higgs = true;
      }
      
      while( abs( genpdgId[gmoid] ) != 5)
      {
	      gmoid = midx[gmoid];
	      if ( gmoid > genpdgId.size() || gmoid < -1 ) break;
	      //cout<<"gmoidx = "<<gmoid<<" "<<endl;
	      //cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;
	      if( abs( genpdgId[gmoid] ) == 5 )
	      {
	        //cout<<"The lepton comes from a quark b."<<endl;
	        wrong_decay = true;
	        break;
	      }
	      if( abs( genpdgId[gmoid] ) == 6 )
	      {	
	        //cout<<"The lepton comes from a quark t."<<endl;
	        decay_top = true;
	      }
	      if( abs( genpdgId[gmoid] ) == 25 )
	      {
	        //cout<<"The lepton comes from a Higgs boson."<<endl;
	        decay_Higgs = true;
	      }	     
	      if( abs( genpdgId[gmoid] ) == 8000001 )
	      {
	        //cout<<"The lepton comes from the T'."<<endl;
	      }
      }
    }

    if( abs( genpdgId[moid] ) == 24 || abs( genpdgId[moid] ) == 15)
    {
      flag = true;
      //cout<<"yeaaaah"<<endl;
    }

    if( i == 0 && decay_W == true && wrong_decay == false && (decay_top == true || decay_Higgs == true))
    {
      flag = true;
      //cout<<"yeaaaah"<<endl;
    }
    if( i == 1 && decay_W == true && wrong_decay == false && decay_top == true && decay_Higgs == true)
    {
      flag = true;
      //cout<<"yeaaaah"<<endl;
    }
    out.emplace_back(flag);
    //cout<<"out = "<<out<<endl;
  }          
  return out;
}

FourVectorVec isLfromtoporHiggs (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx, FourVectorVec &genleptons4vecs)
{
  //cout<<"B2B2B2"<<endl;
  unsigned int nl = pdgId.size();
  //cout<<"nl = "<<nl<<endl;
  bool flag = false;
  bool wrong_decay = false;
  bool decay_W = false;
  bool decay_top = false;
  bool decay_Higgs = false;
  bool top_leading = false;
  bool Higgs_leading = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    wrong_decay = false;
    decay_W = false;
    //cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = lgenidx[i];
    //cout<<"Nature of the particle = "<<pdgId[i]<<" "<<endl;
    //cout<<"genlidx = "<<genlidx<<" "<<endl;
    //cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<" "<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    if ( moid > genpdgId.size() || moid < -1 ) break;
    //cout<<"moidx = "<<moid<<" "<<endl;
    //cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) != 24 && (abs( genpdgId[moid] ) == 13 || abs( genpdgId[moid] ) == 11))
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      //cout<<"moidx = "<<moid<<" "<<endl;
      //cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
      if( abs( genpdgId[moid] ) == 5 )
      {
	      //cout<<"The lepton comes from a quark b."<<endl;
	      wrong_decay = true;
	      break;
      }
    }
    
    if( abs( genpdgId[moid] ) == 24 && wrong_decay == false )
    {
      decay_W = true;
      //cout<<"The lepton comes from a W."<<endl;
      int gmoid = midx[moid];
      //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
      if ( gmoid > genpdgId.size() || gmoid < -1 ) break;
      //cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;
      
      // make sure the lepton isn't coming from a quark b
      if( abs( genpdgId[gmoid] ) == 5 )
      {
	      //cout<<"The lepton comes from a quark b."<<endl;
	      wrong_decay = true;
	      break;
      }
      if( abs( genpdgId[gmoid] ) == 6 )
      {
	      //cout<<"The lepton comes from a quark t."<<endl;
	      decay_top = true;
      }
      if( abs( genpdgId[gmoid] ) == 25 )
      {
	      //cout<<"The lepton comes from a Higgs boson."<<endl;
	      decay_Higgs = true;
      }
      
      while( abs( genpdgId[gmoid] ) != 5)
      {
	      gmoid = midx[gmoid];
	      if ( gmoid > genpdgId.size() || gmoid < -1 ) break;
	      //cout<<"gmoidx = "<<gmoid<<" "<<endl;
	
	      //cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;
	      if( abs( genpdgId[gmoid] ) == 5 )
	      {
	        //cout<<"The lepton comes from a quark b."<<endl;
	        wrong_decay = true;
	        break;
	      }
	      if( abs( genpdgId[gmoid] ) == 6 )
	      {
	        //cout<<"The lepton comes from a quark t."<<endl;
	        decay_top = true;
	      }
	      if( abs( genpdgId[gmoid] ) == 25 )
	      {
	        //cout<<"The lepton comes from a Higgs boson."<<endl;
	        decay_Higgs = true;
	      }	     
	      if( abs( genpdgId[gmoid] ) == 8000001 )
	      {
	        //cout<<"The lepton comes from the T'."<<endl;
	      }
      }
    }
    
    if( i == 0 && decay_W == true && wrong_decay == false && decay_top == true)
    {
      top_leading = true;
      //cout<<"yeaaaah"<<endl;
    }
    if( i == 1 && decay_W == true && wrong_decay == false && decay_top == true && decay_Higgs == true && top_leading == true)
    {
      //cout<<"yeaaaah"<<endl;
    }
    if( i == 0 && decay_W == true && wrong_decay == false && decay_Higgs == true)
    {
      Higgs_leading = true;
      //cout<<"yeaaaah"<<endl;
    }
    if( i == 1 && decay_W == true && wrong_decay == false && decay_top == true && decay_Higgs == true && Higgs_leading == true)
    {
      swap(genleptons4vecs[0],genleptons4vecs[1]);
      //cout<<"yeaaaah"<<endl;
      //cout<<"Swap!"<<endl;
    }
  }          
  return genleptons4vecs;
}

bools isLfromW (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx)
{
  // cout<<"B3B3B3"<<endl;
  bools out;
  unsigned int nl = pdgId.size();
  // cout<<"nl = "<<nl<<endl;
  bool flag = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    // cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = lgenidx[i];
    // cout<<"Nature of the particle = "<<pdgId[i]<<" "<<endl;
    // cout<<"genlidx = "<<genlidx<<" "<<endl;
    // cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<" "<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    //cout<<"moid = "<<moid<<" nature = "<<genpdgId[moid]<<endl;
    //if ( moid > genpdgId.size() || moid < -1 ) break;
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) != 24 && (abs( genpdgId[moid] ) == 13 || abs( genpdgId[moid] ) == 11))
    {
      moid = midx[moid];
      //if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    }
    
    if( abs( genpdgId[moid] ) == 24)
    {
      flag = true;
      // cout<<"yeaaaah"<<endl;
    }

//     // check if the lepton comes from a tau or a Z boson
//     while( abs( genpdgId[moid] ) != 15) // && abs( genpdgId[moid] ) != 23 )
//     {
//       moid = midx[moid];
//       if ( moid > genpdgId.size() || moid < -1 ) break;
//       //cout<<"moidx = "<<moid<<" "<<endl;
//       //cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
//     }
//     //cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;

//     if( abs( genpdgId[moid] ) == 15 || abs( genpdgId[moid] ) == 23 )
//     {
//       flag = true;
//       //cout<<"The lepton comes from a tau or a Z boson."<<endl;
//     }

//     if( abs( genpdgId[moid] ) == 24)
//     {
//       while( abs( genpdgId[moid] ) != 8000001)
//       {
// 	   moid = midx[moid];
// 	   if ( moid > genpdgId.size() || moid < -1 ) break;
// 	   cout<<"moidx = "<<moid<<" "<<endl;
// 	   cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;

// 	   if( abs( genpdgId[moid] ) == 15 || abs( genpdgId[moid] ) == 23 )
// 	   {
// 	     flag = false;
// 	     cout<<"The lepton comes from a tau or a Z boson."<<endl;
// 	     break;
// 	   }
//       }

//       if( abs( genpdgId[moid] ) == 8000001)
//       {
// 	   flag = true;
// 	   cout<<"yeaaaah"<<endl;
//       }
//     }

    out.emplace_back(flag);
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

ints MotherfromL (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx)
{
  //cout<<"B4B4B4"<<endl;
  ints out = {-1, -1};
  unsigned int nl = pdgId.size();
  // cout<<"nl = "<<nl<<endl;
  bool flag = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    // cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = lgenidx[i];
    // cout<<"Nature of the particle = "<<pdgId[i]<<" "<<endl;
    // cout<<"genlidx = "<<genlidx<<" "<<endl;
    // cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<" "<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    // cout<<"moid = "<<moid<<" nature = "<<genpdgId[moid]<<endl;
    if ( moid > genpdgId.size() || moid < -1 ) break;
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while(abs( genpdgId[moid] ) == 13 || abs( genpdgId[moid] ) == 11)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    }
    
    // cout<<"Nature of its mmother = "<<genpdgId[moid]<<" "<<endl;
    out[i] = genpdgId[moid];
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isWfromHiggs (ints &genpdgId, ints &wgenidx, ints &midx, floats &genmass)
{
  // cout<<"B5B5B5"<<endl;
  bools out;
  unsigned int nl = wgenidx.size();
  // cout<<"nl = "<<nl<<endl;
  bool flag = false;
  // bool onshell = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    // onshell = false;
    // cout<<"i = "<<i<<endl;
    
    int genlidx = wgenidx[i];
    // cout<<"genlidx = "<<genlidx<<" "<<endl;
    // cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<", Mass of the Gen particle = "<<genmass[genlidx]<<" "<<endl;
    // if(genmass[genlidx] > 70 && genmass[genlidx] < 90)
    // {
    //   onshell = true;
    // }
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    // cout<<"moid = "<<moid<<" nature = "<<genpdgId[moid]<<endl;
    if ( moid > genpdgId.size() || moid < -1 ) break;
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<", Mass of the mother = "<<genmass[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) == 24)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<", Mass of the mother = "<<genmass[moid]<<" "<<endl;
    }
    
    if( abs( genpdgId[moid] ) == 25)// && onshell == true)
    {
      flag = true;
      // cout<<"yeaaaah"<<endl;
    }

    out.emplace_back(flag);
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isWfromtop (ints &genpdgId, ints &wgenidx, ints &midx, floats &genmass)
{
  // cout<<"B6B6B6"<<endl;
  bools out;
  unsigned int nl = wgenidx.size();
  // cout<<"nl = "<<nl<<endl;
  bool flag = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    //cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = wgenidx[i];
    // cout<<"genlidx = "<<genlidx<<" "<<endl;
    // cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<", Mass of the Gen particle = "<<genmass[genlidx]<<" "<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    // cout<<"moid = "<<moid<<" nature = "<<genpdgId[moid]<<endl;
    if ( moid > genpdgId.size() || moid < -1 ) break;
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<", Mass of the mother = "<<genmass[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) == 24)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<", Mass of the mother = "<<genmass[moid]<<" "<<endl;
    }
    
    if( abs( genpdgId[moid] ) == 6)
    {
      flag = true;
      // cout<<"yeaaaah"<<endl;
    }

    out.emplace_back(flag);
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isLfromtopantitop (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx)//, floats &leptonpt)
{
  // cout<<"B7B7B7"<<endl;
  bools out;
  unsigned int nl = pdgId.size();
  // cout<<"nl = "<<nl<<endl;
  bool flag = false;
  bool wrong_decay = false;
  bool decay_W = false;
  bool decay_top = false;
  bool decay_antitop = false;

  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    wrong_decay = false;
    decay_W = false;
    // cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = lgenidx[i];
    // cout<<"Nature of the particle = "<<pdgId[i]<<" "<<endl;
    //cout<<"Lepton_pt = "<<leptonpt[i]<<" "<<endl;
    // cout<<"genlidx = "<<genlidx<<" "<<endl;
    // cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<" "<<endl;
    //if(pdgId[i] == -genpdgId[genlidx])cout<<"The particle and the Gen particle don't have the same charge."<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    // cout<<"moid = "<<moid<<" nature = "<<genpdgId[moid]<<endl;
    // if ( moid > genpdgId.size() || moid < -1 ) break;
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) != 24 && (abs( genpdgId[moid] ) == 13 || abs( genpdgId[moid] ) == 11))
    {
      moid = midx[moid];
      // if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
      if( abs( genpdgId[moid] ) == 5 )
      {
	      // cout<<"The lepton comes from a quark b."<<endl;
	      wrong_decay = true;
	      break;
      }
    }
    
    if( abs( genpdgId[moid] ) == 24 && wrong_decay == false )
    {
      decay_W = true;
      // cout<<"The lepton comes from a W."<<endl;
      int gmoid = midx[moid];
      //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
      if ( gmoid > genpdgId.size() || gmoid < -1 ) break;
      // cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;
      
      // make sure the lepton isn't coming from a quark b
      if( abs( genpdgId[gmoid] ) == 5 )
      {
	      // cout<<"The lepton comes from a quark b."<<endl;
	      wrong_decay = true;
	      break;
      }
      if( genpdgId[gmoid] == 6 )
      {
	      // cout<<"The lepton comes from a quark t."<<endl;
	      decay_top = true;
      }
      if( genpdgId[gmoid] == -6 )
      {
	      // cout<<"The lepton comes from an antiquark t."<<endl;
	      decay_antitop = true;
      }
      
      while( abs( genpdgId[gmoid] ) != 5)
      {
	      gmoid = midx[gmoid];
	      if ( gmoid > genpdgId.size() || gmoid < -1 ) break;
	      // cout<<"gmoidx = "<<gmoid<<" "<<endl;
	      // cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;
	      if( abs( genpdgId[gmoid] ) == 5 )
	      {
	        //cout<<"The lepton comes from a quark b."<<endl;
	        wrong_decay = true;
	        break;
	      }
	      if( genpdgId[gmoid] == 6 )
	      {	
	        // cout<<"The lepton comes from a quark t."<<endl;
	        decay_top = true;
	      }
	      if( genpdgId[gmoid] == -6 )
	      {
	        // cout<<"The lepton comes from an antiquark t."<<endl;
	        decay_antitop = true;
	      }	     
	      if( abs( genpdgId[gmoid] ) == 8000001 )
	      {
	        // cout<<"The lepton comes from the T'."<<endl;
	      }
      }
    }

    if( i == 0 && decay_W == true && wrong_decay == false && (decay_top == true || decay_antitop == true))
    {
      flag = true;
      // cout<<"yeaaaah"<<endl;
    }
    if( i == 1 && decay_W == true && wrong_decay == false && decay_top == true && decay_antitop == true)
    {
      flag = true;
      // cout<<"yeaaaah"<<endl;
    }
    out.emplace_back(flag);
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

FourVectorVec isLfromtoporantitop (ints &pdgId, ints &genpdgId, ints &lgenidx, ints &midx, FourVectorVec &genleptons4vecs)
{
  cout<<"B8B8B8"<<endl;
  unsigned int nl = pdgId.size();
  cout<<"nl = "<<nl<<endl;
  bool flag = false;
  bool wrong_decay = false;
  bool decay_W = false;
  bool decay_top = false;
  bool decay_antitop = false;
  bool top_leading = false;
  bool antitop_leading = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    wrong_decay = false;
    decay_W = false;
    cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = lgenidx[i];
    cout<<"Nature of the particle = "<<pdgId[i]<<" "<<endl;
    cout<<"genlidx = "<<genlidx<<" "<<endl;
    cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<" "<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    // if ( moid > genpdgId.size() || moid < -1 ) break;
    cout<<"moidx = "<<moid<<" "<<endl;
    cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) != 24 && (abs( genpdgId[moid] ) == 13 || abs( genpdgId[moid] ) == 11))
    {
      moid = midx[moid];
      // if ( moid > genpdgId.size() || moid < -1 ) break;
      cout<<"moidx = "<<moid<<" "<<endl;
      cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
      if( abs( genpdgId[moid] ) == 5 )
      {
	      //cout<<"The lepton comes from a quark b."<<endl;
	      wrong_decay = true;
	      break;
      }
    }
    
    if( abs( genpdgId[moid] ) == 24 && wrong_decay == false )
    {
      decay_W = true;
      cout<<"The lepton comes from a W."<<endl;
      int gmoid = midx[moid];
      //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
      if ( gmoid > genpdgId.size() || gmoid < -1 ) break;
      cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;
      
      // make sure the lepton isn't coming from a quark b
      if( abs( genpdgId[gmoid] ) == 5 )
      {
	      //cout<<"The lepton comes from a quark b."<<endl;
	      wrong_decay = true;
	      break;
      }
      if( genpdgId[gmoid] == 6 )
      {
	      cout<<"The lepton comes from a quark t."<<endl;
	      decay_top = true;
      }
      if( genpdgId[gmoid] == -6 )
      {
	      cout<<"The lepton comes from an antiquark t."<<endl;
	      decay_antitop = true;
      }
      
      while( abs( genpdgId[gmoid] ) != 5)
      {
	      gmoid = midx[gmoid];
	      if ( gmoid > genpdgId.size() || gmoid < -1 ) break;
	      cout<<"gmoidx = "<<gmoid<<" "<<endl;
	
	      cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;
	      if( abs( genpdgId[gmoid] ) == 5 )
	      {
	        //cout<<"The lepton comes from a quark b."<<endl;
	        wrong_decay = true;
	        break;
	      }
	      if( genpdgId[gmoid] == 6 )
	      {
	        cout<<"The lepton comes from a quark t."<<endl;
	        decay_top = true;
	      }
	      if( genpdgId[gmoid] == -6 )
	      {
	        cout<<"The lepton comes from an antiquark t."<<endl;
	        decay_antitop = true;
	      }	     
	      if( abs( genpdgId[gmoid] ) == 8000001 )
	      {
	        //cout<<"The lepton comes from the T'."<<endl;
	      }
      }
    }
    
    if( i == 0 && decay_W == true && wrong_decay == false && decay_top == true)
    {
      top_leading = true;
      cout<<"yeaaaah"<<endl;
    }
    if( i == 1 && decay_W == true && wrong_decay == false && decay_top == true && decay_antitop == true && top_leading == true)
    {
      cout<<"yeaaaah"<<endl;
    }
    if( i == 0 && decay_W == true && wrong_decay == false && decay_antitop == true)
    {
      antitop_leading = true;
      cout<<"yeaaaah"<<endl;
    }
    if( i == 1 && decay_W == true && wrong_decay == false && decay_top == true && decay_antitop == true && antitop_leading == true)
    {
      swap(genleptons4vecs[0],genleptons4vecs[1]);
      cout<<"yeaaaah"<<endl;
      cout<<"Swap!"<<endl;
    }
  }          
  return genleptons4vecs;
}

bools isWfromtoponly (ints &genpdgId, ints &wgenidx, ints &midx, floats &genmass)
{
  cout<<"B9B9B9"<<endl;
  bools out;
  unsigned int nl = wgenidx.size();
  cout<<"nl = "<<nl<<endl;
  bool flag = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = wgenidx[i];
    cout<<"genlidx = "<<genlidx<<" "<<endl;
    cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<", Mass of the Gen particle = "<<genmass[genlidx]<<" "<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    cout<<"moid = "<<moid<<" nature = "<<genpdgId[moid]<<endl;
    if ( moid > genpdgId.size() || moid < -1 ) break;
    cout<<"moidx = "<<moid<<" "<<endl;
    cout<<"Nature of its mother = "<<genpdgId[moid]<<", Mass of the mother = "<<genmass[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) == 24)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      cout<<"moidx = "<<moid<<" "<<endl;
      cout<<"Nature of its mother = "<<genpdgId[moid]<<", Mass of the mother = "<<genmass[moid]<<" "<<endl;
    }
    
    if( genpdgId[moid] == 6)
    {
      flag = true;
      cout<<"yeaaaah"<<endl;
    }

    out.emplace_back(flag);
    cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isWfromantitop (ints &genpdgId, ints &wgenidx, ints &midx, floats &genmass)
{
  cout<<"B10B10B10"<<endl;
  bools out;
  unsigned int nl = wgenidx.size();
  cout<<"nl = "<<nl<<endl;
  bool flag = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < nl; i++)
  {
    flag = false;
    cout<<"i = "<<i<<endl;
    
    // matching the lepton to the same lepton in the Generator tree
    int genlidx = wgenidx[i];
    cout<<"genlidx = "<<genlidx<<" "<<endl;
    cout<<"Nature of the Gen particle = "<<genpdgId[genlidx]<<", Mass of the Gen particle = "<<genmass[genlidx]<<" "<<endl;
    
    // tracking the mother
    int moid = midx[genlidx];
    //just to be sure if it is not out of the range... (otherwise it arises segmentation fault)
    cout<<"moid = "<<moid<<" nature = "<<genpdgId[moid]<<endl;
    if ( moid > genpdgId.size() || moid < -1 ) break;
    cout<<"moidx = "<<moid<<" "<<endl;
    cout<<"Nature of its mother = "<<genpdgId[moid]<<", Mass of the mother = "<<genmass[moid]<<" "<<endl;
    
    // check if the lepton comes from a W
    while( abs( genpdgId[moid] ) == 24)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      cout<<"moidx = "<<moid<<" "<<endl;
      cout<<"Nature of its mother = "<<genpdgId[moid]<<", Mass of the mother = "<<genmass[moid]<<" "<<endl;
    }
    
    if( genpdgId[moid] == -6)
    {
      flag = true;
      cout<<"yeaaaah"<<endl;
    }

    out.emplace_back(flag);
    cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools CheckOverlaps(FourVectorVec &jets, FourVectorVec &leps)
{       
  bools out;
  doubles mindrlepton;
  bool flag = false;
  for (auto ajet: jets)
  {
    flag = false;
    auto mindr = 6.0;
    for (auto alepton: leps)
    {
      auto dr = ROOT::Math::VectorUtil::DeltaR(ajet, alepton);
      if(dr < mindr)
      { 
        mindr = dr;
      }
    }
    if(mindr > 0.4)
    {
      flag = true;
    }
    else
    {
      flag = false;
    }
    out.emplace_back(flag);
  }
  return out;
}

RVec< RVec<int> > distinct_comb_3(int n)
{
  // cout<<"C1C1C1"<<endl;
  // cout <<"n = "<<n<< endl;
  RVec<int> idx1;
  RVec<int> idx2;
  RVec<int> idx3;
  RVec< RVec<int> > res;
  idx1.reserve(n*(n-1)/2);
  idx2.reserve(n*(n-1)/2);
  idx3.reserve(n*(n-1)/2);

  if(n < 3)
  {
    idx1.emplace_back(0);
    idx2.emplace_back(0);
    idx3.emplace_back(0);
  }  

  for(int i = 0; i < n-2; i++)
  {
    for(int j = i+1; j < n-1; j++)
    {
      for(unsigned int k = j+1; k < n; k++)
      {
	      idx1.emplace_back(i);
	      idx2.emplace_back(j);
	      idx3.emplace_back(k);
      }
    }
  }      
  // cout<<"idx1 = "<<idx1<<", idx2 = "<<idx2<<"idx3 = "<<idx3<<endl;
  res.emplace_back(idx1);
  res.emplace_back(idx2);
  res.emplace_back(idx3);
  return res;      
};

FourVectorRVec Threejets (FourVectorRVec jet1, FourVectorRVec jet2, FourVectorRVec jet3, int jet_number, int bjet_number)
{
  // cout<<"C2C2C2"<<endl;  
  FourVectorRVec threejets;
  if(jet_number>=3 && bjet_number>=1)
  {
    threejets = jet1 + jet2 + jet3;
    //cout<<"njet = "<<jet_number<<endl; 
  }
  else
  {
    threejets = jet1 - jet1;
    //cout<<"njet = "<<jet_number<<endl; 
  }
  // cout<<"threejets = "<<threejets<<endl;
  return threejets;
}

FourVectorRVec Twojets (FourVectorRVec jet1, FourVectorRVec jet2, FourVectorRVec jet3, floats &jet1_btag, floats &jet2_btag, floats &jet3_btag, int jet_number, int bjet_number)
{
  // cout<<"C3C3C3"<<endl;  
  FourVectorRVec twojets;
  int ncombi = jet1_btag.size();
  // cout<<"jet1_btag = "<<jet1_btag<<", jet2_btag = "<<jet2_btag<<", jet3_btag = "<<jet3_btag<<endl;
  if(jet_number>=3 && bjet_number>=1)
  {
    // cout<<"pouet"<<endl;
    for(unsigned int i = 0; i < ncombi; i++)
    {
      if (jet3_btag[i] < jet1_btag[i] && jet3_btag[i] < jet2_btag[i])
      {
        // cout<<"jet 1 and jet 2 are chosen for W."<<endl;
        twojets = jet1 + jet2;
      }
      if (jet2_btag[i] < jet1_btag[i] && jet2_btag[i] < jet3_btag[i])
      {
        // cout<<"jet 1 and jet 3 are chosen for W."<<endl;
        twojets = jet1 + jet3;
      }
      if (jet1_btag[i] < jet2_btag[i] && jet1_btag[i] < jet3_btag[i])
      {
        // cout<<"jet 2 and jet 3 are chosen for W."<<endl;
        twojets = jet2 + jet3;
      }
      //cout<<"njet = "<<jet_number<<endl; 
    }
  }
  else
  {
    twojets = jet1 - jet1;
    //cout<<"njet = "<<jet_number<<endl; 
  }
  return twojets;
}

bools isPfromtopHiggs (ints &genpdgId, ints &midx, ints &lgenidx)
{
  // cout<<"D1D1D1"<<endl;
  bools out;
  unsigned int np = genpdgId.size();
  // cout<<"\n np = "<<np<<endl;
  bool flag = false;
  bool lepton = false;
  bool decay_top = false;
  bool decay_Higgs = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < np; i++)
  {
    flag = false;
    decay_top = false;
    decay_Higgs = false;
    lepton = false;
    // cout<<"i = "<<i<<endl;
    // cout<<"Nature of the particle = "<<genpdgId[i]<<" "<<endl;
    if(lgenidx.size()==2)
    {
      if(i == lgenidx[0] || i == lgenidx[1])
      {
        // cout<<"The particle is a selected electron or a muon."<<endl;
        lepton = true;
      }
    }
    
    // tracking the mother
    int moid = midx[i];
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the particle comes from a top or a Higgs
    while(abs(genpdgId[moid]) != 6 && abs(genpdgId[moid]) != 25 && lepton == false)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    }
    
    if(abs(genpdgId[moid]) == 6)
    {
      // cout<<"The particle comes from a quark t."<<endl;
      decay_top = true;
    }
    
    if(abs(genpdgId[moid]) == 25)
    {
      // cout<<"The particle comes from a Higgs boson."<<endl;
      decay_Higgs = true;
    }
    
    if((decay_top == true || decay_Higgs == true) && lepton == false)
    {
      // int gmoid = midx[gmoid];
      //cout<<"gmoidx = "<<gmoid<<" "<<endl;
      //cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;	     
      // if(abs(genpdgId[gmoid]) == 8000001)
      // {
	      //cout<<"The particle comes from the T'."<<endl;
      // }
      flag = true;
      //cout<<"yeaaaah"<<endl;
    }
    out.emplace_back(flag);
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isPfromHiggs (ints &genpdgId, ints &midx, ints &lgenidx)
{
  // cout<<"D2D2D2"<<endl;
  bools out;
  unsigned int np = genpdgId.size();
  // cout<<"np = "<<np<<endl;
  bool flag = false;
  bool lepton = false;
  bool decay_Higgs = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < np; i++)
  {
    flag = false;
    decay_Higgs = false;
    lepton = false;
    // cout<<"i = "<<i<<endl;
    // cout<<"Nature of the particle = "<<genpdgId[i]<<" "<<endl;
    if(lgenidx.size()==2)
    {
      if(i == lgenidx[0] || i == lgenidx[1])
      {
        // cout<<"The particle is a selected electron or a muon."<<endl;
        lepton = true;
      }
    }
    
    // tracking the mother
    int moid = midx[i];
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the particle comes from a top or a Higgs
    while(abs(genpdgId[moid]) != 25 && lepton == false)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    }
      
    if(abs(genpdgId[moid]) == 25)
    {
      // cout<<"The particle comes from a Higgs boson."<<endl;
      decay_Higgs = true;
    }
    
    if(decay_Higgs == true && lepton == false)
    { 
      int gmoid = midx[moid];
      // cout<<"gmoidx = "<<gmoid<<" "<<endl;
      // cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;	     
      if(abs(genpdgId[gmoid]) == 8000001)
      {
      	// cout<<"The particle comes from the T'."<<endl;
      }
      flag = true;
      // cout<<"yeaaaah"<<endl;
    }
    out.emplace_back(flag);
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isPfromtop (ints &genpdgId, ints &midx, ints &lgenidx)
{
  // cout<<"D3D3D3"<<endl;
  bools out;
  unsigned int np = genpdgId.size();
  // cout<<"np = "<<np<<endl;
  bool flag = false;
  bool lepton = false;
  bool decay_top = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < np; i++)
  {
    flag = false;
    decay_top = false;
    lepton = false;
    // cout<<"i = "<<i<<endl;
    // cout<<"Nature of the particle = "<<genpdgId[i]<<" "<<endl;
    if(lgenidx.size()==2)
    {
      if(i == lgenidx[0] || i == lgenidx[1])
      {
        // cout<<"The particle is a selected electron or a muon."<<endl;
        lepton = true;
      }
    }

    // tracking the mother
    int moid = midx[i];
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the particle comes from a top or a Higgs
    while(abs(genpdgId[moid]) != 6 && lepton == false)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    }
    
    if(abs(genpdgId[moid]) == 6)
    {
      //cout<<"The particle comes from a quark t."<<endl;
      decay_top = true;
    }
    
    if(decay_top == true && lepton == false)
    {
      // int gmoid = midx[gmoid];
      // cout<<"gmoidx = "<<gmoid<<" "<<endl;
      // cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;	     
      // if(abs(genpdgId[gmoid]) == 8000001)
      // {
	      //cout<<"The particle comes from the T'."<<endl;
      // }
      flag = true;
      //cout<<"yeaaaah"<<endl;
    }
    out.emplace_back(flag);
    //cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isPfromtoponly (ints &genpdgId, ints &midx, ints &lgenidx)
{
  // cout<<"D4D4D4"<<endl;
  bools out;
  unsigned int np = genpdgId.size();
  // cout<<"np = "<<np<<endl;
  bool flag = false;
  bool lepton = false;
  bool decay_top = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < np; i++)
  {
    flag = false;
    decay_top = false;
    lepton = false;
    // cout<<"i = "<<i<<endl;
    // cout<<"Nature of the particle = "<<genpdgId[i]<<" "<<endl;
    if(lgenidx.size()==2)
    {
      if(i == lgenidx[0] || i == lgenidx[1])
      {
        // cout<<"The particle is a selected electron or a muon."<<endl;
        lepton = true;
      }
    }

    // tracking the mother
    int moid = midx[i];
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the particle comes from a top or a Higgs
    while(genpdgId[moid] != 6 && lepton == false)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    }
    
    if(genpdgId[moid] == 6)
    {
      // cout<<"The particle comes from a quark t."<<endl;
      decay_top = true;
    }
    
    if(decay_top == true && lepton == false)
    {
      // int gmoid = midx[gmoid];
      // cout<<"gmoidx = "<<gmoid<<" "<<endl;
      // cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;	     
      // if(abs(genpdgId[gmoid]) == 8000001)
      // {
	      //cout<<"The particle comes from the T'."<<endl;
      // }
      flag = true;
      // cout<<"yeaaaah"<<endl;
    }
    out.emplace_back(flag);
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isPfromantitop (ints &genpdgId, ints &midx, ints &lgenidx)
{
  // cout<<"D5D5D5"<<endl;
  bools out;
  unsigned int np = genpdgId.size();
  // cout<<"np = "<<np<<endl;
  bool flag = false;
  bool lepton = false;
  bool decay_antitop = false;
  
  // GenParticle loop 
  for(unsigned int i = 0; i < np; i++)
  {
    flag = false;
    decay_antitop = false;
    lepton = false;
    // cout<<"i = "<<i<<endl;
    // cout<<"Nature of the particle = "<<genpdgId[i]<<" "<<endl;
    if(lgenidx.size()==2)
    {
      if(i == lgenidx[0] || i == lgenidx[1])
      {
        // cout<<"The particle is a selected electron or a muon."<<endl;
        lepton = true;
      }
    }

    // tracking the mother
    int moid = midx[i];
    // cout<<"moidx = "<<moid<<" "<<endl;
    // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    
    // check if the particle comes from a top or a Higgs
    while(genpdgId[moid] != -6 && lepton == false)
    {
      moid = midx[moid];
      if ( moid > genpdgId.size() || moid < -1 ) break;
      // cout<<"moidx = "<<moid<<" "<<endl;
      // cout<<"Nature of its mother = "<<genpdgId[moid]<<" "<<endl;
    }
    
    if(genpdgId[moid] == -6)
    {
      // cout<<"The particle comes from an antiquark t."<<endl;
      decay_antitop = true;
    }
    
    if(decay_antitop == true && lepton == false)
    {
      // int gmoid = midx[gmoid];
      // cout<<"gmoidx = "<<gmoid<<" "<<endl;
      // cout<<"Nature of its grandmother = "<<genpdgId[gmoid]<<" "<<endl;	     
      // if(abs(genpdgId[gmoid]) == 8000001)
      // {
	      //cout<<"The particle comes from the T'."<<endl;
      // }
      flag = true;
      // cout<<"yeaaaah"<<endl;
    }
    out.emplace_back(flag);
    // cout<<"out = "<<out<<endl;
  }          
  return out;
}

bools isdR04 (floats &eta1, floats &eta2, floats &phi1, floats &phi2)
{
  // cout<<"EEE"<<endl;
  bools out;
  for( unsigned int i = 0; i < eta1.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    bool isless04 = false;	 
    float mindr = FLT_MAX;
    // cout<<"isless04 = "<<isless04<<", mindr = "<<mindr<<endl;
    for( unsigned int j = 0; j < eta2.size(); j++ )
    {
      // cout<<"j = "<<j<<endl;
      auto dr = ROOT::VecOps::DeltaR(eta1[i],eta2[j],phi1[i],phi2[j]);
      if( dr < mindr ) mindr = dr;
      // cout<<"mindr = "<<mindr<<endl;
    }
    if( mindr < 0.4 ) isless04 = true;
    // if(isless04 == true) cout<<" mindr is lower than 0.4."<<endl;	 
    out.emplace_back(isless04);
    // cout<<"out = "<<out<<endl;
  }
  return out;
}

bools isJetFromgJet (floats &jet_pt, floats &gjet_pt, floats &gjetTp_pt, ints &jet_gidx)//, floats &jet_eta, floats &jet_phi, floats &jet_mass, ints &jet_jetid, floats &gjet_eta, floats &gjet_phi, floats &gjet_mass)
{
  // cout<<"F1F1F1"<<endl;
  bools out;
  for( unsigned int i = 0; i < jet_pt.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    bool isFromTprime = false;
    //index of genJet matched with Jet i
    int gIdx = jet_gidx[i];
    // cout<<"jet_pt[i] = "<<jet_pt[i]<<", gIdx = "<<gIdx<<", gjet_pt[gIdx] = "<<gjet_pt[gIdx]<<endl;
//     cout<<"jet_eta[i] = "<<jet_eta[i]<<", gjet_eta[gIdx] = "<<gjet_eta[gIdx]<<endl;	   
//     cout<<"jet_phi[i] = "<<jet_phi[i]<<", gjet_phi[gIdx] = "<<gjet_phi[gIdx]<<endl;
//     cout<<"jet_mass[i] = "<<jet_mass[i]<<", gjet_mass[gIdx] = "<<gjet_mass[gIdx]<<endl;
//     //cout<<"jet_energy[i] = "<<jet_energy[i]<<", gjet_energy[gIdx] = "<<gjet_energy[gIdx]<<endl;
//     cout<<"jet_jetId[i] = "<<jet_jetid[i]<<endl;
    if ( gIdx <= int(gjet_pt.size()) && gIdx > -1 )
    {
      // cout<<"F11F11F11"<<endl;
      for( unsigned int j = 0; j < gjetTp_pt.size(); j++ )
      {
	      // cout<<"j = "<<j<<endl;
	      // cout<<"gjetTp_pt[j] = "<<gjetTp_pt[j]<<endl;
	      if( abs((gjet_pt[gIdx] - gjetTp_pt[j])) < 0.01 )
	      {
	        isFromTprime = true;
	        // if(gjet_pt[gIdx] != gjetTp_pt[j]) cout<<"gjet_pt[gIdx] - gjetTp_pt[j] = "<<gjet_pt[gIdx] - gjetTp_pt[j]<<endl;
	        // cout<<"yeaaaah"<<endl;
	      }
      }
    }
    out.emplace_back(isFromTprime);
    // cout<<"out = "<<out<<endl;
  }
  return out;
}

bools isJetFromgJet_2 (double &jet_pt, floats &gjet_pt, floats &gjetTp_pt, int &jet_gidx)
{
  //cout<<"F2F2F2"<<endl;
  bools out;
  bool isFromTprime = false;
  int gIdx = jet_gidx;
  if ( gIdx <= int(gjet_pt.size()) && gIdx > -1 )
  {
    for( unsigned int j = 0; j < gjetTp_pt.size(); j++ )
    {
      //cout<<"j = "<<j<<endl;
      //cout<<"jet_pt = "<<jet_pt<<", gjet_pt[gIdx] = "<<gjet_pt[gIdx]<<", gjetTp_pt[j] = "<<gjetTp_pt[j]<<endl;
      if( abs((gjet_pt[gIdx] - gjetTp_pt[j])) < 0.01 )
      {
	isFromTprime = true;
	//cout<<"yeaaaah"<<endl;
      }
    }
  }
  out.emplace_back(isFromTprime);
  //cout<<"out = "<<out<<endl;
  return out;
}

floats DR_1 (floats &eta1, floats &eta2, floats &phi1, floats &phi2)
{
  // cout<<"G1G1G1"<<endl;
  floats out;
  float deltaR = FLT_MAX;
  for( unsigned int i = 0; i < eta1.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    for( unsigned int j = 0; j < eta2.size(); j++ )
    {
      // cout<<"j = "<<j<<endl;
      deltaR = FLT_MAX;
      auto deltaR = ROOT::VecOps::DeltaR(eta1[i],eta2[j],phi1[i],phi2[j]);
	    // cout<<"deltaR = "<<deltaR<<endl;
      out.emplace_back(deltaR);
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

floats DR_2 (double &eta1, floats &eta2, double &phi1, floats &phi2)
{
  // cout<<"G2G2G2"<<endl;
  floats out;
  float deltaR = FLT_MAX;
  for( unsigned int i = 0; i < eta2.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    deltaR = FLT_MAX;
    auto deltaR = ROOT::VecOps::DeltaR(float(eta1),eta2[i],float(phi1),phi2[i]);
    // cout<<"deltaR = "<<deltaR<<endl;
    out.emplace_back(deltaR);
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

int maxDR (double &eta1, floats &eta2, double &phi1, floats &phi2)
{
  //cout<<"G3G3G3"<<endl;
  int out = -1;
  float maxdr = -FLT_MAX;
  for( unsigned int i = 0; i < eta2.size(); i++ )
  {
    //cout<<"i = "<<i<<endl; 
    auto dr = ROOT::VecOps::DeltaR(float(eta1),eta2[i],float(phi1),phi2[i]);
    if( dr > maxdr )
    {
      out = i;
      maxdr = dr;
      //cout<<"maxdr = "<<maxdr<<endl;
    }
  }
  //cout<<"out = "<<out<<endl;
  return out;
}

ints minDR_1 (floats &eta1, floats &eta2, floats &phi1, floats &phi2)
{
  // cout<<"G4G4G4"<<endl;
  ints out = {-1, -1}; 
  ints mout = {-1, -1};
  float mindr = FLT_MAX;
  float mmindr = FLT_MAX;
  // cout<<"nj = "<<eta2.size()<<endl;
  if(eta2.size() < 3)
  {
    for( unsigned int i = 0; i < eta1.size(); i++ )
    {
      // cout<<"i = "<<i<<endl;
      for( unsigned int j = 0; j < eta2.size(); j++ )
      {
        // cout<<"j = "<<j<<endl;
        auto dr = ROOT::VecOps::DeltaR(eta1[i],eta2[j],phi1[i],phi2[j]);
        if( dr < mindr )
        {
	        out = {i, j};
	        mindr = dr;
	        // cout<<"mindr = "<<mindr<<endl;
        }
      }
    }
  }
  else if(eta2.size() >= 3)
  {
    for( unsigned int i = 0; i < eta1.size(); i++ )
    {
      // cout<<"i = "<<i<<endl;
      for( unsigned int j = 0; j < eta2.size(); j++ )
      {
        // cout<<"j = "<<j<<endl;
        auto dr = ROOT::VecOps::DeltaR(eta1[i],eta2[j],phi1[i],phi2[j]);
        if( dr < mindr )
        {
          if( dr < mmindr )
          {
            out = mout;
            mout = {i, j};
            mindr = mmindr;
            mmindr = dr;
            // cout<<"mindr = "<<mindr<<endl;
            // cout<<"mmindr = "<<mmindr<<endl;
          }
          else
          {
            out = {i, j};
            mindr = dr;
            // cout<<"mindr = "<<mindr<<endl;
            // cout<<"mmindr = "<<mmindr<<endl;
          }
        }
      }
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

ints minDR_2 (floats &eta1, floats &eta2, floats &phi1, floats &phi2)
{
  // cout<<"G5G5G5"<<endl;
  ints out = {-1, -1}; 
  float mindr = FLT_MAX;
  for( unsigned int i = 0; i < eta1.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;	 
    for( unsigned int j = i+1; j < eta2.size(); j++ )
    {
      // cout<<"j = "<<j<<endl;
      auto dr = ROOT::VecOps::DeltaR(eta1[i],eta2[j],phi1[i],phi2[j]);
      if( dr < mindr )
      {
	      out = {i, j};
	      mindr = dr;
	      // cout<<"mindr = "<<mindr<<endl;
      }
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

int minDR_3 (double &eta1, floats &eta2, double &phi1, floats &phi2)
{
  //cout<<"G6G6G6"<<endl;
  int out = -1;
  float mindr = FLT_MAX;
  for( unsigned int i = 0; i < eta2.size(); i++ )
  {
    //cout<<"i = "<<i<<endl; 
    auto dr = ROOT::VecOps::DeltaR(float(eta1),eta2[i],float(phi1),phi2[i]);
    if( dr < mindr )
    {
      out = i;
      mindr = dr;
      //cout<<"mindr = "<<mindr<<endl;
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

ints minDR_4 (floats &eta1, floats &eta2, floats &phi1, floats &phi2, floats &jet_pt, double &bjet_pt)
{
  //cout<<"G7G7G7"<<endl;
  ints out = {-1, -1}; 
  float mindr = FLT_MAX;
  for( unsigned int i = 0; i < eta1.size(); i++ )
  {
    //cout<<"i = "<<i<<endl;	 
    for( unsigned int j = i+1; j < eta2.size(); j++ )
    {
      //cout<<"j = "<<j<<endl;
      if(jet_pt[i] != bjet_pt && jet_pt[j] != bjet_pt)
      {
	      auto dr = ROOT::VecOps::DeltaR(eta1[i],eta2[j],phi1[i],phi2[j]);
	      if( dr < mindr )
	      {
	        out = {i, j};
	        mindr = dr;
	        //cout<<"mindr = "<<mindr<<endl;
	      }
      }
    }
  }
  //cout<<"out = "<<out<<endl;
  return out;
}

float minDR_5 (floats &eta1, floats &eta2, floats &phi1, floats &phi2)
{
  // cout<<"G8G8G8"<<endl;
  float out = 100;
  float mindr = FLT_MAX;
  for( unsigned int i = 0; i < eta1.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    for( unsigned int j = 0; j < eta2.size(); j++ )
    {
      // cout<<"j = "<<j<<endl;
      auto dr = ROOT::VecOps::DeltaR(eta1[i],eta2[j],phi1[i],phi2[j]);
      if( dr < mindr )
      {
	      mindr = dr;
        out = mindr;
	      // cout<<"mindr = "<<mindr<<endl;
      }
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

float minDR_6 (double &eta1, floats &eta2, double &phi1, floats &phi2)
{
  // cout<<"G9G9G9"<<endl;
  float out = 100;
  float mindr = FLT_MAX;
  for( unsigned int i = 0; i < eta2.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    auto dr = ROOT::VecOps::DeltaR(float(eta1),eta2[i],float(phi1),phi2[i]);
    if( dr < mindr )
    {
	    mindr = dr;
      out = mindr;
	    // cout<<"mindr = "<<mindr<<endl;
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

int Max (floats &variable)
{
  // cout<<"G10G10G10"<<endl;
  float value = 0.;
  int out = -1;
  for( unsigned int i = 0; i < variable.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    auto temp = variable[i];
    if(temp > value)
    {
      value = temp;
      out = i;
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

ints mass_W (floats jet_pt, floats jet_eta, floats jet_phi, floats jet_mass)
{
  // cout<<"G11G11G11"<<endl;
  const double W_mass = 83.9;
  float value = FLT_MAX;
  FourVectorRVec sumW4Rvecs;
  FourVectorRVec jet1;
  FourVectorRVec jet2;
  ints out = {-1, -1};
  for( unsigned int i = 0; i < jet_pt.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    for( unsigned int j = i+1; j < jet_pt.size(); j++ )
    {
      // cout<<"j = "<<j<<endl;
      jet1 = generate_4Rvec_2(jet_pt[i],jet_eta[i],jet_phi[i],jet_mass[i]);
      jet2 = generate_4Rvec_2(jet_pt[j],jet_eta[j],jet_phi[j],jet_mass[j]);
      sumW4Rvecs = sum_two_4Rvecs(jet1, jet2);
      auto temp = sumW4Rvecs[0].M();
      // cout<<"temp = "<<temp<<endl;
      if(abs(temp - W_mass) < value)
      {
	      value = abs(temp - W_mass);
        out = {i,j};
	      // cout<<"out = "<<out<<endl;
        // cout<<"mass_W = "<<temp<<endl;
      }
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

ints mass_W_offshell (floats jet_pt, floats jet_eta, floats jet_phi, floats jet_mass)
{
  // cout<<"G12G12G12"<<endl;
  float value = FLT_MAX;
  FourVectorRVec sumW4Rvecs;
  FourVectorRVec jet1;
  FourVectorRVec jet2;
  ints out = {-1, -1};
  for( unsigned int i = 0; i < jet_pt.size(); i++ )
  {
    // cout<<"i = "<<i<<endl;
    for( unsigned int j = i+1; j < jet_pt.size(); j++ )
    {
      // cout<<"j = "<<j<<endl;
      jet1 = generate_4Rvec_2(jet_pt[i],jet_eta[i],jet_phi[i],jet_mass[i]);
      jet2 = generate_4Rvec_2(jet_pt[j],jet_eta[j],jet_phi[j],jet_mass[j]);
      sumW4Rvecs = sum_two_4Rvecs(jet1,jet2);
      auto temp = sumW4Rvecs[0].M();
      if(temp < value)
      {
	      value = temp;
        out = {i,j};
	      // cout<<"out = "<<out<<endl;
        // cout<<"mass_W = "<<value<<endl;
      }
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

bool false_mass_dilepton (FourVectorVec &leptons)
{
  // cout<<"H1H1H1"<<endl;
  bool out = false;
  float dilepton_mass = 0;

  for( unsigned int i = 0; i < leptons.size(); i++)
  {
    // cout<<"i = "<<i<<endl;
    for( unsigned int j = i+1; j < leptons.size(); j++)
    {
      // cout<<"j = "<<j<<endl;
      auto dilepton = leptons[i] + leptons[j];
      dilepton_mass = dilepton.M();
      // cout<<"dilepton_mass = "<<dilepton_mass<<endl;
      if(dilepton_mass > 80 && dilepton_mass < 100) out = true;
      break;
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

floats ratio_1 (floats &lepton_pt, floats &jet_pt)
{
  // cout<<"I1I1I1"<<endl;
  floats out;
  float ratio_pt = 0;

  for( unsigned int i = 0; i < lepton_pt.size(); i++)
  {
    // cout<<"i = "<<i<<endl;
    for( unsigned int j = 0; j < jet_pt.size(); j++)
    {
      // cout<<"j = "<<j<<endl;
      float ratio_pt = 0;
      ratio_pt = lepton_pt[i]/jet_pt[j];
      // cout<<"ratio_pt = "<<ratio_pt<<endl;
      out.emplace_back(ratio_pt);
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

floats ratio_2 (double &lepton_pt, floats &jet_pt)
{
  // cout<<"I2I2I2"<<endl;
  floats out;
  float ratio_pt = 0;

  for( unsigned int i = 0; i < jet_pt.size(); i++)
  {
    // cout<<"i = "<<i<<endl;
    float ratio_pt = 0;
    ratio_pt = float(lepton_pt)/jet_pt[i];
    // cout<<"ratio_pt = "<<ratio_pt<<endl;
    out.emplace_back(ratio_pt);
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

float minratio_1 (floats &lepton_pt, floats &jet_pt)
{
  // cout<<"I3I3I3"<<endl;
  float out = 100;
  float minratio_pt = FLT_MAX;

  for( unsigned int i = 0; i < lepton_pt.size(); i++)
  {
    // cout<<"i = "<<i<<endl;
    for( unsigned int j = 0; j < jet_pt.size(); j++)
    {
      // cout<<"j = "<<j<<endl;
      auto ratio_pt = lepton_pt[i]/jet_pt[j];
      if( ratio_pt < minratio_pt )
      {
        minratio_pt = ratio_pt;
        out = minratio_pt;
        // cout<<"minratio_pt = "<<minratio_pt<<endl;
      }
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

float minratio_2 (double &lepton_pt, floats &jet_pt)
{
  // cout<<"I4I4I4"<<endl;
  float out = 100;
  float minratio_pt = FLT_MAX;

  for( unsigned int i = 0; i < jet_pt.size(); i++)
  {
    // cout<<"i = "<<i<<endl;
    auto ratio_pt = float(lepton_pt)/jet_pt[i];
    if( ratio_pt < minratio_pt )
    {
      minratio_pt = ratio_pt;
      out = minratio_pt;
      // cout<<"minratio_pt = "<<minratio_pt<<endl;
    }
  }
  // cout<<"out = "<<out<<endl;
  return out;
}

FourVectorRVec removebjet(FourVectorRVec jet4Rvecs, floats jetpt, double bjetpt)
{
  // cout<<"J1J1J1"<<endl;
  unsigned int nj = jetpt.size();
  unsigned int index = 0;
  bool bjet = false;
  FourVectorRVec temp;
  temp.reserve(1);
  for( unsigned int i = 0; i < nj; i++)
  {
    // cout<<"i = "<<i<<endl;
    if(jetpt[i] == (float) bjetpt)
    {
      // cout<<"The bjet coming from the top has the index "<<i<<"."<<endl;
      index = i;
      bjet = true;
    }
    if(bjet == true && index != i)
    {
      temp[0] = jet4Rvecs[i];
      jet4Rvecs[i] = jet4Rvecs[index];
      jet4Rvecs[index] = temp[0];
      index = i;
    }
  }
  // temp = jet4Rvecs[nj];
  // jet4Rvecs[nj] = jet4Rvecs[index];
  // jet4Rvecs[index] = temp;
  jet4Rvecs.pop_back();
  // cout<<"jet4Rvecs = "<<jet4Rvecs<<endl;

  return jet4Rvecs;
}

RVec< RVec<int> > distinct_comb_2(int nl, int nj)
{
  // cout<<"J3J3J3"<<endl;
  // cout<<"nl = "<<nl<<endl;
  // cout<<"nj = "<<nj<<endl;
  RVec<int> idx1;
  RVec<int> idx2;
  RVec< RVec<int> > res;
  idx1.reserve(nl*(nl-1)/2);
  idx2.reserve(nj*(nj-1)/2);

  if(nl < 2 || nj < 1)
  {
    idx1.emplace_back(0);
    idx2.emplace_back(0);
  }  

  for(unsigned int i = 0; i < nl; i++)
  {
    for(unsigned int j = 0; j < nj; j++)
    {
	    idx1.emplace_back(i);
	    idx2.emplace_back(j);
    }
  }      
  // cout<<"idx1 = "<<idx1<<", idx2 = "<<idx2<<endl;
  res.emplace_back(idx1);
  res.emplace_back(idx2);
  return res;      
};

FourVectorRVec sum_two_4Rvecs (FourVectorRVec vec1, FourVectorRVec vec2)
{
  // cout<<"J4J4J4"<<endl;
  FourVectorRVec vec;
  vec = vec1 + vec2;
  // cout<<"vec = "<<vec<<endl;
  return vec;
};

RVec< RVec<int> > distinct_comb_2_bis(int nj)
{
  // cout<<"J5J5J5"<<endl;
  // cout<<"nj = "<<nj<<endl;
  RVec<int> idx1;
  RVec<int> idx2;
  RVec< RVec<int> > res;
  idx1.reserve(nj*(nj-1)/2);
  idx2.reserve(nj*(nj-1)/2);

  if(nj < 2)
  {
    idx1.emplace_back(0);
    idx2.emplace_back(0);
  }  

  for(unsigned int i = 0; i < nj; i++)
  {
    for(unsigned int j = i+1; j < nj; j++)
    {
	    idx1.emplace_back(i);
	    idx2.emplace_back(j);
    }
  }      
  // cout<<"idx1 = "<<idx1<<", idx2 = "<<idx2<<endl;
  res.emplace_back(idx1);
  res.emplace_back(idx2);
  return res;      
};

FourVectorRVec reconstruction_neutrino(double leptontoppt, double leptontopeta, double leptontopphi, double leptontopmass, double bjetpt, double bjeteta, double bjetphi, double bjetmass, double leptonHiggspt, double leptonHiggseta, double leptonHiggsphi, double leptonHiggsmass, float MET_px, float MET_py)
{
  // cout<<"K1K1K1"<<endl;
  FourVectorRVec vecs;
  const double PI = acos(-1.0);
  const double W_mass = 83.9;
  const double W_width = 10.8;
  const double top_mass = 175.9;
  const double top_width = 17.2;
  double pt = 0.;
  double eta = 0.;
  double phi = 0.;
  double px = 0.;
  double py = 0.;
  // double ptp = 0.;
  // double etap = 0.;
  // double phip = 0.;
  // double pxp = 0.;
  // double pyp = 0.;
  double mass_W_temp = FLT_MAX;
  // double mass_Wp_temp = FLT_MAX;
  double mass_top_temp = FLT_MAX;
  double chi2_W_temp = FLT_MAX;
  double chi2_top_temp = FLT_MAX;
  double mass_W = FLT_MAX;
  double mass_top = FLT_MAX;
  double chi2_W = FLT_MAX;
  double chi2_top = FLT_MAX;
  double nutoppt = 0.;
  double nutopeta = 0.;
  double nutopphi = 0.;
  FourVectorRVec nutop4Rvecs;
  // FourVectorRVec nuHiggs4Rvecs;
  FourVectorRVec sumW4Rvecs;
  // FourVectorRVec sumWp4Rvecs;
  FourVectorRVec sumtop4Rvecs;
  double ptmin = 0.;
  double ptmax = 150.;
  double ptstep = (ptmax-ptmin)/1;
  double etamin = leptontopeta-1.5;
  double etamax = leptontopeta+1.5;
  double etastep = (etamax-etamin)/2;
  // cout<<"etamin = "<<etamin<<", etamax = "<<etamax<<", etastep = "<<etastep<<endl;
  double phimin = leptontopphi-PI/2;
  double phimax = leptontopphi+PI/2;
  double phistep = (phimax-phimin)/1;
  // cout<<"phimin = "<<phimin<<", phimax = "<<phimax<<", phistep = "<<phistep<<endl;
  FourVectorRVec leptontop4Rvecs = generate_4Rvec_2(leptontoppt, leptontopeta, leptontopphi, leptontopmass);
  // FourVectorRVec leptonHiggs4Rvecs = generate_4Rvec_2(leptonHiggspt, leptonHiggseta, leptonHiggsphi, leptonHiggsmass);
  FourVectorRVec bjet4Rvecs = generate_4Rvec_2(bjetpt, bjeteta, bjetphi, bjetmass);

  for(double i = ptmin; i < ptmax; i += ptstep)
  { 
    pt = i;
    cout<<"pt = "<<pt<<endl;
    for(double j = etamin; j < etamax; j += etastep)
    {
      eta = j;
      cout<<"eta = "<<eta<<endl;
      for(double k = phimin; k < phimax; k += phistep)
      {
        phi = k;
        // cout<<"phi = "<<phi<<endl;
        if(phi < -PI)
        {
          phi += 2*PI;
        }
        else if(phi > PI)
        {
          phi -= 2*PI;
        }
        cout<<"phi = "<<phi<<endl;
        nutop4Rvecs = generate_4Rvec_4(pt, eta, phi);
        sumW4Rvecs = sum_two_4Rvecs(leptontop4Rvecs, nutop4Rvecs);
        sumtop4Rvecs = sum_two_4Rvecs(sumW4Rvecs, bjet4Rvecs);
        mass_W_temp = sumW4Rvecs[0].M();
        mass_top_temp = sumtop4Rvecs[0].M();
        chi2_W_temp = ((mass_W_temp-W_mass)/W_width)*((mass_W_temp-W_mass)/W_width);
        chi2_top_temp = ((mass_top_temp-top_mass)/top_width)*((mass_top_temp-top_mass)/top_width);
        // cout<<"mass_W_temp = "<<mass_W_temp<<", chi2_W_temp = "<<chi2_W_temp<<endl;
        // cout<<"mass_top_temp = "<<mass_top_temp<<", chi2_top_temp = "<<chi2_top_temp<<endl;
        // px = pt*cos(phi);
        // py = pt*sin(phi);
        // pxp = MET_px - px;
        // pyp = MET_py - py;
        // ptp = sqrt(pxp*pxp + pyp*pyp);
        // phip = atan2(pyp,pxp);
        // etap = eta_neutrino_Higgs_2(leptontopeta, leptontopphi, leptonHiggseta, leptonHiggsphi, eta, phi, phip);
        // nuHiggs4Rvecs = generate_4Rvec_4(ptp, etap, phip);
        // sumWp4Rvecs = sum_two_4Rvecs(leptonHiggs4Rvecs, nuHiggs4Rvecs);
        // mass_Wp_temp = sumWp4Rvecs[0].M();
        // cout<<"mass_Wp_temp = "<<mass_Wp_temp<<endl;
        // if(abs(mass_W_temp-W_mass) < abs(mass_W-W_mass) && abs(mass_W_temp-W_mass) < 20 && abs(mass_top_temp-top_mass) < abs(mass_top-top_mass) && abs(mass_top_temp-top_mass) < 50)
        // if(chi2_W_temp < chi2_W && chi2_top_temp < chi2_top && abs(mass_W_temp-W_mass) < 20 && abs(mass_top_temp-top_mass) < 50)
        if(chi2_W_temp < chi2_W && chi2_top_temp < chi2_top)
        {
          mass_W = mass_W_temp;
          mass_top = mass_top_temp;
          chi2_W = chi2_W_temp;
          chi2_top = chi2_top_temp;
          nutoppt = pt;
          nutopeta = eta;
          nutopphi = phi;
          vecs = nutop4Rvecs;
          cout<<"chi2_W = "<<chi2_W<<", chi2_top = "<<chi2_top<<endl;
          cout<<"W_pt = "<<sumW4Rvecs[0].Pt()<<", W_eta = "<<sumW4Rvecs[0].Eta()<<", W_phi = "<<sumW4Rvecs[0].Phi()<<", W_mass = "<<mass_W<<endl;
          cout<<"top_pt = "<<sumtop4Rvecs[0].Pt()<<", top_eta = "<<sumtop4Rvecs[0].Eta()<<", top_phi = "<<sumtop4Rvecs[0].Phi()<<", top_mass = "<<mass_top<<endl;
          cout<<"pt = "<<nutoppt<<", eta = "<<nutopeta<<", phi = "<<nutopphi<<endl;
        }
      }
    }
  }
  // cout<<"vecs = "<<vecs<<endl;
  return vecs;
}

double eta_neutrino_Higgs(double leptontopeta, double leptontopphi, double leptonHiggseta, double leptonHiggsphi, floats &nutopeta, floats &nutopphi, floats &nuHiggsphi)
{
  // cout<<"K2K2K2"<<endl;
  double nuHiggseta1;
  double nuHiggseta2;
  double leptondeltaR = ROOT::VecOps::DeltaR((float) leptontopeta, (float) leptonHiggseta, (float) leptontopphi, (float) leptonHiggsphi);
  double nudeltaPhi = abs(ROOT::VecOps::DeltaPhi((float) nutopphi[0], (float) nuHiggsphi[0]));
  double nudeltaEta = sqrt(abs(leptondeltaR*leptondeltaR - nudeltaPhi*nudeltaPhi));
  // cout<<"leptondeltaR = "<<leptondeltaR<<", nudeltaPhi = "<<nudeltaPhi<<", nudeltaEta = "<<nudeltaEta<<endl;
  nuHiggseta1 = nudeltaEta + (double) nutopeta[0];
  nuHiggseta2 = nudeltaEta - (double) nutopeta[0];
  // cout<<"nuHiggseta1 = "<<nuHiggseta1<<", nuHiggseta2 = "<<nuHiggseta2<<endl;

  if(abs(nuHiggseta1) < abs(nuHiggseta2))
  {
    return nuHiggseta1;
  }
  else
  {
    return nuHiggseta2;
  }
}

double eta_neutrino_Higgs_2(double leptontopeta, double leptontopphi, double leptonHiggseta, double leptonHiggsphi, double nutopeta, double nutopphi, double nuHiggsphi)
{
  // cout<<"K3K3K3"<<endl;
  double nuHiggseta1;
  double nuHiggseta2;
  double leptondeltaR = ROOT::VecOps::DeltaR((float) leptontopeta, (float) leptonHiggseta, (float) leptontopphi, (float) leptonHiggsphi);
  double nudeltaPhi = abs(ROOT::VecOps::DeltaPhi((float) nutopphi, (float) nuHiggsphi));
  double nudeltaEta = sqrt(abs(leptondeltaR*leptondeltaR - nudeltaPhi*nudeltaPhi));
  // cout<<"leptondeltaR = "<<leptondeltaR<<", nudeltaPhi = "<<nudeltaPhi<<", nudeltaEta = "<<nudeltaEta<<endl;
  nuHiggseta1 = nudeltaEta + (double) nutopeta;
  nuHiggseta2 = nudeltaEta - (double) nutopeta;
  // cout<<"nuHiggseta1 = "<<nuHiggseta1<<", nuHiggseta2 = "<<nuHiggseta2<<endl;

  if(abs(nuHiggseta1) < abs(nuHiggseta2))
  {
    return nuHiggseta1;
  }
  else
  {
    return nuHiggseta2;
  }
}

// double Basic_Mt2_332_Calculator::mT2Fcn::operator()(const vector<double>& par) const
// {
//   double qT1x = par[0];
//   double qT1y = par[1];
    
//   double qT2x = theExmiss - qT1x;
//   double qT2y = theEymiss - qT1y;
    
//   double ETj1 = sqrt(theMass1*theMass1 +  thePT1x*thePT1x + thePT1y*thePT1y);    // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
//   double ETj2 = sqrt(theMass2*theMass2 +  thePT2x*thePT2x + thePT2y*thePT2y);     // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless! 
    
//   double ETchi1 = sqrt(qT1x*qT1x + qT1y*qT1y + theMchi*theMchi);
//   double ETchi2 = sqrt(qT2x*qT2x + qT2y*qT2y + theMchi*theMchi);
    
//   double mTsq1 = theMass1*theMass1 + theMchi*theMchi + 2.0*(ETj1*ETchi1 - (thePT1x*qT1x + thePT1y*qT1y)); // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
//   double mTsq2 = theMass2*theMass2 + theMchi*theMchi + 2.0*(ETj2*ETchi2 - (thePT2x*qT2x + thePT2y*qT2y)); // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!

//     // cout << "MOO " << sqrt(fmax(mTsq1,mTsq2)) << "\t" << qT1x << "\t" << qT1y << "\t" << 0.5*(qT1x-qT2x) << "\t" << 0.5*(qT1y-qT2y ) << endl;
 
//   return fmax(mTsq1,mTsq2);
    
// }

// double mt2_332_Sq(FourVectorRVec visA, FourVectorRVec visB, float ptmiss, float phimiss, const double mEachInvisible)
// {
//   float pxmiss = ptmiss*cos(phimiss);
//   float pymiss = ptmiss*sin(phimiss);

//   mT2Fcn theFCN(pxmiss, pymiss, mEachInvisible, visA.Px(), visA.Py(), visA.M(), visB.Px(), visB.Py(), visB.M());

//   const double massScale = (ptmiss.Pt() + mEachInvisible + visA.Pt() + visA.M() + visB.Pt() + visB.M())/6.0;
//   // DANG! Try to get rid of Minuit output:
//   //ofstream    DANG_log("/dev/null");
//   //streambuf * DANG_save = cerr.rdbuf();
//   //if (DANG_log.is_open()) {
//   //  cerr.rdbuf(DANG_log.rdbuf());
//   // }

//   double guessx = 0.5*pxmiss;
//   double guessy = 0.5*pymiss;
    
//   MnUserParameters upar;
//   upar.Add("etx", guessx, 0.02*massScale); 
//   upar.Add("ety", guessy, 0.02*massScale);
//   const int highQuality=2;    

//   // Usually migrad produces the best minumum.
//   // But when the minimum is in a fold, migrad can fail badly.
//   // On the fold, simplex does well.  We therefore do both separately
//   // and record the answer of the one that did best.
    
//   // Further to the above notes, it now seems that by choosing the massScale sensibly, and by making the "tolerance" (the second argument to the call that extracts the FunctionMinimum below) much smaller, the simplex algorithm seems to work so well that we don't need the migrad algorithm.  This is good news, as the migrad algorithm produces lots of error output that we can't get rid of.

//   MnSimplex simplex(theFCN, upar, highQuality);
//   FunctionMinimum minS = simplex(0,massScale*0.000001);

//   etxAt = minS.UserState().Value("etx");
//   etyAt = minS.UserState().Value("ety");
//   etxBt = pxmiss-etxAt;
//   etyBt = pymiss-etyAt;
//   //MnMigrad migrad(theFCN, upar, highQuality);
//   //FunctionMinimum minM = migrad(0,massScale*0.000001);
//   //const double best = fmin(minS.Fval(), minM.Fval());
//   const double best = minS.Fval();

//   // DANG! Undoing our attempt to get rid of Minuit output:
//   //if (DANG_log.is_open()) {
//   //  cerr.rdbuf(DANG_save);
//   //}

//   return best;    
// }

// double mt2_332(FourVectorRVec visA, FourVectorRVec visB, float ptmiss, float phimiss, const double mEachInvisible)
// {
//   const double mSq = mt2_332_Sq(visA, visB, ptmiss, phimiss, mEachInvisible);

//   if (mSq>=0)
//   {
//     return sqrt(mSq);
//   }
//   else
//   {
//     // something went wrong -- maybe just rounding errors;
//     return sqrt(fabs(mSq));
//   }
// }

// void solve(double* ETmiss, FourVectorVec b, double* bb, double* lp, double* lm, double mWp, double mWm, double mt, double mtb, double mnu, double mnub, vector<double> *nupx, vector<double> *nupy, vector<double> *nupz, vector<double> *nubpx, vector<double> *nubpy, vector<double> *nubpz, vector<double> *cd_diff, int& cubic_single_root_cmplx)
// floats solve(float metpx, float metpy, floats &be, floats &bpx, floats &bpy, floats &bpz, floats &je, floats &jpx, floats &jpy, floats &jpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &mj, floats &ml, floats &nue, floats &nupx, floats &nupy, floats &nupz, floats &nupt, floats &nueta, floats &nuphi, floats &num, floats &wmasses)
floats solve(float metpx, float metpy, floats &be, floats &bpx, floats &bpy, floats &bpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &ml, floats &nue, floats &nupx, floats &nupy, floats &nupz, floats &nupt, floats &nueta, floats &nuphi, floats &num, floats &wtmass, floats &watmass)
{

  cout<<"L1L1L1"<<endl;
  // // cubic_single_root_cmplx = 0;
  // double radic = be[0]*be[0]-bpx[0]*bpx[0]-bpy[0]*bpy[0]-bpz[0]*bpz[0];
  // double mb = 0.;
  // if (radic > 0.) mb = sqrt(radic);
  // cout<<"mb = "<<mb<<endl;

  // // radic = sqr(bb[0])-sqr(bb[1])-sqr(bb[2])-sqr(bb[3]);
  // // double mbb = 0.;
  // // if (radic > 0.) mbb = sqrt(radic);

  // radic = le[0]*le[0]-lpx[0]*lpx[0]-lpy[0]*lpy[0]-lpz[0]*lpz[0];
  // double mlp = 0.;
  // if (radic > 0.) mlp = sqrt(radic);
  // cout<<"mlp = "<<mlp<<endl;

  // radic = le[1]*le[1]-lpx[1]*lpx[1]-lpy[1]*lpy[1]-lpz[1]*lpz[1];
  // double mlm = 0.;
  // if (radic > 0.) mlm = sqrt(radic);
  // cout<<"mlm = "<<mlm<<endl;

  double epsilon = 0.000001;
  // float wmass = wmasses[0];
  // cout<<"wmass = "<<wmass<<endl;
  float wtopmass = wtmass[0];
  float wantitopmass = watmass[0];
  cout<<"wtopmass = "<<wtopmass<<", wantitopmass = "<<wantitopmass<<endl;

  double a1 = (be[0]+le[0])*(wtopmass*wtopmass-ml[0]*ml[0])-le[0]*(172.69*172.69-mb[0]*mb[0]-ml[0]*ml[0])+2.*be[0]*le[0]*le[0]
               -2.*le[0]*(bpx[0]*lpx[0]+bpy[0]*lpy[0]+bpz[0]*lpz[0]);
  double a2 = 2.*(be[0]*lpx[0]-le[0]*bpx[0]);
  double a3 = 2.*(be[0]*lpy[0]-le[0]*bpy[0]);
  double a4 = 2.*(be[0]*lpz[0]-le[0]*bpz[0]);
  cout<<"a1 = "<<a1<<", a2 = "<<a2<<", a3 = "<<a3<<", a4 = "<<a4<<endl;

  //c00*nupy^2+(c10*nupx+c11)*nupy+(c20*nupx^2+c21*nupx+c22)
  double c22 = (wtopmass*wtopmass-ml[0]*ml[0])*(wtopmass*wtopmass-ml[0]*ml[0])*a4*a4-4.*(le[0]*le[0]-lpz[0]*lpz[0])*a1*a1
               -4.*(wtopmass*wtopmass-ml[0]*ml[0])*lpz[0]*a1*a4; //c1
  double c21 = -8.*(le[0]*le[0]-lpz[0]*lpz[0])*a1*a2+4.*(wtopmass*wtopmass-ml[0]*ml[0])*(lpx[0]*a4*a4-lpz[0]*a2*a4)-8.*lpx[0]*lpz[0]*a1*a4; //c2
  double c11 = -8.*(le[0]*le[0]-lpz[0]*lpz[0])*a1*a3+4.*(wtopmass*wtopmass-ml[0]*ml[0])*(lpy[0]*a4*a4-lpz[0]*a3*a4)-8.*lpy[0]*lpz[0]*a1*a4; //c3
  double c20 = -4.*(le[0]*le[0]-lpx[0]*lpx[0])*a4*a4-4.*(le[0]*le[0]-lpz[0]*lpz[0])*a2*a2-8.*lpx[0]*lpz[0]*a2*a4; //c4
  double c10 = -8.*(le[0]*le[0]-lpz[0]*lpz[0])*a2*a3+8.*lpx[0]*lpy[0]*a4*a4-8.*lpx[0]*lpz[0]*a3*a4-8.*lpy[0]*lpz[0]*a2*a4; //c5
  double c00 = -4.*(le[0]*le[0]-lpy[0]*lpy[0])*a4*a4-4.*(le[0]*le[0]-lpz[0]*lpz[0])*a3*a3-8.*lpy[0]*lpz[0]*a3*a4; //c6
  cout<<"c22 = "<<c22<<", c21 = "<<c21<<", c11 = "<<c11<<", c20 = "<<c20<<", c10 = "<<c10<<", c00 = "<<c00<<endl;

  double b1 = (be[1]+le[1])*(wantitopmass*wantitopmass-ml[1]*ml[1])-le[1]*(172.69*172.69-mb[1]*mb[1]-ml[1]*ml[1])+2.*be[1]*le[1]*le[1]
               -2.*le[1]*(bpx[1]*lpx[1]+bpy[1]*lpy[1]+bpz[1]*lpz[1]);
  double b2 = 2.*(be[1]*lpx[1]-le[1]*bpx[1]);
  double b3 = 2.*(be[1]*lpy[1]-le[1]*bpy[1]);
  double b4 = 2.*(be[1]*lpz[1]-le[1]*bpz[1]);
  cout<<"b1 = "<<b1<<", b2 = "<<b2<<", b3 = "<<b3<<", b4 = "<<b4<<endl;

  double dp22 = (wantitopmass*wantitopmass-ml[1]*ml[1])*(wantitopmass*wantitopmass-ml[1]*ml[1])*b4*b4-4.*(le[1]*le[1]-lpz[1]*lpz[1])*b1*b1
               -4.*(wantitopmass*wantitopmass-ml[1]*ml[1])*lpz[1]*b1*b4; //d1
  double dp21 = -8.*(le[1]*le[1]-lpz[1]*lpz[1])*b1*b2+4.*(wantitopmass*wantitopmass-ml[1]*ml[1])*(lpx[1]*b4*b4-lpz[1]*b2*b4)-8.*lpx[1]*lpz[1]*b1*b4; //d2
  double dp11 = -8.*(le[1]*le[1]-lpz[1]*lpz[1])*b1*b3+4.*(wantitopmass*wantitopmass-ml[1]*ml[1])*(lpy[1]*b4*b4-lpz[1]*b3*b4)-8.*lpy[1]*lpz[1]*b1*b4; //d3
  double dp20 = -4.*(le[1]*le[1]-lpx[1]*lpx[1])*b4*b4-4.*(le[1]*le[1]-lpz[1]*lpz[1])*b2*b2-8.*lpx[1]*lpz[1]*b2*b4; //d4
  double dp10 = -8.*(le[1]*le[1]-lpz[1]*lpz[1])*b2*b3+8.*lpx[1]*lpy[1]*b4*b4-8.*lpx[1]*lpz[1]*b3*b4-8.*lpy[1]*lpz[1]*b2*b4; //d5
  double dp00 = -4.*(le[1]*le[1]-lpy[1]*lpy[1])*b4*b4-4.*(le[1]*le[1]-lpz[1]*lpz[1])*b3*b3-8.*lpy[1]*lpz[1]*b3*b4; //d6
  cout<<"dp22 = "<<dp22<<", dp21 = "<<dp21<<", dp11 = "<<dp11<<", dp20 = "<<dp20<<", dp10 = "<<dp10<<", dp00 = "<<dp00<<endl;

  // double b1 = (je[0]+je[1]+le[1])*(wmass*wmass-ml[1]*ml[1])-le[1]*(125.25*125.25-mj[0]*mj[0]-mj[1]*mj[1]-ml[1]*ml[1])
  //             +2.*le[1]*(le[1]*(je[0]+je[1])+je[0]*je[1])
  //             -2.*le[1]*(lpx[1]*(jpx[0]+jpx[1])+lpy[1]*(jpy[0]+jpy[1])+lpz[1]*(jpz[0]+jpz[1])+jpx[0]*jpx[1]+jpy[0]*jpy[1]+jpz[0]*jpz[1]);
  // double b2 = 2.*((je[0]+je[1])*lpx[1]-le[1]*(jpx[0]+jpx[1]));
  // double b3 = 2.*((je[0]+je[1])*lpy[1]-le[1]*(jpy[0]+jpy[1]));
  // double b4 = 2.*((je[0]+je[1])*lpz[1]-le[1]*(jpz[0]+jpz[1]));
  // cout<<"b1 = "<<b1<<", b2 = "<<b2<<", b3 = "<<b3<<", b4 = "<<b4<<endl;

  // double dp22 = (wmass*wmass-ml[1]*ml[1])*(wmass*wmass-ml[1]*ml[1])*b4*b4-4.*((le[1])*(le[1])-lpz[1]*lpz[1])*b1*b1
  //               -4.*(wmass*wmass-ml[1]*ml[1])*lpz[1]*b1*b4; //d1
  // double dp21 = -8.*(le[1]*le[1]-lpz[1]*lpz[1])*b1*b2+4.*(wmass*wmass-ml[1]*ml[1])*(lpx[1]*b4*b4-lpz[1]*b2*b4)-8.*lpx[1]*lpz[1]*b1*b4; //d2
  // double dp11 = -8.*(le[1]*le[1]-lpz[1]*lpz[1])*b1*b3+4.*(wmass*wmass-ml[1]*ml[1])*(lpy[1]*b4*b4-lpz[1]*b3*b4)-8.*lpy[1]*lpz[1]*b1*b4; //d3
  // double dp20 = -4.*(le[1]*le[1]-lpx[1]*lpx[1])*b4*b4-4.*(le[1]*le[1]-lpz[1]*lpz[1])*b2*b2-8.*lpx[1]*lpz[1]*b2*b4; //d4
  // double dp10 = -8.*(le[1]*le[1]-lpz[1]*lpz[1])*b2*b3+8.*lpx[1]*lpy[1]*b4*b4-8.*lpx[1]*lpz[1]*b3*b4-8.*lpy[1]*lpz[1]*b2*b4; //d5
  // double dp00 = -4.*(le[1]*le[1]-lpy[1]*lpy[1])*b4*b4-4.*(le[1]*le[1]-lpz[1]*lpz[1])*b3*b3-8.*lpy[1]*lpz[1]*b3*b4; //d6
  // cout<<"dp22 = "<<dp22<<", dp21 = "<<dp21<<", dp11 = "<<dp11<<", dp20 = "<<dp20<<", dp10 = "<<dp10<<", dp00 = "<<dp00<<endl;

  //d00*nupy^2+(d10*nupx+d11)*nupy+(d20*nupx^2+d21*nupx+d22)
  double d22 = dp22+metpx*metpx*dp20+metpy*metpy*dp00+metpx*metpy*dp10+metpx*dp21+metpy*dp11;
  double d20 = dp20;
  double d00 = dp00;
  double d10 = dp10;
  double d21 = -dp21-2.*metpx*dp20-metpy*dp10;
  double d11 = -dp11-2.*metpy*dp00-metpx*dp10;
  cout<<"d22 = "<<d22<<", d21 = "<<d21<<", d11 = "<<d11<<", d20 = "<<d20<<", d10 = "<<d10<<", d00 = "<<d00<<endl;

  //                  |c0    d0   |
  // resultant(nupy)= |c1 c0 d1 d0|
  //                  |c2 c1 d2 d1|
  //                  |   c2    d2|
  //
  //expressions c0,1,2,d0,1,2 are polynomials in nupx of degree 2, to be multiplied out below

  vector<double> polx(5);

  //publication formulae
  polx[0] = c00*c00*d22*d22+c11*d22*(c11*d00-c00*d11)+c00*c22*(d11*d11-2.*d00*d22)+c22*d00*(c22*d00-c11*d11); //x^0

  polx[1] = c00*d21*(2.*c00*d22-c11*d11)+c00*d11*(2.*c22*d10+c21*d11)+c22*d00*(2.*c21*d00-c11*d10)-c00*d22*(c11*d10+c10*d11)
            -2.*c00*d00*(c22*d21+c21*d22)-d00*d11*(c11*c21+c10*c22)+c11*d00*(c11*d21+2.*c10*d22); //X^1

  polx[2] = c00*c00*(2.*d22*d20+d21*d21)-c00*d21*(c11*d10+c10*d11)+c11*d20*(c11*d00-c00*d11)+c00*d10*(c22*d10-c10*d22)
            +c00*d11*(2.*c21*d10+c20*d11)+d00*d00*(2.*c22*c20+c21*c21)-2.*c00*d00*(c22*d20+c21*d21+c20*d22)
            +c10*d00*(2.*c11*d21+c10*d22)-d00*d10*(c11*c21+c10*c22)-d00*d11*(c11*c20+c10*c21); //x^2

  polx[3] = c00*d21*(2.*c00*d20-c10*d10)-c00*d20*(c11*d10+c10*d11)+c00*d10*(c21*d10+2.*c20*d11)-2.*c00*d00*(c21*d20+c20*d21)
            +c10*d00*(2.*c11*d20+c10*d21)+c20*d00*(2.*c21*d00-c10*d11)-d00*d10*(c11*c20+c10*c21); //x^3 //sign typo for the sixth term

  polx[4] = c00*c00*d20*d20+c10*d20*(c10*d00-c00*d10)+c20*d10*(c00*d10-c10*d00)+c20*d00*(c20*d00-2.*c00*d20); //x^4
  cout<<"polx[0] = "<<polx[0]<<", polx[1] = "<<polx[1]<<", polx[2] = "<<polx[2]<<", polx[3] = "<<polx[3]<<", polx[4] = "<<polx[4]<<endl;

  vector<double> nupxt;
  int cubic_single_root_cmplx = 0;
  quartic(polx, &nupxt, cubic_single_root_cmplx);

  double c0 = c00;
  double c1, c2;
  double d0 = d00;
  double d1, d2;
  cout<<"size of nupxt = "<<nupxt.size()<<endl;
  for (size_t i=0; i<nupxt.size(); ++i)
  {
    c1 = c10*nupxt[i]+c11;
    c2 = c20*nupxt[i]*nupxt[i]+c21*nupxt[i]+c22;
    d1 = d10*nupxt[i]+d11;
    d2 = d20*nupxt[i]*nupxt[i]+d21*nupxt[i]+d22;
    double denom = c1*d0-c0*d1;
    cout<<"c1 = "<<c1<<", c2 = "<<c2<<", d1 = "<<d1<<", d2 = "<<d2<<", denom = "<<denom<<endl;

    // (*cd_diff).push_back(denom); //should never happen
    if (fabs(denom) < epsilon) continue;

    double lpbz_diff = le[0]*bpz[0]-be[0]*lpz[0];
    // double lmbbz_diff = le[1]*(jpz[1]+jpz[2])-(je[1]+je[2])*lpz[1];
    double lmbbz_diff = le[1]*bpz[1]-be[1]*lpz[1];
    double thisnupx = nupxt[i];
    double thisnupy = (c0*d2-c2*d0)/denom;
    double thisnubpx = metpx-nupxt[i];
    double thisnubpy = metpy-thisnupy;
    double thisnupz, thisnubpz;

    int error=0;
    if (fabs(lpbz_diff)<epsilon)
    { //circumvent singularity
      // error = (be, bpx, bpy, bpz, je, jpx, jpy, jpz, le, lpx, lpy, lpz, mb, mj, ml, thisnupx, thisnupy, &thisnupz, epsilon);
      error = algebraic_pz_1(be, bpx, bpy, bpz, le, lpx, lpy, lpz, mb, ml, thisnupx, thisnupy, &thisnupz, epsilon, wtopmass);
      if (error==-1) continue; //next nupx, nubpx solutions
    }
    else
    {
      thisnupz = (-a1-a2*thisnupx-a3*thisnupy)/a4;
    }

    if (fabs(lmbbz_diff)<epsilon)
    { //circumvent singularity
      // error = algebraic_pz_2(be, bpx, bpy, bpz, je, jpx, jpy, jpz, le, lpx, lpy, lpz, mb, mj, ml, thisnupx, thisnupy, &thisnupz, epsilon);
      error = algebraic_pz_2(be, bpx, bpy, bpz, le, lpx, lpy, lpz, mb, ml, thisnupx, thisnupy, &thisnupz, epsilon, wantitopmass);
      if (error == -1) continue; //next nupx, nubpx solutions
    }
    else
    {
      thisnubpz = (-b1-b2*thisnubpx-b3*thisnubpy)/b4;
    }

    // (*nupx).push_back(thisnupx);
    // (*nupy).push_back(thisnupy);
    // (*nupz).push_back(thisnupz);
    // (*nupx).push_back(thisnubpx);
    // (*nupy).push_back(thisnubpy);
    // (*nupz).push_back(thisnubpz);
    // nupx.push_back(thisnupx);
    // nupy.push_back(thisnupy);
    // nupz.push_back(thisnupz);
    // nupx.push_back(thisnubpx);
    // nupy.push_back(thisnubpy);
    // nupz.push_back(thisnubpz);
    // nupx[2*i] = thisnupx;
    // nupy[2*i] = thisnupy;
    // nupz[2*i] = thisnupz;
    // nupt[2*i] = sqrt(thisnupx*thisnupx + thisnupy*thisnupy);
    // nueta[2*i] = asinh(thisnupz/sqrt(thisnupx*thisnupx + thisnupy*thisnupy));
    // nuphi[2*i] = acos(thisnupx/sqrt(thisnupx*thisnupx + thisnupy*thisnupy));
    // nue[2*i] = sqrt(thisnupx*thisnupx + thisnupy*thisnupy + thisnupz*thisnupz);
    // cout<<"nupx["<<2*i<<"] = "<<nupx[2*i]<<", nupy["<<2*i<<"] = "<<nupy[2*i]<<", nupz["<<2*i<<"] = "<<nupz[2*i]<<", nupt["<<2*i<<"] = "<<nupt[2*i]<<", nueta["<<2*i<<"] = "<<nueta[2*i]<<", nuphi["<<2*i<<"] = "<<nuphi[2*i]<<", nue["<<2*i<<"] = "<<nue[2*i]<<endl;
    // nupx[2*i+1] = thisnubpx;
    // nupy[2*i+1] = thisnubpy;
    // nupz[2*i+1] = thisnubpz;
    // nupt[2*i+1] = sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy);
    // nueta[2*i+1] = asinh(thisnubpz/sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy));
    // nuphi[2*i+1] = acos(thisnubpx/sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy));
    // nue[2*i+1] = sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy + thisnubpz*thisnubpz);
    // cout<<"nupx["<<2*i+1<<"] = "<<nupx[2*i+1]<<", nupy["<<2*i+1<<"] = "<<nupy[2*i+1]<<", nupz["<<2*i+1<<"] = "<<nupz[2*i+1]<<", nupt["<<2*i+1<<"] = "<<nupt[2*i+1]<<", nueta["<<2*i+1<<"] = "<<nueta[2*i+1]<<", nuphi["<<2*i+1<<"] = "<<nuphi[2*i+1]<<", nue["<<2*i+1<<"] = "<<nue[2*i+1]<<endl;
    
    if(i == 0)
    {
      nupx[0] = thisnupx;
      nupy[0] = thisnupy;
      nupz[0] = thisnupz;
      nupt[0] = sqrt(thisnupx*thisnupx + thisnupy*thisnupy);
      nueta[0] = asinh(thisnupz/sqrt(thisnupx*thisnupx + thisnupy*thisnupy));
      nuphi[0] = acos(thisnupx/sqrt(thisnupx*thisnupx + thisnupy*thisnupy));
      nue[0] = sqrt(thisnupx*thisnupx + thisnupy*thisnupy + thisnupz*thisnupz);
      cout<<"nupx["<<0<<"] = "<<nupx[0]<<", nupy["<<0<<"] = "<<nupy[0]<<", nupz["<<0<<"] = "<<nupz[0]<<", nupt["<<0<<"] = "<<nupt[0]<<", nueta["<<0<<"] = "<<nueta[0]<<", nuphi["<<0<<"] = "<<nuphi[0]<<", nue["<<0<<"] = "<<nue[0]<<", num["<<0<<"] = "<<num[0]<<endl;
      nupx[1] = thisnubpx;
      nupy[1] = thisnubpy;
      nupz[1] = thisnubpz;
      nupt[1] = sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy);
      nueta[1] = asinh(thisnubpz/sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy));
      nuphi[1] = acos(thisnubpx/sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy));
      nue[1] = sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy + thisnubpz*thisnubpz);
      cout<<"nupx["<<1<<"] = "<<nupx[1]<<", nupy["<<1<<"] = "<<nupy[1]<<", nupz["<<1<<"] = "<<nupz[1]<<", nupt["<<1<<"] = "<<nupt[1]<<", nueta["<<1<<"] = "<<nueta[1]<<", nuphi["<<1<<"] = "<<nuphi[1]<<", nue["<<1<<"] = "<<nue[1]<<", num["<<1<<"] = "<<num[1]<<endl;
    }
    else if(i > 0)
    {
      nupx.push_back(thisnupx);
      nupy.push_back(thisnupy);
      nupz.push_back(thisnupz);
      nupt.push_back(sqrt(thisnupx*thisnupx + thisnupy*thisnupy));
      nueta.push_back(asinh(thisnupz/sqrt(thisnupx*thisnupx + thisnupy*thisnupy)));
      nuphi.push_back(acos(thisnupx/sqrt(thisnupx*thisnupx + thisnupy*thisnupy)));
      nue.push_back(sqrt(thisnupx*thisnupx + thisnupy*thisnupy + thisnupz*thisnupz));
      num.push_back(thisnupx-thisnupx);
      cout<<"nupx["<<2*i<<"] = "<<nupx[2*i]<<", nupy["<<2*i<<"] = "<<nupy[2*i]<<", nupz["<<2*i<<"] = "<<nupz[2*i]<<", nupt["<<2*i<<"] = "<<nupt[2*i]<<", nueta["<<2*i<<"] = "<<nueta[2*i]<<", nuphi["<<2*i<<"] = "<<nuphi[2*i]<<", nue["<<2*i<<"] = "<<nue[2*i]<<", num["<<2*i<<"] = "<<num[2*i]<<endl;
      nupx.push_back(thisnubpx);
      nupy.push_back(thisnubpy);
      nupz.push_back(thisnubpz);
      nupt.push_back(sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy));
      nueta.push_back(asinh(thisnubpz/sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy)));
      nuphi.push_back(acos(thisnubpx/sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy)));
      nue.push_back(sqrt(thisnubpx*thisnubpx + thisnubpy*thisnubpy + thisnubpz*thisnubpz));
      num.push_back(thisnupx-thisnupx);
      cout<<"nupx["<<2*i+1<<"] = "<<nupx[2*i+1]<<", nupy["<<2*i+1<<"] = "<<nupy[2*i+1]<<", nupz["<<2*i+1<<"] = "<<nupz[2*i+1]<<", nupt["<<2*i+1<<"] = "<<nupt[2*i+1]<<", nueta["<<2*i+1<<"] = "<<nueta[2*i+1]<<", nuphi["<<2*i+1<<"] = "<<nuphi[2*i+1]<<", nue["<<2*i+1<<"] = "<<nue[2*i+1]<<", num["<<2*i+1<<"] = "<<num[2*i+1]<<endl;
    }
  }

  return bpx;
}

void quartic(vector<double> polx, vector<double> *nupx, int& cubic_single_root_cmplx)
{
  cout<<"Quartic equation"<<endl;
  vector<double> polxt(polx.size());
  cubic_single_root_cmplx = 0;
  if (polx[4] == 0.)
    cubic(polx, nupx);
  else
  {
    for (size_t i=0; i<polxt.size(); ++i)
    {
      polxt[i]=polx[i]/polx.back(); //normilize to coefficient of highest order (=nupx^4)
      cout << "polxt[" << i << "]=" << polxt[i] << endl;
    }
    if (polxt[0] == 0.)
    {
      (*nupx).push_back(0.);
      cubic(polxt, nupx);
    }
    else
    {
      double e=polxt[2]-3.*polxt[3]*polxt[3]/8.; //k1
      double f=polxt[1]+polxt[3]*polxt[3]*polxt[3]/8.-polxt[2]*polxt[3]/2.; //k2
      double g=polxt[0]-3.*polxt[3]*polxt[3]*polxt[3]*polxt[3]/256.+polxt[3]*polxt[3]*polxt[2]/16.-polxt[3]*polxt[1]/4.; //k3
      cout << "e = " << e << ", f = " << f << ", g = " << g << endl;
      if (g == 0.)
      {
        cout<<"g = 0"<<endl;
        (*nupx).push_back(-polxt[3]/4.);
        vector<double> polxt2(4);
        polxt2[0]=f;
        polxt2[1]=e;
        polxt2[2]=0.;
        polxt2[3]=1.;
        cubic(polxt2, nupx);
        for (size_t i=1; i<(*nupx).size(); ++i)
        {
          (*nupx)[i]-=polxt[3]/4.;
        }
      }
      else if (f == 0.)
      {
        cout<<"f = 0"<<endl;
        vector<double> polxt2(3);
        polxt2[0]=g;
        polxt2[1]=e;
        polxt2[2]=1.;
        vector<double> polxt3;
        quadratic(polxt2, &polxt3);
        for (size_t i=0; i<polxt3.size(); ++i)
        {
          if (polxt3[i]>=0.)
          {
            (*nupx).push_back(sqrt(polxt3[i]-polxt[3]/4.));
            (*nupx).push_back(-sqrt(polxt3[i]-polxt[3]/4.));
          }
        }
      }
      else
      { //d,f,g != 0, default case
        cout<<"g =/= f =/= 0"<<endl;
        vector<double> polxt2(4);
        polxt2[0]=-f*f;
        polxt2[1]=(e*e-4.*g);
        polxt2[2]=2.*e;
        polxt2[3]=1.;
        vector<double> polxt3;
        cubic(polxt2, &polxt3);
        if (polxt3.size()==1 && polxt3[0]<0.)
        {
          cout <<"Warning: quartic: single cubic h^2 solution is negative! (unique sol=" << polxt3[0] << ") => complex root h" << endl;
          cubic_single_root_cmplx++;
          return;
        }
        cout<<"cubic_single_root_cmplx = "<<cubic_single_root_cmplx<<endl;
        double h;
        cout<<"size of polxt3 = "<<polxt3.size()<<endl;
        for (size_t i=0; i<polxt3.size(); ++i)
        {
          cout << "polxt3[" << i << "]=" << polxt3[i] << endl;
          if (polxt3[i]>0.) h=sqrt(polxt3[i]); //t1
        }
        double j=(e+h*h-f/h)/2.; //t2
        cout<<"h = "<<h<<", j = "<<j<<endl;
        vector<double> polxt4(3);
        polxt4[0]=j;
        polxt4[1]=h;
        polxt4[2]=1.; //first polynomial of second degree
        cout<<"polxt4[0] = "<<polxt4[0]<<", polxt4[1] = "<<polxt4[1]<<", polxt4[2] = "<<polxt4[2]<<endl;
        vector<double> polxt5;
        quadratic(polxt4,&polxt5);
        cout<<"size of polxt5 = "<<polxt5.size()<<endl;
        for (size_t i=0; i<polxt5.size(); ++i)
        {
          (*nupx).push_back(polxt5[i]-polxt[3]/4.);
        }
        polxt4[0]=g/j;
        polxt4[1]=-h;
        polxt4[2]=1.; //second polynomial of second degree
        cout<<"polxt4[0] = "<<polxt4[0]<<", polxt4[1] = "<<polxt4[1]<<", polxt4[2] = "<<polxt4[2]<<endl;
        vector<double> polxt6;
        quadratic(polxt4,&polxt6);
        cout<<"size of polxt6 = "<<polxt6.size()<<endl;
        for (size_t i=0; i<polxt6.size(); ++i)
        {
          (*nupx).push_back(polxt6[i]-polxt[3]/4.);
        }
      }
    }
  }
  return;
}

void cubic(vector<double> polx, vector<double> *nupx)
{
  cout<<"Cubic equation"<<endl;
  const double PI = acos(-1.0);
  if (polx[3] == 0.)
    quadratic(polx, nupx);
  else
  {
    cout<<"polx[3] =/= 0"<<endl;
    double q=(polx[2]*polx[2]-3.*polx[1])/9.;
    double r=(2.*polx[2]*polx[2]*polx[2]-9.*polx[1]*polx[2]+27.*polx[0])/54.;
    cout<<"q = "<<q<<", r = "<<r<<endl;
    if (q == 0.)
    { //single root of degree three
      (*nupx).push_back(-polx[2]/3.);
    }
    else if (q*q*q > r*r)
    { //3 real roots
      cout<<"3 real roots!"<<endl;
      double theta=acos(r/sqrt(q*q*q));
      cout<<"-2.*sqrt(q)*cos(theta/3.)-polx[2]/3. = "<<-2.*sqrt(q)*cos(theta/3.)-polx[2]/3.<<", -2.*sqrt(q)*cos((theta+2.*PI)/3.)-polx[2]/3. = "<<-2.*sqrt(q)*cos((theta+2.*PI)/3.)-polx[2]/3.<<", -2.*sqrt(q)*cos((theta+4.*PI)/3.)-polx[2]/3. ="<<-2.*sqrt(q)*cos((theta+4.*PI)/3.)-polx[2]/3.<<endl;
      (*nupx).push_back(-2.*sqrt(q)*cos(theta/3.)-polx[2]/3.);
      (*nupx).push_back(-2.*sqrt(q)*cos((theta+2.*PI)/3.)-polx[2]/3.);
      (*nupx).push_back(-2.*sqrt(q)*cos((theta+4.*PI)/3.)-polx[2]/3.);
    }
    else
    { // 1 real root
      cout<<"1 real root!"<<endl;
      //double powthrd = pow(sqrt(r*r-q*q*q)+fabs(r),1./3.);
      double radicant = -r+sqrt(r*r-q*q*q);
      double powthrd = pow(fabs(radicant),1./3.);
      //double a = -sign(r)*powthrd;
      double a = TMath::Sign(1,radicant)*powthrd;
      double b = (a!=0.) ? q/a : 0.;
      cout<<"radicant = "<<radicant<<", powthrd = "<<powthrd<<", a = "<<a<<", b = "<<b<<", a+b-polx[2]/3. = "<<a+b-polx[2]/3.<<endl;
      (*nupx).push_back(a+b-polx[2]/3.);
    }
  }
  return;
}

void quadratic(vector<double> polx, vector<double> *nupx)
{
  cout<<"Quadratic equation"<<endl;
  cout<<"polx[1]*polx[1]-4.*polx[2]*polx[0] = "<<polx[1]*polx[1]-4.*polx[2]*polx[0]<<endl;
  if (polx[2] == 0.)
  { //linear equation
    cout<<"Linear equation"<<endl;
    (*nupx).push_back(-polx[0]/polx[1]); // a3*x+a4=0 -> x=-a4/a3
  }
  else if (polx[1]*polx[1]==4.*polx[2]*polx[0])
  {
    cout<<"Single root"<<endl;
    (*nupx).push_back(-0.5*polx[1]/polx[2]); // single root of degree two
  }
  else if (polx[1]*polx[1]>4.*polx[2]*polx[0])
  {
    cout<<"Two roots"<<endl;
    double q=-0.5*(polx[1]+TMath::Sign(1,polx[1])*sqrt(polx[1]*polx[1]-4.*polx[2]*polx[0]));
    (*nupx).push_back(q/polx[2]);
    (*nupx).push_back(polx[0]/q);
  }
  return;
}

// int algebraic_pz_1(floats &be, floats &bpx, floats &bpy, floats &bpz, floats &je, floats &jpx, floats &jpy, floats &jpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &mj, floats &ml, double nupx, double nupy, double* nupz, double epsilon)
// {
//   double mblp=sqrt((be[0]+le[0])*(be[0]+le[0])-(bpx[0]+lpx[0])*(bpx[0]+lpx[0])-(bpy[0]+lpy[0])*(bpy[0]+lpy[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));

//   // a1=a11+a12*nu[1]+a13*nu[2];
//   vector<double> a1(3);
//   a1[0]=(80.5*80.5-ml[0]*ml[0])*lpz[0]*0.5 / (le[0]*le[0]-lpz[0]*lpz[0]);
//   a1[1]=lpx[0]*lpz[0] / (le[0]*le[0]-lpz[0]*lpz[0]);
//   a1[2]=lpy[0]*lpz[0] / (le[0]*le[0]-lpz[0]*lpz[0]);

//   //a2=a21+a22*nu[1]+a23*nu[2]+a24*nu[1]^2+a25**nu[1]*nu[2]+a26*nu[2]^2
//   vector<double> a2(6);
//   a2[0]=(80.5*80.5*80.5*80.5+ml[0]*ml[0]*ml[0]*ml[0]-2.*80.5*80.5*ml[0]*ml[0])*0.25 / (le[0]*le[0]-lpz[0]*lpz[0]);
//   a2[3]=-(le[0]*le[0]-lpx[0]*lpx[0]) / (le[0]*le[0]-lpz[0]*lpz[0]);
//   a2[5]=-(le[0]*le[0]-lpy[0]*lpy[0]) / (le[0]*le[0]-lpz[0]*lpz[0]);
//   a2[1]=(80.5*80.5-ml[0]*ml[0])*lpx[0] / (le[0]*le[0]-lpz[0]*lpz[0]);
//   a2[2]=(80.5*80.5-ml[0]*ml[0])*lpy[0] / (le[0]*le[0]-lpz[0]*lpz[0]);
//   a2[4]=2.*lpx[0]*lpy[0] / (le[0]*le[0]-lpz[0]*lpz[0]);

//   // b1=b11+b12*nu[1]+b13*nu[2];
//   vector<double> b1(3);
//   b1[0]=(172.69*172.69-mblp*mblp)*(bpz[0]+lpz[0])*0.5 / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
//   b1[1]=(bpx[0]+lpx[0])*(bpz[0]+lpz[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
//   b1[2]=(bpy[0]+lpy[0])*(bpz[0]+lpz[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));

//   //b2=b21+b22*nu[1]+b23*nu[2]+b24*nu[1]^2+b25*nu[2]*nu[1]+b26*nu[2]^2
//   vector<double> b2(6);
//   b2[0]=(172.69*172.69*172.69*172.69+mblp*mblp*mblp*mblp-2.*172.69*172.69*mblp*mblp)*0.25 / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
//   b2[3]=-((be[0]+le[0])*(be[0]+le[0])-(bpx[0]+lpx[0])*(bpx[0]+lpx[0])) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
//   b2[5]=-((be[0]+le[0])*(be[0]+le[0])-(bpy[0]+lpy[0])*(bpy[0]+lpy[0])) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
//   b2[1]=(172.69*172.69-mblp*mblp)*(bpx[0]+lpx[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
//   b2[2]=(172.69*172.69-mblp*mblp)*(bpy[0]+lpy[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
//   b2[4]=2.*(bpx[0]+lpx[0])*(bpy[0]+lpy[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));

//   //determine first temporary nupz, nubpz
//   double a1val = evalterm1(&a1, nupx, nupy);
//   double a2val = evalterm2(&a2, nupx, nupy);
//   double b1val = evalterm1(&b1, nupx, nupy);
//   double b2val = evalterm2(&b2, nupx, nupy);

//   vector<double> nupz_a, nupz_b;
//   double radicant=a1val*a1val+a2val;
//   if (radicant>=0.)
//   {
//     nupz_a.push_back(a1val+sqrt(radicant));
//     nupz_a.push_back(a1val-sqrt(radicant));
//   }
//   else if (fabs(radicant)<epsilon)
//   {
//     nupz_a.push_back(a1val);
//   }
//   else
//   { //negative radicant => no solution, should never happen
//          cout << "ttdilepsolve: algebraic_pz: nupz_a radicant=" << radicant
//           << " should never happen!" << endl;
//   }
//   radicant=b1val*b1val+b2val;
//   if (radicant>=0.)
//   {
//     nupz_b.push_back(b1val+sqrt(radicant));
//     nupz_b.push_back(b1val-sqrt(radicant));
//     //cout << "nupz_b solutions= " << b1val+sqrt(radicant) << " and " << b1val-sqrt(radicant) << endl;
//   }
//   else if (fabs(radicant)<epsilon)
//   {
//       nupz_b.push_back(b1val);
//   }
//   else
//   { //negative radicant => no solution, should never happen
//          cout << "ttdilepsolve: algebraic_pz: nupz_b radicant=" << radicant
//           << " should never happen!" << endl;
//   }

//   if (nupz_a.size()==0 || nupz_b.size()==0) return -1; //error

//   double nupzchi;
//   double nupzchimin=fabs(nupz_a[0]-nupz_b[0]);
//   int a_min_ind=0, b_min_ind=0;
//   for (size_t j=0; j<nupz_a.size(); ++j)
//   {
//     for (size_t k=0; k<nupz_b.size(); ++k)
//     {
//       nupzchi = fabs(nupz_a[j]-nupz_b[k]);
//       if (nupzchi < nupzchimin)
//       {
//         nupzchimin = nupzchi;
//         a_min_ind=j;
//         b_min_ind=k;
//       }
//     }
//   }

//   if (nupzchimin<sqrt(epsilon))
//   {
//     *nupz = 0.5*(nupz_a[a_min_ind]+nupz_b[b_min_ind]);
//   }
//   else return -1; //error

//   return 0; //success
// }

// int algebraic_pz_2(floats &be, floats &bpx, floats &bpy, floats &bpz, floats &je, floats &jpx, floats &jpy, floats &jpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &mj, floats &ml, double nupx, double nupy, double* nupz, double epsilon, float wmass)
// {
//   double mblp=sqrt((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpx[1]+jpx[2]+lpx[1])*(jpx[1]+jpx[2]+lpx[1])-(jpy[1]+jpy[2]+lpy[1])*(jpy[1]+jpy[2]+lpy[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));

//   // a1=a11+a12*nu[1]+a13*nu[2];
//   vector<double> a1(3);
//   a1[1]=(wmass*wmass-ml[1]*ml[1])*lpz[1]*0.5 / (le[1]*le[1]-lpz[1]*lpz[1]);
//   a1[1]=lpx[1]*lpz[1] / (le[1]*le[1]-lpz[1]*lpz[1]);
//   a1[2]=lpy[1]*lpz[1] / (le[1]*le[1]-lpz[1]*lpz[1]);

//   //a2=a21+a22*nu[1]+a23*nu[2]+a24*nu[1]^2+a25**nu[1]*nu[2]+a26*nu[2]^2
//   vector<double> a2(6);
//   a2[1]=(wmass*wmass*wmass*wmass+ml[1]*ml[1]*ml[1]*ml[1]-2.*wmass*wmass*ml[1]*ml[1])*0.25 / (le[1]*le[1]-lpz[1]*lpz[1]);
//   a2[3]=-(le[1]*le[1]-lpx[1]*lpx[1]) / (le[1]*le[1]-lpz[1]*lpz[1]);
//   a2[5]=-(le[1]*le[1]-lpy[1]*lpy[1]) / (le[1]*le[1]-lpz[1]*lpz[1]);
//   a2[1]=(wmass*wmass-ml[1]*ml[1])*lpx[1] / (le[1]*le[1]-lpz[1]*lpz[1]);
//   a2[2]=(wmass*wmass-ml[1]*ml[1])*lpy[1] / (le[1]*le[1]-lpz[1]*lpz[1]);
//   a2[4]=2.*lpx[1]*lpy[1] / (le[1]*le[1]-lpz[1]*lpz[1]);

//   // b1=b11+b12*nu[1]+b13*nu[2];
//   vector<double> b1(3);
//   b1[1]=(172.69*172.69-mblp*mblp)*(jpz[1]+jpz[2]+lpz[1])*0.5 / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));
//   b1[1]=(bpx[1]+lpx[1])*(bpz[1]+lpz[1]) / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));
//   b1[2]=(bpy[1]+lpy[1])*(bpz[1]+lpz[1]) / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));

//   //b2=b21+b22*nu[1]+b23*nu[2]+b24*nu[1]^2+b25*nu[2]*nu[1]+b26*nu[2]^2
//   vector<double> b2(6);
//   b2[1]=(172.69*172.69*172.69*172.69+mblp*mblp*mblp*mblp-2.*172.69*172.69*mblp*mblp)*0.25 / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));
//   b2[3]=-((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpx[1]+jpx[2]+lpx[1])*(jpx[1]+jpx[2]+lpx[1])) / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));
//   b2[5]=-((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpy[1]+jpy[2]+lpy[1])*(jpy[1]+jpy[2]+lpy[1])) / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));
//   b2[1]=(172.69*172.69-mblp*mblp)*(jpx[1]+jpx[2]+lpx[1]) / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));
//   b2[2]=(172.69*172.69-mblp*mblp)*(jpy[1]+jpy[2]+lpy[1]) / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));
//   b2[4]=2.*(jpx[1]+jpx[2]+lpx[1])*(jpy[1]+jpy[2]+lpy[1]) / ((je[1]+je[2]+le[1])*(je[1]+je[2]+le[1])-(jpz[1]+jpz[2]+lpz[1])*(jpz[1]+jpz[2]+lpz[1]));

//   //determine first temporary nupz, nubpz
//   double a1val = evalterm1(&a1, nupx, nupy);
//   double a2val = evalterm2(&a2, nupx, nupy);
//   double b1val = evalterm1(&b1, nupx, nupy);
//   double b2val = evalterm2(&b2, nupx, nupy);

//   vector<double> nupz_a, nupz_b;
//   double radicant=a1val*a1val+a2val;
//   if (radicant>=0.)
//   {
//     nupz_a.push_back(a1val+sqrt(radicant));
//     nupz_a.push_back(a1val-sqrt(radicant));
//   }
//   else if (fabs(radicant)<epsilon)
//   {
//     nupz_a.push_back(a1val);
//   }
//   else
//   { //negative radicant => no solution, should never happen
//          cout << "ttdilepsolve: algebraic_pz: nupz_a radicant=" << radicant
//           << " should never happen!" << endl;
//   }
//   radicant=b1val*b1val+b2val;
//   if (radicant>=0.)
//   {
//     nupz_b.push_back(b1val+sqrt(radicant));
//     nupz_b.push_back(b1val-sqrt(radicant));
//     //cout << "nupz_b solutions= " << b1val+sqrt(radicant) << " and " << b1val-sqrt(radicant) << endl;
//   }
//   else if (fabs(radicant)<epsilon)
//   {
//       nupz_b.push_back(b1val);
//   }
//   else
//   { //negative radicant => no solution, should never happen
//          cout << "ttdilepsolve: algebraic_pz: nupz_b radicant=" << radicant
//           << " should never happen!" << endl;
//   }

//   if (nupz_a.size()==0 || nupz_b.size()==0) return -1; //error

//   double nupzchi;
//   double nupzchimin=fabs(nupz_a[0]-nupz_b[0]);
//   int a_min_ind=0, b_min_ind=0;
//   for (size_t j=0; j<nupz_a.size(); ++j)
//   {
//     for (size_t k=0; k<nupz_b.size(); ++k)
//     {
//       nupzchi = fabs(nupz_a[j]-nupz_b[k]);
//       if (nupzchi < nupzchimin)
//       {
//         nupzchimin = nupzchi;
//         a_min_ind=j;
//         b_min_ind=k;
//       }
//     }
//   }

//   if (nupzchimin<sqrt(epsilon))
//   {
//     *nupz = 0.5*(nupz_a[a_min_ind]+nupz_b[b_min_ind]);
//   }
//   else return -1; //error

//   return 0; //success
// }

int algebraic_pz_1(floats &be, floats &bpx, floats &bpy, floats &bpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &ml, double nupx, double nupy, double* nupz, double epsilon, float wmass)
{
  double mblp=sqrt((be[0]+le[0])*(be[0]+le[0])-(bpx[0]+lpx[0])*(bpx[0]+lpx[0])-(bpy[0]+lpy[0])*(bpy[0]+lpy[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));

  // a1=a11+a12*nu[1]+a13*nu[2];
  vector<double> a1(3);
  a1[0]=(wmass*wmass-ml[0]*ml[0])*lpz[0]*0.5 / (le[0]*le[0]-lpz[0]*lpz[0]);
  a1[1]=lpx[0]*lpz[0] / (le[0]*le[0]-lpz[0]*lpz[0]);
  a1[2]=lpy[0]*lpz[0] / (le[0]*le[0]-lpz[0]*lpz[0]);

  //a2=a21+a22*nu[1]+a23*nu[2]+a24*nu[1]^2+a25**nu[1]*nu[2]+a26*nu[2]^2
  vector<double> a2(6);
  a2[0]=(wmass*wmass*wmass*wmass+ml[0]*ml[0]*ml[0]*ml[0]-2.*wmass*wmass*ml[0]*ml[0])*0.25 / (le[0]*le[0]-lpz[0]*lpz[0]);
  a2[3]=-(le[0]*le[0]-lpx[0]*lpx[0]) / (le[0]*le[0]-lpz[0]*lpz[0]);
  a2[5]=-(le[0]*le[0]-lpy[0]*lpy[0]) / (le[0]*le[0]-lpz[0]*lpz[0]);
  a2[1]=(wmass*wmass-ml[0]*ml[0])*lpx[0] / (le[0]*le[0]-lpz[0]*lpz[0]);
  a2[2]=(wmass*wmass-ml[0]*ml[0])*lpy[0] / (le[0]*le[0]-lpz[0]*lpz[0]);
  a2[4]=2.*lpx[0]*lpy[0] / (le[0]*le[0]-lpz[0]*lpz[0]);

  // b1=b11+b12*nu[1]+b13*nu[2];
  vector<double> b1(3);
  b1[0]=(172.69*172.69-mblp*mblp)*(bpz[0]+lpz[0])*0.5 / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
  b1[1]=(bpx[0]+lpx[0])*(bpz[0]+lpz[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
  b1[2]=(bpy[0]+lpy[0])*(bpz[0]+lpz[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));

  //b2=b21+b22*nu[1]+b23*nu[2]+b24*nu[1]^2+b25*nu[2]*nu[1]+b26*nu[2]^2
  vector<double> b2(6);
  b2[0]=(172.69*172.69*172.69*172.69+mblp*mblp*mblp*mblp-2.*172.69*172.69*mblp*mblp)*0.25 / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
  b2[3]=-((be[0]+le[0])*(be[0]+le[0])-(bpx[0]+lpx[0])*(bpx[0]+lpx[0])) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
  b2[5]=-((be[0]+le[0])*(be[0]+le[0])-(bpy[0]+lpy[0])*(bpy[0]+lpy[0])) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
  b2[1]=(172.69*172.69-mblp*mblp)*(bpx[0]+lpx[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
  b2[2]=(172.69*172.69-mblp*mblp)*(bpy[0]+lpy[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));
  b2[4]=2.*(bpx[0]+lpx[0])*(bpy[0]+lpy[0]) / ((be[0]+le[0])*(be[0]+le[0])-(bpz[0]+lpz[0])*(bpz[0]+lpz[0]));

  //determine first temporary nupz, nubpz
  double a1val = evalterm1(&a1, nupx, nupy);
  double a2val = evalterm2(&a2, nupx, nupy);
  double b1val = evalterm1(&b1, nupx, nupy);
  double b2val = evalterm2(&b2, nupx, nupy);

  vector<double> nupz_a, nupz_b;
  double radicant=a1val*a1val+a2val;
  if (radicant>=0.)
  {
    nupz_a.push_back(a1val+sqrt(radicant));
    nupz_a.push_back(a1val-sqrt(radicant));
  }
  else if (fabs(radicant)<epsilon)
  {
    nupz_a.push_back(a1val);
  }
  else
  { //negative radicant => no solution, should never happen
         cout << "ttdilepsolve: algebraic_pz: nupz_a radicant=" << radicant
          << " should never happen!" << endl;
  }
  radicant=b1val*b1val+b2val;
  if (radicant>=0.)
  {
    nupz_b.push_back(b1val+sqrt(radicant));
    nupz_b.push_back(b1val-sqrt(radicant));
    //cout << "nupz_b solutions= " << b1val+sqrt(radicant) << " and " << b1val-sqrt(radicant) << endl;
  }
  else if (fabs(radicant)<epsilon)
  {
      nupz_b.push_back(b1val);
  }
  else
  { //negative radicant => no solution, should never happen
         cout << "ttdilepsolve: algebraic_pz: nupz_b radicant=" << radicant
          << " should never happen!" << endl;
  }

  if (nupz_a.size()==0 || nupz_b.size()==0) return -1; //error

  double nupzchi;
  double nupzchimin=fabs(nupz_a[0]-nupz_b[0]);
  int a_min_ind=0, b_min_ind=0;
  for (size_t j=0; j<nupz_a.size(); ++j)
  {
    for (size_t k=0; k<nupz_b.size(); ++k)
    {
      nupzchi = fabs(nupz_a[j]-nupz_b[k]);
      if (nupzchi < nupzchimin)
      {
        nupzchimin = nupzchi;
        a_min_ind=j;
        b_min_ind=k;
      }
    }
  }

  if (nupzchimin<sqrt(epsilon))
  {
    *nupz = 0.5*(nupz_a[a_min_ind]+nupz_b[b_min_ind]);
  }
  else return -1; //error

  return 0; //success
}

int algebraic_pz_2(floats &be, floats &bpx, floats &bpy, floats &bpz, floats &le, floats &lpx, floats &lpy, floats &lpz, floats &mb, floats &ml, double nupx, double nupy, double* nupz, double epsilon, float wmass)
{
  double mblp=sqrt((be[1]+le[1])*(be[1]+le[1])-(bpx[1]+lpx[1])*(bpx[1]+lpx[1])-(bpy[1]+lpy[1])*(bpy[1]+lpy[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));

  // a1=a11+a12*nu[1]+a13*nu[2];
  vector<double> a1(3);
  a1[0]=(wmass*wmass-ml[1]*ml[1])*lpz[1]*0.5 / (le[1]*le[1]-lpz[1]*lpz[1]);
  a1[1]=lpx[1]*lpz[1] / (le[1]*le[1]-lpz[1]*lpz[1]);
  a1[2]=lpy[1]*lpz[1] / (le[1]*le[1]-lpz[1]*lpz[1]);

  //a2=a21+a22*nu[1]+a23*nu[2]+a24*nu[1]^2+a25**nu[1]*nu[2]+a26*nu[2]^2
  vector<double> a2(6);
  a2[0]=(wmass*wmass*wmass*wmass+ml[1]*ml[1]*ml[1]*ml[1]-2.*wmass*wmass*ml[1]*ml[1])*0.25 / (le[1]*le[1]-lpz[1]*lpz[1]);
  a2[3]=-(le[1]*le[1]-lpx[1]*lpx[1]) / (le[1]*le[1]-lpz[1]*lpz[1]);
  a2[5]=-(le[1]*le[1]-lpy[1]*lpy[1]) / (le[1]*le[1]-lpz[1]*lpz[1]);
  a2[1]=(wmass*wmass-ml[1]*ml[1])*lpx[1] / (le[1]*le[1]-lpz[1]*lpz[1]);
  a2[2]=(wmass*wmass-ml[1]*ml[1])*lpy[1] / (le[1]*le[1]-lpz[1]*lpz[1]);
  a2[4]=2.*lpx[1]*lpy[1] / (le[1]*le[1]-lpz[1]*lpz[1]);

  // b1=b11+b12*nu[1]+b13*nu[2];
  vector<double> b1(3);
  b1[0]=(172.69*172.69-mblp*mblp)*(bpz[1]+lpz[1])*0.5 / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));
  b1[1]=(bpx[1]+lpx[1])*(bpz[1]+lpz[1]) / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));
  b1[2]=(bpy[1]+lpy[1])*(bpz[1]+lpz[1]) / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));

  //b2=b21+b22*nu[1]+b23*nu[2]+b24*nu[1]^2+b25*nu[2]*nu[1]+b26*nu[2]^2
  vector<double> b2(6);
  b2[0]=(172.69*172.69*172.69*172.69+mblp*mblp*mblp*mblp-2.*172.69*172.69*mblp*mblp)*0.25 / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));
  b2[3]=-((be[1]+le[1])*(be[1]+le[1])-(bpx[1]+lpx[1])*(bpx[1]+lpx[1])) / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));
  b2[5]=-((be[1]+le[1])*(be[1]+le[1])-(bpy[1]+lpy[1])*(bpy[1]+lpy[1])) / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));
  b2[1]=(172.69*172.69-mblp*mblp)*(bpx[1]+lpx[1]) / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));
  b2[2]=(172.69*172.69-mblp*mblp)*(bpy[1]+lpy[1]) / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));
  b2[4]=2.*(bpx[1]+lpx[1])*(bpy[1]+lpy[1]) / ((be[1]+le[1])*(be[1]+le[1])-(bpz[1]+lpz[1])*(bpz[1]+lpz[1]));

  //determine first temporary nupz, nubpz
  double a1val = evalterm1(&a1, nupx, nupy);
  double a2val = evalterm2(&a2, nupx, nupy);
  double b1val = evalterm1(&b1, nupx, nupy);
  double b2val = evalterm2(&b2, nupx, nupy);

  vector<double> nupz_a, nupz_b;
  double radicant=a1val*a1val+a2val;
  if (radicant>=0.)
  {
    nupz_a.push_back(a1val+sqrt(radicant));
    nupz_a.push_back(a1val-sqrt(radicant));
  }
  else if (fabs(radicant)<epsilon)
  {
    nupz_a.push_back(a1val);
  }
  else
  { //negative radicant => no solution, should never happen
         cout << "ttdilepsolve: algebraic_pz: nupz_a radicant=" << radicant
          << " should never happen!" << endl;
  }
  radicant=b1val*b1val+b2val;
  if (radicant>=0.)
  {
    nupz_b.push_back(b1val+sqrt(radicant));
    nupz_b.push_back(b1val-sqrt(radicant));
    //cout << "nupz_b solutions= " << b1val+sqrt(radicant) << " and " << b1val-sqrt(radicant) << endl;
  }
  else if (fabs(radicant)<epsilon)
  {
      nupz_b.push_back(b1val);
  }
  else
  { //negative radicant => no solution, should never happen
         cout << "ttdilepsolve: algebraic_pz: nupz_b radicant=" << radicant
          << " should never happen!" << endl;
  }

  if (nupz_a.size()==0 || nupz_b.size()==0) return -1; //error

  double nupzchi;
  double nupzchimin=fabs(nupz_a[0]-nupz_b[0]);
  int a_min_ind=0, b_min_ind=0;
  for (size_t j=0; j<nupz_a.size(); ++j)
  {
    for (size_t k=0; k<nupz_b.size(); ++k)
    {
      nupzchi = fabs(nupz_a[j]-nupz_b[k]);
      if (nupzchi < nupzchimin)
      {
        nupzchimin = nupzchi;
        a_min_ind=j;
        b_min_ind=k;
      }
    }
  }

  if (nupzchimin<sqrt(epsilon))
  {
    *nupz = 0.5*(nupz_a[a_min_ind]+nupz_b[b_min_ind]);
  }
  else return -1; //error

  return 0; //success
}

double evalterm1(vector<double> *a1, double nupx, double nupy)
{
  return (*a1)[0]+(*a1)[1]*nupx+(*a1)[2]*nupy;
}

double evalterm2(vector<double> *a2, double nupx, double nupy)
{
  return (*a2)[0]+(*a2)[1]*nupx+(*a2)[2]*nupy+(*a2)[3]*nupx*nupx+(*a2)[4]*nupx*nupy+(*a2)[5]*nupy*nupy;
}

RVec< RVec<int> > distinct_comb_1_even(int nn)
{
  // cout<<"L2L2L2"<<endl;
  // cout<<"nn = "<<nn<<endl;
  RVec<int> idx1;
  RVec< RVec<int> > res;
  idx1.reserve(nn*(nn-1)/4);

  if(nn < 2)
  {
    idx1.emplace_back(0);
    return res;
  }  

  for(unsigned int i = 0; i < nn; i = i+2)
  {
	    idx1.emplace_back(i);
  }      
  // cout<<"idx1 = "<<idx1<<endl;
  res.emplace_back(idx1);
  return res;      
};

RVec< RVec<int> > distinct_comb_1_odd(int nn)
{
  // cout<<"L3L3L3"<<endl;
  // cout<<"nn = "<<nn<<endl;
  RVec<int> idx1;
  RVec< RVec<int> > res;
  idx1.reserve(nn*(nn-1)/4);

  if(nn < 2)
  {
    idx1.emplace_back(0);
    return res;
  }  

  for(unsigned int i = 1; i < nn; i = i+2)
  {
	    idx1.emplace_back(i);
  }      
  // cout<<"idx1 = "<<idx1<<endl;
  res.emplace_back(idx1);
  return res;      
};

RVec< RVec<int> > distinct_comb_2_even_odd(int nn)
{
  // cout<<"L4L4L4"<<endl;
  // cout<<"nn = "<<nn<<endl;
  RVec<int> idx1;
  RVec<int> idx2;
  RVec< RVec<int> > res;
  idx1.reserve(nn*(nn-1)/4);
  idx2.reserve(nn*(nn-1)/4);

  if(nn < 2)
  {
    idx1.emplace_back(0);
    idx2.emplace_back(0);
    return res;
  }  

  for(unsigned int i = 0; i < nn; i = i+2)
  {
    for(unsigned int j = 1; j < nn; j = j+2)
    {
	    idx1.emplace_back(i);
	    idx2.emplace_back(j);
    }
  }      
  // cout<<"idx1 = "<<idx1<<", idx2 = "<<idx2<<endl;
  res.emplace_back(idx1);
  res.emplace_back(idx2);
  return res;      
};

FourVectorRVec Leptonneutrino (FourVectorRVec lepton, FourVectorRVec neutrino, int nn)
{
  // cout<<"H1H1H1"<<endl;
  FourVectorRVec vecs;
  FourVectorRVec temp;
  vecs.reserve(nn/2);
  // cout<<"lepton = "<<lepton<<endl;

  for(unsigned int i = 0; i < nn/2; i++)
  {
    // cout<<"neutrino["<<i<<"] = "<<neutrino[i]<<endl;
    if(neutrino[i].Pt() == 0)
    {
      vecs.emplace_back(0, 0, 0, 0);
      return vecs;
    }
    temp = lepton + neutrino[i];
    // cout<<"temp = "<<temp<<", Pt = "<<temp[0].Pt()<<", Eta = "<<temp[0].Eta()<<", Phi = "<<temp[0].Phi()<<", Mass = "<<temp[0].M()<<endl;
    vecs.emplace_back(temp[0].Pt(),temp[0].Eta(),temp[0].Phi(),temp[0].M());
  }
  // cout<<"vecs = "<<vecs<<endl;
  return vecs;
};

floats deltaeta(double &leta, floats &neta)
{
  // cout<<"H2H2H2"<<endl;
  floats deta;
  unsigned int nn = neta.size();
  if(neta[0] == 0)
  {
    deta.push_back(0);
    return deta;
  }
  for(unsigned int i = 0; i < nn; i++)
  {
    double delta = abs(float(leta) - neta[i]);
    deta.push_back(delta);
  }
  // cout<<"deta = "<<deta<<endl;
  return deta;
};

floats deltaphi(double &lphi, floats &nphi)
{
  // cout<<"H3H3H3"<<endl;
  floats dphi;
  unsigned int nn = nphi.size();
  if(nphi[0] == 0)
  {
    dphi.push_back(0);
    return dphi;
  }
  for(unsigned int i = 0; i < nn; i++)
  {
    double delta = ROOT::VecOps::DeltaPhi(float(lphi), nphi[i]);
    dphi.push_back(delta);
  }
  // cout<<"dphi = "<<dphi<<endl;
  return dphi;
};

floats deltaR(double &leta, floats &neta, double &lphi, floats &nphi)
{
  // cout<<"H4H4H4"<<endl;
  floats dR;
  unsigned int nn = neta.size();
  if(neta[0] == 0)
  {
    dR.push_back(0);
    return dR;
  }
  for(unsigned int i = 0; i < nn; i++)
  {
    double delta = ROOT::VecOps::DeltaR(float(leta), neta[i], float(lphi), nphi[i]);
    dR.push_back(delta);
  }
  // cout<<"dR = "<<dR<<endl;
  return dR;
};

floats transverse_mass(double &lmass, floats &te, double &lpt, floats &npt, floats &dphi)
{
  // cout<<"H5H5H5"<<endl;
  floats tm;
  unsigned int nn = npt.size();
  if(npt[0] == 0)
  {
    tm.push_back(0);
    return tm;
  }
  for(unsigned int i = 0; i < nn; i++)
  {
    double transmass = sqrt(float(lmass)*float(lmass) + 2*te[i]*te[i] - 2*float(lpt)*npt[i]*cos(dphi[i]));
    tm.push_back(transmass);
  }
  // cout<<"tm = "<<tm<<endl;
  return tm;
};

FourVectorRVec Leptonneutrinobjet (FourVectorRVec bjet, FourVectorRVec leptonneutrino, int nn)
{
  // cout<<"H6H6H6"<<endl;
  FourVectorRVec vecs;
  FourVectorRVec temp;
  vecs.reserve(nn/2);
  // cout<<"bjet = "<<bjet<<endl;

  for(unsigned int i = 0; i < nn/2; i++)
  {
    // cout<<"leptonneutrino["<<i<<"] = "<<leptonneutrino[i]<<endl;
    if(leptonneutrino[i].Pt() == 0)
    {
      vecs.emplace_back(0, 0, 0, 0);
      return vecs;
    }
    temp = bjet + leptonneutrino[i];
    // cout<<"temp = "<<temp<<", Pt = "<<temp[0].Pt()<<", Eta = "<<temp[0].Eta()<<", Phi = "<<temp[0].Phi()<<", Mass = "<<temp[0].M()<<endl;
    vecs.emplace_back(temp[0].Pt(),temp[0].Eta(),temp[0].Phi(),temp[0].M());
  }
  // cout<<"vecs = "<<vecs<<endl;
  return vecs;
};

floats transverse_mass_bjet(double &bmass, floats &lnmass, floats &te, double &jpt, floats &lnpt, floats &dphi)
{
  // cout<<"H7H7H7"<<endl;
  floats tm;
  unsigned int nn = lnpt.size();
  if(lnpt[0] == 0)
  {
    tm.push_back(0);
    return tm;
  }
  for(unsigned int i = 0; i < nn; i++)
  {
    double transmass = sqrt(float(bmass)*float(bmass) + lnmass[i]*lnmass[i] + 2*te[i]*te[i] - 2*float(jpt)*lnpt[i]*cos(dphi[i]));
    tm.push_back(transmass);
  }
  // cout<<"tm = "<<tm<<endl;
  return tm;
};