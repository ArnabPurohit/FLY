/*
 * TprimeAnalyser.cpp
 *
 *  Created on: May 6, 2022
 *      Author: suyong
 *      Developper: cdozen
 */

#include "TprimeAnalyser.h"
#include "utility.h"
#include "Math/GenVector/VectorUtil.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/TProcessExecutor.hxx>
#include <TStopwatch.h>
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include "correction.h"
#include "TCanvas.h"
// #include "Mt2/Basic_Mt2_332_Calculator.h"
// #include "Mt2/Mt2Units.h"
// #include "Minuit2/MnUserParameterState.h"
// #include "Minuit2/MnMigrad.h"
// #include "Minuit2/MnSimplex.h"
// #include "Minuit2/MnScan.h"
// #include "Minuit2/FunctionMinimum.h"
#include <fstream>
using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;
using floats =  ROOT::VecOps::RVec<float>;
using ints =  ROOT::VecOps::RVec<int>;
using FourVector = ROOT::Math::PtEtaPhiMVector;
using FourVectorRVec = ROOT::VecOps::RVec<FourVector>;
using FourVectorVec = std::vector<FourVector>;
using correction::CorrectionSet;
// using ROOT::Minuit2::MnUserParameters;
// using ROOT::Minuit2::MnMigrad;
// using ROOT::Minuit2::MnSimplex;
// using ROOT::Minuit2::MnScan;
// using ROOT::Minuit2::FunctionMinimum;

TprimeAnalyser::TprimeAnalyser(TTree *t, std::string outfilename, std::string cuts, std::string Region, std::string year)
:NanoAODAnalyzerrdframe(t, outfilename)
{
    //initiliaze the HLT names in your analyzer class
    HLT2018Names= {"HLT_IsoMu24","HLT_Ele32_WPTight_Gsf","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"};
	HLT2017Names= {"HLT_IsoMu27","HLT_Ele35_WPTight_Gsf","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"};
	HLT2016Names= {"Name1","Name2"};
}

// Define your cuts here
void TprimeAnalyser::defineCuts(std::string cuts)
{
	if (debug)
	{
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
	}

	auto Nentry = _rlm.Count();
	//This is how you can express a range of the first 100 entries
	// _rlm = _rlm.Range(0, 10000);
	//auto Nentry_100 = _rlm.Count();
	std::cout<< "-------------------------------------------------------------------" << std::endl;
	cout << "Usage of ranges:\n"
	     << " - All entries: " << *Nentry << endl;
		//<< " - Entries from 0 to 100: " << *Nentry_100 << endl;
	std::cout<< "-------------------------------------------------------------------" << std::endl;

    // std::cout<<"cuts = "<<cuts<<std::endl;
    std::string delimiter = "\n";
    std::string Cut0_bis = cuts.substr(0, cuts.find(delimiter));
    // std::cout<<"Cut0_bis = "<<Cut0_bis<<std::endl;
    std::string Cut0 = Form("%s && (%s)",Cut0_bis.c_str(),setHLT().c_str());
    // std::cout<<"Cut0 = "<<Cut0<<std::endl;
    addCuts(Cut0.c_str(),"0");
    cuts.erase(0, cuts.find(delimiter) + delimiter.length());
    std::string Cut00 = cuts.substr(0, cuts.find(delimiter));
    // std::cout<<"Cut00 = "<<Cut00<<std::endl;
    addCuts(Cut00.c_str(),"00");
    cuts.erase(0, cuts.find(delimiter) + delimiter.length());
    std::string Cut000 = cuts.substr(0, cuts.find(delimiter));
    // std::cout<<"Cut000 = "<<Cut000<<std::endl;
    addCuts(Cut000.c_str(),"000");
    cuts.erase(0, cuts.find(delimiter) + delimiter.length());
    std::string Cut0000 = cuts.substr(0, cuts.find(delimiter));
    // std::cout<<"Cut0000 = "<<Cut0000<<std::endl;
    addCuts(Cut0000.c_str(),"0000");
    cuts.erase(0, cuts.find(delimiter) + delimiter.length());
    std::string Cut00000 = cuts.substr(0, cuts.find(delimiter));
    // std::cout<<"Cut00000 = "<<Cut00000<<std::endl;
    addCuts(Cut00000.c_str(),"00000");
    // cuts.erase(0, cuts.find(delimiter) + delimiter.length());
    // std::string Cut000000 = cuts.substr(0, cuts.find(delimiter));
    // // std::cout<<"Cut000000 = "<<Cut000000<<std::endl;
    // addCuts(Cut000000.c_str(),"000000");
    
}

//===============================Find Good Electrons===========================================//
//: Define Good Electrons in rdata frame
//=============================================================================================//

void TprimeAnalyser::selectElectrons(std::string Region)
{
    cout << "select good electrons" << endl;
    if (debug)
    {
    std::cout<< "================================//=================================" << std::endl;
    std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
    std::cout<< "================================//=================================" << std::endl;
    }

    _rlm = _rlm.Define("goodElecIDlooseTTSL", ElectronID(2));
    _rlm = _rlm.Define("goodElectronslooseTTSL","goodElecIDlooseTTSL && Electron_pt > 25 && (abs(Electron_eta) < 1.442 || (abs(Electron_eta) > 1.566 && abs(Electron_eta) < 2.5)) && Electron_miniPFRelIso_all < 0.10 && Electron_sip3d < 2 && Electron_tightCharge == 2")
               .Define("Selected_electron_looseTTSL_pt","Electron_pt[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_number","int(Selected_electron_looseTTSL_pt.size())")
               .Define("Selected_electron_looseTTSL_leading_pt","(Selected_electron_looseTTSL_number > 0) ? Selected_electron_looseTTSL_pt[0]:-100")
               .Define("Selected_electron_looseTTSL_subleading_pt","(Selected_electron_looseTTSL_number > 1) ? Selected_electron_looseTTSL_pt[1]:-100")
               .Define("Selected_electron_looseTTSL_sum_two_electrons_pt","(Selected_electron_looseTTSL_number > 1) ? Selected_electron_looseTTSL_pt[0] + Selected_electron_looseTTSL_pt[1]:-100")
               .Define("Selected_electron_looseTTSL_sum_all_electrons_pt","Sum(Selected_electron_looseTTSL_pt)")
               .Define("Selected_electron_looseTTSL_eta","Electron_eta[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_leading_eta","(Selected_electron_looseTTSL_number > 0) ? Selected_electron_looseTTSL_eta[0]:-10")
               .Define("Selected_electron_looseTTSL_subleading_eta","(Selected_electron_looseTTSL_number > 1) ? Selected_electron_looseTTSL_eta[1]:-10")
               .Define("Selected_electron_looseTTSL_phi","Electron_phi[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_mass","Electron_mass[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_transverse_energy","sqrt(Selected_electron_looseTTSL_pt*Selected_electron_looseTTSL_pt + Selected_electron_looseTTSL_mass*Selected_electron_looseTTSL_mass)")
               .Define("Selected_electron_looseTTSL_charge","Electron_charge[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_charge_sum","(Selected_electron_looseTTSL_number > 1) ? Selected_electron_looseTTSL_charge[0] + Selected_electron_looseTTSL_charge[1]:-100")
               .Define("Selected_electron_looseTTSL_deltaeta","(Selected_electron_looseTTSL_number > 1) ? abs(Selected_electron_looseTTSL_eta[1] - Selected_electron_looseTTSL_eta[0]):-10")
               .Define("Selected_electron_looseTTSL_deltaphi","(Selected_electron_looseTTSL_number > 1) ? abs(ROOT::VecOps::DeltaPhi(Selected_electron_looseTTSL_phi[0],Selected_electron_looseTTSL_phi[1])):-10")
               .Define("Selected_electron_looseTTSL_deltaR","(Selected_electron_looseTTSL_number > 1) ? ROOT::VecOps::DeltaR(Selected_electron_looseTTSL_eta[0],Selected_electron_looseTTSL_eta[1],Selected_electron_looseTTSL_phi[0],Selected_electron_looseTTSL_phi[1]):-10")
               .Define("Selected_electron_looseTTSL_miniPFRelIso_all","Electron_miniPFRelIso_all[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_miniPFRelIso_chg","Electron_miniPFRelIso_chg[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_pfRelIso03_all","Electron_pfRelIso03_all[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_pfRelIso03_chg","Electron_pfRelIso03_chg[goodElectronslooseTTSL]");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_electron_looseTTSL_genPartIdx","Electron_genPartIdx[goodElectronslooseTTSL]")
                   .Define("Selected_electron_looseTTSL_genPartpdgId","ROOT::VecOps::Take(GenPart_pdgId,Selected_electron_looseTTSL_genPartIdx)")
                   .Define("Selected_electron_looseTTSL_genPartIdxMother","ROOT::VecOps::Take(GenPart_genPartIdxMother,Selected_electron_looseTTSL_genPartIdx)")
                   .Define("Selected_electron_looseTTSL_genPartpdgIdMother","ROOT::VecOps::Take(GenPart_pdgId,Selected_electron_looseTTSL_genPartIdxMother)")
                   .Define("Selected_electron_looseTTSL_genPartFlav","Electron_genPartFlav[goodElectronslooseTTSL]");
    }
    _rlm = _rlm.Define("Selected_electron_looseTTSL_pdgId","Electron_pdgId[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_ip3d","Electron_ip3d[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_sip3d","Electron_sip3d[goodElectronslooseTTSL]")
               .Define("Selected_electron_looseTTSL_idx", ::good_idx,{"goodElectronslooseTTSL"});
    
    _rlm = _rlm.Define("goodElecID", ElectronID(4));
    _rlm = _rlm.Define("goodElectrons","goodElecID && Electron_pt > 25 && (abs(Electron_eta) < 1.442 || (abs(Electron_eta) > 1.566 && abs(Electron_eta) < 2.5)) && Electron_miniPFRelIso_all < 0.05 && Electron_sip3d < 2 && Electron_tightCharge == 2") // && Electron_tightCharge == 2
               .Define("Selected_electron_pt","Electron_pt[goodElectrons]")
               .Define("Selected_electron_number","int(Selected_electron_pt.size())")
               .Define("Selected_electron_leading_pt","(Selected_electron_number > 0) ? Selected_electron_pt[0]:-100")
               .Define("Selected_electron_subleading_pt","(Selected_electron_number > 1) ? Selected_electron_pt[1]:-100")
               .Define("Selected_electron_sum_two_electrons_pt","(Selected_electron_number > 1) ? Selected_electron_pt[0] + Selected_electron_pt[1]:-100")
               .Define("Selected_electron_sum_all_electrons_pt","Sum(Selected_electron_pt)")
               .Define("Selected_electron_eta","Electron_eta[goodElectrons]")
               .Define("Selected_electron_leading_eta","(Selected_electron_number > 0) ? Selected_electron_eta[0]:-10")
               .Define("Selected_electron_subleading_eta","(Selected_electron_number > 1) ? Selected_electron_eta[1]:-10")
               .Define("Selected_electron_phi","Electron_phi[goodElectrons]")
               .Define("Selected_electron_mass","Electron_mass[goodElectrons]")
               .Define("Selected_electron_transverse_energy","sqrt(Selected_electron_pt*Selected_electron_pt + Selected_electron_mass*Selected_electron_mass)")
               .Define("Selected_electron_charge","Electron_charge[goodElectrons]")
               .Define("Selected_electron_charge_sum","(Selected_electron_number > 1) ? Selected_electron_charge[0] + Selected_electron_charge[1]:-100")
               .Define("Selected_electron_deltaeta","(Selected_electron_number > 1) ? abs(Selected_electron_eta[1] - Selected_electron_eta[0]):-10")
               .Define("Selected_electron_deltaphi","(Selected_electron_number > 1) ? abs(ROOT::VecOps::DeltaPhi(Selected_electron_phi[0],Selected_electron_phi[1])):-10")
               .Define("Selected_electron_deltaR","(Selected_electron_number > 1) ? ROOT::VecOps::DeltaR(Selected_electron_eta[0],Selected_electron_eta[1],Selected_electron_phi[0],Selected_electron_phi[1]):-10")
               .Define("Selected_electron_miniPFRelIso_all","Electron_miniPFRelIso_all[goodElectrons]")
               .Define("Selected_electron_miniPFRelIso_chg","Electron_miniPFRelIso_chg[goodElectrons]")
               .Define("Selected_electron_pfRelIso03_all","Electron_pfRelIso03_all[goodElectrons]")
               .Define("Selected_electron_pfRelIso03_chg","Electron_pfRelIso03_chg[goodElectrons]")
               .Define("Selected_electron_HEEPconditions","Electron_vidNestedWPBitmapHEEP[goodElectrons]");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_electron_genPartIdx","Electron_genPartIdx[goodElectrons]")
                   .Define("Selected_electron_genPartpdgId","ROOT::VecOps::Take(GenPart_pdgId,Selected_electron_genPartIdx)")
                   .Define("Selected_electron_genPartIdxMother","ROOT::VecOps::Take(GenPart_genPartIdxMother,Selected_electron_genPartIdx)")
                   .Define("Selected_electron_genPartpdgIdMother","ROOT::VecOps::Take(GenPart_pdgId,Selected_electron_genPartIdxMother)")
                   .Define("Selected_electron_genPartFlav","Electron_genPartFlav[goodElectrons]");
    }
    _rlm = _rlm.Define("Selected_electron_pdgId","Electron_pdgId[goodElectrons]")
               .Define("Selected_electron_ip3d","Electron_ip3d[goodElectrons]")
               .Define("Selected_electron_sip3d","Electron_sip3d[goodElectrons]")
               .Define("Selected_electron_idx", ::good_idx,{"goodElectrons"});

    _rlm = _rlm.Define("elec4vecs", ::generate_4vec,{"Selected_electron_pt","Selected_electron_eta","Selected_electron_phi","Selected_electron_mass"})
               .Define("Vectorial_sum_two_electrons","(Selected_electron_number > 1) ? elec4vecs[0] + elec4vecs[1]:elec4vecs[0]")
               .Define("Vectorial_sum_two_electrons_pt","(Selected_electron_number > 1) ? Vectorial_sum_two_electrons.Pt():-100")
               .Define("Vectorial_sum_two_electrons_eta","(Selected_electron_number > 1) ? Vectorial_sum_two_electrons.Eta():-10")
               .Define("Vectorial_sum_two_electrons_phi","(Selected_electron_number > 1) ? Vectorial_sum_two_electrons.Phi():-10")
               .Define("Vectorial_sum_two_electrons_mass","(Selected_electron_number > 1) ? Vectorial_sum_two_electrons.M():-100");

}

//===============================Find Good Muons===============================================//
//: Define Good Muons in rdata frame
//=============================================================================================//
void TprimeAnalyser::selectMuons(std::string Region)
{

    cout << "select good muons" << endl;
    if (debug)
    {
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

    _rlm = _rlm.Define("goodMuonIDloose", MuonID(2));
    _rlm = _rlm.Define("goodMuonsloose","goodMuonIDloose && Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_miniPFRelIso_all < 0.40")
               .Define("Selected_muon_loose_pt","Muon_pt[goodMuonsloose]")
               .Define("Selected_muon_loose_eta","Muon_eta[goodMuonsloose]")
               .Define("Selected_muon_loose_phi","Muon_phi[goodMuonsloose]")
               .Define("Selected_muon_loose_mass","Muon_mass[goodMuonsloose]")
               .Define("Selected_muon_loose_number","int(Selected_muon_loose_pt.size())")
               .Define("muon4vecs_loose", ::generate_4vec,{"Selected_muon_loose_pt","Selected_muon_loose_eta","Selected_muon_loose_phi","Selected_muon_loose_mass"});

    _rlm = _rlm.Define("goodMuonIDlooseTTSL", MuonID(4));
    _rlm = _rlm.Define("goodMuonslooseTTSL","goodMuonsloose && goodMuonIDlooseTTSL && Muon_pt > 20 && abs(Muon_eta) < 2.4 && Muon_miniPFRelIso_all < 0.40 && Muon_sip3d < 3 && Muon_tightCharge == 2")
               .Define("Selected_muon_looseTTSL_pt","Muon_pt[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_number","int(Selected_muon_looseTTSL_pt.size())")
               .Define("Selected_muon_looseTTSL_leading_pt","(Selected_muon_looseTTSL_number > 0) ? Selected_muon_looseTTSL_pt[0]:-100")
               .Define("Selected_muon_looseTTSL_subleading_pt","(Selected_muon_looseTTSL_number > 1) ? Selected_muon_looseTTSL_pt[1]:-100")
               .Define("Selected_muon_looseTTSL_sum_two_muons_pt","(Selected_muon_looseTTSL_number>1) ? Selected_muon_looseTTSL_pt[0] + Selected_muon_looseTTSL_pt[1]:-100")
               .Define("Selected_muon_looseTTSL_sum_all_muons_pt","Sum(Selected_muon_looseTTSL_pt)")
               .Define("Selected_muon_looseTTSL_eta","Muon_eta[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_leading_eta","(Selected_muon_looseTTSL_number > 0) ? Selected_muon_looseTTSL_eta[0]:-10")
               .Define("Selected_muon_looseTTSL_subleading_eta","(Selected_muon_looseTTSL_number > 1) ? Selected_muon_looseTTSL_eta[1]:-10")
               .Define("Selected_muon_looseTTSL_phi","Muon_phi[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_mass","Muon_mass[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_transverse_energy","sqrt(Selected_muon_looseTTSL_pt*Selected_muon_looseTTSL_pt + Selected_muon_looseTTSL_mass*Selected_muon_looseTTSL_mass)")
               .Define("Selected_muon_looseTTSL_charge","Muon_charge[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_charge_sum","(Selected_muon_looseTTSL_number > 1) ? Selected_muon_looseTTSL_charge[0] + Selected_muon_looseTTSL_charge[1]:-100")
               .Define("Selected_muon_looseTTSL_deltaeta","(Selected_muon_looseTTSL_number > 1) ? abs(Selected_muon_looseTTSL_eta[1] - Selected_muon_looseTTSL_eta[0]):-10")
               .Define("Selected_muon_looseTTSL_deltaphi","(Selected_muon_looseTTSL_number > 1) ? abs(ROOT::VecOps::DeltaPhi(Selected_muon_looseTTSL_phi[0],Selected_muon_looseTTSL_phi[1])):-10")
               .Define("Selected_muon_looseTTSL_deltaR","(Selected_muon_looseTTSL_number >1) ? ROOT::VecOps::DeltaR(Selected_muon_looseTTSL_eta[0],Selected_muon_looseTTSL_eta[1],Selected_muon_looseTTSL_phi[0],Selected_muon_looseTTSL_phi[1]):-10")
               .Define("Selected_muon_looseTTSL_miniPFRelIso_all","Muon_miniPFRelIso_all[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_miniPFRelIso_chg","Muon_miniPFRelIso_chg[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_pfRelIso03_all","Muon_pfRelIso03_all[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_pfRelIso03_chg","Muon_pfRelIso03_chg[goodMuonslooseTTSL]");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_muon_looseTTSL_genPartIdx","Muon_genPartIdx[goodMuonslooseTTSL]")
                   .Define("Selected_muon_looseTTSL_genPartpdgId","ROOT::VecOps::Take(GenPart_pdgId,Selected_muon_looseTTSL_genPartIdx)")
                   .Define("Selected_muon_looseTTSL_genPartIdxMother","ROOT::VecOps::Take(GenPart_genPartIdxMother,Selected_muon_looseTTSL_genPartIdx)")
                   .Define("Selected_muon_looseTTSL_genPartpdgIdMother","ROOT::VecOps::Take(GenPart_pdgId,Selected_muon_looseTTSL_genPartIdxMother)")
                   .Define("Selected_muon_looseTTSL_genPartFlav","Muon_genPartFlav[goodMuonslooseTTSL]");
    }
    _rlm = _rlm.Define("Selected_muon_looseTTSL_pdgId","Muon_pdgId[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_ip3d","Muon_ip3d[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_sip3d","Muon_sip3d[goodMuonslooseTTSL]")
               .Define("Selected_muon_looseTTSL_idx", ::good_idx,{"goodMuonslooseTTSL"});

    _rlm = _rlm.Define("goodMuonID", MuonID(4));
    _rlm = _rlm.Define("goodMuons","goodMuonsloose && goodMuonID && Muon_pt > 20 && abs(Muon_eta) < 2.4 && Muon_miniPFRelIso_all < 0.05 && Muon_sip3d < 3 && Muon_tightCharge == 2") // && Muon_tightCharge == 2
               .Define("Selected_muon_pt","Muon_pt[goodMuons]")
               .Define("Selected_muon_number","int(Selected_muon_pt.size())")
               .Define("Selected_muon_leading_pt","(Selected_muon_number > 0) ? Selected_muon_pt[0]:-100")
               .Define("Selected_muon_subleading_pt","(Selected_muon_number > 1) ? Selected_muon_pt[1]:-100")
               .Define("Selected_muon_sum_two_muons_pt","(Selected_muon_number > 1) ? Selected_muon_pt[0] + Selected_muon_pt[1]:-100")
               .Define("Selected_muon_sum_all_muons_pt","Sum(Selected_muon_pt)")
               .Define("Selected_muon_eta","Muon_eta[goodMuons]")
               .Define("Selected_muon_leading_eta","(Selected_muon_number > 0) ? Selected_muon_eta[0]:-10")
               .Define("Selected_muon_subleading_eta","(Selected_muon_number > 1) ? Selected_muon_eta[1]:-10")
               .Define("Selected_muon_phi","Muon_phi[goodMuons]")
               .Define("Selected_muon_mass","Muon_mass[goodMuons]")
               .Define("Selected_muon_transverse_energy","sqrt(Selected_muon_pt*Selected_muon_pt + Selected_muon_mass*Selected_muon_mass)")
               .Define("Selected_muon_charge","Muon_charge[goodMuons]")
               .Define("Selected_muon_charge_sum","(Selected_muon_number > 1) ? Selected_muon_charge[0] + Selected_muon_charge[1]:-100")
               .Define("Selected_muon_deltaeta","(Selected_muon_number > 1) ? abs(Selected_muon_eta[1] - Selected_muon_eta[0]):-10")
               .Define("Selected_muon_deltaphi","(Selected_muon_number > 1) ? abs(ROOT::VecOps::DeltaPhi(Selected_muon_phi[0],Selected_muon_phi[1])):-10")
               .Define("Selected_muon_deltaR","(Selected_muon_number > 1) ? ROOT::VecOps::DeltaR(Selected_muon_eta[0],Selected_muon_eta[1],Selected_muon_phi[0],Selected_muon_phi[1]):-10")
               .Define("Selected_muon_miniPFRelIso_all","Muon_miniPFRelIso_all[goodMuons]")
               .Define("Selected_muon_miniPFRelIso_chg","Muon_miniPFRelIso_chg[goodMuons]")
               .Define("Selected_muon_pfRelIso03_all","Muon_pfRelIso03_all[goodMuons]")
               .Define("Selected_muon_pfRelIso03_chg","Muon_pfRelIso03_chg[goodMuons]")
               .Define("Selected_muon_pfRelIso04_all","Muon_pfRelIso04_all[goodMuons]");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_muon_genPartIdx","Muon_genPartIdx[goodMuons]")
                   .Define("Selected_muon_genPartpdgId","ROOT::VecOps::Take(GenPart_pdgId,Selected_muon_genPartIdx)")
                   .Define("Selected_muon_genPartIdxMother","ROOT::VecOps::Take(GenPart_genPartIdxMother,Selected_muon_genPartIdx)")
                   .Define("Selected_muon_genPartpdgIdMother","ROOT::VecOps::Take(GenPart_pdgId,Selected_muon_genPartIdxMother)")
                   .Define("Selected_muon_genPartFlav","Muon_genPartFlav[goodMuons]");
    }
    _rlm = _rlm.Define("Selected_muon_pdgId","Muon_pdgId[goodMuons]")
               .Define("Selected_muon_ip3d","Muon_ip3d[goodMuons]")
               .Define("Selected_muon_sip3d","Muon_sip3d[goodMuons]")
               .Define("Selected_muon_idx", ::good_idx,{"goodMuons"});

    _rlm = _rlm.Define("muon4vecs", ::generate_4vec,{"Selected_muon_pt","Selected_muon_eta","Selected_muon_phi","Selected_muon_mass"})
               .Define("Vectorial_sum_two_muons","(Selected_muon_number > 1) ? muon4vecs[0] + muon4vecs[1]:muon4vecs[0]")
               .Define("Vectorial_sum_two_muons_pt","(Selected_muon_number > 1) ? Vectorial_sum_two_muons.Pt():-100")
               .Define("Vectorial_sum_two_muons_eta","(Selected_muon_number > 1) ? Vectorial_sum_two_muons.Eta():-10")
               .Define("Vectorial_sum_two_muons_phi","(Selected_muon_number > 1) ? Vectorial_sum_two_muons.Phi():-10")
               .Define("Vectorial_sum_two_muons_mass","(Selected_muon_number > 1) ? Vectorial_sum_two_muons.M():-100");

}

void TprimeAnalyser::selectLeptons(std::string Region)
{
    _rlm = _rlm.Define("Selected_electron_muon_looseTTSL_pt","ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_pt,Selected_muon_looseTTSL_pt)")
               .Define("Selected_electron_muon_looseTTSL_pt_sort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Selected_electron_muon_looseTTSL_pt))")
               .Define("Selected_lepton_looseTTSL_pt","ROOT::VecOps::Take(Selected_electron_muon_looseTTSL_pt, Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_number","int(Selected_lepton_looseTTSL_pt.size())")
               .Define("Selected_lepton_looseTTSL_leading_pt","(Selected_lepton_looseTTSL_number > 0) ? Selected_lepton_looseTTSL_pt[0]:-100")
               .Define("Selected_lepton_looseTTSL_subleading_pt","(Selected_lepton_looseTTSL_number > 1) ? Selected_lepton_looseTTSL_pt[1]:-100")
               .Define("Selected_lepton_looseTTSL_sum_two_leptons_pt","(Selected_lepton_looseTTSL_number > 1) ? Selected_lepton_looseTTSL_pt[0] + Selected_lepton_looseTTSL_pt[1]:-100")
               .Define("Selected_lepton_looseTTSL_sum_all_leptons_pt","Sum(Selected_lepton_looseTTSL_pt)")
               .Define("Selected_lepton_looseTTSL_eta","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_eta,Selected_muon_looseTTSL_eta),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_leading_eta","(Selected_lepton_looseTTSL_number > 0) ? Selected_lepton_looseTTSL_eta[0]:-10")
               .Define("Selected_lepton_looseTTSL_subleading_eta","(Selected_lepton_looseTTSL_number > 1) ? Selected_lepton_looseTTSL_eta[1]:-10")
               .Define("Selected_lepton_looseTTSL_phi","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_phi,Selected_muon_looseTTSL_phi),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_mass","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_mass,Selected_muon_looseTTSL_mass),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_transverse_energy","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_transverse_energy,Selected_muon_looseTTSL_transverse_energy),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_charge","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_charge,Selected_muon_looseTTSL_charge),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_charge_sum","(Selected_lepton_looseTTSL_number > 1) ? Selected_lepton_looseTTSL_charge[0] + Selected_lepton_looseTTSL_charge[1]:-100")
               .Define("Selected_lepton_looseTTSL_deltaeta","(Selected_lepton_looseTTSL_number > 1) ? abs(Selected_lepton_looseTTSL_eta[1] - Selected_lepton_looseTTSL_eta[0]):-10")
               .Define("Selected_lepton_looseTTSL_deltaphi","(Selected_lepton_looseTTSL_number > 1) ? abs(ROOT::VecOps::DeltaPhi(Selected_lepton_looseTTSL_phi[0],Selected_lepton_looseTTSL_phi[1])):-10")
               .Define("Selected_lepton_looseTTSL_deltaR","(Selected_lepton_looseTTSL_number > 1) ? ROOT::VecOps::DeltaR(Selected_lepton_looseTTSL_eta[0],Selected_lepton_looseTTSL_eta[1],Selected_lepton_looseTTSL_phi[0],Selected_lepton_looseTTSL_phi[1]):-10")
               .Define("Selected_lepton_looseTTSL_miniPFRelIso_all","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_miniPFRelIso_all,Selected_muon_looseTTSL_miniPFRelIso_all),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_miniPFRelIso_chg","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_miniPFRelIso_chg,Selected_muon_looseTTSL_miniPFRelIso_chg),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_pfRelIso03_all","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_pfRelIso03_all,Selected_muon_looseTTSL_pfRelIso03_all),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_pfRelIso03_chg","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_pfRelIso03_chg,Selected_muon_looseTTSL_pfRelIso03_chg),Selected_electron_muon_looseTTSL_pt_sort)");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_lepton_looseTTSL_genPartIdx","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_genPartIdx,Selected_muon_looseTTSL_genPartIdx),Selected_electron_muon_looseTTSL_pt_sort)")
                   .Define("Selected_lepton_looseTTSL_genPartIdxMother","ROOT::VecOps::Take(GenPart_genPartIdxMother,Selected_lepton_looseTTSL_genPartIdx)")
                   .Define("Selected_lepton_looseTTSL_genPartpdgId","ROOT::VecOps::Take(GenPart_pdgId,Selected_lepton_looseTTSL_genPartIdx)")
                   .Define("Selected_lepton_looseTTSL_genPartpdgIdMother","ROOT::VecOps::Take(GenPart_pdgId,Selected_lepton_looseTTSL_genPartIdxMother)")
                   .Define("Selected_lepton_looseTTSL_genPartFlav","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_genPartFlav,Selected_muon_looseTTSL_genPartFlav),Selected_electron_muon_looseTTSL_pt_sort)");
    }
    _rlm = _rlm.Define("Selected_lepton_looseTTSL_pdgId","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_pdgId,Selected_muon_looseTTSL_pdgId),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton1_looseTTSL_pdgId","(Selected_lepton_looseTTSL_number > 1) ? abs(Selected_lepton_looseTTSL_pdgId[0]):0")
               .Define("Selected_lepton2_looseTTSL_pdgId","(Selected_lepton_looseTTSL_number > 1) ? abs(Selected_lepton_looseTTSL_pdgId[1]):0")
               .Define("Selected_lepton12_looseTTSL_pdgId","(Selected_lepton_looseTTSL_number > 1) ? abs(Selected_lepton_looseTTSL_pdgId[0]) + abs(Selected_lepton_looseTTSL_pdgId[1]):0")
               .Define("Selected_lepton_looseTTSL_ip3d","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_ip3d,Selected_muon_looseTTSL_ip3d),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_sip3d","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_sip3d,Selected_muon_looseTTSL_sip3d),Selected_electron_muon_looseTTSL_pt_sort)")
               .Define("Selected_lepton_looseTTSL_idx","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_looseTTSL_idx,Selected_muon_looseTTSL_idx),Selected_electron_muon_looseTTSL_pt_sort)");

    _rlm = _rlm.Define("lep4vecs_looseTTSL", ::generate_4vec,{"Selected_lepton_looseTTSL_pt","Selected_lepton_looseTTSL_eta","Selected_lepton_looseTTSL_phi","Selected_lepton_looseTTSL_mass"})
               .Define("lep4Rvecs_looseTTSL", ::generate_4Rvec,{"Selected_lepton_looseTTSL_pt","Selected_lepton_looseTTSL_eta","Selected_lepton_looseTTSL_phi","Selected_lepton_looseTTSL_mass"})
               .Define("Vectorial_sum_two_leptons_looseTTSL","(Selected_lepton_looseTTSL_number > 1) ? lep4vecs_looseTTSL[0] + lep4vecs_looseTTSL[1]:lep4vecs_looseTTSL[0]")
               .Define("Vectorial_sum_two_leptons_looseTTSL_pt","(Selected_lepton_looseTTSL_number > 1) ? (lep4vecs_looseTTSL[0] + lep4vecs_looseTTSL[1]).Pt():-100")
               .Define("Vectorial_sum_two_leptons_looseTTSL_eta","(Selected_lepton_looseTTSL_number > 1) ? (lep4vecs_looseTTSL[0] + lep4vecs_looseTTSL[1]).Eta():-10")
               .Define("Vectorial_sum_two_leptons_looseTTSL_phi","(Selected_lepton_looseTTSL_number > 1) ? (lep4vecs_looseTTSL[0] + lep4vecs_looseTTSL[1]).Phi():-10")
               .Define("Vectorial_sum_two_leptons_looseTTSL_mass","(Selected_lepton_looseTTSL_number > 1) ? (lep4vecs_looseTTSL[0] + lep4vecs_looseTTSL[1]).M():-100");

    _rlm = _rlm.Define("Selected_electron_muon_pt","ROOT::VecOps::Concatenate(Selected_electron_pt,Selected_muon_pt)")
               .Define("Selected_electron_muon_pt_sort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Selected_electron_muon_pt))")
               .Define("Selected_lepton_pt","ROOT::VecOps::Take(Selected_electron_muon_pt, Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_number","int(Selected_lepton_pt.size())")
               .Define("Selected_lepton_leading_pt","(Selected_lepton_number > 0) ? Selected_lepton_pt[0]:-100")
               .Define("Selected_lepton_subleading_pt","(Selected_lepton_number > 1) ? Selected_lepton_pt[1]:-100")
               .Define("Selected_lepton_sum_two_leptons_pt","(Selected_lepton_number > 1) ? Selected_lepton_pt[0] + Selected_lepton_pt[1]:-100")
               .Define("Selected_lepton_sum_all_leptons_pt","Sum(Selected_lepton_pt)")
               .Define("Selected_lepton_eta","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_eta,Selected_muon_eta),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_leading_eta","(Selected_lepton_number > 0) ? Selected_lepton_eta[0]:-10")
               .Define("Selected_lepton_subleading_eta","(Selected_lepton_number > 1) ? Selected_lepton_eta[1]:-10")
               .Define("Selected_lepton_phi","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_phi,Selected_muon_phi),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_mass","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_mass,Selected_muon_mass),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_transverse_energy","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_transverse_energy,Selected_muon_transverse_energy),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_charge","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_charge,Selected_muon_charge),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_charge_sum","(Selected_lepton_number > 1) ? Selected_lepton_charge[0] + Selected_lepton_charge[1]:-100")
               .Define("Selected_lepton_deltaeta","(Selected_lepton_number > 1) ? abs(Selected_lepton_eta[1] - Selected_lepton_eta[0]):-10")
               .Define("Selected_lepton_deltaphi","(Selected_lepton_number > 1) ? abs(ROOT::VecOps::DeltaPhi(Selected_lepton_phi[0],Selected_lepton_phi[1])):-10")
               .Define("Selected_lepton_deltaR","(Selected_lepton_number > 1) ? ROOT::VecOps::DeltaR(Selected_lepton_eta[0],Selected_lepton_eta[1],Selected_lepton_phi[0],Selected_lepton_phi[1]):-10")
               .Define("Selected_lepton_miniPFRelIso_all","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_miniPFRelIso_all,Selected_muon_miniPFRelIso_all),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_miniPFRelIso_chg","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_miniPFRelIso_chg,Selected_muon_miniPFRelIso_chg),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_pfRelIso03_all","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_pfRelIso03_all,Selected_muon_pfRelIso03_all),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_pfRelIso03_chg","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_pfRelIso03_chg,Selected_muon_pfRelIso03_chg),Selected_electron_muon_pt_sort)");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_lepton_genPartIdx","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_genPartIdx,Selected_muon_genPartIdx),Selected_electron_muon_pt_sort)")
                   .Define("Selected_lepton_genPartIdxMother","ROOT::VecOps::Take(GenPart_genPartIdxMother,Selected_lepton_genPartIdx)")
                   .Define("Selected_lepton_genPartpdgId","ROOT::VecOps::Take(GenPart_pdgId,Selected_lepton_genPartIdx)")
                   .Define("Selected_lepton_genPartpdgIdMother","ROOT::VecOps::Take(GenPart_pdgId,Selected_lepton_genPartIdxMother)")
                   .Define("Selected_lepton_genPartFlav","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_genPartFlav,Selected_muon_genPartFlav),Selected_electron_muon_pt_sort)");
    }
    _rlm = _rlm.Define("Selected_lepton_pdgId","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_pdgId,Selected_muon_pdgId),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton1_pdgId","(Selected_lepton_number > 0) ? abs(Selected_lepton_pdgId[0]):-1000")
               .Define("Selected_lepton2_pdgId","(Selected_lepton_number > 1) ? abs(Selected_lepton_pdgId[1]):-1000")
               .Define("Selected_lepton12_pdgId","(Selected_lepton_number > 1) ? abs(Selected_lepton_pdgId[0]) + abs(Selected_lepton_pdgId[1]):-1000")
               .Define("Selected_lepton_ip3d","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_ip3d,Selected_muon_ip3d),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_sip3d","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_sip3d,Selected_muon_sip3d),Selected_electron_muon_pt_sort)")
               .Define("Selected_lepton_idx","ROOT::VecOps::Take(ROOT::VecOps::Concatenate(Selected_electron_idx,Selected_muon_idx),Selected_electron_muon_pt_sort)");

    _rlm = _rlm.Define("lep4vecs", ::generate_4vec,{"Selected_lepton_pt","Selected_lepton_eta","Selected_lepton_phi","Selected_lepton_mass"})
               .Define("lep4Rvecs", ::generate_4Rvec,{"Selected_lepton_pt","Selected_lepton_eta","Selected_lepton_phi","Selected_lepton_mass"})
               .Define("Selected_lepton_energy","ROOT::VecOps::Map(lep4Rvecs, FourVecEnergy)")
               .Define("Selected_lepton_px","ROOT::VecOps::Map(lep4Rvecs, FourVecPx)")
               .Define("Selected_lepton_py","ROOT::VecOps::Map(lep4Rvecs, FourVecPy)")
               .Define("Selected_lepton_pz","ROOT::VecOps::Map(lep4Rvecs, FourVecPz)")
               .Define("Vectorial_sum_two_leptons","(Selected_lepton_number > 1) ? lep4vecs[0] + lep4vecs[1]:lep4vecs[0]")
               .Define("Vectorial_sum_two_leptons_pt","(Selected_lepton_number > 1) ? Vectorial_sum_two_leptons.Pt():-100")
               .Define("Vectorial_sum_two_leptons_eta","(Selected_lepton_number > 1) ? Vectorial_sum_two_leptons.Eta():-10")
               .Define("Vectorial_sum_two_leptons_phi","(Selected_lepton_number > 1) ? Vectorial_sum_two_leptons.Phi():-10")
               .Define("Vectorial_sum_two_leptons_mass","(Selected_lepton_number > 1) ? Vectorial_sum_two_leptons.M():-100")
               .Define("Mass_dilepton_false",::false_mass_dilepton,{"lep4vecs"}); 
}

void TprimeAnalyser::selectMET()
{
    if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

//   _rlm = _rlm.Define("goodMET","MET_pt > 0")
//              .Define("Selected_MET_pt","MET_pt[goodMET]")
//              .Define("Selected_MET_sum_pt_unclustered","MET_sumPtUnclustered[goodMET]")
//              .Define("Selected_MET_phi","MET_phi[goodMET]")
//              .Define("Selected_MET_sum_transverse_energy","MET_sumEt[goodMET]");

    _rlm = _rlm.Define("MET_px","MET_pt*cos(MET_phi)")
               .Define("MET_py","MET_pt*sin(MET_phi)");

}

void TprimeAnalyser::selectReconstructedLeptons()
{
    _rlm = _rlm.Define("sel_leptongentopHiggs",::isLfromtopHiggs,{"Selected_lepton_pdgId","GenPart_pdgId","Selected_lepton_genPartIdx","GenPart_genPartIdxMother"})//, "Selected_lepton_pt"})
               .Define("LeptonFromtopHiggs_pt","Selected_lepton_pt[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_leading_pt","LeptonFromtopHiggs_pt[0]")
               .Define("LeptonFromtopHiggs_subleading_pt","LeptonFromtopHiggs_pt[1]")
               .Define("LeptonFromtopHiggs_sum_two_muons_pt","LeptonFromtopHiggs_pt[0] + LeptonFromtopHiggs_pt[1]")
               .Define("LeptonFromtopHiggs_sum_all_muons_pt","Sum(LeptonFromtopHiggs_pt)")
               .Define("LeptonFromtopHiggs_eta","Selected_lepton_eta[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_phi","Selected_lepton_phi[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_mass","Selected_lepton_mass[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_transverse_energy","sqrt(LeptonFromtopHiggs_pt*LeptonFromtopHiggs_pt + LeptonFromtopHiggs_mass*LeptonFromtopHiggs_mass)")
               .Define("LeptonFromtopHiggs_charge","Selected_lepton_charge[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_charge_sum","LeptonFromtopHiggs_charge[0] + LeptonFromtopHiggs_charge[1]")
               .Define("LeptonFromtopHiggs_number","int(LeptonFromtopHiggs_pt.size())")
               .Define("LeptonFromtopHiggs_deltaeta","abs(LeptonFromtopHiggs_eta[1] - LeptonFromtopHiggs_eta[0])")
               .Define("LeptonFromtopHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtopHiggs_phi[0],LeptonFromtopHiggs_phi[1]))")
	           .Define("LeptonFromtopHiggs_deltaR","ROOT::VecOps::DeltaR(LeptonFromtopHiggs_eta[0],LeptonFromtopHiggs_eta[1],LeptonFromtopHiggs_phi[0],LeptonFromtopHiggs_phi[1])")
               .Define("LeptonFromtopHiggs_miniPFRelIso_all","Selected_lepton_miniPFRelIso_all[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_miniPFRelIso_chg","Selected_lepton_miniPFRelIso_chg[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_pfRelIso03_all","Selected_lepton_pfRelIso03_all[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_pfRelIso03_chg","Selected_lepton_pfRelIso03_chg[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_genPartIdx","Selected_lepton_genPartIdx[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_genPartFlav","Selected_lepton_genPartFlav[sel_leptongentopHiggs]")
               .Define("LeptonFromtopHiggs_pdgId","Selected_lepton_pdgId[sel_leptongentopHiggs]")
               .Define("p4_LeptonFromtopHiggs",::generate_4vec,{"LeptonFromtopHiggs_pt","LeptonFromtopHiggs_eta","LeptonFromtopHiggs_phi","LeptonFromtopHiggs_mass"});

    _rlm = _rlm.Define("sel_leptongentoporHiggs",::isLfromtoporHiggs,{"LeptonFromtopHiggs_pdgId","GenPart_pdgId","LeptonFromtopHiggs_genPartIdx","GenPart_genPartIdxMother","p4_LeptonFromtopHiggs"})
               .Define("LeptonFromtop4vecs","sel_leptongentoporHiggs[0]")
               .Define("LeptonFromtop_pt","LeptonFromtop4vecs.Pt()")
               .Define("LeptonFromtop_eta","LeptonFromtop4vecs.Eta()")
               .Define("LeptonFromtop_phi","LeptonFromtop4vecs.Phi()")
               .Define("LeptonFromtop_mass","LeptonFromtop4vecs.M()")
               .Define("LeptonFromtop_transverse_energy","sqrt(LeptonFromtop_pt*LeptonFromtop_pt + LeptonFromtop_mass*LeptonFromtop_mass)")
               .Define("LeptonFromtop_energy","LeptonFromtop4vecs.E()")
               .Define("LeptonFromtop_px","LeptonFromtop4vecs.Px()")
               .Define("LeptonFromtop_py","LeptonFromtop4vecs.Py()")
               .Define("LeptonFromtop_pz","LeptonFromtop4vecs.Pz()")
               .Define("LeptonFromHiggs4vecs","sel_leptongentoporHiggs[1]")
               .Define("LeptonFromHiggs_pt","LeptonFromHiggs4vecs.Pt()")
               .Define("LeptonFromHiggs_eta","LeptonFromHiggs4vecs.Eta()")
               .Define("LeptonFromHiggs_phi","LeptonFromHiggs4vecs.Phi()")
               .Define("LeptonFromHiggs_mass","LeptonFromHiggs4vecs.M()")
               .Define("LeptonFromHiggs_transverse_energy","sqrt(LeptonFromHiggs_pt*LeptonFromHiggs_pt + LeptonFromHiggs_mass*LeptonFromHiggs_mass)")
               .Define("LeptonFromHiggs_energy","LeptonFromHiggs4vecs.E()")
               .Define("LeptonFromHiggs_px","LeptonFromHiggs4vecs.Px()")
               .Define("LeptonFromHiggs_py","LeptonFromHiggs4vecs.Py()")
               .Define("LeptonFromHiggs_pz","LeptonFromHiggs4vecs.Pz()")
               .Define("LeptonFromtopHiggs_ordered_pt",::concatenate,{"LeptonFromtop_pt","LeptonFromHiggs_pt"})
               .Define("LeptonFromtopHiggs_ordered_px",::concatenate,{"LeptonFromtop_px","LeptonFromHiggs_px"})
               .Define("LeptonFromtopHiggs_ordered_py",::concatenate,{"LeptonFromtop_py","LeptonFromHiggs_py"})
               .Define("LeptonFromtopHiggs_ordered_pz",::concatenate,{"LeptonFromtop_pz","LeptonFromHiggs_pz"})
               .Define("LeptonFromtopHiggs_ordered_energy",::concatenate,{"LeptonFromtop_energy","LeptonFromHiggs_energy"})
               .Define("LeptonFromtopHiggs_ordered_mass",::concatenate,{"LeptonFromtop_mass","LeptonFromHiggs_mass"})
               .Define("LeptonFromtop_LeptonFromHiggs","LeptonFromtop4vecs + LeptonFromHiggs4vecs")
               .Define("LeptonFromtop_LeptonFromHiggs_pt","LeptonFromtop_LeptonFromHiggs.Pt()")
               .Define("LeptonFromtop_LeptonFromHiggs_eta","LeptonFromtop_LeptonFromHiggs.Eta()")
               .Define("LeptonFromtop_LeptonFromHiggs_phi","LeptonFromtop_LeptonFromHiggs.Phi()")
               .Define("LeptonFromtop_LeptonFromHiggs_mass","LeptonFromtop_LeptonFromHiggs.M()")
               .Define("LeptonFromtop_LeptonFromHiggs_transverse_energy","sqrt(LeptonFromtop_LeptonFromHiggs_pt*LeptonFromtop_LeptonFromHiggs_pt + LeptonFromtop_LeptonFromHiggs_mass*LeptonFromtop_LeptonFromHiggs_mass)")
               .Define("LeptonFromtop_LeptonFromHiggs_deltaeta","abs(LeptonFromtop_eta - LeptonFromHiggs_eta)")
               .Define("LeptonFromtop_LeptonFromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtop_phi,LeptonFromHiggs_phi))")
               .Define("LeptonFromtop_LeptonFromHiggs_deltaR","ROOT::VecOps::DeltaR(LeptonFromtop_eta,LeptonFromHiggs_eta,LeptonFromtop_phi,LeptonFromHiggs_phi)");

    // _rlm = _rlm.Define("sel_leptongenW",::isLfromW,{"Selected_lepton_pdgId","GenPart_pdgId","Selected_lepton_genPartIdx","GenPart_genPartIdxMother"})
    //            .Define("LeptonFromW_pt","Selected_lepton_pt[sel_leptongenW]")
    //            .Define("LeptonFromW_leading_pt","LeptonFromW_pt[0]")
    //            .Define("LeptonFromW_subleading_pt","LeptonFromW_pt[1]")
    //            .Define("LeptonFromW_sum_two_muons_pt","LeptonFromW_pt[0] + LeptonFromW_pt[1]")
    //            .Define("LeptonFromW_sum_all_muons_pt","Sum(LeptonFromW_pt)")
    //            .Define("LeptonFromW_eta","Selected_lepton_eta[sel_leptongenW]")
    //            .Define("LeptonFromW_phi","Selected_lepton_phi[sel_leptongenW]")
    //            .Define("LeptonFromW_mass","Selected_lepton_mass[sel_leptongenW]")
    //            .Define("LeptonFromW_charge","Selected_lepton_charge[sel_leptongenW]")
    //            .Define("LeptonFromW_charge_sum","LeptonFromW_charge[0] + LeptonFromW_charge[1]")
    //            .Define("LeptonFromW_number","int(LeptonFromW_pt.size())")
    //            .Define("LeptonFromW_deltaeta","abs(LeptonFromW_eta[1] - LeptonFromW_eta[0])")
    //            .Define("LeptonFromW_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromW_phi[0],LeptonFromW_phi[1]))")
    // 	          .Define("LeptonFromW_deltaR","ROOT::VecOps::DeltaR(LeptonFromW_eta[0],LeptonFromW_eta[1],LeptonFromW_phi[0],LeptonFromW_phi[1])")
    //            .Define("LeptonFromW_genPartIdx","Selected_lepton_genPartIdx[sel_leptongenW]")
    //            .Define("LeptonFromW_genPartFlav","Selected_lepton_genPartFlav[sel_leptongenW]")
    //            .Define("LeptonFromW_pdgId","Selected_lepton_pdgId[sel_leptongenW]")
    //            .Define("p4_LeptonFromW",::generate_4vec,{"LeptonFromW_pt","LeptonFromW_eta","LeptonFromW_phi","LeptonFromW_mass"});

    _rlm = _rlm.Define("sel_muongenW",::isLfromW,{"Selected_muon_pdgId","GenPart_pdgId","Selected_muon_genPartIdx","GenPart_genPartIdxMother"})
               .Define("MuonFromW_pt","Selected_muon_pt[sel_muongenW]")
               .Define("MuonFromW_charge","Selected_muon_charge[sel_muongenW]")
               .Define("MuonFromW_number","int(MuonFromW_pt.size())");

    _rlm = _rlm.Define("sel_electrongenW",::isLfromW,{"Selected_electron_pdgId","GenPart_pdgId","Selected_electron_genPartIdx","GenPart_genPartIdxMother"})
               .Define("ElectronFromW_pt","Selected_electron_pt[sel_electrongenW]")
               .Define("ElectronFromW_charge","Selected_electron_charge[sel_electrongenW]")
               .Define("ElectronFromW_number","int(ElectronFromW_pt.size())");

    _rlm = _rlm.Define("sel_leptongenMother",::MotherfromL,{"Selected_lepton_pdgId","GenPart_pdgId","Selected_lepton_genPartIdx","GenPart_genPartIdxMother"})
               .Define("MotherfromLepton_pdgId","sel_leptongenMother")
               .Define("MotherfromLepton1_pdgId","abs(sel_leptongenMother[0])")     
               .Define("MotherfromLepton2_pdgId","abs(sel_leptongenMother[1])")
               .Define("MotherfromLepton12_pdgId","abs(sel_leptongenMother[0]) + abs(sel_leptongenMother[1])");

    // _rlm = _rlm.Define("GenLeptons","abs(GenPart_pdgId) == 13 || abs(GenPart_pdgId) == 11")
    //            .Define("MotherGenLeptons","GenPart_genPartIdxMother[GenLeptons]")
    //            .Define("MotherGenLeptonsPdgId","ROOT::VecOps::Take(GenPart_pdgId, MotherGenLeptons)")
    //            .Define("IsFromW","abs(MotherGenLeptonsPdgId) == 24");

    // _rlm = _rlm.Define("GenMuons","abs(GenPart_pdgId) == 13")
    //            .Define("GenMuonsPdgId","GenPart_pdgId[GenMuons]")
    //            .Define("MotherGenMuons","GenPart_genPartIdxMother[GenMuons]")
    //            .Define("MotherGenMuonsPdgId","ROOT::VecOps::Take(GenPart_pdgId, MotherGenMuons)")
    //            .Define("IsMuonFromW","abs(MotherGenMuonsPdgId) == 24");

    // _rlm = _rlm.Define("GenElectrons","abs(GenPart_pdgId) == 11")
    //            .Define("GenElectronsPdgId","GenPart_pdgId[GenElectrons]")
    //            .Define("MotherGenElectrons","GenPart_genPartIdxMother[GenElectrons]")
    //            .Define("MotherGenElectronsPdgId","ROOT::VecOps::Take(GenPart_pdgId, MotherGenElectrons)")
    //            .Define("IsElectronFromW","abs(MotherGenElectronsPdgId) == 24");

    // _rlm = _rlm.Define("GenW","abs(GenPart_pdgId) == 24")
    //            .Define("MotherGenW","GenPart_genPartIdxMother[GenW]")
    //            .Define("MotherGenWPdgId","ROOT::VecOps::Take(GenPart_pdgId, MotherGenW)")
    //            .Define("IsFromHiggs","abs(MotherGenWPdgId) == 25")
    //            .Define("IsFromtop","abs(MotherGenWPdgId) == 6");

    _rlm = _rlm.Define("GenNeutrinos","abs(GenPart_pdgId) == 14 || abs(GenPart_pdgId) == 12")
               .Define("GenNeutrinosPdgId","GenPart_pdgId[GenNeutrinos]")
               .Define("GenNeutrinos_pt","GenPart_pt[GenNeutrinos]")
               .Define("GenNeutrinos_eta","GenPart_eta[GenNeutrinos]")
               .Define("GenNeutrinos_phi","GenPart_phi[GenNeutrinos]")
               .Define("MotherGenNeutrinos","GenPart_genPartIdxMother[GenNeutrinos]")
               .Define("MotherGenNeutrinosPdgId","ROOT::VecOps::Take(GenPart_pdgId, MotherGenNeutrinos)")
               .Define("IsNeutrinoFromW","abs(MotherGenNeutrinosPdgId) == 24")
               .Define("GenNeutrinoFromWPdgId","GenNeutrinosPdgId[IsNeutrinoFromW]")
               .Define("GenNeutrinoFromW_pt","GenNeutrinos_pt[IsNeutrinoFromW]")
               .Define("GenNeutrinoFromW_eta","GenNeutrinos_eta[IsNeutrinoFromW]")
               .Define("GenNeutrinoFromW_phi","GenNeutrinos_phi[IsNeutrinoFromW]")
               .Define("GenNeutrinoFromW_number","int(GenNeutrinoFromW_pt.size())")
               .Define("MotherNeutrinoFromW","MotherGenNeutrinos[IsNeutrinoFromW]")
               .Define("MotherNeutrinoFromWPdgId","ROOT::VecOps::Take(GenPart_pdgId, MotherNeutrinoFromW)")
               .Define("MotherNeutrinoFromWMass","ROOT::VecOps::Take(GenPart_mass, MotherNeutrinoFromW)")
               .Define("sel_NeutrinoFromWgenHiggs",::isWfromHiggs,{"GenPart_pdgId","MotherNeutrinoFromW","GenPart_genPartIdxMother","GenPart_mass"})
               .Define("sel_NeutrinoFromWgentop",::isWfromtop,{"GenPart_pdgId","MotherNeutrinoFromW","GenPart_genPartIdxMother","GenPart_mass"})
               .Define("GenNeutrinoFromtop_pt","GenNeutrinoFromW_pt[sel_NeutrinoFromWgentop]")
               .Define("GenNeutrinoFromtop_eta","GenNeutrinoFromW_eta[sel_NeutrinoFromWgentop]")
               .Define("GenNeutrinoFromtop_phi","GenNeutrinoFromW_phi[sel_NeutrinoFromWgentop]")
               .Define("GenNeutrinoFromtop_px","GenNeutrinoFromtop_pt*cos(GenNeutrinoFromtop_phi)")
               .Define("GenNeutrinoFromtop_py","GenNeutrinoFromtop_pt*sin(GenNeutrinoFromtop_phi)")
               .Define("p4_GenNeutrinoFromtop",::generate_4vec_3,{"GenNeutrinoFromtop_pt","GenNeutrinoFromtop_eta","GenNeutrinoFromtop_phi"})
               .Define("p4_GenNeutrinoFromtop_first","p4_GenNeutrinoFromtop[0]")
               .Define("GenNeutrinoFromtop_first_pt","p4_GenNeutrinoFromtop_first.Pt()")
               .Define("GenNeutrinoFromtop_first_eta","p4_GenNeutrinoFromtop_first.Eta()")
               .Define("GenNeutrinoFromtop_first_phi","p4_GenNeutrinoFromtop_first.Phi()")
               .Define("GenNeutrinoFromHiggs_pt","GenNeutrinoFromW_pt[sel_NeutrinoFromWgenHiggs]")
               .Define("GenNeutrinoFromHiggs_eta","GenNeutrinoFromW_eta[sel_NeutrinoFromWgenHiggs]")
               .Define("GenNeutrinoFromHiggs_phi","GenNeutrinoFromW_phi[sel_NeutrinoFromWgenHiggs]")
               .Define("GenNeutrinoFromHiggs_px","GenNeutrinoFromHiggs_pt*cos(GenNeutrinoFromHiggs_phi)")
               .Define("GenNeutrinoFromHiggs_py","GenNeutrinoFromHiggs_pt*sin(GenNeutrinoFromHiggs_phi)")
               .Define("p4_GenNeutrinoFromHiggs",::generate_4vec_3,{"GenNeutrinoFromHiggs_pt","GenNeutrinoFromHiggs_eta","GenNeutrinoFromHiggs_phi"})
               .Define("p4_GenNeutrinoFromHiggs_first","p4_GenNeutrinoFromHiggs[0]")
               .Define("GenNeutrinoFromHiggs_first_pt","p4_GenNeutrinoFromHiggs_first.Pt()")
               .Define("GenNeutrinoFromHiggs_first_eta","p4_GenNeutrinoFromHiggs_first.Eta()")
               .Define("GenNeutrinoFromHiggs_first_phi","p4_GenNeutrinoFromHiggs_first.Phi()")
               .Define("MotherNeutrinoFromHiggsPdgid","MotherNeutrinoFromWPdgId[sel_NeutrinoFromWgenHiggs]")
               .Define("MotherNeutrinoFromHiggsMass","MotherNeutrinoFromWMass[sel_NeutrinoFromWgenHiggs]")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs","p4_GenNeutrinoFromtop_first + p4_GenNeutrinoFromHiggs_first")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs_pt","NeutrinoFromtop_NeutrinoFromHiggs.Pt()")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs_eta","NeutrinoFromtop_NeutrinoFromHiggs.Eta()")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs_phi","NeutrinoFromtop_NeutrinoFromHiggs.Phi()")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs_mass","NeutrinoFromtop_NeutrinoFromHiggs.M()")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs_transverse_energy","sqrt(NeutrinoFromtop_NeutrinoFromHiggs_pt*NeutrinoFromtop_NeutrinoFromHiggs_pt + NeutrinoFromtop_NeutrinoFromHiggs_mass*NeutrinoFromtop_NeutrinoFromHiggs_mass)")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs_deltaeta","abs(GenNeutrinoFromtop_first_eta - GenNeutrinoFromHiggs_first_eta)")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(GenNeutrinoFromtop_first_phi,GenNeutrinoFromHiggs_first_phi))")
               .Define("NeutrinoFromtop_NeutrinoFromHiggs_deltaR","ROOT::VecOps::DeltaR(GenNeutrinoFromtop_first_eta,GenNeutrinoFromHiggs_first_eta,GenNeutrinoFromtop_first_phi,GenNeutrinoFromHiggs_first_phi)");

    _rlm = _rlm.Define("NeutrinoFromtop_LeptonFromtop","p4_GenNeutrinoFromtop_first + LeptonFromtop4vecs")
               .Define("NeutrinoFromtop_LeptonFromtop_pt","NeutrinoFromtop_LeptonFromtop.Pt()")
               .Define("NeutrinoFromtop_LeptonFromtop_eta","NeutrinoFromtop_LeptonFromtop.Eta()")
               .Define("NeutrinoFromtop_LeptonFromtop_phi","NeutrinoFromtop_LeptonFromtop.Phi()")
               .Define("NeutrinoFromtop_LeptonFromtop_mass","NeutrinoFromtop_LeptonFromtop.M()")
               .Define("NeutrinoFromtop_LeptonFromtop_transverse_energy","sqrt(NeutrinoFromtop_LeptonFromtop_pt*NeutrinoFromtop_LeptonFromtop_pt + NeutrinoFromtop_LeptonFromtop_mass*NeutrinoFromtop_LeptonFromtop_mass)")
               .Define("NeutrinoFromtop_LeptonFromtop_deltaeta","abs(GenNeutrinoFromtop_first_eta - LeptonFromtop_eta)")
               .Define("NeutrinoFromtop_LeptonFromtop_deltaphi","abs(ROOT::VecOps::DeltaPhi(GenNeutrinoFromtop_first_phi,LeptonFromtop_phi))")
               .Define("NeutrinoFromtop_LeptonFromtop_deltaR","ROOT::VecOps::DeltaR(GenNeutrinoFromtop_first_eta,LeptonFromtop_eta,GenNeutrinoFromtop_first_phi,LeptonFromtop_phi)")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs","p4_GenNeutrinoFromHiggs_first + LeptonFromHiggs4vecs")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs_pt","NeutrinoFromHiggs_LeptonFromHiggs.Pt()")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs_eta","NeutrinoFromHiggs_LeptonFromHiggs.Eta()")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs_phi","NeutrinoFromHiggs_LeptonFromHiggs.Phi()")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs_mass","NeutrinoFromHiggs_LeptonFromHiggs.M()")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs_transverse_energy","sqrt(NeutrinoFromHiggs_LeptonFromHiggs_pt*NeutrinoFromHiggs_LeptonFromHiggs_pt + NeutrinoFromHiggs_LeptonFromHiggs_mass*NeutrinoFromHiggs_LeptonFromHiggs_mass)")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs_deltaeta","abs(GenNeutrinoFromHiggs_first_eta - LeptonFromHiggs_eta)")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(GenNeutrinoFromHiggs_first_phi,LeptonFromHiggs_phi))")
               .Define("NeutrinoFromHiggs_LeptonFromHiggs_deltaR","ROOT::VecOps::DeltaR(GenNeutrinoFromHiggs_first_eta,LeptonFromHiggs_eta,GenNeutrinoFromHiggs_first_phi,LeptonFromHiggs_phi)");

    // _rlm = _rlm.Define("GenLeptons","abs(GenPart_pdgId) == 13 || abs(GenPart_pdgId) == 11")
    //            .Define("GenLeptonsPdgId","GenPart_pdgId[GenLeptons]")
    //            .Define("GenLeptons_pt","GenPart_pt[GenLeptons]")
    //            .Define("GenLeptons_eta","GenPart_eta[GenLeptons]")
    //            .Define("GenLeptons_phi","GenPart_phi[GenLeptons]");

}

//=================================Select Jets=================================================//
//check the twiki page :    https://twiki.cern.ch/twiki/bin/view/CMS/JetID
//to find jetId working points for the purpose of  your analysis.
    //jetId==2 means: pass tight ID, fail tightLepVeto
    //jetId==6 means: pass tight ID and tightLepVeto ID.
//=============================================================================================//

void TprimeAnalyser::selectJets(std::string Region)
{

    cout << "select good jets" << endl;
    if (debug)
    {
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

    _rlm = _rlm.Define("goodJets","Jet_pt > 30 && abs(Jet_eta) < 4.5 && Jet_jetId >= 6")
               .Define("Selected_jet_pt","Jet_pt[goodJets]")
               .Define("Selected_jet_leading_pt","(Selected_jet_pt.size()>0) ? Selected_jet_pt[0]:-100")
               .Define("Selected_jet_subleading_pt","(Selected_jet_pt.size()>1) ? Selected_jet_pt[1]:-100")
               .Define("Selected_jet_subsubleading_pt","(Selected_jet_pt.size()>2) ? Selected_jet_pt[2]:-100")
               .Define("Selected_jet_Ht","Sum(Selected_jet_pt)")
               .Define("Selected_jet_eta","Jet_eta[goodJets]")
               .Define("Selected_jet_phi","Jet_phi[goodJets]")
               .Define("Selected_jet_mass","Jet_mass[goodJets]")
               .Define("Selected_jet_number","int(Selected_jet_pt.size())")
               .Define("Selected_jet_btagDeepB","Jet_btagDeepB[goodJets]")
               .Define("Selected_jet_btagDeepFlavB","Jet_btagDeepFlavB[goodJets]");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_jet_genJetIdx","Jet_genJetIdx[goodJets]")
                   .Define("Selected_jet_partonFlavour","Jet_partonFlavour[goodJets]")
                   .Define("Selected_jet_hadronFlavour","Jet_hadronFlavour[goodJets]");
    }
    _rlm = _rlm.Define("Selected_jet_jetId","Jet_jetId[goodJets]")      
               .Define("jet4vecs", ::generate_4vec,{"Selected_jet_pt","Selected_jet_eta","Selected_jet_phi","Selected_jet_mass"});

    _rlm = _rlm.Define("goodJets_btag","goodJets && abs(Jet_eta) < 2.5 && Jet_btagDeepFlavB > 0.2783") //Jet_btagDeepB > 0.1208") DeepCSV loose = 0.1208, DeepJet loose = 0.0490, DeepJet medium = 0.2783
               .Define("Selected_bjet_pt","Jet_pt[goodJets_btag]")
               .Define("Selected_bjet_leading_pt","(Selected_bjet_pt.size()>0) ? Selected_bjet_pt[0]:-100")
               .Define("Selected_bjet_subleading_pt","(Selected_bjet_pt.size()>1) ? Selected_bjet_pt[1]:-100")
               .Define("Selected_bjet_subsubleading_pt","(Selected_bjet_pt.size()>2) ? Selected_bjet_pt[2]:-100")
               .Define("Selected_bjet_sum_all_jets_pt","Sum(Selected_bjet_pt)")
               .Define("Selected_bjet_eta","Jet_eta[goodJets_btag]")
               .Define("Selected_bjet_phi","Jet_phi[goodJets_btag]")
               .Define("Selected_bjet_mass","Jet_mass[goodJets_btag]")
               .Define("Selected_bjet_number","int(Selected_bjet_pt.size())")
               .Define("Selected_bjet_btagDeepB","Jet_btagDeepB[goodJets_btag]")
               .Define("Selected_bjet_btagDeepFlavB","Jet_btagDeepFlavB[goodJets_btag]");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_bjet_genJetIdx","Jet_genJetIdx[goodJets_btag]")
                   .Define("Selected_bjet_partonFlavour","Jet_partonFlavour[goodJets_btag]")
                   .Define("Selected_bjet_hadronFlavour","Jet_hadronFlavour[goodJets_btag]");
    }
    _rlm = _rlm.Define("Selected_bjet_jetId","Jet_jetId[goodJets_btag]")  
               .Define("bjet4vecs", ::generate_4vec,{"Selected_bjet_pt","Selected_bjet_eta","Selected_bjet_phi","Selected_bjet_mass"});

    if(Region!="CR_TT_SL")
    {
        _rlm = _rlm.Define("Lepton_Jet_overlap",::CheckOverlaps,{"jet4vecs","lep4vecs"});
    }
    if(Region=="CR_TT_SL")
    {
        _rlm = _rlm.Define("Lepton_Jet_overlap",::CheckOverlaps,{"jet4vecs","lep4vecs_looseTTSL"});
    }
    _rlm = _rlm.Define("Selected_clean_jet_pt","Selected_jet_pt[Lepton_Jet_overlap]")
               .Define("Selected_clean_jet_leading_pt","(Selected_clean_jet_pt.size()>0) ? Selected_clean_jet_pt[0]:-100")
               .Define("Selected_clean_jet_subleading_pt","(Selected_clean_jet_pt.size()>1) ? Selected_clean_jet_pt[1]:-100")
               .Define("Selected_clean_jet_subsubleading_pt","(Selected_clean_jet_pt.size()>2) ? Selected_clean_jet_pt[2]:-100")
               .Define("Selected_clean_jet_Ht","Sum(Selected_clean_jet_pt)")
               .Define("Selected_clean_jet_eta","Selected_jet_eta[Lepton_Jet_overlap]")
               .Define("Selected_clean_jet_phi","Selected_jet_phi[Lepton_Jet_overlap]")
               .Define("Selected_clean_jet_mass","Selected_jet_mass[Lepton_Jet_overlap]")
               .Define("Selected_clean_jet_transverse_energy","sqrt(Selected_clean_jet_pt*Selected_clean_jet_pt + Selected_clean_jet_mass*Selected_clean_jet_mass)")
               .Define("Selected_clean_jet_number","int(Selected_clean_jet_pt.size())")
               .Define("Selected_clean_jet_btagDeepB","Selected_jet_btagDeepB[Lepton_Jet_overlap]")
               .Define("Selected_clean_jet_btagDeepFlavB","Selected_jet_btagDeepFlavB[Lepton_Jet_overlap]");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_clean_jet_genJetIdx","Selected_jet_genJetIdx[Lepton_Jet_overlap]")
                   .Define("Selected_clean_jet_partonFlavour","Selected_jet_partonFlavour[Lepton_Jet_overlap]")
                   .Define("Selected_clean_jet_hadronFlavour","Selected_jet_hadronFlavour[Lepton_Jet_overlap]");
    }
    _rlm = _rlm.Define("Selected_clean_jet_jetId","Selected_jet_jetId[Lepton_Jet_overlap]")      
               .Define("cleanjet4vecs", ::generate_4vec,{"Selected_clean_jet_pt","Selected_clean_jet_eta","Selected_clean_jet_phi","Selected_clean_jet_mass"})
               .Define("cleanjet4Rvecs",::generate_4Rvec,{"Selected_clean_jet_pt","Selected_clean_jet_eta","Selected_clean_jet_phi","Selected_clean_jet_mass"})
               .Define("Selected_clean_jet_energy","ROOT::VecOps::Map(cleanjet4Rvecs, FourVecEnergy)")
               .Define("Selected_clean_jet_px","ROOT::VecOps::Map(cleanjet4Rvecs, FourVecPx)")
               .Define("Selected_clean_jet_py","ROOT::VecOps::Map(cleanjet4Rvecs, FourVecPy)")
               .Define("Selected_clean_jet_pz","ROOT::VecOps::Map(cleanjet4Rvecs, FourVecPz)");

    if(Region!="CR_TT_SL")
    {
        _rlm = _rlm.Define("Lepton_bJet_overlap",::CheckOverlaps,{"bjet4vecs","lep4vecs"});
    }
    if(Region=="CR_TT_SL")
    {
        _rlm = _rlm.Define("Lepton_bJet_overlap",::CheckOverlaps,{"bjet4vecs","lep4vecs_looseTTSL"});
    }
    _rlm = _rlm.Define("Selected_clean_bjet_pt","Selected_bjet_pt[Lepton_bJet_overlap]")
               .Define("Selected_clean_bjet_leading_pt","(Selected_clean_bjet_pt.size()>0) ? Selected_clean_bjet_pt[0]:-100")
               .Define("Selected_clean_bjet_subleading_pt","(Selected_clean_bjet_pt.size()>1) ? Selected_clean_bjet_pt[1]:-100")
               .Define("Selected_clean_bjet_subsubleading_pt","(Selected_clean_bjet_pt.size()>2) ? Selected_clean_bjet_pt[2]:-100")
               .Define("Selected_clean_bjet_sum_all_jets_pt","Sum(Selected_clean_bjet_pt)")
               .Define("Selected_clean_bjet_eta","Selected_bjet_eta[Lepton_bJet_overlap]")
               .Define("Selected_clean_bjet_phi","Selected_bjet_phi[Lepton_bJet_overlap]")
               .Define("Selected_clean_bjet_mass","Selected_bjet_mass[Lepton_bJet_overlap]")
               .Define("Selected_clean_bjet_transverse_energy","sqrt(Selected_clean_bjet_pt*Selected_clean_bjet_pt + Selected_clean_bjet_mass*Selected_clean_bjet_mass)")
               .Define("Selected_clean_bjet_number","int(Selected_clean_bjet_pt.size())")
               .Define("Selected_clean_bjet_btagDeepB","Selected_bjet_btagDeepB[Lepton_bJet_overlap]")
               .Define("Selected_clean_bjet_btagDeepFlavB","Selected_bjet_btagDeepFlavB[Lepton_bJet_overlap]");
    if(!_isData)
    {
        _rlm = _rlm.Define("Selected_clean_bjet_genJetIdx","Selected_bjet_genJetIdx[Lepton_bJet_overlap]")
                   .Define("Selected_clean_bjet_partonFlavour","Selected_bjet_partonFlavour[Lepton_bJet_overlap]")
                   .Define("Selected_clean_bjet_hadronFlavour","Selected_bjet_hadronFlavour[Lepton_bJet_overlap]");
    }
    _rlm = _rlm.Define("Selected_clean_bjet_jetId","Selected_bjet_jetId[Lepton_bJet_overlap]")  
               .Define("cleanbjet4vecs", ::generate_4vec,{"Selected_clean_bjet_pt","Selected_clean_bjet_eta","Selected_clean_bjet_phi","Selected_clean_bjet_mass"})
               .Define("cleanbjet4Rvecs",::generate_4Rvec,{"Selected_clean_bjet_pt","Selected_clean_bjet_eta","Selected_clean_bjet_phi","Selected_clean_bjet_mass"})
               .Define("Selected_clean_bjet_energy","ROOT::VecOps::Map(cleanbjet4Rvecs, FourVecEnergy)")
               .Define("Selected_clean_bjet_px","ROOT::VecOps::Map(cleanbjet4Rvecs, FourVecPx)")
               .Define("Selected_clean_bjet_py","ROOT::VecOps::Map(cleanbjet4Rvecs, FourVecPy)")
               .Define("Selected_clean_bjet_pz","ROOT::VecOps::Map(cleanbjet4Rvecs, FourVecPz)");

    _rlm = _rlm.Define("dis_comb",::distinct_comb_3,{"Selected_clean_jet_number"})
               .Define("Selected_clean_jet1","ROOT::VecOps::Take(cleanjet4Rvecs, dis_comb[0])")
               .Define("Selected_clean_jet1_btagDeepFlavB","ROOT::VecOps::Take(Selected_clean_jet_btagDeepFlavB, dis_comb[0])")
               .Define("Selected_clean_jet2","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? ROOT::VecOps::Take(cleanjet4Rvecs, dis_comb[1]):ROOT::VecOps::Take(cleanjet4Rvecs, dis_comb[0])")
               .Define("Selected_clean_jet2_btagDeepFlavB","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? ROOT::VecOps::Take(Selected_clean_jet_btagDeepFlavB, dis_comb[1]):ROOT::VecOps::Take(Selected_clean_jet_btagDeepFlavB, dis_comb[0])")
               .Define("Selected_clean_jet3","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? ROOT::VecOps::Take(cleanjet4Rvecs, dis_comb[2]):ROOT::VecOps::Take(cleanjet4Rvecs, dis_comb[0])")
               .Define("Selected_clean_jet3_btagDeepFlavB","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? ROOT::VecOps::Take(Selected_clean_jet_btagDeepFlavB, dis_comb[2]):ROOT::VecOps::Take(Selected_clean_jet_btagDeepFlavB, dis_comb[0])")
               .Define("Selected_clean_jet_condition","Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1")
               .Define("Selected_clean_jet_condition_btagDeepFlavB","Selected_clean_jet1_btagDeepFlavB > 0.2783 || Selected_clean_jet2_btagDeepFlavB > 0.2783 || Selected_clean_jet3_btagDeepFlavB > 0.2783")
               .Define("Vectorial_sum_three_clean_jets",::Threejets,{"Selected_clean_jet1","Selected_clean_jet2","Selected_clean_jet3","Selected_clean_jet_number","Selected_clean_bjet_number"})
               .Define("Vectorial_sum_three_clean_jets_pt","ROOT::VecOps::Map(Vectorial_sum_three_clean_jets, FourVecPt)*Selected_clean_jet_condition*Selected_clean_jet_condition_btagDeepFlavB - 100*(1-Selected_clean_jet_condition) - 100*(Selected_clean_jet_condition*(Selected_clean_jet_condition-Selected_clean_jet_condition_btagDeepFlavB))")
               .Define("Vectorial_sum_three_clean_jets_eta","ROOT::VecOps::Map(Vectorial_sum_three_clean_jets, FourVecEta)*Selected_clean_jet_condition*Selected_clean_jet_condition_btagDeepFlavB - 10*(1-Selected_clean_jet_condition) - 10*(Selected_clean_jet_condition*(Selected_clean_jet_condition-Selected_clean_jet_condition_btagDeepFlavB))")
               .Define("Vectorial_sum_three_clean_jets_phi","ROOT::VecOps::Map(Vectorial_sum_three_clean_jets, FourVecPhi)*Selected_clean_jet_condition*Selected_clean_jet_condition_btagDeepFlavB - 10*(1-Selected_clean_jet_condition) - 10*(Selected_clean_jet_condition*(Selected_clean_jet_condition-Selected_clean_jet_condition_btagDeepFlavB))")
               .Define("Vectorial_sum_three_clean_jets_mass","ROOT::VecOps::Map(Vectorial_sum_three_clean_jets, FourVecMass)*Selected_clean_jet_condition*Selected_clean_jet_condition_btagDeepFlavB + 175.9*(1-Selected_clean_jet_condition) + 4000*(Selected_clean_jet_condition*(Selected_clean_jet_condition-Selected_clean_jet_condition_btagDeepFlavB))") 
               .Define("Vectorial_sum_three_clean_jets_mass_offset","Vectorial_sum_three_clean_jets_mass - 175.9")
               .Define("Vectorial_sum_three_clean_jets_mass_min","Min(abs(Vectorial_sum_three_clean_jets_mass_offset))")
               .Define("Vectorial_sum_two_clean_jets",::Twojets,{"Selected_clean_jet1","Selected_clean_jet2","Selected_clean_jet3","Selected_clean_jet1_btagDeepFlavB","Selected_clean_jet2_btagDeepFlavB","Selected_clean_jet3_btagDeepFlavB","Selected_clean_jet_number","Selected_clean_bjet_number"})
               .Define("Vectorial_sum_two_clean_jets_mass","ROOT::VecOps::Map(Vectorial_sum_two_clean_jets, FourVecMass)*Selected_clean_jet_condition*Selected_clean_jet_condition_btagDeepFlavB + 83.9*(1-Selected_clean_jet_condition) + 4000*(Selected_clean_jet_condition*(Selected_clean_jet_condition-Selected_clean_jet_condition_btagDeepFlavB))") 
               .Define("Vectorial_sum_two_clean_jets_mass_offset","Vectorial_sum_two_clean_jets_mass - 83.9")
               .Define("Vectorial_sum_two_clean_jets_mass_min","Min(abs(Vectorial_sum_two_clean_jets_mass_offset))");

}

void TprimeAnalyser::selectReconstructedJets()
{
    // _rlm = _rlm.Define("sel_partgentopHiggs",::isPfromtopHiggs,{"GenPart_pdgId","GenPart_genPartIdxMother","Selected_lepton_genPartIdx"})
    //            .Define("PartFromtopHiggs_eta","GenPart_eta[sel_partgentopHiggs]")
    //            .Define("PartFromtopHiggs_phi","GenPart_phi[sel_partgentopHiggs]");

    // _rlm = _rlm.Define("sel_jgentopHiggs",::isdR04,{"GenJet_eta","PartFromtopHiggs_eta","GenJet_phi","PartFromtopHiggs_phi"})
    //            .Define("JFromtopHiggs_pt","GenJet_pt[sel_jgentopHiggs]");

    // _rlm = _rlm.Define("sel_jetgentopHiggs",::isJetFromgJet,{"Selected_clean_jet_pt","GenJet_pt","JFromtopHiggs_pt","Selected_clean_jet_genJetIdx"})//,"Selected_clean_jet_eta","Selected_clean_jet_phi","Selected_clean_jet_mass","Selected_clean_jet_jetId","GenJet_eta","GenJet_phi","GenJet_mass"})
    //            .Define("JetFromtopHiggs_pt","Selected_clean_jet_pt[sel_jetgentopHiggs]")
    //            .Define("JetFromtopHiggs_leading_pt","(JetFromtopHiggs_pt.size()>0) ? JetFromtopHiggs_pt[0]:-100")
    //            .Define("JetFromtopHiggs_subleading_pt","(JetFromtopHiggs_pt.size()>1) ? JetFromtopHiggs_pt[1]:-100")
    //            .Define("JetFromtopHiggs_subsubleading_pt","(JetFromtopHiggs_pt.size()>2) ? JetFromtopHiggs_pt[2]:-100")          
    //            .Define("JetFromtopHiggs_Ht","Sum(JetFromtopHiggs_pt)")
    //            .Define("JetFromtopHiggs_eta","Selected_clean_jet_eta[sel_jetgentopHiggs]")
    //            .Define("JetFromtopHiggs_phi","Selected_clean_jet_phi[sel_jetgentopHiggs]")
    //            .Define("JetFromtopHiggs_mass","Selected_clean_jet_mass[sel_jetgentopHiggs]")
    //            .Define("JetFromtopHiggs_transverse_energy","sqrt(JetFromtopHiggs_pt*JetFromtopHiggs_pt + JetFromtopHiggs_mass*JetFromtopHiggs_mass)")
    //            .Define("JetFromtopHiggs_number","int(JetFromtopHiggs_pt.size())")
    //            .Define("JetFromtopHiggs_btagDeepB","Selected_clean_jet_btagDeepB[sel_jetgentopHiggs]")
    //            .Define("JetFromtopHiggs_btagDeepFlavB","Selected_clean_jet_btagDeepFlavB[sel_jetgentopHiggs]")
    //            .Define("JetFromtopHiggs_genJetIdx","Selected_clean_jet_genJetIdx[sel_jetgentopHiggs]")
    //            .Define("JetFromtopHiggs_partonFlavour","Selected_clean_jet_partonFlavour[sel_jetgentopHiggs]")
    //            .Define("JetFromtopHiggs_hadronFlavour","Selected_clean_jet_hadronFlavour[sel_jetgentopHiggs]")
    //            .Define("p4_JetFromtopHiggs",::generate_4vec,{"JetFromtopHiggs_pt","JetFromtopHiggs_eta","JetFromtopHiggs_phi","JetFromtopHiggs_mass"})
    //            .Define("Rp4_JetFromtopHiggs",::generate_4Rvec,{"JetFromtopHiggs_pt","JetFromtopHiggs_eta","JetFromtopHiggs_phi","JetFromtopHiggs_mass"})
    //            .Define("JetFromtopHiggs_energy","ROOT::VecOps::Map(Rp4_JetFromtopHiggs, FourVecEnergy)")
    //            .Define("JetFromtopHiggs_px","ROOT::VecOps::Map(Rp4_JetFromtopHiggs, FourVecPx)")
    //            .Define("JetFromtopHiggs_py","ROOT::VecOps::Map(Rp4_JetFromtopHiggs, FourVecPy)")
    //            .Define("JetFromtopHiggs_pz","ROOT::VecOps::Map(Rp4_JetFromtopHiggs, FourVecPz)");

    // _rlm = _rlm.Define("sel_bjetgentopHiggs",::isJetFromgJet,{"Selected_clean_bjet_pt","GenJet_pt","JFromtopHiggs_pt","Selected_clean_bjet_genJetIdx"})//,"Selected_clean_bjet_eta","Selected_clean_bjet_phi","Selected_clean_bjet_mass","Selected_clean_bjet_jetId","GenJet_eta","GenJet_phi","GenJet_mass"})
    //            .Define("bJetFromtopHiggs_pt","Selected_clean_bjet_pt[sel_bjetgentopHiggs]")
    //            .Define("bJetFromtopHiggs_leading_pt","(bJetFromtopHiggs_pt.size()>0) ? bJetFromtopHiggs_pt[0]:-100")
    //            .Define("bJetFromtopHiggs_subleading_pt","(bJetFromtopHiggs_pt.size()>1) ? bJetFromtopHiggs_pt[1]:-100")
    //            .Define("bJetFromtopHiggs_subsubleading_pt","(bJetFromtopHiggs_pt.size()>2) ? bJetFromtopHiggs_pt[2]:-100")
    //            .Define("bJetFromtopHiggs_Ht","Sum(bJetFromtopHiggs_pt)")
    //            .Define("bJetFromtopHiggs_eta","Selected_clean_bjet_eta[sel_bjetgentopHiggs]")
    //            .Define("bJetFromtopHiggs_phi","Selected_clean_bjet_phi[sel_bjetgentopHiggs]")
    //            .Define("bJetFromtopHiggs_mass","Selected_clean_bjet_mass[sel_bjetgentopHiggs]")
    //            .Define("bJetFromtopHiggs_transverse_energy","sqrt(bJetFromtopHiggs_pt*bJetFromtopHiggs_pt + bJetFromtopHiggs_mass*bJetFromtopHiggs_mass)")
    //            .Define("bJetFromtopHiggs_number","int(bJetFromtopHiggs_pt.size())")
    //            .Define("bJetFromtopHiggs_btagDeepB","Selected_clean_bjet_btagDeepB[sel_bjetgentopHiggs]")
    //            .Define("bJetFromtopHiggs_btagDeepFlavB","Selected_clean_bjet_btagDeepFlavB[sel_bjetgentopHiggs]")
    //            .Define("bJetFromtopHiggs_genJetIdx","Selected_clean_bjet_genJetIdx[sel_bjetgentopHiggs]")
    //            .Define("bJetFromtopHiggs_partonFlavour","Selected_clean_bjet_partonFlavour[sel_bjetgentopHiggs]")
    //            .Define("bJetFromtopHiggs_hadronFlavour","Selected_clean_bjet_hadronFlavour[sel_bjetgentopHiggs]")
    //            .Define("p4_bJetFromtopHiggs",::generate_4vec,{"bJetFromtopHiggs_pt","bJetFromtopHiggs_eta","bJetFromtopHiggs_phi","bJetFromtopHiggs_mass"})
    //            .Define("Rp4_bJetFromtopHiggs",::generate_4Rvec,{"bJetFromtopHiggs_pt","bJetFromtopHiggs_eta","bJetFromtopHiggs_phi","bJetFromtopHiggs_mass"})
    //            .Define("bJetFromtopHiggs_energy","ROOT::VecOps::Map(Rp4_bJetFromtopHiggs, FourVecEnergy)")
    //            .Define("bJetFromtopHiggs_px","ROOT::VecOps::Map(Rp4_bJetFromtopHiggs, FourVecPx)")
    //            .Define("bJetFromtopHiggs_py","ROOT::VecOps::Map(Rp4_bJetFromtopHiggs, FourVecPy)")
    //            .Define("bJetFromtopHiggs_pz","ROOT::VecOps::Map(Rp4_bJetFromtopHiggs, FourVecPz)");

    _rlm = _rlm.Define("sel_partgenHiggs",::isPfromHiggs,{"GenPart_pdgId","GenPart_genPartIdxMother","Selected_lepton_genPartIdx"})
               .Define("PartFromHiggs_eta","GenPart_eta[sel_partgenHiggs]")
               .Define("PartFromHiggs_phi","GenPart_phi[sel_partgenHiggs]");

    _rlm = _rlm.Define("sel_jgenHiggs",::isdR04,{"GenJet_eta","PartFromHiggs_eta","GenJet_phi","PartFromHiggs_phi"})
               .Define("JFromHiggs_pt","GenJet_pt[sel_jgenHiggs]");

    _rlm = _rlm.Define("sel_jetgenHiggs",::isJetFromgJet,{"Selected_clean_jet_pt","GenJet_pt","JFromHiggs_pt","Selected_clean_jet_genJetIdx"})//,"Selected_clean_jet_eta","Selected_clean_jet_phi","Selected_clean_jet_mass","Selected_clean_jet_jetId","GenJet_eta","GenJet_phi","GenJet_mass"})
               .Define("JetFromHiggs_pt","Selected_clean_jet_pt[sel_jetgenHiggs]")
               .Define("JetFromHiggs_leading_pt","(JetFromHiggs_pt.size()>0) ? JetFromHiggs_pt[0]:-100")
               .Define("JetFromHiggs_subleading_pt","(JetFromHiggs_pt.size()>1) ? JetFromHiggs_pt[1]:-100")
               .Define("JetFromHiggs_subsubleading_pt","(JetFromHiggs_pt.size()>2) ? JetFromHiggs_pt[2]:-100")
               .Define("JetFromHiggs_Ht","Sum(JetFromHiggs_pt)")
               .Define("JetFromHiggs_eta","Selected_clean_jet_eta[sel_jetgenHiggs]")
               .Define("JetFromHiggs_phi","Selected_clean_jet_phi[sel_jetgenHiggs]")
               .Define("JetFromHiggs_mass","Selected_clean_jet_mass[sel_jetgenHiggs]")
               .Define("JetFromHiggs_transverse_energy","sqrt(JetFromHiggs_pt*JetFromHiggs_pt + JetFromHiggs_mass*JetFromHiggs_mass)")
               .Define("JetFromHiggs_number","int(JetFromHiggs_pt.size())")
               .Define("JetFromHiggs_btagDeepB","Selected_clean_jet_btagDeepB[sel_jetgenHiggs]")
               .Define("JetFromHiggs_btagDeepFlavB","Selected_clean_jet_btagDeepFlavB[sel_jetgenHiggs]")
               .Define("JetFromHiggs_genJetIdx","Selected_clean_jet_genJetIdx[sel_jetgenHiggs]")
               .Define("JetFromHiggs_partonFlavour","Selected_clean_jet_partonFlavour[sel_jetgenHiggs]")
               .Define("JetFromHiggs_hadronFlavour","Selected_clean_jet_hadronFlavour[sel_jetgenHiggs]")
               .Define("p4_JetFromHiggs",::generate_4vec,{"JetFromHiggs_pt","JetFromHiggs_eta","JetFromHiggs_phi","JetFromHiggs_mass"})
               .Define("Rp4_JetFromHiggs",::generate_4Rvec,{"JetFromHiggs_pt","JetFromHiggs_eta","JetFromHiggs_phi","JetFromHiggs_mass"})
               .Define("JetFromHiggs_energy","ROOT::VecOps::Map(Rp4_JetFromHiggs, FourVecEnergy)")
               .Define("JetFromHiggs_px","ROOT::VecOps::Map(Rp4_JetFromHiggs, FourVecPx)")
               .Define("JetFromHiggs_py","ROOT::VecOps::Map(Rp4_JetFromHiggs, FourVecPy)")
               .Define("JetFromHiggs_pz","ROOT::VecOps::Map(Rp4_JetFromHiggs, FourVecPz)");

    _rlm = _rlm.Define("sel_partgentop",::isPfromtop,{"GenPart_pdgId","GenPart_genPartIdxMother","Selected_lepton_genPartIdx"})
               .Define("PartFromtop_eta","GenPart_eta[sel_partgentop]")
               .Define("PartFromtop_phi","GenPart_phi[sel_partgentop]");

    _rlm = _rlm.Define("sel_jgentop",::isdR04,{"GenJet_eta","PartFromtop_eta","GenJet_phi","PartFromtop_phi"})
               .Define("JFromtop_pt","GenJet_pt[sel_jgentop]");

    _rlm = _rlm.Define("sel_bjetgentop",::isJetFromgJet,{"Selected_clean_bjet_pt","GenJet_pt","JFromtop_pt","Selected_clean_bjet_genJetIdx"})//,"Selected_clean_bjet_eta","Selected_clean_bjet_phi","Selected_clean_bjet_mass","Selected_clean_bjet_jetId","GenJet_eta","GenJet_phi","GenJet_mass"})
               .Define("bJetFromtop_pt","Selected_clean_bjet_pt[sel_bjetgentop]")
               .Define("bJetFromtop_leading_pt","(bJetFromtop_pt.size()>0) ? bJetFromtop_pt[0]:-100")
               .Define("bJetFromtop_subleading_pt","(bJetFromtop_pt.size()>1) ? bJetFromtop_pt[1]:-100")
               .Define("bJetFromtop_subsubleading_pt","(bJetFromtop_pt.size()>2) ? bJetFromtop_pt[2]:-100")
               .Define("bJetFromtop_Ht","Sum(bJetFromtop_pt)")
               .Define("bJetFromtop_eta","Selected_clean_bjet_eta[sel_bjetgentop]")
               .Define("bJetFromtop_phi","Selected_clean_bjet_phi[sel_bjetgentop]")
               .Define("bJetFromtop_mass","Selected_clean_bjet_mass[sel_bjetgentop]")
               .Define("bJetFromtop_transverse_energy","sqrt(bJetFromtop_pt*bJetFromtop_pt + bJetFromtop_mass*bJetFromtop_mass)")
               .Define("bJetFromtop_number","int(bJetFromtop_pt.size())")
               .Define("bJetFromtop_btagDeepB","Selected_clean_bjet_btagDeepB[sel_bjetgentop]")
               .Define("bJetFromtop_btagDeepFlavB","Selected_clean_bjet_btagDeepFlavB[sel_bjetgentop]")
               .Define("bJetFromtop_genJetIdx","Selected_clean_bjet_genJetIdx[sel_bjetgentop]")
               .Define("bJetFromtop_partonFlavour","Selected_clean_bjet_partonFlavour[sel_bjetgentop]")
               .Define("bJetFromtop_hadronFlavour","Selected_clean_bjet_hadronFlavour[sel_bjetgentop]")
               .Define("p4_bJetFromtop",::generate_4vec,{"bJetFromtop_pt","bJetFromtop_eta","bJetFromtop_phi","bJetFromtop_mass"})
               .Define("Rp4_bJetFromtop",::generate_4Rvec,{"bJetFromtop_pt","bJetFromtop_eta","bJetFromtop_phi","bJetFromtop_mass"})
            //    .Define("bJetFromtop_energy","ROOT::VecOps::Map(Rp4_bJetFromtop, FourVecEnergy)")
               .Define("bJetFromtop_px","ROOT::VecOps::Map(Rp4_bJetFromtop, FourVecPx)")
               .Define("bJetFromtop_py","ROOT::VecOps::Map(Rp4_bJetFromtop, FourVecPy)")
               .Define("bJetFromtop_pz","ROOT::VecOps::Map(Rp4_bJetFromtop, FourVecPz)");

    // _rlm = _rlm.Define("sel_partgentop",::isPfromtoponly,{"GenPart_pdgId","GenPart_genPartIdxMother","Selected_lepton_genPartIdx"})
    //            .Define("PartFromtop_eta","GenPart_eta[sel_partgentop]")
    //            .Define("PartFromtop_phi","GenPart_phi[sel_partgentop]");

    // _rlm = _rlm.Define("sel_jgentop",::isdR04,{"GenJet_eta","PartFromtop_eta","GenJet_phi","PartFromtop_phi"})
    //            .Define("JFromtop_pt","GenJet_pt[sel_jgentop]");

    // _rlm = _rlm.Define("sel_bjetgentop",::isJetFromgJet,{"Selected_clean_bjet_pt","GenJet_pt","JFromtop_pt","Selected_clean_bjet_genJetIdx"})//,"Selected_clean_bjet_eta","Selected_clean_bjet_phi","Selected_clean_bjet_mass","Selected_clean_bjet_jetId","GenJet_eta","GenJet_phi","GenJet_mass"})
    //            .Define("bJetFromtop_pt","Selected_clean_bjet_pt[sel_bjetgentop]")
    //            .Define("bJetFromtop_leading_pt","(bJetFromtop_pt.size()>0) ? bJetFromtop_pt[0]:-100")
    //            .Define("bJetFromtop_subleading_pt","(bJetFromtop_pt.size()>1) ? bJetFromtop_pt[1]:-100")
    //            .Define("bJetFromtop_subsubleading_pt","(bJetFromtop_pt.size()>2) ? bJetFromtop_pt[2]:-100")
    //            .Define("bJetFromtop_Ht","Sum(bJetFromtop_pt)")
    //            .Define("bJetFromtop_eta","Selected_clean_bjet_eta[sel_bjetgentop]")
    //            .Define("bJetFromtop_phi","Selected_clean_bjet_phi[sel_bjetgentop]")
    //            .Define("bJetFromtop_mass","Selected_clean_bjet_mass[sel_bjetgentop]")
    //            .Define("bJetFromtop_transverse_energy","sqrt(bJetFromtop_pt*bJetFromtop_pt + bJetFromtop_mass*bJetFromtop_mass)")
    //            .Define("bJetFromtop_number","int(bJetFromtop_pt.size())")
    //            .Define("bJetFromtop_btagDeepB","Selected_clean_bjet_btagDeepB[sel_bjetgentop]")
    //            .Define("bJetFromtop_btagDeepFlavB","Selected_clean_bjet_btagDeepFlavB[sel_bjetgentop]")
    //            .Define("bJetFromtop_genJetIdx","Selected_clean_bjet_genJetIdx[sel_bjetgentop]")
    //            .Define("bJetFromtop_partonFlavour","Selected_clean_bjet_partonFlavour[sel_bjetgentop]")
    //            .Define("bJetFromtop_hadronFlavour","Selected_clean_bjet_hadronFlavour[sel_bjetgentop]")
    //            .Define("p4_bJetFromtop",::generate_4vec,{"bJetFromtop_pt","bJetFromtop_eta","bJetFromtop_phi","bJetFromtop_mass"})
    //            .Define("p4_bJetFromtop_first","p4_bJetFromtop[0]")
    //            .Define("bJetFromtop_first_pt","p4_bJetFromtop_first.Pt()")
    //            .Define("bJetFromtop_first_eta","p4_bJetFromtop_first.Eta()")
    //            .Define("bJetFromtop_first_phi","p4_bJetFromtop_first.Phi()")
    //            .Define("bJetFromtop_first_mass","p4_bJetFromtop_first.M()")
    //            .Define("bJetFromtop_first_energy","p4_bJetFromtop_first.E()")
    //            .Define("bJetFromtop_first_px","p4_bJetFromtop_first.Px()")
    //            .Define("bJetFromtop_first_py","p4_bJetFromtop_first.Py()")
    //            .Define("bJetFromtop_first_pz","p4_bJetFromtop_first.Pz()")
    //            .Define("Rp4_bJetFromtop_first",::generate_4Rvec_2,{"bJetFromtop_first_pt","bJetFromtop_first_eta","bJetFromtop_first_phi","bJetFromtop_first_mass"});

    // _rlm = _rlm.Define("GenJets_pt","GenJet_pt")
    //            .Define("GenJets_eta","GenJet_eta")
    //            .Define("GenJets_phi","GenJet_phi")
    //            .Define("GenJets_mass","GenJet_mass")
    //            .Define("GenJets_hadronFlavour","GenJet_hadronFlavour")
    //            .Define("GenJets_partonFlavour","GenJet_partonFlavour")
    //            .Define("GenJets_light","abs(GenJets_partonFlavour) < 5")
    //            .Define("GenJets_light_pt","GenJet_pt[GenJets_light]")
    //            .Define("GenJets_light_eta","GenJet_eta[GenJets_light]")
    //            .Define("GenJets_light_phi","GenJet_phi[GenJets_light]")
    //            .Define("GenJets_light_mass","GenJet_mass[GenJets_light]")
    //            .Define("GenJets_light_hadronFlavour","GenJet_hadronFlavour[GenJets_light]")
    //            .Define("GenJets_light_partonFlavour","GenJet_partonFlavour[GenJets_light]");

}

void TprimeAnalyser::calculateEvWeight(std::string Region, std::string year)
{		
  int _case = 1;
  std::vector<std::string> Jets_vars_names = {"Selected_clean_bjet_hadronFlavour", "Selected_clean_bjet_eta", "Selected_clean_bjet_pt"};  
  if(_case !=1)
  {
    Jets_vars_names.emplace_back("Selected_clean_bjet_btagDeepFlavB");
  }
  std::string output_btag_column_name = "btag_SF_";
  _rlm = calculateBTagSF(_rlm, Jets_vars_names, _case, output_btag_column_name);

  std::vector<std::string> Muon_vars_names;
  if(Region!="CR_TT_SL")
  {
    Muon_vars_names = {"Selected_muon_eta", "Selected_muon_pt"};
  }
  else
  {
    Muon_vars_names = {"Selected_muon_looseTTSL_eta", "Selected_muon_looseTTSL_pt"};
  }
  std::string output_mu_column_name = "muon_SF_";
  _rlm = calculateMuSF(_rlm, Muon_vars_names, output_mu_column_name);

  std::vector<std::string> Electron_vars_names;
  if(Region!="CR_TT_SL")
  {
    Electron_vars_names = {"Selected_electron_eta", "Selected_electron_pt"};
  }
  else
  {
    Electron_vars_names = {"Selected_electron_looseTTSL_eta", "Selected_electron_looseTTSL_pt"};
  }
  std::string output_ele_column_name = "ele_SF_";
  _rlm = calculateEleSF(_rlm, Electron_vars_names, output_ele_column_name);

  _rlm = applyPrefiringWeight(_rlm);
  //   Total event Weight:
  //_rlm = _rlm.Define("evWeight", " pugenWeight * btag_SF_central * muon_SF_central * ele_SF_central"); 
  _rlm = _rlm.Define("evWeight", " pugenWeight * prefiring_SF_central * btag_SF_central * muon_SF_central * ele_SF_central"); 
  // _rlm = _rlm.Define("evWeight", " pugenWeight * btag_SF_central "); 

}

//=============================define variables==================================================//
void TprimeAnalyser::defineMoreVars(std::string Region, std::string year)
{
    if (debug)
    {
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

    //combination of particles
    if(Region=="CR_TT_SL")
    {
        _rlm = _rlm.Define("St","Sum(Selected_electron_looseTTSL_pt) + Sum(Selected_electron_pt) + Sum(Selected_muon_looseTTSL_pt) + Sum(Selected_clean_jet_pt)")
                   .Define("St_with_MET","Sum(Selected_electron_looseTTSL_pt) + Sum(Selected_electron_pt) + Sum(Selected_muon_looseTTSL_pt) + Sum(Selected_clean_jet_pt) + MET_pt");
    }
    if(Region!="CR_TT_SL")
    {
        _rlm = _rlm.Define("St","Sum(Selected_electron_pt) + Sum(Selected_muon_pt) + Sum(Selected_clean_jet_pt)")
                   .Define("St_with_MET","Sum(Selected_electron_pt) + Sum(Selected_muon_pt) + Sum(Selected_clean_jet_pt) + MET_pt");
    }

    // _rlm = _rlm.Define("Selected_clean_jet_btagDeepFlavB_condition","Selected_clean_jet_btagDeepFlavB[0] > 0.2783 || Selected_clean_jet_btagDeepFlavB[1] > 0.2783 || Selected_clean_jet_btagDeepFlavB[2] > 0.2783")
    //            .Define("Selected_clean_third_jet","(Selected_clean_jet_btagDeepFlavB_condition == 1) ? cleanjet4vecs[2] : cleanbjet4vecs[0]")
    //            .Define("sumparticles4vecs","lep4vecs[0] + lep4vecs[1] + cleanjet4vecs[0] + cleanjet4vecs[1] + Selected_clean_third_jet")
    //            .Define("Vectorial_sum_particles_pt","sumparticles4vecs.Pt()")
    //            .Define("Vectorial_sum_particles_eta","sumparticles4vecs.Eta()")
    //            .Define("Vectorial_sum_particles_phi","sumparticles4vecs.Phi()")
    //            .Define("Vectorial_sum_particles_mass","sumparticles4vecs.M()");

    // _rlm = _rlm.Define("deltaR_Selected_leading_lepton_Selected_clean_leading_jet","ROOT::VecOps::DeltaR(Selected_lepton_eta[0],Selected_clean_jet_eta[0],Selected_lepton_phi[0],Selected_clean_jet_eta[0])")
    //            .Define("deltaR_Selected_leading_lepton_Selected_clean_subleading_jet","ROOT::VecOps::DeltaR(Selected_lepton_eta[0],Selected_clean_jet_eta[1],Selected_lepton_phi[0],Selected_clean_jet_eta[1])")
    //            .Define("deltaR_Selected_leading_lepton_Selected_clean_subsubleading_jet","ROOT::VecOps::DeltaR(Selected_lepton_eta[0],Selected_clean_jet_eta[2],Selected_lepton_phi[0],Selected_clean_jet_eta[2])")
    //            .Define("deltaR_Selected_leading_lepton_Selected_clean_leading_bjet","ROOT::VecOps::DeltaR(Selected_lepton_eta[0],Selected_clean_bjet_eta[0],Selected_lepton_phi[0],Selected_clean_bjet_eta[0])")
    //            .Define("deltaR_Selected_subleading_lepton_Selected_clean_leading_jet","ROOT::VecOps::DeltaR(Selected_lepton_eta[1],Selected_clean_jet_eta[0],Selected_lepton_phi[1],Selected_clean_jet_eta[0])")
    //            .Define("deltaR_Selected_subleading_lepton_Selected_clean_subleading_jet","ROOT::VecOps::DeltaR(Selected_lepton_eta[1],Selected_clean_jet_eta[1],Selected_lepton_phi[1],Selected_clean_jet_eta[1])")
    //            .Define("deltaR_Selected_subleading_lepton_Selected_clean_subsubleading_jet","ROOT::VecOps::DeltaR(Selected_lepton_eta[1],Selected_clean_jet_eta[2],Selected_lepton_phi[1],Selected_clean_jet_eta[2])")
    //            .Define("deltaR_Selected_subleading_lepton_Selected_clean_leading_bjet","ROOT::VecOps::DeltaR(Selected_lepton_eta[1],Selected_clean_bjet_eta[0],Selected_lepton_phi[1],Selected_clean_bjet_eta[0])");

    // _rlm = _rlm.Define("Selected_leading_lepton_Pt_over_Selected_clean_leading_jet_Pt","Selected_lepton_leading_pt/Selected_clean_jet_leading_pt")
    //            .Define("Selected_leading_lepton_Pt_over_Selected_clean_subleading_jet_Pt","Selected_lepton_leading_pt/Selected_clean_jet_subleading_pt")
    //            .Define("Selected_leading_lepton_Pt_over_Selected_clean_subsubleading_jet_Pt","Selected_lepton_leading_pt/Selected_clean_jet_subsubleading_pt")
    //            .Define("Selected_leading_lepton_Pt_over_Selected_clean_leading_bjet_Pt","Selected_lepton_leading_pt/Selected_clean_bjet_leading_pt")
    //            .Define("Selected_subleading_lepton_Pt_over_Selected_clean_leading_jet_Pt","Selected_lepton_subleading_pt/Selected_clean_jet_leading_pt")
    //            .Define("Selected_subleading_lepton_Pt_over_Selected_clean_subleading_jet_Pt","Selected_lepton_subleading_pt/Selected_clean_jet_subleading_pt")
    //            .Define("Selected_subleading_lepton_Pt_over_Selected_clean_subsubleading_jet_Pt","Selected_lepton_subleading_pt/Selected_clean_jet_subsubleading_pt")
    //            .Define("Selected_subleading_lepton_Pt_over_Selected_clean_leading_bjet_Pt","Selected_lepton_subleading_pt/Selected_clean_bjet_leading_pt");

    // _rlm = _rlm.Define("deltaR_Selected_lepton_Selected_clean_jet",::DR_1,{"Selected_lepton_eta","Selected_clean_jet_eta","Selected_lepton_phi","Selected_clean_jet_phi"})
    //            .Define("deltaR_Selected_lepton_Selected_clean_bjet",::DR_1,{"Selected_lepton_eta","Selected_clean_bjet_eta","Selected_lepton_phi","Selected_clean_bjet_phi"});
    //            .Define("deltaR_Selected_muon_Selected_clean_jet",::DR_1,{"Selected_muon_eta","Selected_clean_jet_eta","Selected_muon_phi","Selected_clean_jet_phi"})
    //            .Define("deltaR_Selected_muon_Selected_clean_bjet",::DR_1,{"Selected_muon_eta","Selected_clean_bjet_eta","Selected_muon_phi","Selected_clean_bjet_phi"});

    // _rlm = _rlm.Define("Pt_lepton_over_Pt_clean_jet",::ratio_1,{"Selected_lepton_pt","Selected_clean_jet_pt"})
    //            .Define("Pt_lepton_over_Pt_clean_bjet",::ratio_1,{"Selected_lepton_pt","Selected_clean_bjet_pt"});
    //            .Define("Pt_muon_over_Pt_clean_jet",::ratio_1,{"Selected_muon_pt","Selected_clean_jet_pt"})
    //            .Define("Pt_muon_over_Pt_clean_bjet",::ratio_1,{"Selected_muon_pt","Selected_clean_bjet_pt"});

    // _rlm = _rlm.Define("Min_deltaR_Selected_lepton_Selected_clean_jet",::minDR_5,{"Selected_lepton_eta","Selected_clean_jet_eta","Selected_lepton_phi","Selected_clean_jet_phi"})
    //            .Define("Min_deltaR_Selected_lepton_Selected_clean_bjet",::minDR_5,{"Selected_lepton_eta","Selected_clean_bjet_eta","Selected_lepton_phi","Selected_clean_bjet_phi"});
    //            .Define("Min_deltaR_Selected_muon_Selected_clean_jet",::minDR_5,{"Selected_muon_eta","Selected_clean_jet_eta","Selected_muon_phi","Selected_clean_jet_phi"})
    //            .Define("Min_deltaR_Selected_muon_Selected_clean_bjet",::minDR_5,{"Selected_muon_eta","Selected_clean_bjet_eta","Selected_muon_phi","Selected_clean_bjet_phi"});

    // _rlm = _rlm.Define("Min_Pt_lepton_over_Pt_clean_jet",::minratio_1,{"Selected_lepton_pt","Selected_clean_jet_pt"})
    //            .Define("Min_Pt_lepton_over_Pt_clean_bjet",::minratio_1,{"Selected_lepton_pt","Selected_clean_bjet_pt"});
    //            .Define("Min_Pt_muon_over_Pt_clean_jet",::minratio_1,{"Selected_muon_pt","Selected_clean_jet_pt"})
    //            .Define("Min_Pt_muon_over_Pt_clean_bjet",::minratio_1,{"Selected_muon_pt","Selected_clean_bjet_pt"});

    if(Region=="CR_TT_SL")
    {
        _rlm = _rlm.Define("dis_comb_leptonbjet",::distinct_comb_2,{"Selected_lepton_looseTTSL_number","Selected_clean_bjet_number"})
                   .Define("Selected_lepton1","ROOT::VecOps::Take(lep4Rvecs_looseTTSL, dis_comb_leptonbjet[0])")
                   .Define("Selected_lepton1_pt","ROOT::VecOps::Take(Selected_lepton_looseTTSL_pt, dis_comb_leptonbjet[0])");
    }
    if(Region!="CR_TT_SL")
    {
        _rlm = _rlm.Define("dis_comb_leptonbjet",::distinct_comb_2,{"Selected_lepton_number","Selected_clean_bjet_number"})
                   .Define("Selected_lepton1","ROOT::VecOps::Take(lep4Rvecs, dis_comb_leptonbjet[0])")
                   .Define("Selected_lepton1_pt","ROOT::VecOps::Take(Selected_lepton_pt, dis_comb_leptonbjet[0])");
    }
    _rlm = _rlm.Define("Selected_clean_bjet1","ROOT::VecOps::Take(cleanbjet4Rvecs, dis_comb_leptonbjet[1])")
               .Define("Selected_clean_bjet1_pt","ROOT::VecOps::Take(Selected_clean_bjet_pt, dis_comb_leptonbjet[1])")
               .Define("LeptonbjetFromtop4Rvecs",::sum_two_4Rvecs,{"Selected_lepton1","Selected_clean_bjet1"})
               .Define("LeptonbjetFromtop_pt","ROOT::VecOps::Map(LeptonbjetFromtop4Rvecs, FourVecPt)")
               .Define("LeptonbjetFromtop_eta","ROOT::VecOps::Map(LeptonbjetFromtop4Rvecs, FourVecEta)")
               .Define("LeptonbjetFromtop_phi","ROOT::VecOps::Map(LeptonbjetFromtop4Rvecs, FourVecPhi)")
               .Define("LeptonbjetFromtop_mass","ROOT::VecOps::Map(LeptonbjetFromtop4Rvecs, FourVecMass)")
               .Define("LeptonbjetFromtop_energy","ROOT::VecOps::Map(LeptonbjetFromtop4Rvecs, FourVecEnergy)")
               .Define("LeptonbjetFromtop_transverse_energy","sqrt(LeptonbjetFromtop_pt*LeptonbjetFromtop_pt + LeptonbjetFromtop_mass*LeptonbjetFromtop_mass)");
    if(Region=="CR_TT_SL")
    {
        _rlm = _rlm.Define("LeptonbjetFromtop_deltaeta","abs(ROOT::VecOps::Take(Selected_lepton_looseTTSL_eta, dis_comb_leptonbjet[0]) - ROOT::VecOps::Take(Selected_clean_bjet_eta, dis_comb_leptonbjet[1]))")
                   .Define("LeptonbjetFromtop_deltaphi","abs(ROOT::VecOps::DeltaPhi(ROOT::VecOps::Take(Selected_lepton_looseTTSL_phi, dis_comb_leptonbjet[0]),ROOT::VecOps::Take(Selected_clean_bjet_phi, dis_comb_leptonbjet[1])))")
                   .Define("LeptonbjetFromtop_deltaR","ROOT::VecOps::DeltaR(ROOT::VecOps::Take(Selected_lepton_looseTTSL_eta, dis_comb_leptonbjet[0]),ROOT::VecOps::Take(Selected_clean_bjet_eta, dis_comb_leptonbjet[1]),ROOT::VecOps::Take(Selected_lepton_looseTTSL_phi, dis_comb_leptonbjet[0]),ROOT::VecOps::Take(Selected_clean_bjet_phi, dis_comb_leptonbjet[1]))")
                   .Define("LeptonbjetFromtop_transverse_mass","sqrt(ROOT::VecOps::Take(Selected_lepton_looseTTSL_mass, dis_comb_leptonbjet[0])*ROOT::VecOps::Take(Selected_lepton_looseTTSL_mass, dis_comb_leptonbjet[0]) + ROOT::VecOps::Take(Selected_clean_bjet_mass, dis_comb_leptonbjet[1])*ROOT::VecOps::Take(Selected_clean_bjet_mass, dis_comb_leptonbjet[1]) + 2*ROOT::VecOps::Take(Selected_lepton_looseTTSL_transverse_energy, dis_comb_leptonbjet[0])*ROOT::VecOps::Take(Selected_clean_bjet_transverse_energy, dis_comb_leptonbjet[1]) - 2*ROOT::VecOps::Take(Selected_lepton_looseTTSL_pt, dis_comb_leptonbjet[0])*ROOT::VecOps::Take(Selected_clean_bjet_pt, dis_comb_leptonbjet[1])*cos(LeptonbjetFromtop_deltaphi))");
    }
    if(Region!="CR_TT_SL")
    {
        _rlm = _rlm.Define("LeptonbjetFromtop_deltaeta","abs(ROOT::VecOps::Take(Selected_lepton_eta, dis_comb_leptonbjet[0]) - ROOT::VecOps::Take(Selected_clean_bjet_eta, dis_comb_leptonbjet[1]))")
                   .Define("LeptonbjetFromtop_deltaphi","abs(ROOT::VecOps::DeltaPhi(ROOT::VecOps::Take(Selected_lepton_phi, dis_comb_leptonbjet[0]),ROOT::VecOps::Take(Selected_clean_bjet_phi, dis_comb_leptonbjet[1])))")
                   .Define("LeptonbjetFromtop_deltaR","ROOT::VecOps::DeltaR(ROOT::VecOps::Take(Selected_lepton_eta, dis_comb_leptonbjet[0]),ROOT::VecOps::Take(Selected_clean_bjet_eta, dis_comb_leptonbjet[1]),ROOT::VecOps::Take(Selected_lepton_phi, dis_comb_leptonbjet[0]),ROOT::VecOps::Take(Selected_clean_bjet_phi, dis_comb_leptonbjet[1]))")
                   .Define("LeptonbjetFromtop_transverse_mass","sqrt(ROOT::VecOps::Take(Selected_lepton_mass, dis_comb_leptonbjet[0])*ROOT::VecOps::Take(Selected_lepton_mass, dis_comb_leptonbjet[0]) + ROOT::VecOps::Take(Selected_clean_bjet_mass, dis_comb_leptonbjet[1])*ROOT::VecOps::Take(Selected_clean_bjet_mass, dis_comb_leptonbjet[1]) + 2*ROOT::VecOps::Take(Selected_lepton_transverse_energy, dis_comb_leptonbjet[0])*ROOT::VecOps::Take(Selected_clean_bjet_transverse_energy, dis_comb_leptonbjet[1]) - 2*ROOT::VecOps::Take(Selected_lepton_pt, dis_comb_leptonbjet[0])*ROOT::VecOps::Take(Selected_clean_bjet_pt, dis_comb_leptonbjet[1])*cos(LeptonbjetFromtop_deltaphi))");
    }

    if(Region=="CR_TT_SL")
    {
        _rlm = _rlm.Define("Min_deltaR_bJetFromtop_leptonFromtop",::minDR_1,{"Selected_lepton_looseTTSL_eta","Selected_clean_bjet_eta","Selected_lepton_looseTTSL_phi","Selected_clean_bjet_phi"});
    }
    if(Region!="CR_TT_SL")
    {
        _rlm = _rlm.Define("Min_deltaR_bJetFromtop_leptonFromtop",::minDR_1,{"Selected_lepton_eta","Selected_clean_bjet_eta","Selected_lepton_phi","Selected_clean_bjet_phi"});
    }
    _rlm = _rlm.Define("Lepton2Fromtop","lep4vecs[Min_deltaR_bJetFromtop_leptonFromtop[0]]")
               .Define("Lepton2Fromtop_pt","Lepton2Fromtop.Pt()")
               .Define("Lepton2Fromtop_eta","Lepton2Fromtop.Eta()")
               .Define("Lepton2Fromtop_phi","Lepton2Fromtop.Phi()")
               .Define("Lepton2Fromtop_mass","Lepton2Fromtop.M()")
               .Define("Lepton2Fromtop_transverse_energy","sqrt(Lepton2Fromtop_pt*Lepton2Fromtop_pt + Lepton2Fromtop_mass*Lepton2Fromtop_mass)")
               .Define("bJet1Fromtop","bjet4vecs[Min_deltaR_bJetFromtop_leptonFromtop[1]]")
               .Define("bJet1Fromtop_pt","bJet1Fromtop.Pt()")
               .Define("bJet1Fromtop_eta","bJet1Fromtop.Eta()")
               .Define("bJet1Fromtop_phi","bJet1Fromtop.Phi()")
               .Define("bJet1Fromtop_mass","bJet1Fromtop.M()")
               .Define("bJet1Fromtop_transverse_energy","sqrt(bJet1Fromtop_pt*bJet1Fromtop_pt + bJet1Fromtop_mass*bJet1Fromtop_mass)")
    	       .Define("Lepton1FromHiggs","(Min_deltaR_bJetFromtop_leptonFromtop[0] == 0) ? lep4vecs[1]:lep4vecs[0]")
               .Define("Lepton1FromHiggs_pt","Lepton1FromHiggs.Pt()")
               .Define("Lepton1FromHiggs_eta","Lepton1FromHiggs.Eta()")
               .Define("Lepton1FromHiggs_phi","Lepton1FromHiggs.Phi()")
    	       .Define("Lepton1FromHiggs_mass","Lepton1FromHiggs.M()")
               .Define("Lepton1FromHiggs_transverse_energy","sqrt(Lepton1FromHiggs_pt*Lepton1FromHiggs_pt + Lepton1FromHiggs_mass*Lepton1FromHiggs_mass)");

    _rlm = _rlm.Define("cleanjet4Rvecs_withoutbjet",::removebjet,{"cleanjet4Rvecs","Selected_clean_jet_pt","bJet1Fromtop_pt"})
               .Define("cleanjet4Rvecs_withoutbjet_pt","ROOT::VecOps::Map(cleanjet4Rvecs_withoutbjet, FourVecPt)")
               .Define("cleanjet4Rvecs_withoutbjet_eta","ROOT::VecOps::Map(cleanjet4Rvecs_withoutbjet, FourVecEta)")
               .Define("cleanjet4Rvecs_withoutbjet_phi","ROOT::VecOps::Map(cleanjet4Rvecs_withoutbjet, FourVecPhi)")
               .Define("cleanjet4Rvecs_withoutbjet_mass","ROOT::VecOps::Map(cleanjet4Rvecs_withoutbjet, FourVecMass)")
               .Define("cleanjet4Rvecs_withoutbjet_transverse_energy","sqrt(cleanjet4Rvecs_withoutbjet_pt*cleanjet4Rvecs_withoutbjet_pt + cleanjet4Rvecs_withoutbjet_mass*cleanjet4Rvecs_withoutbjet_mass)")
               .Define("cleanjet4Rvecs_withoutbjet_number","int(cleanjet4Rvecs_withoutbjet_pt.size())")
               .Define("cleanjet4vecs_withoutbjet",::generate_4vec,{"cleanjet4Rvecs_withoutbjet_pt","cleanjet4Rvecs_withoutbjet_eta","cleanjet4Rvecs_withoutbjet_phi","cleanjet4Rvecs_withoutbjet_mass"});

    _rlm = _rlm.Define("dis_comb_jets_withoutbjet",::distinct_comb_2_bis,{"cleanjet4Rvecs_withoutbjet_number"})
               .Define("Selected_jet1","ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet, dis_comb_jets_withoutbjet[0])")
               .Define("Selected_jet1_pt","ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_pt, dis_comb_jets_withoutbjet[0])")
               .Define("Selected_jet2","ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet, dis_comb_jets_withoutbjet[1])")
               .Define("Selected_jet2_pt","ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_pt, dis_comb_jets_withoutbjet[1])")
               .Define("JetsFromHiggs4Rvecs",::sum_two_4Rvecs,{"Selected_jet1","Selected_jet2"})
               .Define("JetsFromHiggs_pt","ROOT::VecOps::Map(JetsFromHiggs4Rvecs, FourVecPt)")
               .Define("JetsFromHiggs_eta","ROOT::VecOps::Map(JetsFromHiggs4Rvecs, FourVecEta)")
               .Define("JetsFromHiggs_phi","ROOT::VecOps::Map(JetsFromHiggs4Rvecs, FourVecPhi)")
               .Define("JetsFromHiggs_mass","ROOT::VecOps::Map(JetsFromHiggs4Rvecs, FourVecMass)")
               .Define("JetsFromHiggs_energy","ROOT::VecOps::Map(JetsFromHiggs4Rvecs, FourVecEnergy)")
               .Define("JetsFromHiggs_transverse_energy","sqrt(JetsFromHiggs_pt*JetsFromHiggs_pt + JetsFromHiggs_mass*JetsFromHiggs_mass)")
               .Define("JetsFromHiggs_deltaeta","abs(ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_eta, dis_comb_jets_withoutbjet[0]) - ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_eta, dis_comb_jets_withoutbjet[1]))")
               .Define("JetsFromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_phi, dis_comb_jets_withoutbjet[0]),ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_phi, dis_comb_jets_withoutbjet[1])))")
               .Define("JetsFromHiggs_deltaR","ROOT::VecOps::DeltaR(ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_eta, dis_comb_jets_withoutbjet[0]),ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_eta, dis_comb_jets_withoutbjet[1]),ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_phi, dis_comb_jets_withoutbjet[0]),ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_phi, dis_comb_jets_withoutbjet[1]))")
               .Define("JetsFromHiggs_transverse_mass","sqrt(ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_mass, dis_comb_jets_withoutbjet[0])*ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_mass, dis_comb_jets_withoutbjet[0]) + ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_mass, dis_comb_jets_withoutbjet[1])*ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_mass, dis_comb_jets_withoutbjet[1]) + 2*ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_transverse_energy, dis_comb_jets_withoutbjet[0])*ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_transverse_energy, dis_comb_jets_withoutbjet[1]) - 2*ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_pt, dis_comb_jets_withoutbjet[0])*ROOT::VecOps::Take(cleanjet4Rvecs_withoutbjet_pt, dis_comb_jets_withoutbjet[1])*cos(JetsFromHiggs_deltaphi))");

    _rlm = _rlm.Define("Min_deltaR_JetsFromHiggs",::minDR_2,{"cleanjet4Rvecs_withoutbjet_eta","cleanjet4Rvecs_withoutbjet_eta","cleanjet4Rvecs_withoutbjet_phi","cleanjet4Rvecs_withoutbjet_phi"})
    	       .Define("Jet1FromHiggs","cleanjet4vecs_withoutbjet[Min_deltaR_JetsFromHiggs[0]]")
               .Define("Jet1FromHiggs_pt","Jet1FromHiggs.Pt()")
               .Define("Jet1FromHiggs_eta","Jet1FromHiggs.Eta()")
               .Define("Jet1FromHiggs_phi","Jet1FromHiggs.Phi()")
               .Define("Jet1FromHiggs_mass","Jet1FromHiggs.M()")
               .Define("Jet1FromHiggs_transverse_energy","sqrt(Jet1FromHiggs_pt*Jet1FromHiggs_pt + Jet1FromHiggs_mass*Jet1FromHiggs_mass)")
               .Define("Jet2FromHiggs","cleanjet4vecs_withoutbjet[Min_deltaR_JetsFromHiggs[1]]")
               .Define("Jet2FromHiggs_pt","Jet2FromHiggs.Pt()")
               .Define("Jet2FromHiggs_eta","Jet2FromHiggs.Eta()")
               .Define("Jet2FromHiggs_phi","Jet2FromHiggs.Phi()")
               .Define("Jet2FromHiggs_mass","Jet2FromHiggs.M()")
               .Define("Jet2FromHiggs_transverse_energy","sqrt(Jet2FromHiggs_pt*Jet2FromHiggs_pt + Jet2FromHiggs_mass*Jet2FromHiggs_mass)")
               .Define("Jet1FromHiggs_Jet2FromHiggs","Jet1FromHiggs + Jet2FromHiggs")
               .Define("Jet1FromHiggs_Jet2FromHiggs_pt","Jet1FromHiggs_Jet2FromHiggs.Pt()")
               .Define("Jet1FromHiggs_Jet2FromHiggs_eta","Jet1FromHiggs_Jet2FromHiggs.Eta()")
               .Define("Jet1FromHiggs_Jet2FromHiggs_phi","Jet1FromHiggs_Jet2FromHiggs.Phi()")
               .Define("Jet1FromHiggs_Jet2FromHiggs_mass","Jet1FromHiggs_Jet2FromHiggs.M()")
               .Define("Jet1FromHiggs_Jet2FromHiggs_transverse_energy","sqrt(Jet1FromHiggs_Jet2FromHiggs_pt*Jet1FromHiggs_Jet2FromHiggs_pt + Jet1FromHiggs_Jet2FromHiggs_mass*Jet1FromHiggs_Jet2FromHiggs_mass)")
               .Define("Jet1FromHiggs_Jet2FromHiggs_deltaeta","abs(Jet2FromHiggs_eta - Jet1FromHiggs_eta)")
               .Define("Jet1FromHiggs_Jet2FromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(Jet1FromHiggs_phi,Jet2FromHiggs_phi))")
               .Define("Jet1FromHiggs_Jet2FromHiggs_deltaR","ROOT::VecOps::DeltaR(Jet1FromHiggs_eta,Jet2FromHiggs_eta,Jet1FromHiggs_phi,Jet2FromHiggs_phi)")
               .Define("Jet1FromHiggs_Jet2FromHiggs_transverse_mass","sqrt(Jet1FromHiggs_mass*Jet1FromHiggs_mass + Jet2FromHiggs_mass*Jet2FromHiggs_mass + 2*Jet1FromHiggs_transverse_energy*Jet2FromHiggs_transverse_energy - 2*Jet1FromHiggs_pt*Jet2FromHiggs_pt*cos(Jet1FromHiggs_Jet2FromHiggs_deltaphi))");

    // _rlm = _rlm.Define("Min_deltaR_JetsFromHiggs",::minDR_4,{"Selected_clean_jet_eta","Selected_clean_jet_eta","Selected_clean_jet_phi","Selected_clean_jet_phi","Selected_clean_jet_pt","bJet1Fromtop_pt"})
    //     	      .Define("Jet1FromHiggs","jet4vecs[Min_deltaR_JetsFromHiggs[0]]")
    //            .Define("Jet1FromHiggs_pt","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet1FromHiggs.Pt():-100")
    //            .Define("Jet1FromHiggs_eta","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet1FromHiggs.Eta():-10")
    //            .Define("Jet1FromHiggs_phi","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet1FromHiggs.Phi():-10")
    //            .Define("Jet1FromHiggs_mass","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet1FromHiggs.M():-10")
    //            .Define("Jet1FromHiggs_transverse_energy","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? sqrt(Jet1FromHiggs_pt*Jet1FromHiggs_pt + Jet1FromHiggs_mass*Jet1FromHiggs_mass):-100")
    //            .Define("Jet2FromHiggs","jet4vecs[Min_deltaR_JetsFromHiggs[1]]")
    //            .Define("Jet2FromHiggs_pt","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet2FromHiggs.Pt():-100")
    //            .Define("Jet2FromHiggs_eta","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet2FromHiggs.Eta():-10")
    //            .Define("Jet2FromHiggs_phi","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet2FromHiggs.Phi():-10")
    //            .Define("Jet2FromHiggs_mass","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet2FromHiggs.M():-100")
    //            .Define("Jet2FromHiggs_transverse_energy","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? sqrt(Jet2FromHiggs_pt*Jet2FromHiggs_pt + Jet2FromHiggs_mass*Jet2FromHiggs_mass):-100")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs","Jet1FromHiggs + Jet2FromHiggs")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_pt","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet1FromHiggs_Jet2FromHiggs.Pt():-100")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_eta","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet1FromHiggs_Jet2FromHiggs.Eta():-10")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_phi","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet1FromHiggs_Jet2FromHiggs.Phi():-10")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_mass","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Jet1FromHiggs_Jet2FromHiggs.M():-10")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_transverse_energy","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? sqrt(Jet1FromHiggs_Jet2FromHiggs_pt*Jet1FromHiggs_Jet2FromHiggs_pt + Jet1FromHiggs_Jet2FromHiggs_mass*Jet1FromHiggs_Jet2FromHiggs_mass):-100")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_deltaeta","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? abs(Jet2FromHiggs_eta - Jet1FromHiggs_eta):-10")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_deltaphi","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? abs(ROOT::VecOps::DeltaPhi(Jet1FromHiggs_phi,Jet2FromHiggs_phi)):-10")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_deltaR","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? ROOT::VecOps::DeltaR(Jet1FromHiggs_eta,Jet2FromHiggs_eta,Jet1FromHiggs_phi,Jet2FromHiggs_phi):-10")
    //            .Define("Jet1FromHiggs_Jet2FromHiggs_transverse_mass","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? sqrt(Jet1FromHiggs_mass*Jet1FromHiggs_mass + Jet2FromHiggs_mass*Jet2FromHiggs_mass + 2*Jet1FromHiggs_transverse_energy*Jet2FromHiggs_transverse_energy - 2*Jet1FromHiggs_pt*Jet2FromHiggs_pt*cos(Jet1FromHiggs_Jet2FromHiggs_deltaphi)):-100");

    // _rlm = _rlm.Define("Min_deltaR_LeptonFromHiggs",::minDR_3,{"Jet1FromHiggs_Jet2FromHiggs_eta","Selected_lepton_eta","Jet1FromHiggs_Jet2FromHiggs_phi","Selected_lepton_phi"})
    // 	          .Define("Lepton1FromHiggs","lep4vecs[Min_deltaR_LeptonFromHiggs]")
    //            .Define("Lepton1FromHiggs_pt","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton1FromHiggs.Pt():-100")
    //            .Define("Lepton1FromHiggs_eta","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton1FromHiggs.Eta():-10")
    //            .Define("Lepton1FromHiggs_phi","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton1FromHiggs.Phi():-10")
    //            .Define("Lepton1FromHiggs_mass","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton1FromHiggs.M():-100")
    // 	          .Define("Lepton2Fromtop","(Min_deltaR_LeptonFromHiggs == 0) ? lep4vecs[1]:lep4vecs[0]")
    //            .Define("Lepton2Fromtop_pt","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton2Fromtop.Pt():-100")
    //            .Define("Lepton2Fromtop_eta","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton2Fromtop.Eta():-10")
    //            .Define("Lepton2Fromtop_phi","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton2Fromtop.Phi():-10")
    //            .Define("Lepton2Fromtop_mass","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton2Fromtop.M():-100");

    // _rlm = _rlm.Define("Min_deltaR_bJetFromtop",::minDR_3,{"Lepton2Fromtop_eta","Selected_clean_bjet_eta","Lepton2Fromtop_phi","Selected_clean_bjet_phi"})
    //            .Define("bJet1Fromtop","bjet4vecs[Min_deltaR_bJetFromtop]")
    //            .Define("bJet1Fromtop_pt","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? bJet1Fromtop.Pt():-100")
    //            .Define("bJet1Fromtop_eta","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? bJet1Fromtop.Eta():-10")
    //            .Define("bJet1Fromtop_phi","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? bJet1Fromtop.Phi():-10")
    //            .Define("bJet1Fromtop_mass","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? bJet1Fromtop.M():-100");

    // _rlm = _rlm.Define("bJet1Fromtop_genJetIdx","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Selected_clean_bjet_genJetIdx[Min_deltaR_bJetFromtop[1]]:-1")
    //            .Define("sel_bjetgentop_gencheck",::isJetFromgJet_2,{"bJet1Fromtop_pt","GenJet_pt","JFromtop_pt","bJet1Fromtop_genJetIdx"})
    //            .Define("Jet1FromHiggs_genJetIdx","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Selected_clean_jet_genJetIdx[Min_deltaR_JetsFromHiggs[0]]:-1")
    //            .Define("sel_jet1genHiggs_gencheck",::isJetFromgJet_2,{"Jet1FromHiggs_pt","GenJet_pt","JFromHiggs_pt","Jet1FromHiggs_genJetIdx"})
    //            .Define("Jet2FromHiggs_genJetIdx","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Selected_clean_jet_genJetIdx[Min_deltaR_JetsFromHiggs[1]]:-1")
    // 	          .Define("sel_jet2genHiggs_gencheck",::isJetFromgJet_2,{"Jet2FromHiggs_pt","GenJet_pt","JFromHiggs_pt","Jet2FromHiggs_genJetIdx"})
    //            .Define("sel_lepton1genHiggs_gencheck","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton1FromHiggs_pt == LeptonFromHiggs_pt:0")
    //            .Define("sel_lepton2gentop_gencheck","(Selected_clean_jet_number >=3 && Selected_clean_bjet_number >= 1) ? Lepton2Fromtop_pt == LeptonFromtop_pt:0");

    _rlm = _rlm.Define("Lepton2Fromtop_bJet1Fromtop","Lepton2Fromtop + bJet1Fromtop")
               .Define("Lepton2Fromtop_bJet1Fromtop_pt","Lepton2Fromtop_bJet1Fromtop.Pt()")      
               .Define("Lepton2Fromtop_bJet1Fromtop_eta","Lepton2Fromtop_bJet1Fromtop.Eta()") 
               .Define("Lepton2Fromtop_bJet1Fromtop_phi","Lepton2Fromtop_bJet1Fromtop.Phi()") 
               .Define("Lepton2Fromtop_bJet1Fromtop_mass","Lepton2Fromtop_bJet1Fromtop.M()") 
               .Define("Lepton2Fromtop_bJet1Fromtop_transverse_energy","sqrt(Lepton2Fromtop_bJet1Fromtop_pt*Lepton2Fromtop_bJet1Fromtop_pt + Lepton2Fromtop_bJet1Fromtop_mass*Lepton2Fromtop_bJet1Fromtop_mass)")
               .Define("Lepton2Fromtop_bJet1Fromtop_deltaeta","abs(Lepton2Fromtop_eta - bJet1Fromtop_eta)")
               .Define("Lepton2Fromtop_bJet1Fromtop_deltaphi","abs(ROOT::VecOps::DeltaPhi(Lepton2Fromtop_phi,bJet1Fromtop_phi))")
               .Define("Lepton2Fromtop_bJet1Fromtop_deltaR","ROOT::VecOps::DeltaR(Lepton2Fromtop_eta,bJet1Fromtop_eta,Lepton2Fromtop_phi,bJet1Fromtop_phi)")
               .Define("Lepton2Fromtop_bJet1Fromtop_transverse_mass","sqrt(Lepton2Fromtop_mass*Lepton2Fromtop_mass + bJet1Fromtop_mass*bJet1Fromtop_mass + 2*Lepton2Fromtop_transverse_energy*bJet1Fromtop_transverse_energy - 2*Lepton2Fromtop_pt*bJet1Fromtop_pt*cos(Lepton2Fromtop_bJet1Fromtop_deltaphi))")
               .Define("Lepton2Fromtop_Jet1FromHiggs","Lepton2Fromtop + Jet1FromHiggs")
               .Define("Lepton2Fromtop_Jet1FromHiggs_pt","Lepton2Fromtop_Jet1FromHiggs.Pt()")      
               .Define("Lepton2Fromtop_Jet1FromHiggs_eta","Lepton2Fromtop_Jet1FromHiggs.Eta()") 
               .Define("Lepton2Fromtop_Jet1FromHiggs_phi","Lepton2Fromtop_Jet1FromHiggs.Phi()") 
               .Define("Lepton2Fromtop_Jet1FromHiggs_mass","Lepton2Fromtop_Jet1FromHiggs.M()") 
               .Define("Lepton2Fromtop_Jet1FromHiggs_transverse_energy","sqrt(Lepton2Fromtop_Jet1FromHiggs_pt*Lepton2Fromtop_Jet1FromHiggs_pt + Lepton2Fromtop_Jet1FromHiggs_mass*Lepton2Fromtop_Jet1FromHiggs_mass)")
               .Define("Lepton2Fromtop_Jet1FromHiggs_deltaeta","abs(Lepton2Fromtop_eta - Jet1FromHiggs_eta)")
               .Define("Lepton2Fromtop_Jet1FromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(Lepton2Fromtop_phi,Jet1FromHiggs_phi))")
               .Define("Lepton2Fromtop_Jet1FromHiggs_deltaR","ROOT::VecOps::DeltaR(Lepton2Fromtop_eta,Jet1FromHiggs_eta,Lepton2Fromtop_phi,Jet1FromHiggs_phi)")
               .Define("Lepton2Fromtop_Jet1FromHiggs_transverse_mass","sqrt(Lepton2Fromtop_mass*Lepton2Fromtop_mass + Jet1FromHiggs_mass*Jet1FromHiggs_mass + 2*Lepton2Fromtop_transverse_energy*Jet1FromHiggs_transverse_energy - 2*Lepton2Fromtop_pt*Jet1FromHiggs_pt*cos(Lepton2Fromtop_Jet1FromHiggs_deltaphi))")
               .Define("Lepton2Fromtop_Jet2FromHiggs","Lepton2Fromtop + Jet2FromHiggs")
               .Define("Lepton2Fromtop_Jet2FromHiggs_pt","Lepton2Fromtop_Jet2FromHiggs.Pt()")      
               .Define("Lepton2Fromtop_Jet2FromHiggs_eta","Lepton2Fromtop_Jet2FromHiggs.Eta()") 
               .Define("Lepton2Fromtop_Jet2FromHiggs_phi","Lepton2Fromtop_Jet2FromHiggs.Phi()") 
               .Define("Lepton2Fromtop_Jet2FromHiggs_mass","Lepton2Fromtop_Jet2FromHiggs.M()") 
               .Define("Lepton2Fromtop_Jet2FromHiggs_transverse_energy","sqrt(Lepton2Fromtop_Jet2FromHiggs_pt*Lepton2Fromtop_Jet2FromHiggs_pt + Lepton2Fromtop_Jet2FromHiggs_mass*Lepton2Fromtop_Jet2FromHiggs_mass)")
               .Define("Lepton2Fromtop_Jet2FromHiggs_deltaeta","abs(Lepton2Fromtop_eta - Jet2FromHiggs_eta)")
               .Define("Lepton2Fromtop_Jet2FromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(Lepton2Fromtop_phi,Jet2FromHiggs_phi))")
               .Define("Lepton2Fromtop_Jet2FromHiggs_deltaR","ROOT::VecOps::DeltaR(Lepton2Fromtop_eta,Jet2FromHiggs_eta,Lepton2Fromtop_phi,Jet2FromHiggs_phi)")
               .Define("Lepton2Fromtop_Jet2FromHiggs_transverse_mass","sqrt(Lepton2Fromtop_mass*Lepton2Fromtop_mass + Jet2FromHiggs_mass*Jet2FromHiggs_mass + 2*Lepton2Fromtop_transverse_energy*Jet2FromHiggs_transverse_energy - 2*Lepton2Fromtop_pt*Jet2FromHiggs_pt*cos(Lepton2Fromtop_Jet2FromHiggs_deltaphi))")
               .Define("Lepton1FromHiggs_bJet1Fromtop","Lepton1FromHiggs + bJet1Fromtop")
               .Define("Lepton1FromHiggs_bJet1Fromtop_pt","Lepton1FromHiggs_bJet1Fromtop.Pt()")      
               .Define("Lepton1FromHiggs_bJet1Fromtop_eta","Lepton1FromHiggs_bJet1Fromtop.Eta()") 
               .Define("Lepton1FromHiggs_bJet1Fromtop_phi","Lepton1FromHiggs_bJet1Fromtop.Phi()") 
               .Define("Lepton1FromHiggs_bJet1Fromtop_mass","Lepton1FromHiggs_bJet1Fromtop.M()")
               .Define("Lepton1FromHiggs_bJet1Fromtop_transverse_energy","sqrt(Lepton1FromHiggs_bJet1Fromtop_pt*Lepton1FromHiggs_bJet1Fromtop_pt + Lepton1FromHiggs_bJet1Fromtop_mass*Lepton1FromHiggs_bJet1Fromtop_mass)")
               .Define("Lepton1FromHiggs_bJet1Fromtop_deltaeta","abs(Lepton1FromHiggs_eta - bJet1Fromtop_eta)")
               .Define("Lepton1FromHiggs_bJet1Fromtop_deltaphi","abs(ROOT::VecOps::DeltaPhi(Lepton1FromHiggs_phi,bJet1Fromtop_phi))")
               .Define("Lepton1FromHiggs_bJet1Fromtop_deltaR","ROOT::VecOps::DeltaR(Lepton1FromHiggs_eta,bJet1Fromtop_eta,Lepton1FromHiggs_phi,bJet1Fromtop_phi)")
               .Define("Lepton1FromHiggs_bJet1Fromtop_transverse_mass","sqrt(Lepton1FromHiggs_mass*Lepton1FromHiggs_mass + bJet1Fromtop_mass*bJet1Fromtop_mass + 2*Lepton1FromHiggs_transverse_energy*bJet1Fromtop_transverse_energy - 2*Lepton1FromHiggs_pt*bJet1Fromtop_pt*cos(Lepton1FromHiggs_bJet1Fromtop_deltaphi))")
               .Define("Lepton1FromHiggs_Jet1FromHiggs","Lepton1FromHiggs + Jet1FromHiggs")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_pt","Lepton1FromHiggs_Jet1FromHiggs.Pt()")      
               .Define("Lepton1FromHiggs_Jet1FromHiggs_eta","Lepton1FromHiggs_Jet1FromHiggs.Eta()") 
               .Define("Lepton1FromHiggs_Jet1FromHiggs_phi","Lepton1FromHiggs_Jet1FromHiggs.Phi()") 
               .Define("Lepton1FromHiggs_Jet1FromHiggs_mass","Lepton1FromHiggs_Jet1FromHiggs.M()") 
               .Define("Lepton1FromHiggs_Jet1FromHiggs_transverse_energy","sqrt(Lepton1FromHiggs_Jet1FromHiggs_pt*Lepton1FromHiggs_Jet1FromHiggs_pt + Lepton1FromHiggs_Jet1FromHiggs_mass*Lepton1FromHiggs_Jet1FromHiggs_mass)")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_deltaeta","abs(Lepton1FromHiggs_eta - Jet1FromHiggs_eta)")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(Lepton1FromHiggs_phi,Jet1FromHiggs_phi))")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_deltaR","ROOT::VecOps::DeltaR(Lepton1FromHiggs_eta,Jet1FromHiggs_eta,Lepton1FromHiggs_phi,Jet1FromHiggs_phi)")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_transverse_mass","sqrt(Lepton1FromHiggs_mass*Lepton1FromHiggs_mass + Jet1FromHiggs_mass*Jet1FromHiggs_mass + 2*Lepton1FromHiggs_transverse_energy*Jet1FromHiggs_transverse_energy - 2*Lepton1FromHiggs_pt*Jet1FromHiggs_pt*cos(Lepton1FromHiggs_Jet1FromHiggs_deltaphi))")
               .Define("Lepton1FromHiggs_Jet2FromHiggs","Lepton1FromHiggs + Jet2FromHiggs")
               .Define("Lepton1FromHiggs_Jet2FromHiggs_pt","Lepton1FromHiggs_Jet2FromHiggs.Pt()")      
               .Define("Lepton1FromHiggs_Jet2FromHiggs_eta","Lepton1FromHiggs_Jet2FromHiggs.Eta()") 
               .Define("Lepton1FromHiggs_Jet2FromHiggs_phi","Lepton1FromHiggs_Jet2FromHiggs.Phi()") 
               .Define("Lepton1FromHiggs_Jet2FromHiggs_mass","Lepton1FromHiggs_Jet2FromHiggs.M()") 
               .Define("Lepton1FromHiggs_Jet2FromHiggs_transverse_energy","sqrt(Lepton1FromHiggs_Jet2FromHiggs_pt*Lepton1FromHiggs_Jet2FromHiggs_pt + Lepton1FromHiggs_Jet2FromHiggs_mass*Lepton1FromHiggs_Jet2FromHiggs_mass)")
               .Define("Lepton1FromHiggs_Jet2FromHiggs_deltaeta","abs(Lepton1FromHiggs_eta - Jet2FromHiggs_eta)")
               .Define("Lepton1FromHiggs_Jet2FromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(Lepton1FromHiggs_phi,Jet2FromHiggs_phi))")
               .Define("Lepton1FromHiggs_Jet2FromHiggs_deltaR","ROOT::VecOps::DeltaR(Lepton1FromHiggs_eta,Jet2FromHiggs_eta,Lepton1FromHiggs_phi,Jet2FromHiggs_phi)")
               .Define("Lepton1FromHiggs_Jet2FromHiggs_transverse_mass","sqrt(Lepton1FromHiggs_mass*Lepton1FromHiggs_mass + Jet2FromHiggs_mass*Jet2FromHiggs_mass + 2*Lepton1FromHiggs_transverse_energy*Jet2FromHiggs_transverse_energy - 2*Lepton1FromHiggs_pt*Jet2FromHiggs_pt*cos(Lepton1FromHiggs_Jet2FromHiggs_deltaphi))")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs","Lepton1FromHiggs + Jet1FromHiggs_Jet2FromHiggs")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_pt","Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs.Pt()")      
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_eta","Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs.Eta()") 
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_phi","Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs.Phi()") 
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_mass","Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs.M()") 
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_transverse_energy","sqrt(Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_pt*Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_pt + Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_mass*Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_mass)")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_deltaeta","abs(Lepton1FromHiggs_eta - Jet1FromHiggs_Jet2FromHiggs_eta)")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(Lepton1FromHiggs_phi,Jet1FromHiggs_Jet2FromHiggs_phi))")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_deltaR","ROOT::VecOps::DeltaR(Lepton1FromHiggs_eta,Jet1FromHiggs_Jet2FromHiggs_eta,Lepton1FromHiggs_phi,Jet1FromHiggs_Jet2FromHiggs_phi)")
               .Define("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_transverse_mass","sqrt(Lepton1FromHiggs_mass*Lepton1FromHiggs_mass + Jet1FromHiggs_Jet2FromHiggs_mass*Jet1FromHiggs_Jet2FromHiggs_mass + 2*Lepton1FromHiggs_transverse_energy*Jet1FromHiggs_Jet2FromHiggs_transverse_energy - 2*Lepton1FromHiggs_pt*Jet1FromHiggs_Jet2FromHiggs_pt*cos(Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_deltaphi))")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs","Lepton2Fromtop_bJet1Fromtop + Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_pt","L2bJ1Fromtop_L1J1J2FromHiggs.Pt()")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_eta","L2bJ1Fromtop_L1J1J2FromHiggs.Eta()")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_phi","L2bJ1Fromtop_L1J1J2FromHiggs.Phi()")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_mass","L2bJ1Fromtop_L1J1J2FromHiggs.M()")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_px","L2bJ1Fromtop_L1J1J2FromHiggs_pt*cos(L2bJ1Fromtop_L1J1J2FromHiggs_phi)")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_py","L2bJ1Fromtop_L1J1J2FromHiggs_pt*sin(L2bJ1Fromtop_L1J1J2FromHiggs_phi)")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_transverse_energy","sqrt(L2bJ1Fromtop_L1J1J2FromHiggs_pt*L2bJ1Fromtop_L1J1J2FromHiggs_pt + L2bJ1Fromtop_L1J1J2FromHiggs_mass*L2bJ1Fromtop_L1J1J2FromHiggs_mass)")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_deltaeta","abs(Lepton2Fromtop_bJet1Fromtop_eta - Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_eta)")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(Lepton2Fromtop_bJet1Fromtop_phi,Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_phi))")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_deltaR","ROOT::VecOps::DeltaR(Lepton2Fromtop_bJet1Fromtop_eta,Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_eta,Lepton2Fromtop_bJet1Fromtop_phi,Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_phi)")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_transverse_mass","sqrt(Lepton2Fromtop_bJet1Fromtop_mass*Lepton2Fromtop_bJet1Fromtop_mass + Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_mass*Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_mass + 2*Lepton2Fromtop_bJet1Fromtop_transverse_energy*Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_transverse_energy - 2*Lepton2Fromtop_bJet1Fromtop_pt*Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_pt*cos(L2bJ1Fromtop_L1J1J2FromHiggs_deltaphi))")
               .Define("L2bJ1Fromtop_L1J1J2FromHiggs_Relative_St","(Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_pt + Lepton2Fromtop_bJet1Fromtop_pt) / St")
               .Define("Tprime_transverse_mass","sqrt((L2bJ1Fromtop_L1J1J2FromHiggs_transverse_energy)*(L2bJ1Fromtop_L1J1J2FromHiggs_transverse_energy) - ((L2bJ1Fromtop_L1J1J2FromHiggs_px )*(L2bJ1Fromtop_L1J1J2FromHiggs_px) + (L2bJ1Fromtop_L1J1J2FromHiggs_py)*(L2bJ1Fromtop_L1J1J2FromHiggs_py)))");
            //    .Define("Tprime_transverse_mass","sqrt(L2bJ1Fromtop_L1J1J2FromHiggs_mass*L2bJ1Fromtop_L1J1J2FromHiggs_mass + 2*L2bJ1Fromtop_L1J1J2FromHiggs_transverse_energy*MET_pt - 2*Lepton2Fromtop_bJet1Fromtop_pt*MET_pt*cos(abs(ROOT::VecOps::DeltaPhi(L2bJ1Fromtop_L1J1J2FromHiggs_phi,MET_phi))))");

    //================================Store variables in tree=======================================//
    // define variables that you want to store
    //==============================================================================================//
    addVartoStore("run");
    addVartoStore("luminosityBlock");
    addVartoStore("event");
    // addVartoStore("evWeight.*");
    addVartoStore("evWeight");
    if(!_isData)
    {
        addVartoStore("genWeight"); // doesn't work for data
	    addVartoStore("pugenWeight");// doesn't work for data
	    // addVartoStore("puWeight");// doesn't work for data
        addVartoStore("btag_SF_case1");
        addVartoStore("btagWeight_case1_central");
        if(year == "2018")
        {
            addVartoStore("muonID_SF");
            addVartoStore("muon_iso_weight");
            addVartoStore("evWeight_MuonIDISO");
        }
    }

    //electron
    addVartoStore("nElectron");
    addVartoStore("Electron_cutBased");
    addVartoStore("Electron_pt");
    addVartoStore("Electron_eta");    
    addVartoStore("Electron_sip3d");
    addVartoStore("Electron_miniPFRelIso_all");
    addVartoStore("goodElectronslooseTTSL");
    addVartoStore("Selected_electron_looseTTSL_pt");
    addVartoStore("Selected_electron_looseTTSL_leading_pt");
    addVartoStore("Selected_electron_looseTTSL_subleading_pt");
    addVartoStore("Selected_electron_looseTTSL_sum_two_electrons_pt");
    addVartoStore("Selected_electron_looseTTSL_sum_all_electrons_pt");
    addVartoStore("Selected_electron_looseTTSL_eta");
    addVartoStore("Selected_electron_looseTTSL_leading_eta");
    addVartoStore("Selected_electron_looseTTSL_subleading_eta");
    addVartoStore("Selected_electron_looseTTSL_phi");
    addVartoStore("Selected_electron_looseTTSL_mass");
    addVartoStore("Selected_electron_looseTTSL_transverse_energy");
    addVartoStore("Selected_electron_looseTTSL_charge");
    addVartoStore("Selected_electron_looseTTSL_charge_sum");
    addVartoStore("Selected_electron_looseTTSL_number");
    addVartoStore("Selected_electron_looseTTSL_deltaeta");
    addVartoStore("Selected_electron_looseTTSL_deltaphi");
    addVartoStore("Selected_electron_looseTTSL_deltaR");
    addVartoStore("Selected_electron_looseTTSL_miniPFRelIso_all");
    addVartoStore("Selected_electron_looseTTSL_miniPFRelIso_chg");
    addVartoStore("Selected_electron_looseTTSL_pfRelIso03_all");
    addVartoStore("Selected_electron_looseTTSL_pfRelIso03_chg");
    if(!_isData)
    {
        addVartoStore("Selected_electron_looseTTSL_genPartIdx");
        addVartoStore("Selected_electron_looseTTSL_genPartpdgId");
        addVartoStore("Selected_electron_looseTTSL_genPartIdxMother");
        addVartoStore("Selected_electron_looseTTSL_genPartpdgIdMother");
        addVartoStore("Selected_electron_looseTTSL_genPartFlav");
    }
    addVartoStore("Selected_electron_looseTTSL_pdgId");
    addVartoStore("Selected_electron_looseTTSL_ip3d");
    addVartoStore("Selected_electron_looseTTSL_sip3d");
    addVartoStore("Selected_electron_looseTTSL_idx");
    addVartoStore("goodElectrons");
    addVartoStore("Selected_electron_pt");
    addVartoStore("Selected_electron_leading_pt");
    addVartoStore("Selected_electron_subleading_pt");
    addVartoStore("Selected_electron_sum_two_electrons_pt");
    addVartoStore("Selected_electron_sum_all_electrons_pt");
    addVartoStore("Selected_electron_eta");
    addVartoStore("Selected_electron_leading_eta");
    addVartoStore("Selected_electron_subleading_eta");
    addVartoStore("Selected_electron_phi");
    addVartoStore("Selected_electron_mass");
    addVartoStore("Selected_electron_transverse_energy");
    addVartoStore("Selected_electron_charge");
    addVartoStore("Selected_electron_charge_sum");
    addVartoStore("Selected_electron_number");
    addVartoStore("Selected_electron_deltaeta");
    addVartoStore("Selected_electron_deltaphi");
    addVartoStore("Selected_electron_deltaR");
    addVartoStore("Selected_electron_miniPFRelIso_all");
    addVartoStore("Selected_electron_miniPFRelIso_chg");
    addVartoStore("Selected_electron_pfRelIso03_all");
    addVartoStore("Selected_electron_pfRelIso03_chg");
    addVartoStore("Selected_electron_HEEPconditions");
    if(!_isData)
    {
        addVartoStore("Selected_electron_genPartIdx");
        addVartoStore("Selected_electron_genPartpdgId");
        addVartoStore("Selected_electron_genPartIdxMother");
        addVartoStore("Selected_electron_genPartpdgIdMother");
        addVartoStore("Selected_electron_genPartFlav");
    }
    addVartoStore("Selected_electron_pdgId");
    addVartoStore("Selected_electron_ip3d");
    addVartoStore("Selected_electron_sip3d");
    addVartoStore("Selected_electron_idx");
    // addVartoStore("elec4vecs");
    // addVartoStore("Vectorial_sum_two_electrons");
    // addVartoStore("Vectorial_sum_two_electrons_pt");
    // addVartoStore("Vectorial_sum_two_electrons_eta");
    // addVartoStore("Vectorial_sum_two_electrons_phi");
    // addVartoStore("Vectorial_sum_two_electrons_mass");
    
    //muon
    addVartoStore("nMuon");
    addVartoStore("Muon_pt");
    addVartoStore("Muon_eta");
    addVartoStore("Muon_looseId");
    addVartoStore("Muon_mediumId");
    addVartoStore("Muon_tightId");
    addVartoStore("Muon_miniPFRelIso_all");
    addVartoStore("Muon_sip3d");
    addVartoStore("goodMuonsloose");
    addVartoStore("Selected_muon_loose_pt");
    addVartoStore("Selected_muon_loose_eta");
    addVartoStore("Selected_muon_loose_phi");
    addVartoStore("Selected_muon_loose_mass");
    addVartoStore("Selected_muon_loose_number");
    // addVartoStore("muon4vecs_loose");
    addVartoStore("goodMuonslooseTTSL");
    addVartoStore("Selected_muon_looseTTSL_pt");
    addVartoStore("Selected_muon_looseTTSL_leading_pt");
    addVartoStore("Selected_muon_looseTTSL_subleading_pt");
    addVartoStore("Selected_muon_looseTTSL_sum_two_muons_pt");
    addVartoStore("Selected_muon_looseTTSL_sum_all_muons_pt");   
    addVartoStore("Selected_muon_looseTTSL_eta");
    addVartoStore("Selected_muon_looseTTSL_leading_eta");
    addVartoStore("Selected_muon_looseTTSL_subleading_eta");
    addVartoStore("Selected_muon_looseTTSL_phi");
    addVartoStore("Selected_muon_looseTTSL_mass");
    addVartoStore("Selected_muon_looseTTSL_transverse_energy");
    addVartoStore("Selected_muon_looseTTSL_charge");
    addVartoStore("Selected_muon_looseTTSL_charge_sum");
    addVartoStore("Selected_muon_looseTTSL_number");
    addVartoStore("Selected_muon_looseTTSL_deltaeta");
    addVartoStore("Selected_muon_looseTTSL_deltaphi");
    addVartoStore("Selected_muon_looseTTSL_deltaR");
    addVartoStore("Selected_muon_looseTTSL_miniPFRelIso_all");
    addVartoStore("Selected_muon_looseTTSL_miniPFRelIso_chg");
    addVartoStore("Selected_muon_looseTTSL_pfRelIso03_all");
    addVartoStore("Selected_muon_looseTTSL_pfRelIso03_chg");
    if(!_isData)
    {
        addVartoStore("Selected_muon_looseTTSL_genPartIdx");
        addVartoStore("Selected_muon_looseTTSL_genPartpdgId");
        addVartoStore("Selected_muon_looseTTSL_genPartIdxMother");
        addVartoStore("Selected_muon_looseTTSL_genPartpdgIdMother");
        addVartoStore("Selected_muon_looseTTSL_genPartFlav");
    }
    addVartoStore("Selected_muon_looseTTSL_pdgId");
    addVartoStore("Selected_muon_looseTTSL_ip3d");
    addVartoStore("Selected_muon_looseTTSL_sip3d");
    addVartoStore("Selected_muon_looseTTSL_idx");
    addVartoStore("goodMuons");
    addVartoStore("Selected_muon_pt");
    addVartoStore("Selected_muon_leading_pt");
    addVartoStore("Selected_muon_subleading_pt");
    addVartoStore("Selected_muon_sum_two_muons_pt");
    addVartoStore("Selected_muon_sum_all_muons_pt");   
    addVartoStore("Selected_muon_eta");
    addVartoStore("Selected_muon_leading_eta");
    addVartoStore("Selected_muon_subleading_eta");
    addVartoStore("Selected_muon_phi");
    addVartoStore("Selected_muon_mass");
    addVartoStore("Selected_muon_transverse_energy");
    addVartoStore("Selected_muon_charge");
    addVartoStore("Selected_muon_charge_sum");
    addVartoStore("Selected_muon_number");
    addVartoStore("Selected_muon_deltaeta");
    addVartoStore("Selected_muon_deltaphi");
    addVartoStore("Selected_muon_deltaR");
    addVartoStore("Selected_muon_miniPFRelIso_all");
    addVartoStore("Selected_muon_miniPFRelIso_chg");
    addVartoStore("Selected_muon_pfRelIso03_all");
    addVartoStore("Selected_muon_pfRelIso03_chg");
    addVartoStore("Selected_muon_pfRelIso04_all");
    if(!_isData)
    {
        addVartoStore("Selected_muon_genPartIdx");
        addVartoStore("Selected_muon_genPartpdgId");
        addVartoStore("Selected_muon_genPartIdxMother");
        addVartoStore("Selected_muon_genPartpdgIdMother");
        addVartoStore("Selected_muon_genPartFlav");
    }
    addVartoStore("Selected_muon_pdgId");
    addVartoStore("Selected_muon_ip3d");
    addVartoStore("Selected_muon_sip3d");
    addVartoStore("Selected_muon_idx");
    // addVartoStore("muon4vecs");
    // addVartoStore("Vectorial_sum_two_muons");
    // addVartoStore("Vectorial_sum_two_muons_pt");
    // addVartoStore("Vectorial_sum_two_muons_eta");
    // addVartoStore("Vectorial_sum_two_muons_phi");
    // addVartoStore("Vectorial_sum_two_muons_mass");

    //lepton
    addVartoStore("Selected_lepton_looseTTSL_pt");
    addVartoStore("Selected_lepton_looseTTSL_leading_pt");
    addVartoStore("Selected_lepton_looseTTSL_subleading_pt");
    addVartoStore("Selected_lepton_looseTTSL_sum_two_leptons_pt");
    addVartoStore("Selected_lepton_looseTTSL_sum_all_leptons_pt");
    addVartoStore("Selected_lepton_looseTTSL_eta");
    addVartoStore("Selected_lepton_looseTTSL_leading_eta");
    addVartoStore("Selected_lepton_looseTTSL_subleading_eta");
    addVartoStore("Selected_lepton_looseTTSL_phi");
    addVartoStore("Selected_lepton_looseTTSL_mass");
    addVartoStore("Selected_lepton_looseTTSL_transverse_energy");
    addVartoStore("Selected_lepton_looseTTSL_charge");
    addVartoStore("Selected_lepton_looseTTSL_charge_sum");
    addVartoStore("Selected_lepton_looseTTSL_number");
    addVartoStore("Selected_lepton_looseTTSL_deltaeta");
    addVartoStore("Selected_lepton_looseTTSL_deltaphi");
    addVartoStore("Selected_lepton_looseTTSL_deltaR");
    addVartoStore("Selected_lepton_looseTTSL_miniPFRelIso_all");
    addVartoStore("Selected_lepton_looseTTSL_miniPFRelIso_chg");
    addVartoStore("Selected_lepton_looseTTSL_pfRelIso03_all");
    addVartoStore("Selected_lepton_looseTTSL_pfRelIso03_chg");
    if(!_isData)
    {
        addVartoStore("Selected_lepton_looseTTSL_genPartIdx");
        addVartoStore("Selected_lepton_looseTTSL_genPartpdgId");
        addVartoStore("Selected_lepton_looseTTSL_genPartIdxMother");
        addVartoStore("Selected_lepton_looseTTSL_genPartpdgIdMother");
        addVartoStore("Selected_lepton_looseTTSL_genPartFlav");
    }
    addVartoStore("Selected_lepton_looseTTSL_pdgId");
    addVartoStore("Selected_lepton1_looseTTSL_pdgId");
    addVartoStore("Selected_lepton2_looseTTSL_pdgId");
    addVartoStore("Selected_lepton12_looseTTSL_pdgId");
    addVartoStore("Selected_lepton_looseTTSL_ip3d");
    addVartoStore("Selected_lepton_looseTTSL_sip3d");
    addVartoStore("Selected_lepton_looseTTSL_idx");
    addVartoStore("lep4vecs_looseTTSL");
    // addVartoStore("lep4Rvecs_looseTTSL");
    // addVartoStore("Vectorial_sum_two_leptons_looseTTSL");
    addVartoStore("Vectorial_sum_two_leptons_looseTTSL_pt");
    addVartoStore("Vectorial_sum_two_leptons_looseTTSL_eta");
    addVartoStore("Vectorial_sum_two_leptons_looseTTSL_phi");
    addVartoStore("Vectorial_sum_two_leptons_looseTTSL_mass");
    addVartoStore("Selected_lepton_pt");
    addVartoStore("Selected_lepton_leading_pt");
    addVartoStore("Selected_lepton_subleading_pt");
    addVartoStore("Selected_lepton_sum_two_leptons_pt");
    addVartoStore("Selected_lepton_sum_all_leptons_pt");
    addVartoStore("Selected_lepton_eta");
    addVartoStore("Selected_lepton_leading_eta");
    addVartoStore("Selected_lepton_subleading_eta");
    addVartoStore("Selected_lepton_phi");
    addVartoStore("Selected_lepton_mass");
    addVartoStore("Selected_lepton_transverse_energy");
    addVartoStore("Selected_lepton_charge");
    addVartoStore("Selected_lepton_charge_sum");
    addVartoStore("Selected_lepton_number");
    addVartoStore("Selected_lepton_deltaeta");
    addVartoStore("Selected_lepton_deltaphi");
    addVartoStore("Selected_lepton_deltaR");
    addVartoStore("Selected_lepton_miniPFRelIso_all");
    addVartoStore("Selected_lepton_miniPFRelIso_chg");
    addVartoStore("Selected_lepton_pfRelIso03_all");
    addVartoStore("Selected_lepton_pfRelIso03_chg");
    if(!_isData)
    {
        addVartoStore("Selected_lepton_genPartIdx");
        addVartoStore("Selected_lepton_genPartpdgId");
        addVartoStore("Selected_lepton_genPartIdxMother");
        addVartoStore("Selected_lepton_genPartpdgIdMother");
        addVartoStore("Selected_lepton_genPartFlav");
    }
    addVartoStore("Selected_lepton_pdgId");
    addVartoStore("Selected_lepton1_pdgId");
    addVartoStore("Selected_lepton2_pdgId");
    addVartoStore("Selected_lepton12_pdgId");
    addVartoStore("Selected_lepton_ip3d");
    addVartoStore("Selected_lepton_sip3d");
    addVartoStore("Selected_lepton_idx");
    addVartoStore("lep4vecs");
    // addVartoStore("lep4Rvecs");
    addVartoStore("Selected_lepton_energy");
    addVartoStore("Selected_lepton_px");
    addVartoStore("Selected_lepton_py");
    addVartoStore("Selected_lepton_pz");
    addVartoStore("Vectorial_sum_two_leptons");
    addVartoStore("Vectorial_sum_two_leptons_pt");
    addVartoStore("Vectorial_sum_two_leptons_eta");
    addVartoStore("Vectorial_sum_two_leptons_phi");
    addVartoStore("Vectorial_sum_two_leptons_mass");
    addVartoStore("Mass_dilepton_false");

    //MET
    // addVartoStore("goodMET");
    // addVartoStore("Selected_MET_pt");
    // addVartoStore("Selected_MET_sum_pt_unclustered");
    // addVartoStore("Selected_MET_phi");
    // addVartoStore("Selected_MET_sum_transverse_energy");
    addVartoStore("MET_pt");
    addVartoStore("MET_sumPtUnclustered");
    addVartoStore("MET_phi");
    addVartoStore("MET_sumEt");
    addVartoStore("MET_px");
    addVartoStore("MET_py");
    
    //jets
    addVartoStore("Selected_jet_pt");
    addVartoStore("Selected_jet_leading_pt");
    addVartoStore("Selected_jet_subleading_pt");
    addVartoStore("Selected_jet_subsubleading_pt");
    addVartoStore("Selected_jet_Ht");
    addVartoStore("Selected_jet_eta");
    addVartoStore("Selected_jet_phi");
    addVartoStore("Selected_jet_mass");
    addVartoStore("Selected_jet_number");
    addVartoStore("Selected_jet_btagDeepB");
    addVartoStore("Selected_jet_btagDeepFlavB");
    if(!_isData)
    {
        addVartoStore("Selected_jet_genJetIdx");
        addVartoStore("Selected_jet_partonFlavour");
        addVartoStore("Selected_jet_hadronFlavour");
    }
    addVartoStore("Selected_jet_jetId");
    addVartoStore("jet4vecs");
    addVartoStore("Selected_bjet_pt");
    addVartoStore("Selected_bjet_leading_pt");
    addVartoStore("Selected_bjet_subleading_pt");
    addVartoStore("Selected_bjet_subsubleading_pt");
    addVartoStore("Selected_bjet_sum_all_jets_pt");
    addVartoStore("Selected_bjet_eta");
    addVartoStore("Selected_bjet_phi");
    addVartoStore("Selected_bjet_mass");
    addVartoStore("Selected_bjet_number");
    addVartoStore("Selected_bjet_btagDeepB");
    addVartoStore("Selected_bjet_btagDeepFlavB");
    if(!_isData)
    {
        addVartoStore("Selected_bjet_genJetIdx");
        addVartoStore("Selected_bjet_partonFlavour");
        addVartoStore("Selected_bjet_hadronFlavour");
    }
    addVartoStore("Selected_bjet_jetId");
    addVartoStore("bjet4vecs");
    addVartoStore("Lepton_Jet_overlap");
    addVartoStore("Selected_clean_jet_pt");
    addVartoStore("Selected_clean_jet_leading_pt");
    addVartoStore("Selected_clean_jet_subleading_pt");
    addVartoStore("Selected_clean_jet_subsubleading_pt");
    addVartoStore("Selected_clean_jet_Ht");
    addVartoStore("Selected_clean_jet_eta");
    addVartoStore("Selected_clean_jet_phi");
    addVartoStore("Selected_clean_jet_mass");
    addVartoStore("Selected_clean_jet_transverse_energy");
    addVartoStore("Selected_clean_jet_number");
    addVartoStore("Selected_clean_jet_btagDeepB");
    addVartoStore("Selected_clean_jet_btagDeepFlavB");
    if(!_isData)
    {
        addVartoStore("Selected_clean_jet_genJetIdx");
        addVartoStore("Selected_clean_jet_partonFlavour");
        addVartoStore("Selected_clean_jet_hadronFlavour");
    }
    addVartoStore("Selected_clean_jet_jetId");
    addVartoStore("cleanjet4vecs");
    // addVartoStore("cleanjet4Rvecs");
    addVartoStore("Selected_clean_jet_energy");
    addVartoStore("Selected_clean_jet_px");
    addVartoStore("Selected_clean_jet_py");
    addVartoStore("Selected_clean_jet_pz");
    addVartoStore("Lepton_bJet_overlap");
    addVartoStore("Selected_clean_bjet_pt");
    addVartoStore("Selected_clean_bjet_leading_pt");
    addVartoStore("Selected_clean_bjet_subleading_pt");
    addVartoStore("Selected_clean_bjet_subsubleading_pt");
    addVartoStore("Selected_clean_bjet_sum_all_jets_pt");
    addVartoStore("Selected_clean_bjet_eta");
    addVartoStore("Selected_clean_bjet_phi");
    addVartoStore("Selected_clean_bjet_mass");
    addVartoStore("Selected_clean_bjet_transverse_energy");
    addVartoStore("Selected_clean_bjet_number");
    addVartoStore("Selected_clean_bjet_btagDeepB");
    addVartoStore("Selected_clean_bjet_btagDeepFlavB");
    if(!_isData)
    {
        addVartoStore("Selected_clean_bjet_genJetIdx");
        addVartoStore("Selected_clean_bjet_partonFlavour");
        addVartoStore("Selected_clean_bjet_hadronFlavour");
    }
    addVartoStore("Selected_clean_bjet_jetId");
    addVartoStore("cleanbjet4vecs");
    // addVartoStore("cleanbjet4Rvecs");
    addVartoStore("Selected_clean_bjet_energy");
    addVartoStore("Selected_clean_bjet_px");
    addVartoStore("Selected_clean_bjet_py");
    addVartoStore("Selected_clean_bjet_pz");
    // addVartoStore("dis_comb");
    // addVartoStore("Selected_clean_jet1");
    addVartoStore("Selected_clean_jet1_btagDeepFlavB");
    // addVartoStore("Selected_clean_jet2");
    addVartoStore("Selected_clean_jet2_btagDeepFlavB");
    // addVartoStore("Selected_clean_jet3");
    addVartoStore("Selected_clean_jet3_btagDeepFlavB");
    addVartoStore("Selected_clean_jet_condition");
    addVartoStore("Selected_clean_jet_condition_btagDeepFlavB");
    // addVartoStore("Vectorial_sum_three_clean_jets");
    addVartoStore("Vectorial_sum_three_clean_jets_pt");
    addVartoStore("Vectorial_sum_three_clean_jets_eta");
    addVartoStore("Vectorial_sum_three_clean_jets_phi");
    addVartoStore("Vectorial_sum_three_clean_jets_mass");
    addVartoStore("Vectorial_sum_three_clean_jets_mass_offset");
    addVartoStore("Vectorial_sum_three_clean_jets_mass_min");
    // // addVartoStore("Vectorial_sum_two_clean_jets");
    // addVartoStore("Vectorial_sum_two_clean_jets_mass");
    // addVartoStore("Vectorial_sum_two_clean_jets_mass_offset");
    // addVartoStore("Vectorial_sum_two_clean_jets_mass_min");

    //combination of particles
    addVartoStore("St");
    addVartoStore("St_with_MET");
    // addVartoStore("Selected_clean_jet_btagDeepFlavB_condition");
    // addVartoStore("Selected_clean_third_jet");
    // addVartoStore("sumparticles4vecs");
    // addVartoStore("Vectorial_sum_particles_pt");
    // addVartoStore("Vectorial_sum_particles_eta");
    // addVartoStore("Vectorial_sum_particles_phi");
    // addVartoStore("Vectorial_sum_particles_mass");
    // addVartoStore("deltaR_Selected_leading_lepton_Selected_clean_leading_jet");
    // addVartoStore("deltaR_Selected_leading_lepton_Selected_clean_subleading_jet");
    // addVartoStore("deltaR_Selected_leading_lepton_Selected_clean_subsubleading_jet");
    // addVartoStore("deltaR_Selected_leading_lepton_Selected_clean_leading_bjet");
    // addVartoStore("deltaR_Selected_subleading_lepton_Selected_clean_leading_jet");
    // addVartoStore("deltaR_Selected_subleading_lepton_Selected_clean_subleading_jet");
    // addVartoStore("deltaR_Selected_subleading_lepton_Selected_clean_subsubleading_jet");
    // addVartoStore("deltaR_Selected_subleading_lepton_Selected_clean_leading_bjet");
    // addVartoStore("Selected_leading_lepton_Pt_over_Selected_clean_leading_jet_Pt");
    // addVartoStore("Selected_leading_lepton_Pt_over_Selected_clean_subleading_jet_Pt");
    // addVartoStore("Selected_leading_lepton_Pt_over_Selected_clean_subsubleading_jet_Pt");
    // addVartoStore("Selected_leading_lepton_Pt_over_Selected_clean_leading_bjet_Pt");
    // addVartoStore("Selected_subleading_lepton_Pt_over_Selected_clean_leading_jet_Pt");
    // addVartoStore("Selected_subleading_lepton_Pt_over_Selected_clean_subleading_jet_Pt");
    // addVartoStore("Selected_subleading_lepton_Pt_over_Selected_clean_subsubleading_jet_Pt");
    // addVartoStore("Selected_subleading_lepton_Pt_over_Selected_clean_leading_bjet_Pt");
    // addVartoStore("deltaR_Selected_lepton_Selected_clean_jet");
    // addVartoStore("deltaR_Selected_lepton_Selected_clean_bjet");
    // addVartoStore("deltaR_Selected_muon_Selected_clean_jet");
    // addVartoStore("deltaR_Selected_muon_Selected_clean_bjet");
    // addVartoStore("Pt_lepton_over_Pt_clean_jet");
    // addVartoStore("Pt_lepton_over_Pt_clean_bjet");
    // addVartoStore("Pt_muon_over_Pt_clean_jet");
    // addVartoStore("Pt_muon_over_Pt_clean_bjet");
    // addVartoStore("Min_deltaR_Selected_lepton_Selected_clean_jet");
    // addVartoStore("Min_deltaR_Selected_lepton_Selected_clean_bjet");
    // addVartoStore("Min_deltaR_Selected_muon_Selected_clean_jet");
    // addVartoStore("Min_deltaR_Selected_muon_Selected_clean_bjet");
    // addVartoStore("Min_Pt_lepton_over_Pt_clean_jet");
    // addVartoStore("Min_Pt_lepton_over_Pt_clean_bjet");
    // addVartoStore("Min_Pt_muon_over_Pt_clean_jet");
    // addVartoStore("Min_Pt_muon_over_Pt_clean_bjet");
    // addVartoStore("dis_comb_leptonbjet");
    // addVartoStore("Selected_lepton1");
    addVartoStore("Selected_lepton1_pt");
    // addVartoStore("Selected_clean_bjet1");
    addVartoStore("Selected_clean_bjet1_pt");
    // addVartoStore("LeptonbjetFromtop4Rvecs");
    addVartoStore("LeptonbjetFromtop_pt");
    addVartoStore("LeptonbjetFromtop_eta");
    addVartoStore("LeptonbjetFromtop_phi");
    addVartoStore("LeptonbjetFromtop_mass");
    addVartoStore("LeptonbjetFromtop_energy");
    addVartoStore("LeptonbjetFromtop_transverse_energy");
    addVartoStore("LeptonbjetFromtop_deltaeta");
    addVartoStore("LeptonbjetFromtop_deltaphi");
    addVartoStore("LeptonbjetFromtop_deltaR");
    addVartoStore("LeptonbjetFromtop_transverse_mass");
    addVartoStore("Min_deltaR_bJetFromtop_leptonFromtop");
    addVartoStore("Lepton1FromHiggs");
    addVartoStore("Lepton1FromHiggs_pt");
    addVartoStore("Lepton1FromHiggs_eta");
    addVartoStore("Lepton1FromHiggs_phi");
    addVartoStore("Lepton1FromHiggs_mass");
    addVartoStore("Lepton1FromHiggs_transverse_energy");
    addVartoStore("Lepton2Fromtop");
    addVartoStore("Lepton2Fromtop_pt");
    addVartoStore("Lepton2Fromtop_eta");
    addVartoStore("Lepton2Fromtop_phi");
    addVartoStore("Lepton2Fromtop_mass");
    addVartoStore("Lepton2Fromtop_transverse_energy");
    addVartoStore("bJet1Fromtop");
    addVartoStore("bJet1Fromtop_pt");
    addVartoStore("bJet1Fromtop_eta");
    addVartoStore("bJet1Fromtop_phi");
    addVartoStore("bJet1Fromtop_mass");
    addVartoStore("bJet1Fromtop_transverse_energy");
    // addVartoStore("Min_deltaR_bJetFromtop");
    // addVartoStore("cleanjet4Rvecs_withoutbjet");
    addVartoStore("cleanjet4Rvecs_withoutbjet_pt");
    addVartoStore("cleanjet4Rvecs_withoutbjet_eta");
    addVartoStore("cleanjet4Rvecs_withoutbjet_phi");
    addVartoStore("cleanjet4Rvecs_withoutbjet_mass");
    addVartoStore("cleanjet4Rvecs_withoutbjet_transverse_energy");
    addVartoStore("cleanjet4Rvecs_withoutbjet_number");
    addVartoStore("cleanjet4vecs_withoutbjet");
    // addVartoStore("dis_comb_jets_withoutbjet");
    // addVartoStore("Selected_jet1");
    addVartoStore("Selected_jet1_pt");
    // addVartoStore("Selected_jet2");
    addVartoStore("Selected_jet2_pt");
    // addVartoStore("JetsFromHiggs4Rvecs");
    addVartoStore("JetsFromHiggs_pt");
    addVartoStore("JetsFromHiggs_eta");
    addVartoStore("JetsFromHiggs_phi");
    addVartoStore("JetsFromHiggs_mass");
    addVartoStore("JetsFromHiggs_energy");
    addVartoStore("JetsFromHiggs_transverse_energy");
    addVartoStore("JetsFromHiggs_deltaeta");
    addVartoStore("JetsFromHiggs_deltaphi");
    addVartoStore("JetsFromHiggs_deltaR");
    addVartoStore("JetsFromHiggs_transverse_mass");
    addVartoStore("Min_deltaR_JetsFromHiggs");
    // addVartoStore("Jet1FromHiggs");
    addVartoStore("Jet1FromHiggs_pt");
    addVartoStore("Jet1FromHiggs_eta");
    addVartoStore("Jet1FromHiggs_phi");
    addVartoStore("Jet1FromHiggs_mass");
    addVartoStore("Jet1FromHiggs_transverse_energy");
    // addVartoStore("Jet2FromHiggs");
    addVartoStore("Jet2FromHiggs_pt");
    addVartoStore("Jet2FromHiggs_eta");
    addVartoStore("Jet2FromHiggs_phi");
    addVartoStore("Jet2FromHiggs_mass");
    addVartoStore("Jet2FromHiggs_transverse_energy");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_pt");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_eta");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_phi");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_mass");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_transverse_energy");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_deltaeta");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_deltaphi");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_deltaR");
    addVartoStore("Jet1FromHiggs_Jet2FromHiggs_transverse_mass");

    // // addVartoStore("bJet1Fromtop_genJetIdx");
    // // addVartoStore("sel_bjetgentop_gencheck");
    // addVartoStore("Jet1FromHiggs_genJetIdx");
    // // addVartoStore("sel_jet1genHiggs_gencheck");
    // addVartoStore("Jet2FromHiggs_genJetIdx");
    // // addVartoStore("sel_jet2genHiggs_gencheck");
    // addVartoStore("sel_lepton1genHiggs_gencheck");
    // addVartoStore("sel_lepton2gentop_gencheck");

    addVartoStore("Lepton2Fromtop_bJet1Fromtop");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_pt");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_eta");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_phi");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_mass");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_transverse_energy");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_deltaeta");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_deltaphi");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_deltaR");
    addVartoStore("Lepton2Fromtop_bJet1Fromtop_transverse_mass");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_pt");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_eta");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_phi");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_mass");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_transverse_energy");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_deltaeta");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_deltaphi");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_deltaR");
    addVartoStore("Lepton2Fromtop_Jet1FromHiggs_transverse_mass");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_pt");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_eta");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_phi");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_mass");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_transverse_energy");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_deltaeta");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_deltaphi");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_deltaR");
    addVartoStore("Lepton2Fromtop_Jet2FromHiggs_transverse_mass");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_pt");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_eta");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_phi");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_mass");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_transverse_energy");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_deltaeta");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_deltaphi");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_deltaR");
    addVartoStore("Lepton1FromHiggs_bJet1Fromtop_transverse_mass");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_pt");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_eta");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_phi");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_mass");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_transverse_energy");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_deltaeta");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_deltaphi");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_deltaR");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_transverse_mass");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_pt");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_eta");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_phi");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_mass");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_transverse_energy");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_deltaeta");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_deltaphi");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_deltaR");
    addVartoStore("Lepton1FromHiggs_Jet2FromHiggs_transverse_mass");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_pt");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_eta");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_phi");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_mass");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_transverse_energy");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_deltaeta");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_deltaphi");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_deltaR");
    addVartoStore("Lepton1FromHiggs_Jet1FromHiggs_Jet2FromHiggs_transverse_mass");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_pt");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_eta");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_phi");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_px");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_py");
    if(!_isData)
    {
        addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_mass");
    }
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_transverse_energy");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_deltaeta");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_deltaphi");
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_deltaR");
    if(!_isData)
    {
        addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_transverse_mass");
    }
    addVartoStore("L2bJ1Fromtop_L1J1J2FromHiggs_Relative_St");
    addVartoStore("Tprime_transverse_mass");
}

void TprimeAnalyser::defineMoreReconstructedVars()
{
    // _rlm = _rlm.Define("GenParticles_mass","GenPart_mass")
    //            .Define("GenParticles_pdgId","GenPart_pdgId")
    //            .Define("MotherGenParticles","GenPart_genPartIdxMother")
    //            .Define("MotherGenParticlesPdgId","ROOT::VecOps::Take(GenPart_pdgId,MotherGenParticles)");

    // _rlm = _rlm.Define("deltaR_Selected_leptonfromtop_Selected_clean_jet",::DR_2,{"LeptonFromtop_eta","Selected_clean_jet_eta","LeptonFromtop_phi","Selected_clean_jet_phi"})
    //            .Define("deltaR_Selected_leptonfromtop_Selected_clean_bjet",::DR_2,{"LeptonFromtop_eta","Selected_clean_bjet_eta","LeptonFromtop_phi","Selected_clean_bjet_phi"})
    //            .Define("deltaR_Selected_leptonfromHiggs_Selected_clean_jet",::DR_2,{"LeptonFromHiggs_eta","Selected_clean_jet_eta","LeptonFromHiggs_phi","Selected_clean_jet_phi"})
    //            .Define("deltaR_Selected_leptonfromHiggs_Selected_clean_bjet",::DR_2,{"LeptonFromHiggs_eta","Selected_clean_bjet_eta","LeptonFromHiggs_phi","Selected_clean_bjet_phi"});

    // _rlm = _rlm.Define("Pt_leptonfromtop_over_Pt_clean_jet",::ratio_2,{"LeptonFromtop_pt","Selected_clean_jet_pt"})
    //            .Define("Pt_leptonfromtop_over_Pt_clean_bjet",::ratio_2,{"LeptonFromtop_pt","Selected_clean_bjet_pt"})
    //            .Define("Pt_leptonfromHiggs_over_Pt_clean_jet",::ratio_2,{"LeptonFromHiggs_pt","Selected_clean_jet_pt"})
    //            .Define("Pt_leptonfromHiggs_over_Pt_clean_bjet",::ratio_2,{"LeptonFromHiggs_pt","Selected_clean_bjet_pt"});

    // _rlm = _rlm.Define("Min_deltaR_Selected_leptonfromtop_Selected_clean_jet",::minDR_6,{"LeptonFromtop_eta","Selected_clean_jet_eta","LeptonFromtop_phi","Selected_clean_jet_phi"})
    //            .Define("Min_deltaR_Selected_leptonfromtop_Selected_clean_bjet",::minDR_6,{"LeptonFromtop_eta","Selected_clean_bjet_eta","LeptonFromtop_phi","Selected_clean_bjet_phi"})
    //            .Define("Min_deltaR_Selected_leptonfromHiggs_Selected_clean_jet",::minDR_6,{"LeptonFromHiggs_eta","Selected_clean_jet_eta","LeptonFromHiggs_phi","Selected_clean_jet_phi"})
    //            .Define("Min_deltaR_Selected_leptonfromHiggs_Selected_clean_bjet",::minDR_6,{"LeptonFromHiggs_eta","Selected_clean_bjet_eta","LeptonFromHiggs_phi","Selected_clean_bjet_phi"});

    // _rlm = _rlm.Define("Min_Pt_leptonfromtop_over_Pt_clean_jet","LeptonFromtop_pt/Selected_clean_jet_pt[0]")
    //            .Define("Min_Pt_leptonfromtop_over_Pt_clean_bjet","LeptonFromtop_pt/Selected_clean_bjet_pt[0]")
    //            .Define("Min_Pt_leptonfromHiggs_over_Pt_clean_jet","LeptonFromHiggs_pt/Selected_clean_jet_pt[0]")
    //            .Define("Min_Pt_leptonfromHiggs_over_Pt_clean_bjet","LeptonFromHiggs_pt/Selected_clean_bjet_pt[0]");

    _rlm = _rlm.Define("Max_deltaR_top",::maxDR,{"LeptonFromtop_eta","bJetFromtop_eta","LeptonFromtop_phi","bJetFromtop_phi"})
               .Define("bJetFromtop1","p4_bJetFromtop[Max_deltaR_top]")
               .Define("bJetFromtop1_pt","bJetFromtop1.Pt()")
               .Define("bJetFromtop1_eta","bJetFromtop1.Eta()")
               .Define("bJetFromtop1_phi","bJetFromtop1.Phi()")
               .Define("bJetFromtop1_mass","bJetFromtop1.M()")
               .Define("bJetFromtop1_transverse_energy","sqrt(bJetFromtop1_pt*bJetFromtop1_pt + bJetFromtop1_mass*bJetFromtop1_mass)");

    _rlm = _rlm.Define("Min_deltaR_Higgs",::minDR_2,{"JetFromHiggs_eta","JetFromHiggs_eta","JetFromHiggs_phi","JetFromHiggs_phi"})
               .Define("JetFromHiggs1","p4_JetFromHiggs[Min_deltaR_Higgs[0]]")
               .Define("JetFromHiggs1_pt","JetFromHiggs1.Pt()")
               .Define("JetFromHiggs1_eta","JetFromHiggs1.Eta()")
               .Define("JetFromHiggs1_phi","JetFromHiggs1.Phi()")
               .Define("JetFromHiggs1_mass","JetFromHiggs1.M()")
               .Define("JetFromHiggs1_transverse_energy","sqrt(JetFromHiggs1_pt*JetFromHiggs1_pt + JetFromHiggs1_mass*JetFromHiggs1_mass)")
               .Define("JetFromHiggs2","p4_JetFromHiggs[Min_deltaR_Higgs[1]]")
               .Define("JetFromHiggs2_pt","JetFromHiggs2.Pt()")
               .Define("JetFromHiggs2_eta","JetFromHiggs2.Eta()")
               .Define("JetFromHiggs2_phi","JetFromHiggs2.Phi()")
               .Define("JetFromHiggs2_mass","JetFromHiggs2.M()")
               .Define("JetFromHiggs2_transverse_energy","sqrt(JetFromHiggs2_pt*JetFromHiggs2_pt + JetFromHiggs2_mass*JetFromHiggs2_mass)");

    // _rlm = _rlm.Define("Min_deltaR_Higgs",::minDR_2,{"JetFromtopHiggs_eta","JetFromtopHiggs_eta","JetFromtopHiggs_phi","JetFromtopHiggs_phi","JetFromtopHiggs_pt","bJetFromtop1_pt"})
    //            .Define("JetFromHiggs1","p4_JetFromtopHiggs[Min_deltaR_Higgs[0]]")
    //            .Define("JetFromHiggs1_pt","JetFromHiggs1.Pt()")
    //            .Define("JetFromHiggs1_eta","JetFromHiggs1.Eta()")
    //            .Define("JetFromHiggs1_phi","JetFromHiggs1.Phi()")
    //            .Define("JetFromHiggs1_mass","JetFromHiggs1.M()")
    //            .Define("JetFromHiggs2","p4_JetFromtopHiggs[Min_deltaR_Higgs[1]]")
    //            .Define("JetFromHiggs2_pt","JetFromHiggs2.Pt()")
    //            .Define("JetFromHiggs2_eta","JetFromHiggs2.Eta()")
    //            .Define("JetFromHiggs2_phi","JetFromHiggs2.Phi()")
    //            .Define("JetFromHiggs2_mass","JetFromHiggs2.M()");

    _rlm = _rlm.Define("LeptonFromtop_bJetFromtop1","LeptonFromtop4vecs + bJetFromtop1")
               .Define("LeptonFromtop_bJetFromtop1_pt","LeptonFromtop_bJetFromtop1.Pt()")
               .Define("LeptonFromtop_bJetFromtop1_eta","LeptonFromtop_bJetFromtop1.Eta()")
               .Define("LeptonFromtop_bJetFromtop1_phi","LeptonFromtop_bJetFromtop1.Phi()")
               .Define("LeptonFromtop_bJetFromtop1_mass","LeptonFromtop_bJetFromtop1.M()")
               .Define("LeptonFromtop_bJetFromtop1_transverse_energy","sqrt(LeptonFromtop_bJetFromtop1_pt*LeptonFromtop_bJetFromtop1_pt + LeptonFromtop_bJetFromtop1_mass*LeptonFromtop_bJetFromtop1_mass)")
               .Define("LeptonFromtop_bJetFromtop1_deltaeta","abs(LeptonFromtop_eta - bJetFromtop1_eta)")
               .Define("LeptonFromtop_bJetFromtop1_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtop_phi,bJetFromtop1_phi))")
               .Define("LeptonFromtop_bJetFromtop1_deltaR","ROOT::VecOps::DeltaR(LeptonFromtop_eta,bJetFromtop1_eta,LeptonFromtop_phi,bJetFromtop1_phi)")
               .Define("LeptonFromtop_bJetFromtop1_transverse_mass","sqrt(LeptonFromtop_mass*LeptonFromtop_mass + bJetFromtop1_mass*bJetFromtop1_mass + 2*LeptonFromtop_transverse_energy*bJetFromtop1_transverse_energy - 2*LeptonFromtop_pt*bJetFromtop1_pt*cos(LeptonFromtop_bJetFromtop1_deltaphi))")
    //            .Define("WFromtop_transverse_mass","sqrt(LeptonFromtop_pt*MET_pt*(1 - cos(LeptonFromtop_phi-(MET_phi+LeptonFromtop_phi-LeptonFromHiggs_phi)/2)))")
               .Define("LeptonFromtop_JetFromHiggs1","LeptonFromtop4vecs + JetFromHiggs1")
               .Define("LeptonFromtop_JetFromHiggs1_pt","LeptonFromtop_JetFromHiggs1.Pt()")
               .Define("LeptonFromtop_JetFromHiggs1_eta","LeptonFromtop_JetFromHiggs1.Eta()")
               .Define("LeptonFromtop_JetFromHiggs1_phi","LeptonFromtop_JetFromHiggs1.Phi()")
               .Define("LeptonFromtop_JetFromHiggs1_mass","LeptonFromtop_JetFromHiggs1.M()")
               .Define("LeptonFromtop_JetFromHiggs1_transverse_energy","sqrt(LeptonFromtop_JetFromHiggs1_pt*LeptonFromtop_JetFromHiggs1_pt + LeptonFromtop_JetFromHiggs1_mass*LeptonFromtop_JetFromHiggs1_mass)")
               .Define("LeptonFromtop_JetFromHiggs1_deltaeta","abs(LeptonFromtop_eta - JetFromHiggs1_eta)")
               .Define("LeptonFromtop_JetFromHiggs1_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtop_phi,JetFromHiggs1_phi))")
               .Define("LeptonFromtop_JetFromHiggs1_deltaR","ROOT::VecOps::DeltaR(LeptonFromtop_eta,JetFromHiggs1_eta,LeptonFromtop_phi,JetFromHiggs1_phi)")
               .Define("LeptonFromtop_JetFromHiggs1_transverse_mass","sqrt(LeptonFromtop_mass*LeptonFromtop_mass + JetFromHiggs1_mass*JetFromHiggs1_mass + 2*LeptonFromtop_transverse_energy*JetFromHiggs1_transverse_energy - 2*LeptonFromtop_pt*JetFromHiggs1_pt*cos(LeptonFromtop_JetFromHiggs1_deltaphi))")
               .Define("LeptonFromtop_JetFromHiggs2","LeptonFromtop4vecs + JetFromHiggs2")
               .Define("LeptonFromtop_JetFromHiggs2_pt","LeptonFromtop_JetFromHiggs2.Pt()")
               .Define("LeptonFromtop_JetFromHiggs2_eta","LeptonFromtop_JetFromHiggs2.Eta()")
               .Define("LeptonFromtop_JetFromHiggs2_phi","LeptonFromtop_JetFromHiggs2.Phi()")
               .Define("LeptonFromtop_JetFromHiggs2_mass","LeptonFromtop_JetFromHiggs2.M()")
               .Define("LeptonFromtop_JetFromHiggs2_transverse_energy","sqrt(LeptonFromtop_JetFromHiggs2_pt*LeptonFromtop_JetFromHiggs2_pt + LeptonFromtop_JetFromHiggs2_mass*LeptonFromtop_JetFromHiggs2_mass)")
               .Define("LeptonFromtop_JetFromHiggs2_deltaeta","abs(LeptonFromtop_eta - JetFromHiggs2_eta)")
               .Define("LeptonFromtop_JetFromHiggs2_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtop_phi,JetFromHiggs2_phi))")
               .Define("LeptonFromtop_JetFromHiggs2_deltaR","ROOT::VecOps::DeltaR(LeptonFromtop_eta,JetFromHiggs2_eta,LeptonFromtop_phi,JetFromHiggs2_phi)")
               .Define("LeptonFromtop_JetFromHiggs2_transverse_mass","sqrt(LeptonFromtop_mass*LeptonFromtop_mass + JetFromHiggs2_mass*JetFromHiggs2_mass + 2*LeptonFromtop_transverse_energy*JetFromHiggs2_transverse_energy - 2*LeptonFromtop_pt*JetFromHiggs2_pt*cos(LeptonFromtop_JetFromHiggs2_deltaphi))")
               .Define("LeptonFromHiggs_bJetFromtop1","LeptonFromHiggs4vecs + bJetFromtop1")
               .Define("LeptonFromHiggs_bJetFromtop1_pt","LeptonFromHiggs_bJetFromtop1.Pt()")
               .Define("LeptonFromHiggs_bJetFromtop1_eta","LeptonFromHiggs_bJetFromtop1.Eta()")
               .Define("LeptonFromHiggs_bJetFromtop1_phi","LeptonFromHiggs_bJetFromtop1.Phi()")
               .Define("LeptonFromHiggs_bJetFromtop1_mass","LeptonFromHiggs_bJetFromtop1.M()")
               .Define("LeptonFromHiggs_bJetFromtop1_transverse_energy","sqrt(LeptonFromHiggs_bJetFromtop1_pt*LeptonFromHiggs_bJetFromtop1_pt + LeptonFromHiggs_bJetFromtop1_mass*LeptonFromHiggs_bJetFromtop1_mass)")
               .Define("LeptonFromHiggs_bJetFromtop1_deltaeta","abs(LeptonFromHiggs_eta - bJetFromtop1_eta)")
               .Define("LeptonFromHiggs_bJetFromtop1_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromHiggs_phi,bJetFromtop1_phi))")
               .Define("LeptonFromHiggs_bJetFromtop1_deltaR","ROOT::VecOps::DeltaR(LeptonFromHiggs_eta,bJetFromtop1_eta,LeptonFromHiggs_phi,bJetFromtop1_phi)")
               .Define("LeptonFromHiggs_bJetFromtop1_transverse_mass","sqrt(LeptonFromHiggs_mass*LeptonFromHiggs_mass + bJetFromtop1_mass*bJetFromtop1_mass + 2*LeptonFromHiggs_transverse_energy*bJetFromtop1_transverse_energy - 2*LeptonFromHiggs_pt*bJetFromtop1_pt*cos(LeptonFromHiggs_bJetFromtop1_deltaphi))")
               .Define("LeptonFromHiggs_JetFromHiggs1","LeptonFromHiggs4vecs + JetFromHiggs1")
               .Define("LeptonFromHiggs_JetFromHiggs1_pt","LeptonFromHiggs_JetFromHiggs1.Pt()")
               .Define("LeptonFromHiggs_JetFromHiggs1_eta","LeptonFromHiggs_JetFromHiggs1.Eta()")
               .Define("LeptonFromHiggs_JetFromHiggs1_phi","LeptonFromHiggs_JetFromHiggs1.Phi()")
               .Define("LeptonFromHiggs_JetFromHiggs1_mass","LeptonFromHiggs_JetFromHiggs1.M()")
               .Define("LeptonFromHiggs_JetFromHiggs1_transverse_energy","sqrt(LeptonFromHiggs_JetFromHiggs1_pt*LeptonFromHiggs_JetFromHiggs1_pt + LeptonFromHiggs_JetFromHiggs1_mass*LeptonFromHiggs_JetFromHiggs1_mass)")
               .Define("LeptonFromHiggs_JetFromHiggs1_deltaeta","abs(LeptonFromHiggs_eta - JetFromHiggs1_eta)")
               .Define("LeptonFromHiggs_JetFromHiggs1_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromHiggs_phi,JetFromHiggs1_phi))")
               .Define("LeptonFromHiggs_JetFromHiggs1_deltaR","ROOT::VecOps::DeltaR(LeptonFromHiggs_eta,JetFromHiggs1_eta,LeptonFromHiggs_phi,JetFromHiggs1_phi)")
               .Define("LeptonFromHiggs_JetFromHiggs1_transverse_mass","sqrt(LeptonFromHiggs_mass*LeptonFromHiggs_mass + JetFromHiggs1_mass*JetFromHiggs1_mass + 2*LeptonFromHiggs_transverse_energy*JetFromHiggs1_transverse_energy - 2*LeptonFromHiggs_pt*JetFromHiggs1_pt*cos(LeptonFromHiggs_JetFromHiggs1_deltaphi))")
               .Define("LeptonFromHiggs_JetFromHiggs2","LeptonFromHiggs4vecs + JetFromHiggs2")
               .Define("LeptonFromHiggs_JetFromHiggs2_pt","LeptonFromHiggs_JetFromHiggs2.Pt()")
               .Define("LeptonFromHiggs_JetFromHiggs2_eta","LeptonFromHiggs_JetFromHiggs2.Eta()")
               .Define("LeptonFromHiggs_JetFromHiggs2_phi","LeptonFromHiggs_JetFromHiggs2.Phi()")
               .Define("LeptonFromHiggs_JetFromHiggs2_mass","LeptonFromHiggs_JetFromHiggs2.M()")
               .Define("LeptonFromHiggs_JetFromHiggs2_transverse_energy","sqrt(LeptonFromHiggs_JetFromHiggs2_pt*LeptonFromHiggs_JetFromHiggs2_pt + LeptonFromHiggs_JetFromHiggs2_mass*LeptonFromHiggs_JetFromHiggs2_mass)")
               .Define("LeptonFromHiggs_JetFromHiggs2_deltaeta","abs(LeptonFromHiggs_eta - JetFromHiggs2_eta)")
               .Define("LeptonFromHiggs_JetFromHiggs2_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromHiggs_phi,JetFromHiggs2_phi))")
               .Define("LeptonFromHiggs_JetFromHiggs2_deltaR","ROOT::VecOps::DeltaR(LeptonFromHiggs_eta,JetFromHiggs2_eta,LeptonFromHiggs_phi,JetFromHiggs2_phi)")
               .Define("LeptonFromHiggs_JetFromHiggs2_transverse_mass","sqrt(LeptonFromHiggs_mass*LeptonFromHiggs_mass + JetFromHiggs2_mass*JetFromHiggs2_mass + 2*LeptonFromHiggs_transverse_energy*JetFromHiggs2_transverse_energy - 2*LeptonFromHiggs_pt*JetFromHiggs2_pt*cos(LeptonFromHiggs_JetFromHiggs2_deltaphi))")
               .Define("JetFromHiggs1_JetFromHiggs2","JetFromHiggs1 + JetFromHiggs2")
               .Define("JetFromHiggs1_JetFromHiggs2_pt","JetFromHiggs1_JetFromHiggs2.Pt()")
               .Define("JetFromHiggs1_JetFromHiggs2_eta","JetFromHiggs1_JetFromHiggs2.Eta()")
               .Define("JetFromHiggs1_JetFromHiggs2_phi","JetFromHiggs1_JetFromHiggs2.Phi()")
               .Define("JetFromHiggs1_JetFromHiggs2_mass","JetFromHiggs1_JetFromHiggs2.M()")
               .Define("JetFromHiggs1_JetFromHiggs2_transverse_energy","sqrt(JetFromHiggs1_JetFromHiggs2_pt*JetFromHiggs1_JetFromHiggs2_pt + JetFromHiggs1_JetFromHiggs2_mass*JetFromHiggs1_JetFromHiggs2_mass)")
               .Define("JetFromHiggs1_JetFromHiggs2_deltaeta","abs(JetFromHiggs2_eta - JetFromHiggs1_eta)")
               .Define("JetFromHiggs1_JetFromHiggs2_deltaphi","abs(ROOT::VecOps::DeltaPhi(JetFromHiggs1_phi,JetFromHiggs2_phi))")
               .Define("JetFromHiggs1_JetFromHiggs2_deltaR","ROOT::VecOps::DeltaR(JetFromHiggs1_eta,JetFromHiggs2_eta,JetFromHiggs1_phi,JetFromHiggs2_phi)")
               .Define("JetFromHiggs1_JetFromHiggs2_transverse_mass","sqrt(JetFromHiggs1_mass*JetFromHiggs1_mass + JetFromHiggs2_mass*JetFromHiggs2_mass + 2*JetFromHiggs1_transverse_energy*JetFromHiggs2_transverse_energy - 2*JetFromHiggs1_pt*JetFromHiggs2_pt*cos(JetFromHiggs1_JetFromHiggs2_deltaphi))")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2","LeptonFromHiggs4vecs + JetFromHiggs1_JetFromHiggs2")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_pt","LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2.Pt()")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_eta","LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2.Eta()")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_phi","LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2.Phi()")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_mass","LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2.M()")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_energy","sqrt(LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_pt*LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_pt + LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_mass*LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_mass)")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaeta","abs(LeptonFromHiggs_eta - JetFromHiggs1_JetFromHiggs2_eta)")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromHiggs_phi,JetFromHiggs1_JetFromHiggs2_phi))")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaR","ROOT::VecOps::DeltaR(LeptonFromHiggs_eta,JetFromHiggs1_JetFromHiggs2_eta,LeptonFromHiggs_phi,JetFromHiggs1_JetFromHiggs2_phi)")
               .Define("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_mass","sqrt(LeptonFromHiggs_mass*LeptonFromHiggs_mass + JetFromHiggs1_JetFromHiggs2_mass*JetFromHiggs1_JetFromHiggs2_mass + 2*LeptonFromHiggs_transverse_energy*JetFromHiggs1_JetFromHiggs2_transverse_energy - 2*LeptonFromHiggs_pt*LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_pt*cos(LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaphi))")
               .Define("LbJFromtop1_LJJFromHiggs12","LeptonFromtop_bJetFromtop1 + LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2")
               .Define("LbJFromtop1_LJJFromHiggs12_pt","LbJFromtop1_LJJFromHiggs12.Pt()")
               .Define("LbJFromtop1_LJJFromHiggs12_eta","LbJFromtop1_LJJFromHiggs12.Eta()")
               .Define("LbJFromtop1_LJJFromHiggs12_phi","LbJFromtop1_LJJFromHiggs12.Phi()")
               .Define("LbJFromtop1_LJJFromHiggs12_mass","LbJFromtop1_LJJFromHiggs12.M()")
               .Define("LbJFromtop1_LJJFromHiggs12_transverse_energy","sqrt(LbJFromtop1_LJJFromHiggs12_pt*LbJFromtop1_LJJFromHiggs12_pt + LbJFromtop1_LJJFromHiggs12_mass*LbJFromtop1_LJJFromHiggs12_mass)")
               .Define("LbJFromtop1_LJJFromHiggs12_deltaeta","abs(LeptonFromtop_bJetFromtop1_eta - LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_eta)")
               .Define("LbJFromtop1_LJJFromHiggs12_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtop_bJetFromtop1_phi,LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_phi))")
               .Define("LbJFromtop1_LJJFromHiggs12_deltaR","ROOT::VecOps::DeltaR(LeptonFromtop_bJetFromtop1_eta,LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_eta,LeptonFromtop_bJetFromtop1_phi,LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_phi)")
               .Define("LbJFromtop1_LJJFromHiggs12_transverse_mass","sqrt(LeptonFromtop_bJetFromtop1_mass*LeptonFromtop_bJetFromtop1_mass + LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_mass*LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_mass + 2*LeptonFromtop_bJetFromtop1_transverse_energy*LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_energy - 2*LeptonFromtop_bJetFromtop1_pt*LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_pt*cos(LbJFromtop1_LJJFromHiggs12_deltaphi))");
               
    _rlm = _rlm.Define("LeptonFromtop_NeutrinoFromtop","LeptonFromtop4vecs + p4_GenNeutrinoFromtop_first")
               .Define("LeptonFromtop_NeutrinoFromtop_pt","LeptonFromtop_NeutrinoFromtop.Pt()")      
               .Define("LeptonFromtop_NeutrinoFromtop_eta","LeptonFromtop_NeutrinoFromtop.Eta()") 
               .Define("LeptonFromtop_NeutrinoFromtop_phi","LeptonFromtop_NeutrinoFromtop.Phi()") 
               .Define("LeptonFromtop_NeutrinoFromtop_mass","LeptonFromtop_NeutrinoFromtop.M()") 
               .Define("LeptonFromtop_NeutrinoFromtop_transverse_energy","sqrt(LeptonFromtop_NeutrinoFromtop_pt*LeptonFromtop_NeutrinoFromtop_pt + LeptonFromtop_NeutrinoFromtop_mass*LeptonFromtop_NeutrinoFromtop_mass)")
               .Define("LeptonFromtop_NeutrinoFromtop_deltaeta","abs(LeptonFromtop_eta - GenNeutrinoFromtop_first_eta)")
               .Define("LeptonFromtop_NeutrinoFromtop_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtop_phi,GenNeutrinoFromtop_first_phi))")
               .Define("LeptonFromtop_NeutrinoFromtop_deltaR","ROOT::VecOps::DeltaR(LeptonFromtop_eta,GenNeutrinoFromtop_first_eta,LeptonFromtop_phi,GenNeutrinoFromtop_first_phi)")
               .Define("LeptonFromtop_NeutrinoFromtop_transverse_mass","sqrt(2*LeptonFromtop_pt*GenNeutrinoFromtop_first_pt*(1 - cos(LeptonFromtop_NeutrinoFromtop_deltaphi)))")
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1","LeptonFromtop_NeutrinoFromtop + bJetFromtop1")
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_pt","LeptonFromtop_NeutrinoFromtop_bJetFromtop1.Pt()")      
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_eta","LeptonFromtop_NeutrinoFromtop_bJetFromtop1.Eta()") 
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_phi","LeptonFromtop_NeutrinoFromtop_bJetFromtop1.Phi()") 
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_mass","LeptonFromtop_NeutrinoFromtop_bJetFromtop1.M()") 
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_transverse_energy","sqrt(LeptonFromtop_NeutrinoFromtop_bJetFromtop1_pt*LeptonFromtop_NeutrinoFromtop_bJetFromtop1_pt + LeptonFromtop_NeutrinoFromtop_bJetFromtop1_mass*LeptonFromtop_NeutrinoFromtop_bJetFromtop1_mass)")
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_deltaeta","abs(LeptonFromtop_NeutrinoFromtop_eta - bJetFromtop1_eta)")
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtop_NeutrinoFromtop_phi,bJetFromtop1_phi))")
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_deltaR","ROOT::VecOps::DeltaR(LeptonFromtop_NeutrinoFromtop_eta,bJetFromtop1_eta,LeptonFromtop_NeutrinoFromtop_phi,bJetFromtop1_phi)")
               .Define("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_transverse_mass","sqrt(LeptonFromtop_NeutrinoFromtop_mass*LeptonFromtop_NeutrinoFromtop_mass + bJetFromtop1_mass*bJetFromtop1_mass + 2*LeptonFromtop_NeutrinoFromtop_transverse_energy*bJetFromtop1_transverse_energy - 2*LeptonFromtop_NeutrinoFromtop_pt*bJetFromtop1_pt*cos(LeptonFromtop_NeutrinoFromtop_bJetFromtop1_deltaphi))")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs","LeptonFromHiggs4vecs + p4_GenNeutrinoFromHiggs_first")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_pt","LeptonFromHiggs_NeutrinoFromHiggs.Pt()")      
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_eta","LeptonFromHiggs_NeutrinoFromHiggs.Eta()") 
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_phi","LeptonFromHiggs_NeutrinoFromHiggs.Phi()") 
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_mass","LeptonFromHiggs_NeutrinoFromHiggs.M()")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_transverse_energy","sqrt(LeptonFromHiggs_NeutrinoFromHiggs_pt*LeptonFromHiggs_NeutrinoFromHiggs_pt + LeptonFromHiggs_NeutrinoFromHiggs_mass*LeptonFromHiggs_NeutrinoFromHiggs_mass)")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_deltaeta","abs(LeptonFromHiggs_eta - GenNeutrinoFromHiggs_first_eta)")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromHiggs_phi,GenNeutrinoFromHiggs_first_phi))")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_deltaR","ROOT::VecOps::DeltaR(LeptonFromHiggs_eta,GenNeutrinoFromHiggs_first_eta,LeptonFromHiggs_phi,GenNeutrinoFromHiggs_first_phi)")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_transverse_mass","sqrt(2*LeptonFromHiggs_pt*GenNeutrinoFromHiggs_first_pt*(1 - cos(LeptonFromHiggs_NeutrinoFromHiggs_deltaphi)))")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2","LeptonFromHiggs_NeutrinoFromHiggs + JetFromHiggs1_JetFromHiggs2")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_pt","LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2.Pt()")      
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_eta","LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2.Eta()") 
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_phi","LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2.Phi()") 
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_mass","LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2.M()") 
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_energy","sqrt(LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_pt*LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_pt + LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_mass*LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_mass)")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaeta","abs(LeptonFromHiggs_NeutrinoFromHiggs_eta - JetFromHiggs1_JetFromHiggs2_eta)")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromHiggs_NeutrinoFromHiggs_phi,JetFromHiggs1_JetFromHiggs2_phi))")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaR","ROOT::VecOps::DeltaR(LeptonFromHiggs_NeutrinoFromHiggs_eta,JetFromHiggs1_JetFromHiggs2_eta,LeptonFromHiggs_NeutrinoFromHiggs_phi,JetFromHiggs1_JetFromHiggs2_phi)")
               .Define("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_mass","sqrt(LeptonFromHiggs_NeutrinoFromHiggs_mass*LeptonFromHiggs_NeutrinoFromHiggs_mass + JetFromHiggs1_JetFromHiggs2_mass*JetFromHiggs1_JetFromHiggs2_mass + 2*LeptonFromHiggs_NeutrinoFromHiggs_transverse_energy*JetFromHiggs1_JetFromHiggs2_transverse_energy - 2*LeptonFromHiggs_NeutrinoFromHiggs_pt*JetFromHiggs1_JetFromHiggs2_pt*cos(LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaphi))")
               .Define("LNbJFromtop1_LNJJFromHiggs12","LeptonFromtop_NeutrinoFromtop_bJetFromtop1 + LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2")
               .Define("LNbJFromtop1_LNJJFromHiggs12_pt","LNbJFromtop1_LNJJFromHiggs12.Pt()")
               .Define("LNbJFromtop1_LNJJFromHiggs12_eta","LNbJFromtop1_LNJJFromHiggs12.Eta()")
               .Define("LNbJFromtop1_LNJJFromHiggs12_phi","LNbJFromtop1_LNJJFromHiggs12.Phi()")
               .Define("LNbJFromtop1_LNJJFromHiggs12_mass","LNbJFromtop1_LNJJFromHiggs12.M()")
               .Define("LNbJFromtop1_LNJJFromHiggs12_transverse_energy","sqrt(LNbJFromtop1_LNJJFromHiggs12_pt*LNbJFromtop1_LNJJFromHiggs12_pt + LNbJFromtop1_LNJJFromHiggs12_mass*LNbJFromtop1_LNJJFromHiggs12_mass)")
               .Define("LNbJFromtop1_LNJJFromHiggs12_deltaeta","abs(LeptonFromtop_NeutrinoFromtop_bJetFromtop1_eta - LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_eta)")
               .Define("LNbJFromtop1_LNJJFromHiggs12_deltaphi","abs(ROOT::VecOps::DeltaPhi(LeptonFromtop_NeutrinoFromtop_bJetFromtop1_phi,LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_phi))")
               .Define("LNbJFromtop1_LNJJFromHiggs12_deltaR","ROOT::VecOps::DeltaR(LeptonFromtop_NeutrinoFromtop_bJetFromtop1_eta,LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_eta,LeptonFromtop_NeutrinoFromtop_bJetFromtop1_phi,LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_phi)")
               .Define("LNbJFromtop1_LNJJFromHiggs12_transverse_mass","sqrt(LeptonFromtop_NeutrinoFromtop_bJetFromtop1_mass*LeptonFromtop_NeutrinoFromtop_bJetFromtop1_mass + LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_mass*LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_mass + 2*LeptonFromtop_NeutrinoFromtop_bJetFromtop1_transverse_energy*LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_energy - 2*LeptonFromtop_NeutrinoFromtop_bJetFromtop1_pt*LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_pt*cos(LNbJFromtop1_LNJJFromHiggs12_deltaphi))")
               .Define("LNbJFromtop1_LNJJFromHiggs12_Relative_St","(LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_pt + LeptonFromtop_NeutrinoFromtop_bJetFromtop1_pt) / St_with_MET");

    // _rlm = _rlm.Define("METFromNeutrinos_pt","sqrt((NeutrinoFromtop_pt*cos(NeutrinoFromtop_phi)+NeutrinoFromHiggs_pt*cos(NeutrinoFromHiggs_phi))*(NeutrinoFromtop_pt*cos(NeutrinoFromtop_phi)+NeutrinoFromHiggs_pt*cos(NeutrinoFromHiggs_phi)) + (NeutrinoFromtop_pt*sin(NeutrinoFromtop_phi)+NeutrinoFromHiggs_pt*sin(NeutrinoFromHiggs_phi))*(NeutrinoFromtop_pt*sin(NeutrinoFromtop_phi)+NeutrinoFromHiggs_pt*sin(NeutrinoFromHiggs_phi)))")
    //            .Define("METFromNeutrinos_phi","tan((NeutrinoFromtop_pt*cos(NeutrinoFromtop_phi)+NeutrinoFromHiggs_pt*cos(NeutrinoFromHiggs_phi)) / (NeutrinoFromtop_pt*sin(NeutrinoFromtop_phi)+NeutrinoFromHiggs_pt*sin(NeutrinoFromHiggs_phi)))")
    //            .Define("LeptonFromtop_NeutrinoFromtop_transverse_mass","sqrt(2*LeptonFromtop_pt*NeutrinoFromtop_pt*(1-cos(LeptonFromtop_phi-NeutrinoFromtop_phi)))")
    //            .Define("LeptonFromHiggs_NeutrinoFromHiggs_transverse_mass","sqrt(2*LeptonFromHiggs_pt*NeutrinoFromHiggs_pt*(1-cos(LeptonFromHiggs_phi-NeutrinoFromHiggs_phi)))");

    // addVartoStore("GenParticles_mass");
    // addVartoStore("GenParticles_pdgId");
    // addVartoStore("MotherGenParticles");
    // addVartoStore("MotherGenParticlesPdgId");
    // addVartoStore("GenJets_pt");
    // addVartoStore("GenJets_eta");
    // addVartoStore("GenJets_phi");
    // addVartoStore("GenJets_mass");
    // addVartoStore("GenJets_hadronFlavour");
    // addVartoStore("GenJets_partonFlavour");
    // addVartoStore("GenJets_light");
    // addVartoStore("GenJets_light_pt");
    // addVartoStore("GenJets_light_eta");
    // addVartoStore("GenJets_light_phi");
    // addVartoStore("GenJets_light_mass");
    // addVartoStore("GenJets_light_hadronFlavour");
    // addVartoStore("GenJets_light_partonFlavour");

    // //reconstructed leptons
    // // addVartoStore("sel_leptongentopHiggs");
    // addVartoStore("LeptonFromtopHiggs_pt");
    // addVartoStore("LeptonFromtopHiggs_leading_pt");
    // addVartoStore("LeptonFromtopHiggs_subleading_pt");
    // addVartoStore("LeptonFromtopHiggs_sum_two_muons_pt");
    // addVartoStore("LeptonFromtopHiggs_sum_all_muons_pt");  
    // addVartoStore("LeptonFromtopHiggs_eta");
    // addVartoStore("LeptonFromtopHiggs_phi");
    // addVartoStore("LeptonFromtopHiggs_mass");
    // addVartoStore("LeptonFromtopHiggs_transverse_energy");
    // addVartoStore("LeptonFromtopHiggs_charge");
    // addVartoStore("LeptonFromtopHiggs_charge_sum");
    // addVartoStore("LeptonFromtopHiggs_number");
    // addVartoStore("LeptonFromtopHiggs_deltaeta");
    // addVartoStore("LeptonFromtopHiggs_deltaphi");
    // addVartoStore("LeptonFromtopHiggs_deltaR");
    // addVartoStore("LeptonFromtopHiggs_miniPFRelIso_all");
    // addVartoStore("LeptonFromtopHiggs_miniPFRelIso_chg");
    // addVartoStore("LeptonFromtopHiggs_pfRelIso03_all");
    // addVartoStore("LeptonFromtopHiggs_pfRelIso03_chg");
    // addVartoStore("LeptonFromtopHiggs_genPartIdx");
    // addVartoStore("LeptonFromtopHiggs_genPartFlav");
    // addVartoStore("LeptonFromtopHiggs_pdgId"); 
    // addVartoStore("p4_LeptonFromtopHiggs");
    // // addVartoStore("sel_leptongentoporHiggs");
    // addVartoStore("LeptonFromtop4vecs");
    // addVartoStore("LeptonFromtop_pt");
    // addVartoStore("LeptonFromtop_eta");
    // addVartoStore("LeptonFromtop_phi");
    // addVartoStore("LeptonFromtop_mass");
    // addVartoStore("LeptonFromtop_transverse_energy");
    // addVartoStore("LeptonFromtop_energy");
    // addVartoStore("LeptonFromtop_px");
    // addVartoStore("LeptonFromtop_py");
    // addVartoStore("LeptonFromtop_pz");
    // addVartoStore("LeptonFromHiggs4vecs");
    // addVartoStore("LeptonFromHiggs_pt");
    // addVartoStore("LeptonFromHiggs_eta");
    // addVartoStore("LeptonFromHiggs_phi");
    // addVartoStore("LeptonFromHiggs_mass");
    // addVartoStore("LeptonFromHiggs_transverse_energy");
    // addVartoStore("LeptonFromHiggs_energy");
    // addVartoStore("LeptonFromHiggs_px");
    // addVartoStore("LeptonFromHiggs_py");
    // addVartoStore("LeptonFromHiggs_pz");
    // addVartoStore("LeptonFromtopHiggs_ordered_pt");
    // addVartoStore("LeptonFromtopHiggs_ordered_px");
    // addVartoStore("LeptonFromtopHiggs_ordered_py");
    // addVartoStore("LeptonFromtopHiggs_ordered_pz");
    // addVartoStore("LeptonFromtopHiggs_ordered_energy");
    // addVartoStore("LeptonFromtopHiggs_ordered_mass");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs_pt");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs_eta");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs_phi");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs_mass");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs_transverse_energy");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs_deltaeta");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs_deltaphi");
    // addVartoStore("LeptonFromtop_LeptonFromHiggs_deltaR");
    // // addVartoStore("sel_leptongenW");
    // addVartoStore("LeptonFromW_pt");
    // addVartoStore("LeptonFromW_leading_pt");
    // addVartoStore("LeptonFromW_subleading_pt");
    // addVartoStore("LeptonFromW_sum_two_muons_pt");
    // addVartoStore("LeptonFromW_sum_all_muons_pt");  
    // addVartoStore("LeptonFromW_eta");
    // addVartoStore("LeptonFromW_phi");
    // addVartoStore("LeptonFromW_mass");
    // addVartoStore("LeptonFromW_charge");
    // addVartoStore("LeptonFromW_charge_sum");
    // addVartoStore("LeptonFromW_number");
    // addVartoStore("LeptonFromW_deltaeta");
    // addVartoStore("LeptonFromW_deltaphi");
    // addVartoStore("LeptonFromW_deltaR");
    // addVartoStore("LeptonFromW_genPartIdx");
    // addVartoStore("LeptonFromW_genPartFlav");
    // addVartoStore("LeptonFromW_pdgId"); 
    // addVartoStore("p4_LeptonFromW");
    // addVartoStore("sel_muongenW");
    addVartoStore("MuonFromW_pt");
    addVartoStore("MuonFromW_charge");
    addVartoStore("MuonFromW_number");
    // addVartoStore("sel_electrongenW");
    addVartoStore("ElectronFromW_pt");
    addVartoStore("ElectronFromW_charge");
    addVartoStore("ElectronFromW_number");
    // addVartoStore("sel_leptongenMother");
    addVartoStore("MotherfromLepton_pdgId");
    addVartoStore("MotherfromLepton1_pdgId");
    addVartoStore("MotherfromLepton2_pdgId");
    addVartoStore("MotherfromLepton12_pdgId");
    // addVartoStore("GenLeptons");
    // addVartoStore("MotherGenLeptons");
    // addVartoStore("MotherGenLeptonsPdgId");
    // addVartoStore("IsFromW");
    // addVartoStore("GenMuons");
    // addVartoStore("GenMuonsPdgId");
    // addVartoStore("MotherGenMuons");
    // addVartoStore("MotherGenMuonsPdgId");
    // addVartoStore("IsMuonFromW");
    // addVartoStore("GenElectrons");
    // addVartoStore("GenElectronsPdgId");
    // addVartoStore("MotherGenElectrons");
    // addVartoStore("MotherGenElectronsPdgId");
    // addVartoStore("IsElectronFromW");
    // addVartoStore("GenW");
    // addVartoStore("MotherGenW");
    // addVartoStore("MotherGenWPdgId");
    // addVartoStore("IsFromHiggs");
    // addVartoStore("IsFromtop");
    // addVartoStore("GenLeptons");
    // addVartoStore("GenLeptonsPdgId");
    // addVartoStore("GenLeptons_pt");
    // addVartoStore("GenLeptons_eta");
    // addVartoStore("GenLeptons_phi");

    // // reconstructed neutrinos
    // addVartoStore("GenNeutrinos");
    // addVartoStore("GenNeutrinosPdgId");
    // addVartoStore("GenNeutrinos_pt");
    // addVartoStore("GenNeutrinos_eta");
    // addVartoStore("GenNeutrinos_phi");
    // addVartoStore("MotherGenNeutrinos");
    // addVartoStore("MotherGenNeutrinosPdgId");
    // addVartoStore("IsNeutrinoFromW");
    // addVartoStore("GenNeutrinoFromWPdgId");
    // addVartoStore("GenNeutrinoFromW_pt");
    // addVartoStore("GenNeutrinoFromW_eta");
    // addVartoStore("GenNeutrinoFromW_phi");
    // addVartoStore("GenNeutrinoFromW_number");
    // addVartoStore("MotherNeutrinoFromW");
    // addVartoStore("MotherNeutrinoFromWPdgId");
    // addVartoStore("MotherNeutrinoFromWMass");
    // // addVartoStore("sel_NeutrinoFromWgenHiggs");
    // // addVartoStore("sel_NeutrinoFromWgentop");
    // addVartoStore("GenNeutrinoFromtop_pt");
    // addVartoStore("GenNeutrinoFromtop_eta");
    // addVartoStore("GenNeutrinoFromtop_phi");
    // addVartoStore("GenNeutrinoFromtop_px");
    // addVartoStore("GenNeutrinoFromtop_py");
    // addVartoStore("p4_GenNeutrinoFromtop");
    // addVartoStore("p4_GenNeutrinoFromtop_first");
    // addVartoStore("GenNeutrinoFromtop_first_pt");
    // addVartoStore("GenNeutrinoFromtop_first_eta");
    // addVartoStore("GenNeutrinoFromtop_first_phi");
    // addVartoStore("GenNeutrinoFromHiggs_pt");
    // addVartoStore("GenNeutrinoFromHiggs_eta");
    // addVartoStore("GenNeutrinoFromHiggs_phi");
    // addVartoStore("GenNeutrinoFromHiggs_px");
    // addVartoStore("GenNeutrinoFromHiggs_py");
    // addVartoStore("p4_GenNeutrinoFromHiggs");
    // addVartoStore("p4_GenNeutrinoFromHiggs_first");
    // addVartoStore("GenNeutrinoFromHiggs_first_pt");
    // addVartoStore("GenNeutrinoFromHiggs_first_eta");
    // addVartoStore("GenNeutrinoFromHiggs_first_phi");
    // addVartoStore("MotherNeutrinoFromHiggsPdgid");
    // addVartoStore("MotherNeutrinoFromHiggsMass");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs_pt");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs_eta");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs_phi");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs_mass");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs_transverse_energy");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs_deltaeta");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs_deltaphi");
    // addVartoStore("NeutrinoFromtop_NeutrinoFromHiggs_deltaR");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop_pt");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop_eta");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop_phi");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop_mass");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop_transverse_energy");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop_deltaeta");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop_deltaphi");
    // addVartoStore("NeutrinoFromtop_LeptonFromtop_deltaR");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs_pt");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs_eta");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs_phi");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs_mass");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs_transverse_energy");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs_deltaeta");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs_deltaphi");
    // addVartoStore("NeutrinoFromHiggs_LeptonFromHiggs_deltaR");

    // // reconstructed jets
    // addVartoStore("JFromtopHiggs_pt");
    // addVartoStore("sel_jetgentopHiggs");
    // addVartoStore("JetFromtopHiggs_pt");
    // addVartoStore("JetFromtopHiggs_leading_pt");
    // addVartoStore("JetFromtopHiggs_subleading_pt");
    // addVartoStore("JetFromtopHiggs_subsubleading_pt");
    // addVartoStore("JetFromtopHiggs_Ht");
    // addVartoStore("JetFromtopHiggs_eta");
    // addVartoStore("JetFromtopHiggs_phi");
    // addVartoStore("JetFromtopHiggs_mass");
    // addVartoStore("JetFromtopHiggs_transverse_energy");
    // addVartoStore("JetFromtopHiggs_number");
    // addVartoStore("JetFromtopHiggs_btagDeepB");
    // addVartoStore("JetFromtopHiggs_btagDeepFlavB");
    // addVartoStore("JetFromtopHiggs_genJetIdx");
    // addVartoStore("JetFromtopHiggs_partonFlavour");
    // addVartoStore("JetFromtopHiggs_hadronFlavour");
    // addVartoStore("p4_JetFromtopHiggs");
    // addVartoStore("Rp4_JetFromtopHiggs");
    // addVartoStore("JetFromtopHiggs_energy");
    // addVartoStore("JetFromtopHiggs_px");
    // addVartoStore("JetFromtopHiggs_py");
    // addVartoStore("JetFromtopHiggs_pz");
    // addVartoStore("sel_bjetgentopHiggs");
    // addVartoStore("bJetFromtopHiggs_pt");
    // addVartoStore("bJetFromtopHiggs_leading_pt");
    // addVartoStore("bJetFromtopHiggs_subleading_pt");
    // addVartoStore("bJetFromtopHiggs_subsubleading_pt");
    // addVartoStore("bJetFromtopHiggs_Ht");
    // addVartoStore("bJetFromtopHiggs_eta");
    // addVartoStore("bJetFromtopHiggs_phi");
    // addVartoStore("bJetFromtopHiggs_mass");
    // addVartoStore("bJetFromtopHiggs_transverse_energy");
    // addVartoStore("bJetFromtopHiggs_number");
    // addVartoStore("bJetFromtopHiggs_btagDeepB");
    // addVartoStore("bJetFromtopHiggs_btagDeepFlavB");
    // addVartoStore("bJetFromtopHiggs_genJetIdx");
    // addVartoStore("bJetFromtopHiggs_partonFlavour");
    // addVartoStore("bJetFromtopHiggs_hadronFlavour");
    // addVartoStore("p4_bJetFromtopHiggs");
    // addVartoStore("Rp4_bJetFromtopHiggs");
    // addVartoStore("bJetFromtopHiggs_energy");
    // addVartoStore("bJetFromtopHiggs_px");
    // addVartoStore("bJetFromtopHiggs_py");
    // addVartoStore("bJetFromtopHiggs_pz");
    // // addVartoStore("JFromHiggs_pt");
    // // addVartoStore("sel_jetgenHiggs");
    // addVartoStore("JetFromHiggs_pt");
    // addVartoStore("JetFromHiggs_leading_pt");
    // addVartoStore("JetFromHiggs_subleading_pt");
    // addVartoStore("JetFromHiggs_subsubleading_pt");
    // addVartoStore("JetFromHiggs_Ht");
    // addVartoStore("JetFromHiggs_eta");
    // addVartoStore("JetFromHiggs_phi");
    // addVartoStore("JetFromHiggs_mass");
    // addVartoStore("JetFromHiggs_number");
    // addVartoStore("JetFromHiggs_transverse_energy");
    // addVartoStore("JetFromHiggs_btagDeepB");
    // addVartoStore("JetFromHiggs_btagDeepFlavB");
    // addVartoStore("JetFromHiggs_genJetIdx");
    // addVartoStore("JetFromHiggs_partonFlavour");
    // addVartoStore("JetFromHiggs_hadronFlavour");
    // // addVartoStore("p4_JetFromHiggs");
    // // addVartoStore("Rp4_JetFromHiggs");
    // addVartoStore("JetFromHiggs_energy");
    // addVartoStore("JetFromHiggs_px");
    // addVartoStore("JetFromHiggs_py");
    // addVartoStore("JetFromHiggs_pz");
    // // addVartoStore("JFromtop_pt");
    // // addVartoStore("sel_bjetgentop");
    // addVartoStore("bJetFromtop_pt");
    // addVartoStore("bJetFromtop_leading_pt");
    // addVartoStore("bJetFromtop_subleading_pt");
    // addVartoStore("bJetFromtop_subsubleading_pt");
    // addVartoStore("bJetFromtop_Ht");
    // addVartoStore("bJetFromtop_eta");
    // addVartoStore("bJetFromtop_phi");
    // addVartoStore("bJetFromtop_mass");
    // addVartoStore("bJetFromtop_transverse_energy");
    // addVartoStore("bJetFromtop_number");
    // addVartoStore("bJetFromtop_btagDeepB");
    // addVartoStore("bJetFromtop_btagDeepFlavB");
    // addVartoStore("bJetFromtop_genJetIdx");
    // addVartoStore("bJetFromtop_partonFlavour");
    // addVartoStore("bJetFromtop_hadronFlavour");
    // // // addVartoStore("p4_bJetFromtop");
    // addVartoStore("JetFromtop_px");
    // addVartoStore("JetFromtop_py");
    // addVartoStore("JetFromtop_pz");
    // // // addVartoStore("p4_bJetFromtop_first");
    // // addVartoStore("bJetFromtop_first_pt");
    // // addVartoStore("bJetFromtop_first_eta");
    // // addVartoStore("bJetFromtop_first_phi");
    // // addVartoStore("bJetFromtop_first_mass");
    // // addVartoStore("bJetFromtop_first_energy");
    // // addVartoStore("bJetFromtop_first_px");
    // // addVartoStore("bJetFromtop_first_py");
    // // addVartoStore("bJetFromtop_first_pz");
    // // // addVartoStore("Rp4_bJetFromtop_first");
    // addVartoStore("Max_deltaR_top");
    // // addVartoStore("Min_deltaR_top");
    // addVartoStore("bJetFromtop1");
    // addVartoStore("bJetFromtop1_pt");
    // addVartoStore("bJetFromtop1_eta");
    // addVartoStore("bJetFromtop1_phi");
    // addVartoStore("bJetFromtop1_mass");
    // addVartoStore("bJetFromtop1_transverse_energy");   
    // // addVartoStore("Min_deltaR_Higgs");
    // addVartoStore("JetFromHiggs1");
    // addVartoStore("JetFromHiggs1_pt");
    // addVartoStore("JetFromHiggs1_eta");
    // addVartoStore("JetFromHiggs1_phi");
    // addVartoStore("JetFromHiggs1_mass");
    // addVartoStore("JetFromHiggs1_transverse_energy");
    // addVartoStore("JetFromHiggs2");
    // addVartoStore("JetFromHiggs2_pt");
    // addVartoStore("JetFromHiggs2_eta");
    // addVartoStore("JetFromHiggs2_phi");
    // addVartoStore("JetFromHiggs2_mass");
    // addVartoStore("JetFromHiggs2_transverse_energy");

    // //combination of reconstructed particles
    // // addVartoStore("deltaR_Selected_leptonfromtop_Selected_clean_jet");
    // // addVartoStore("deltaR_Selected_leptonfromtop_Selected_clean_bjet");
    // // addVartoStore("deltaR_Selected_leptonfromHiggs_Selected_clean_jet");
    // // addVartoStore("deltaR_Selected_leptonfromHiggs_Selected_clean_bjet");
    // // addVartoStore("Pt_leptonfromtop_over_Pt_clean_jet");
    // // addVartoStore("Pt_leptonfromtop_over_Pt_clean_bjet");
    // // addVartoStore("Pt_leptonfromHiggs_over_Pt_clean_jet");
    // // addVartoStore("Pt_leptonfromHiggs_over_Pt_clean_bjet");
    // // addVartoStore("Min_deltaR_Selected_leptonfromtop_Selected_clean_jet");
    // // addVartoStore("Min_deltaR_Selected_leptonfromtop_Selected_clean_bjet");
    // // addVartoStore("Min_deltaR_Selected_leptonfromHiggs_Selected_clean_jet");
    // // addVartoStore("Min_deltaR_Selected_leptonfromHiggs_Selected_clean_bjet");
    // // addVartoStore("Min_Pt_leptonfromtop_over_Pt_clean_jet");
    // // addVartoStore("Min_Pt_leptonfromtop_over_Pt_clean_bjet");
    // // addVartoStore("Min_Pt_leptonfromHiggs_over_Pt_clean_jet");
    // // addVartoStore("Min_Pt_leptonfromHiggs_over_Pt_clean_bjet");
    // addVartoStore("LeptonFromtop_bJetFromtop1");
    // addVartoStore("LeptonFromtop_bJetFromtop1_pt");
    // addVartoStore("LeptonFromtop_bJetFromtop1_eta");
    // addVartoStore("LeptonFromtop_bJetFromtop1_phi");
    // addVartoStore("LeptonFromtop_bJetFromtop1_mass");
    // addVartoStore("LeptonFromtop_bJetFromtop1_transverse_energy");
    // addVartoStore("LeptonFromtop_bJetFromtop1_deltaeta");
    // addVartoStore("LeptonFromtop_bJetFromtop1_deltaphi");
    // addVartoStore("LeptonFromtop_bJetFromtop1_deltaR");
    // addVartoStore("LeptonFromtop_bJetFromtop1_transverse_mass");
    // // addVartoStore("WFromtop_transverse_mass");
    // addVartoStore("LeptonFromtop_JetFromHiggs1");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_pt");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_eta");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_phi");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_mass");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_transverse_energy");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_deltaeta");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_deltaphi");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_deltaR");
    // addVartoStore("LeptonFromtop_JetFromHiggs1_transverse_mass");
    // addVartoStore("LeptonFromtop_JetFromHiggs2");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_pt");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_eta");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_phi");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_mass");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_transverse_energy");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_deltaeta");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_deltaphi");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_deltaR");
    // addVartoStore("LeptonFromtop_JetFromHiggs2_transverse_mass");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_pt");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_eta");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_phi");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_mass");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_transverse_energy");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_deltaeta");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_deltaphi");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_deltaR");
    // addVartoStore("LeptonFromHiggs_bJetFromtop1_transverse_mass");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_pt");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_eta");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_phi");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_mass");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_transverse_energy");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_deltaeta");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_deltaphi");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_deltaR");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_transverse_mass");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_pt");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_eta");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_phi");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_mass");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_transverse_energy");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_deltaeta");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_deltaphi");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_deltaR"); 
    // addVartoStore("LeptonFromHiggs_JetFromHiggs2_transverse_mass");     
    // addVartoStore("JetFromHiggs1_JetFromHiggs2");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_pt");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_eta");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_phi");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_mass");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_transverse_energy");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_deltaeta");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_deltaphi");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_deltaR");
    // addVartoStore("JetFromHiggs1_JetFromHiggs2_transverse_mass");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_pt");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_eta");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_phi");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_mass");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_energy");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaeta");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaphi");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaR");
    // addVartoStore("LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_mass");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_pt");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_eta");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_phi");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_mass");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_transverse_energy");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_deltaeta");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_deltaphi");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_deltaR");
    // addVartoStore("LbJFromtop1_LJJFromHiggs12_transverse_mass");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_pt");      
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_eta"); 
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_phi"); 
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_mass"); 
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_transverse_energy");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_deltaeta");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_deltaphi");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_deltaR");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_transverse_mass");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_pt");      
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_eta"); 
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_phi"); 
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_mass"); 
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_transverse_energy");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_deltaeta");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_deltaphi");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_deltaR");
    // addVartoStore("LeptonFromtop_NeutrinoFromtop_bJetFromtop1_transverse_mass");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_pt");      
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_eta"); 
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_phi"); 
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_mass");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_transverse_energy");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_deltaeta");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_deltaphi");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_deltaR");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_transverse_mass");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_pt");      
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_eta"); 
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_phi"); 
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_mass"); 
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_transverse_energy");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaeta");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaphi");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaR");
    // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_JetFromHiggs1_JetFromHiggs1_transverse_mass");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_pt");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_eta");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_phi");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_mass");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_transverse_energy");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_deltaeta");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_deltaphi");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_deltaR");
    // addVartoStore("LNbJFromtop1_LNJJFromHiggs12_transverse_mass");

    // // addVartoStore("METFromNeutrinos_pt");
    // // addVartoStore("METFromNeutrinos_phi");
    // // addVartoStore("LeptonFromtop_NeutrinoFromtop_transverse_mass");
    // // addVartoStore("LeptonFromHiggs_NeutrinoFromHiggs_transverse_mass");

}

void TprimeAnalyser::bookHists(std::string Region)
{
    if (debug)
    {
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

    add1DHist( {"hnevents","Number of Events; Number of Events;Events;", 2, -0.5, 1.5}, "one","evWeight","");
    if(!_isData)
    {
        add1DHist( {"genweight","Genweight of Events; Genweight of Events;Events;", 40, -400.0, 400.0}, "genWeight","one","");
    }
    
    ////=== electrons =====//
    if(Region=="CR_TT_SL")
    {
        add1DHist( {"Number_Electrons_Loose","Number of loose electrons; Number of loose electrons;Events;", 5, 0.0, 5.0}, "Selected_electron_looseTTSL_number","evWeight","0");
        add1DHist( {"Pt_Electrons_Loose","Pt of loose electrons; Pt of the loose electrons (GeV);Events;", 20, 0.0, 400.0}, "Selected_electron_looseTTSL_pt","evWeight","");
        add1DHist( {"Leading_Pt_Electrons_Loose","Pt of leading loose electron; Pt of the leading loose electron (GeV);Events;", 20, 0.0, 400.0}, "Selected_electron_looseTTSL_leading_pt","evWeight","");
        add1DHist( {"Subleading_Pt_Electrons_Loose","Pt of subleading loose electron; Pt of the subleading loose electron (GeV);Events;", 10, 0.0, 200.0}, "Selected_electron_looseTTSL_subleading_pt","evWeight","");
        add1DHist( {"Eta_Electrons_Loose","Eta of loose electrons; Eta of the loose electrons;Events;", 30, -3.0, 3.0}, "Selected_electron_looseTTSL_eta","evWeight","");
        add1DHist( {"Leading_Eta_Electrons_Loose","Eta of loose leading electron; Eta of the leading loose electron (GeV);Events;", 30, -3.0, 3.0}, "Selected_electron_looseTTSL_leading_eta","evWeight","");
        add1DHist( {"Subleading_Eta_Electrons_Loose","Eta of loose subleading electron; Eta of the subleading loose electron (GeV);Events;", 30, -3.0, 3.0}, "Selected_electron_looseTTSL_subleading_eta","evWeight","");
        add1DHist( {"Phi_Electrons_Loose","Phi of loose electrons; Phi of the loose electrons;Events;", 35, -3.5, 3.5}, "Selected_electron_looseTTSL_phi","evWeight","");
        add1DHist( {"Charge_Electrons_Loose","Charge of loose electrons; Charge of the loose electrons;Events;", 5, -2.0, 3.0}, "Selected_electron_looseTTSL_charge","evWeight","");
        add1DHist( {"miniPFRelIso_all_Electrons_Loose","miniPFRelIso_all isolation of loose electrons; miniPFRel_all isolation of the loose electrons;Events;", 40, 0.0, 0.4}, "Selected_electron_looseTTSL_miniPFRelIso_all","evWeight","");
        add1DHist( {"Sip_Electrons_Loose","Impact parameter significance of loose electrons; Impact parameter significance of the loose electrons wrt PV;Events;", 10, 0.0, 2.0}, "Selected_electron_looseTTSL_sip3d","evWeight","");
    }

    add1DHist( {"Number_Electrons","Number of  electrons; Number of electrons;Events;", 5, 0.0, 5.0}, "Selected_electron_number","evWeight","0");
    add1DHist( {"Pt_Electrons","Pt of electrons; Pt of the electrons (GeV);Events;", 20, 0.0, 400.0}, "Selected_electron_pt","evWeight","");
    add1DHist( {"Leading_Pt_Electrons","Pt of leading electron; Pt of the leading electron (GeV);Events;", 20, 0.0, 400.0}, "Selected_electron_leading_pt","evWeight","");
    add1DHist( {"Subleading_Pt_Electrons","Pt of subleading electron; Pt of the subleading electron (GeV);Events;", 10, 0.0, 200.0}, "Selected_electron_subleading_pt","evWeight","");
    add1DHist( {"Eta_Electrons","Eta of electrons; Eta of the electrons;Events;", 30, -3.0, 3.0}, "Selected_electron_eta","evWeight","");
    add1DHist( {"Leading_Eta_Electrons","Eta of leading electron; Eta of the leading electron (GeV);Events;", 30, -3.0, 3.0}, "Selected_electron_leading_eta","evWeight","");
    add1DHist( {"Subleading_Eta_Electrons","Eta of subleading electron; Eta of the subleading electron (GeV);Events;", 30, -3.0, 3.0}, "Selected_electron_subleading_eta","evWeight","");
    add1DHist( {"Phi_Electrons","Phi of electrons; Phi of the electrons;Events;", 35, -3.5, 3.5}, "Selected_electron_phi","evWeight","");
    add1DHist( {"Charge_Electrons","Charge of electrons; Charge of the electrons;Events;", 5, -2.0, 3.0}, "Selected_electron_charge","evWeight","");
    add1DHist( {"miniPFRelIso_all_Electrons","miniPFRelIso_all isolation of electrons; miniPFRel_all isolation of the electrons;Events;", 40, 0.0, 0.4}, "Selected_electron_miniPFRelIso_all","evWeight","");
    add1DHist( {"Sip_Electrons","Impact parameter significance of electrons; Impact parameter significance of the electrons wrt PV;Events;", 10, 0.0, 2.0}, "Selected_electron_sip3d","evWeight","");

    ////=== muons =====//
    if(Region=="CR_TT_SL")
    {
        add1DHist( {"Number_Muons_Loose","Number of loose muons; Number of loose muons;Events;", 5, 0.0, 5.0}, "Selected_muon_looseTTSL_number","evWeight","0");
        add1DHist( {"Pt_Muons_Loose","Pt of loose muons; Pt of the loose muons (GeV);Events;", 20, 0.0, 400.0}, "Selected_muon_looseTTSL_pt","evWeight","");
        add1DHist( {"Leading_Pt_Muons_Loose","Pt of leading loose muon; Pt of the leading loose muon (GeV);Events;", 20, 0.0, 400.0}, "Selected_muon_looseTTSL_leading_pt","evWeight","");
        add1DHist( {"Subleading_Pt_Muons_Loose","Pt of subleading loose muon; Pt of the subleading loose muon (GeV);Events;", 10, 0.0, 200.0}, "Selected_muon_looseTTSL_subleading_pt","evWeight","");
        add1DHist( {"Eta_Muons_Loose","Eta of loose muons; Eta of the loose muons;Events;", 30, -3.0, 3.0}, "Selected_muon_looseTTSL_eta","evWeight","");
        add1DHist( {"Leading_Eta_Muons_Loose","Eta of leading loose muon; Eta of the leading loose muon (GeV);Events;", 30, -3.0, 3.0}, "Selected_muon_looseTTSL_leading_eta","evWeight","");
        add1DHist( {"Subleading_Eta_Muons_Loose","Eta of subleading loose muon; Eta of the subleading loose muon (GeV);Events;", 30, -3.0, 3.0}, "Selected_muon_looseTTSL_subleading_eta","evWeight","");
        add1DHist( {"Phi_Muons_Loose","Phi of loose muons; Phi of the loose muons;Events;", 35, -3.5, 3.5}, "Selected_muon_looseTTSL_phi","evWeight","");
        add1DHist( {"Charge_Muons_Loose","Charge of loose muons; Charge of the loose muons;Events;", 5, -2.0, 3.0}, "Selected_muon_looseTTSL_charge","evWeight","");
        add1DHist( {"miniPFRelIso_all_Muons_Loose","miniPFRelIso_all isolation of loose muons; miniPFRel_all isolation of the loose muons;Events;", 40, 0.0, 0.4}, "Selected_muon_looseTTSL_miniPFRelIso_all","evWeight","");
        add1DHist( {"Sip_Muons_Loose","Impact parameter significance of loose muons; Impact parameter significance of the muons wrt PV;Events;", 15, 0.0, 3.0}, "Selected_muon_looseTTSL_sip3d","evWeight","");
    }
    
    add1DHist( {"Number_Muons","Number of muons; Number of muons;Events;", 5, 0.0, 5.0}, "Selected_muon_number","evWeight","0");
    add1DHist( {"Pt_Muons","Pt of muons; Pt of the muons (GeV);Events;", 20, 0.0, 400.0}, "Selected_muon_pt","evWeight","");
    add1DHist( {"Leading_Pt_Muons","Pt of leading muon; Pt of the leading muon (GeV);Events;", 20, 0.0, 400.0}, "Selected_muon_leading_pt","evWeight","");
    add1DHist( {"Subleading_Pt_Muons","Pt of subleading muon; Pt of the subleading muon (GeV);Events;", 10, 0.0, 200.0}, "Selected_muon_subleading_pt","evWeight","");
    add1DHist( {"Eta_Muons","Eta of muons; Eta of the muons;Events;", 30, -3.0, 3.0}, "Selected_muon_eta","evWeight","");
    add1DHist( {"Leading_Eta_Muons","Eta of leading muon; Eta of the leading muon (GeV);Events;", 30, -3.0, 3.0}, "Selected_muon_leading_eta","evWeight","");
    add1DHist( {"Subleading_Eta_Muons","Eta of subleading muon; Eta of the subleading muon (GeV);Events;", 30, -3.0, 3.0}, "Selected_muon_subleading_eta","evWeight","");
    add1DHist( {"Phi_Muons","Phi of muons; Phi of the muons;Events;", 35, -3.5, 3.5}, "Selected_muon_phi","evWeight","");
    add1DHist( {"Charge_Muons","Charge of muons; Charge of the muons;Events;", 5, -2.0, 3.0}, "Selected_muon_charge","evWeight","");
    add1DHist( {"miniPFRelIso_all_Muons","miniPFRelIso_all isolation of muons; miniPFRel_all isolation of the muons;Events;", 40, 0.0, 0.4}, "Selected_muon_miniPFRelIso_all","evWeight","");
    add1DHist( {"Sip_Muons","Impact parameter significance of muons; Impact parameter significance of the muons wrt PV;Events;", 15, 0.0, 3.0}, "Selected_muon_sip3d","evWeight","");

    ////=== leptons =====//
    if(Region=="CR_TT_SL")
    {
        add1DHist( {"Number_Leptons_Loose","Number of loose leptons; Number of loose leptons;Events;", 5, 0.0, 5.0}, "Selected_lepton_looseTTSL_number","evWeight","0");
        add1DHist( {"Pt_Leptons_Loose","Pt of loose leptons; Pt of the loose leptons (GeV);Events;", 20, 0.0, 400.0}, "Selected_lepton_looseTTSL_pt","evWeight","");
        add1DHist( {"Leading_Pt_Leptons_Loose","Pt of leading loose lepton; Pt of the leading loose lepton (GeV);Events;", 20, 0.0, 400.0}, "Selected_lepton_looseTTSL_leading_pt","evWeight","0");
        add1DHist( {"Subleading_Pt_Leptons_Loose","Pt of subleading loose lepton; Pt of the subleading loose lepton (GeV);Events;", 10, 0.0, 200.0}, "Selected_lepton_looseTTSL_subleading_pt","evWeight","0");
        add1DHist( {"Sum_Pt_Two_Leptons_Loose","Sum of Pt of two loose leptons; Sum of the Pt of the two loose leptons (GeV);Events;", 30, 0.0, 600.0}, "Selected_lepton_looseTTSL_sum_two_leptons_pt","evWeight","0");
        add1DHist( {"Eta_Leptons_Loose","Eta of loose leptons; Eta of the loose leptons;Events;", 30, -3.0, 3.0}, "Selected_lepton_looseTTSL_eta","evWeight","");
        add1DHist( {"Leading_Eta_Leptons_Loose","Eta of leading loose lepton; Eta of the leading loose lepton (GeV);Events;", 30, -3.0, 3.0}, "Selected_lepton_looseTTSL_leading_eta","evWeight","0");
        add1DHist( {"Subleading_Eta_Leptons_Loose","Eta of subleading loose lepton; Eta of the subleading loose lepton (GeV);Events;", 30, -3.0, 3.0}, "Selected_lepton_looseTTSL_subleading_eta","evWeight","0");
        add1DHist( {"Phi_Leptons_Loose","Phi of loose leptons; Phi of the loose leptons;Events;", 35, -3.5, 3.5}, "Selected_lepton_looseTTSL_phi","evWeight","");
        add1DHist( {"Charge_Leptons_Loose","Charge of loose leptons; Charge of the loose leptons;Events;", 5, -2.0, 3.0}, "Selected_lepton_looseTTSL_charge","evWeight","");
        add1DHist( {"Deltaeta_Leptons_Loose","Delta eta between 2 loose leptons; Delta eta between the two loose leptons;Events;", 20, 0.0, 4.0}, "Selected_lepton_looseTTSL_deltaeta","evWeight","0");
        add1DHist( {"Deltaphi_Leptons_Loose","Delta phi between 2 loose leptons; Delta phi between the two loose leptons;Events;", 20, 0.0, 4.0}, "Selected_lepton_looseTTSL_deltaphi","evWeight","0");
        add1DHist( {"DeltaR_Leptons_Loose","Delta R between 2 loose leptons; Delta R between the two loose leptons;Events;", 30, 0.0, 6.0}, "Selected_lepton_looseTTSL_deltaR","evWeight","0");
        add1DHist( {"miniPFRelIso_all_Leptons_Loose","miniPFRelIso_all isolation of loose leptons; miniPFRel_all isolation of the loose leptons;Events;", 40, 0.0, 0.4}, "Selected_lepton_looseTTSL_miniPFRelIso_all","evWeight","");
        add1DHist( {"Sip_Leptons_Loose","Impact parameter significance of loose leptons; Impact parameter significance of the loose leptons wrt PV;Events;", 15, 0.0, 3.0}, "Selected_lepton_looseTTSL_sip3d","evWeight","");
        add1DHist( {"Pt_Sum_Two_Leptons_Loose","Pt of sum of two loose leptons; Pt of the sum of the two loose leptons (GeV);Events;", 30, 0.0, 600.0}, "Vectorial_sum_two_leptons_looseTTSL_pt","evWeight","0");
        add1DHist( {"Mass_Sum_Two_Leptons_Loose","Invariant mass of two loose leptons; Invariant mass of the two loose leptons (GeV);Events;", 20, 0.0, 800.0}, "Vectorial_sum_two_leptons_looseTTSL_mass","evWeight","0");
    }
    
    add1DHist( {"Number_Leptons","Number of leptons; Number of leptons;Events;", 5, 0.0, 5.0}, "Selected_lepton_number","evWeight","0");
    add1DHist( {"Pt_Leptons","Pt of leptons; Pt of the leptons (GeV);Events;", 20, 0.0, 400.0}, "Selected_lepton_pt","evWeight","");
    add1DHist( {"Leading_Pt_Leptons","Pt of leading lepton; Pt of the leading lepton (GeV);Events;", 20, 0.0, 400.0}, "Selected_lepton_leading_pt","evWeight","0");
    add1DHist( {"Subleading_Pt_Leptons","Pt of subleading lepton; Pt of the subleading lepton (GeV);Events;", 10, 0.0, 200.0}, "Selected_lepton_subleading_pt","evWeight","0");
    add1DHist( {"Sum_Pt_Two_Leptons","Sum of Pt of two leptons; Sum of the Pt of the two leptons (GeV);Events;", 30, 0.0, 600.0}, "Selected_lepton_sum_two_leptons_pt","evWeight","0");
    add1DHist( {"Eta_Leptons","Eta of leptons; Eta of the leptons;Events;", 30, -3.0, 3.0}, "Selected_lepton_eta","evWeight","");
    add1DHist( {"Leading_Eta_Leptons","Eta of leading lepton; Eta of the leading lepton (GeV);Events;", 30, -3.0, 3.0}, "Selected_lepton_leading_eta","evWeight","0");
    add1DHist( {"Subleading_Eta_Leptons","Eta of subleading lepton; Eta of the subleading lepton (GeV);Events;", 30, -3.0, 3.0}, "Selected_lepton_subleading_eta","evWeight","0");
    add1DHist( {"Phi_Leptons","Phi of leptons; Phi of the leptons;Events;", 35, -3.5, 3.5}, "Selected_lepton_phi","evWeight","");
    add1DHist( {"Charge_Leptons","Charge of leptons; Charge of the leptons;Events;", 5, -2.0, 3.0}, "Selected_lepton_charge","evWeight","");
    add1DHist( {"Deltaeta_Leptons","Delta eta between 2 leptons; Delta eta between the two leptons;Events;", 20, 0.0, 4.0}, "Selected_lepton_deltaeta","evWeight","0");
    add1DHist( {"Deltaphi_Leptons","Delta phi between 2 leptons; Delta phi between the two leptons;Events;", 20, 0.0, 4.0}, "Selected_lepton_deltaphi","evWeight","0");
    add1DHist( {"DeltaR_Leptons","Delta R between 2 leptons; Delta R between the two leptons;Events;", 30, 0.0, 6.0}, "Selected_lepton_deltaR","evWeight","0");
    add1DHist( {"miniPFRelIso_all_Leptons","miniPFRelIso_all isolation of leptons; miniPFRel_all isolation of the leptons;Events;", 40, 0.0, 0.4}, "Selected_lepton_miniPFRelIso_all","evWeight","");
    add1DHist( {"Sip_Leptons","Impact parameter significance of leptons; Impact parameter significance of the leptons wrt PV;Events;", 15, 0.0, 3.0}, "Selected_lepton_sip3d","evWeight","");
    add1DHist( {"Pt_Sum_Two_Leptons","Pt of sum of two leptons; Pt of the sum of the two leptons (GeV);Events;", 30, 0.0, 600.0}, "Vectorial_sum_two_leptons_pt","evWeight","0");
    add1DHist( {"Mass_Sum_Two_Leptons","Invariant mass of two leptons; Invariant mass of the two leptons (GeV);Events;", 20, 0.0, 800.0}, "Vectorial_sum_two_leptons_mass","evWeight","0");

    ////=== jets =====//
    add1DHist( {"Number_Jets","Number of jets; Number of jets;Events;", 10, 0.0, 10.0}, "Selected_clean_jet_number","evWeight","0");
    add1DHist( {"Pt_Jets","Pt of jets; Pt of the jets (GeV);Events;", 25, 0.0, 500.0}, "Selected_clean_jet_pt","evWeight","");
    add1DHist( {"Leading_Pt_Jets","Pt of leading jet; Pt of the leading jet (GeV);Events;", 25, 0.0, 500.0}, "Selected_clean_jet_leading_pt","evWeight","0");
    add1DHist( {"Subleading_Pt_Jets","Pt of subleading jet; Pt of the subleading jet (GeV);Events;", 20, 0.0, 400.0}, "Selected_clean_jet_subleading_pt","evWeight","0");
    add1DHist( {"Subsubleading_Pt_Jets","Pt of subsubleading jet; Pt of the subsubleading jet (GeV);Events;", 15, 0.0, 300.0}, "Selected_clean_jet_subsubleading_pt","evWeight","0");
    add1DHist( {"Ht_Jets","Ht; Ht (GeV);Events;", 30, 0.0, 1200.0}, "Selected_clean_jet_Ht","evWeight","");
    add1DHist( {"Eta_Jets","Eta of jets; Eta of the jets;Events;", 50, -5.0, 5.0}, "Selected_clean_jet_eta","evWeight","");
    add1DHist( {"Phi_Jets","Phi of jets; Phi of the jets;Events;", 35, -3.5, 3.5}, "Selected_clean_jet_phi","evWeight","");
    add1DHist( {"Mass_Jets","Mass of jets; Mass of the jets (GeV);Events;", 10, 0.0, 50.0}, "Selected_clean_jet_mass","evWeight","");
    add1DHist( {"DeepB_Jets","DeepCSV of jets; DeepCSV of the jets;Events;", 20, 0.0, 1.0}, "Selected_clean_jet_btagDeepB","evWeight","");
    add1DHist( {"DeepFlavB_Jets","DeepJet of jets; DeepJet of the jets;Events;", 20, 0.0, 1.0}, "Selected_clean_jet_btagDeepFlavB","evWeight","");
    add1DHist( {"Number_bJets","Number of bjets; Number of bjets;Events;", 10, 0.0, 10.0}, "Selected_clean_bjet_number","evWeight","0");
    add1DHist( {"Pt_bJets","Pt of bjets; Pt of the bjets (GeV);Events;", 25, 0.0, 500.0}, "Selected_clean_bjet_pt","evWeight","");
    add1DHist( {"Leading_Pt_bJets","Pt of leading bjet; Pt of the leading bjet (GeV);Events;", 25, 0.0, 500.0}, "Selected_clean_bjet_leading_pt","evWeight","0");
    add1DHist( {"Subleading_Pt_bJets","Pt of subleading bjet; Pt of the subleading bjet (GeV);Events;", 20, 0.0, 400.0}, "Selected_clean_bjet_subleading_pt","evWeight","0");
    add1DHist( {"Subsubleading_Pt_bJets","Pt of subsubleading bjet; Pt of the subsubleading bjet (GeV);Events;", 15, 0.0, 300.0}, "Selected_clean_bjet_subsubleading_pt","evWeight","0");
    add1DHist( {"Eta_bJets","Eta of bjets; Eta of the bjets;Events;", 30, -3.0, 3.0}, "Selected_clean_bjet_eta","evWeight","");
    add1DHist( {"Phi_bJets","Phi of bjets; Phi of the bjets;Events;", 35, -3.5, 3.5}, "Selected_clean_bjet_phi","evWeight","");
    add1DHist( {"Mass_bJets","Mass of bjets; Mass of the bjets (GeV);Events;", 10, 0.0, 50.0}, "Selected_clean_bjet_mass","evWeight","");
    add1DHist( {"DeepB_bJets","DeepCSV of bjets; DeepCSV of the bjets;Events;", 20, 0.0, 1.0}, "Selected_clean_bjet_btagDeepB","evWeight","");
    add1DHist( {"DeepFlavB_bJets","DeepJet of bjets; DeepJet of the bjets;Events;", 20, 0.0, 1.0}, "Selected_clean_bjet_btagDeepFlavB","evWeight","");
    add1DHist( {"Mass_Sum_Three_Jets","Invariant mass of three jets; Invariant mass(3j) (GeV);Events;", 25, 0.0, 3000.0}, "Vectorial_sum_three_clean_jets_mass","evWeight","0");
    add1DHist( {"Mass_Sum_Three_Jets_offset","Invariant mass of three jets (offset); Invariant mass (3j) - top mass (GeV);Events;", 25, -200.0, 2800.0}, "Vectorial_sum_three_clean_jets_mass_offset","evWeight","0");
    add1DHist( {"Mass_Sum_Three_Jets_min","Invariant mass of three jets (min); Minimum invariant mass (3j) (GeV);Events;", 40, 0.0, 400.0}, "Vectorial_sum_three_clean_jets_mass_min","evWeight","0");

    ////=== others ====//
    add1DHist( {"St","St; St (GeV);Events;", 25, 0.0, 2000.0}, "St","evWeight","0");
    add1DHist( {"St_with_MET","St with MET; St with MET (GeV);Events;", 25, 0.0, 2000.0}, "St_with_MET","evWeight","0");
    // add1DHist( {"Mass_Sum_All_Particles","Invariant mass of leptons + jets; Invariant mass(leptons+jets) (GeV);Events;", 25, 0.0, 3000.0}, "Vectorial_sum_particles_mass","evWeight","0");
    // add1DHist( {"DeltaR_Leptons_Jets","DeltaR between leptons and jets; DeltaR between the leptons and the jets;Events;", 60, 0.0, 6.0}, "deltaR_Selected_lepton_Selected_clean_jet","evWeight","0");
    // add1DHist( {"DeltaR_Leptons_bJets","DeltaR between leptons and bjets; DeltaR between the leptons and the bjets;Events;", 60, 0.0, 6.0}, "deltaR_Selected_lepton_Selected_clean_bjet","evWeight","0");
    // add1DHist( {"DeltaR_Muons_Jets","DeltaR between muons and jets; DeltaR between the muons and the jets;Events;", 60, 0.0, 6.0}, "deltaR_Selected_muon_Selected_clean_jet","evWeight","0");
    // add1DHist( {"DeltaR_Muons_bJets","DeltaR between muons and bjets; DeltaR between the muons and the bjets;Events;", 60, 0.0, 6.0}, "deltaR_Selected_muon_Selected_clean_bjet","evWeight","0");
    // add1DHist( {"Ratio_Pt_Leptons_Pt_Jets","Ratio between Pt of leptons and Pt of jets; Ratio between the Pt of leptons and the Pt of jets;Events;", 50, 0.0, 10.0}, "Pt_lepton_over_Pt_clean_jet","evWeight","0");
    // add1DHist( {"Ratio_Pt_Leptons_Pt_bJets","Ratio between Pt of leptons and Pt of bjets; Ratio between the Pt of leptons and the Pt of bjets;Events;", 50, 0.0, 10.0}, "Pt_lepton_over_Pt_clean_bjet","evWeight","0");
    // add1DHist( {"Ratio_Pt_Muons_Pt_Jets","Ratio between Pt of muons and Pt of jets; Ratio between the Pt of muons and the Pt of jets;Events;", 50, 0.0, 10.0}, "Pt_muon_over_Pt_clean_jet","evWeight","0");
    // add1DHist( {"Ratio_Pt_Muons_Pt_bJets","Ratio between Pt of muons and Pt of bjets; Ratio between the Pt of muons and the Pt of bjets;Events;", 50, 0.0, 10.0}, "Pt_muon_over_Pt_clean_bjet","evWeight","0");
    // add1DHist( {"Min_DeltaR_Leptons_Jets","Minimum of DeltaR between leptons and jets; Minimum of DeltaR between the leptons and the jets;Events;", 40, 0.0, 4.0}, "Min_deltaR_Selected_lepton_Selected_clean_jet","evWeight","0");
    // add1DHist( {"Min_DeltaR_Leptons_bJets","Minimum of DeltaR between leptons and bjets; Minimum of DeltaR between the leptons and the bjets;Events;", 40, 0.0, 4.0}, "Min_deltaR_Selected_lepton_Selected_clean_bjet","evWeight","0");
    // add1DHist( {"Min_DeltaR_Muons_Jets","Minimum of DeltaR between muons and jets; Minimum of DeltaR between the muons and the jets;Events;", 40, 0.0, 4.0}, "Min_deltaR_Selected_muon_Selected_clean_jet","evWeight","0");
    // add1DHist( {"Min_DeltaR_Muons_bJets","Minimum of DeltaR between muons and bjets; Minimum of DeltaR between the muons and the bjets;Events;", 40, 0.0, 4.0}, "Min_deltaR_Selected_muon_Selected_clean_bjet","evWeight","0");
    // add1DHist( {"Min_Ratio_Pt_Leptons_Pt_Jets","Minimum of Ratio between Pt of leptons and Pt of jets; Minimum of Ratio between the Pt of leptons and the Pt of jets;Events;", 50, 0.0, 5.0}, "Min_Pt_lepton_over_Pt_clean_jet","evWeight","0");
    // add1DHist( {"Min_Ratio_Pt_Leptons_Pt_bJets","Minimum of Ratio between Pt of leptons and Pt of bjets; Minimum of Ratio between the Pt of leptons and the Pt of bjets;Events;", 50, 0.0, 5.0}, "Min_Pt_lepton_over_Pt_clean_bjet","evWeight","0");
    // add1DHist( {"Min_Ratio_Pt_Muons_Pt_Jets","Minimum of Ratio between Pt of muons and Pt of jets; Minimum of Ratio between the Pt of muons and the Pt of jets;Events;", 50, 0.0, 5.0}, "Min_Pt_muon_over_Pt_clean_jet","evWeight","0");
    // add1DHist( {"Min_Ratio_Pt_Muons_Pt_bJets","Minimum of Ratio between Pt of muons and Pt of bjets; Minimum of Ratio between the Pt of muons and the Pt of bjets;Events;", 50, 0.0, 5.0}, "Min_Pt_muon_over_Pt_clean_bjet","evWeight","0");
    if(!_isData)
    {
        add1DHist( {"Tprime_mass","Mass of Tprime; Mass of the Tprime (GeV);Events;", 20, 0.0, 1600.0}, "L2bJ1Fromtop_L1J1J2FromHiggs_mass","evWeight","0");
        add1DHist( {"Tprime_transverse_mass","Transverse mass of Tprime; Transverse mass of the Tprime (GeV);Events;", 20, 0.0, 1600.0}, "L2bJ1Fromtop_L1J1J2FromHiggs_transverse_mass","evWeight","0");
    }

    add2DHist( {"Tprime_transverse_mass_Sum_Pt_Two_Leptons","Transverse mass of Tprime * Sum of Pt of two leptons; Transverse mass of Tprime (GeV); Sum of the Pt of the two leptons (GeV);", 20, 0.0, 1600.0, 50, 0.0, 600.0}, "L2bJ1Fromtop_L1J1J2FromHiggs_transverse_mass","Selected_lepton_sum_two_leptons_pt","evWeight","0");
    add2DHist( {"Tprime_transverse_mass_Mass_Sum_Three_Jets_min","Transverse mass of Tprime * Invariant mass of three jets (min); Transverse mass of Tprime (GeV); Minimum invariant mass (3j) (GeV);", 20, 0.0, 1600.0, 40, 0.0, 400.0}, "L2bJ1Fromtop_L1J1J2FromHiggs_transverse_mass","Vectorial_sum_three_clean_jets_mass_min","evWeight","0");

}

void TprimeAnalyser::bookReconstructedHists()
{
    // ////=== reconstructed leptons ====//
    // add1DHist( {"Pt_LeptonFromtop","Pt of lepton from top; Pt of the lepton coming from top (GeV);Events;", 20, 0.0, 400.0}, "LeptonFromtop_pt","evWeight","0");
    // add1DHist( {"Eta_LeptonFromtop","Eta of lepton from top; Eta of the lepton coming from top;Events;", 30, -3.0, 3.0}, "LeptonFromtop_eta","evWeight","0");
    // add1DHist( {"Phi_LeptonFromtop","Phi of lepton from top; Phi of the lepton coming from top;Events;", 35, -3.5, 3.5}, "LeptonFromtop_phi","evWeight","0");
    // add1DHist( {"Pt_LeptonFromHiggs","Pt of lepton from Higgs; Pt of the lepton coming from Higgs (GeV);Events;", 20, 0.0, 400.0}, "LeptonFromHiggs_pt","evWeight","0");
    // add1DHist( {"Eta_LeptonFromHiggs","Eta of lepton from Higgs; Eta of the lepton coming from Higgs;Events;", 100, -3.0, 3.0}, "LeptonFromHiggs_eta","evWeight","0");
    // add1DHist( {"Phi_LeptonFromHiggs","Phi of lepton from Higgs; Phi of the lepton coming from Higgs;Events;", 35, -3.5, 3.5}, "LeptonFromHiggs_phi","evWeight","0");
    // add1DHist( {"Lepton_pdgId","Nature of lepton; Nature of the leptons;Events;", 50, -25, 25}, "LeptonFromW_pdgId","evWeight","");
    add1DHist( {"MotherFromLepton1_pdgId","Nature of mother of first lepton; Nature of the mother of the first lepton;Events;", 50, 0, 50}, "MotherfromLepton1_pdgId","evWeight","0");
    add1DHist( {"MotherFromLepton2_pdgId","Nature of mother of second lepton; Nature of the mother of the second lepton;Events;", 50, 0, 50}, "MotherfromLepton2_pdgId","evWeight","0");
    add1DHist( {"MotherFromLepton12_pdgId","Sum of nature of mothers of leptons; Sum of the nature of the mothers of the two leptons;Events;", 50, 0, 50}, "MotherfromLepton12_pdgId","evWeight","0");
    // add1DHist( {"Pt_NeutrinoFromtop","Pt of neutrino from top; Pt of the neutrino coming from top (GeV);Events;", 20, 0.0, 400.0}, "NeutrinoFromtop_pt","evWeight","0");
    // add1DHist( {"Phi_NeutrinoFromtop","Phi of neutrino from top; Phi of the neutrino coming from top;Events;", 35, -3.5, 3.5}, "NeutrinoFromtop_phi","evWeight","0");
    // add1DHist( {"Pt_NeutrinoFromHiggs","Pt of neutrino from Higgs; Pt of the neutrino coming from Higgs (GeV);Events;", 20, 0.0, 400.0}, "NeutrinoFromHiggs_pt","evWeight","0");
    // add1DHist( {"Phi_NeutrinoFromHiggs","Phi of neutrino from Higgs; Phi of the neutrino coming from Higgs;Events;", 35, -3.5, 3.5}, "NeutrinoFromHiggs_phi","evWeight","0");

    // ////=== reconstructed jets ====//
    // add1DHist( {"Pt_bJetFromtop","Pt of bjet from top; Pt of the bjet coming from top (GeV);Events;", 80, 0.0, 400.0}, "bJetFromtop1_pt","one","0");
    // add1DHist( {"Eta_bJetFromtop","Eta of bjet from top; Eta of the bjet coming from top;Events;", 100, -3.0, 3.0}, "bJetFromtop1_eta","one","0");
    // add1DHist( {"Phi_bJetFromtop","Phi of bjet from top; Phi of the bjet coming from top;Events;", 100, -3.5, 3.5}, "bJetFromtop1_phi","one","0");
    // add1DHist( {"Mass_bJetFromtop","Mass of bjet from top; Mass of the bjet coming from top (GeV);Events;", 50, 0.0, 50.0}, "bJetFromtop1_mass","one","0"); 
    // add1DHist( {"Pt_JetFromHiggs1","Pt of first jet from Higgs; Pt of the first jet coming from Higgs (GeV);Events;", 80, 0.0, 400.0}, "JetFromHiggs1_pt","one","0");
    // add1DHist( {"Eta_JetFromHiggs1","Eta of first jet from Higgs; Eta of the first jet coming from Higgs;Events;", 100, -3.0, 3.0}, "JetFromHiggs1_eta","one","0");
    // add1DHist( {"Phi_JetFromHiggs1","Phi of fist jet from Higgs; Phi of the first jet coming from Higgs;Events;", 100, -3.5, 3.5}, "JetFromHiggs1_phi","one","0");
    // add1DHist( {"Mass_JetFromHiggs1","Mass of first jet from Higgs; Mass of the first jet coming from Higgs (GeV);Events;", 10, 0.0, 50.0}, "JetFromHiggs1_mass","one","0");
    // add1DHist( {"Pt_JetFromHiggs2","Pt of second jet from Higgs; Pt of the second jet coming from Higgs (GeV);Events;", 40, 0.0, 200.0}, "JetFromHiggs2_pt","one","0");
    // add1DHist( {"Eta_JetFromHiggs2","Eta of second jet from Higgs; Eta of the second jet coming from Higgs;Events;", 100, -3.0, 3.0}, "JetFromHiggs2_eta","one","0");
    // add1DHist( {"Phi_JetFromHiggs2","Phi of second jet from Higgs; Phi of the second jet coming from Higgs;Events;", 100, -3.5, 3.5}, "JetFromHiggs2_phi","one","0");
    // add1DHist( {"Mass_JetFromHiggs2","Mass of second jet from Higgs; Mass of the second jet coming from Higgs (GeV);Events;", 10, 0.0, 30.0}, "JetFromHiggs2_mass","one","0");

    ////=== combination of reconstructed particles ====//
    // add1DHist( {"DeltaR_LeptonFromtop_Jets","DeltaR between lepton coming from top and jets; DeltaR between the lepton coming from top and the jets;Events;", 60, 0.0, 6.0}, "deltaR_Selected_leptonfromtop_Selected_clean_jet","one","0");
    // add1DHist( {"DeltaR_LeptonFromtop_bJets","DeltaR between lepton coming from top and bjets; DeltaR between the lepton coming from top and the bjets;Events;", 60, 0.0, 6.0}, "deltaR_Selected_leptonfromtop_Selected_clean_bjet","one","0");
    // add1DHist( {"DeltaR_LeptonFromHiggs_Jets","DeltaR between lepton coming from Higgs and jets; DeltaR between the lepton coming from Higgs and the jets;Events;", 60, 0.0, 6.0}, "deltaR_Selected_leptonfromHiggs_Selected_clean_jet","one","0");
    // add1DHist( {"DeltaR_LeptonFromHiggs_bJets","DeltaR between lepton coming from Higgs and bjets; DeltaR between the lepton coming from Higgs and the bjets;Events;", 60, 0.0, 6.0}, "deltaR_Selected_leptonfromHiggs_Selected_clean_bjet","one","0");
    // add1DHist( {"Ratio_Pt_LeptonFromtop_Pt_Jets","Ratio between Pt of lepton coming from top and Pt of jets; Ratio between the Pt of lepton coming from top and the Pt of jets;Events;", 50, 0.0, 10.0}, "Pt_leptonfromtop_over_Pt_clean_jet","one","0");
    // add1DHist( {"Ratio_Pt_LeptonFromtop_Pt_bJets","Ratio between Pt of lepton coming from top and Pt of bjets; Ratio between the Pt of lepton coming from top and the Pt of bjets;Events;", 50, 0.0, 10.0}, "Pt_leptonfromtop_over_Pt_clean_bjet","one","0");
    // add1DHist( {"Ratio_Pt_LeptonFromHiggs_Pt_Jets","Ratio between Pt of lepton coming from Higgs and Pt of jets; Ratio between the Pt of lepton coming from Higgs and the Pt of jets;Events;", 50, 0.0, 10.0}, "Pt_leptonfromHiggs_over_Pt_clean_jet","one","0");
    // add1DHist( {"Ratio_Pt_LeptonFromHiggs_Pt_bJets","Ratio between Pt of lepton coming from Higgs and Pt of bjets; Ratio between the Pt of lepton coming from Higgs and the Pt of bjets;Events;", 50, 0.0, 10.0}, "Pt_leptonfromHiggs_over_Pt_clean_bjet","one","0");
    // add1DHist( {"Min_DeltaR_LeptonFromtop_Jets","Minimum of DeltaR between lepton coming from top and jets; Minimum of DeltaR between the lepton coming from top and the jets;Events;", 40, 0.0, 4.0}, "Min_deltaR_Selected_leptonfromtop_Selected_clean_jet","one","0");
    // add1DHist( {"Min_DeltaR_LeptonFromtop_bJets","Minimum of DeltaR between lepton coming from top and bjets; Minimum of DeltaR between the lepton coming from top and the bjets;Events;", 40, 0.0, 4.0}, "Min_deltaR_Selected_leptonfromtop_Selected_clean_bjet","one","0");
    // add1DHist( {"Min_DeltaR_LeptonFromHiggs_Jets","Minimum of DeltaR between lepton coming from Higgs and jets; Minimum of DeltaR between the lepton coming from Higgs and the jets;Events;", 40, 0.0, 4.0}, "Min_deltaR_Selected_leptonfromHiggs_Selected_clean_jet","one","0");
    // add1DHist( {"Min_DeltaR_LeptonFromHiggs_bJets","Minimum of DeltaR between lepton coming from Higgs and bjets; Minimum of DeltaR between the lepton coming from Higgs and the bjets;Events;", 40, 0.0, 4.0}, "Min_deltaR_Selected_leptonfromHiggs_Selected_clean_bjet","one","0");
    // add1DHist( {"Min_Ratio_Pt_LeptonFromtop_Pt_Jets","Minimum of Ratio between Pt of lepton coming from top and Pt of jets; Minimum of Ratio between the Pt of lepton coming from top and the Pt of jets;Events;", 50, 0.0, 5.0}, "Min_Pt_leptonfromtop_over_Pt_clean_jet","one","0");
    // add1DHist( {"Min_Ratio_Pt_LeptonFromtop_Pt_bJets","Minimum of Ratio between Pt of lepton coming from top and Pt of bjets; Minimum of Ratio between the Pt of lepton coming from top and the Pt of bjets;Events;", 50, 0.0, 5.0}, "Min_Pt_leptonfromtop_over_Pt_clean_bjet","one","0");
    // add1DHist( {"Min_Ratio_Pt_LeptonFromHiggs_Pt_Jets","Minimum of Ratio between Pt of lepton coming from Higgs and Pt of jets; Minimum of Ratio between the Pt of lepton coming from Higgs and the Pt of jets;Events;", 50, 0.0, 5.0}, "Min_Pt_leptonfromHiggs_over_Pt_clean_jet","one","0");
    // add1DHist( {"Min_Ratio_Pt_LeptonFromHiggs_Pt_bJets","Minimum of Ratio between Pt of lepton coming from Higgs and Pt of bjets; Minimum of Ratio between the Pt of lepton coming from Higgs and the Pt of bjets;Events;", 50, 0.0, 5.0}, "Min_Pt_leptonfromHiggs_over_Pt_clean_bjet","one","0");

    // add1DHist( {"Pt_LeptonFromtop_bJetFromtop","Pt of lepton from top and bjet from top; Pt(lepton_top,bjet_top) (GeV);Events;", 100, 0.0, 500.0}, "LeptonFromtop_bJetFromtop1_pt","one","0");
    // add1DHist( {"Eta_LeptonFromtop_bJetFromtop","Eta of lepton from top and bjet from top; Eta(lepton_top,bjet_top);Events;", 100, -3.0, 3.0}, "LeptonFromtop_bJetFromtop1_eta","one","0");
    // add1DHist( {"Phi_LeptonFromtop_bJetFromtop","Phi of lepton from top and bjet from top; Phi(lepton_top,bjet_top);Events;", 100, -3.5, 3.5}, "LeptonFromtop_bJetFromtop1_phi","one","0");
    // add1DHist( {"Mass_LeptonFromtop_bJetFromtop","Invariant mass of lepton from top and bjet from top; Invariant mass(lepton_top,bjet_top) (GeV);Events;", 50, 0.0, 200.0}, "LeptonFromtop_bJetFromtop1_mass","one","0");
    // add1DHist( {"Deltaeta_LeptonFromtop_bJetFromtop","Delta eta of lepton from top and bjet from top; Delta_eta(lepton_top,bjet_top);Events;", 30, 0.0, 4.0}, "LeptonFromtop_bJetFromtop1_deltaeta","one","0");    
    // add1DHist( {"Deltaphi_LeptonFromtop_bJetFromtop","Delta phi of lepton from top and bjet from top; Delta_phi(lepton_top,bjet_top);Events;", 30, 0.0, 4.0}, "LeptonFromtop_bJetFromtop1_deltaphi","one","0");    
    // add1DHist( {"DeltaR_LeptonFromtop_bJetFromtop","Delta R of lepton from top and bjet from top; Delta_R(lepton_top,bjet_top);Events;", 30, 0.0, 6.0}, "LeptonFromtop_bJetFromtop1_deltaR","one","0");
    // add1DHist( {"Pt_LeptonFromtop_JetFromHiggs1","Pt of lepton from top and first jet from Higgs; Pt(lepton_top,jet1_Higgs) (GeV);Events;", 100, 0.0, 500.0}, "LeptonFromtop_JetFromHiggs1_pt","one","0");
    // add1DHist( {"Eta_LeptonFromtop_JetFromHiggs1","Eta of lepton from top and first jet from Higgs; Eta(lepton_top,jet1_Higgs);Events;", 100, -3.0, 3.0}, "LeptonFromtop_JetFromHiggs1_eta","one","0");
    // add1DHist( {"Phi_LeptonFromtop_JetFromHiggs1","Phi of lepton from top and first jet from Higgs; Phi(lepton_top,jet1_Higgs);Events;", 100, -3.5, 3.5}, "LeptonFromtop_JetFromHiggs1_phi","one","0");
    // add1DHist( {"Mass_LeptonFromtop_JetFromHiggs1","Invariant mass of lepton from top and first jet from Higgs; Invariant mass(lepton_top,jet1_Higgs) (GeV);Events;", 125, 0.0, 500.0}, "LeptonFromtop_JetFromHiggs1_mass","one","0");
    // add1DHist( {"Deltaeta_LeptonFromtop_JetFromHiggs1","Delta eta of lepton from top and first jet from Higgs; Delta_eta(lepton_top,jet1_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromtop_JetFromHiggs1_deltaeta","one","0");    
    // add1DHist( {"Deltaphi_LeptonFromtop_JetFromHiggs1","Delta phi of lepton from top and first jet from Higgs; Delta_phi(lepton_top,jet1_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromtop_JetFromHiggs1_deltaphi","one","0");    
    // add1DHist( {"DeltaR_LeptonFromtop_JetFromHiggs1","Delta R of lepton from top and first jet from Higgs; Delta_R(lepton_top,jet1_Higgs);Events;", 30, 0.0, 6.0}, "LeptonFromtop_JetFromHiggs1_deltaR","one","0");
    // add1DHist( {"Pt_LeptonFromtop_JetFromHiggs2","Pt of lepton from top and second jet from Higgs; Pt(lepton_top,jet2_Higgs) (GeV);Events;", 40, 0.0, 200.0}, "LeptonFromtop_JetFromHiggs2_pt","one","0");
    // add1DHist( {"Eta_LeptonFromtop_JetFromHiggs2","Eta of lepton from top and second jet from Higgs; Eta(lepton_top,jet2_Higgs);Events;", 100, -3.0, 3.0}, "LeptonFromtop_JetFromHiggs2_eta","one","0");
    // add1DHist( {"Phi_LeptonFromtop_JetFromHiggs2","Phi of lepton from top and second jet from Higgs; Phi(lepton_top,jet2_Higgs);Events;", 100, -3.5, 3.5}, "LeptonFromtop_JetFromHiggs2_phi","one","0");
    // add1DHist( {"Mass_LeptonFromtop_JetFromHiggs2","Invariant mass of lepton from top and second jet from Higgs; Invariant mass(lepton_top,jet2_Higgs) (GeV);Events;", 125, 0.0, 500.0}, "LeptonFromtop_JetFromHiggs2_mass","one","0");
    // add1DHist( {"Deltaeta_LeptonFromtop_JetFromHiggs2","Delta eta of lepton from top and second jet from Higgs; Delta_eta(lepton_top,jet2_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromtop_JetFromHiggs2_deltaeta","one","0");    
    // add1DHist( {"Deltaphi_LeptonFromtop_JetFromHiggs2","Delta phi of lepton from top and second jet from Higgs; Delta_phi(lepton_top,jet2_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromtop_JetFromHiggs2_deltaphi","one","0");    
    // add1DHist( {"DeltaR_LeptonFromtop_JetFromHiggs2","Delta R of lepton from top and second jet from Higgs; Delta_R(lepton_top,jet2_Higgs);Events;", 30, 0.0, 6.0}, "LeptonFromtop_JetFromHiggs2_deltaR","one","0");
    // add1DHist( {"Pt_LeptonFromHiggs_bJetFromtop","Pt of lepton from Higgs and bjet from top; Pt(lepton_Higgs,bjet_top) (GeV);Events;", 80, 0.0, 400.0}, "LeptonFromHiggs_bJetFromtop1_pt","one","0");
    // add1DHist( {"Eta_LeptonFromHiggs_bJetFromtop","Eta of lepton from Higgs and bjet from top; Eta(lepton_Higgs,bjet_top);Events;", 100, -3.0, 3.0}, "LeptonFromHiggs_bJetFromtop1_eta","one","0");
    // add1DHist( {"Phi_LeptonFromHiggs_bJetFromtop","Phi of lepton from Higgs and bjet from top; Phi(lepton_Higgs,bjet_top);Events;", 100, -3.5, 3.5}, "LeptonFromHiggs_bJetFromtop1_phi","one","0");
    // add1DHist( {"Mass_LeptonFromHiggs_bJetFromtop","Invariant mass of lepton from Higgs and bjet from top; Invariant mass(lepton_Higgs,bjet_top) (GeV);Events;", 125, 0.0, 500.0}, "LeptonFromHiggs_bJetFromtop1_mass","one","0");
    // add1DHist( {"Deltaeta_LeptonFromHiggs_bJetFromtop","Delta eta of lepton from Higgs and bjet from top; Delta_eta(lepton_Higgs,bjet_top);Events;", 30, 0.0, 4.0}, "LeptonFromHiggs_bJetFromtop1_deltaeta","one","0");    
    // add1DHist( {"Deltaphi_LeptonFromHiggs_bJetFromtop","Delta phi of lepton from Higgs and bjet from top; Delta_phi(lepton_Higgs,bjet_top);Events;", 30, 0.0, 4.0}, "LeptonFromHiggs_bJetFromtop1_deltaphi","one","0");    
    // add1DHist( {"DeltaR_LeptonFromHiggs_bJetFromtop","Delta R of lepton from Higgs and bjet from top; Delta_R(lepton_Higgs,bjet_top);Events;", 30, 0.0, 6.0}, "LeptonFromHiggs_bJetFromtop1_deltaR","one","0");
    // add1DHist( {"Pt_LeptonFromHiggs_JetFromHiggs1","Pt of lepton from Higgs and first jet from Higgs; Pt(lepton_Higgs,jet1_Higgs) (GeV);Events;", 100, 0.0, 500.0}, "LeptonFromHiggs_JetFromHiggs1_pt","one","0");
    // add1DHist( {"Eta_LeptonFromHiggs_JetFromHiggs1","Eta of lepton from Higgs and first jet from Higgs; Eta(lepton_Higgs,jet1_Higgs);Events;", 100, -3.0, 3.0}, "LeptonFromHiggs_JetFromHiggs1_eta","one","0");
    // add1DHist( {"Phi_LeptonFromHiggs_JetFromHiggs1","Phi of lepton from Higgs and first jet from Higgs; Phi(lepton_Higgs,jet1_Higgs);Events;", 100, -3.5, 3.5}, "LeptonFromHiggs_JetFromHiggs1_phi","one","0");
    // add1DHist( {"Mass_LeptonFromHiggs_JetFromHiggs1","Invariant mass of lepton from Higgs and first jet from Higgs; Invariant mass(lepton_Higgs,jet1_Higgs) (GeV);Events;", 40, 0.0, 150.0}, "LeptonFromHiggs_JetFromHiggs1_mass","one","0");
    // add1DHist( {"Deltaeta_LeptonFromHiggs_JetFromHiggs1","Delta eta of lepton from Higgs and first jet from Higgs; Delta_eta(lepton_Higgs,jet1_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromHiggs_JetFromHiggs1_deltaeta","one","0");    
    // add1DHist( {"Deltaphi_LeptonFromHiggs_JetFromHiggs1","Delta phi of lepton from Higgs and first jet from Higgs; Delta_phi(lepton_Higgs,jet1_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromHiggs_JetFromHiggs1_deltaphi","one","0");    
    // add1DHist( {"DeltaR_LeptonFromHiggs_JetFromHiggs1","Delta R of lepton from Higgs and first jet from Higgs; Delta_R(lepton_Higgs,jet1_Higgs);Events;", 30, 0.0, 6.0}, "LeptonFromHiggs_JetFromHiggs1_deltaR","one","0");
    // add1DHist( {"Pt_LeptonFromHiggs_JetFromHiggs2","Pt of lepton from Higgs and second jet from Higgs; Pt(lepton_Higgs,jet2_Higgs) (GeV);Events;", 60, 0.0, 300.0}, "LeptonFromHiggs_JetFromHiggs2_pt","one","0");
    // add1DHist( {"Eta_LeptonFromHiggs_JetFromHiggs2","Eta of lepton from Higgs and second jet from Higgs; Eta(lepton_Higgs,jet2_Higgs);Events;", 100, -3.0, 3.0}, "LeptonFromHiggs_JetFromHiggs2_eta","one","0");
    // add1DHist( {"Phi_LeptonFromHiggs_JetFromHiggs2","Phi of lepton from Higgs and second jet from Higgs; Phi(lepton_Higgs,jet2_Higgs);Events;", 100, -3.5, 3.5}, "LeptonFromHiggs_JetFromHiggs2_phi","one","0");
    // add1DHist( {"Mass_LeptonFromHiggs_JetFromHiggs2","Invariant mass of lepton from Higgs and second jet from Higgs; Invariant mass(lepton_Higgs,jet2_Higgs) (GeV);Events;", 25, 0.0, 100.0}, "LeptonFromHiggs_JetFromHiggs2_mass","one","0");
    // add1DHist( {"Deltaeta_LeptonFromHiggs_JetFromHiggs2","Delta eta of lepton from Higgs and second jet from Higgs; Delta_eta(lepton_Higgs,jet2_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromHiggs_JetFromHiggs2_deltaeta","one","0");    
    // add1DHist( {"Deltaphi_LeptonFromHiggs_JetFromHiggs2","Delta phi of lepton from Higgs and second jet from Higgs; Delta_phi(lepton_Higgs,jet2_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromHiggs_JetFromHiggs2_deltaphi","one","0");    
    // add1DHist( {"DeltaR_LeptonFromHiggs_JetFromHiggs2","Delta R of lepton from Higgs and second jet from Higgs; Delta_R(lepton_Higgs,jet2_Higgs);Events;", 30, 0.0, 6.0}, "LeptonFromHiggs_JetFromHiggs2_deltaR","one","0");
    // add1DHist( {"Pt_JetFromHiggs1_JetFromHiggs2","Pt of first and second jets from Higgs; Pt(jet1_Higgs,jet2_Higgs) (GeV);Events;", 100, 0.0, 500.0}, "JetFromHiggs1_JetFromHiggs2_pt","one","0");
    // add1DHist( {"Eta_JetFromHiggs1_JetFromHiggs2","Eta of first and second jets from Higgs; Eta(jet1_Higgs,jet2_Higgs);Events;", 100, -3.0, 3.0}, "JetFromHiggs1_JetFromHiggs2_eta","one","0");
    // add1DHist( {"Phi_JetFromHiggs1_JetFromHiggs2","Phi of first and second jets from Higgs; Phi(jet1_Higgs,jet2_Higgs);Events;", 100, -3.5, 3.5}, "JetFromHiggs1_JetFromHiggs2_phi","one","0");
    // add1DHist( {"Mass_JetFromHiggs1_JetFromHiggs2","Invariant mass of first and second jets from Higgs; Invariant mass(jet1_Higgs,jet2_Higgs) (GeV);Events;", 50, 0.0, 200.0}, "JetFromHiggs1_JetFromHiggs2_mass","one","0");
    // add1DHist( {"Deltaeta_JetFromHiggs1_JetFromHiggs2","Delta eta of first and second jets from Higgs; Delta_eta(jet1_Higgs,jet2_Higgs);Events;", 30, 0.0, 4.0}, "JetFromHiggs1_JetFromHiggs2_deltaeta","one","0");    
    // add1DHist( {"Deltaphi_JetFromHiggs1_JetFromHiggs2","Delta phi of first and second jets from Higgs; Delta_phi(jet1_Higgs,jet2_Higgs);Events;", 30, 0.0, 4.0}, "JetFromHiggs1_JetFromHiggs2_deltaphi","one","0");    
    // add1DHist( {"DeltaR_JetFromHiggs1_JetFromHiggs2","Delta R of first and second jets from Higgs; Delta_R(jet1_Higgs,jet2_Higgs);Events;", 30, 0.0, 6.0}, "JetFromHiggs1_JetFromHiggs2_deltaR","one","0");
    // add1DHist( {"Pt_LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2","Pt of lepton and first and second jets from Higgs; Pt(lepton_Higgs,jet12_Higgs) (GeV);Events;", 120, 0.0, 600.0}, "LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_pt","one","0");
    // add1DHist( {"Eta_LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2","Eta of lepton and first and second jets from Higgs; Eta(lepton_Higgs,jet12_Higgs);Events;", 100, -3.0, 3.0}, "LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_eta","one","0");
    // add1DHist( {"Phi_LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2","Phi of lepton and first and second jets from Higgs; Phi(lepton_Higgs,jet12_Higgs);Events;", 100, -3.5, 3.5}, "LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_phi","one","0");
    // add1DHist( {"Mass_LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2","Invariant mass of lepton and first and second jets from Higgs; Invariant mass(lepton_Higgs,jet12_Higgs) (GeV);Events;", 50, 0.0, 200.0}, "LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_mass","one","0");
    // add1DHist( {"Deltaeta_LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2","Delta eta of lepton and first and second jets from Higgs; Delta_eta(lepton_Higgs,jet12_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaeta","one","0");    
    // add1DHist( {"Deltaphi_LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2","Delta phi of lepton and first and second jets from Higgs; Delta_phi(lepton_Higgs,jet12_Higgs);Events;", 30, 0.0, 4.0}, "LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaphi","one","0");    
    // add1DHist( {"DeltaR_LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2","Delta R of lepton and first and second jets from Higgs; Delta_R(lepton_Higgs,jet12_Higgs);Events;", 30, 0.0, 6.0}, "LeptonFromHiggs_JetFromHiggs1_JetFromHiggs2_deltaR","one","0");
    // add1DHist( {"Tprime_mass_Gen","Mass of Gen Tprime; Mass of the Gen Tprime (GeV);Events;", 20, 0.0, 1600.0}, "LbJFromtop1_LJJFromHiggs12_mass","evWeight","0");
    // add1DHist( {"Tprime_transverse_mass_Gen","Transverse mass of Gen Tprime; Transverse mass of the Gen Tprime (GeV);Events;", 20, 0.0, 1600.0}, "LbJFromtop1_LJJFromHiggs12_transverse_mass","evWeight","0");

}

void TprimeAnalyser::setTree(TTree *t, std::string outfilename, std::string cuts, std::string Region, std::string year)
{
	if (debug)
	{
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
	}
	
	_rd = ROOT::RDataFrame(*t);
	_rlm = RNode(_rd);
	_outfilename = outfilename;
	_hist1dinfovector.clear();
	_hist2dinfovector.clear();
	_th1dhistos.clear();
	_th2dhistos.clear();
	_varstostore.clear();
	// _hist1dinfovector.clear();
	_selections.clear();

	this->setupAnalysis(cuts, Region, year);
}
//================================Selected Object Definitions====================================//

void TprimeAnalyser::setupObjects(std::string Region, std::string year)
{
  // Object selection will be defined in sequence.
  // Selected objects will be stored in new vectors.
  selectElectrons(Region);
  selectMuons(Region);
  selectLeptons(Region);
  selectMET();
  if(!_isData)
  {
    selectReconstructedLeptons();
  }
  selectJets(Region);
  if(!_isData)
  {
    selectReconstructedJets();
  }
  if(!_isData){
    this->calculateEvWeight(Region, year);
  }
}

void TprimeAnalyser::setupAnalysis(std::string cuts, std::string Region, std::string year)
{
    if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }
	
 	cout<<"year===="<< _year<< "==runtype=== " <<  _runtype <<endl;

    //==========================================event/gen/ weights==========================================//
    // Event weight for data it's always one. For MC, it depends on the sign
    //=====================================================================================================//
   
	_rlm = _rlm.Define("one", "1.0");
    // _rlm = _rlm.Define("evWeight", "1.0");
	
	if (_isData && !isDefined("evWeight"))
	{
		_rlm = _rlm.Define("evWeight", [](){
				return 1.0;
			}, {} );
	}
    if(!_isData ) // Only use genWeight
    {
        //_rlm = _rlm.Define("evWeight", "genWeight");
        
       //std::cout<<"Using evWeight = genWeight"<<std::endl;

        auto sumgenweight = _rd.Sum("genWeight");
        string sumofgenweight = Form("%f",*sumgenweight);
        _rlm = _rlm.Define("genEventSumw",sumofgenweight.c_str());
        std::cout<<"Sum of genWeights = "<<sumofgenweight.c_str()<<std::endl;
	}

  defineCuts(cuts);
  defineMoreVars(Region, year);
  if(!_isData)
  {
    defineMoreReconstructedVars();
  }
  bookHists(Region);
  if(!_isData)
  {
    bookReconstructedHists();
  }
  setupCuts_and_Hists();
  setupTree();
}
