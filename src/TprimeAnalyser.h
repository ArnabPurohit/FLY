/*
 * TprimeAnalyser.h
 *
 *  Created on: May 6, 2022
 *      Author: suyong
 *		Developper: cdozen
 */

#ifndef BASEANALYSER_H_
#define BASEANALYSER_H_

#include "NanoAODAnalyzerrdframe.h"

class TprimeAnalyser: public NanoAODAnalyzerrdframe
{

  	public:
  	TprimeAnalyser(TTree *t, std::string outfilename, std::string cuts, std::string Region, std::string year);
  	void defineCuts(std::string cuts);		//define a series of cuts from defined variables only. you must implement this in your subclassed analysis 
  	void defineMoreVars(std::string Region, std::string year); 	//define higher-level variables from basic ones, you must implement this in your subclassed analysis code
  	void defineMoreReconstructedVars();
  	void bookHists(std::string Region); 		//book histograms, you must implement this in your subclassed analysis code
  	void bookReconstructedHists();
  	void setTree(TTree *t, std::string outfilename, std::string cuts, std::string Region, std::string year);
  	void setupObjects(std::string Region, std::string year);
  	void setupAnalysis(std::string cuts, std::string Region, std::string year);
  	// object selectors
  	void selectElectrons(std::string Region);
  	void selectMuons(std::string Region);
  	void selectLeptons(std::string Region);
  	void selectMET();
  	void selectReconstructedLeptons();
  	void selectJets(std::string Region);
  	void calculateEvWeight(std::string Region, std::string year);
  	void selectReconstructedJets();

 	bool debug = true;
	bool _jsonOK;
	string _outfilename;
	string _putag;

	TFile *_outrootfile;
	vector<string> _outrootfilenames;

};

#endif /* BASEANALYSER_H_ */
