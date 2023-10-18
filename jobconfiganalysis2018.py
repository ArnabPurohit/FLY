"""
File contains job options 
"""

# options for Nanoaodrdframe
config = {
        # tree name of input file(s)
        'intreename': "Events",
        # tree name of output file(s) it cannot be the same as the input tree name or it'll crash
        'outtreename': "outputTree",
        #data year (2016,2017,2018)
        'year': 2018,
        # is ReReco or Ultra Legacy
        'runtype': 'UL',
        'datatype': -1,


        #for correction
        # good json file 
        'goodjson': 'data/Legacy_RunII/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt', #UL_2018

        # pileup weight for MC 
        'pileupfname': 'data/LUM/2018_UL/puWeights.json',

        'pileuptag': 'Collisions18_UltraLegacy_goldenJSON',

        # json filename for BTV correction 
        'btvfname': 'data/BTV/2018_UL/btagging.json',
        # BTV correction type: //to use btvtype from the json file for the btag SFs
        # 'btvtype': 'deepJet_mujets', #for fixed wp : case1 
        'btvtype': 'deepJet_comb', #for fixed wp : case2
        #'btvtype': 'deepJet_shape', #for central shape : case3
        
        'muon_fname': 'data/MUON/2018_UL/muon_Z.json', 
        'muontype': 'NUM_TightID_DEN_TrackerMuons',#'Medium ID UL scale factor',
        # 'muontype': 'NUM_TightRelIso_DEN_MediumID',#'Medium ISO UL scale factor',

        # json file name for JERC
        'jercfname': 'data/JERC/jetUL18_jerc.json',

        # conbined correction type for jets
        'jerctag': 'Summer19UL18_V5_MC_L1L2L3Res_AK4PFchs', 
        # JER/jet resolution : Summer19UL18_JRV2_MC_PtResolution_AK4PFchs
        # jet uncertainty 
        'jercunctag': 'Summer19UL18_V5_MC_Total_AK4PFchs',
        # JER/jetresolution scale factor: Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs
        
        }

# processing options
procflags = {
        # how many jobs?
        'split': 'Max',
        # if False, one output file per input file, if True then one output file for everything
        'allinone': False, 
        # if True then skip existing analyzed files
        'skipold': True,
        # travel through the subdirectories and their subdirecties when processing.
        # becareful not to mix MC and real DATA in them.
        'recursive': True,
        # if False then only selected branches which is done in the .cpp file will be saved
        'saveallbranches': False,
        #How many input files?
        'nrootfiles': 10000,
        }