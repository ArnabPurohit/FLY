import os
import datetime

def dir_checker(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)

year='2018'
pwd = 'list_datasets/UL' + year + '/'
Samples = [

    'List_Data_DoubleMuon_UL2018A',
    'List_Data_DoubleMuon_UL2018B',
    'List_Data_DoubleMuon_UL2018C',
    'List_Data_DoubleMuon_UL2018D',
    # 'List_Data_SingleMuon_UL2018A',
    # 'List_Data_SingleMuon_UL2018B',
    # 'List_Data_SingleMuon_UL2018C',
    # 'List_Data_SingleMuon_UL2018D',
    # 'List_Data_MuonEG_UL2018A',
    # 'List_Data_MuonEG_UL2018B',
    # 'List_Data_MuonEG_UL2018C',
    # 'List_Data_MuonEG_UL2018D',
    # 'List_Data_EGamma_UL2018A',
    # 'List_Data_EGamma_UL2018B',
    # 'List_Data_EGamma_UL2018C',
    # 'List_Data_EGamma_UL2018D',
]

pwdToCuts = 'list_cuts/'
Cuts = [

    'Dimuon_Signal_data',
    'Dimuon_CR_TT_2L_data',
    'Dimuon_CR_TT_SL_data',
    'Dimuon_CR_ttX_data',
    # 'Onemuononeelectron_Signal_data',
    # 'Onemuononeelectron_CR_TT_2L_data',
    # 'Onemuononeelectron_CR_TT_SL_data',
    # 'Onemuononeelectron_CR_ttX_data',
    # 'Dielectron_Signal_data',
    # 'Dielectron_CR_TT_2L_data',
    # 'Dielectron_CR_TT_SL_data',
    # 'Dielectron_CR_ttX_data',
    # 'Total_Signal_data',
    # 'Total_CR_TT_2L_data',
    # 'Total_CR_TT_2L_data',
    # 'Total_CR_TT_SL_data',
    # 'Total_CR_ttX_data',
    
]

files = (
    (Samples, pwd, Cuts, pwdToCuts),
    )


tasks = 'All'
cpus_per_task = '1'
copyInstance = True

########## 0 = MC ; 1 = DATA #############
n = 0

for sample in files[n][0]:

    if len(files[n][0]) == 0: break

    analysisOutput = 'analyzed/' + str(datetime.datetime.now().strftime('%Y_%m_%d_%H')) + '/'
    dir_checker(analysisOutput)

    if copyInstance: os.system('cp src/TprimeAnalyser.cpp ' + analysisOutput + 'TprimeAnalyser_COPY.cpp')

    analysisOutput = 'analyzed/' + str(datetime.datetime.now().strftime('%Y_%m_%d_%H')) + '/UL' + year + '/'
    dir_checker(analysisOutput)

    analysisOutput = 'analyzed/' + str(datetime.datetime.now().strftime('%Y_%m_%d_%H')) + '/UL' + year + '/' + sample.replace('List_','') + '/'
    dir_checker(analysisOutput)

    analysisOutput_temp = analysisOutput

    analysisOutput_merged = 'analyzed/' + str(datetime.datetime.now().strftime('%Y_%m_%d_%H')) + '/Merged/'
    dir_checker(analysisOutput_merged)

    sample += '.txt'

    for cuts in files[n][2]:
        analysisOutput = analysisOutput_temp + sample.replace('List_','').replace('.txt','') + '_' + cuts + '/'
        dir_checker(analysisOutput)
        dir_checker(analysisOutput + 'slurmOutputs')

        Region = 'Signal'
        if('CR_TT_2L' in cuts):
            Region = 'CR_TT_2L'
        elif('CR_TT_SL' in cuts):
            Region = 'CR_TT_SL'
        elif('CR_ttX' in cuts):
            Region = 'CR_ttX'
        print(Region)

        cuts += '.txt'

        if copyInstance: os.system('cp list_cuts/' + cuts + ' ' + analysisOutput + cuts.replace('.txt','_COPY.txt'))

        if tasks == 'All':
            with open(files[n][1] + sample, 'r') as fToRead:
                ntasks = str(len(fToRead.readlines()))
        else: ntasks = tasks

        slurmConfigs = [

            ' --cpus-per-task=' + cpus_per_task,
            ' --ntasks=' + ntasks,
            ' --job-name=' + sample.replace('List_','').replace('.txt','') + '_' + cuts.replace('.txt',''),
            ' --output=' + analysisOutput + 'slurmOutputs/output.out',
            ' --error=' + analysisOutput + 'slurmOutputs/errors.err'
        ]

        args = [
            ' ' + files[n][1] + sample, #inputFile
            ' ' + analysisOutput, #outputFiles
            ' ' + files[n][3] + cuts, #list of Cuts
            ' ' + str(Region), #which Region is selected (Signal, CR_TT_2L, CR_TT_SL, CR_ttX)
            ' ' + year,
        ]

        os.system('sbatch' + ''.join(slurmConfigs) + ' submitBatchJob.sh' + ''.join(args))