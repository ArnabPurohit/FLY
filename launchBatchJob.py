import os
import datetime

def dir_checker(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)

year='2017'
pwdToMC = 'list_datasets/UL' + year + '/'
MCSamples = [

    # 'List_Signal600',
    # 'List_Signal625',
    # 'List_Signal650',
    # 'List_Signal675',
    'List_Signal700',
    # 'List_Signal800',
    # 'List_Signal900',
    # 'List_Signal1000',
    # 'List_Signal1100',
    # 'List_Signal1200',
    # 'List_Signal1300',
    # 'List_Signal1400',
    # 'List_Signal1500',
    # 'List_Signal1600',
    # 'List_Signal1700',
    # 'List_Signal1800',
    # 'List_TT_2L',
    # 'List_TT_SL',
    # 'List_DY_50',
    # 'List_DY_10_50',
    # 'List_ST_tW_antitop',
    # 'List_ST_tW_top',
    # 'List_ttW',
    # 'List_ttH',
    # 'List_ttZ',
    # 'List_WW',
    # 'List_WW_DS',
    # 'List_WWW',
    # 'List_WJets',
    # 'List_WH_plus',
    # 'List_WH_minus',
    # 'List_WZ',
    # 'List_ZH',
    # 'List_ZZ',
    # 'List_TT_Had',
    # 'List_ST_tW_All_antitop',
    # 'List_ST_tW_All_top',
    # 'List_TTJets_amcatnlo',
    # 'List_TTJets_madgraph',
    # 'List_TTJets_SingleLeptfromTbar',
    # 'List_TTJets_SingleLeptfromT',
    # 'List_TT_Mtop175',
]

pwdToCuts = 'list_cuts/'
Cuts = [

    # 'Dimuon_Signal',
    # 'Dimuon_CR_TT_2L',
    # 'Dimuon_CR_TT_SL',
    # 'Dimuon_CR_ttX',
    # 'Onemuononeelectron_Signal',
    # 'Onemuononeelectron_CR_TT_2L',
    # 'Onemuononeelectron_CR_TT_SL',
    # 'Onemuononeelectron_CR_ttX',
    # 'Dielectron_Signal',
    # 'Dielectron_CR_TT_2L',
    # 'Dielectron_CR_TT_SL',
    # 'Dielectron_CR_ttX',
    'Total_Signal',
    # 'Total_CR_TT_2L',
    # 'Total_CR_TT_SL',
    # 'Total_CR_ttX',
    # 'Total_Gen_Signal',

]

files = (
    (MCSamples, pwdToMC, Cuts, pwdToCuts),
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
        if '2016' in year:
            cut = cut + sample.replace('.txt','_2016.txt')

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