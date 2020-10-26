
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_muons_2018_ZprimeCT14_M30000_6000toInfpythia'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/ZprimeToMuMuCT14-30000_M-6000ToInf_13TeV-pythia8/minxi-RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v8-59bbcc9933f31ad81e7726585f9ffceb/USER'
config.Data.inputDBS = 'phys03'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_muons_2018_ZprimeCT14_M30000_6000toInfpythia'
config.Data.outLFNDirBase = '/store/user/minxi/'
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
#config.General.instance = 'preprod' 
config.Site.whitelist = ["T2_US_*"]
config.Site.blacklist = ['T2_US_Caltech']
config.Site.storageSite = 'T2_US_Purdue'
#config.JobType.maxMemoryMB  = 4000
config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

