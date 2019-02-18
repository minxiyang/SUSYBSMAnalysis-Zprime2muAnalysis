
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_electrons_2016_ADDGravToLL_LambdaT10000_M1700'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/ADDGravToLL_LambdaT-10000_M-1700_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_electrons_2016_ADDGravToLL_LambdaT10000_M1700'
config.Data.outLFNDirBase = '/store/user/zhangfa/ADD2016MC'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
#config.Site.whitelist = ["T2_IT_Bari"]
#config.Site.storageSite = 'T2_US_Purdue'
config.Site.whitelist = ["T2_US_Caltech", "T2_US_USCD", "T2_US_Nebraska"]
config.Site.storageSite = 'T3_US_FNALLPC'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

