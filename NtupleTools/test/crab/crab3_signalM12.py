from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SUSYGluGluToHToAA_AToBB_AToTauTau_M-12'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 300
config.JobType.maxMemoryMB = 2500
config.JobType.psetName = 'make_ntuples_cfg_crab.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['out.root']
config.JobType.inputFiles = ['Summer19UL18_V5_DATA.db', 'Summer19UL18_V5_MC.db', 'Summer19UL17_RunBCDEF_V5_DATA.db', 'Summer19UL17_V5_MC.db', 'Summer19UL17_RunBCDEF_V5_DATA.db', 'Summer19UL16_V7_MC.db', 'Summer19UL16_RunBCDEFGH_Combined_V7_DATA.db', 'Summer19UL16APV_V7_MC.db', 'RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt']

config.Data.inputDataset = '/SUSYGluGluToHToAA_AToBB_AToTauTau_M-12_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/pdas/'

config.Site.storageSite = 'T2_US_Wisconsin'
