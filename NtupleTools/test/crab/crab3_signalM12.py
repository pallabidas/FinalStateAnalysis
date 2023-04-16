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
config.JobType.inputFiles = ['Autumn18_RunABCD_V19_DATA.db', 'Autumn18_V19_MC.db', 'Fall17_17Nov2017_V32_94X_DATA.db', 'Fall17_17Nov2017_V32_94X_MC.db', 'Summer16_07Aug2017All_V11_DATA.db', 'Summer16_07Aug2017_V11_MC.db', 'Summer16_23Sep2016AllV4_DATA.db', 'Summer16_23Sep2016V4_MC.db', 'RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt']

config.Data.inputDataset = '/SUSYGluGluToHToAA_AToBB_AToTauTau_M-12_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/pdas/'

config.Site.storageSite = 'T2_US_Wisconsin'
