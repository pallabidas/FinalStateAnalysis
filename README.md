FinalStateAnalysis Package Description (miniAOD_dev version)
============================================================

The Final State Analysis (FSA) package is a CMSSW analysis framework.  
The package contains a complete implementatation to build a POG-approved 
PAT tuple, and utilities for generating plain ROOT ntuples from the PAT tuple.

**Documentation:** https://github.com/uwcms/FinalStateAnalysis/wiki


Installation
------------

Current recommended UL CMSSW version: ``CMSSW_10_6_30``.

Get release area:

```bash
  cmsrel CMSSW_10_6_30
  cd CMSSW_10_6_30/src
  # Setup your CMSSW environment
  cmsenv
  git cms-init
```

Checkout the FinalStateAnalysis repository:

```bash
  git clone -b UL_10_6_30 https://github.com/pallabidas/FinalStateAnalysis.git
  git cms-merge-topic cms-egamma:EgammaPostRecoTools
  git cms-addpkg EgammaAnalysis/ElectronTools
  rm -rf EgammaAnalysis/ElectronTools/data
  git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
  git cms-addpkg  RecoJets/JetProducers
  git clone -b 94X_weights_DYJets_inc_v2 git@github.com:cms-jet/PUjetID.git PUJetIDweights/
  cp PUJetIDweights/weights/pileupJetId_{94,102}X_Eta* $CMSSW_BASE/src/RecoJets/JetProducers/data/
  
  cd FinalStateAnalysis
  source environment.sh
```

Compile:

```bash
  USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable -Wno-error=unused-function -Wno-error=sign-compare -Wno-error=reorder -Wno-error=delete-non-virtual-dtor -fpermissive -std=c++17" scram b -j 12
```

Run:

```bash
  cd FinalStateAnalysis/NtupleTools/test/
  cmsRun make_ntuples_cfg.py htt=1 era="2018" isMC=1 isEmbedded=0 skipMET=1 maxEvents=200 paramFile=../python/parameters/ztt.py runningLocal=1 fullJES=1 metShift=1 inputFiles=root://cms-xrd-global.cern.ch://store/mc/RunIIAutumn18MiniAOD/SUSYGluGluToHToAA_AToBB_AToTauTau_M-12_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/A734E146-A24C-7549-AF2E-66B2E809DBA0.root
```
