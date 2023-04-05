FinalStateAnalysis Package Description (miniAOD_dev version)
============================================================

The Final State Analysis (FSA) package is a CMSSW analysis framework.  
The package contains a complete implementatation to build a POG-approved 
PAT tuple, and utilities for generating plain ROOT ntuples from the PAT tuple.

**Documentation:** https://github.com/uwcms/FinalStateAnalysis/wiki


Installation
------------

```bash
cmsrel CMSSW_10_2_22
cd CMSSW_10_2_22/src/
cmsenv
git cms-init
git clone --recursive -b miniAOD_10_2_22 git@github.com:pallabidas/FinalStateAnalysis.git
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data -rf
git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
git cms-addpkg RecoMET/METFilters
git cms-addpkg  RecoJets/JetProducers
git clone -b 94X_weights_DYJets_inc_v2 git@github.com:cms-jet/PUjetID.git PUJetIDweights/
cp PUJetIDweights/weights/pileupJetId_{94,102}X_Eta* $CMSSW_BASE/src/RecoJets/JetProducers/data/
git cms-merge-topic -u alefisico:PUID_102X
source $CMSSW_BASE/src/FinalStateAnalysis/environment.sh
USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable -Wno-error=unused-function -Wno-error=sign-compare -Wno-error=reorder -Wno-error=delete-non-virtual-dtor -fpermissive -std=c++17" scram b -j 12

```
