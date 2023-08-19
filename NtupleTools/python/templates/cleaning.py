'''

Ntuple branch template sets for applying cleaning and extra object vetoes

Each string is transformed into an expression on a FinalStateEvent object.

{object} should be replaced by an expression which evaluates to a pat::Muon
i.e. daughter(1) or somesuch.

Author: Evan K. Friis

'''

from FinalStateAnalysis.Utilities.cfgtools import addargs

# Vetos on extra stuff in the event
vetos = addargs(
    #MUON VETOS

    muVetoZTTp001dxyzR0 = 'vetoMuons(0.0, "pt > 10 & abs(eta) < 2.4 & ( ( pfIsolationR04().sumChargedHadronPt + max( pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5 * pfIsolationR04().sumPUPt, 0.0)) / pt() ) < 0.3 & isMediumMuon > 0 & abs( userFloat(\'ipDXY\') ) < 0.045 & abs( userFloat(\'dz\') ) < 0.2").size()',
    dimuonVeto = 'vetoSecondMuon(0.15,"pt > 15 & abs(eta) < 2.4 & ( ( pfIsolationR04().sumChargedHadronPt + max( pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5 * pfIsolationR04().sumPUPt, 0.0)) / pt() ) < 0.3 & isGlobalMuon > 0 & isTrackerMuon > 0 & isPFMuon > 0 & abs( userFloat(\'ipDXY\') ) < 0.045 & abs( userFloat(\'dz\') ) < 0.2").size()',
    
    ##ELECTRON VETOS

    eVetoZTTp001dxyzR0 = 'vetoElectrons(0.0, "pt > 10 & abs(eta) < 2.5 & ( pfIsolationVariables().sumChargedHadronPt + max(0.0,pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt - userFloat(\'rho_fastjet\')*userFloat(\'EffectiveArea\'))) / pt() < 0.3 & electronID(\'mvaEleID-Fall17-noIso-V2-wp90\') > 0 & passConversionVeto() > 0 & abs( userFloat(\'ipDXY\') ) < 0.045 & abs( userFloat(\'dz\') ) < 0.2").size()',
    dielectronVeto = 'vetoSecondElectron(0.15, "pt > 15 & abs(eta) < 2.5 & ( pfIsolationVariables().sumChargedHadronPt + max(0.0,pfIsolationVariables().sumNeutralHadronEt + pfIsolationVariables().sumPhotonEt - userFloat(\'rho_fastjet\')*userFloat(\'EffectiveArea\'))) / pt() < 0.3 & electronID(\'cutBasedElectronID-Fall17-94X-V2-veto\') > 0 & abs( userFloat(\'ipDXY\') ) < 0.045 & abs( userFloat(\'dz\') ) < 0.2").size()',

    #(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & 
    #B-JET Vetos
    bjetDeepCSVVeto20Loose_2016_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.2217").size()',
    bjetDeepCSVVeto20Medium_2016_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.6321").size()',
    bjetDeepCSVVeto20Tight_2016_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.8953").size()',
    bjetDeepCSVVeto20Medium_2016_DR0 = 'vetoJets(0.0, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.6321").size()',

    bjetDeepCSVVeto20Loose_2017_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1522").size()',
    bjetDeepCSVVeto20Medium_2017_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.4941").size()',
    bjetDeepCSVVeto20Tight_2017_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.8001").size()',
    bjetDeepCSVVeto20Medium_2017_DR0 = 'vetoJets(0.0, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.4941").size()',

    bjetDeepCSVVeto20Loose_2018_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241").size()',
    bjetDeepCSVVeto20Medium_2018_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.4184").size()',
    bjetDeepCSVVeto20Tight_2018_DR0p5 = 'vetoJets(0.5, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.7527").size()',
    bjetDeepCSVVeto20Medium_2018_DR0 = 'vetoJets(0.0, "pt > 25 & abs(eta) < 2.4 & userFloat(\'idTight\') > 0.5 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.4184").size()',


    # Leading and subleading BTagged Jets (ordered by btag score)
    deepcsvb1_pt = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(0)',
    deepcsvb1_eta = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(1)',
    deepcsvb1_phi = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(2)',
    deepcsvb1_m = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(3)',
    deepcsvb1_btagscore = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(4)',
    deepcsvb1_hadronflavour = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(5)',
    deepcsvb2_pt = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(6)',
    deepcsvb2_eta = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(7)',
    deepcsvb2_phi = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(8)',
    deepcsvb2_m = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(9)',
    deepcsvb2_btagscore = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(10)',
    deepcsvb2_hadronflavour = 'deepCSVJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepCSVJetTags:probb\') + bDiscriminator(\'pfDeepCSVJetTags:probbb\')) > 0.1241", 0.5).at(11)',

    deepflavourb1_pt = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(0)',
    deepflavourb1_eta = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(1)',
    deepflavourb1_phi = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(2)',
    deepflavourb1_m = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(3)',
    deepflavourb1_btagscore = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(4)',
    deepflavourb1_hadronflavour = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(5)',
    deepflavourb2_pt = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(6)',
    deepflavourb2_eta = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(7)',
    deepflavourb2_phi = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(8)',
    deepflavourb2_m = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(9)',
    deepflavourb2_btagscore = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(10)',
    deepflavourb2_hadronflavour = 'deepFlavourJetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & (bDiscriminator(\'pfDeepFlavourJetTags:probb\') + bDiscriminator(\'pfDeepFlavourJetTags:probbb\') + bDiscriminator(\'pfDeepFlavourJetTags:problepb\')) > 0.0494 & userFloat(\'matchesGenBB\')", 0.5).at(11)',

    mergedb_pt = 'jetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & userFloat(\'matchesGenBB\') > 0.", 0.5).at(0)',
    mergedb_eta = 'jetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & userFloat(\'matchesGenBB\') > 0.", 0.5).at(1)',
    mergedb_phi = 'jetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & userFloat(\'matchesGenBB\') > 0.", 0.5).at(2)',
    mergedb_deepcsvscore = 'jetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & userFloat(\'matchesGenBB\') > 0.", 0.5).at(3)',
    mergedb_deepflavourscore = 'jetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & userFloat(\'matchesGenBB\') > 0.", 0.5).at(4)',
    mergedb_hadronflavour = 'jetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & userFloat(\'matchesGenBB\') > 0.", 0.5).at(5)',
    mergedb_pnetScore = 'jetVariables("pt > 20 & userFloat(\'idTight\') > 0.5 & abs(eta) < 2.4 & userFloat(\'matchesGenBB\') > 0.", 0.5).at(6)',

    ## Leading and sublead jets
    ### Comment for 2017
    j1pt = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(0)',
    j1eta = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(1)',
    j1phi = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(2)',
    j1_deepcsvscore = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(3)',
    j1_deepflavourscore = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(4)',
    j1_hadronflavour = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(5)',
    j1_pnetScore = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(6)',
    j2pt = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(7)',
    j2eta = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(8)',
    j2phi = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(9)',
    j2_deepcsvscore = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(10)',
    j2_deepflavourscore = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(11)',
    j2_hadronflavour = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(12)',
    j2_pnetScore = 'jetVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(13)',

    jetVeto20 = 'vetoJets(0.5, "pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)").size()',
    jetVeto30 = 'vetoJets(0.5, "pt > 30 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)").size()',
    jetVeto30_JetEnUp = 'vetoJets(0.5, "userCand(\'jes+\').pt > 30 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)").size()',
    jetVeto30_JetEnDown = 'vetoJets(0.5, "userCand(\'jes-\').pt > 30 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)").size()',
    vbfMass = 'vbfVariables("pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).mass',


    ## Uncomment for 2017 Noisy jets
    #j1pt = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(0)',
    #j1eta = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(1)',
    #j1phi = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(2)',
    #j1_deepcsvscore = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(3)',
    #j1_deepflavourscore = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(4)',
    #j1_hadronflavour = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(5)',
    #j2pt = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(6)',
    #j2eta = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(7)',
    #j2phi = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(8)',
    #j2_deepcsvscore = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(9)',
    #j2_deepflavourscore = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(10)',
    #j2_hadronflavour = 'jetVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).at(11)',
    #jetVeto20 = 'vetoJets(0.5, "(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)").size()',
    #jetVeto30 = 'vetoJets(0.5, "(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 30 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)").size()',
    #jetVeto30_JetEnUp = 'vetoJets(0.5, "(userCand(\'jes+\').pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & userCand(\'jes+\').pt > 30 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)").size()',
    #jetVeto30_JetEnDown = 'vetoJets(0.5, "(userCand(\'jes-\').pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & userCand(\'jes-\').pt > 30 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)").size()',
    #vbfMass = 'vbfVariables("(pt > 50 | (abs(eta) < 2.65 | abs(eta) > 3.139)) & pt > 20 & abs(eta) < 4.7 & userFloat(\'idTight\') > 0.5 & (userInt(\'pileupJetId:fullId\')>3 | pt>50)", 0.5).mass',


)

overlaps = addargs(
   
)
