# Embed IDs for jets
import FWCore.ParameterSet.Config as cms
import os

def preJets(process, jSrc, jupSrc, jdownSrc, vSrc, metSrc,mSrc, eSrc, **kwargs):
    postfix = kwargs.pop('postfix','')
    jType = kwargs.pop('jType','AK4PFchs')
    doBTag = kwargs.pop('doBTag',False)
    doFullJESUnc = kwargs.pop('doFullJESUnc',False)
    runningLocal = kwargs.pop('runningLocal',False)


    mod = cms.EDProducer(
        "MiniAODJetIdEmbedder",
        src = cms.InputTag(jSrc),
        genParticles = cms.InputTag("prunedGenParticles"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
        preprocess_json = cms.FileInPath("FinalStateAnalysis/PatTools/data/preprocess.json"),
        model_path = cms.FileInPath("FinalStateAnalysis/PatTools/data/model.onnx"),
    )
    modName = 'miniPatJets{0}'.format(postfix)
    setattr(process,modName,mod)
    jSrc = modName

    pathName = 'runMiniAODJetEmbedding{0}'.format(postfix)
    setattr(process,pathName,cms.Path(getattr(process,modName)))
    process.schedule.append(getattr(process,pathName))

    if doFullJESUnc :
        if runningLocal : fName = "../../NtupleTools/data/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt" 
        else :
            fName = "RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt"

        modName = 'miniAODJetFullSystematicsEmbedding{0}'.format(postfix)
        mod = cms.EDProducer(
	    "MiniAODJetFullSystematicsEmbedder",
            src = cms.InputTag(jSrc),
            srcMET=cms.InputTag(metSrc),
            corrLabel = cms.string(jType),
            fName = cms.string(fName)
        )
        jSrc = modName
        setattr(process,modName,mod)

        pathName = 'jetFullSystematicsEmbedding{0}'.format(postfix)
        path = cms.Path(getattr(process,modName))
        setattr(process,pathName,path)

        print modName+" for  MET?" 
 
        process.schedule.append(getattr(process,pathName))


    print jSrc 

    modName = 'miniAODJetSystematicsEmbedding{0}'.format(postfix)
    mod = cms.EDProducer(
	"MiniAODJetSystematicsEmbedder",
        src = cms.InputTag(jSrc),
        corrLabel = cms.string(jType)
    )
    jSrc = modName
    setattr(process,modName,mod)
    pathName = 'jetSystematicsEmbedding{0}'.format(postfix)
    path = cms.Path(getattr(process,modName))
    setattr(process,pathName,path)
    process.schedule.append(getattr(process,pathName))

    modName = 'miniAODJERSystematicsEmbedding{0}'.format(postfix)
    mod = cms.EDProducer(
        "MiniAODJERSystematicsEmbedder",
        src = cms.InputTag(jSrc),
        up = cms.InputTag(jupSrc),
        down = cms.InputTag(jdownSrc)
    )
    jSrc = modName
    setattr(process,modName,mod)
    pathName = 'jerSystematicsEmbedding{0}'.format(postfix)
    path = cms.Path(getattr(process,modName))
    setattr(process,pathName,path)
    process.schedule.append(getattr(process,pathName))

    return jSrc

