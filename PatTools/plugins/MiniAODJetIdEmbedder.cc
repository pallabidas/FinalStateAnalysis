/*
 * Embed PF Jet IDs (see https://twiki.cern.ch/twiki/bin/view/CMS/JetID)
 * into pat::Jets
 *
 * Author: Evan K. Friis, UW Madison
 */


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"

class MiniAODJetIdEmbedder : public edm::EDProducer {
  public:
    MiniAODJetIdEmbedder(const edm::ParameterSet& pset);
    virtual ~MiniAODJetIdEmbedder(){}
    void produce(edm::Event& evt, const edm::EventSetup& es);
  private:
    edm::EDGetTokenT<edm::View<pat::Jet> > srcToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
};

MiniAODJetIdEmbedder::MiniAODJetIdEmbedder(const edm::ParameterSet& pset) {
  srcToken_ = consumes<edm::View<pat::Jet> >(pset.getParameter<edm::InputTag>("src"));
  genParticlesToken_ = consumes<edm::View<reco::GenParticle> >(pset.getParameter<edm::InputTag>("genParticles"));
  produces<pat::JetCollection>();
}

void MiniAODJetIdEmbedder::produce(edm::Event& evt, const edm::EventSetup& es) {
  std::unique_ptr<pat::JetCollection> output(new pat::JetCollection);

  edm::Handle<edm::View<pat::Jet> > input;
  evt.getByToken(srcToken_, input);

  output->reserve(input->size());

  edm::Handle<edm::View<reco::GenParticle> > genPart;
  evt.getByToken(genParticlesToken_, genPart);
  TLorentzVector v1;
  TLorentzVector v2;
  //TLorentzVector higgs;
  for ( auto& genp : *genPart ) {
    //if(genp.pdgId() == 25 && genp.status() == 62) higgs.SetPtEtaPhiM(gen.pt(), genp.eta(), genp.phi(), 125.09);
    if(genp.pdgId() == 5 && genp.status() == 23) v1.SetPtEtaPhiM(genp.pt(), genp.eta(), genp.phi(), 4.18);
    if(genp.pdgId() == -5 && genp.status() == 23) v2.SetPtEtaPhiM(genp.pt(), genp.eta(), genp.phi(), 4.18);
  }
  size_t matchedindex = 99;
  float dR_genBs = v1.DeltaR(v2);
  if(dR_genBs < 0.8){
    TLorentzVector diB = v1 + v2;
    size_t jetindex = 0;
    edm::View<pat::Jet>::const_iterator jet;
    float deltar = 0.4;
    const edm::View<pat::Jet> & jets = *input;
    for(jet = jets.begin(); jet != jets.end(); ++jet){
      if(reco::deltaR(jet->eta(), jet->phi(), diB.Eta(), diB.Phi()) < deltar){
        deltar = reco::deltaR(jet->eta(), jet->phi(), diB.Eta(), diB.Phi());
        matchedindex = jetindex;
      }
      jetindex++;
    }
  }

  for (size_t i = 0; i < input->size(); ++i) {
    pat::Jet jet = input->at(i);
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID
    bool loose = true;
    bool tight = true;
    bool tightLepVeto = true;
    if (std::abs(jet.eta()) <= 2.6) {
      if (jet.neutralHadronEnergyFraction() >= 0.99) {
        loose = false;
      }
      if (jet.neutralHadronEnergyFraction() >= 0.90) {
        tight = false;
        tightLepVeto = false;
      }
      if (jet.neutralEmEnergyFraction() >= 0.99) {
        loose = false;
      }
      if (jet.neutralEmEnergyFraction() >= 0.90) {
        tight = false;
        tightLepVeto = false;
      }
      if (jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1){
        loose = false;
        tight = false;
        tightLepVeto = false;
      }
      if (jet.muonEnergyFraction() >= 0.8)
        {
          tightLepVeto = false;
        }
      if (jet.chargedHadronEnergyFraction() == 0) {
        loose = false;
        tight = false;
        tightLepVeto = false;
      }
      if (jet.chargedHadronMultiplicity() == 0) {
        loose = false;
        tight = false;
        tightLepVeto = false;
      }
      if (jet.chargedEmEnergyFraction() >= 0.80) {
        tightLepVeto = false;
      }
    }

    if (std::abs(jet.eta()) >2.6 && std::abs(jet.eta()) <= 2.7) {
      if (jet.neutralHadronEnergyFraction() >= 0.99) {
        loose = false;
      }
      if (jet.neutralHadronEnergyFraction() >= 0.90) {
        tight = false;
        tightLepVeto = false;
      }
      if (jet.neutralEmEnergyFraction() >= 0.99) {
        loose = false;
        tight = false;
        tightLepVeto = false;
      }
      if (jet.muonEnergyFraction() >= 0.8){
        tightLepVeto = false;
      }
      if (jet.chargedHadronMultiplicity() == 0) {
        loose = false;
        tight = false;
        tightLepVeto = false;
      }
      if (jet.chargedEmEnergyFraction() >= 0.80) {
        tightLepVeto = false;
      }
    }

    if (std::abs(jet.eta()) > 2.7 && std::abs(jet.eta()) <= 3.0) {
      if (jet.neutralEmEnergyFraction() >= 0.99 or jet.neutralEmEnergyFraction()<=0.02) {
        loose = false;
        tight = false;
      }
      if (jet.neutralMultiplicity()<=2) {
        loose = false;
        tight = false;
      }
    }

    if (std::abs(jet.eta()) > 3.0) {
      if (jet.neutralEmEnergyFraction() >= 0.90) {
        loose = false;
        tight = false;
      }
      if (jet.neutralHadronEnergyFraction() <= 0.2) {
        tight = false;
        tightLepVeto = false;
      }
      if (jet.neutralMultiplicity()<=10) {
        loose = false;
        tight = false;
      }
    }
    jet.addUserFloat("idLoose", loose);
    jet.addUserFloat("idTight", tight);
    jet.addUserFloat("idTightLepVeto", tightLepVeto);

    // Pileup discriminant
    bool passPU = true;
    float jpumva = jet.userFloat("pileupJetId:fullDiscriminant");
    if(jet.pt() > 20)
      {
        if(fabs(jet.eta()) > 3.)
          {
            if(jpumva <= -0.45) passPU = false;
          }
        else if(fabs(jet.eta()) > 2.75)
          {
            if(jpumva <= -0.55) passPU = false;
          }
        else if(fabs(jet.eta()) > 2.5)
          {
            if(jpumva <= -0.6) passPU = false;
          }
        else if(jpumva <= -0.63) passPU = false;
      }
    else
      {
        if(fabs(jet.eta()) > 3.)
          {
            if(jpumva <= -0.95) passPU = false;
          }
        else if(fabs(jet.eta()) > 2.75)
          {
            if(jpumva <= -0.94) passPU = false;
          }
        else if(fabs(jet.eta()) > 2.5)
          {
            if(jpumva <= -0.96) passPU = false;
          }
        else if(jpumva <= -0.95) passPU = false;
      }

    jet.addUserFloat("puID", float(passPU));


    // Add matching to merged gen bs
    bool passMatching = false;
    if(i == matchedindex) passMatching = true;

    jet.addUserFloat("matchesGenBB", float(passMatching));
    output->push_back(jet);
  }

  evt.put(std::move(output));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MiniAODJetIdEmbedder);
