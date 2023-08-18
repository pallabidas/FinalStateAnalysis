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
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <nlohmann/json.hpp>

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace cms::Ort;
using namespace btagbtvdeep;

class MiniAODJetIdEmbedder : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
  public:
    MiniAODJetIdEmbedder(const edm::ParameterSet& pset, const ONNXRuntime *);
    virtual ~MiniAODJetIdEmbedder(){}
    void produce(edm::Event& evt, const edm::EventSetup& es);
    static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &);
    static void globalEndJob(const ONNXRuntime *);

  private:
    edm::EDGetTokenT<edm::View<pat::Jet> > srcToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> svToken_;

    std::vector<std::string> input_names_; // names of each input group - the ordering is important!
    std::vector<std::vector<int64_t>> input_shapes_; // shapes of each input group (-1 for dynamic axis)
    std::vector<unsigned> input_sizes_; // total length of each input vector
    std::unordered_map<std::string, PreprocessParams> prep_info_map_; // preprocessing info for each input group
    FloatArrays data_; // each stream hosts its own data

    const reco::Vertex *pv_ = nullptr;
    std::vector<float> jet_pfcand_px;
    std::vector<float> jet_pfcand_py;
    std::vector<float> jet_pfcand_pz;
    std::vector<float> jet_pfcand_energy;
    std::vector<float> jet_pfcand_pt;
    std::vector<float> jet_pfcand_pt_log;
    std::vector<float> jet_pfcand_e_log;
    std::vector<float> jet_pfcand_phirel;
    std::vector<float> jet_pfcand_etarel;
    std::vector<float> jet_pfcand_abseta;
    std::vector<float> jet_pfcand_puppiw;
    std::vector<float> jet_pfcand_charge;
    std::vector<float> jet_pfcand_isEl;
    std::vector<float> jet_pfcand_isMu;
    std::vector<float> jet_pfcand_isChargedHad;
    std::vector<float> jet_pfcand_isGamma;
    std::vector<float> jet_pfcand_isNeutralHad;
    std::vector<float> jet_pfcand_hcalFrac;
    std::vector<float> jet_pfcand_hcalFracCalib;
    std::vector<float> jet_pfcand_VTX_ass;
    std::vector<float> jet_pfcand_fromPV;
    std::vector<float> jet_pfcand_lostInnerHits;
    std::vector<float> jet_pfcand_trackHighPurity;
    std::vector<float> jet_pfcand_dz;
    std::vector<float> jet_pfcand_dzsig;
    std::vector<float> jet_pfcand_dxy;
    std::vector<float> jet_pfcand_dxysig;
    std::vector<float> jet_pfcand_normchi2;
    std::vector<float> jet_pfcand_quality;
    std::vector<float> jet_pfcand_btagEtaRel;
    std::vector<float> jet_pfcand_btagPtRatio;
    std::vector<float> jet_pfcand_btagPParRatio;
    std::vector<float> jet_pfcand_btagSip2dVal;
    std::vector<float> jet_pfcand_btagSip2dSig;
    std::vector<float> jet_pfcand_btagSip3dVal;
    std::vector<float> jet_pfcand_btagSip3dSig;
    std::vector<float> jet_pfcand_btagJetDistVal;
    std::vector<float> jet_sv_pt;
    std::vector<float> jet_sv_pt_log;
    std::vector<float> jet_sv_ptrel;
    std::vector<float> jet_sv_ptrel_log;
    std::vector<float> jet_sv_eta;
    std::vector<float> jet_sv_phi;
    std::vector<float> jet_sv_mass;
    std::vector<float> jet_sv_energy;
    std::vector<float> jet_sv_energy_log;
    std::vector<float> jet_sv_erel;
    std::vector<float> jet_sv_erel_log;
    std::vector<float> jet_sv_deta;
    std::vector<float> jet_sv_dphi;
    std::vector<float> jet_sv_chi2;
    std::vector<float> jet_sv_dxy;
    std::vector<float> jet_sv_dxysig;
    std::vector<float> jet_sv_d3d;
    std::vector<float> jet_sv_d3dsig;
    std::vector<float> jet_sv_ntrack;

    const std::vector<std::string>   pf_points_{"jet_pfcand_etarel",
                                                "jet_pfcand_phirel"};
    const std::vector<std::string> pf_features_{"jet_pfcand_pt_log",
                                                "jet_pfcand_pt_log",
                                                "jet_pfcand_e_log",
                                                "jet_pfcand_etarel",
                                                "jet_pfcand_phirel",
                                                "jet_pfcand_abseta",
                                                "jet_pfcand_charge",
                                                "jet_pfcand_VTX_ass",
                                                "jet_pfcand_lostInnerHits",
                                                "jet_pfcand_normchi2",
                                                "jet_pfcand_quality",
                                                "jet_pfcand_dz",
                                                "jet_pfcand_dzsig",
                                                "jet_pfcand_dxy",
                                                "jet_pfcand_dxysig",
                                                "jet_pfcand_btagEtaRel",
                                                "jet_pfcand_btagPtRatio",
                                                "jet_pfcand_btagPParRatio",
                                                "jet_pfcand_btagSip3dVal",
                                                "jet_pfcand_btagSip3dSig",
                                                "jet_pfcand_btagJetDistVal",
                                                "jet_pfcand_puppiw",
                                                "jet_pfcand_isEl",
                                                "jet_pfcand_isMu"};
    const std::vector<std::string>     pf_mask_{"pfcand_mask"};
    const std::vector<std::string>   sv_points_{"jet_sv_eta",
                                                "jet_sv_phi"};
    const std::vector<std::string> sv_features_{"jet_sv_pt_log",
                                                "jet_sv_mass",
                                                "jet_sv_eta",
                                                "jet_sv_phi",
                                                "jet_sv_ntrack",
                                                "jet_sv_chi2",
                                                "jet_sv_dxy",
                                                "jet_sv_dxysig",
                                                "jet_sv_d3d",
                                                "jet_sv_d3dsig"};
    const std::vector<std::string>     sv_mask_{"sv_mask"};
};

int center_norm_pad(const std::vector<float> &input,
                    float center,
                    float norm_factor,
                    unsigned min_length,
                    unsigned max_length,
                    std::vector<float> &datavec,
                    int startval,
                    float pad_value,
                    float replace_inf_value,
                    float min,
                    float max) {
  assert(min <= pad_value && pad_value <= max);
  assert(min_length <= max_length);

  unsigned target_length = std::clamp((unsigned)input.size(), min_length, max_length);
  for (unsigned i = 0; i < target_length; ++i) {
    if (i < input.size()) {
      datavec[i + startval] = std::clamp((btagbtvdeep::catch_infs(input[i], replace_inf_value) - center) * norm_factor, min, max);
    } else {
      datavec[i + startval] = pad_value;
    }
  }
  return target_length;
}

MiniAODJetIdEmbedder::MiniAODJetIdEmbedder(const edm::ParameterSet& pset, const ONNXRuntime *cache) {
  srcToken_ = consumes<edm::View<pat::Jet>>(pset.getParameter<edm::InputTag>("src"));
  genParticlesToken_ = consumes<edm::View<reco::GenParticle>>(pset.getParameter<edm::InputTag>("genParticles"));
  vtxToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertices"));
  svToken_ = consumes<edm::View<reco::VertexCompositePtrCandidate>>(pset.getParameter<edm::InputTag>("secondaryVertices"));
  auto json_path = pset.getParameter<edm::FileInPath>("preprocess_json");

  input_names_.clear(); input_shapes_.clear();
  std::ifstream ifs(edm::FileInPath(json_path).fullPath());
  nlohmann::json js = nlohmann::json::parse(ifs);
  js.at("input_names").get_to(input_names_);
  for (const auto &group_name : input_names_) {
    const auto &group_pset = js.at(group_name);
    auto &prep_params = prep_info_map_[group_name];
    group_pset.at("var_names").get_to(prep_params.var_names);
    if (group_pset.contains("var_length")) {
      prep_params.min_length = group_pset.at("var_length");
      prep_params.max_length = prep_params.min_length;
      input_shapes_.push_back({1, (int64_t)prep_params.var_names.size(), prep_params.min_length});
    } else {
      prep_params.min_length = group_pset.at("min_length");
      prep_params.max_length = group_pset.at("max_length");
      input_shapes_.push_back({1, (int64_t)prep_params.var_names.size(), -1});
    }
    const auto &var_info_pset = group_pset.at("var_infos");
    for (const auto &var_name : prep_params.var_names) {
      const auto &var_pset = var_info_pset.at(var_name);
      double median = var_pset.at("median");
      double norm_factor = var_pset.at("norm_factor");
      double replace_inf_value = var_pset.at("replace_inf_value");
      double lower_bound = var_pset.at("lower_bound");
      double upper_bound = var_pset.at("upper_bound");
      double pad = var_pset.contains("pad") ? double(var_pset.at("pad")) : 0;
      prep_params.var_info_map[var_name] =
          PreprocessParams::VarInfo(median, norm_factor, replace_inf_value, lower_bound, upper_bound, pad);
    }
    if (&data_ != nullptr) {
      const auto &len = input_sizes_.emplace_back(prep_params.max_length * prep_params.var_names.size());
      data_.emplace_back(len, 0);
    }
  }

  produces<pat::JetCollection>();
}

std::unique_ptr<ONNXRuntime> MiniAODJetIdEmbedder::initializeGlobalCache(const edm::ParameterSet &iConfig) {
  return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void MiniAODJetIdEmbedder::globalEndJob(const ONNXRuntime *cache) {}

void MiniAODJetIdEmbedder::produce(edm::Event& evt, const edm::EventSetup& es) {
  std::unique_ptr<pat::JetCollection> output(new pat::JetCollection);

  edm::Handle<edm::View<pat::Jet> > input;
  evt.getByToken(srcToken_, input);

  output->reserve(input->size());

  edm::Handle<reco::VertexCollection> vtxHandle;
  evt.getByToken(vtxToken_, vtxHandle);
  pv_ = &vtxHandle->at(0);

  edm::Handle<edm::View<reco::VertexCompositePtrCandidate> > svHandle;
  evt.getByToken(svToken_, svHandle);

  edm::ESHandle<TransientTrackBuilder> track_builder_;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);
  TrackInfoBuilder trackinfo(track_builder_);

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

    // Add ParticleNet score
    float PNetScore = -1.;
    if(jet.pt() > 10. && abs(jet.eta()) < 2.4) {

      // create jet features
      DeepBoostedJetFeatures features;
      for (const auto &name : pf_points_)
        features.add(name);
      for (const auto &name : pf_features_)
        features.add(name);
      for (const auto &name : pf_mask_)
        features.add(name);
      for (const auto &name : sv_points_)
        features.add(name);
      for (const auto &name : sv_features_)
        features.add(name);
      for (const auto &name : sv_mask_)
        features.add(name);

      // build trackinfo
      math::XYZVector jet_dir = jet.momentum().Unit();
      GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());

      // reserve space
      int nConstituents = jet.numberOfSourceCandidatePtrs();
      for (const auto &name : pf_points_)
        features.reserve(name, nConstituents);
      for (const auto &name : pf_features_)
        features.reserve(name, nConstituents);
      for (const auto &name : pf_mask_)
        features.reserve(name, nConstituents);

      // fill pfcand features
      float etasign = jet.eta() > 0 ? 1 : -1;
      for(int k = 0; k < nConstituents; k++){
        reco::CandidatePtr pfcand = jet.sourceCandidatePtr(k);
        const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*pfcand));
        features.fill("pfcand_mask", 1);
        features.fill("jet_pfcand_pt_log", btagbtvdeep::catch_infs(std::log(packed_cand->pt()), -99));
        features.fill("jet_pfcand_e_log", btagbtvdeep::catch_infs(std::log(packed_cand->energy()), -99));
        features.fill("jet_pfcand_phirel", reco::deltaPhi(*packed_cand, jet));
        features.fill("jet_pfcand_etarel", etasign * (packed_cand->eta() - jet.eta()));
        features.fill("jet_pfcand_abseta", std::abs(packed_cand->eta()));
        features.fill("jet_pfcand_puppiw", packed_cand->puppiWeight());
        features.fill("jet_pfcand_charge", packed_cand->charge());
        features.fill("jet_pfcand_isEl", std::abs(packed_cand->pdgId()) == 11);
        features.fill("jet_pfcand_isMu", std::abs(packed_cand->pdgId()) == 13);
        features.fill("jet_pfcand_VTX_ass", packed_cand->pvAssociationQuality());
        features.fill("jet_pfcand_lostInnerHits", packed_cand->lostInnerHits());
        // impact parameters
        features.fill("jet_pfcand_dz", btagbtvdeep::catch_infs(packed_cand->dz()));
        features.fill("jet_pfcand_dzsig", packed_cand->bestTrack() ? btagbtvdeep::catch_infs(packed_cand->dz() / packed_cand->dzError()) : 0);
        features.fill("jet_pfcand_dxy", btagbtvdeep::catch_infs(packed_cand->dxy()));
        features.fill("jet_pfcand_dxysig", packed_cand->bestTrack() ? btagbtvdeep::catch_infs(packed_cand->dxy() / packed_cand->dxyError()) : 0);
        // track info
        if (packed_cand->bestTrack()) {
          const auto* trk = packed_cand->bestTrack();
          features.fill("jet_pfcand_normchi2", btagbtvdeep::catch_infs(trk->normalizedChi2()));
          features.fill("jet_pfcand_quality", trk->qualityMask());
  
          trackinfo.buildTrackInfo(&(*pfcand), jet_dir, jet_ref_track_dir, *pv_);
          features.fill("jet_pfcand_btagEtaRel", trackinfo.getTrackEtaRel());
          features.fill("jet_pfcand_btagPtRatio", trackinfo.getTrackPtRatio());
          features.fill("jet_pfcand_btagPParRatio", trackinfo.getTrackPParRatio());
          features.fill("jet_pfcand_btagSip3dVal", trackinfo.getTrackSip3dVal());
          features.fill("jet_pfcand_btagSip3dSig", trackinfo.getTrackSip3dSig());
          features.fill("jet_pfcand_btagJetDistVal", trackinfo.getTrackJetDistVal());
        }
        else {
          features.fill("jet_pfcand_normchi2", 999);
          features.fill("jet_pfcand_quality", 0);
          features.fill("jet_pfcand_btagEtaRel", 0);
          features.fill("jet_pfcand_btagPtRatio", 0);
          features.fill("jet_pfcand_btagPParRatio", 0);
          features.fill("jet_pfcand_btagSip3dVal", 0);
          features.fill("jet_pfcand_btagSip3dSig", 0);
          features.fill("jet_pfcand_btagJetDistVal", 0);
        }
      }

      // find SVs associated with the jet
      std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
      for(const auto &sv : *svHandle){
        if(reco::deltaR(jet, sv) < 0.4){
          jetSVs.push_back(&sv);
        }
      }

      // sort by dxy significance
      std::sort(jetSVs.begin(), jetSVs.end(), [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) { return btagbtvdeep::sv_vertex_comparator(*sva, *svb, *pv_); });

      // reserve space
      int nsv = jetSVs.size();
      for (const auto &name : sv_points_)
        features.reserve(name, nsv);
      for (const auto &name : sv_features_)
        features.reserve(name, nsv);
      for (const auto &name : sv_mask_)
        features.reserve(name, nsv);

      // fill associated secondary vertices information
      for (const auto *jetsv : jetSVs) {
        features.fill("sv_mask", 1);
        features.fill("jet_sv_pt_log", std::log(jetsv->pt()));
        features.fill("jet_sv_eta", jetsv->eta());
        features.fill("jet_sv_phi", jetsv->phi());
        features.fill("jet_sv_mass", jetsv->mass());
        features.fill("jet_sv_chi2", jetsv->vertexChi2());
        const auto &dxy = btagbtvdeep::vertexDxy(*jetsv, *pv_);
        features.fill("jet_sv_dxy", dxy.value());
        features.fill("jet_sv_dxysig", dxy.significance());
        const auto &d3d = btagbtvdeep::vertexD3d(*jetsv, *pv_);
        features.fill("jet_sv_d3d", d3d.value());
        features.fill("jet_sv_d3dsig", d3d.significance());
        features.fill("jet_sv_ntrack", jetsv->numberOfDaughters());
      }

      features.check_consistency(pf_points_);
      features.check_consistency(pf_features_);
      features.check_consistency(pf_mask_);
      features.check_consistency(sv_points_);
      features.check_consistency(sv_features_);
      features.check_consistency(sv_mask_);
 
      // make input
      // following: https://github.com/cms-sw/cmssw/blob/4e6da5e/RecoBTag/ONNXRuntime/plugins/BoostedJetONNXJetTagsProducer.cc#L196
      for (unsigned igroup = 0; igroup < input_names_.size(); ++igroup) {
        const auto &group_name = input_names_[igroup];
        const auto &prep_params = prep_info_map_.at(group_name);
        auto &group_values = data_[igroup];
        group_values.resize(input_sizes_[igroup]);
        // first reset group_values to 0
        std::fill(group_values.begin(), group_values.end(), 0);
        unsigned curr_pos = 0;
        // transform/pad
        for (unsigned i = 0; i < prep_params.var_names.size(); ++i) {
          const auto &varname = prep_params.var_names[i];
          const auto &raw_value = features.get(varname);
          const auto &info = prep_params.info(varname);
          int insize = center_norm_pad(raw_value,
                                       info.center,
                                       info.norm_factor,
                                       prep_params.min_length,
                                       prep_params.max_length,
                                       group_values,
                                       curr_pos,
                                       info.pad,
                                       info.replace_inf_value,
                                       info.lower_bound,
                                       info.upper_bound);
          curr_pos += insize;
          if (i == 0 && (!input_shapes_.empty())) {
            input_shapes_[igroup][2] = insize;
          }
        }
        group_values.resize(curr_pos);
      }

      // get output
      std::vector<float> outputs = globalCache()->run(input_names_, data_, input_shapes_)[0];
      PNetScore = outputs[0];
    }

    jet.addUserFloat("PNetScore", float(PNetScore));
    output->push_back(jet);
  }

  evt.put(std::move(output));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MiniAODJetIdEmbedder);
