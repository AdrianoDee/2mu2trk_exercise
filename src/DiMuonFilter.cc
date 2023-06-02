// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormat includes
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"


class DiMuonFilter : public edm::stream::EDProducer<>  {
 public:
  explicit DiMuonFilter(const edm::ParameterSet&);
  ~DiMuonFilter() override {};
  UInt_t isTriggerMatched(const pat::CompositeCandidate *);
 private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT<std::vector<pat::CompositeCandidate>> theOnias_;
  StringCutObjectSelector<reco::Candidate, true> SingleMuonSelection_;
  StringCutObjectSelector<reco::Candidate, true> DiMuonSelection_;
  bool do_trigger_match_;
  std::vector<std::string> HLTFilters_;
};

DiMuonFilter::DiMuonFilter(const edm::ParameterSet& iConfig):
  theOnias_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("OniaTag"))),
  SingleMuonSelection_(iConfig.existsAs<std::string>("singlemuonSelection") ? iConfig.getParameter<std::string>("singlemuonSelection") : ""),
  DiMuonSelection_(iConfig.existsAs<std::string>("dimuonSelection") ? iConfig.getParameter<std::string>("dimuonSelection") : ""),
  do_trigger_match_(iConfig.getParameter<bool>("do_trigger_match")),
  HLTFilters_(iConfig.getParameter<std::vector<std::string>>("HLTFilters"))
{
  produces<pat::CompositeCandidateCollection>();
}

UInt_t DiMuonFilter::isTriggerMatched(const pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  UInt_t matched = 0;  // if no list is given, is not matched

// if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
     const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
     const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
     if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
  }
  return matched;
}

// ------------ method called to produce the data  ------------
void DiMuonFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::unique_ptr<pat::CompositeCandidateCollection> mumuOutput(new pat::CompositeCandidateCollection);
  edm::Handle<pat::CompositeCandidateCollection> onias_;
  iEvent.getByToken(theOnias_, onias_);
  if (onias_.isValid() && !onias_->empty()) {
    const pat::CompositeCandidate *ionia = nullptr;
    for (size_t ii = 0, nn=onias_->size(); ii < nn; ii++ ) {
       ionia = &(onias_->at(ii));
       if (ionia && DiMuonSelection_(*ionia) &&
           SingleMuonSelection_(*ionia->daughter("muon1")) &&
           SingleMuonSelection_(*ionia->daughter("muon2")) &&
           ( !do_trigger_match_ || isTriggerMatched(ionia))
          ) mumuOutput->push_back(*ionia);
    }
  }
  iEvent.put(std::move(mumuOutput));
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonFilter);
