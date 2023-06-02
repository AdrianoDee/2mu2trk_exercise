#ifndef __OniaPseudoTrackTrackProducer_h_
#define __OniaPseudoTrackTrackProducer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include <vector>
class OniaPseudoTrackTrackProducer : public edm::stream::EDProducer<>  {

 public:
  explicit OniaPseudoTrackTrackProducer(const edm::ParameterSet& ps);

 private:

  virtual void produce(edm::Event& event, const edm::EventSetup& esetup) override;
  virtual void endJob();

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> OniaCollection_;
  const edm::EDGetTokenT<reco::BeamSpot> BeamSpot_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> TrakCollection_;
  std::vector<double> OniaMassCuts_;
  std::vector<double> CandidateMassCuts_;
  const double Track1Mass_;
  const double Track2Mass_;
  const double ConstraintMass_;

	edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> builderToken_;

  const bool IsTheSame(const reco::Track& tk, const pat::Muon& mu);
  const bool IsTheSame(const reco::Track& tk1, const reco::Track& tk2);

  const RefCountedKinematicParticle FitPhoton(const reco::Track &tk0, const reco::Track &tk1, edm::ESHandle<TransientTrackBuilder> &theB);

  int candidates;
  int nevents;
  std::vector<reco::Track> kaons;
  std::vector<int> kaonsPV;
  std::vector<int> refVtx;
  std::vector<int> pvAssocQ;
  std::vector<float> dzAssocPV;

  GreaterByPt<reco::Track> PtComparator;

  template<typename T>
  struct GreaterByVProb {
         typedef T first_argument_type;
         typedef T second_argument_type;
         bool operator()( const T & t1, const T & t2 ) const { return t1.userFloat("vProb") > t2.userFloat("vProb"); }
  };
};

#endif
