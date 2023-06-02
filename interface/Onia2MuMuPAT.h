#ifndef _Onia2MuMuPAT_h
#define _Onia2MuMuPAT_h

// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

template<typename T>
struct GreaterByVProb {
  bool operator()( const T & t1, const T & t2 ) const {
    return t1.userFloat("vProb") > t2.userFloat("vProb");
  }
};


//
// class decleration
//

class Onia2MuMuPAT : public edm::stream::EDProducer<> {
 public:
  explicit Onia2MuMuPAT(const edm::ParameterSet&);
  ~Onia2MuMuPAT() override;

 private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  bool isAbHadron(int pdgID);
  bool isAMixedbHadron(int pdgID, int momPdgID);
  std::pair<int, float> findJpsiMCInfo(reco::GenParticleRef genJpsi);

  // ----------member data ---------------------------
 private:

  edm::EDGetTokenT<edm::View<pat::Muon>> muons_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  edm::EDGetTokenT<reco::TrackCollection> revtxtrks_;
  edm::EDGetTokenT<reco::BeamSpot> revtxbs_;
  StringCutObjectSelector<pat::Muon> higherPuritySelection_;
  StringCutObjectSelector<pat::Muon> lowerPuritySelection_;
  StringCutObjectSelector<reco::Candidate, true> dimuonSelection_;
  bool addCommonVertex_, addMuonlessPrimaryVertex_;
  bool resolveAmbiguity_;
  bool addMCTruth_;

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBuilderToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> theFieldToken_;

  GreaterByVProb<pat::CompositeCandidate> vPComparator_;

  InvariantMassFromVertex massCalculator;

  int nMuMu = 0;
  int nMuMuM = 0;
  int nEvents = 0;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//

#endif
