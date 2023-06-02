#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>
#include <sstream>

class OniaRecoTrackTrackRootupler : public edm::one::EDAnalyzer<> {
   public:
      explicit OniaRecoTrackTrackRootupler(const edm::ParameterSet&);
      ~OniaRecoTrackTrackRootupler() override {};
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void analyze(const edm::Event&, const edm::EventSetup&) override;

  std::string file_name;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> TheCandidateLabel;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> TheUpsLabel;
  const edm::EDGetTokenT<reco::VertexCollection>            ThePrimaryVertexLabel;
  const edm::EDGetTokenT<edm::TriggerResults>               TheTriggerResultLabel;
  const edm::EDGetTokenT<reco::GenParticleCollection>       TheGenParticleLabel;
  const double Track1Mass_;
  const double Track2Mass_;
  const double DimuonMass_;
  const int  candidate_pdgid_, onia_pdgid_, ditrack_pdgid_, track1_pdgid_, track2_pdgid_;
  const bool isMC_,OnlyBest_;

  UInt_t run, event, nCandPerEvent, numPrimaryVertices, trigger;

  TLorentzVector candidate_p4;
  TLorentzVector dimuon_p4;
  TLorentzVector track1_p4;
  TLorentzVector track2_p4;
  TLorentzVector muonp_p4;
  TLorentzVector muonn_p4;
  TLorentzVector ditrack_p4;

  Int_t    candidate_charge, track1_nvsh, track1_nvph, track2_nvsh, track2_nvph, gen_candidate_charge;
  Double_t candidate_vMass, candidate_vProb,  candidate_vChi2, candidate_cosAlpha, candidate_ctauPV, candidate_ctauErrPV;
  Double_t candidate_cosAlpha3D, candidate_lxy, candidate_lxyErr, candidate_lxyz, candidate_lxyzErr;

  Double_t thePrimaryV_X, thePrimaryV_Y, thePrimaryV_Z, TheDecayVertex_X, TheDecayVertex_Y, TheDecayVertex_Z, thePrimaryV_2D_position, thePrimaryV_3D_position, TheDecayVertex_2D_position, TheDecayVertex_3D_position, TheVertexDistance_2D, TheVertexDistance_3D;

  Double_t track1_d0, track1_d0Err, track1_dz, track1_dxy, track1_dzErr, track1_dxyErr, track1_charge, track1_dzAssocPV;
  Double_t track2_d0, track2_d0Err, track2_dz, track2_dxy, track2_dzErr, track2_dxyErr, track2_charge, track2_dzAssocPV;
  Double_t dimuon_vProb, dimuon_vChi2, dimuon_DCA, dimuon_ctauPV, dimuon_ctauErrPV, dimuon_cosAlpha, dimuon_nSigma;
  Double_t track2_dRdimuon, track1_dRdimuon, ditrack_dRdimuon;

  Int_t track1_PV, track2_PV, track1_refVtx, track2_refVtx, track1_pvAssocQ, track2_pvAssocQ;

  Int_t dimuon_vertexWeight, iPVwithmuons;

  Double_t invm1Skk;
  Double_t invm2Skk;

  TLorentzVector gen_candidate_p4;
  Int_t          gen_candidate_pdgId;
  TLorentzVector gen_dimuon_p4;
  Int_t          gen_onia_pdgId;
  TLorentzVector gen_ditrack_p4;
  TLorentzVector gen_track1_p4;
  Int_t          gen_track1_pdgid;
  TLorentzVector gen_track2_p4;
  Int_t          gen_track2_pdgid;
  TLorentzVector gen_muonp_p4;
  TLorentzVector gen_muonn_p4;

  TLorentzVector track1;
  TLorentzVector track2;

  TLorentzVector ups_p4, muonP_p4, muonN_p4;
  Double_t ups_vMass, ups_vertexWeight, ups_vProb, ups_vChi2, ups_DCA, ups_ctauPV, ups_ctauErrPV, ups_cosAlpha;
  Double_t ups_lxyPV, ups_lxyErrPV, ups_ctauBS, ups_ctauErrBS, ups_lxyBS, ups_lxyErrBS;
  Double_t mu1_pt, mu1_ptErr, mu1_d0, mu1_d0Err, mu1_dz, mu1_dzErr, mu1_dxy, mu1_dxyErr, mu2_pt, mu2_ptErr, mu2_d0, mu2_d0Err, mu2_dz, mu2_dzErr, mu2_dxy, mu2_dxyErr;
  Int_t mu1_nvsh, mu1_nvph, mu2_nvsh, mu2_nvph, iPVwithmuons_ups, mu1_charge, mu2_charge;

  TTree* TheTree;
  TTree* UpsTree;
  bool is_dimuon_;
};

static const Double_t ups1SMass =  9.46030;
static const Double_t ups2SMass = 10.02326;

OniaRecoTrackTrackRootupler::OniaRecoTrackTrackRootupler(const edm::ParameterSet& iConfig):
TheCandidateLabel(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("TheCandidates"))),
TheUpsLabel(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("TheUps"))),
ThePrimaryVertexLabel(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag > ("PrimaryVertices"))),
TheTriggerResultLabel(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag > ("TriggerResults"))),
TheGenParticleLabel(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
Track1Mass_(iConfig.getParameter<double>("Track1Mass")),
Track2Mass_(iConfig.getParameter<double>("Track2Mass")),
DimuonMass_(iConfig.getParameter<double>("DimuonMass")),
candidate_pdgid_(iConfig.getParameter<uint32_t>("candidate_pdgid")),
onia_pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
ditrack_pdgid_(iConfig.getParameter<uint32_t>("ditrack_pdgid")),
track1_pdgid_(iConfig.getParameter<int32_t>("track1_pdgid")),
track2_pdgid_(iConfig.getParameter<int32_t>("track2_pdgid")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest"))
{
        is_dimuon_ = (onia_pdgid_ == 333 || onia_pdgid_ == 443 || onia_pdgid_ == 100443 || onia_pdgid_ == 553 || onia_pdgid_ == 100553 || onia_pdgid_ == 200553);
	edm::Service<TFileService> fs;
        TheTree = fs->make<TTree>("CandidateTree","CandidateTree");

        TheTree->Branch("run",                &run,                "run/I");
        TheTree->Branch("event",              &event,              "event/I");
        TheTree->Branch("nCandPerEvent", &nCandPerEvent, "nCandPerEvent/I");
        TheTree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        TheTree->Branch("trigger",            &trigger,            "trigger/I");

        TheTree->Branch("candidate_p4","TLorentzVector", &candidate_p4);
        TheTree->Branch("track1_p4",   "TLorentzVector", &track1_p4);
        TheTree->Branch("track2_p4",   "TLorentzVector", &track2_p4);
        TheTree->Branch("ditrack_p4",   "TLorentzVector", &ditrack_p4);
        TheTree->Branch("dimuon_p4",   "TLorentzVector", &dimuon_p4);

        TheTree->Branch("muonp_p4",    "TLorentzVector", &muonp_p4);
        TheTree->Branch("muonn_p4",    "TLorentzVector", &muonn_p4);

        TheTree->Branch("invm1Skk",      &invm1Skk,          "invm1Skk/D");
        TheTree->Branch("invm2Skk",      &invm2Skk,          "invm2Skk/D");

        TheTree->Branch("iPVwithmuons",        &iPVwithmuons,        "iPVwithmuons/I");

        TheTree->Branch("dimuon_vertexWeight",        &dimuon_vertexWeight,        "dimuon_vertexWeight/I");
        TheTree->Branch("dimuon_vProb",        &dimuon_vProb,        "dimuon_vProb/D");
        TheTree->Branch("dimuon_vNChi2",       &dimuon_vChi2,        "dimuon_vNChi2/D");
        TheTree->Branch("dimuon_DCA",          &dimuon_DCA,          "dimuon_DCA/D");
        TheTree->Branch("dimuon_ctauPV",       &dimuon_ctauPV,       "dimuon_ctauPV/D");
        TheTree->Branch("dimuon_ctauErrPV",    &dimuon_ctauErrPV,    "dimuon_ctauErrPV/D");
        TheTree->Branch("dimuon_cosAlpha",     &dimuon_cosAlpha,     "dimuon_cosAlpha/D");
        TheTree->Branch("dimuon_nSigma",       &dimuon_nSigma,       "dimuon_nSigma/D");

        TheTree->Branch("candidate_vMass",      &candidate_vMass,        "candidate_vMass/D");
        TheTree->Branch("candidate_vProb",      &candidate_vProb,        "candidate_vProb/D");
        TheTree->Branch("candidate_vChi2",      &candidate_vChi2,        "candidate_vChi2/D");
        TheTree->Branch("candidate_cosAlpha",   &candidate_cosAlpha,     "candidate_cosAlpha/D");
        TheTree->Branch("candidate_ctauPV",     &candidate_ctauPV,       "candidate_ctauPV/D");
        TheTree->Branch("candidate_ctauErrPV",  &candidate_ctauErrPV,    "candidate_ctauErrPV/D");
        TheTree->Branch("candidate_charge",     &candidate_charge,       "candidate_charge/I");
        TheTree->Branch("candidate_lxy",        &candidate_lxy,          "candidate_lxy/D");
        TheTree->Branch("candidate_lxyErr",     &candidate_lxyErr,       "candidate_lxyErr/D");
        TheTree->Branch("candidate_lxyz",       &candidate_lxyz,         "candidate_lxyz/D");
        TheTree->Branch("candidate_lxyzErr",    &candidate_lxyzErr,      "candidate_lxyzErr/D");

        TheTree->Branch("thePrimaryV_X",      &thePrimaryV_X,        "thePrimaryV_X/D");
        TheTree->Branch("thePrimaryV_Y",      &thePrimaryV_Y,        "thePrimaryV_Y/D");
        TheTree->Branch("thePrimaryV_Z",      &thePrimaryV_Z,        "thePrimaryV_Z/D");
        TheTree->Branch("TheDecayVertex_X",      &TheDecayVertex_X,        "TheDecayVertex_X/D");
        TheTree->Branch("TheDecayVertex_Y",      &TheDecayVertex_Y,        "TheDecayVertex_Y/D");
        TheTree->Branch("TheDecayVertex_Z",      &TheDecayVertex_Z,        "TheDecayVertex_Z/D");
        TheTree->Branch("thePrimaryV_2D_position",      &thePrimaryV_2D_position,        "thePrimaryV_2D_position/D");
        TheTree->Branch("thePrimaryV_3D_position",      &thePrimaryV_3D_position,        "thePrimaryV_3D_position/D");
        TheTree->Branch("TheDecayVertex_2D_position",      &TheDecayVertex_2D_position,        "TheDecayVertex_2D_position/D");
        TheTree->Branch("TheDecayVertex_3D_position",      &TheDecayVertex_3D_position,        "TheDecayVertex_3D_position/D");
        TheTree->Branch("TheVertexDistance_2D",      &TheVertexDistance_2D,        "TheVertexDistance_2D/D");
        TheTree->Branch("TheVertexDistance_3D",      &TheVertexDistance_3D,        "TheVertexDistance_3D/D");

        TheTree->Branch("track1_d0",    &track1_d0,    "track1_d0/D");
        TheTree->Branch("track1_d0Err", &track1_d0Err, "track1_d0Err/D");
        TheTree->Branch("track1_dz",    &track1_dz,    "track1_dz/D");
        TheTree->Branch("track1_dzErr",    &track1_dzErr,    "track1_dzErr/D");
        TheTree->Branch("track1_dxy",   &track1_dxy,   "track1_dxy/D");
        TheTree->Branch("track1_dxyErr",   &track1_dxyErr,   "track1_dxyErr/D");
        TheTree->Branch("track1_nvsh",  &track1_nvsh,  "track1_nvsh/I");
        TheTree->Branch("track1_nvph",  &track1_nvph,  "track1_nvph/I");
        TheTree->Branch("track1_dRdimuon",  &track1_dRdimuon,  "track1_dRdimuon/D");
        TheTree->Branch("track1_charge",  &track1_charge,  "track1_charge/D");
        TheTree->Branch("track1_PV",  &track1_PV,  "track1_PV/I");
        TheTree->Branch("track1_refVtx",  &track1_refVtx,  "track1_refVtx/I");
        TheTree->Branch("track1_pvAssocQ",  &track1_pvAssocQ,  "track1_pvAssocQ/I");
        TheTree->Branch("track1_dzAssocPV",    &track1_dzAssocPV,    "track1_dzAssocPV/D");

        TheTree->Branch("track2_d0",    &track2_d0,    "track2_d0/D");
        TheTree->Branch("track2_d0Err", &track2_d0Err, "track2_d0Err/D");
        TheTree->Branch("track2_dz",    &track2_dz,    "track2_dz/D");
        TheTree->Branch("track2_dzErr",    &track2_dzErr,    "track2_dzErr/D");
        TheTree->Branch("track2_dxy",   &track2_dxy,   "track2_dxy/D");
        TheTree->Branch("track2_dxyErr",   &track2_dxyErr,   "track2_dxyErr/D");
        TheTree->Branch("track2_nvsh",  &track2_nvsh,  "track2_nvsh/I");
        TheTree->Branch("track2_nvph",  &track2_nvph,  "track2_nvph/I");
        TheTree->Branch("track2_dRdimuon",  &track2_dRdimuon,  "track2_dRdimuon/D");
        TheTree->Branch("track2_charge",  &track2_charge,  "track2_charge/D");
        TheTree->Branch("track2_PV",  &track2_PV,  "track2_PV/I");
        TheTree->Branch("track2_refVtx",  &track2_refVtx,  "track2_refVtx/I");
        TheTree->Branch("track2_pvAssocQ",  &track2_pvAssocQ,  "track2_pvAssocQ/I");
        TheTree->Branch("track2_dzAssocPV",    &track2_dzAssocPV,    "track2_dzAssocPV/D");

        TheTree->Branch("ditrack_dRdimuon",  &ditrack_dRdimuon,  "ditrack_dRdimuon/D");

	if(isMC_)
	  {
      TheTree->Branch("gen_candidate_pdgId", &gen_candidate_pdgId, "gen_candidate_pdgId/I");
      TheTree->Branch("gen_onia_pdgId",      &gen_onia_pdgId,      "gen_onia_pdgId/I");
	  TheTree->Branch("gen_candidate_p4","TLorentzVector", &gen_candidate_p4);
      TheTree->Branch("gen_dimuon_p4",   "TLorentzVector", &gen_dimuon_p4);
      if (ditrack_pdgid_) TheTree->Branch("gen_ditrack_p4",   "TLorentzVector", &gen_ditrack_p4);
	    TheTree->Branch("gen_track1_p4",   "TLorentzVector", &gen_track1_p4);
      TheTree->Branch("gen_track1_pdgid", &gen_track1_pdgid, "gen_track1_pdgid/I");
      TheTree->Branch("gen_track2_p4",   "TLorentzVector", &gen_track2_p4);
      TheTree->Branch("gen_track2_pdgid", &gen_track2_pdgid, "gen_track2_pdgid/I");
      TheTree->Branch("gen_muonp_p4",    "TLorentzVector", &gen_muonp_p4);
      TheTree->Branch("gen_muonn_p4",    "TLorentzVector", &gen_muonn_p4);
      TheTree->Branch("gen_candidate_charge",     &gen_candidate_charge,       "gen_candidate_charge/I");
	  }

        UpsTree = fs->make<TTree>("UpsTree","UpsTree");

        UpsTree->Branch("run",                &run,                "run/I");
        UpsTree->Branch("event",              &event,              "event/I");
        UpsTree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        UpsTree->Branch("trigger",            &trigger,            "trigger/I");
        UpsTree->Branch("ups_p4",   "TLorentzVector", &ups_p4);
        UpsTree->Branch("muonP_p4",    "TLorentzVector", &muonP_p4);
        UpsTree->Branch("muonN_p4",    "TLorentzVector", &muonN_p4);

        UpsTree->Branch("iPVwithmuons_ups",        &iPVwithmuons_ups,        "iPVwithmuons_ups/I");

        UpsTree->Branch("ups_vertexWeight",        &ups_vertexWeight,        "ups_vertexWeight/D");
        UpsTree->Branch("ups_vProb",        &ups_vProb,        "ups_vProb/D");
        UpsTree->Branch("ups_vMass",        &ups_vMass,        "ups_vMass/D");
        UpsTree->Branch("ups_vNChi2",       &ups_vChi2,        "ups_vNChi2/D");
        UpsTree->Branch("ups_DCA",          &ups_DCA,          "ups_DCA/D");
        UpsTree->Branch("ups_ctauPV",       &ups_ctauPV,       "ups_ctauPV/D");
        UpsTree->Branch("ups_ctauErrPV",    &ups_ctauErrPV,    "ups_ctauErrPV/D");
        UpsTree->Branch("ups_lxyPV",        &ups_lxyPV,          "ups_lxyPV/D");
        UpsTree->Branch("ups_lxyErrPV",     &ups_lxyErrPV,       "ups_lxyErrPV/D");
        UpsTree->Branch("ups_cosAlpha",     &ups_cosAlpha,     "ups_cosAlpha/D");
        UpsTree->Branch("ups_ctauBS",       &ups_ctauBS,       "ups_ctauBS/D");
        UpsTree->Branch("ups_ctauErrBS",    &ups_ctauErrBS,    "ups_ctauErrBS/D");
        UpsTree->Branch("ups_lxyBS",        &ups_lxyBS,          "ups_lxyBS/D");
        UpsTree->Branch("ups_lxyErrBS",     &ups_lxyErrBS,       "ups_lxyErrBS/D");

        UpsTree->Branch("mu1_pt",    &mu1_pt,    "mu1_pt/D");
        UpsTree->Branch("mu1_ptErr",    &mu1_ptErr,    "mu1_ptErr/D");
        UpsTree->Branch("mu1_d0",    &mu1_d0,    "mu1_d0/D");
        UpsTree->Branch("mu1_d0Err", &mu1_d0Err, "mu1_d0Err/D");
        UpsTree->Branch("mu1_dz",    &mu1_dz,    "mu1_dz/D");
        UpsTree->Branch("mu1_dzErr",    &mu1_dzErr,    "mu1_dzErr/D");
        UpsTree->Branch("mu1_dxy",   &mu1_dxy,   "mu1_dxy/D");
        UpsTree->Branch("mu1_dxyErr",   &mu1_dxyErr,   "mu1_dxyErr/D");
        UpsTree->Branch("mu1_nvsh",  &mu1_nvsh,  "mu1_nvsh/I");
        UpsTree->Branch("mu1_nvph",  &mu1_nvph,  "mu1_nvph/I");
        UpsTree->Branch("mu1_charge",  &mu1_charge,  "mu1_charge/I");

        UpsTree->Branch("mu2_pt",    &mu2_pt,    "mu2_pt/D");
        UpsTree->Branch("mu2_ptErr",    &mu2_ptErr,    "mu2_ptErr/D");
        UpsTree->Branch("mu2_d0",    &mu2_d0,    "mu2_d0/D");
        UpsTree->Branch("mu2_d0Err", &mu2_d0Err, "mu2_d0Err/D");
        UpsTree->Branch("mu2_dz",    &mu2_dz,    "mu2_dz/D");
        UpsTree->Branch("mu2_dzErr",    &mu2_dzErr,    "mu2_dzErr/D");
        UpsTree->Branch("mu2_dxy",   &mu2_dxy,   "mu2_dxy/D");
        UpsTree->Branch("mu2_dxyErr",   &mu2_dxyErr,   "mu2_dxyErr/D");
        UpsTree->Branch("mu2_nvsh",  &mu2_nvsh,  "mu2_nvsh/I");
        UpsTree->Branch("mu2_nvph",  &mu2_nvph,  "mu2_nvph/I");
        UpsTree->Branch("mu2_charge",  &mu2_charge,  "mu2_charge/I");

}

void OniaRecoTrackTrackRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;

  edm::Handle < pat::CompositeCandidateCollection >TheCandidates;
  iEvent.getByToken(TheCandidateLabel, TheCandidates);

  edm::Handle < pat::CompositeCandidateCollection >TheUps;
  iEvent.getByToken(TheUpsLabel, TheUps);

  edm::Handle < reco::VertexCollection  >ThePrimaryVertices;
  iEvent.getByToken(ThePrimaryVertexLabel, ThePrimaryVertices);

  edm::Handle < edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken(TheTriggerResultLabel, triggerResults_handle);

  numPrimaryVertices = ThePrimaryVertices->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

  gen_ditrack_p4.SetPtEtaPhiM(0,0,0,0);

  if ( isMC_ ) {
    edm::Handle<reco::GenParticleCollection> GenParticles;
    iEvent.getByToken(TheGenParticleLabel, GenParticles);
    int foundit = 0;
    //int gen_track1_pdgid_ = 0;
    //int gen_track2_pdgid_ = 0;
    gen_candidate_pdgId = 0;
    if (GenParticles.isValid() ) {
      for ( reco::GenParticleCollection::const_iterator itParticle = GenParticles->begin(); itParticle != GenParticles->end(); ++itParticle ) {
        int pdgId = itParticle->pdgId();
        if ( abs(pdgId) ==  candidate_pdgid_ ) {
          //const reco::Candidate* gen_y2s = itParticle
          gen_candidate_charge = itParticle->charge();
          gen_candidate_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass());
          gen_candidate_pdgId = pdgId;
          foundit++;
          for (uint i = 0; i < itParticle->numberOfDaughters(); ++i) {
            const reco::Candidate* b = itParticle->daughter(i);
            int bpdgid = b->pdgId();
            if ( abs(bpdgid) == onia_pdgid_ && b->status() == 2 ) {
              gen_onia_pdgId = bpdgid;
              const reco::Candidate* d = nullptr;
              d = b;
              if (d) {
                   gen_dimuon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                   foundit++;
                   for (uint j = 0; j < d->numberOfDaughters(); ++j) {
                     const reco::Candidate* p = d->daughter(j);
                     if ( p->pdgId() == 13 && p->status() == 1 ) {
                        gen_muonp_p4.SetPtEtaPhiM(p->pt(),p->eta(),p->phi(),p->mass());
                        foundit++;
                     }
                     if ( p->pdgId() == -13 && p->status() == 1 ) {
                        gen_muonn_p4.SetPtEtaPhiM(p->pt(),p->eta(),p->phi(),p->mass());
                        foundit++;
                     }
                   }
              } // d
            }
            if (ditrack_pdgid_ && abs(bpdgid) == ditrack_pdgid_) {
              gen_ditrack_p4.SetPtEtaPhiM(b->pt(),b->eta(),b->phi(),b->mass());
              for (uint k = 0; k < b->numberOfDaughters(); ++k) {
                const reco::Candidate* p = b->daughter(k);
                if ( p->pdgId() == track1_pdgid_ && p->status() == 1 ) {
                   //std::cout<<" in dipion - pion 1 "<<std::endl;
                   //track1.SetPtEtaPhiM(p->pt(),p->eta(),p->phi(),p->mass());
                   //gen_track1_pdgid_ = track1_pdgid_;
                   gen_track1_p4.SetPtEtaPhiM(p->pt(),p->eta(),p->phi(),p->mass());
                   gen_track1_pdgid = track1_pdgid_;
                   foundit++;
                }
                if ( p->pdgId() == track2_pdgid_ && p->status() == 1 ) {
                   //std::cout<<" in dipion - pion 2 "<<std::endl;
                   //track2.SetPtEtaPhiM(p->pt(),p->eta(),p->phi(),p->mass());
                   //gen_track2_pdgid_ = track2_pdgid_;
                   gen_track2_p4.SetPtEtaPhiM(p->pt(),p->eta(),p->phi(),p->mass());
                   gen_track2_pdgid = track2_pdgid_;
                   foundit++;
                }
              }
            } else {
              if ( bpdgid == track1_pdgid_ && b->status() == 1 ) {
                //std::cout<<" pion 1 "<<std::endl;
                //track1.SetPtEtaPhiM(b->pt(),b->eta(),b->phi(),b->mass());
                //gen_track1_pdgid_ = track1_pdgid_;
                gen_track1_p4.SetPtEtaPhiM(b->pt(),b->eta(),b->phi(),b->mass());
                gen_track1_pdgid = track1_pdgid_;
                foundit++;
              }
              if ( bpdgid == track2_pdgid_ && b->status() == 1 ) {
                //std::cout<<" pion 2 "<<std::endl;
                //track2.SetPtEtaPhiM(b->pt(),b->eta(),b->phi(),b->mass());
                //gen_track2_pdgid_ = track2_pdgid_;
                gen_track2_p4.SetPtEtaPhiM(b->pt(),b->eta(),b->phi(),b->mass());
                gen_track2_pdgid = track2_pdgid_;
                foundit++;
              }
            }
          }
          if ( foundit == 6 ) break;
          else {
            foundit = 0;
            gen_candidate_pdgId = 0;
          }
        } // if ( abs(
      }   // for ( reco
    }
    if (!gen_candidate_pdgId) {std::cout << "OniaRecoTrackTrackRootupler: didn't find the given decay " << run << "," << event << std::endl;} /*else {std::cout << "I found it yupppiiiii " << run << "," << event << std::endl; }*/
  } // end if isMC

   trigger = 0;
   if (triggerResults_handle.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      unsigned int NTRIGGERS = 10;
      std::string TriggersToTest[NTRIGGERS] = {
        "HLT_Dimuon20_Jpsi_Barrel_Seagulls","HLT_Dimuon25_Jpsi",
        "HLT_Dimuon10_PsiPrime_Barrel_Seagulls","HLT_Dimuon18_PsiPrime",
        "HLT_Dimuon10_Upsilon_Barrel_Seagulls","HLT_Dimuon12_Upsilon_eta1p5",
        "HLT_Dimuon14_Phi_Barrel_Seagulls","HLT_Dimuon12_Upsilon_y1p4",
        "HLT_Dimuon8_Upsilon_Barrel","HLT_Dimuon13_Upsilon"
      };

      for (unsigned int i = 0; i < NTRIGGERS; i++) {
         for (int version = 1; version < 19; version++) {
            std::stringstream ss;
            ss << TriggersToTest[i] << "_v" << version;
            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
            if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
               trigger += (1<<i);
               break;
            }
         }
      }
    } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

// grabbing candidate information. Notice we are just keeping combinations with succesfull vertex fit
  if (TheCandidates.isValid() && !TheCandidates->empty()) {
    pat::CompositeCandidate TheCandidate_;
    nCandPerEvent = TheCandidates->size();
    for (unsigned int i=0; i< TheCandidates->size(); i++){
      TheCandidate_       = TheCandidates->at(i);
      candidate_vMass     = TheCandidate_.userFloat("vMass");
      candidate_vProb     = TheCandidate_.userFloat("vProb");
      candidate_vChi2     = TheCandidate_.userFloat("vChi2");
      candidate_cosAlpha  = TheCandidate_.userFloat("cosAlpha");
      candidate_ctauPV    = TheCandidate_.userFloat("ctauPV");
      candidate_ctauErrPV = TheCandidate_.userFloat("ctauErrPV");
      candidate_charge    = TheCandidate_.charge();
      candidate_lxy       = TheCandidate_.userFloat("lxy");
      candidate_lxyErr    = TheCandidate_.userFloat("lxyErr");
      candidate_lxyz      = TheCandidate_.userFloat("lxyz");
      candidate_lxyzErr   = TheCandidate_.userFloat("lxyzErr");

      thePrimaryV_X = TheCandidate_.userFloat("thePrimaryV_X");
      thePrimaryV_Y = TheCandidate_.userFloat("thePrimaryV_Y");
      thePrimaryV_Z = TheCandidate_.userFloat("thePrimaryV_Z");
      TheDecayVertex_X = TheCandidate_.userFloat("TheDecayVertex_X");
      TheDecayVertex_Y = TheCandidate_.userFloat("TheDecayVertex_Y");
      TheDecayVertex_Z = TheCandidate_.userFloat("TheDecayVertex_Z");
      thePrimaryV_2D_position = TheCandidate_.userFloat("thePrimaryV_2D_position");
      thePrimaryV_3D_position = TheCandidate_.userFloat("thePrimaryV_3D_position");
      TheDecayVertex_2D_position = TheCandidate_.userFloat("TheDecayVertex_2D_position");
      TheDecayVertex_3D_position = TheCandidate_.userFloat("TheDecayVertex_3D_position");
      TheVertexDistance_2D = TheCandidate_.userFloat("TheVertexDistance_2D");
      TheVertexDistance_3D = TheCandidate_.userFloat("TheVertexDistance_3D");

      ditrack_dRdimuon    = TheCandidate_.userFloat("ditrack_dRdimuon");
      track1_dRdimuon    = TheCandidate_.userFloat("track1_dRdimuon");
      track2_dRdimuon    = TheCandidate_.userFloat("track2_dRdimuon");
      track1_PV  = TheCandidate_.userInt("track1_PV");
      track2_PV  = TheCandidate_.userInt("track2_PV");
      track1_refVtx  = TheCandidate_.userInt("track1_refVtx");
      track2_refVtx  = TheCandidate_.userInt("track2_refVtx");
      track1_pvAssocQ  = TheCandidate_.userInt("track1_pvAssocQ");
      track2_pvAssocQ  = TheCandidate_.userInt("track2_pvAssocQ");

      track1_dzAssocPV = TheCandidate_.userFloat("track1_dzAssocPV");
      track2_dzAssocPV = TheCandidate_.userFloat("track2_dzAssocPV");

      const pat::CompositeCandidate *TheDimuon_ = nullptr;
      const reco::Vertex *ThePrimaryV_ = nullptr;

      TheDimuon_ = dynamic_cast <pat::CompositeCandidate *>(TheCandidate_.daughter("onia"));
      ThePrimaryV_ = TheDimuon_->userData<reco::Vertex>("PVwithmuons");

      const reco::Track *TheTrack1_ =  dynamic_cast <const reco::Track *>(TheCandidate_.userData<reco::Track>("track1"));
      const reco::Track *TheTrack2_ =  dynamic_cast <const reco::Track *>(TheCandidate_.userData<reco::Track>("track2"));

      candidate_p4.SetPtEtaPhiM(TheCandidate_.pt(),TheCandidate_.eta(),TheCandidate_.phi(),TheCandidate_.mass());
      dimuon_p4.SetPtEtaPhiM(TheDimuon_->pt(),TheDimuon_->eta(),TheDimuon_->phi(),TheDimuon_->mass());
      track1_p4.SetPtEtaPhiM(TheTrack1_->pt(),TheTrack1_->eta(),TheTrack1_->phi(),Track1Mass_);
      track2_p4.SetPtEtaPhiM(TheTrack2_->pt(),TheTrack2_->eta(),TheTrack2_->phi(),Track2Mass_);
      ditrack_p4 = track1_p4 + track2_p4;

      typedef math::XYZPoint Point;
      Point pv_(ThePrimaryV_->x(),ThePrimaryV_->y(),ThePrimaryV_->z());

      track1_d0      = TheTrack1_->d0();
      track1_d0Err   = TheTrack1_->d0Error();
      track1_dz      = TheTrack1_->dz(pv_);
      track1_dzErr   = TheTrack1_->dzError();
      track1_dxy     = TheTrack1_->dxy(pv_);
      track1_dxyErr  = TheTrack1_->dxyError();
      track1_nvsh    = TheTrack1_->hitPattern().numberOfValidStripHits();
      track1_nvph    = TheTrack1_->hitPattern().numberOfValidPixelHits();
      track1_charge  = TheTrack1_->charge();

      track2_d0      = TheTrack2_->d0();
      track2_d0Err   = TheTrack2_->d0Error();
      track2_dz      = TheTrack2_->dz(pv_);
      track2_dzErr   = TheTrack2_->dzError();
      track2_dxy     = TheTrack2_->dxy(pv_);
      track2_dxyErr  = TheTrack2_->dxyError();
      track2_nvsh    = TheTrack2_->hitPattern().numberOfValidStripHits();
      track2_nvph    = TheTrack2_->hitPattern().numberOfValidPixelHits();
      track2_charge  = TheTrack2_->charge();

      // iPVwithmuons = TheDimuon_->userInt("iPV");

      dimuon_vertexWeight = TheDimuon_->userFloat("vertexWeight");
      dimuon_vProb        = TheDimuon_->userFloat("vProb");
      dimuon_vChi2        = TheDimuon_->userFloat("vNChi2");
      dimuon_DCA          = TheDimuon_->userFloat("DCA");
      dimuon_ctauPV       = TheDimuon_->userFloat("ppdlPV");
      dimuon_ctauErrPV    = TheDimuon_->userFloat("ppdlErrPV");
      dimuon_cosAlpha     = TheDimuon_->userFloat("cosAlpha");

      const reco::Candidate::LorentzVector vP = TheDimuon_->daughter("muon1")->p4();
      const reco::Candidate::LorentzVector vM = TheDimuon_->daughter("muon2")->p4();
      if (TheDimuon_->daughter("muon1")->charge() > 0) {
      	 muonp_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
      	 muonn_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
      } else {
         muonn_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
         muonp_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
      }

      double QValue = candidate_p4.M() - dimuon_p4.M();
      invm1Skk = QValue + ups1SMass;
      invm2Skk = QValue + ups2SMass;

      TheTree->Fill();
      if (OnlyBest_) break;
    }
  } else std::cout<< "No candidate information " << run << "," << event <<std::endl;

  if (TheUps.isValid() && !TheUps->empty()) {
    pat::CompositeCandidate TheUps_;
    for (unsigned int i=0; i< TheUps->size(); i++){
      TheUps_       = TheUps->at(i);

      ups_p4.SetPtEtaPhiM(TheUps_.pt(),TheUps_.eta(),TheUps_.phi(),TheUps_.mass());

      ups_vertexWeight = TheUps_.userFloat("vertexWeight");
      ups_vProb        = TheUps_.userFloat("vProb");
      ups_vMass        = TheUps_.userFloat("vMass");
      ups_vChi2        = TheUps_.userFloat("vNChi2");
      ups_DCA          = TheUps_.userFloat("DCA");
      ups_ctauPV       = TheUps_.userFloat("ppdlPV");
      ups_ctauErrPV    = TheUps_.userFloat("ppdlErrPV");
      ups_cosAlpha     = TheUps_.userFloat("cosAlpha");

      ups_lxyPV        = TheUps_.userFloat("lxyPV");
      ups_lxyErrPV     = TheUps_.userFloat("lxyErrPV");
      ups_ctauBS       = TheUps_.userFloat("ppdlBS");
      ups_ctauErrBS    = TheUps_.userFloat("ppdlErrBS");
      ups_lxyBS        = TheUps_.userFloat("lxyBS");
      ups_lxyErrBS     = TheUps_.userFloat("lxyErrBS");

      // iPVwithmuons_ups = TheUps_.userInt("iPV");

      const reco::Candidate::LorentzVector muP = TheUps_.daughter("muon1")->p4();
      const reco::Candidate::LorentzVector muM = TheUps_.daughter("muon2")->p4();
      if (TheUps_.daughter("muon1")->charge() > 0) {
         muonP_p4.SetPtEtaPhiM(muP.pt(), muP.eta(), muP.phi(), muP.mass());
         muonN_p4.SetPtEtaPhiM(muM.pt(), muM.eta(), muM.phi(), muM.mass());
      } else {
         muonN_p4.SetPtEtaPhiM(muP.pt(), muP.eta(), muP.phi(), muP.mass());
         muonP_p4.SetPtEtaPhiM(muM.pt(), muM.eta(), muM.phi(), muM.mass());
      }

      //double testPt1 = 0.;
      if (TheUps_.userInt("mu1_charge") > 0) {
        //testPt1 = TheUps_.userFloat("mu1_pt");
        mu1_pt      = TheUps_.userFloat("mu1_pt");
        mu1_ptErr   = TheUps_.userFloat("mu1_ptErr");
        mu1_d0      = TheUps_.userFloat("mu1_d0");
        mu1_d0Err   = TheUps_.userFloat("mu1_d0Err");
        mu1_dz      = TheUps_.userFloat("mu1_dz");
        mu1_dzErr   = TheUps_.userFloat("mu1_dzErr");
        mu1_dxy     = TheUps_.userFloat("mu1_dxy");
        mu1_dxyErr  = TheUps_.userFloat("mu1_dxyErr");
        mu1_nvsh    = TheUps_.userInt("mu1_nvsh");
        mu1_nvph    = TheUps_.userInt("mu1_nvph");
        mu1_charge  = TheUps_.userInt("mu1_charge");

        mu2_pt      = TheUps_.userFloat("mu2_pt");
        mu2_ptErr   = TheUps_.userFloat("mu2_ptErr");
        mu2_d0      = TheUps_.userFloat("mu2_d0");
        mu2_d0Err   = TheUps_.userFloat("mu2_d0Err");
        mu2_dz      = TheUps_.userFloat("mu2_dz");
        mu2_dzErr   = TheUps_.userFloat("mu2_dzErr");
        mu2_dxy     = TheUps_.userFloat("mu2_dxy");
        mu2_dxyErr  = TheUps_.userFloat("mu2_dxyErr");
        mu2_nvsh    = TheUps_.userInt("mu2_nvsh");
        mu2_nvph    = TheUps_.userInt("mu2_nvph");
        mu2_charge  = TheUps_.userInt("mu2_charge");
      } else {
        //testPt1 = TheUps_.userFloat("mu2_pt");
        mu1_pt      = TheUps_.userFloat("mu2_pt");
        mu1_ptErr   = TheUps_.userFloat("mu2_ptErr");
        mu1_d0      = TheUps_.userFloat("mu2_d0");
        mu1_d0Err   = TheUps_.userFloat("mu2_d0Err");
        mu1_dz      = TheUps_.userFloat("mu2_dz");
        mu1_dzErr   = TheUps_.userFloat("mu2_dzErr");
        mu1_dxy     = TheUps_.userFloat("mu2_dxy");
        mu1_dxyErr  = TheUps_.userFloat("mu2_dxyErr");
        mu1_nvsh    = TheUps_.userInt("mu2_nvsh");
        mu1_nvph    = TheUps_.userInt("mu2_nvph");
        mu1_charge  = TheUps_.userInt("mu2_charge");

        mu2_pt      = TheUps_.userFloat("mu1_pt");
        mu2_ptErr   = TheUps_.userFloat("mu1_ptErr");
        mu2_d0      = TheUps_.userFloat("mu1_d0");
        mu2_d0Err   = TheUps_.userFloat("mu1_d0Err");
        mu2_dz      = TheUps_.userFloat("mu1_dz");
        mu2_dzErr   = TheUps_.userFloat("mu1_dzErr");
        mu2_dxy     = TheUps_.userFloat("mu1_dxy");
        mu2_dxyErr  = TheUps_.userFloat("mu1_dxyErr");
        mu2_nvsh    = TheUps_.userInt("mu1_nvsh");
        mu2_nvph    = TheUps_.userInt("mu1_nvph");
        mu2_charge  = TheUps_.userInt("mu1_charge");
      }

      //std::cout<<" ===> pt = "<<muonP_p4.Pt()<<" single one = "<<testPt1<<std::endl;

      UpsTree->Fill();
    }
  } else std::cout<< "No Upsilon candidate information " << run << "," << event <<std::endl;



}

void OniaRecoTrackTrackRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(OniaRecoTrackTrackRootupler);
