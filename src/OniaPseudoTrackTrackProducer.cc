#include "MuMuTrkTrk/MuMuTrkTrk/interface/OniaPseudoTrackTrackProducer.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TVector3.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

OniaPseudoTrackTrackProducer::OniaPseudoTrackTrackProducer(const edm::ParameterSet& ps):
  OniaCollection_(consumes(ps.getParameter<edm::InputTag>("Onia"))),
  BeamSpot_(consumes(ps.getParameter<edm::InputTag>("BeamSpot"))),
  TrakCollection_(consumes(ps.getParameter<edm::InputTag>("Track"))),
  OniaMassCuts_(ps.getParameter<std::vector<double>>("OniaMassCuts")),
  CandidateMassCuts_(ps.getParameter<std::vector<double>>("CandidateMassCuts")),
  Track1Mass_(ps.getParameter<double>("Track1Mass")),
  Track2Mass_(ps.getParameter<double>("Track2Mass")),
  ConstraintMass_(ps.getParameter<double>("ConstraintMass"))
{
  builderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  produces<pat::CompositeCandidateCollection>();
  candidates = 0;
  nevents = 0;
}

void OniaPseudoTrackTrackProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::unique_ptr<pat::CompositeCandidateCollection> TheCandidateColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> Onias;
  event.getByToken(OniaCollection_,Onias);

  edm::Handle<pat::PackedCandidateCollection> Tracks;
  event.getByToken(TrakCollection_,Tracks);

  edm::Handle<reco::BeamSpot> theBeamSpot;
  event.getByToken(BeamSpot_, theBeamSpot);


  // edm::ESHandle<TransientTrackBuilder> theB;
  auto const &theB = esetup.getData(builderToken_);
  // esetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  for (pat::CompositeCandidateCollection::const_iterator TheOnia = Onias->begin(); TheOnia != Onias->end(); ++TheOnia) {
     float oniaM = TheOnia->mass();
     if ( oniaM < OniaMassCuts_[1]  && oniaM > OniaMassCuts_[0] ) {
       const pat::Muon *pmu1, *pmu2;
       const reco::Vertex *thePrimaryV;
       RefCountedKinematicParticle  fittedPhoton=nullptr;
       pmu1 = dynamic_cast<const pat::Muon*>(TheOnia->daughter("muon1"));
       pmu2 = dynamic_cast<const pat::Muon*>(TheOnia->daughter("muon2"));
       thePrimaryV = TheOnia->userData<reco::Vertex>("PVwithmuons");
       kaons.clear();
       kaonsPV.clear();
       refVtx.clear();
       pvAssocQ.clear();
       dzAssocPV.clear();
       //std::cout<<"track size = "<<Tracks->size()<<std::endl;

       if (Tracks.isValid() && !Tracks->empty()) {
         for (std::vector<pat::PackedCandidate>::const_iterator pp = Tracks->begin(); pp!= Tracks->end(); ++pp) {

           //if (!pp->trackHighPurity() || !pp->fromPV() || std::abs(pp->pdgId())!=321 || !pp->hasTrackDetails()) continue;

           if (!pp->trackHighPurity() || !pp->hasTrackDetails()) continue;

           //std::cout<<" ######### I have kaons for you ######### "<<std::endl;
           
           //if ( pp->vertexRef().key()!=THEKEY) continue; //This reinforce the fromPV request above
           const reco::Track* TheTrack = &pp->pseudoTrack();
           //if (TheTrack->pt()<0.5 || std::abs(TheTrack->eta())>2.5 || TheTrack->charge()==0 ) continue;

           if (TheTrack->charge()==0) continue;

           //std::cout<<" ######### K charge is not ZERO ######### "<<std::endl;
           if ( IsTheSame(*TheTrack,*pmu1) || IsTheSame(*TheTrack,*pmu2) ) continue;
           //std::cout<<" ######### Kaons does not match with Muons ######### "<<std::endl;

           kaons.push_back(*TheTrack);
           kaonsPV.push_back(pp->fromPV());
           refVtx.push_back(pp->vertexRef().key());
           pvAssocQ.push_back(pp->pvAssociationQuality());
           dzAssocPV.push_back(pp->dzAssociatedPV());
           //std::cout<<" ######### PV Association Quality = "<< pp->pvAssociationQuality() <<" ######### "<<" PV id = " << pp->fromPV() <<" ######### "<<" Ref. Vertex i = " << pp->vertexRef().key() <<" ######### "<<std::endl;
           //std::cout<< " dzAssociatedPV = " << pp->dzAssociatedPV() << std::endl;
         }
       }

       //std::sort(kaons.begin(),kaons.end(),PtComparator);

       for (unsigned int ii=0; ii<kaons.size();++ii) {
         reco::Track* Track1 = &kaons.at(ii);
         for (unsigned int jj=ii+1; jj<kaons.size();++jj) {
            reco::Track* Track2 = &kaons.at(jj);
            if (Track1 == Track2) continue;

            //if(Track1->charge() * Track2->charge() >= 0) continue; // Right sign
            //if(Track1->charge() * Track2->charge() <= 0) continue; // Wrong sign

            //std::cout<<" charge = "<<Track1->charge() * Track2->charge()<<std::endl;
            float deltaR_track1 = std::sqrt(reco::deltaR2(*Track1,*TheOnia));
            float deltaR_track2 = std::sqrt(reco::deltaR2(*Track2,*TheOnia));

            double trk1_e = sqrt(Track1Mass_*Track1Mass_ + Track1->p()*Track1->p());
            double trk2_e = sqrt(Track2Mass_*Track2Mass_ + Track2->p()*Track2->p());
            reco::Candidate::LorentzVector Trk1_p4(Track1->px(),Track1->py(),Track1->pz(),trk1_e);
            reco::Candidate::LorentzVector Trk2_p4(Track2->px(),Track2->py(),Track2->pz(),trk2_e);
            reco::Candidate::LorentzVector ditrack_p4 = Trk1_p4 + Trk2_p4;
            if (ditrack_p4.M()<0.98 || ditrack_p4.M()>1.06) continue;

            //std::cout<<" ######### You have Kaons in phi mass range ######### "<<std::endl;
            float deltaR_ditrack = std::sqrt(reco::deltaR2(ditrack_p4,*TheOnia));
            pat::CompositeCandidate TheCandidate;
            TheCandidate.addDaughter(*TheOnia,"onia");
            TheCandidate.addUserData<reco::Track>( "track1", *Track1 );
            TheCandidate.addUserData<reco::Track>( "track2", *Track2 );

            TheCandidate.addUserFloat("ditrack_dRdimuon",deltaR_ditrack);
            TheCandidate.addUserFloat("track1_dRdimuon",deltaR_track1);
            TheCandidate.addUserFloat("track2_dRdimuon",deltaR_track2);
            TheCandidate.addUserInt("track1_PV",kaonsPV[ii]);
            TheCandidate.addUserInt("track2_PV",kaonsPV[jj]);
            TheCandidate.addUserInt("track1_refVtx",refVtx[ii]);
            TheCandidate.addUserInt("track2_refVtx",refVtx[jj]);
            TheCandidate.addUserInt("track1_pvAssocQ",pvAssocQ[ii]);
            TheCandidate.addUserInt("track2_pvAssocQ",pvAssocQ[jj]);
            TheCandidate.addUserFloat("track1_dzAssocPV",dzAssocPV[ii]);
            TheCandidate.addUserFloat("track2_dzAssocPV",dzAssocPV[jj]);

            TheCandidate.setCharge(Track1->charge()+Track2->charge());
            reco::Candidate::LorentzVector vCandidate = TheOnia->p4() + ditrack_p4;
            TheCandidate.setP4(vCandidate);
            if ( TheCandidate.mass() > CandidateMassCuts_[1] || TheCandidate.mass() < CandidateMassCuts_[0]) continue;
            //std::cout<<" ######### You have some possible candidates ######### "<<std::endl;
            std::vector<reco::TransientTrack> MuMuTk;
            MuMuTk.push_back(theB.build(*pmu1->innerTrack()));
            MuMuTk.push_back(theB.build(*pmu2->innerTrack()));
            MuMuTk.push_back(theB.build(*Track1));
            MuMuTk.push_back(theB.build(*Track2));

            KinematicParticleFactoryFromTransientTrack pFactory;
            const ParticleMass muMass(0.1056583);
            float muSigma = muMass*1E-6;
            const ParticleMass tk1Mass(Track1Mass_);
            float tk1Sigma = tk1Mass*1E-6;
            const ParticleMass tk2Mass(Track2Mass_);
            float tk2Sigma = tk2Mass*1E-6;

            ParticleMass mass_(ConstraintMass_);
            std::vector<RefCountedKinematicParticle> allDaughters_;
            allDaughters_.push_back(pFactory.particle (MuMuTk[0], muMass, float(0), float(0), muSigma));
            allDaughters_.push_back(pFactory.particle (MuMuTk[1], muMass, float(0), float(0), muSigma));
            allDaughters_.push_back(pFactory.particle (MuMuTk[2], tk1Mass, float(0), float(0), tk1Sigma));
            allDaughters_.push_back(pFactory.particle (MuMuTk[3], tk2Mass, float(0), float(0), tk2Sigma));

            KinematicConstrainedVertexFitter constVertexFitter;
            MultiTrackKinematicConstraint *onia_mtc = new  TwoTrackMassKinematicConstraint(mass_);
            RefCountedKinematicTree TheParticleTree = constVertexFitter.fit(allDaughters_,onia_mtc);

            TVector3 pvtx(thePrimaryV->position().x(),thePrimaryV->position().y(),0);
            TVector3 pvtx3D(thePrimaryV->position().x(),thePrimaryV->position().y(),thePrimaryV->position().z());

            if (TheParticleTree->isEmpty()) continue;
            TheParticleTree->movePointerToTheTop();
            RefCountedKinematicParticle TheParticle = TheParticleTree->currentParticle();
            RefCountedKinematicVertex TheDecayVertex = TheParticleTree->currentDecayVertex();
            if (!TheParticle->currentState().isValid()) continue;
            //std::cout<<" ######### TheDecayVertex = "<<TheDecayVertex->position()<<" ######### "<<std::endl;
            if (TheParticle->currentState().mass() > CandidateMassCuts_[1] || TheParticle->currentState().mass() < CandidateMassCuts_[0]) continue;
            double Prob = ChiSquaredProbability((double)(TheDecayVertex->chiSquared()),(double)(TheDecayVertex->degreesOfFreedom()));
            if(Prob<0.005) continue;

            //std::cout<<" ######### You have CANDIDATES ######### "<<std::endl;

            TheCandidate.addUserFloat("vMass",TheParticle->currentState().mass());
            TheCandidate.addUserFloat("vChi2",TheDecayVertex->chiSquared());

            TVector3 vtx(TheDecayVertex->position().x(),TheDecayVertex->position().y(),0);
            TVector3 pperp(TheParticle->currentState().kinematicParameters().momentum().x(), TheParticle->currentState().kinematicParameters().momentum().y(), 0);
            AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
            TVector3 vdiff = vtx - pvtx;
            double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
            GlobalError v1e = (reco::Vertex(*TheDecayVertex)).error();
            GlobalError v2e = thePrimaryV->error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix() + v2e.matrix();

            float lxy = vdiff.Perp();
            ROOT::Math::SVector<double, 3> vDiff; // needed by Similarity method
            vDiff[0] = vdiff.x(); vDiff[1] = vdiff.y(); vDiff[2] = 0; // needed by Similarity method
            float lxyErr = sqrt(ROOT::Math::Similarity(vDiff,vXYe)) / vdiff.Perp();

            VertexDistanceXY vdistXY;
            Measurement1D distXY = vdistXY.distance(reco::Vertex(*TheDecayVertex),* thePrimaryV);
            double ctauPV = distXY.value() * cosAlpha * TheParticle->currentState().mass() / pperp.Perp();
            double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe)) * TheParticle->currentState().mass() / pperp.Perp2();

            TVector3 vtx3D(TheDecayVertex->position().x(),TheDecayVertex->position().y(),TheDecayVertex->position().z());
            //std::cout<<" ######### TheDecayVertex 3D = "<<vtx3D.Mag()<<" ######### PVwithmuons 3D = "<<pvtx3D.Mag()<<" ######### "<<std::endl;
            //double Pos = sqrt(TheDecayVertex->position().x()*TheDecayVertex->position().x() + TheDecayVertex->position().y()*TheDecayVertex->position().y() + TheDecayVertex->position().z()*TheDecayVertex->position().z());
            //std::cout<<" ######### TheDecayVertex 3D manuel = "<<Pos<<" #########TheDecayVertex 3D = "<<vtx3D.Mag()<<" ######### "<<std::endl;
            /*double pVx = TheDecayVertex->position().x() - thePrimaryV->position().x();
            double pVy = TheDecayVertex->position().y() - thePrimaryV->position().y();
            double pVz = TheDecayVertex->position().z() - thePrimaryV->position().z();

            double pV3D = sqrt(pVx*pVx + pVy*pVy + pVz*pVz);*/

            VertexDistance3D vdistXYZ;
            Measurement1D distXYZ = vdistXYZ.distance(reco::Vertex(*TheDecayVertex),* thePrimaryV);

            TVector3 pperp3D(TheParticle->currentState().kinematicParameters().momentum().x(), TheParticle->currentState().kinematicParameters().momentum().y(), TheParticle->currentState().kinematicParameters().momentum().z());
            TVector3 vdiff3D = vtx3D - pvtx3D;
            //double dis = vtx3D.Mag()
            //std::cout<<" ######### Vertex distance 3V = "<<vdiff3D.Mag()<<" ######### Vertex distance 3D = "<<distXYZ.value()<<" ######### Vertex distance 3D manuel = "<<pV3D<<" ######### "<<std::endl;
            double cosAlpha3D = vdiff3D.Dot(pperp3D) / (vdiff3D.Mag() * pperp3D.Mag());

            float lxyz = vdiff3D.Mag();
            ROOT::Math::SVector<double, 3> vDiff3D; // needed by Similarity method
            vDiff3D[0] = vdiff3D.x(); vDiff3D[1] = vdiff3D.y(); vDiff3D[2] = vdiff3D.z(); // needed by Similarity method
            float lxyzErr = sqrt(ROOT::Math::Similarity(vDiff3D,vXYe)) / vdiff3D.Mag();

            TheCandidate.addUserInt("vStatus",1);
            TheCandidate.addUserFloat("vProb",Prob);
            TheCandidate.addUserFloat("cosAlpha",cosAlpha);
            TheCandidate.addUserFloat("cosAlpha3D",cosAlpha3D);
            TheCandidate.addUserFloat("ctauPV",ctauPV);
            TheCandidate.addUserFloat("ctauErrPV",ctauErrPV);
            TheCandidate.addUserFloat("lxy",lxy);
            TheCandidate.addUserFloat("lxyErr",lxyErr);
            TheCandidate.addUserFloat("lxyz",lxyz);
            TheCandidate.addUserFloat("lxyzErr",lxyzErr);
            TheCandidate.addUserFloat("thePrimaryV_X",thePrimaryV->position().x());
            TheCandidate.addUserFloat("thePrimaryV_Y",thePrimaryV->position().y());
            TheCandidate.addUserFloat("thePrimaryV_Z",thePrimaryV->position().z());
            TheCandidate.addUserFloat("TheDecayVertex_X",TheDecayVertex->position().x());
            TheCandidate.addUserFloat("TheDecayVertex_Y",TheDecayVertex->position().y());
            TheCandidate.addUserFloat("TheDecayVertex_Z",TheDecayVertex->position().z());
            TheCandidate.addUserFloat("thePrimaryV_2D_position",pvtx.Mag());
            TheCandidate.addUserFloat("thePrimaryV_3D_position",pvtx3D.Mag());
            TheCandidate.addUserFloat("TheDecayVertex_2D_position",vtx.Mag());
            TheCandidate.addUserFloat("TheDecayVertex_3D_position",vtx3D.Mag());
            TheCandidate.addUserFloat("TheVertexDistance_2D",distXY.value());
            TheCandidate.addUserFloat("TheVertexDistance_3D",distXYZ.value());

            TheCandidate.setVertex(reco::Vertex(*TheDecayVertex).position());
            TheCandidateColl->push_back(TheCandidate);
            candidates++;
         } //Second track

       } //First

     }
  }
  OniaPseudoTrackTrackProducer::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  std::sort(TheCandidateColl->begin(),TheCandidateColl->end(),vPComparator);
  event.put(std::move(TheCandidateColl));
  nevents++;
}

void OniaPseudoTrackTrackProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "OniaPseudoTrackTrackProducer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " X_b candidates." << std::endl;
  std::cout << "###########################" << std::endl;
 /* std::cout << "dz max  " << dz << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "dxy max  " << dxy << std::endl;
  std::cout << "###########################" << std::endl;
  */
  //std::cout<< " PV ---- n0 = " << fpv0 << ", n1 = " << fpv1 << ", n2 = " << fpv2 << ", n3 = " << fpv3 << ", n23 = " << fpv23 << std::endl;
  //std::cout << "###########################" << std::endl;
}

const bool OniaPseudoTrackTrackProducer::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  return std::abs(mu.eta()-tk.eta()) < 0.02 && std::abs(mu.p()-tk.p()) < 0.02;
}

const bool OniaPseudoTrackTrackProducer::IsTheSame(const reco::Track& tk1, const reco::Track& tk2){
  return std::abs(tk1.eta()-tk2.eta()) < 0.02 && std::abs(tk1.p()-tk2.p()) < 0.02;
}

DEFINE_FWK_MODULE(OniaPseudoTrackTrackProducer);
