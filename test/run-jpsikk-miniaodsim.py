from glob import glob
import FWCore.ParameterSet.Config as cms
process = cms.Process('MuMuTrkTrk')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v2', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))
process.options.numberOfConcurrentLuminosityBlocks = 1

#from Ponia.Onia.run2_UL_miniAODv2_files import UL2018MC

import FWCore.ParameterSet.VarParsing as VarParsing

input_filename = [#"/store/relval/CMSSW_13_0_0_pre3/RelValBsToJpsiPhi_mumuKK_14TeV/MINIAODSIM/130X_mcRun3_2022_realistic_v2-v1/00000/d26ae0bf-3a06-4098-9c24-c29452079aa0.root"]
    "file:/lustre/home/adrianodif/analyses/das_exercise/cmssw/d26ae0bf-3a06-4098-9c24-c29452079aa0.root"]
ouput_filename = 'mumukk.root'

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))

process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.options.numberOfThreads=cms.untracked.uint32(16)
process.options.numberOfStreams=cms.untracked.uint32(16)

process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuons'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && (pt > 2.)'
   ),
   filter = cms.bool(True)
)

process.load("Ponia.Onia.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('oniaSelectedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("")
process.onia2MuMuPAT.dimuonSelection=cms.string("2.5 < mass && mass < 3.5")
process.onia2MuMuPAT.addMCTruth = cms.bool(False)

triggers = [
'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
]

hltpathsV = cms.vstring([h + '_v*' for h in triggers ])

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = hltpathsV,
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.Onia2MuMuFiltered = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("onia2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.5 < mass && mass < 3.5 && charge==0 && userFloat('vProb') > 0.01"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = hltpathsV,
)

process.DiMuonCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFiltered"),
    minNumber = cms.uint32(1),
)

process.OniaPseudoTrackTrackCandidateProducer = cms.EDProducer('OniaPseudoTrackTrackProducer',
    Onia = cms.InputTag("Onia2MuMuFiltered"),
    BeamSpot = cms.InputTag('offlineBeamSpot'),
    Track = cms.InputTag("packedPFCandidates"),
    OniaMassCuts = cms.vdouble(2.9,3.3), #J
    CandidateMassCuts = cms.vdouble(4.0,6.0),
    Track1Mass = cms.double(0.493677),#kaons
    Track2Mass = cms.double(0.493677),#kaons
    ConstraintMass = cms.double(3.096916),#J/Psi
)

process.candidateCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("OniaPseudoTrackTrackCandidateProducer"),
    minNumber = cms.uint32(1),
)

process.jpsitrktrkSequence = cms.Sequence(
                                   process.triggerSelection *
                                   process.oniaSelectedMuons *
                                   process.onia2MuMuPAT*
                                   process.Onia2MuMuFiltered *
                                   process.DiMuonCounter *
                                   process.OniaPseudoTrackTrackCandidateProducer *
                                   process.candidateCounter
                                   )

process.rootuple = cms.EDAnalyzer('OniaRecoTrackTrackRootupler',
                          TheCandidates = cms.InputTag("OniaPseudoTrackTrackCandidateProducer"),
                          TheUps = cms.InputTag("Onia2MuMuFiltered"),
                          PrimaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          GenParticles = cms.InputTag("prunedGenParticles"),
                          Track1Mass = cms.double(0.493677),
                          Track2Mass = cms.double(0.493677),
                          DimuonMass = cms.double(3.096916),#J/Psi
                          candidate_pdgid = cms.uint32(531),
                          onia_pdgid = cms.uint32(443),
                          ditrack_pdgid = cms.uint32(333),
                          track1_pdgid = cms.int32(321),
                          track2_pdgid = cms.int32(-321),
                          isMC = cms.bool(True),
                          OnlyBest = cms.bool(False)
)

process.p = cms.Path(process.jpsitrktrkSequence*process.rootuple)
