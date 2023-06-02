# µµTrkTrk reconstruction in MiniAOD 

## The setup
This package is mean to be run using Ultra Legacy MINIAODv2

* Setup in `13_1_0` 

```
scram p -n cmssw CMSSW_13_1_0
cd cmssw/src/
cmsenv
git clone git@github.com:AdrianoDee/2mu2trk_exercise.git -b jpsiphi MuMuTrkTrk/MuMuTrkTrk/
scram b -j 8
```

To run the example on `/RelValBsToJpsiPhi_mumuKK_14TeV/CMSSW_13_0_0_pre3-130X_mcRun3_2022_realistic_v2-v1/MINIAODSIM`

If you haven't, set up your certificate to access the data:
```
voms-proxy-init -rfc -voms cms -valid 192:00
```

Then simply run the config
```
cmsRun MuMuTrkTrk/MuMuTrkTrk/test/run-jpsikk-miniaodsim.py
```

----

## The Decay

The decay we want to reconstruct a $B^0_s$ meson candidate decaying as:

$$B_s^0 \to J/\psi(\to \mu\mu)\phi(\to KK)  $$

a bottom strange meson ($s\bar{b}$). You can find further information in the [PDG](https://pdg.lbl.gov/2023/web/viewer.html?file=../tables/rpp2023-tab-mesons-bottom-strange.pdf).

![Bs meson table from pdg](images/bs_meson.png "Bs meson table from pdg")

We can sketch the decay topology as follows:

![BsJpsiPhi](images/bs_jps_phi.png "BsJpsiPhi")

the $B_s^0$ candidate, before decaying, flies (given the non negligible lifetime) w.r.t the collision area (beam spot) and it decays in:
- $J/\psi$ meson, a $c\bar{c}$ bound state with a mass of $\sim 3096 MeV$(see the [PDG](https://pdg.lbl.gov/2023/web/viewer.html?file=../tables/rpp2023-tab-mesons-c-cbar.pdf) for further details) 
- $\phi$ meson ($s\bar{s}$) a mass of $\sim 1020 MeV$(see the [PDG](https://pdg.lbl.gov/2023/tables/contents_tables_mesons.html ) for further details).

The $J/psi$ and the $\phi$ then furhterly decay into a pair of opposite sign muons and a pair of opposite sign kaons, almost promptly w.r.t. the $B^0_s$ meson.

## The Analyzer

## The Config

We select the global tag suitable for the dataset we have
```python
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v2', '')
```
this is something you can usually get from the dataset name itself.

Next step let's configure the inputs and the outputs with the `PoolSource` module and the `TFileService` (used to output a simple `ROOT` file).

```python
input_filename = ["/store/relval/CMSSW_13_0_0_pre3/RelValBsToJpsiPhi_mumuKK_14TeV/MINIAODSIM/130X_mcRun3_2022_realistic_v2-v1/00000/d26ae0bf-3a06-4098-9c24-c29452079aa0.root"]
ouput_filename = 'mumukk.root'
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
```

Then we may select the events that has fired a specific trigger:

```python
triggers = [
'HLT_DoubleMu4_JpsiTrkTrk_Displaced', ##Run3 trigger!
]
```
this is a trigger dedicated to $J/\psi + 2 Tracks$ displaced topologies (exactly what we need!). The $4$ there is the $p_T$ threshold on the $2\mu$ system. So we add to the process a module to filter the events:

```
hltpathsV = cms.vstring([h + '_v*' for h in triggers ])
process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = hltpathsV,
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

```

In `MINIAOD` the muons are stored under the `slimmedMuons` collection of `pat::Muons`. Instead of using all of them we may run a first selection based on a combination of cuts with the `PATMuonSelector` module: 

```python
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
```
In this case, e.g., we are selecting muons:
- `innerTrack.hitPattern.pixelLayersWithMeasurement>0`, whose track has at least one pixel on at least one layer;
- `innerTrack.hitPattern.pixelLayersWithMeasurement>0`, whose track has at least one hit on at least five layers;
- `innerTrack.quality(\"highPurity\")`, an high purity track
- `pt > 2.`, etc...


Now it's the moment to build the µµ candidate:

```python
### µµ building
process.load("MuMuTrkTrk.MuMuTrkTrk.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('oniaSelectedMuons') #using as input the selected muons
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices') #the PVs
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot') #the Bs
process.onia2MuMuPAT.higherPuritySelection=cms.string("") 
process.onia2MuMuPAT.lowerPuritySelection=cms.string("")
process.onia2MuMuPAT.dimuonSelection=cms.string("2.5 < mass && mass < 3.5") #mass window cuts
process.onia2MuMuPAT.addMCTruth = cms.bool(False) ## not interested for the moment
```


```python
process.DiMuonCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFiltered"),
    minNumber = cms.uint32(1),
)
```


```python

process.OniaPseudoTrackTrackCandidateProducer = cms.EDProducer('OniaPseudoTrackTrackProducer',
    Onia = cms.InputTag("Onia2MuMuFiltered"),
    BeamSpot = cms.InputTag('offlineBeamSpot'),
    Track = cms.InputTag("packedPFCandidates"),
    OniaMassCuts = cms.vdouble(2.9,3.3), 
    CandidateMassCuts = cms.vdouble(4.0,6.0),
    Track1Mass = cms.double(0.493677),#kaons
    Track2Mass = cms.double(0.493677),#kaons
    ConstraintMass = cms.double(3.096916),#J/Psi
)

```

## The Notebook

The `bs_decay.ipynb` notebook under `MuMuTrkTrk/MuMuTrkTrk/test/`

It's easier if you follow the instructions directly there. If you want you may run this notebook in [SWAN](https://swan.web.cern.ch/swan/) from which you will be able to acces any area you have acces on `eos` (e.g. your `lxplus` home, in my case `/afs/cern.ch/user/a/adiflori`).


 <a href="https://cern.ch/swanserver/cgi-bin/go?projurl=https://raw.githubusercontent.com/AdrianoDee/2mu2trk_exercise/main/test/bs_decay.ipynb" target="_blank">
                            <img src="https://swanserver.web.cern.ch/swanserver/images/badge_swan_white_150.png" alt="Open in SWAN" style="height:1.5em">
                        </a>

## To you!

As an excercise you may try to modify the analyzer chain and the python config in order to reconstruct another decay with a similar daughters:

$$\psi(2S) \to J/\psi(\to \mu\mu)\pi\pi  $$

![Psi(2s)toJpsiPiPi](images/psi2s.png "Psi(2s)toJpsiPiPi")

There are some differences:

1. in this case the decay is __prompt__, meaning we expect all the candidates to come from a PV;
2. there is no constraint on the $\pi\pi$ system mass since they are non-resonant;
3. the $\psi(2S)$ has a lower mass w.r.t. the $B^0_s$

You can try to use the dataset

`/RelValPsi2SToJPsiPiPi_14TeV/CMSSW_13_0_0_pre3-130X_mcRun3_2022_realistic_v2-v1/MINIAODSIM`

that has the same GT.

(Suggestion: also the `HLT` you select will need to be changed, if you don't know which to use, simply turn off the selection.)

Another possibility is to check the PU datasets for $B^0_s$ such as:

`/RelValBsToJpsiPhi_mumuKK_14TeV/CMSSW_13_1_0_pre4-PU_131X_mcRun3_2022_realistic_v3_2023_BPH-v1/MINIAODSIM`

and see what happens. Beware! The GT is different!
