# UpsPhi-Analysis

This package is mean to be run using Ultra Legacy MINIAODv2

* Setup: (it has being tested on 10_6_27)

```
export SCRAM_ARCH=slc7_amd64_gcc700
scram p -n CMSSW_10627_ykk CMSSW_10_6_27
cd CMSSW_10627_ykk/src/
cmsenv
git clone git@github.com:AdrianoDee/YPhi_Analysis_UL_MINIAODv2.git -b jpsiphi Ponia/Onia/
scram b

```

* Run: (use your favorite input sample)

```
voms-proxy-init -rfc -voms cms -valid 192:00
cmsRun Ponia/Onia/test/run-chib3p2upskk-miniaodsim.py (for chi_b1(3P)->Y(1S)Phi reconstruction using UL MC samples)
```
