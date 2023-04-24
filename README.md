# EgammaAnalysis-TnPTreeProducer

Package of the EGamma group to produce Tag-and-Probe trees

## :fire: Dedicated changes to HToEleEleGamma analysis
---
The official TnPTreeProducer is modified to make Tnp Ntuples for measuring the DiPhoton30_22(18) trigger scale facters used for HToEleEleGamma analysis. The main changes are as follows.

- `plugins/MergedElectronMvaProducer.cc`: To add the variables needed for merged electron ID.
- `python/egmTreesSetup_cff.py`: additional process are add to require tag electron match to seed filters of DiPhoton30_22 HLT for unseed leg Ntuples. (`process.tagEleLeadMatch`)

The filters for DiPhoton30_22 HLT
- Seeded leg 
    - hltEG30LIso60CaloId15b35eHE12R9Id50b80eEcalIsoLastFilter
    - hltEG30LR9Id85b90eHE12R9Id50b80eR9IdLastFilter
- Unseeded leg
    - hltEG22R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter (2017 and 2018)
    - hltEG22Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter (2017 and 2018)
    - hltEG18R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter (2016)
    - hltEG18Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter (2016)

The way to add filters is referred to `flashgg`:
   - https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Validation/test/makeHltTreePhotons.py
   - https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Validation/python/treeMakerOptionsHLT_cfi.py

Bacause the variables needed for merged electron ID contains the variable which reuires EGamma post processing after reconstruction, the `EgammaPostRecoTools` need to be set up. To setup the TnPTreeProducer used for HToEleEleGamma analysis, please use the following commands.
```bash
# install the TnPTreeProducer
cmsrel CMSSW_10_6_13
cd CMSSW_10_6_13/src
cmsenv
git clone -b RunIIfinal https://github.com/chw1207/EgammaAnalysis-TnPTreeProducer.git EgammaAnalysis/TnPTreeProducer

# install the EgammaPostRecoTools
git cms-addpkg RecoEgamma/EgammaTools
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools

# compile
scram b -j8
```

To submit the crab jods, please executes
```bash
source /cvmfs/cms.cern.ch/common/crab-setup.sh
python MyTnpCrabSubmit.py
```
---
## Overview of branches

| Branch                                     | release            | tnpEleIDs          | tnpPhoIDs          | tnpEleTrig         | tnpEleReco         | purpose                                |
| ------------------------------------------ | ------------------ |:------------------:|:------------------:|:------------------:|:------------------:|:--------------------------------------:|
|                                            |                    | *miniAOD*          | *miniAOD*          | *miniAOD*          | *AOD*              |                                        |
| [RunIIfinal](../../tree/RunIIfinal)        | CMSSW\_10\_2       | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | Run II analysis                        |
| [RunIIfinal](../../tree/RunIIfinal)        | CMSSW\_10\_6       | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | Run II analysis using ultra-legacy     |
| [CMSSW\_11\_X\_Y](../../tree/CMSSW_11_X_Y) | CMSSW\_11          | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :white_check_mark: | Development for Run III (experimental) |

Note: because of a dataformat CMSSW\_10\_6 can only be used for ultra-legacy samples, and CMSSW\_10\_2 should be used for the rereco samples.

## Available tuples
### ReReco 2016, 2017 and 2018
If you do not need changes to the default code, you can simply use existing flat tag and probe trees, avalaible for both 2016, 2017 and 2018 (RunIIfinal branch):

```bash
ls /eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-06-09/*/merged/
```

These inlcude the tnpEleTrig, tnpEleIDs and tnpPhoIDs trees produced with the RunIIfinal branch.
*Main change with respect to the 2020-02-28 production is the inclusion of some additional branches, e.g. the leptonMva's*

### ReReco 2016, 2017 and 2018 - L1 matched
In case you need L1 matching for the measurement of doubleEle HLT triggers, you can use the tnpEleTrig trees found in:

```bash
ls /eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-03-03/*/merged/*L1matched.root
```

### UL2017 and UL2018
For ultra-legacy  we have tnpEleTrig, tnpEleIDs and tnpPhoIDs trees available at:
```
ls /eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2018/merged
ls /eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2017/merged
ls /eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016postVFP/merged
ls /eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged
```


## To produce new tuples
### 1a. Install for rereco (CMSSW\_10\_2\_X with X=10 or higher, works for 2016, 2017 and 2018 data/MC)

```bash
cmsrel CMSSW_10_2_22
cd CMSSW_10_2_22/src
cmsenv
git clone -b RunIIfinal git@github.com:cms-egamma/EgammaAnalysis-TnPTreeProducer.git EgammaAnalysis/TnPTreeProducer
scram b -j8
```

### 1b. Install for ultra-legacy (CMSSW\_10\_6\_X, works for UL2017 and UL2018 data/MC)

```bash
cmsrel CMSSW_10_6_13
cd CMSSW_10_6_13/src
cmsenv
git clone -b RunIIfinal git@github.com:cms-egamma/EgammaAnalysis-TnPTreeProducer.git EgammaAnalysis/TnPTreeProducer
scram b -j8
```

### 2. Try-out
You can find the cmsRun executable in EgammaAnalysis/TnPTreeProducer/python:
```bash
cd EgammaAnalysis/TnPTreeProducer/python/
cmsRun TnPTreeProducer_cfg.py isMC=True doTrigger=True era=UL2018
```
Check [TnPTreeProducer\_cfg.py](python/TnPTreeProducer_cfg.py) for all available options. Update the code if you need to implement custom-made recipes.

Test files can be defined in [python/etc/tnpInputTestFiles\_cff.py](python/etc/tnpInputTestFiles_cff.py)
If you update the code, you can use the ./runTests.py script in the test directory to check for new differences in the 2016, 2017 and 2018 test files.

### 3. Submit jobs
Check in EgammaAnalysis/TnPTreeProducer//crab the tnpCrabSubmit.py script to submit your jobs using crab

```bash
source /cvmfs/cms.cern.ch/common/crab-setup.sh
```

## To make a pull request to this repository
1. On github fork the package https://github.com/cms-analysis/EgammaAnalysis-TnPTreeProducer
2. Add the remote
```bash
git remote add username-push git@github.com:username/EgammaAnalysis-TnPTreeProducer.git
```
3. push commits to fork and then standard pull request process
```bash
git push username-push branchname
```

## Adding new workingpoints
You can add new electron workingpoints in [python/egmElectronIDModules\_cff.py](python/egmElectronIDModules_cff.py) and new photon workingpoints
in [python/egmPhotonIDModules\_cff.py](python/egmPhotonIDModules_cff.py). Each new workingpoint added in these python config fragments will
add a new "passing<WP>" boolean in the electron and photon trees respectively. Of course, one can also choose to simply add a variable in
[python/egmTreesContent\_cff.py](python/egmTreesContent\_cff.py), which might be preferred for MVA variables when you want to have the
flexibility to explore different workingpoints: you can simply put a cut on these variable in the egm\_tnp\_analysis package.

## Description of variables

Description of some of variables in the output tree is given [here](VariablesInfo.md).
