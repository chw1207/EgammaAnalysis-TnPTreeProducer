import FWCore.ParameterSet.Config as cms

#
# Sequence to add merged electron MVA
#
def mergedMvaSequence(process, options, tnpVars):
    process.mergedMva = cms.EDProducer("MergedElectronMvaProducer",
        # weightFileEB    = cms.FileInPath("EgammaAnalysis/TnPTreeProducer/data/Models/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EB/XGB/XGB_modelXGB.txt"),
        # weightFileEE    = cms.FileInPath("EgammaAnalysis/TnPTreeProducer/data/Models/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EE/XGB/XGB_modelXGB.txt"),
        probes          = cms.InputTag("slimmedElectrons"),
        gsfTracks       = cms.InputTag("reducedEgamma:reducedGsfTracks"),
        rhoInputTag     = cms.InputTag("fixedGridRhoFastjetAll"),
        photons         = cms.InputTag("slimmedPhotons"), # used for Hgg preselection (without conversion safe electron veto)
        rhoAllInputTag  = cms.InputTag("fixedGridRhoAll") # used for Hgg preselection (rho correction for isolation)
    )
    mergedMva_sequence = cms.Sequence(process.mergedMva)
    
    #
    # Adding the new variables to the trees
    #
    newVariables = {
        # "el_mergedMva"      : cms.InputTag("mergedMva:mergedMVA"),
        "el_pterr"          : cms.InputTag("mergedMva:pterr"),
        "el_ntks"           : cms.InputTag("mergedMva:ntks"), # number of associated Gsf tracks
        "el_tksdr"          : cms.InputTag("mergedMva:tksdr"),
        "el_tksPtRatio"     : cms.InputTag("mergedMva:tksPtRatio"),
        "el_tksRelPtRatio"  : cms.InputTag("mergedMva:tksRelPtRatio"),
        "el_EoverPInv"      : cms.InputTag("mergedMva:EoverPInv"),
        "el_isHggPresel"    : cms.InputTag("mergedMva:isHggPresel")
    }
    for i, j in newVariables.iteritems():
        setattr(tnpVars.CommonStuffForGsfElectronProbe.variables, i, j)
        
    return mergedMva_sequence