#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "EgammaAnalysis/TnPTreeProducer/plugins/WriteValueMap.h"
#include "EgammaAnalysis/TnPTreeProducer/plugins/isolations.h"

#include "TMath.h"
#include "TLorentzVector.h"
// #include <xgboost/c_api.h>

// #include "XGBoostCMSSW/XGBoostInterface/interface/XGBReader.h"

/*
    #define safe_xgboost(call) { \
        int err = (call);\
        if (err != 0) {\
            cms::Exception exception("XGBoost error");\
            exception << __FILE__ << ":" << __LINE__ << " in " << __FUNCTION__ \
            << " :" << XGBGetLastError();\
            throw exception;\
        } \
    }
*/

class MergedElectronMvaProducer: public edm::EDProducer{
public:
    explicit MergedElectronMvaProducer(const edm::ParameterSet &iConfig);
    virtual ~MergedElectronMvaProducer();

    virtual void beginJob();
    virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override;

private:
    // std::string weightFileNameEB_;
    // std::string weightFileNameEE_;
    edm::EDGetTokenT<std::vector<pat::Electron>> probesToken_;
    edm::EDGetTokenT<edm::View<reco::GsfTrack>> gsfTracksToken_;
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<edm::View<pat::Photon>> photonsToken_;
    edm::EDGetTokenT<double> rhoAllToken_;

    // std::vector<XGBReader*> readers = {NULL, NULL};
};


MergedElectronMvaProducer::MergedElectronMvaProducer(const edm::ParameterSet &iConfig):
    // weightFileNameEB_(iConfig.getParameter<edm::FileInPath>("weightFileEB").fullPath()),
    // weightFileNameEE_(iConfig.getParameter<edm::FileInPath>("weightFileEE").fullPath()),
    probesToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("probes"))),
    gsfTracksToken_(consumes<edm::View<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("gsfTracks"))),
    rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInputTag"))),
    photonsToken_(consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
    rhoAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAllInputTag"))){
    
    // merged electron ID MVA
    // produces<edm::ValueMap<float>>("mergedMVA");
    
    // addtional gsf track related variables and egm post calibration 
    produces<edm::ValueMap<float>>("ntks");
    produces<edm::ValueMap<float>>("tksdr");
    produces<edm::ValueMap<float>>("tksPtRatio");
    produces<edm::ValueMap<float>>("tksRelPtRatio");
    produces<edm::ValueMap<float>>("pterr");
    produces<edm::ValueMap<float>>("EoverPInv");
    produces<edm::ValueMap<float>>("isHggPresel");
}


MergedElectronMvaProducer::~MergedElectronMvaProducer(){
    // for (int i = 0; i < 2; i++){
    //     // delete readers[i];
    //     safe_xgboost(XGBoosterFree(booster));
    // }
    // const char* dummy = weightFileNameEB_.c_str();
}


void MergedElectronMvaProducer::beginJob(){
    // const char* fName = weightFileNameEB_.c_str();
    // safe_xgboost(XGBoosterCreate(NULL, 0, &booster));
    //     // safe_xgboost(XGBoosterSetParam(booster, "nthread", "1"));
    // safe_xgboost(XGBoosterLoadModel(booster, fName));
    
    // for (int i = 0; i < 2; i++){
    //     const char* fName = (i == 0) ? weightFileNameEB_.c_str() : weightFileNameEE_.c_str();
    //     std::cout << "Reading weight file: " << fName << std::endl;
    //     safe_xgboost(XGBoosterCreate(NULL, 0, &boosters[i]));
    //     // safe_xgboost(XGBoosterSetParam(boosters[i], "nthread", "1"));
    //     safe_xgboost(XGBoosterLoadModel(boosters[i], fName));
    // }
}


void MergedElectronMvaProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup){
    edm::Handle<std::vector<pat::Electron>> probes;
    iEvent.getByToken(probesToken_, probes);

    edm::Handle<edm::View<reco::GsfTrack>> gsfTracks;
    iEvent.getByToken(gsfTracksToken_, gsfTracks);

    edm::Handle<double> rhoH;
    iEvent.getByToken(rhoToken_, rhoH);

    edm::Handle<edm::View<pat::Photon>> photons;
    iEvent.getByToken(photonsToken_, photons);

    edm::Handle<double> rhoAll;
    iEvent.getByToken(rhoAllToken_, rhoAll);

    std::vector<float> isHggPreselVals;  

    // BoosterHandle booster;
    // safe_xgboost(XGBoosterCreate(NULL, 0, &booster));
    //     // safe_xgboost(XGBoosterSetParam(booster, "nthread", "1"));
    // safe_xgboost(XGBoosterLoadModel(booster, weightFileNameEB_.c_str()));

    // Hgg pre-selection
    // https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Taggers/python/flashggPreselectedDiPhotons_cfi.py
    for (auto probe = probes->begin(); probe != probes->end(); ++probe){
        // find the photon which shares the same sc as electron
        int phoIdx = -1;
        for (size_t idx = 0; idx < photons->size(); ++idx){
            auto photon = photons->ptrAt(idx);
            float ele_sc_eta = probe->superCluster()->eta();
            float pho_sc_eta = photon->superCluster()->eta();
            if (ele_sc_eta == pho_sc_eta){
                phoIdx = idx;
                break;
            }
        }

        float isHgg = 0;
        if (phoIdx != -1){
            // matched photon should pass Hgg pre-selection
            auto photon_matched = photons->ptrAt(phoIdx);

            float phoSCEta = photon_matched->superCluster()->eta();
            float phoEffArea = 0.13212;
            if (fabs(phoSCEta) > 0 && fabs(phoSCEta) < 1.5) 
                phoEffArea = 0.16544;

            float rhocorr = *(rhoAll.product()) * phoEffArea;
            float phoPFPhoIso_corr = TMath::Max(photon_matched->photonIso() - rhocorr, (float) 0.);

            bool isEB = (fabs(phoSCEta) < 1.4442);
            bool isEE = (fabs(phoSCEta) > 1.566 && fabs(phoSCEta) < 2.5);

            float phoR9Full5x5 = photon_matched->full5x5_r9();
            float phoTrkIsoHollowConeDR03 = photon_matched->trkSumPtHollowConeDR03();
            float phoSigmaIEtaIEtaFull5x5 = photon_matched->full5x5_sigmaIetaIeta();
            bool isHR9_EB = isEB && phoR9Full5x5 >  0.85;
            bool isLR9_EB = isEB && phoR9Full5x5 <= 0.85 && phoR9Full5x5 > 0.5 && phoPFPhoIso_corr < 4 && phoTrkIsoHollowConeDR03 < 6 && phoSigmaIEtaIEtaFull5x5 < 0.015;
            bool isHR9_EE = isEE && phoR9Full5x5 >  0.9;
            bool isLR9_EE = isEE && phoR9Full5x5 <= 0.9  && phoR9Full5x5 > 0.8 && phoPFPhoIso_corr < 4 && phoTrkIsoHollowConeDR03 < 6 && phoSigmaIEtaIEtaFull5x5 < 0.035;

            // cuts here mimic the miniAOD photon cuts and the non-category based trigger cuts
            bool isAOD = photon_matched->hadTowOverEm() < 0.08 && (phoR9Full5x5 > 0.8 || photon_matched->chargedHadronIso() < 20. || (photon_matched->chargedHadronIso()/photon_matched->et() < 0.3));

            // Hgg preselection without passElectronVeto
            isHgg = (isAOD && (isHR9_EB || isLR9_EB || isHR9_EE || isLR9_EE)) ? 1. : 0.;
        }
        isHggPreselVals.push_back(isHgg);
        // std::cout << isHgg << std::endl;
    }

    // a tricky way to find the ambiguousTracks for the pat::Electron from gsftracks
    // first of all, check whether the gsf track is the main gsf track of a electron 
    std::vector<int> isMainGSF;  
    for (auto ig = gsfTracks->begin(); ig != gsfTracks->end(); ++ig){
        int ismain = 0;
        for (auto probe = probes->begin(); probe != probes->end(); ++probe){
            float main_tk_eta = probe->gsfTrack()->eta();
            float main_tk_phi = probe->gsfTrack()->phi();
            float tk_eta = ig->eta();
            float tk_phi = ig->phi();
            if ((tk_eta == main_tk_eta) && (tk_phi == main_tk_phi)){
                ismain = 1;
                break;
            }
        }
        isMainGSF.push_back(ismain);
    }  
      
    // then extract the idex of the main gsf tracks in gsfTracks
    std::vector<int> mainGsfIdx;
    for (size_t i = 0; i < isMainGSF.size(); ++i){
        if (isMainGSF[i] == 1)
            mainGsfIdx.push_back(i);
    } 
    if (probes->size() != mainGsfIdx.size()){ // this should not happen
        cms::Exception exception("# of matched main Gsf tracks != # of electrons");
        exception << __FILE__ << ":" << __LINE__ << " in " << __FUNCTION__
                  << ": Number of electrons: " << probes->size() << ", Number of matched Gsf tracks: " << mainGsfIdx.size();
        throw exception;
    }
    
    // loop over the electron and add the associated tracks
    std::vector<float> ntksVals;
    std::vector<float> tksdrVals;
    std::vector<float> tksPtRatioVals;
    std::vector<float> tksRelPtRatioVals;
    std::vector<float> pterrVals;
    std::vector<float> EoverPInvVals;
    for (auto probe = probes->begin(); probe != probes->end(); ++probe){
        std::pair<int, int> results(-1, -1);
        for (size_t j = 0; j < mainGsfIdx.size(); j++){
            auto gsfTrack = gsfTracks->ptrAt(mainGsfIdx[j]);
            float main_tk_eta = probe->gsfTrack()->eta();
            float main_tk_phi = probe->gsfTrack()->phi();
            float tk_eta = gsfTrack->eta();
            float tk_phi = gsfTrack->phi();
            if ((tk_eta == main_tk_eta) && (tk_phi == main_tk_phi)){
                results.first = j; // mainGsfIdx index number
                results.second = mainGsfIdx[j]; // gsf track index number 
            }
        }

        // start from the main track of an electron, 
        // and keep picking the rest track as the ambiguous track for the electron until meeting the other main track
        std::vector<int> associatedGsf; // the indecies of gsf tracks associated with the electron.
        if (results.second == mainGsfIdx.back()){ // the last main gsf track in this event
            for (size_t k = results.second; k < isMainGSF.size(); ++k){
                associatedGsf.push_back(k);
            }
        }
        else{
            for (int k = results.second; k < mainGsfIdx[results.first+1]; ++k){
                associatedGsf.push_back(k);
            }
        }

        int nAssTrks = associatedGsf.size();
        ntksVals.push_back(nAssTrks);
    
        auto main_tk = gsfTracks->ptrAt(associatedGsf[0]);
        TLorentzVector main_tk_v;
        main_tk_v.SetPtEtaPhiM(main_tk->pt(), main_tk->eta(), main_tk->phi(), 0.000511); // assum it is electron 
        if (nAssTrks > 1){
            // find the cloest tk to main tk to be the sub tk
            float tmp_dr = 999.;
            int el_subtk_idx = -1;
            for (size_t iat = 1; iat < associatedGsf.size(); ++iat){ // 0 is main tk
                auto tmp_sub_tk = gsfTracks->ptrAt(associatedGsf[iat]);
                float dr = deltaR(main_tk->eta(), main_tk->phi(), tmp_sub_tk->eta(), tmp_sub_tk->phi());
                if (dr < tmp_dr){
                tmp_dr = dr;
                el_subtk_idx = associatedGsf[iat];
                }
            }

            auto sub_tk = gsfTracks->ptrAt(el_subtk_idx);
            TLorentzVector sub_tk_v;
            sub_tk_v.SetPtEtaPhiM(sub_tk->pt(), sub_tk->eta(), sub_tk->phi(), 0.000511); // assume it is electron 
            
            tksdrVals.push_back(deltaR(main_tk->eta(), main_tk->phi(), sub_tk->eta(), sub_tk->phi()));
            tksPtRatioVals.push_back(sub_tk->pt()/main_tk->pt());
            tksRelPtRatioVals.push_back((sub_tk_v+main_tk_v).Pt()/probe->superCluster()->rawEnergy());
        }
        else{
            tksdrVals.push_back(-999.);
            tksPtRatioVals.push_back(-999.);
            tksRelPtRatioVals.push_back(main_tk->pt()/probe->superCluster()->rawEnergy());
        }

        pterrVals.push_back(probe->userFloat("ecalTrkEnergyErrPostCorr")*probe->pt()/probe->p());

        if (probe->ecalEnergy() == 0)   
            EoverPInvVals.push_back(1e30);
        else if (!TMath::Finite(probe->ecalEnergy()))  
            EoverPInvVals.push_back(1e30);
        else  
            EoverPInvVals.push_back((1.0 - probe->eSuperClusterOverP())/probe->ecalEnergy());
    }
      
    // add the merged electron ID MVA
    // std::vector<float> mergedMVAVals;
    // for (size_t i = 0; i < probes->size(); i++){
    //     float mva_score = -999.;
    //     if (ntksVals.at(i) > 1){
    //         pat::Electron probe = probes->at(i);

    //         int iBE = (probe.superCluster()->eta() < 1.479) ? 0 : 1;
    //         reco::GsfElectron::PflowIsolationVariables pfIso = probe.pfIsolationVariables();
    //         std::vector<float> features = {
    //             (float) *(rhoH.product()),
    //             (float) probe.superCluster()->eta(),
    //             (float) probe.superCluster()->rawEnergy(),
    //             (float) probe.deltaEtaSuperClusterTrackAtVtx(),
    //             (float) probe.deltaPhiSuperClusterTrackAtVtx(),
    //             (float) pterrVals.at(i),
    //             (float) probe.hcalOverEcal(),
    //             (float) probe.eSuperClusterOverP(),
    //             (float) probe.eEleClusterOverPout(),
    //             (float) EoverPInvVals.at(i),
    //             (float) probe.superCluster()->etaWidth(),
    //             (float) probe.superCluster()->phiWidth(),
    //             (float) probe.full5x5_sigmaIetaIeta(),
    //             (float) probe.full5x5_sigmaIphiIphi(),
    //             (float) probe.full5x5_r9(),
    //             (float) probe.fbrem(),
    //             (float) pfIso.sumChargedHadronPt,
    //             (float) pfIso.sumPhotonEt,
    //             (float) pfIso.sumNeutralHadronEt,
    //             (float) tksPtRatioVals.at(i),
    //             (float) tksdrVals.at(i),
    //             (float) tksRelPtRatioVals.at(i)
    //         };
    //         for (size_t jf = 0; jf < features.size(); jf++){
    //             float var_tmp = features[jf]; // Nan value check
    //             features[jf] = (TMath::Finite(var_tmp)) ? (float) var_tmp : -99999.;
    //         }

    //         // // prediction
    //         // DMatrixHandle dpred;
    //         // bst_ulong out_shape = 0;  
    //         // float const *out_results;  
    //         // safe_xgboost(XGDMatrixCreateFromMat((float*)features.data(), 1, features.size(), -99999., &dpred));
    //         // safe_xgboost(XGBoosterPredict(booster, dpred, 0, 0, &out_shape, &out_results));

    //         // // std::vector<float> scores = readers[iBE]->Compute(features);
    //         // mva_score = out_results[0];
    //         // std::cout << mva_score << std::endl;
    //         // safe_xgboost(XGDMatrixFree(dpred));
    //     }
    //     // std::cout << (1.0 - probe.eSuperClusterOverP())/probe.ecalEnergy() << std::endl;
    //     mergedMVAVals.push_back(mva_score);
    // }
    // safe_xgboost(XGBoosterFree(booster));

    // writeValueMap(iEvent, probes, mergedMVAVals,        "mergedMVA");
    writeValueMap(iEvent, probes, ntksVals,             "ntks");
    writeValueMap(iEvent, probes, tksdrVals,            "tksdr");
    writeValueMap(iEvent, probes, tksPtRatioVals,       "tksPtRatio");
    writeValueMap(iEvent, probes, tksRelPtRatioVals,    "tksRelPtRatio");
    writeValueMap(iEvent, probes, pterrVals,            "pterr");
    writeValueMap(iEvent, probes, EoverPInvVals,        "EoverPInv");
    writeValueMap(iEvent, probes, isHggPreselVals,      "isHggPresel");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MergedElectronMvaProducer);
