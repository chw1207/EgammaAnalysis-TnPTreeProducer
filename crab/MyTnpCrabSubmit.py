import os
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from CRABClient.UserUtilities import config
from httplib import HTTPException
from multiprocessing import Process

# source /cvmfs/cms.cern.ch/common/crab-setup.sh
# python MyTnpCrabSubmit.py

def submit(cfg_):
    try:
        crabCommand("submit", config=cfg_)
    except HTTPException as hte:
        print("\033[91m"+"[ERROR] Failed submitting task: %s" %(hte.headers)+"\033[0m")
    except ClientException as cle:
        print("\033[91m"+"[ERROR] Failed submitting task: %s" %(cle)+"\033[0m")
        

def getLumiMask(era_):
    if era_ == "UL2016preVFP":
        return "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
    elif era_ == "UL2016postVFP":
        return "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
    elif era_ == "UL2017":
        return "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
    elif era_ == "UL2018":
        return "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
    else:
        raise ValueError("Unknown era: {}. [UL2016preVFP, UL2016postVFP, UL2017, UL2018]".format(era_))


def main(requestName_, sample_, era_):
    isMC = "SIM" in sample_
    
    config.General.requestName = requestName_
    config.Data.inputDataset   = sample_
    config.Data.outLFNDirBase  = "{}/{}/{}/".format(mainOutputDir, era_, "mc" if isMC else "data")
    config.Data.splitting      = "FileBased" if isMC else "LumiBased"
    config.Data.lumiMask       = None if isMC else getLumiMask(era_)
    config.Data.unitsPerJob    = 5 if isMC else 25
    config.JobType.pyCfgParams = defaultArgs + ["isMC=True" if isMC else "isMC=False", "era={}".format(era_)]
        
    print(config)
    p = Process(target=submit, args=(config, ))
    p.start()
    p.join()
    print("")
    
    
if __name__ == "__main__":
    # submitVersion = "2023-04-03" # add some date here # -> UnSeeded leg
    submitVersion = "2023-04-11" # add some date here # -> Seeded leg
    # defaultArgs   = ["doEleID=True", "doPhoID=False", "doTrigger=True"] # -> UnSeeded leg
    defaultArgs   = ["doEleID=True", "doPhoID=False", "doTrigger=True", "DiphoLeg=LEAD"] # -> Seeded leg
    mainOutputDir = "/store/user/{}/hdalitz/tnpTuples/{}".format(os.environ["USER"], submitVersion)
    
    config = config()
    config.General.transferLogs            = False
    config.General.workArea                = "crab_{}".format(submitVersion)
    config.JobType.pluginName              = "Analysis"
    config.JobType.psetName                = "../python/TnPTreeProducer_cfg.py"
    config.JobType.sendExternalFolder      = True
    config.JobType.allowUndistributedCMSSW = True
    config.Data.inputDBS                   = "global"
    config.Data.publication                = False
    config.Data.allowNonValidInputDataset  = True
    config.Site.storageSite                = "T2_TW_NCHC"
    
    # UL2016preVFP
    main("UL2016preVFP_Run2016B", "/SingleElectron/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD", "UL2016preVFP")
    main("UL2016preVFP_Run2016C", "/SingleElectron/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD",      "UL2016preVFP")
    main("UL2016preVFP_Run2016D", "/SingleElectron/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD",      "UL2016preVFP")
    main("UL2016preVFP_Run2016E", "/SingleElectron/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD",      "UL2016preVFP")
    main("UL2016preVFP_Run2016F", "/SingleElectron/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD",      "UL2016preVFP")
    main("UL2016preVFP_DY_NLO", "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM", "UL2016preVFP")
    main("UL2016preVFP_DY_LO",  "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",  "UL2016preVFP")
    
    # UL2016postVFP
    main("UL2016postVFP_Run2016F", "/SingleElectron/Run2016F-21Feb2020_UL2016-v1/MINIAOD", "UL2016postVFP")
    main("UL2016postVFP_Run2016G", "/SingleElectron/Run2016G-21Feb2020_UL2016-v1/MINIAOD", "UL2016postVFP")
    main("UL2016postVFP_Run2016H", "/SingleElectron/Run2016H-21Feb2020_UL2016-v2/MINIAOD", "UL2016postVFP")
    main("UL2016postVFP_DY_NLO", "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM", "UL2016postVFP")
    main("UL2016postVFP_DY_LO",  "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM",  "UL2016postVFP")
    
    # UL2017
    main("UL2017_Run2017B", "/SingleElectron/Run2017B-09Aug2019_UL2017-v1/MINIAOD",                 "UL2017")
    main("UL2017_Run2017C", "/SingleElectron/Run2017C-09Aug2019_UL2017-v1/MINIAOD",                 "UL2017")
    main("UL2017_Run2017D", "/SingleElectron/Run2017D-09Aug2019_UL2017-v1/MINIAOD",                 "UL2017")
    main("UL2017_Run2017E", "/SingleElectron/Run2017E-09Aug2019_UL2017-v1/MINIAOD",                 "UL2017")
    main("UL2017_Run2017F", "/SingleElectron/Run2017F-09Aug2019_UL2017_EcalRecovery-v1/MINIAOD",    "UL2017")
    main("UL2017_DY_NLO", "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6-v2/MINIAODSIM", "UL2017")
    main("UL2017_DY_LO",  "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6-v2/MINIAODSIM",  "UL2017")
    
    # UL2018
    main("UL2018_Run2018A", "/EGamma/Run2018A-12Nov2019_UL2018-v2/MINIAOD", "UL2018")
    main("UL2018_Run2018B", "/EGamma/Run2018B-12Nov2019_UL2018-v2/MINIAOD", "UL2018")
    main("UL2018_Run2018C", "/EGamma/Run2018C-12Nov2019_UL2018-v2/MINIAOD", "UL2018")
    main("UL2018_Run2018D", "/EGamma/Run2018D-12Nov2019_UL2018-v4/MINIAOD", "UL2018")
    main("UL2018_DY_NLO", "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM", "UL2018")
    main("UL2018_DY_LO",  "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM",  "UL2018")
    