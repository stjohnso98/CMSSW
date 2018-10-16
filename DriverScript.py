import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
process = cms.Process("Demo")

#dataset names
DYMulist = FileUtils.loadListFromFile('/afs/cern.ch/work/s/stjohnso/opend/opend/CMSSW_5_3_32/src/calotowerFt/DemoAnalyzer/CMS_MonteCarlo2012_Summer12_DR53X_DYToMuMu_M-20_CT10_8TeV-powheg-pythia6_AODSIM_PU_S10_START53_V19-v1_00000_file_index.txt')
readFilesDYMu = cms.untracked.vstring(*DYMulist)
DYElelist = FileUtils.loadListFromFile('/afs/cern.ch/work/s/stjohnso/opend/opend/CMSSW_5_3_32/src/calotowerFt/DemoAnalyzer/CMS_MonteCarlo2012_Summer12_DR53X_DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_AODSIM_NewG4Phys_PU_RD1_START53_V7N-v1_00000_file_index.txt')
readFilesDYEle = cms.untracked.vstring(*DYElelist)
DYQCDlist = FileUtils.loadListFromFile('/afs/cern.ch/work/s/stjohnso/opend/opend/CMSSW_5_3_32/src/calotowerFt/DemoAnalyzer/CMS_MonteCarlo2012_Summer12_DR53X_TTbar_8TeV-Madspin_aMCatNLO-herwig_AODSIM_PU_S10_START53_V19-v2_00000_file_index.txt')
readFilesDYQCD = cms.untracked.vstring(*DYQCDlist)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000000) )

process.source = cms.Source("PoolSource",fileNames = readFilesDYQCD)
#Creating Tree
process.ExportTree = cms.EDAnalyzer("RecoTree"

)
#Calling C++ Script
process.demo = cms.EDAnalyzer('ElectronTree',

                                         minTracks=cms.untracked.uint32(0),
                              MC=cms.bool(True)
                              
                              )
#Creating Output File
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('TTBar.root')
                                   )


process.p = cms.Path(process.demo)
