import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticleCandidates_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')


#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
'file:/cmsrm/pc18/pandolf/CMSSW_3_7_0_patch2/src/JetMETCorrections/GammaJet/test/PD.root'
)

)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
#    #wantSummary = cms.untracked.bool(True)
#)

process.MessageLogger.cerr.FwkReport.reportEvery = 10


#############   Include the jet corrections ##########
#from JetMETCorrections.Configuration.JetCorrectionEra_cff import *
#JetCorrectionEra.era = 'Summer09_7TeV_ReReco332' # applies to L2 & L3 only
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")


#monster track event cleaning
process.monster = cms.EDFilter(
   "FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.2)
)





process.myanalysis = cms.EDAnalyzer("HZZlljjAnalyzer",
    debug = cms.bool(False),
    recoProducer = cms.string('ecalRecHit'),
    MCTruthCollection = cms.untracked.InputTag("source"),
    genMet = cms.untracked.InputTag("genMetTrue"),
    met = cms.untracked.InputTag("met"),
    tracks = cms.untracked.InputTag("generalTracks"),
    Electronsrc = cms.untracked.InputTag("gsfelectrons"),
    Muonsrc = cms.untracked.InputTag("muons"),
    recoCollection = cms.string('EcalRecHitsEB'),
    JetAlgo = cms.string('ak5'),
    hbhits = cms.untracked.InputTag("hbhereco"),
    TriggerTag = cms.untracked.InputTag("TriggerResults::HLT"),
    vertices = cms.untracked.InputTag("offlinePrimaryVertices"),
    genjetptthr = cms.double(5.),
    jetptthr = cms.double(4.),
    genjetnmin = cms.int32(10),
    jetnmin = cms.int32(10),
)


process.p = cms.Path(process.monster*process.myanalysis)
