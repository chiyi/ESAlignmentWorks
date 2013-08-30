import FWCore.ParameterSet.Config as cms

process = cms.Process("ESSKIM")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'GR_R_311_V1::All'
process.GlobalTag.globaltag = 'GR_P_V14::All'
process.newTKAlignment = cms.ESSource("PoolDBESSource",
                                        process.CondDBSetup,
                                        connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT'),
                                        timetype = cms.string('runnumber'),
                                        toGet = cms.VPSet(cms.PSet(
                                                record = cms.string('TrackerAlignmentRcd'),
                                                tag = cms.string('TrackerAlignment_GR10_v4_offline')
                                                ))
                                        )
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "newTKAlignment")

process.newGlobalPosition = cms.ESSource("PoolDBESSource",
                                          process.CondDBSetup,
                                          connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT'),
                                          toGet= cms.VPSet(cms.PSet(record = cms.string("GlobalPositionRcd"),
                                                                     tag = cms.string('GlobalAlignment_TkRotMuonFromLastIovV2_offline'))
                                                           )
                                         )
process.es_prefer_GlobalPositionDB = cms.ESPrefer("PoolDBESSource", "newGlobalPosition")

 
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

#'rfio://castorcms/?svcClass=cmscafuser&path=/castor/cern.ch/cms/store/data/Run2011A/Jet/RECO/PromptReco-v1/000/161/078/C8206724-0856-E011-A8DD-0030487CD906.root'
'file:/data/chiyi/ESAlignment/CMSSW_4_1_2/src/DATAFiles/JetPD2011A_PromptReco-v1/DAACC74C-DF57-E011-A45D-001D09F2B30B.root'

  )
)

process.myFilter = cms.EDFilter('ESHitSkimLoose')

process.myfilter = cms.Sequence(process.myFilter)

process.myPath = cms.Path(process.myfilter)

process.mySelection = cms.PSet(
 SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('myPath')
 )
)


process.out = cms.OutputModule("PoolOutputModule",
     process.mySelection,
     outputCommands = cms.untracked.vstring('drop *',
       'keep *_ecalPreshowerRecHit_*_*',
       'keep *_generalTracks_*_*',
       'keep *_offlineBeamSpot_*_*',
       'keep *_siPixelClusters_*_*',
       'keep *_siStripClusters_*_*',
      ),
     fileName = cms.untracked.string('JET2011A_ESSkim_P1.root')
)
 
process.p = cms.EndPath(process.out)


process.schedule = cms.Schedule(process.myPath,process.p)
