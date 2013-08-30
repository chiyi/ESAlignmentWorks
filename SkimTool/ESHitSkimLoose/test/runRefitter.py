import FWCore.ParameterSet.Config as cms

process = cms.Process("REFIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
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


process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter.constraint = ""
#process.TrackRefitter.src = "doConstraint"
 
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'rfio:///castor/cern.ch/cms/store/data/Run2011A/Jet/RECO/PromptReco-v4/000/165/180/F4B900FE-1182-E011-AC8E-0030487CD14E.root'
#'rfio:///castor/cern.ch/cms/store/data/Run2011A/MinimumBias/ALCARECO/TkAlMinBias-v1/000/161/008/B237A638-9B55-E011-B210-0030487CD6E8.root'
#'rfio://castor/cern.ch/cms/store/data/Run2011A/MinimumBias/ALCARECO/TkAlMinBias-v1/000/161/312/CE20632D-8B59-E011-9A9C-0030487CD6D8.root'
'file:JET2011A_ESSkim_P1.root'

  )
)

process.out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring('keep *'),
     fileName = cms.untracked.string('ESSkim_wRefitter.root')
)

process.p = cms.EndPath(process.out)


#process.p1 = cms.Path(process.doConstraint * process.TrackRefitter} 
process.p1 = cms.Path(process.TrackRefitter) 

process.schedule = cms.Schedule( process.p1 , process.p)

