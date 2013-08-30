import FWCore.ParameterSet.Config as cms

process = cms.Process("REFIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'FT_R_53_V6::All'
#process.GlobalTag.globaltag = 'GR_P_V42::All'
process.GlobalTag.globaltag = 'GR_R_53_V18::All'
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)



#process.newTKAlignment = cms.ESSource("PoolDBESSource",
#                                        process.CondDBSetup,
#                                        connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT'),
#                                        timetype = cms.string('runnumber'),
#                                        toGet = cms.VPSet(cms.PSet(
#                                                record = cms.string('TrackerAlignmentRcd'),
#                                                tag = cms.string('TrackerAlignment_v9a_offline')
#                                                ))
#                                        )
#process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "newTKAlignment")

#process.newTKAlignmentError = cms.ESSource("PoolDBESSource",
#                                        process.CondDBSetup,
#                                        connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ALIGNMENT'),
#                                        timetype = cms.string('runnumber'),
#                                        toGet = cms.VPSet(cms.PSet(
#                                                record = cms.string('TrackerAlignmentErrorRcd'),
#                                                tag = cms.string('TrackerAlignmentErrors_v3_offline')
#                                                ))
#                                        )
#process.es_prefer_trackerAlignmentError = cms.ESPrefer("PoolDBESSource", "newTKAlignmentError")
#
#
#process.newTrackerSurfaceDeformation = cms.ESSource("PoolDBESSource",
#                                          process.CondDBSetup,
#                                          connect = cms.string('frontier://FrontierProd/CMS_COND_310X_ALIGN'),
#                                          toGet= cms.VPSet(cms.PSet(record = cms.string("TrackerSurfaceDeformationRcd"),
#                                                                     tag = cms.string('TrackerSurfaceDeformations_v3_offline'))
#                                                           )
#                                         )
#process.es_prefer_wTrackerSurfaceDeformation = cms.ESPrefer("PoolDBESSource", "newTrackerSurfaceDeformation")


#process.newTKAlignment = cms.ESSource("PoolDBESSource",
#                                        process.CondDBSetup,
#                                        connect = cms.string('frontier://FrontierPrep/CMS_COND_ALIGNMENT'),
#                                        timetype = cms.string('runnumber'),
#                                        toGet = cms.VPSet(cms.PSet(
#                                                record = cms.string('TrackerAlignmentRcd'),
#                                                tag = cms.string('TrackerAlignment_Run2012D_offline')
#                                                ))
#                                        )
#process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "newTKAlignment")


process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter.constraint = ""
#process.TrackRefitter.src = "doConstraint"

 
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

 'root://eoscms//eos/cms/store/group/alca_ecalcalib/ESalignment/5_3_2_p4/small_Run2012D_bkup/00DE8F3C-890B-E211-8285-001D09F2AD4D.root',

  )
)

#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_190456-196531_8TeV_29Jun2012ReReco_Collisions12_JSON.txt').getVLuminosityBlockRange()



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
       'keep *_TrackRefitter_*_*',
       'keep *_offlineBeamSpot_*_*',
       'keep *_siPixelClusters_*_*',
       'keep *_siStripClusters_*_*'
      ),
     fileName = cms.untracked.string('file:ESSkim_wRefit.root')
)
 
process.p = cms.EndPath(process.out)
process.p1 = cms.Path(process.TrackRefitter)

process.schedule = cms.Schedule(process.myPath,process.p)

