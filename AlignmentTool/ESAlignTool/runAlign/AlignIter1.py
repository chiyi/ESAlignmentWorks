import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

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


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    firstEvent = cms.untracked.uint32(1),
	fileNames = cms.untracked.vstring(

'file:test.root'

	)

)

process.load("AlignmentTool.ESAlignTool.esaligntool_cfi")
#process.esAlignTool.DrawMagField = cms.bool(True)
#process.esAlignTool.PrintPosition = cms.bool(True)
process.esAlignTool.Overwrite_RotationMatrix_fromGeometry = cms.bool(True)
process.esAlignTool.fromRefitter = cms.bool(True)

process.esAlignTool.Cal_ESorigin_from_Geometry = cms.bool(True)
process.esAlignTool.withRotation = cms.bool(True)
process.esAlignTool.IterN = cms.uint32(1)
process.esAlignTool.e_xxlimit = cms.double(1.)
process.esAlignTool.e_yylimit = cms.double(1.)
process.esAlignTool.e_yxlimit = cms.double(1.)
process.esAlignTool.winlimit = cms.double(3.)

process.p = cms.Path(process.esAlignTool)

