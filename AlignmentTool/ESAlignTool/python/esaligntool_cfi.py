import FWCore.ParameterSet.Config as cms

esAlignTool = cms.EDAnalyzer('ESAlignTool',

    DrawMagField = cms.bool(False),
    PrintPosition = cms.bool(False),
    withRotation = cms.bool(False),
    Cal_ESorigin_from_Geometry = cms.bool(True),
    Cal_ESaxes_from_Geometry = cms.bool(True),
    Selected_idee = cms.uint32(0),
    Selected_RUNmin = cms.int32(0),
    Selected_RUNmax = cms.int32(0),
    PrintMatrix = cms.bool(False),
    ReSetRfromOutside = cms.bool(False),
    Overwrite_RotationMatrix_fromGeometry = cms.bool(False),
 
    IterN = cms.uint32(0),

    Iter1_ESpFdX = cms.double(0),Iter1_ESpFdY = cms.double(0),Iter1_ESpFdZ = cms.double(0),Iter1_ESpFdAlpha = cms.double(0),Iter1_ESpFdBeta = cms.double(0),Iter1_ESpFdGamma = cms.double(0),
    Iter1_ESpRdX = cms.double(0),Iter1_ESpRdY = cms.double(0),Iter1_ESpRdZ = cms.double(0),Iter1_ESpRdAlpha = cms.double(0),Iter1_ESpRdBeta = cms.double(0),Iter1_ESpRdGamma = cms.double(0),
    Iter1_ESmFdX = cms.double(0),Iter1_ESmFdY = cms.double(0),Iter1_ESmFdZ = cms.double(0),Iter1_ESmFdAlpha = cms.double(0),Iter1_ESmFdBeta = cms.double(0),Iter1_ESmFdGamma = cms.double(0),
    Iter1_ESmRdX = cms.double(0),Iter1_ESmRdY = cms.double(0),Iter1_ESmRdZ = cms.double(0),Iter1_ESmRdAlpha = cms.double(0),Iter1_ESmRdBeta = cms.double(0),Iter1_ESmRdGamma = cms.double(0),

    Iter2_ESpFdX = cms.double(0),Iter2_ESpFdY = cms.double(0),Iter2_ESpFdZ = cms.double(0),Iter2_ESpFdAlpha = cms.double(0),Iter2_ESpFdBeta = cms.double(0),Iter2_ESpFdGamma = cms.double(0),
    Iter2_ESpRdX = cms.double(0),Iter2_ESpRdY = cms.double(0),Iter2_ESpRdZ = cms.double(0),Iter2_ESpRdAlpha = cms.double(0),Iter2_ESpRdBeta = cms.double(0),Iter2_ESpRdGamma = cms.double(0),
    Iter2_ESmFdX = cms.double(0),Iter2_ESmFdY = cms.double(0),Iter2_ESmFdZ = cms.double(0),Iter2_ESmFdAlpha = cms.double(0),Iter2_ESmFdBeta = cms.double(0),Iter2_ESmFdGamma = cms.double(0),
    Iter2_ESmRdX = cms.double(0),Iter2_ESmRdY = cms.double(0),Iter2_ESmRdZ = cms.double(0),Iter2_ESmRdAlpha = cms.double(0),Iter2_ESmRdBeta = cms.double(0),Iter2_ESmRdGamma = cms.double(0),

    Iter3_ESpFdX = cms.double(0),Iter3_ESpFdY = cms.double(0),Iter3_ESpFdZ = cms.double(0),Iter3_ESpFdAlpha = cms.double(0),Iter3_ESpFdBeta = cms.double(0),Iter3_ESpFdGamma = cms.double(0),
    Iter3_ESpRdX = cms.double(0),Iter3_ESpRdY = cms.double(0),Iter3_ESpRdZ = cms.double(0),Iter3_ESpRdAlpha = cms.double(0),Iter3_ESpRdBeta = cms.double(0),Iter3_ESpRdGamma = cms.double(0),
    Iter3_ESmFdX = cms.double(0),Iter3_ESmFdY = cms.double(0),Iter3_ESmFdZ = cms.double(0),Iter3_ESmFdAlpha = cms.double(0),Iter3_ESmFdBeta = cms.double(0),Iter3_ESmFdGamma = cms.double(0),
    Iter3_ESmRdX = cms.double(0),Iter3_ESmRdY = cms.double(0),Iter3_ESmRdZ = cms.double(0),Iter3_ESmRdAlpha = cms.double(0),Iter3_ESmRdBeta = cms.double(0),Iter3_ESmRdGamma = cms.double(0),

    Iter4_ESpFdX = cms.double(0),Iter4_ESpFdY = cms.double(0),Iter4_ESpFdZ = cms.double(0),Iter4_ESpFdAlpha = cms.double(0),Iter4_ESpFdBeta = cms.double(0),Iter4_ESpFdGamma = cms.double(0),
    Iter4_ESpRdX = cms.double(0),Iter4_ESpRdY = cms.double(0),Iter4_ESpRdZ = cms.double(0),Iter4_ESpRdAlpha = cms.double(0),Iter4_ESpRdBeta = cms.double(0),Iter4_ESpRdGamma = cms.double(0),
    Iter4_ESmFdX = cms.double(0),Iter4_ESmFdY = cms.double(0),Iter4_ESmFdZ = cms.double(0),Iter4_ESmFdAlpha = cms.double(0),Iter4_ESmFdBeta = cms.double(0),Iter4_ESmFdGamma = cms.double(0),
    Iter4_ESmRdX = cms.double(0),Iter4_ESmRdY = cms.double(0),Iter4_ESmRdZ = cms.double(0),Iter4_ESmRdAlpha = cms.double(0),Iter4_ESmRdBeta = cms.double(0),Iter4_ESmRdGamma = cms.double(0),

    Iter5_ESpFdX = cms.double(0),Iter5_ESpFdY = cms.double(0),Iter5_ESpFdZ = cms.double(0),Iter5_ESpFdAlpha = cms.double(0),Iter5_ESpFdBeta = cms.double(0),Iter5_ESpFdGamma = cms.double(0),
    Iter5_ESpRdX = cms.double(0),Iter5_ESpRdY = cms.double(0),Iter5_ESpRdZ = cms.double(0),Iter5_ESpRdAlpha = cms.double(0),Iter5_ESpRdBeta = cms.double(0),Iter5_ESpRdGamma = cms.double(0),
    Iter5_ESmFdX = cms.double(0),Iter5_ESmFdY = cms.double(0),Iter5_ESmFdZ = cms.double(0),Iter5_ESmFdAlpha = cms.double(0),Iter5_ESmFdBeta = cms.double(0),Iter5_ESmFdGamma = cms.double(0),
    Iter5_ESmRdX = cms.double(0),Iter5_ESmRdY = cms.double(0),Iter5_ESmRdZ = cms.double(0),Iter5_ESmRdAlpha = cms.double(0),Iter5_ESmRdBeta = cms.double(0),Iter5_ESmRdGamma = cms.double(0),

    Iter6_ESpFdX = cms.double(0),Iter6_ESpFdY = cms.double(0),Iter6_ESpFdZ = cms.double(0),Iter6_ESpFdAlpha = cms.double(0),Iter6_ESpFdBeta = cms.double(0),Iter6_ESpFdGamma = cms.double(0),
    Iter6_ESpRdX = cms.double(0),Iter6_ESpRdY = cms.double(0),Iter6_ESpRdZ = cms.double(0),Iter6_ESpRdAlpha = cms.double(0),Iter6_ESpRdBeta = cms.double(0),Iter6_ESpRdGamma = cms.double(0),
    Iter6_ESmFdX = cms.double(0),Iter6_ESmFdY = cms.double(0),Iter6_ESmFdZ = cms.double(0),Iter6_ESmFdAlpha = cms.double(0),Iter6_ESmFdBeta = cms.double(0),Iter6_ESmFdGamma = cms.double(0),
    Iter6_ESmRdX = cms.double(0),Iter6_ESmRdY = cms.double(0),Iter6_ESmRdZ = cms.double(0),Iter6_ESmRdAlpha = cms.double(0),Iter6_ESmRdBeta = cms.double(0),Iter6_ESmRdGamma = cms.double(0),

    Iter7_ESpFdX = cms.double(0),Iter7_ESpFdY = cms.double(0),Iter7_ESpFdZ = cms.double(0),Iter7_ESpFdAlpha = cms.double(0),Iter7_ESpFdBeta = cms.double(0),Iter7_ESpFdGamma = cms.double(0),
    Iter7_ESpRdX = cms.double(0),Iter7_ESpRdY = cms.double(0),Iter7_ESpRdZ = cms.double(0),Iter7_ESpRdAlpha = cms.double(0),Iter7_ESpRdBeta = cms.double(0),Iter7_ESpRdGamma = cms.double(0),
    Iter7_ESmFdX = cms.double(0),Iter7_ESmFdY = cms.double(0),Iter7_ESmFdZ = cms.double(0),Iter7_ESmFdAlpha = cms.double(0),Iter7_ESmFdBeta = cms.double(0),Iter7_ESmFdGamma = cms.double(0),
    Iter7_ESmRdX = cms.double(0),Iter7_ESmRdY = cms.double(0),Iter7_ESmRdZ = cms.double(0),Iter7_ESmRdAlpha = cms.double(0),Iter7_ESmRdBeta = cms.double(0),Iter7_ESmRdGamma = cms.double(0),

    Iter8_ESpFdX = cms.double(0),Iter8_ESpFdY = cms.double(0),Iter8_ESpFdZ = cms.double(0),Iter8_ESpFdAlpha = cms.double(0),Iter8_ESpFdBeta = cms.double(0),Iter8_ESpFdGamma = cms.double(0),
    Iter8_ESpRdX = cms.double(0),Iter8_ESpRdY = cms.double(0),Iter8_ESpRdZ = cms.double(0),Iter8_ESpRdAlpha = cms.double(0),Iter8_ESpRdBeta = cms.double(0),Iter8_ESpRdGamma = cms.double(0),
    Iter8_ESmFdX = cms.double(0),Iter8_ESmFdY = cms.double(0),Iter8_ESmFdZ = cms.double(0),Iter8_ESmFdAlpha = cms.double(0),Iter8_ESmFdBeta = cms.double(0),Iter8_ESmFdGamma = cms.double(0),
    Iter8_ESmRdX = cms.double(0),Iter8_ESmRdY = cms.double(0),Iter8_ESmRdZ = cms.double(0),Iter8_ESmRdAlpha = cms.double(0),Iter8_ESmRdBeta = cms.double(0),Iter8_ESmRdGamma = cms.double(0),

    Iter9_ESpFdX = cms.double(0),Iter9_ESpFdY = cms.double(0),Iter9_ESpFdZ = cms.double(0),Iter9_ESpFdAlpha = cms.double(0),Iter9_ESpFdBeta = cms.double(0),Iter9_ESpFdGamma = cms.double(0),
    Iter9_ESpRdX = cms.double(0),Iter9_ESpRdY = cms.double(0),Iter9_ESpRdZ = cms.double(0),Iter9_ESpRdAlpha = cms.double(0),Iter9_ESpRdBeta = cms.double(0),Iter9_ESpRdGamma = cms.double(0),
    Iter9_ESmFdX = cms.double(0),Iter9_ESmFdY = cms.double(0),Iter9_ESmFdZ = cms.double(0),Iter9_ESmFdAlpha = cms.double(0),Iter9_ESmFdBeta = cms.double(0),Iter9_ESmFdGamma = cms.double(0),
    Iter9_ESmRdX = cms.double(0),Iter9_ESmRdY = cms.double(0),Iter9_ESmRdZ = cms.double(0),Iter9_ESmRdAlpha = cms.double(0),Iter9_ESmRdBeta = cms.double(0),Iter9_ESmRdGamma = cms.double(0),

    Iter10_ESpFdX = cms.double(0),Iter10_ESpFdY = cms.double(0),Iter10_ESpFdZ = cms.double(0),Iter10_ESpFdAlpha = cms.double(0),Iter10_ESpFdBeta = cms.double(0),Iter10_ESpFdGamma = cms.double(0),
    Iter10_ESpRdX = cms.double(0),Iter10_ESpRdY = cms.double(0),Iter10_ESpRdZ = cms.double(0),Iter10_ESpRdAlpha = cms.double(0),Iter10_ESpRdBeta = cms.double(0),Iter10_ESpRdGamma = cms.double(0),
    Iter10_ESmFdX = cms.double(0),Iter10_ESmFdY = cms.double(0),Iter10_ESmFdZ = cms.double(0),Iter10_ESmFdAlpha = cms.double(0),Iter10_ESmFdBeta = cms.double(0),Iter10_ESmFdGamma = cms.double(0),
    Iter10_ESmRdX = cms.double(0),Iter10_ESmRdY = cms.double(0),Iter10_ESmRdZ = cms.double(0),Iter10_ESmRdAlpha = cms.double(0),Iter10_ESmRdBeta = cms.double(0),Iter10_ESmRdGamma = cms.double(0)


)
