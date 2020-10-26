import FWCore.ParameterSet.Config as cms

process = cms.Process("Zprime2muAnalysis")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('dummyFile')
)
process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(
        0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
        0.058939, 0.125497
    ),
    HFdepthOneParameterB = cms.vdouble(
        -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
        0.000425, 0.000209
    ),
    HFdepthTwoParameterA = cms.vdouble(
        0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
        0.051579, 0.086593
    ),
    HFdepthTwoParameterB = cms.vdouble(
        -2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
        0.000157, -3e-06
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.mvaEleID_Fall17_iso_V1_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'pt < 10. && abs(superCluster.eta) < 0.800', 
        'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt < 10. && abs(superCluster.eta) >= 1.479', 
        'pt >= 10. && abs(superCluster.eta) < 0.800', 
        'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt >= 10. && abs(superCluster.eta) >= 1.479'
    ),
    mvaName = cms.string('ElectronMVAEstimatorRun2'),
    mvaTag = cms.string('Fall17IsoV1'),
    nCategories = cms.int32(6),
    variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Fall17V1Variables.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml.gz'
    )
)

process.mvaEleID_Fall17_iso_V2_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'pt < 10. && abs(superCluster.eta) < 0.800', 
        'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt < 10. && abs(superCluster.eta) >= 1.479', 
        'pt >= 10. && abs(superCluster.eta) < 0.800', 
        'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt >= 10. && abs(superCluster.eta) >= 1.479'
    ),
    mvaName = cms.string('ElectronMVAEstimatorRun2'),
    mvaTag = cms.string('Fall17IsoV2'),
    nCategories = cms.int32(6),
    variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB1_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB2_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EE_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB1_10.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB2_10.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EE_10.weights.xml.gz'
    )
)

process.mvaEleID_Fall17_noIso_V1_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'pt < 10. && abs(superCluster.eta) < 0.800', 
        'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt < 10. && abs(superCluster.eta) >= 1.479', 
        'pt >= 10. && abs(superCluster.eta) < 0.800', 
        'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt >= 10. && abs(superCluster.eta) >= 1.479'
    ),
    mvaName = cms.string('ElectronMVAEstimatorRun2'),
    mvaTag = cms.string('Fall17NoIsoV1'),
    nCategories = cms.int32(6),
    variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Fall17V1Variables.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml.gz'
    )
)

process.mvaEleID_Fall17_noIso_V2_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'pt < 10. && abs(superCluster.eta) < 0.800', 
        'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt < 10. && abs(superCluster.eta) >= 1.479', 
        'pt >= 10. && abs(superCluster.eta) < 0.800', 
        'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt >= 10. && abs(superCluster.eta) >= 1.479'
    ),
    mvaName = cms.string('ElectronMVAEstimatorRun2'),
    mvaTag = cms.string('Fall17NoIsoV2'),
    nCategories = cms.int32(6),
    variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB1_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB2_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EE_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB1_10.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB2_10.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EE_10.weights.xml.gz'
    )
)

process.mvaEleID_Spring16_GeneralPurpose_V1_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'abs(superCluster.eta) < 0.800', 
        'abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'abs(superCluster.eta) >= 1.479'
    ),
    mvaName = cms.string('ElectronMVAEstimatorRun2'),
    mvaTag = cms.string('Spring16GeneralPurposeV1'),
    nCategories = cms.int32(3),
    variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml.gz'
    )
)

process.mvaEleID_Spring16_HZZ_V1_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'pt < 10. && abs(superCluster.eta) < 0.800', 
        'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt < 10. && abs(superCluster.eta) >= 1.479', 
        'pt >= 10. && abs(superCluster.eta) < 0.800', 
        'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
        'pt >= 10. && abs(superCluster.eta) >= 1.479'
    ),
    mvaName = cms.string('ElectronMVAEstimatorRun2'),
    mvaTag = cms.string('Spring16HZZV1'),
    nCategories = cms.int32(6),
    variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml.gz', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml.gz'
    )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.trkIsol03CfgV2 = cms.PSet(
    barrelCuts = cms.PSet(
        algosToReject = cms.vstring(),
        allowedQualities = cms.vstring(),
        maxDPtPt = cms.double(0.1),
        maxDR = cms.double(0.3),
        maxDZ = cms.double(0.1),
        minDEta = cms.double(0.005),
        minDR = cms.double(0.0),
        minHits = cms.int32(8),
        minPixelHits = cms.int32(1),
        minPt = cms.double(1.0)
    ),
    endcapCuts = cms.PSet(
        algosToReject = cms.vstring(),
        allowedQualities = cms.vstring(),
        maxDPtPt = cms.double(0.1),
        maxDR = cms.double(0.3),
        maxDZ = cms.double(0.5),
        minDEta = cms.double(0.005),
        minDR = cms.double(0.0),
        minHits = cms.int32(8),
        minPixelHits = cms.int32(1),
        minPt = cms.double(1.0)
    )
)

process.mvaConfigsForEleProducer = cms.VPSet(
    cms.PSet(
        categoryCuts = cms.vstring(
            'pt < 10. && abs(superCluster.eta) < 0.800', 
            'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt < 10. && abs(superCluster.eta) >= 1.479', 
            'pt >= 10. && abs(superCluster.eta) < 0.800', 
            'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt >= 10. && abs(superCluster.eta) >= 1.479'
        ),
        mvaName = cms.string('ElectronMVAEstimatorRun2'),
        mvaTag = cms.string('Spring16HZZV1'),
        nCategories = cms.int32(6),
        variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml.gz'
        )
    ), 
    cms.PSet(
        categoryCuts = cms.vstring(
            'abs(superCluster.eta) < 0.800', 
            'abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'abs(superCluster.eta) >= 1.479'
        ),
        mvaName = cms.string('ElectronMVAEstimatorRun2'),
        mvaTag = cms.string('Spring16GeneralPurposeV1'),
        nCategories = cms.int32(3),
        variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml.gz'
        )
    ), 
    cms.PSet(
        categoryCuts = cms.vstring(
            'pt < 10. && abs(superCluster.eta) < 0.800', 
            'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt < 10. && abs(superCluster.eta) >= 1.479', 
            'pt >= 10. && abs(superCluster.eta) < 0.800', 
            'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt >= 10. && abs(superCluster.eta) >= 1.479'
        ),
        mvaName = cms.string('ElectronMVAEstimatorRun2'),
        mvaTag = cms.string('Fall17NoIsoV1'),
        nCategories = cms.int32(6),
        variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Fall17V1Variables.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml.gz'
        )
    ), 
    cms.PSet(
        categoryCuts = cms.vstring(
            'pt < 10. && abs(superCluster.eta) < 0.800', 
            'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt < 10. && abs(superCluster.eta) >= 1.479', 
            'pt >= 10. && abs(superCluster.eta) < 0.800', 
            'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt >= 10. && abs(superCluster.eta) >= 1.479'
        ),
        mvaName = cms.string('ElectronMVAEstimatorRun2'),
        mvaTag = cms.string('Fall17IsoV1'),
        nCategories = cms.int32(6),
        variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Fall17V1Variables.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml.gz'
        )
    ), 
    cms.PSet(
        categoryCuts = cms.vstring(
            'pt < 10. && abs(superCluster.eta) < 0.800', 
            'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt < 10. && abs(superCluster.eta) >= 1.479', 
            'pt >= 10. && abs(superCluster.eta) < 0.800', 
            'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt >= 10. && abs(superCluster.eta) >= 1.479'
        ),
        mvaName = cms.string('ElectronMVAEstimatorRun2'),
        mvaTag = cms.string('Fall17NoIsoV2'),
        nCategories = cms.int32(6),
        variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB1_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB2_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EE_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB1_10.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB2_10.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EE_10.weights.xml.gz'
        )
    ), 
    cms.PSet(
        categoryCuts = cms.vstring(
            'pt < 10. && abs(superCluster.eta) < 0.800', 
            'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt < 10. && abs(superCluster.eta) >= 1.479', 
            'pt >= 10. && abs(superCluster.eta) < 0.800', 
            'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
            'pt >= 10. && abs(superCluster.eta) >= 1.479'
        ),
        mvaName = cms.string('ElectronMVAEstimatorRun2'),
        mvaTag = cms.string('Fall17IsoV2'),
        nCategories = cms.int32(6),
        variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB1_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB2_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EE_5.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB1_10.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB2_10.weights.xml.gz', 
            'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EE_10.weights.xml.gz'
        )
    )
)

process.Our2017Leptons = cms.EDProducer("Zprime2muLeptonProducer_miniAOD",
    bits = cms.InputTag("TriggerResults","","HLT"),
    electron_id = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70Bitmap"),
    electron_muon_veto_dR = cms.double(-1),
    electron_src = cms.InputTag("slimmedElectrons"),
    gen_src = cms.InputTag("prunedGenParticles"),
    hlt_filter_ele = cms.vstring('hltDiEle33CaloIdLMWPMS2UnseededFilter'),
    l1_filter_ele = cms.string('hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter'),
    muon_cuts = cms.string('isGlobalMuon && isTrackerMuon && pt > 53 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && (globalTrack.hitPattern.numberOfValidMuonHits > 0 || tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0 ) && ( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))'),
    muon_photon_match_src = cms.InputTag("muonPhotonMatchMiniAOD"),
    muon_src = cms.InputTag("slimmedMuons"),
    muon_srcSecond = cms.InputTag("slimmedMuons"),
    muon_track_for_momentum = cms.string('TunePNew'),
    muon_track_for_momentum_CSC = cms.string('Inner'),
    prescaled_trigger_filters = cms.vstring('hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q'),
    prescaled_trigger_path_names = cms.vstring('Mu27'),
    prescales = cms.InputTag("patTrigger"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    trigger_filters = cms.vstring(
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q', 
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q', 
        'hltL3fL1sMu25f0TkFiltered100Q'
    ),
    trigger_match_max_dR = cms.double(0.2),
    trigger_path_full_names = cms.vstring(
        'HLT_Mu50_v*', 
        'HLT_OldMu100_v*', 
        'HLT_TkMu100_v*'
    ),
    trigger_path_names = cms.vstring(
        'Mu50', 
        'OldMu100', 
        'TkMu100'
    ),
    trigger_summary = cms.InputTag("slimmedPatTrigger"),
    vtx_src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.Our2017MuPrescaledCommonLeptons = cms.EDProducer("Zprime2muLeptonProducer_miniAOD",
    bits = cms.InputTag("TriggerResults","","HLT"),
    electron_id = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70Bitmap"),
    electron_muon_veto_dR = cms.double(-1),
    electron_src = cms.InputTag("slimmedElectrons"),
    gen_src = cms.InputTag("prunedGenParticles"),
    hlt_filter_ele = cms.vstring('hltDiEle33CaloIdLMWPMS2UnseededFilter'),
    l1_filter_ele = cms.string('hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter'),
    muon_cuts = cms.string('isGlobalMuon && isTrackerMuon && pt > 27 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && (globalTrack.hitPattern.numberOfValidMuonHits > 0 || tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0 ) && ( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))'),
    muon_photon_match_src = cms.InputTag("muonPhotonMatchMiniAOD"),
    muon_src = cms.InputTag("slimmedMuons"),
    muon_srcSecond = cms.InputTag("slimmedMuons"),
    muon_track_for_momentum = cms.string('TunePNew'),
    muon_track_for_momentum_CSC = cms.string('Inner'),
    prescaled_trigger_filters = cms.vstring('hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q'),
    prescaled_trigger_path_names = cms.vstring('Mu27'),
    prescales = cms.InputTag("patTrigger"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    trigger_filters = cms.vstring(
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q', 
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q', 
        'hltL3fL1sMu25f0TkFiltered100Q'
    ),
    trigger_match_max_dR = cms.double(0.2),
    trigger_path_full_names = cms.vstring(
        'HLT_Mu50_v*', 
        'HLT_OldMu100_v*', 
        'HLT_TkMu100_v*'
    ),
    trigger_path_names = cms.vstring(
        'Mu50', 
        'OldMu100', 
        'TkMu100'
    ),
    trigger_summary = cms.InputTag("slimmedPatTrigger"),
    vtx_src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.Our2017MuPrescaledCommonMuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCompositeCandidatePicker",
    back_to_back_cos_angle_min = cms.double(-0.9998),
    cut = cms.string(''),
    do_remove_overlap = cms.bool(True),
    dpt_over_pt_max = cms.double(0.3),
    max_candidates = cms.uint32(1),
    prefer_Z = cms.bool(True),
    sort_by_pt = cms.bool(True),
    src = cms.InputTag("allOur2017MuPrescaledCommonMuonsPlusMuonsMinus"),
    vertex_chi2_max = cms.double(20)
)


process.Our2017MuPrescaledLeptons = cms.EDProducer("Zprime2muLeptonProducer_miniAOD",
    bits = cms.InputTag("TriggerResults","","HLT"),
    electron_id = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70Bitmap"),
    electron_muon_veto_dR = cms.double(-1),
    electron_src = cms.InputTag("slimmedElectrons"),
    gen_src = cms.InputTag("prunedGenParticles"),
    hlt_filter_ele = cms.vstring('hltDiEle33CaloIdLMWPMS2UnseededFilter'),
    l1_filter_ele = cms.string('hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter'),
    muon_cuts = cms.string('isGlobalMuon && isTrackerMuon && pt > 27 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && (globalTrack.hitPattern.numberOfValidMuonHits > 0 || tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0 ) && ( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))'),
    muon_photon_match_src = cms.InputTag("muonPhotonMatchMiniAOD"),
    muon_src = cms.InputTag("slimmedMuons"),
    muon_srcSecond = cms.InputTag("slimmedMuons"),
    muon_track_for_momentum = cms.string('TunePNew'),
    muon_track_for_momentum_CSC = cms.string('Inner'),
    prescaled_trigger_filters = cms.vstring('hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q'),
    prescaled_trigger_path_names = cms.vstring('Mu27'),
    prescales = cms.InputTag("patTrigger"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    trigger_filters = cms.vstring(
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q', 
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q', 
        'hltL3fL1sMu25f0TkFiltered100Q'
    ),
    trigger_match_max_dR = cms.double(0.2),
    trigger_path_full_names = cms.vstring(
        'HLT_Mu50_v*', 
        'HLT_OldMu100_v*', 
        'HLT_TkMu100_v*'
    ),
    trigger_path_names = cms.vstring(
        'Mu50', 
        'OldMu100', 
        'TkMu100'
    ),
    trigger_summary = cms.InputTag("slimmedPatTrigger"),
    vtx_src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.Our2017MuPrescaledMuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCompositeCandidatePicker",
    back_to_back_cos_angle_min = cms.double(-0.9998),
    cut = cms.string(''),
    do_remove_overlap = cms.bool(True),
    dpt_over_pt_max = cms.double(0.3),
    max_candidates = cms.uint32(1),
    prefer_Z = cms.bool(True),
    sort_by_pt = cms.bool(True),
    src = cms.InputTag("allOur2017MuPrescaledMuonsPlusMuonsMinus"),
    vertex_chi2_max = cms.double(20)
)


process.Our2017MuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCompositeCandidatePicker",
    back_to_back_cos_angle_min = cms.double(-0.9998),
    cut = cms.string(''),
    do_remove_overlap = cms.bool(True),
    dpt_over_pt_max = cms.double(0.3),
    max_candidates = cms.uint32(1),
    prefer_Z = cms.bool(True),
    sort_by_pt = cms.bool(True),
    src = cms.InputTag("allOur2017MuonsPlusMuonsMinus"),
    vertex_chi2_max = cms.double(20)
)


process.Our2018Leptons = cms.EDProducer("Zprime2muLeptonProducer_miniAOD",
    bits = cms.InputTag("TriggerResults","","HLT"),
    electron_id = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70Bitmap"),
    electron_muon_veto_dR = cms.double(-1),
    electron_src = cms.InputTag("slimmedElectrons"),
    gen_src = cms.InputTag("prunedGenParticles"),
    hlt_filter_ele = cms.vstring('hltDiEle33CaloIdLMWPMS2UnseededFilter'),
    l1_filter_ele = cms.string('hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter'),
    muon_cuts = cms.string('isGlobalMuon && isTrackerMuon && pt > 53 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    muon_photon_match_src = cms.InputTag("muonPhotonMatchMiniAOD"),
    muon_src = cms.InputTag("slimmedMuons"),
    muon_srcSecond = cms.InputTag("slimmedMuons"),
    muon_track_for_momentum = cms.string('TunePNew'),
    muon_track_for_momentum_CSC = cms.string('Inner'),
    prescaled_trigger_filters = cms.vstring('hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q'),
    prescaled_trigger_path_names = cms.vstring('Mu27'),
    prescales = cms.InputTag("patTrigger"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    trigger_filters = cms.vstring(
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q', 
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q', 
        'hltL3fL1sMu25f0TkFiltered100Q'
    ),
    trigger_match_max_dR = cms.double(0.2),
    trigger_path_full_names = cms.vstring(
        'HLT_Mu50_v*', 
        'HLT_OldMu100_v*', 
        'HLT_TkMu100_v*'
    ),
    trigger_path_names = cms.vstring(
        'Mu50', 
        'OldMu100', 
        'TkMu100'
    ),
    trigger_summary = cms.InputTag("slimmedPatTrigger"),
    vtx_src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.Our2018MuPrescaledCommonLeptons = cms.EDProducer("Zprime2muLeptonProducer_miniAOD",
    bits = cms.InputTag("TriggerResults","","HLT"),
    electron_id = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70Bitmap"),
    electron_muon_veto_dR = cms.double(-1),
    electron_src = cms.InputTag("slimmedElectrons"),
    gen_src = cms.InputTag("prunedGenParticles"),
    hlt_filter_ele = cms.vstring('hltDiEle33CaloIdLMWPMS2UnseededFilter'),
    l1_filter_ele = cms.string('hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter'),
    muon_cuts = cms.string('isGlobalMuon && isTrackerMuon && pt > 27 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    muon_photon_match_src = cms.InputTag("muonPhotonMatchMiniAOD"),
    muon_src = cms.InputTag("slimmedMuons"),
    muon_srcSecond = cms.InputTag("slimmedMuons"),
    muon_track_for_momentum = cms.string('TunePNew'),
    muon_track_for_momentum_CSC = cms.string('Inner'),
    prescaled_trigger_filters = cms.vstring('hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q'),
    prescaled_trigger_path_names = cms.vstring('Mu27'),
    prescales = cms.InputTag("patTrigger"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    trigger_filters = cms.vstring(
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q', 
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q', 
        'hltL3fL1sMu25f0TkFiltered100Q'
    ),
    trigger_match_max_dR = cms.double(0.2),
    trigger_path_full_names = cms.vstring(
        'HLT_Mu50_v*', 
        'HLT_OldMu100_v*', 
        'HLT_TkMu100_v*'
    ),
    trigger_path_names = cms.vstring(
        'Mu50', 
        'OldMu100', 
        'TkMu100'
    ),
    trigger_summary = cms.InputTag("slimmedPatTrigger"),
    vtx_src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.Our2018MuPrescaledCommonMuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCompositeCandidatePicker",
    back_to_back_cos_angle_min = cms.double(-0.9998),
    cut = cms.string(''),
    do_remove_overlap = cms.bool(True),
    dpt_over_pt_max = cms.double(0.3),
    max_candidates = cms.uint32(1),
    prefer_Z = cms.bool(True),
    sort_by_pt = cms.bool(True),
    src = cms.InputTag("allOur2018MuPrescaledCommonMuonsPlusMuonsMinus"),
    vertex_chi2_max = cms.double(20)
)


process.Our2018MuPrescaledLeptons = cms.EDProducer("Zprime2muLeptonProducer_miniAOD",
    bits = cms.InputTag("TriggerResults","","HLT"),
    electron_id = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70Bitmap"),
    electron_muon_veto_dR = cms.double(-1),
    electron_src = cms.InputTag("slimmedElectrons"),
    gen_src = cms.InputTag("prunedGenParticles"),
    hlt_filter_ele = cms.vstring('hltDiEle33CaloIdLMWPMS2UnseededFilter'),
    l1_filter_ele = cms.string('hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter'),
    muon_cuts = cms.string('isGlobalMuon && isTrackerMuon && pt > 27 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    muon_photon_match_src = cms.InputTag("muonPhotonMatchMiniAOD"),
    muon_src = cms.InputTag("slimmedMuons"),
    muon_srcSecond = cms.InputTag("slimmedMuons"),
    muon_track_for_momentum = cms.string('TunePNew'),
    muon_track_for_momentum_CSC = cms.string('Inner'),
    prescaled_trigger_filters = cms.vstring('hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q'),
    prescaled_trigger_path_names = cms.vstring('Mu27'),
    prescales = cms.InputTag("patTrigger"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    trigger_filters = cms.vstring(
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q', 
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q', 
        'hltL3fL1sMu25f0TkFiltered100Q'
    ),
    trigger_match_max_dR = cms.double(0.2),
    trigger_path_full_names = cms.vstring(
        'HLT_Mu50_v*', 
        'HLT_OldMu100_v*', 
        'HLT_TkMu100_v*'
    ),
    trigger_path_names = cms.vstring(
        'Mu50', 
        'OldMu100', 
        'TkMu100'
    ),
    trigger_summary = cms.InputTag("slimmedPatTrigger"),
    vtx_src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.Our2018MuPrescaledMuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCompositeCandidatePicker",
    back_to_back_cos_angle_min = cms.double(-0.9998),
    cut = cms.string(''),
    do_remove_overlap = cms.bool(True),
    dpt_over_pt_max = cms.double(0.3),
    max_candidates = cms.uint32(1),
    prefer_Z = cms.bool(True),
    sort_by_pt = cms.bool(True),
    src = cms.InputTag("allOur2018MuPrescaledMuonsPlusMuonsMinus"),
    vertex_chi2_max = cms.double(20)
)


process.Our2018MuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCompositeCandidatePicker",
    back_to_back_cos_angle_min = cms.double(-0.9998),
    cut = cms.string(''),
    do_remove_overlap = cms.bool(True),
    dpt_over_pt_max = cms.double(0.3),
    max_candidates = cms.uint32(1),
    prefer_Z = cms.bool(True),
    sort_by_pt = cms.bool(True),
    src = cms.InputTag("allOur2018MuonsPlusMuonsMinus"),
    vertex_chi2_max = cms.double(20)
)


process.allDimuons = cms.EDProducer("Zprime2muCombiner",
    cut = cms.string(''),
    decay = cms.string('leptons:muons@+ leptons:muons@-'),
    loose_cut = cms.string('isGlobalMuon && isTrackerMuon && pt > 53 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    tight_cut = cms.string('userFloat("Mu50_TriggerMatchPt")>=50 || userFloat("OldMu100_TriggerMatchPt")>=100 || userFloat("TkMu100_TriggerMatchPt")>=100 ')
)


process.allOur2017MuPrescaledCommonMuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCombiner",
    cut = cms.string('daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    decay = cms.string('Our2017MuPrescaledCommonLeptons:muons@+ Our2017MuPrescaledCommonLeptons:muons@-'),
    loose_cut = cms.string('isGlobalMuon && isTrackerMuon && pt > 27 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && (globalTrack.hitPattern.numberOfValidMuonHits > 0 || tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0 ) && ( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))'),
    tight_cut = cms.string('userFloat("prescaledMu27_TriggerMatchPt")>=27 ')
)


process.allOur2017MuPrescaledMuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCombiner",
    cut = cms.string('daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    decay = cms.string('Our2017MuPrescaledLeptons:muons@+ Our2017MuPrescaledLeptons:muons@-'),
    loose_cut = cms.string('isGlobalMuon && isTrackerMuon && pt > 27 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && (globalTrack.hitPattern.numberOfValidMuonHits > 0 || tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0 ) && ( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))'),
    tight_cut = cms.string('userFloat("prescaledMu27_TriggerMatchPt")>=27 ')
)


process.allOur2017MuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCombiner",
    cut = cms.string('daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    decay = cms.string('Our2017Leptons:muons@+ Our2017Leptons:muons@-'),
    loose_cut = cms.string('isGlobalMuon && isTrackerMuon && pt > 53 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && (globalTrack.hitPattern.numberOfValidMuonHits > 0 || tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0 ) && ( numberOfMatchedStations > 1 || (numberOfMatchedStations == 1 && !(stationMask == 1 || stationMask == 16)) || (numberOfMatchedStations == 1 && (stationMask == 1 || stationMask == 16) && numberOfMatchedRPCLayers > 2))'),
    tight_cut = cms.string('userFloat("Mu50_TriggerMatchPt")>=50 || userFloat("OldMu100_TriggerMatchPt")>=100 || userFloat("TkMu100_TriggerMatchPt")>=100 ')
)


process.allOur2018MuPrescaledCommonMuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCombiner",
    cut = cms.string('daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    decay = cms.string('Our2018MuPrescaledCommonLeptons:muons@+ Our2018MuPrescaledCommonLeptons:muons@-'),
    loose_cut = cms.string('isGlobalMuon && isTrackerMuon && pt > 27 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    tight_cut = cms.string('userFloat("prescaledMu27_TriggerMatchPt")>=27 ')
)


process.allOur2018MuPrescaledMuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCombiner",
    cut = cms.string('daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    decay = cms.string('Our2018MuPrescaledLeptons:muons@+ Our2018MuPrescaledLeptons:muons@-'),
    loose_cut = cms.string('isGlobalMuon && isTrackerMuon && pt > 27 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    tight_cut = cms.string('userFloat("prescaledMu27_TriggerMatchPt")>=27 ')
)


process.allOur2018MuonsPlusMuonsMinus = cms.EDProducer("Zprime2muCombiner",
    cut = cms.string('daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    decay = cms.string('Our2018Leptons:muons@+ Our2018Leptons:muons@-'),
    loose_cut = cms.string('isGlobalMuon && isTrackerMuon && pt > 53 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    tight_cut = cms.string('userFloat("Mu50_TriggerMatchPt")>=50 || userFloat("OldMu100_TriggerMatchPt")>=100 || userFloat("TkMu100_TriggerMatchPt")>=100 ')
)


process.dimuons = cms.EDProducer("Zprime2muCompositeCandidatePicker",
    back_to_back_cos_angle_min = cms.double(-0.9998),
    cut = cms.string(''),
    do_remove_overlap = cms.bool(True),
    dpt_over_pt_max = cms.double(0.3),
    max_candidates = cms.uint32(1),
    prefer_Z = cms.bool(True),
    sort_by_pt = cms.bool(True),
    src = cms.InputTag("allDimuons"),
    vertex_chi2_max = cms.double(20)
)


process.egmGsfElectronIDs = cms.EDProducer("VersionedGsfElectronIdProducer",
    physicsObjectIDs = cms.VPSet(cms.PSet(
        idDefinition = cms.PSet(
            cutFlow = cms.VPSet(
                cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(35.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    allowedEtaRanges = cms.VPSet(
                        cms.PSet(
                            maxEta = cms.double(1.4442),
                            minEta = cms.double(0.0)
                        ), 
                        cms.PSet(
                            maxEta = cms.double(2.5),
                            minEta = cms.double(1.566)
                        )
                    ),
                    cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False),
                    useAbsEta = cms.bool(True)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleDEtaInSeedCut'),
                    dEtaInSeedCutValueEB = cms.double(0.004),
                    dEtaInSeedCutValueEE = cms.double(0.006),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleDPhiInCut'),
                    dPhiInCutValueEB = cms.double(0.06),
                    dPhiInCutValueEE = cms.double(0.06),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaWithSatCut'),
                    isIgnored = cms.bool(False),
                    maxNrSatCrysIn5x5EB = cms.int32(0),
                    maxNrSatCrysIn5x5EE = cms.int32(0),
                    maxSigmaIEtaIEtaEB = cms.double(9999),
                    maxSigmaIEtaIEtaEE = cms.double(0.03),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    cutName = cms.string('GsfEleFull5x5E2x5OverE5x5WithSatCut'),
                    isIgnored = cms.bool(False),
                    maxNrSatCrysIn5x5EB = cms.int32(0),
                    maxNrSatCrysIn5x5EE = cms.int32(0),
                    minE1x5OverE5x5EB = cms.double(0.83),
                    minE1x5OverE5x5EE = cms.double(-1.0),
                    minE2x5OverE5x5EB = cms.double(0.94),
                    minE2x5OverE5x5EE = cms.double(-1.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    constTermEB = cms.double(1.0),
                    constTermEE = cms.double(5),
                    cutName = cms.string('GsfEleHadronicOverEMLinearCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False),
                    slopeStartEB = cms.double(0.0),
                    slopeStartEE = cms.double(0.0),
                    slopeTermEB = cms.double(0.05),
                    slopeTermEE = cms.double(0.05)
                ), 
                cms.PSet(
                    constTermEB = cms.double(5.0),
                    constTermEE = cms.double(5.0),
                    cutName = cms.string('GsfEleValueMapIsoRhoCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag(""),
                    rhoEAEB = cms.double(0.0),
                    rhoEAEE = cms.double(0.0),
                    rhoEtStartEB = cms.double(999999.0),
                    rhoEtStartEE = cms.double(999999.0),
                    slopeStartEB = cms.double(0.0),
                    slopeStartEE = cms.double(0.0),
                    slopeTermEB = cms.double(0.0),
                    slopeTermEE = cms.double(0.0),
                    value = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso")
                ), 
                cms.PSet(
                    constTermEB = cms.double(2.0),
                    constTermEE = cms.double(2.5),
                    cutName = cms.string('GsfEleEmHadD1IsoRhoCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                    rhoConstant = cms.double(0.28),
                    slopeStartEB = cms.double(0.0),
                    slopeStartEE = cms.double(50.0),
                    slopeTermEB = cms.double(0.03),
                    slopeTermEE = cms.double(0.03)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleDxyCut'),
                    dxyCutValueEB = cms.double(0.02),
                    dxyCutValueEE = cms.double(0.05),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    vertexSrc = cms.InputTag("offlinePrimaryVertices"),
                    vertexSrcMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleMissingHitsCut'),
                    isIgnored = cms.bool(False),
                    maxMissingHitsEB = cms.uint32(1),
                    maxMissingHitsEE = cms.uint32(1),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleEcalDrivenCut'),
                    ecalDrivenEB = cms.int32(1),
                    ecalDrivenEE = cms.int32(1),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                )
            ),
            idName = cms.string('heepElectronID-HEEPV70'),
            isPOGApproved = cms.untracked.bool(True)
        ),
        idMD5 = cms.string('49b6b60e9f16727f241eb34b9d345a8f'),
        isPOGApproved = cms.untracked.bool(True)
    )),
    physicsObjectSrc = cms.InputTag("slimmedElectrons")
)


process.electronMVAValueMapProducer = cms.EDProducer("ElectronMVAValueMapProducer",
    mvaConfigurations = cms.VPSet(
        cms.PSet(
            categoryCuts = cms.vstring(
                'pt < 10. && abs(superCluster.eta) < 0.800', 
                'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt < 10. && abs(superCluster.eta) >= 1.479', 
                'pt >= 10. && abs(superCluster.eta) < 0.800', 
                'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt >= 10. && abs(superCluster.eta) >= 1.479'
            ),
            mvaName = cms.string('ElectronMVAEstimatorRun2'),
            mvaTag = cms.string('Spring16HZZV1'),
            nCategories = cms.int32(6),
            variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml.gz'
            )
        ), 
        cms.PSet(
            categoryCuts = cms.vstring(
                'abs(superCluster.eta) < 0.800', 
                'abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'abs(superCluster.eta) >= 1.479'
            ),
            mvaName = cms.string('ElectronMVAEstimatorRun2'),
            mvaTag = cms.string('Spring16GeneralPurposeV1'),
            nCategories = cms.int32(3),
            variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml.gz'
            )
        ), 
        cms.PSet(
            categoryCuts = cms.vstring(
                'pt < 10. && abs(superCluster.eta) < 0.800', 
                'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt < 10. && abs(superCluster.eta) >= 1.479', 
                'pt >= 10. && abs(superCluster.eta) < 0.800', 
                'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt >= 10. && abs(superCluster.eta) >= 1.479'
            ),
            mvaName = cms.string('ElectronMVAEstimatorRun2'),
            mvaTag = cms.string('Fall17NoIsoV1'),
            nCategories = cms.int32(6),
            variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Fall17V1Variables.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml.gz'
            )
        ), 
        cms.PSet(
            categoryCuts = cms.vstring(
                'pt < 10. && abs(superCluster.eta) < 0.800', 
                'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt < 10. && abs(superCluster.eta) >= 1.479', 
                'pt >= 10. && abs(superCluster.eta) < 0.800', 
                'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt >= 10. && abs(superCluster.eta) >= 1.479'
            ),
            mvaName = cms.string('ElectronMVAEstimatorRun2'),
            mvaTag = cms.string('Fall17IsoV1'),
            nCategories = cms.int32(6),
            variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Fall17V1Variables.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml.gz'
            )
        ), 
        cms.PSet(
            categoryCuts = cms.vstring(
                'pt < 10. && abs(superCluster.eta) < 0.800', 
                'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt < 10. && abs(superCluster.eta) >= 1.479', 
                'pt >= 10. && abs(superCluster.eta) < 0.800', 
                'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt >= 10. && abs(superCluster.eta) >= 1.479'
            ),
            mvaName = cms.string('ElectronMVAEstimatorRun2'),
            mvaTag = cms.string('Fall17NoIsoV2'),
            nCategories = cms.int32(6),
            variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB1_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB2_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EE_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB1_10.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EB2_10.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17NoIsoV2/EE_10.weights.xml.gz'
            )
        ), 
        cms.PSet(
            categoryCuts = cms.vstring(
                'pt < 10. && abs(superCluster.eta) < 0.800', 
                'pt < 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt < 10. && abs(superCluster.eta) >= 1.479', 
                'pt >= 10. && abs(superCluster.eta) < 0.800', 
                'pt >= 10. && abs(superCluster.eta) >= 0.800 && abs(superCluster.eta) < 1.479', 
                'pt >= 10. && abs(superCluster.eta) >= 1.479'
            ),
            mvaName = cms.string('ElectronMVAEstimatorRun2'),
            mvaTag = cms.string('Fall17IsoV2'),
            nCategories = cms.int32(6),
            variableDefinition = cms.string('RecoEgamma/ElectronIdentification/data/ElectronMVAEstimatorRun2Variables.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB1_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB2_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EE_5.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB1_10.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EB2_10.weights.xml.gz', 
                'RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Fall17IsoV2/EE_10.weights.xml.gz'
            )
        )
    ),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("slimmedElectrons","","@skipCurrentProcess")
)


process.electronMVAVariableHelper = cms.EDProducer("GsfElectronMVAVariableHelper",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    beamSpotMiniAOD = cms.InputTag("offlineBeamSpot"),
    conversions = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    vertexCollectionMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.electronRegressionValueMapProducer = cms.EDProducer("ElectronRegressionValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"),
    useFull5x5 = cms.bool(False)
)


process.heepIDVarValueMaps = cms.EDProducer("ElectronHEEPIDValueMapProducer",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    candVetosAOD = cms.vstring(
        'ELES', 
        'NONE', 
        'NONELES'
    ),
    candVetosMiniAOD = cms.vstring(
        'ELES', 
        'NONE', 
        'NONELES'
    ),
    candsAOD = cms.VInputTag("packedCandsForTkIso", "lostTracksForTkIso", "lostTracksForTkIso:eleTracks"),
    candsMiniAOD = cms.VInputTag("packedPFCandidates", "lostTracks", "lostTracks:eleTracks"),
    dataFormat = cms.int32(0),
    ebRecHitsAOD = cms.InputTag("reducedEcalRecHitsEB"),
    ebRecHitsMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeRecHitsAOD = cms.InputTag("reducedEcalRecHitsEB"),
    eeRecHitsMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    elesAOD = cms.InputTag("gedGsfElectrons"),
    elesMiniAOD = cms.InputTag("slimmedElectrons"),
    trkIsoConfig = cms.PSet(
        barrelCuts = cms.PSet(
            algosToReject = cms.vstring(),
            allowedQualities = cms.vstring(),
            maxDPtPt = cms.double(0.1),
            maxDR = cms.double(0.3),
            maxDZ = cms.double(0.1),
            minDEta = cms.double(0.005),
            minDR = cms.double(0.0),
            minHits = cms.int32(8),
            minPixelHits = cms.int32(1),
            minPt = cms.double(1.0)
        ),
        endcapCuts = cms.PSet(
            algosToReject = cms.vstring(),
            allowedQualities = cms.vstring(),
            maxDPtPt = cms.double(0.1),
            maxDR = cms.double(0.3),
            maxDZ = cms.double(0.5),
            minDEta = cms.double(0.005),
            minDR = cms.double(0.0),
            minHits = cms.int32(8),
            minPixelHits = cms.int32(1),
            minPt = cms.double(1.0)
        )
    )
)


process.leptons = cms.EDProducer("Zprime2muLeptonProducer",
    electron_cuts = cms.string(''),
    electron_muon_veto_dR = cms.double(-1),
    electron_src = cms.InputTag("cleanPatElectrons"),
    electron_srcSecond = cms.InputTag("cleanPatElectrons"),
    muon_cuts = cms.string('isGlobalMuon && isTrackerMuon && pt > 53 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    muon_photon_match_src = cms.InputTag("muonPhotonMatch"),
    muon_src = cms.InputTag("cleanPatMuonsTriggerMatch"),
    muon_srcSecond = cms.InputTag("cleanPatMuonsTriggerMatch"),
    muon_track_for_momentum = cms.string('TunePNew'),
    muon_track_for_momentum_CSC = cms.string('Inner'),
    trigger_match_max_dR = cms.double(0.2),
    trigger_summary_src = cms.InputTag("hltTriggerSummaryAOD","","HLT")
)


process.leptonsMini = cms.EDProducer("Zprime2muLeptonProducer_miniAOD",
    bits = cms.InputTag("TriggerResults","","HLT"),
    electron_id = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70Bitmap"),
    electron_muon_veto_dR = cms.double(-1),
    electron_src = cms.InputTag("slimmedElectrons"),
    gen_src = cms.InputTag("prunedGenParticles"),
    hlt_filter_ele = cms.vstring('hltDiEle33CaloIdLMWPMS2UnseededFilter'),
    l1_filter_ele = cms.string('hltEGL1SingleAndDoubleEGNonIsoOrWithEG26WithJetAndTauFilter'),
    muon_cuts = cms.string('isGlobalMuon && isTrackerMuon && pt > 53 && abs(eta) < 2.4 && abs(dB) < 0.2 && isolationR03.sumPt / innerTrack.pt < 0.10 && globalTrack.hitPattern.trackerLayersWithMeasurement > 5 && globalTrack.hitPattern.numberOfValidPixelHits >= 1 && ( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) ) && (( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))'),
    muon_photon_match_src = cms.InputTag("muonPhotonMatchMiniAOD"),
    muon_src = cms.InputTag("slimmedMuons"),
    muon_srcSecond = cms.InputTag("slimmedMuons"),
    muon_track_for_momentum = cms.string('TunePNew'),
    muon_track_for_momentum_CSC = cms.string('Inner'),
    prescaled_trigger_filters = cms.vstring(),
    prescaled_trigger_path_names = cms.vstring(),
    prescales = cms.InputTag("patTrigger"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    trigger_filters = cms.vstring(),
    trigger_match_max_dR = cms.double(0.2),
    trigger_path_full_names = cms.vstring(),
    trigger_path_names = cms.vstring(),
    trigger_summary = cms.InputTag("slimmedPatTrigger"),
    vtx_src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.muonPhotonMatch = cms.EDProducer("TrivialDeltaRViewMatcher",
    distMin = cms.double(0.1),
    matched = cms.InputTag("cleanPatPhotons"),
    src = cms.InputTag("cleanPatMuonsTriggerMatch")
)


process.muonPhotonMatchMiniAOD = cms.EDProducer("TrivialDeltaRViewMatcher",
    distMin = cms.double(0.1),
    matched = cms.InputTag("slimmedPhotons"),
    src = cms.InputTag("slimmedMuons")
)


process.prunedMCLeptons = cms.EDProducer("GenParticlePruner",
    select = cms.vstring(
        'drop *', 
        'keep status == 3', 
        '++keep abs(pdgId) == 13 && (status == 1 || status == 8)', 
        '++keep abs(pdgId) == 11 && status == 1', 
        'keep++ pdgId == 23 || pdgId == 32 || pdgId == 39 || pdgId == 5000039'
    ),
    src = cms.InputTag("prunedGenParticles")
)


process.PrescaleToCommon = cms.EDFilter("PrescaleToCommon_miniAOD",
    L1Prescale_max_src = cms.InputTag("patTrigger","l1max","PAT"),
    L1Prescale_min_src = cms.InputTag("patTrigger","l1min","PAT"),
    Prescale_src = cms.InputTag("patTrigger","","PAT"),
    TriggerResults_src = cms.InputTag("TriggerResults","","HLT"),
    assume_simulation_has_prescale_1 = cms.bool(True),
    debugInfo = cms.bool(False),
    hlt_src = cms.InputTag("TriggerResults","","HLT"),
    no_common = cms.bool(False),
    overall_prescale = cms.int32(500),
    trigger_paths = cms.vstring(
        'HLT_Mu27_v12', 
        'HLT_Mu27_v13'
    )
)


process.PrescaleToCommonMiniAOD = cms.EDFilter("PrescaleToCommon_miniAOD",
    L1Prescale_max_src = cms.InputTag("patTrigger","l1max","PAT"),
    L1Prescale_min_src = cms.InputTag("patTrigger","l1min","PAT"),
    Prescale_src = cms.InputTag("patTrigger","","PAT"),
    TriggerResults_src = cms.InputTag("TriggerResults","","HLT"),
    assume_simulation_has_prescale_1 = cms.bool(True),
    debugInfo = cms.bool(False),
    hlt_src = cms.InputTag("TriggerResults","","HLT"),
    no_common = cms.bool(False),
    overall_prescale = cms.int32(500),
    trigger_paths = cms.vstring(
        'HLT_Mu27_v12', 
        'HLT_Mu27_v13'
    )
)


process.beamHaloMiniAOD = cms.EDFilter("METFilterMiniAOD",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    flag = cms.string('Flag_globalTightHalo2016Filter'),
    src = cms.InputTag("TriggerResults","","RECO")
)


process.defaultSelector = cms.EDFilter("METFilterMiniAOD",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    flag = cms.string('Flag_globalTightHalo2016Filter'),
    src = cms.InputTag("TriggerResults","","RECO")
)


process.dileptonPreseletor = cms.EDFilter("DileptonPreselector",
    muons = cms.InputTag("slimmedMuons"),
    nMuons = cms.double(2),
    ptCut = cms.double(40)
)


process.ecalDeadCellPrimitiveMiniAOD = cms.EDFilter("METFilterMiniAOD",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    flag = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
    src = cms.InputTag("TriggerResults","","RECO")
)


process.eeBadScMiniAOD = cms.EDFilter("METFilterMiniAOD",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    flag = cms.string('Flag_eeBadScFilter'),
    src = cms.InputTag("TriggerResults","","RECO")
)


process.goodDataFilter = cms.EDFilter("HLTHighLevel",
    HLTPaths = cms.vstring('goodDataPrimaryVertexFilter'),
    TriggerResultsTag = cms.InputTag("TriggerResults","","PAT"),
    andOr = cms.bool(False),
    eventSetupPathsKey = cms.string(''),
    throw = cms.bool(True)
)


process.hbheIsoNoiseMiniAOD = cms.EDFilter("METFilterMiniAOD",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    flag = cms.string('Flag_HBHENoiseIsoFilter'),
    src = cms.InputTag("TriggerResults","","RECO")
)


process.hbheNoiseMiniAOD = cms.EDFilter("METFilterMiniAOD",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    flag = cms.string('Flag_HBHENoiseFilter'),
    src = cms.InputTag("TriggerResults","","RECO")
)


process.hltPhysicsDeclared = cms.EDFilter("HLTPhysicsDeclared",
    L1GtReadoutRecordTag = cms.InputTag("gtDigis"),
    invert = cms.bool(False)
)


process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)


process.primaryVertex = cms.EDFilter("GoodVertexFilter",
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2),
    minimumNDOF = cms.uint32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.primaryVertexMiniAOD = cms.EDFilter("METFilterMiniAOD",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    flag = cms.string('Flag_goodVertices'),
    src = cms.InputTag("TriggerResults","","PAT")
)


process.EventCounter = cms.EDAnalyzer("EventCounter",
    genInfoTag = cms.InputTag("generator")
)


process.Our2017MuPrescaledCommonMuonsPlusMuonsMinusHistos = cms.EDAnalyzer("Zprime2muHistosFromPAT",
    LHEInfo = cms.InputTag("source"),
    beamspot_src = cms.InputTag("offlineBeamSpot"),
    dilepton_src = cms.InputTag("Our2017MuPrescaledCommonMuonsPlusMuonsMinus"),
    doElectrons = cms.bool(False),
    hardInteraction = cms.PSet(
        allowFakeResonance = cms.bool(True),
        doingElectrons = cms.bool(False),
        matchTaus = cms.bool(True),
        resonanceIds = cms.vint32(32, 23, 39, 5000039),
        shutUp = cms.bool(True),
        src = cms.InputTag("prunedGenParticles")
    ),
    lepton_src = cms.InputTag("Our2017MuPrescaledCommonLeptons","muons"),
    leptonsFromDileptons = cms.bool(True),
    lrWeightProducer = cms.PSet(
        Lambda = cms.int32(16000),
        calculate = cms.bool(False),
        doingElectrons = cms.bool(False),
        doingLR = cms.bool(True),
        interference = cms.int32(1),
        src = cms.InputTag("prunedMCLeptons")
    ),
    prefireWeights = cms.InputTag("prefiringweight","nonPrefiringProb"),
    pu_src = cms.InputTag("slimmedAddPileupInfo"),
    pu_weights = cms.vstring(),
    useMadgraphWeight = cms.bool(True),
    useTTBarWeight = cms.bool(False),
    use_bs_and_pv = cms.bool(True),
    usekFactor = cms.bool(True),
    vertex_src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    year = cms.int32(2018)
)


process.Our2017MuPrescaledMuonsPlusMuonsMinusHistos = cms.EDAnalyzer("Zprime2muHistosFromPAT",
    LHEInfo = cms.InputTag("source"),
    beamspot_src = cms.InputTag("offlineBeamSpot"),
    dilepton_src = cms.InputTag("Our2017MuPrescaledMuonsPlusMuonsMinus"),
    doElectrons = cms.bool(False),
    hardInteraction = cms.PSet(
        allowFakeResonance = cms.bool(True),
        doingElectrons = cms.bool(False),
        matchTaus = cms.bool(True),
        resonanceIds = cms.vint32(32, 23, 39, 5000039),
        shutUp = cms.bool(True),
        src = cms.InputTag("prunedGenParticles")
    ),
    lepton_src = cms.InputTag("Our2017MuPrescaledLeptons","muons"),
    leptonsFromDileptons = cms.bool(True),
    lrWeightProducer = cms.PSet(
        Lambda = cms.int32(16000),
        calculate = cms.bool(False),
        doingElectrons = cms.bool(False),
        doingLR = cms.bool(True),
        interference = cms.int32(1),
        src = cms.InputTag("prunedMCLeptons")
    ),
    prefireWeights = cms.InputTag("prefiringweight","nonPrefiringProb"),
    pu_src = cms.InputTag("slimmedAddPileupInfo"),
    pu_weights = cms.vstring(),
    useMadgraphWeight = cms.bool(True),
    useTTBarWeight = cms.bool(False),
    use_bs_and_pv = cms.bool(True),
    usekFactor = cms.bool(True),
    vertex_src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    year = cms.int32(2018)
)


process.Our2017MuonsPlusMuonsMinusHistos = cms.EDAnalyzer("Zprime2muHistosFromPAT",
    LHEInfo = cms.InputTag("source"),
    beamspot_src = cms.InputTag("offlineBeamSpot"),
    dilepton_src = cms.InputTag("Our2017MuonsPlusMuonsMinus"),
    doElectrons = cms.bool(False),
    hardInteraction = cms.PSet(
        allowFakeResonance = cms.bool(True),
        doingElectrons = cms.bool(False),
        matchTaus = cms.bool(True),
        resonanceIds = cms.vint32(32, 23, 39, 5000039),
        shutUp = cms.bool(True),
        src = cms.InputTag("prunedGenParticles")
    ),
    lepton_src = cms.InputTag("Our2017Leptons","muons"),
    leptonsFromDileptons = cms.bool(True),
    lrWeightProducer = cms.PSet(
        Lambda = cms.int32(16000),
        calculate = cms.bool(False),
        doingElectrons = cms.bool(False),
        doingLR = cms.bool(True),
        interference = cms.int32(1),
        src = cms.InputTag("prunedMCLeptons")
    ),
    prefireWeights = cms.InputTag("prefiringweight","nonPrefiringProb"),
    pu_src = cms.InputTag("slimmedAddPileupInfo"),
    pu_weights = cms.vstring(),
    useMadgraphWeight = cms.bool(True),
    useTTBarWeight = cms.bool(False),
    use_bs_and_pv = cms.bool(True),
    usekFactor = cms.bool(True),
    vertex_src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    year = cms.int32(2018)
)


process.Our2018MuPrescaledCommonMuonsPlusMuonsMinusHistos = cms.EDAnalyzer("Zprime2muHistosFromPAT",
    LHEInfo = cms.InputTag("source"),
    beamspot_src = cms.InputTag("offlineBeamSpot"),
    dilepton_src = cms.InputTag("Our2018MuPrescaledCommonMuonsPlusMuonsMinus"),
    doElectrons = cms.bool(False),
    hardInteraction = cms.PSet(
        allowFakeResonance = cms.bool(True),
        doingElectrons = cms.bool(False),
        matchTaus = cms.bool(True),
        resonanceIds = cms.vint32(32, 23, 39, 5000039),
        shutUp = cms.bool(True),
        src = cms.InputTag("prunedGenParticles")
    ),
    lepton_src = cms.InputTag("Our2018MuPrescaledCommonLeptons","muons"),
    leptonsFromDileptons = cms.bool(True),
    lrWeightProducer = cms.PSet(
        Lambda = cms.int32(16000),
        calculate = cms.bool(False),
        doingElectrons = cms.bool(False),
        doingLR = cms.bool(True),
        interference = cms.int32(1),
        src = cms.InputTag("prunedMCLeptons")
    ),
    prefireWeights = cms.InputTag("prefiringweight","nonPrefiringProb"),
    pu_src = cms.InputTag("slimmedAddPileupInfo"),
    pu_weights = cms.vstring(),
    useMadgraphWeight = cms.bool(True),
    useTTBarWeight = cms.bool(False),
    use_bs_and_pv = cms.bool(True),
    usekFactor = cms.bool(True),
    vertex_src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    year = cms.int32(2018)
)


process.Our2018MuPrescaledMuonsPlusMuonsMinusHistos = cms.EDAnalyzer("Zprime2muHistosFromPAT",
    LHEInfo = cms.InputTag("source"),
    beamspot_src = cms.InputTag("offlineBeamSpot"),
    dilepton_src = cms.InputTag("Our2018MuPrescaledMuonsPlusMuonsMinus"),
    doElectrons = cms.bool(False),
    hardInteraction = cms.PSet(
        allowFakeResonance = cms.bool(True),
        doingElectrons = cms.bool(False),
        matchTaus = cms.bool(True),
        resonanceIds = cms.vint32(32, 23, 39, 5000039),
        shutUp = cms.bool(True),
        src = cms.InputTag("prunedGenParticles")
    ),
    lepton_src = cms.InputTag("Our2018MuPrescaledLeptons","muons"),
    leptonsFromDileptons = cms.bool(True),
    lrWeightProducer = cms.PSet(
        Lambda = cms.int32(16000),
        calculate = cms.bool(False),
        doingElectrons = cms.bool(False),
        doingLR = cms.bool(True),
        interference = cms.int32(1),
        src = cms.InputTag("prunedMCLeptons")
    ),
    prefireWeights = cms.InputTag("prefiringweight","nonPrefiringProb"),
    pu_src = cms.InputTag("slimmedAddPileupInfo"),
    pu_weights = cms.vstring(),
    useMadgraphWeight = cms.bool(True),
    useTTBarWeight = cms.bool(False),
    use_bs_and_pv = cms.bool(True),
    usekFactor = cms.bool(True),
    vertex_src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    year = cms.int32(2018)
)


process.Our2018MuonsPlusMuonsMinusHistos = cms.EDAnalyzer("Zprime2muHistosFromPAT",
    LHEInfo = cms.InputTag("source"),
    beamspot_src = cms.InputTag("offlineBeamSpot"),
    dilepton_src = cms.InputTag("Our2018MuonsPlusMuonsMinus"),
    doElectrons = cms.bool(False),
    hardInteraction = cms.PSet(
        allowFakeResonance = cms.bool(True),
        doingElectrons = cms.bool(False),
        matchTaus = cms.bool(True),
        resonanceIds = cms.vint32(32, 23, 39, 5000039),
        shutUp = cms.bool(True),
        src = cms.InputTag("prunedGenParticles")
    ),
    lepton_src = cms.InputTag("Our2018Leptons","muons"),
    leptonsFromDileptons = cms.bool(True),
    lrWeightProducer = cms.PSet(
        Lambda = cms.int32(16000),
        calculate = cms.bool(False),
        doingElectrons = cms.bool(False),
        doingLR = cms.bool(True),
        interference = cms.int32(1),
        src = cms.InputTag("prunedMCLeptons")
    ),
    prefireWeights = cms.InputTag("prefiringweight","nonPrefiringProb"),
    pu_src = cms.InputTag("slimmedAddPileupInfo"),
    pu_weights = cms.vstring(),
    useMadgraphWeight = cms.bool(True),
    useTTBarWeight = cms.bool(False),
    use_bs_and_pv = cms.bool(True),
    usekFactor = cms.bool(True),
    vertex_src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    year = cms.int32(2018)
)


process.SimpleNtupler = cms.EDAnalyzer("SimpleNtupler_miniAOD",
    TriggerResults_src = cms.InputTag("TriggerResults","","PAT"),
    beamspot_src = cms.InputTag("offlineBeamSpot"),
    dimu_src = cms.InputTag("Our2018MuonsPlusMuonsMinus"),
    doElectrons = cms.bool(False),
    genEventInfo = cms.untracked.InputTag("generator"),
    hardInteraction = cms.PSet(
        allowFakeResonance = cms.bool(True),
        doingElectrons = cms.bool(False),
        matchTaus = cms.bool(True),
        resonanceIds = cms.vint32(32, 23, 39, 5000039),
        shutUp = cms.bool(True),
        src = cms.InputTag("prunedGenParticles")
    ),
    jet_src = cms.InputTag("slimmedJets"),
    metFilter = cms.VInputTag(cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), cms.InputTag("Flag_globalTightHalo2016Filter")),
    met_src = cms.InputTag("slimmedMETs"),
    vertices_src = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring(
        'FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'
    ),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1000)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring(
        'warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'
    ),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    PrescaleToCommon = cms.PSet(
        initialSeed = cms.untracked.uint32(1219)
    ),
    PrescaleToCommonMiniAOD = cms.PSet(
        initialSeed = cms.untracked.uint32(1219)
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zp2mu_histos.root')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring(
        'HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER'
    )
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(-1),
    version = cms.string('Parabolic')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(-1)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelQualityFromDbRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        )
    )
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGainRcd')
        ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )
    ),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiStripDetVOffRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )
    ),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('94X_mc2017_realistic_v17'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HcalTimeSlewEP = cms.ESSource("HcalTimeSlewEP",
    appendToDataLabel = cms.string('HBHE'),
    timeSlewParametersM2 = cms.VPSet(
        cms.PSet(
            slope = cms.double(-3.178648),
            tmax = cms.double(16.0),
            tzero = cms.double(23.960177)
        ), 
        cms.PSet(
            slope = cms.double(-1.556668),
            tmax = cms.double(10.0),
            tzero = cms.double(13.307784)
        ), 
        cms.PSet(
            slope = cms.double(-1.075824),
            tmax = cms.double(6.25),
            tzero = cms.double(9.109694)
        )
    ),
    timeSlewParametersM3 = cms.VPSet(
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(15.5),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-3.2),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(32.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        )
    )
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(20.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(
            0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
            0.058939, 0.125497
        ),
        HFdepthOneParameterB = cms.vdouble(
            -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
            0.000425, 0.000209
        ),
        HFdepthTwoParameterA = cms.vdouble(
            0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
            0.051579, 0.086593
        ),
        HFdepthTwoParameterB = cms.vdouble(
            -2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
            0.000157, -3e-06
        )
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(36000)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )
    ),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(8)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(206),
        zsThreshold = cms.int32(16)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(9)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(206),
        zsThreshold = cms.int32(16)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201),
        zsThreshold = cms.int32(24)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(False),
    useHEUpgrade = cms.bool(False),
    useHFUpgrade = cms.bool(False),
    useHOUpgrade = cms.bool(True),
    useIeta18depth1 = cms.bool(True),
    useLayer0Weight = cms.bool(False)
)


process.prefer("es_hardcode")

process.egmGsfElectronIDTask = cms.Task(process.egmGsfElectronIDs, process.electronMVAValueMapProducer, process.electronMVAVariableHelper, process.electronRegressionValueMapProducer)


process.Zprime2muAnalysisSequence_MiniAOD = cms.Sequence(process.muonPhotonMatchMiniAOD+process.leptonsMini+process.allDimuons+process.dimuons)


process.Zprime2muAnalysisSequence = cms.Sequence(process.muonPhotonMatch+process.leptons+process.allDimuons+process.dimuons)


process.egmGsfElectronIDSequence = cms.Sequence(cms.Task(process.heepIDVarValueMaps), process.egmGsfElectronIDTask)


process.pathOur2017 = cms.Path(process.primaryVertexMiniAOD+process.EventCounter+process.dileptonPreseletor+process.muonPhotonMatchMiniAOD+process.egmGsfElectronIDSequence+process.Our2017Leptons+process.allOur2017MuonsPlusMuonsMinus+process.Our2017MuonsPlusMuonsMinus+process.Our2017MuonsPlusMuonsMinusHistos)


process.pathOur2018 = cms.Path(process.primaryVertexMiniAOD+process.EventCounter+process.dileptonPreseletor+process.muonPhotonMatchMiniAOD+process.egmGsfElectronIDSequence+process.Our2018Leptons+process.allOur2018MuonsPlusMuonsMinus+process.Our2018MuonsPlusMuonsMinus+process.Our2018MuonsPlusMuonsMinusHistos+process.prunedMCLeptons+process.SimpleNtupler)


process.pathOur2018MuPrescaled = cms.Path(process.primaryVertexMiniAOD+process.EventCounter+process.dileptonPreseletor+process.muonPhotonMatchMiniAOD+process.egmGsfElectronIDSequence+process.Our2018MuPrescaledLeptons+process.allOur2018MuPrescaledMuonsPlusMuonsMinus+process.Our2018MuPrescaledMuonsPlusMuonsMinus+process.Our2018MuPrescaledMuonsPlusMuonsMinusHistos)


process.pathOur2017MuPrescaledCommon = cms.Path(process.PrescaleToCommon+process.primaryVertexMiniAOD+process.EventCounter+process.dileptonPreseletor+process.muonPhotonMatchMiniAOD+process.egmGsfElectronIDSequence+process.Our2017MuPrescaledCommonLeptons+process.allOur2017MuPrescaledCommonMuonsPlusMuonsMinus+process.Our2017MuPrescaledCommonMuonsPlusMuonsMinus+process.Our2017MuPrescaledCommonMuonsPlusMuonsMinusHistos)


process.pathOur2018MuPrescaledCommon = cms.Path(process.PrescaleToCommon+process.primaryVertexMiniAOD+process.EventCounter+process.dileptonPreseletor+process.muonPhotonMatchMiniAOD+process.egmGsfElectronIDSequence+process.Our2018MuPrescaledCommonLeptons+process.allOur2018MuPrescaledCommonMuonsPlusMuonsMinus+process.Our2018MuPrescaledCommonMuonsPlusMuonsMinus+process.Our2018MuPrescaledCommonMuonsPlusMuonsMinusHistos)


process.pathOur2017MuPrescaled = cms.Path(process.primaryVertexMiniAOD+process.EventCounter+process.dileptonPreseletor+process.muonPhotonMatchMiniAOD+process.egmGsfElectronIDSequence+process.Our2017MuPrescaledLeptons+process.allOur2017MuPrescaledMuonsPlusMuonsMinus+process.Our2017MuPrescaledMuonsPlusMuonsMinus+process.Our2017MuPrescaledMuonsPlusMuonsMinusHistos)


