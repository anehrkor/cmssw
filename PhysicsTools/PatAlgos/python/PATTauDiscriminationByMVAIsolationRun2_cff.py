import FWCore.ParameterSet.Config as cms

discriminationByIsolationMVArun2v1raw = cms.EDProducer("PATTauDiscriminationByMVAIsolationRun2",

    # tau collection to discriminate
    PFTauProducer = cms.InputTag('pfTauProducer'), # what is pfTauProducer? should this be replaced by slimmedTaus?

    # Require leading pion ensures that:
    #  1) these is at least one track above threshold (0.5 GeV) in the signal cone
    #  2) a track OR a pi-zero in the signal cone has pT > 5 GeV
    Prediscriminants = requireLeadTrack,
    loadMVAfromDB = cms.bool(True),
    mvaName = cms.string("tauIdMVAnewDMwLT"),
    mvaOpt = cms.string("newDMwLT"),
    
    srcChargedIsoPtSum = cms.InputTag('chargedIsoPtSum'),
    srcNeutralIsoPtSum = cms.InputTag('neutralIsoPtSum'),
    srcPUcorrPtSum = cms.InputTag('puCorrPtSum'),
    srcPhotonPtSumOutsideSignalCone = cms.InputTag('photonPtSumOutsideSignalCone'),
    srcFootprintCorrection = cms.InputTag('footprintCorrection') 
)

#discriminationByIsolationMVArun2v1VLoose = recoTauDiscriminantCutMultiplexer.clone(
#    PFTauProducer = cms.InputTag('pfTauProducer'),    
#    Prediscriminants = requireLeadTrack,
#    toMultiplex = cms.InputTag('discriminationByIsolationMVArun2v1raw'),
#    key = cms.InputTag('discriminationByIsolationMVArun2v1raw:category'),
#    loadMVAfromDB = cms.bool(True),
#    mapping = cms.VPSet(
#        cms.PSet(
#            category = cms.uint32(0),
#            cut = cms.string("newDMwLTEff80"),
#            variable = cms.string("pt"),
#        )
#    )
#)
#discriminationByIsolationMVArun2v1Loose = discriminationByIsolationMVArun2v1VLoose.clone()
#discriminationByIsolationMVArun2v1Loose.mapping[0].cut = cms.string("newDMwLTEff70")
#discriminationByIsolationMVArun2v1Medium = discriminationByIsolationMVArun2v1VLoose.clone()
#discriminationByIsolationMVArun2v1Medium.mapping[0].cut = cms.string("newDMwLTEff60")
#discriminationByIsolationMVArun2v1Tight = discriminationByIsolationMVArun2v1VLoose.clone()
#discriminationByIsolationMVArun2v1Tight.mapping[0].cut = cms.string("newDMwLTEff50")
#discriminationByIsolationMVArun2v1VTight = discriminationByIsolationMVArun2v1VLoose.clone()
#discriminationByIsolationMVArun2v1VTight.mapping[0].cut = cms.string("newDMwLTEff40")

mvaIsolation2SeqRun2 = cms.Sequence(
   discriminationByIsolationMVArun2v1raw
#   + discriminationByIsolationMVArun2v1VLoose
#   + discriminationByIsolationMVArun2v1Loose
#   + discriminationByIsolationMVArun2v1Medium
#   + discriminationByIsolationMVArun2v1Tight
#   + discriminationByIsolationMVArun2v1VTight
)