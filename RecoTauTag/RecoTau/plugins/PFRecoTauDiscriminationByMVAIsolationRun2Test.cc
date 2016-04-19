#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminationByMVAIsolationRun2T.h"

typedef RecoTauDiscriminationByMVAIsolationRun2T<reco::PFTauRef,reco::PFTauDiscriminator,PFTauDiscriminationProducerBase> PFRecoTauDiscriminationByMVAIsolationRun2Test;

DEFINE_FWK_MODULE(PFRecoTauDiscriminationByMVAIsolationRun2Test);
