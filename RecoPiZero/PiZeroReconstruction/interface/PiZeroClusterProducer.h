#ifndef RecoPiZero_PiZeroReconstruction_PiZeroClusterProducer_h
#define RecoPiZero_PiZeroReconstruction_PiZeroClusterProducer_h

#include <vector>

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"

#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"

#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

class PiZeroClusterProducer {
public:

	PiZeroClusterProducer(edm::Handle<EBRecHitCollection> ebHandle, edm::Handle<EERecHitCollection> eeHandle, edm::Handle<ESRecHitCollection> esHandle,
			CaloSubdetectorTopology *estopology, const EcalPreshowerGeometry *esGeometry, const CaloGeometry* geometry,
			double SeedEnergyEB=0, double SeedEnergyEE=0, double SeedEnergyES=0, double barrelClusterPt=0, double endcapClusterPt=0);
	~PiZeroClusterProducer();

	// main functions to build clusters from rec hits
	void buildEBClusters(std::vector<reco::CaloCluster> &ebclusters, const EcalChannelStatus &channelstatus);
	void buildEEClusters(std::vector<reco::CaloCluster> &eseeclusters, std::vector<reco::CaloCluster> &eseeclusters_tot, const EcalChannelStatus &channelstatus);

	// getter & setter
	void setSeedEnergyThresholds(double thresholdEB=0, double thresholdEE=0, double thresholdES=0);
	void setClusterPtThresholds(double thresholdEB=0, double thresholdEE=0);

	// helper functions
	bool checkStatusOfEcalRecHit(const EcalChannelStatus &channelstatus, const EcalRecHit &rechit);
	void convxtalid(int &nphi, int &neta);
	int diff_neta_s(int neta1, int neta2);
	int diff_nphi_s(int nphi1, int nphi2);
	reco::PreshowerCluster makeOnePreshowerCluster(int stripwindow, ESDetId *strip);
	void findESRoad(int stripwindow, ESDetId strip, EcalPreshowerNavigator theESNav, int plane);

private:

    // event handles
    edm::Handle<EBRecHitCollection> ebHandle_;
    edm::Handle<EERecHitCollection> eeHandle_;
    edm::Handle<ESRecHitCollection> esHandle_;

    // geometry variables
    CaloTopology *ebtopology_;
    CaloTopology *eetopology_;
    CaloSubdetectorTopology *estopology_;
    const EcalPreshowerGeometry *esGeometry_;
    const CaloGeometry* geometry_;

    // todo: needed here?
    std::vector<ESDetId> esroad_2d;

    // values obtained from Luca's code (PreshowerTools)
    const double preshowerMIP = 8.108e-05;
    const double preshowerGamma = 0.024;
    const double preshowerCalibPlaneX = 1.0;
    const double preshowerCalibPlaneY = 0.7;
    const int preshowerClusterWindowSize = 15;

    // shower shape variables - what do they even mean?
    bool param_LogWeighted_;
    float param_T0_barl_;
    float param_T0_endc_;
    float param_T0_endcES_;
    float param_W0_;
    float param_X0_;

    // threshold cuts
    double SeedEnergyEB_;
    double SeedEnergyEE_;
    double SeedEnergyES_;
    double barrelClusterPt_;
    double endcapClusterPt_;
};

#endif
