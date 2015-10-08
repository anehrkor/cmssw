#include "RecoPiZero/PiZeroReconstruction/interface/PiZeroClusterProducer.h"
// it needs to be the hardcoded topologies for some reason...
// if the non hardcoded topologies are used we get segfaults.
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapHardcodedTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

PiZeroClusterProducer::PiZeroClusterProducer(edm::Handle<EBRecHitCollection> ebHandle, edm::Handle<EERecHitCollection> eeHandle, edm::Handle<ESRecHitCollection> esHandle,
		CaloSubdetectorTopology *estopology, const EcalPreshowerGeometry *esGeometry, const CaloGeometry* geometry,
		double SeedEnergyEB, double SeedEnergyEE, double SeedEnergyES, double barrelClusterPt, double endcapClusterPt)
{
	ebHandle_ = ebHandle;
	eeHandle_ = eeHandle;
	esHandle_ = esHandle;

    ebtopology_ = new CaloTopology();
    EcalBarrelHardcodedTopology* ebhtopology = new EcalBarrelHardcodedTopology();
    ebtopology_->setSubdetTopology(DetId::Ecal,EcalBarrel,ebhtopology);
    eetopology_ = new CaloTopology();
    EcalEndcapHardcodedTopology* eehtopology = new EcalEndcapHardcodedTopology();
    eetopology_->setSubdetTopology(DetId::Ecal,EcalEndcap,eehtopology);

    estopology_ = estopology;
    esGeometry_ = esGeometry;
    geometry_ = geometry;

    SeedEnergyEB_ = SeedEnergyEB;
    SeedEnergyEE_ = SeedEnergyEE;
    SeedEnergyES_ = SeedEnergyES;
    barrelClusterPt_ = barrelClusterPt;
    endcapClusterPt_ = endcapClusterPt;

    // set shower shape variables (from Luca's code)
    param_LogWeighted_ = true;
    param_T0_barl_     = 7.4;
    param_T0_endc_     = 3.1;
    param_T0_endcES_   = 1.2;
    param_W0_          = 4.2;
    param_X0_          = 0.89;
}

PiZeroClusterProducer::~PiZeroClusterProducer() {
	delete ebtopology_;
	delete eetopology_;
}

void
PiZeroClusterProducer::setSeedEnergyThresholds(double thresholdEB, double thresholdEE, double thresholdES) {
	SeedEnergyEB_ = thresholdEB;
	SeedEnergyEE_ = thresholdEE;
	SeedEnergyES_ = thresholdES;
}

void
PiZeroClusterProducer::setClusterPtThresholds(double thresholdEB, double thresholdEE) {
	barrelClusterPt_ = thresholdEB;
	endcapClusterPt_ = thresholdEE;
}

bool
PiZeroClusterProducer::checkStatusOfEcalRecHit(const EcalChannelStatus &channelstatus, const EcalRecHit &rechit) {
	int status = int(channelstatus[rechit.id().rawId()].getStatusCode());
	if(status > 0) return false;
	return true;
}

void
PiZeroClusterProducer::convxtalid(int &nphi, int &neta) {
	if(neta > 0) neta -= 1;
	if(nphi > 359) nphi = nphi - 360;
}

int
PiZeroClusterProducer::diff_neta_s(int neta1, int neta2) {
	return (neta1-neta2);
}

int
PiZeroClusterProducer::diff_nphi_s(int nphi1, int nphi2) {
	if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2)))
		return nphi1-nphi2;
	else
	{
		int mdiff = 360 - abs(nphi1-nphi2);
		if(nphi1 > nphi2)
			mdiff = -mdiff;
		return mdiff;
	}
}

void
PiZeroClusterProducer::findESRoad(int stripwindow, ESDetId strip, EcalPreshowerNavigator theESNav, int plane) {
	if(strip == ESDetId(0)) return;

	ESDetId next;
	theESNav.setHome(strip);

	esroad_2d.push_back(strip);

	if(plane == 1)
	{
		int neast = 0;
		while((next = theESNav.east()) != ESDetId(0) && next != strip)
		{
			esroad_2d.push_back(next);
			++neast;
			if(neast == stripwindow) break;
		}
		int nwest = 0;
		theESNav.home();
		while((next = theESNav.west()) != ESDetId(0) && next != strip)
		{
			esroad_2d.push_back(next);
			++nwest;
			if(nwest == stripwindow) break;
		}
	}
	else if(plane == 2)
	{
		int nnorth = 0;
		while((next = theESNav.north()) != ESDetId(0) && next != strip)
		{
			esroad_2d.push_back(next);
			++nnorth;
			if(nnorth == stripwindow) break;
		}
		int nsouth = 0;
		theESNav.home();
		while((next = theESNav.south()) != ESDetId(0) && next != strip)
		{
			esroad_2d.push_back(next);
			++nsouth;
			if(nsouth == stripwindow) break;
		}
	}

	theESNav.home();
}

reco::PreshowerCluster
PiZeroClusterProducer::makeOnePreshowerCluster(int stripwindow, ESDetId *strip) {
	reco::PreshowerCluster finalCluster;

	std::map<DetId, EcalRecHit> RecHitsMap;
	for(ESRecHitCollection::const_iterator itRecHit = esHandle_->begin(); itRecHit != esHandle_->end(); itRecHit++)
	{
		RecHitsMap.insert(std::make_pair(itRecHit->id(),*itRecHit));
	}

	std::set<DetId> used_strips;

	esroad_2d.clear();

	int plane = strip->plane();

	// Collection of cluster strips
	EcalRecHitCollection clusterRecHits;
	// Map of strips for position calculation
	std::map<DetId, EcalRecHit> RecHitsPos;

	EcalPreshowerNavigator navigator(*strip, estopology_);
	navigator.setHome(*strip);

	findESRoad(stripwindow,*strip,navigator,plane);

	if(plane == 1)
	{
		ESDetId strip_north = navigator.north();
		findESRoad(stripwindow,strip_north,navigator,plane);
		navigator.home();
		ESDetId strip_south = navigator.south();
		findESRoad(stripwindow,strip_south,navigator,plane);
		navigator.home();
	}
	if(plane == 2)
	{
		ESDetId strip_west = navigator.west();
		findESRoad(stripwindow,strip_west,navigator,plane);
		navigator.home();
		ESDetId strip_east = navigator.east();
		findESRoad(stripwindow,strip_east,navigator,plane);
		navigator.home();
	}

	// Start clustering from strip with max Energy in the road
	double Emax = 0;
	bool found = false;
	std::map<DetId, EcalRecHit>::const_iterator itmax;
	// loop over strips
	for(std::vector<ESDetId>::const_iterator itId = esroad_2d.begin(); itId != esroad_2d.end(); itId++)
	{
		std::map<DetId, EcalRecHit>::const_iterator itStrip = RecHitsMap.find(*itId);
		if((used_strips.find(itStrip->first) != used_strips.end()) || itStrip == RecHitsMap.end() || itStrip->second.energy() <= 0) continue;

		if(itStrip->second.energy() > Emax)
		{
			Emax = itStrip->second.energy();
			found = true;
			itmax = itStrip;
		}
	}

	if(!found)
		return finalCluster;

	// first, save the hottest strip
	clusterRecHits.push_back(itmax->second);
	RecHitsPos.insert(std::make_pair(itmax->first,itmax->second));
	used_strips.insert(itmax->first);

	// find position of adjacent strips
	ESDetId next, strip1, strip2;
	navigator.setHome(itmax->first);
	ESDetId startES = itmax->first;

	if(plane == 1)
	{
		int nadjacent_east = 0;
		while((next = navigator.east()) != ESDetId(0) && next != startES && nadjacent_east < 2)
		{
			++nadjacent_east;
			std::map<DetId, EcalRecHit>::const_iterator itStrip = RecHitsMap.find(next);
			if((used_strips.find(itStrip->first) != used_strips.end()) || itStrip == RecHitsMap.end() || itStrip->second.energy() <= 0) continue;
			clusterRecHits.push_back(itStrip->second);
			if(nadjacent_east == 1) strip1 = next;
			used_strips.insert(itStrip->first);
		}
		navigator.home();
		int nadjacent_west = 0;
		while((next = navigator.west()) != ESDetId(0) && next != startES && nadjacent_west < 2)
		{
			++nadjacent_west;
			std::map<DetId, EcalRecHit>::const_iterator itStrip = RecHitsMap.find(next);
			if((used_strips.find(itStrip->first) != used_strips.end()) || itStrip == RecHitsMap.end() || itStrip->second.energy() <= 0) continue;
			clusterRecHits.push_back(itStrip->second);
			if(nadjacent_west == 1) strip2 = next;
			used_strips.insert(itStrip->first);
		}
		navigator.home();
	}
	else if(plane == 2)
	{
		int nadjacent_north = 0;
		while((next = navigator.north()) != ESDetId(0) && next != startES && nadjacent_north < 2)
		{
			++nadjacent_north;
			std::map<DetId, EcalRecHit>::const_iterator itStrip = RecHitsMap.find(next);
			if((used_strips.find(itStrip->first) != used_strips.end()) || itStrip == RecHitsMap.end() || itStrip->second.energy() <= 0) continue;
			clusterRecHits.push_back(itStrip->second);
			if(nadjacent_north == 1) strip1 = next;
			used_strips.insert(itStrip->first);
		}
		navigator.home();
		int nadjacent_south = 0;
		while((next = navigator.south()) != ESDetId(0) && next != startES && nadjacent_south < 2)
		{
			++nadjacent_south;
			std::map<DetId, EcalRecHit>::const_iterator itStrip = RecHitsMap.find(next);
			if((used_strips.find(itStrip->first) != used_strips.end()) || itStrip == RecHitsMap.end() || itStrip->second.energy() <= 0) continue;
			clusterRecHits.push_back(itStrip->second);
			if(nadjacent_south == 1) strip2 = next;
			used_strips.insert(itStrip->first);
		}
	}
	else
		return finalCluster;

	if(strip1 != ESDetId(0))
		RecHitsPos.insert(std::make_pair(RecHitsMap.find(strip1)->first,RecHitsMap.find(strip1)->second));
	if(strip2 != ESDetId(0))
		RecHitsPos.insert(std::make_pair(RecHitsMap.find(strip2)->first,RecHitsMap.find(strip2)->second));

	double energyPos = 0;
	double xPos(0), yPos(0), zPos(0);

	for(std::map<DetId, EcalRecHit>::const_iterator itPos = RecHitsPos.begin(); itPos != RecHitsPos.end(); itPos++)
	{
		energyPos += itPos->second.energy();
		GlobalPoint position = geometry_->getPosition(itPos->first);
		xPos += itPos->second.energy() * position.x();
		yPos += itPos->second.energy() * position.y();
		zPos += itPos->second.energy() * position.z();
	}

	if(energyPos > 0)
	{
		xPos /= energyPos;
		yPos /= energyPos;
		zPos /= energyPos;
	}

	double energyCluster = 0;
	for(EcalRecHitCollection::const_iterator itCluster = clusterRecHits.begin(); itCluster != clusterRecHits.end(); itCluster++)
	{
		energyCluster += itCluster->energy();
	}

	std::vector<std::pair<DetId, float> > usedHits;
	reco::PreshowerCluster output(energyCluster, math::XYZPoint(xPos,yPos,zPos), usedHits, plane);

	return output;
}

void
PiZeroClusterProducer::buildEBClusters(std::vector<reco::CaloCluster> &ebclusters, const EcalChannelStatus &channelstatus) {
	std::vector<EcalRecHit> ebseeds;
	std::set<EBDetId> detIdsUsed;

	// find seeds in ECAL using threshold cut on crystal energy
	for(EBRecHitCollection::const_iterator itRecHit = ebHandle_->begin(); itRecHit != ebHandle_->end(); itRecHit++)
	{
		if(itRecHit->energy() > SeedEnergyEB_)
			ebseeds.push_back(*itRecHit);
	}

	// sort seeds by energy
	sort(ebseeds.begin(),ebseeds.end(),EcalRecHitLess());

	// loop over seeds and make clusters
	for(std::vector<EcalRecHit>::const_iterator itSeed = ebseeds.begin(); itSeed != ebseeds.end(); itSeed++)
	{

		EBDetId seed_id = itSeed->id();

		// if seed already used, skip it
		if(detIdsUsed.count(seed_id) != 0) continue;

		// find 3x3 matrix of crystals
		std::vector<DetId> vCluster = ebtopology_->getWindow(seed_id,3,3);
		// used in position calculation
		std::vector<std::pair<DetId,double> > clustersUsed;

		// crystals actually used after removing those already used
		std::vector<const EcalRecHit*> RecHitsIn3x3;

		double simpleEnergy = 0;
		double positionEnergy = 0;

		// actually make the 3x3 clusters
		for(std::vector<DetId>::const_iterator itDetId = vCluster.begin(); itDetId != vCluster.end(); itDetId++)
		{
			EBDetId thisId = (*itDetId);
			// skip crystal if already used
			if(detIdsUsed.count(thisId) != 0) continue;

			// find rec hit. if none found, skip
			EBRecHitCollection::const_iterator itCrystal = ebHandle_->find(thisId);
			if(itCrystal == ebHandle_->end()) continue;
			RecHitsIn3x3.push_back(&(*itCrystal));
			clustersUsed.push_back(std::make_pair(*itDetId,1.));

			simpleEnergy += itCrystal->energy();
			if(itCrystal->energy() > 0) positionEnergy += itCrystal->energy();
		}

		if(simpleEnergy <= 0) continue;

		int seed_ieta = seed_id.ieta();
		int seed_iphi = seed_id.iphi();
		convxtalid(seed_iphi,seed_ieta);

		double energy3x3 = 0;
		std::vector<std::pair<DetId,float> > energyFractions;

		double xCluster(0), yCluster(0), zCluster(0);
		double clusterWeight = 0;

		// calculate shower depth (don't really understand this part)
		float maxDepth = param_X0_ * (param_T0_barl_ + log(positionEnergy));
		float maxToFront;
		const CaloCellGeometry* cell = geometry_->getGeometry(seed_id);
		GlobalPoint position = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(0.);
		maxToFront = position.mag();

		bool allRecHitsGood = true;

		// compute energy and position of cluster
		for(std::vector<const EcalRecHit*>::const_iterator itRecHits3x3 = RecHitsIn3x3.begin(); itRecHits3x3 != RecHitsIn3x3.end(); itRecHits3x3++)
		{
			EBDetId thisId = (*itRecHits3x3)->id();
			if(!checkStatusOfEcalRecHit(channelstatus,*(*itRecHits3x3)))
				allRecHitsGood = false;

			int ieta = thisId.ieta();
			int iphi = thisId.iphi();
			convxtalid(iphi,ieta);

			double energy = (*itRecHits3x3)->energy();// todo: calibration coefficient not available without Luca's code ...
			int dx = diff_neta_s(seed_iphi,iphi);
			int dy = diff_nphi_s(seed_ieta,ieta);

			if(abs(dx) <= 1 && abs(dy) <= 1)
			{
				energy3x3 += energy;
				energyFractions.push_back(std::make_pair(thisId, energy));
				detIdsUsed.insert(thisId);
			}

			if(energy > 0)
			{
				double weight = std::max((double)0, param_W0_ + log(energy / positionEnergy));
				double position_geo;
				const CaloCellGeometry* cell = geometry_->getGeometry(thisId);
				GlobalPoint position = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(0.);
				position_geo = position.mag();
				double depth = maxDepth + maxToFront - position_geo;
				GlobalPoint thisPosition = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);

				xCluster = weight * thisPosition.x();
				yCluster = weight * thisPosition.y();
				zCluster = weight * thisPosition.z();
				clusterWeight += weight;
			}
		}

		if(!allRecHitsGood) continue;

		math::XYZPoint clusterPosition(xCluster/clusterWeight,yCluster/clusterWeight,zCluster/clusterWeight);

		// todo: containment corrections not available without Luca's code ...

		double clusterPt = energy3x3/cosh(clusterPosition.eta());
		if(clusterPt < barrelClusterPt_) continue;

		ebclusters.push_back(reco::CaloCluster(energy3x3, clusterPosition, reco::CaloID(reco::CaloID::DET_ECAL_BARREL), energyFractions, reco::CaloCluster::undefined, seed_id));
	}
}

void
PiZeroClusterProducer::buildEEClusters(std::vector<reco::CaloCluster> &eseeclusters, std::vector<reco::CaloCluster> &eseeclusters_tot, const EcalChannelStatus &channelstatus) {
	std::vector<EcalRecHit> eeseeds;
	std::set<EEDetId> detIdsUsed;
	std::vector<reco::CaloCluster> eeclusters;

	// find seeds in ECAL using threshold cut on crystal energy
	for(EERecHitCollection::const_iterator itRecHit = eeHandle_->begin(); itRecHit != eeHandle_->end(); itRecHit++)
	{
		if(itRecHit->energy() > SeedEnergyEE_)
			eeseeds.push_back(*itRecHit);
	}

	// sort seeds by energy
	sort(eeseeds.begin(),eeseeds.end(),EcalRecHitLess());

	// loop over seeds and make clusters
	for(std::vector<EcalRecHit>::const_iterator itSeed = eeseeds.begin(); itSeed != eeseeds.end(); itSeed++)
	{
		EEDetId seed_id = itSeed->id();

		// if seed already used, skip it
		if(detIdsUsed.count(seed_id) != 0) continue;

		// find 3x3 matrix of crystals
		std::vector<DetId> vCluster = eetopology_->getWindow(seed_id,3,3);
		// used in position calculation
		std::vector<std::pair<DetId,double> > clustersUsed;

		// crystals actually used after removing those already used
		std::vector<const EcalRecHit*> RecHitsIn3x3;

		double simpleEnergy = 0;
		double positionEnergy = 0;

		// actually make the 3x3 clusters
		for(std::vector<DetId>::const_iterator itDetId = vCluster.begin(); itDetId != vCluster.end(); itDetId++)
		{
			EEDetId thisId = (*itDetId);
			// skip crystal if already used
			if(detIdsUsed.count(thisId) != 0) continue;

			// find rec hit. if none found, skip
			EERecHitCollection::const_iterator itCrystal = eeHandle_->find(thisId);
			if(itCrystal == eeHandle_->end()) continue;
			RecHitsIn3x3.push_back(&(*itCrystal));
			clustersUsed.push_back(std::make_pair(*itDetId,1.));

			simpleEnergy += itCrystal->energy();
			if(itCrystal->energy() > 0) positionEnergy += itCrystal->energy();
		}

		if(simpleEnergy <= 0) continue;

		int seed_ix = seed_id.ix();
		int seed_iy = seed_id.iy();

		double energy3x3 = 0;
		std::vector<std::pair<DetId,float> > energyFractions;

		double xCluster(0), yCluster(0), zCluster(0);
		double clusterWeight = 0;

		// calculate shower depth (don't really understand this part)
		float maxDepth = param_X0_ * (param_T0_barl_ + log(positionEnergy));
		float maxToFront;
		const CaloCellGeometry* cell = geometry_->getGeometry(seed_id);
		GlobalPoint position = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(0.);
		maxToFront = position.mag();

		bool allRecHitsGood = true;

		// compute energy and position of cluster
		for(std::vector<const EcalRecHit*>::const_iterator itRecHits3x3 = RecHitsIn3x3.begin(); itRecHits3x3 != RecHitsIn3x3.end(); itRecHits3x3++)
		{
			EEDetId thisId = (*itRecHits3x3)->id();
			if(!checkStatusOfEcalRecHit(channelstatus,*(*itRecHits3x3)))
				allRecHitsGood = false;

			int ix = thisId.ix();
			int iy = thisId.iy();

			double energy = (*itRecHits3x3)->energy();// todo: calibration coefficient not available without Luca's code ...
			int dx = seed_ix - ix;
			int dy = seed_iy - iy;

			if(abs(dx) <= 1 && abs(dy) <= 1)
			{
				energy3x3 += energy;
				energyFractions.push_back(std::make_pair(thisId,energy));
				detIdsUsed.insert(thisId);
			}

			if(energy > 0)
			{
				double weight = std::max((double)0, param_W0_ + log(energy / positionEnergy));
				double position_geo;
				const CaloCellGeometry* cell = geometry_->getGeometry(thisId);
				GlobalPoint position = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(0.);
				position_geo = position.mag();
				double depth = maxDepth + maxToFront - position_geo;
				GlobalPoint thisPosition = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);

				xCluster = weight * thisPosition.x();
				yCluster = weight * thisPosition.y();
				zCluster = weight * thisPosition.z();
				clusterWeight += weight;
			}
		}

		if(!allRecHitsGood) continue;

		math::XYZPoint clusterPosition(xCluster/clusterWeight,yCluster/clusterWeight,zCluster/clusterWeight);

		// todo: containment corrections not available without Luca's code ...

		double clusterPt = energy3x3/cosh(clusterPosition.eta());
		if(clusterPt < endcapClusterPt_) continue;

		eeclusters.push_back(reco::CaloCluster(energy3x3, clusterPosition, reco::CaloID(reco::CaloID::DET_ECAL_ENDCAP), energyFractions, reco::CaloCluster::undefined, seed_id));
	}

	// loop over eeclusters to find matches with preshower
	for(std::vector<reco::CaloCluster>::const_iterator iteeclusters = eeclusters.begin(); iteeclusters != eeclusters.end(); iteeclusters++)
	{
		if(fabs(iteeclusters->position().Eta()) > 1.7 && fabs(iteeclusters->position().Eta()) < 2.55)
		{
			const GlobalPoint point(iteeclusters->position().x(),iteeclusters->position().y(),iteeclusters->position().z());

			// preshower consists of two orthogonal layers
			DetId tmp1 = esGeometry_->getClosestCellInPlane(point,1);
			DetId tmp2 = esGeometry_->getClosestCellInPlane(point,2);

			if(tmp1.rawId() != 0 && tmp2.rawId() != 0)
			{
				ESDetId tmp1_conversion = tmp1;
				ESDetId tmp2_conversion = tmp2;

				reco::PreshowerCluster preshowerclusterp1 = makeOnePreshowerCluster(preshowerClusterWindowSize, &tmp1_conversion);
				reco::PreshowerCluster preshowerclusterp2 = makeOnePreshowerCluster(preshowerClusterWindowSize, &tmp2_conversion);

				double energy1 = preshowerclusterp1.energy();
				double energy2 = preshowerclusterp2.energy();
				double tempEnergy = iteeclusters->energy();
				// GeV to #MIPs
				energy1 = energy1 / preshowerMIP;
				energy2 = energy2 / preshowerMIP;

				if(energy1 > 2 || energy2 > 2)
				{
					double deltaE = preshowerGamma * (preshowerCalibPlaneX * energy1 + preshowerCalibPlaneY * energy2);
					tempEnergy += deltaE;

					// todo: containment corrections not available without Luca's code ...

					eseeclusters.push_back(reco::CaloCluster(tempEnergy, iteeclusters->position(), reco::CaloID(reco::CaloID::DET_ECAL_ENDCAP), iteeclusters->hitsAndFractions(), reco::CaloCluster::undefined, iteeclusters->seed()));

					double DZ2 = (preshowerclusterp2.z() - preshowerclusterp1.z()) / 2.;
					GlobalPoint positionCluster(preshowerclusterp1.x() * (1. + DZ2 / preshowerclusterp1.z()), preshowerclusterp2.y() * (1. + DZ2 / preshowerclusterp2.z()), (preshowerclusterp1.z() + preshowerclusterp2.z()) / 2.);

					if(fabs(preshowerclusterp1.z()) > 30 && fabs(preshowerclusterp2.z()) > 30)
					{
						math::XYZPoint positionClusterMath(positionCluster.x(),positionCluster.y(),positionCluster.z());
						eseeclusters_tot.push_back(reco::CaloCluster(tempEnergy, positionClusterMath, reco::CaloID(reco::CaloID::DET_ECAL_ENDCAP), iteeclusters->hitsAndFractions(), reco::CaloCluster::undefined, iteeclusters->seed()));
					}
				}
			}
		}
		else
			eseeclusters_tot.push_back(reco::CaloCluster(iteeclusters->energy(), iteeclusters->position(), reco::CaloID(reco::CaloID::DET_ECAL_ENDCAP), iteeclusters->hitsAndFractions(), reco::CaloCluster::undefined, iteeclusters->seed()));
	}
}
