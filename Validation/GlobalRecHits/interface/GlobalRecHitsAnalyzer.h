#ifndef GlobalRecHitsAnalyzer_h
#define GlobalRecHitsAnalyzer_h

/** \class GlobalHitsProducer
 *  
 *  Class to fill PGlobalRecHit object to be inserted into data stream 
 *  containing information about various sub-systems in global coordinates 
 *  with full geometry
 *
 *  $Date: 2008/03/13 21:20:31 $
 *  $Revision: 1.6 $
 *  \author M. Strang SUNY-Buffalo
 */

// framework & common header files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

//DQM services
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Provenance/interface/Provenance.h"
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// ecal calorimeter info
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

// hcal calorimeter info
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HFDataFrame.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalSourcePositionData.h"

// silicon strip info
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h" 
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h" 
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h" 
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h" 
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h" 
#include "DataFormats/Common/interface/OwnVector.h" 

// silicon pixel info
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

// muon DT info
#include "DataFormats/DTDigi/interface/DTDigi.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Validation/DTRecHits/interface/DTHitQualityUtils.h"

// muon CSC info
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigi.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"

// muon RPC info
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"

// event info
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h" 

// general info 
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h" 
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h" 

#include <iostream>
#include <stdlib.h>
#include <string>
#include <memory>
#include <vector>
#include <map>
#include <math.h>

#include "TString.h"
#include "DQMServices/Core/interface/MonitorElement.h"


class GlobalRecHitsAnalyzer : public edm::EDAnalyzer
{

 public:

  typedef std::map<uint32_t,float,std::less<uint32_t> > MapType;

  explicit GlobalRecHitsAnalyzer(const edm::ParameterSet&);
  virtual ~GlobalRecHitsAnalyzer();
  virtual void beginJob(const edm::EventSetup&);
  virtual void endJob();  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
 private:

  // production related methods
  void fillECal(const edm::Event&, const edm::EventSetup&);
  //void storeECal(PGlobalRecHit&);
  void fillHCal(const edm::Event&, const edm::EventSetup&);
  //void storeHCal(PGlobalRecHit&);
  void fillTrk(const edm::Event&, const edm::EventSetup&);
  //void storeTrk(PGlobalRecHit&);
  void fillMuon(const edm::Event&, const edm::EventSetup&);
  //void storeMuon(PGlobalRecHit&);  

  //void clear();

 private:

  //  parameter information
  std::string fName;
  int verbosity;
  int frequency;
  std::string label;
  bool getAllProvenances;
  bool printProvenanceInfo;

  DQMStore *dbe;

  // Electromagnetic info
  // ECal info
 
  MonitorElement *mehEcaln[3];
  MonitorElement *mehEcalRes[3];

  edm::InputTag ECalEBSrc_;
  edm::InputTag ECalUncalEBSrc_;
  edm::InputTag ECalEESrc_;
  edm::InputTag ECalUncalEESrc_;
  edm::InputTag ECalESSrc_;

  // HCal info

  MonitorElement *mehHcaln[4];
  MonitorElement *mehHcalRes[4];

  edm::InputTag HCalSrc_;

  // Tracker info
  // SiStrip
  
  MonitorElement *mehSiStripn[19];
  MonitorElement *mehSiStripResX[19];
  MonitorElement *mehSiStripResY[19];

  edm::InputTag SiStripSrc_;

  std::vector<PSimHit> matched;
  std::pair<LocalPoint,LocalVector> 
    projectHit( const PSimHit& hit,
		const StripGeomDetUnit* stripDet,
		const BoundPlane& plane);
  edm::ParameterSet conf_;

  // SiPxl

  MonitorElement *mehSiPixeln[7];
  MonitorElement *mehSiPixelResX[7];
  MonitorElement *mehSiPixelResY[7];

  edm::InputTag SiPxlSrc_;

  // Muon info
  // DT

  MonitorElement *mehDtMuonn;
  MonitorElement *mehCSCn;
  MonitorElement *mehRPCn;
  MonitorElement *mehDtMuonRes;
  MonitorElement *mehCSCResRDPhi;
  MonitorElement *mehRPCResX;

  edm::InputTag MuDTSrc_;
  edm::InputTag MuDTSimSrc_;

  // Return a map between DTRecHit1DPair and wireId
  std::map<DTWireId, std::vector<DTRecHit1DPair> >
    map1DRecHitsPerWire(const DTRecHitCollection* dt1DRecHitPairs);
  
  // Compute SimHit distance from wire (cm)
  float simHitDistFromWire(const DTLayer* layer,
			   DTWireId wireId,
			   const PSimHit& hit);
  
  // Find the RecHit closest to the muon SimHit
  template  <typename type>
    const type* 
    findBestRecHit(const DTLayer* layer,
		   DTWireId wireId,
		   const std::vector<type>& recHits,
		   const float simHitDist);
  
  // Compute the distance from wire (cm) of a hits in a DTRecHit1DPair
  float recHitDistFromWire(const DTRecHit1DPair& hitPair, 
			   const DTLayer* layer);
  // Compute the distance from wire (cm) of a hits in a DTRecHit1D
  float recHitDistFromWire(const DTRecHit1D& recHit, const DTLayer* layer);
    
  // Does the real job
  template  <typename type>
    int compute(const DTGeometry *dtGeom,
		 std::map<DTWireId, std::vector<PSimHit> > simHitsPerWire,
		 std::map<DTWireId, std::vector<type> > recHitsPerWire,
		 int step);

  // CSC
  //Defined above....

  edm::InputTag MuCSCSrc_;

  std::map<int, edm::PSimHitContainer> theMap;
  void plotResolution(const PSimHit &simHit, const CSCRecHit2D &recHit,
		      const CSCLayer *layer, int chamberType);

  // RPC

  //Defined above...

  edm::InputTag MuRPCSrc_;
  edm::InputTag MuRPCSimSrc_;

  // private statistics information
  unsigned int count;

}; // end class declaration

#endif

#ifndef GlobalHitMap
#define GlobalHitMap

// geometry mapping
static const int dTrk             = 1;
static const int sdPxlBrl         = 1;
static const int sdPxlFwd         = 2;
static const int sdSiTIB          = 3;
static const int sdSiTID          = 4;
static const int sdSiTOB          = 5;
static const int sdSiTEC          = 6;

static const int dMuon            = 2;
static const int sdMuonDT         = 1;
static const int sdMuonCSC        = 2;
static const int sdMuonRPC        = 3;
static const int sdMuonRPCRgnBrl  = 0;
static const int sdMuonRPCRgnFwdp = 1;
static const int sdMuonRPCRgnFwdn = -1;

static const int dEcal            = 3;
static const int sdEcalBrl        = 1;
static const int sdEcalFwd        = 2;
static const int sdEcalPS         = 3;
static const int sdEcalTT         = 4;
static const int sdEcalLPnD       = 5;

static const int dHcal            = 4;
static const int sdHcalEmpty      = 0;
static const int sdHcalBrl        = 1;
static const int sdHcalEC         = 2;
static const int sdHcalOut        = 3;
static const int sdHcalFwd        = 4;
static const int sdHcalTT         = 5;
static const int sdHcalCalib      = 6;
static const int sdHcalCompst     = 7;

#endif //PGlobalRecHitsProducer_h
