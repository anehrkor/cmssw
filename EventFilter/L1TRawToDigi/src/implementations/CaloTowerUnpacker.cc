#include "FWCore/Framework/interface/MakerMacros.h"

#include "EventFilter/L1TRawToDigi/interface/UnpackerFactory.h"

#include "L1TCollections.h"

namespace l1t {
   class CaloTowerUnpacker : public BaseUnpacker {
      public:
         virtual bool unpack(const unsigned block_id, const unsigned size, const unsigned char *data, UnpackerCollections *coll) override;
   };

   class CaloTowerUnpackerFactory : public BaseUnpackerFactory {
      public:
         virtual std::vector<UnpackerItem> create(const unsigned& fw, const int fedid) override;
   };
}

// Implementation

namespace l1t {
   bool
   CaloTowerUnpacker::unpack(const unsigned block_id, const unsigned size, const unsigned char *data, UnpackerCollections *coll)
   {

     LogDebug("L1T") << "Block ID  = " << block_id << " size = " << size;

     int nBX = int(ceil(size/44.)); // Since there are two Rx links per block with 2*28 slices in barrel and endcap + 2*13 for upgraded HF 

     // Find the first and last BXs
     int firstBX = -(std::ceil((double)nBX/2.)-1);
     int lastBX;
     if (nBX % 2 == 0) {
       lastBX = std::ceil((double)nBX/2.)+1;
     } else {
       lastBX = std::ceil((double)nBX/2.);
     }

     auto res_ = static_cast<L1TCollections*>(coll)->getTowers();
     res_->setBXRange(firstBX, lastBX);

     LogDebug("L1T") << "nBX = " << nBX << " first BX = " << firstBX << " lastBX = " << lastBX;

     // Initialise index
     int unsigned i = 0;

     // Link number is block_ID / 2
     unsigned link = block_id/2;
     
     // Also need link number rounded down to even number
     unsigned link_phi = (link % 2 == 0) ? link : (link -1);

     // Loop over multiple BX and fill towers collection
     for (int bx=firstBX; bx<lastBX; bx++){

       for (unsigned frame=1; frame<42 && frame<(size+1); frame++){

	 uint32_t raw_data = pop(data,i); // pop advances the index i internally

         if ((raw_data & 0xFFFF) != 0) {

           l1t::CaloTower tower1 = l1t::CaloTower();
    
           // First calo tower is in the LSW with phi
           tower1.setHwPt(raw_data & 0x1FF);
           tower1.setHwQual((raw_data >> 12) & 0xF);
           tower1.setHwEtRatio((raw_data >>9) & 0x7);
           tower1.setHwPhi(link_phi+1); // iPhi starts at 1
	 
           if (link % 2==0) { // Even number links carry Eta+
             tower1.setHwEta(frame); // iEta starts at 1
           } else { // Odd number links carry Eta-
             tower1.setHwEta(-1*frame);
           }
	 
           LogDebug("L1T") << "Tower 1: Eta " << tower1.hwEta() 
                           << " phi " << tower1.hwPhi() 
                           << " pT " << tower1.hwPt() 
                           << " frame " << frame 
                           << " qual " << tower1.hwQual() 
                           << " EtRatio " << tower1.hwEtRatio();

           res_->push_back(bx,tower1);
         }

         if (((raw_data >> 16)& 0xFFFF) != 0) {

           // Second calo tower is in the MSW with phi+1
           l1t::CaloTower tower2 = l1t::CaloTower();
	 
           tower2.setHwPt((raw_data >> 16) & 0x1FF);
           tower2.setHwQual((raw_data >> 28 ) & 0xF);
           tower2.setHwEtRatio((raw_data >> 25) & 0x7);
           tower2.setHwPhi(link_phi+2);

           if (link % 2==0) {
             tower2.setHwEta(frame);
           } else {
             tower2.setHwEta(-1*frame);
           }
	 
           LogDebug("L1T") << "Tower 2: Eta " << tower2.hwEta()
                           << " phi " << tower2.hwPhi()
                           << " pT " << tower2.hwPt()
                           << " frame " << frame
                           << " qual " << tower2.hwQual()
                           << " EtRatio " << tower2.hwEtRatio();

           res_->push_back(bx,tower2);
	 }
       }
     }
     
     return true;

  }

   std::vector<UnpackerItem>
   CaloTowerUnpackerFactory::create(const unsigned& fw, const int fedid)
   {

     // This unpacker is only appropriate for the Main Processor input (FED ID=2). Anything else should not be unpacked.
     
     if (fedid==2){

       std::vector<UnpackerItem> towersMap;
    
       // Map all even number links, which are Rx links and need unpacking to the same instance of the CaloTowerUnpacker
       // which receives the block_ID and can convert this to phi

       auto unpacker = std::shared_ptr<BaseUnpacker>(new CaloTowerUnpacker());

       for (int link = 0; link < 144; link++){
         if (link % 2 == 0) towersMap.push_back(std::make_pair(link, unpacker)); 
       }
     
       return towersMap;

     } else {
       
       return {};

     }

   };
};

DEFINE_L1TUNPACKER(l1t::CaloTowerUnpackerFactory);
