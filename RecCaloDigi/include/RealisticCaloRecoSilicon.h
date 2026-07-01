#ifndef REALISTICCALORECOSILICON_H
#define REALISTICCALORECOSILICON_H 1

#include "RealisticCaloReco.h"

/** === RealisticCaloRecoSilicon Processor === <br>
    realistic reconstruction of silicon calorimeter hits
    D.Jeans 02/2016.

    24 March 2016: removed gap corrections - to be put into separate processor

*/

struct RealisticCaloRecoSilicon final : RealisticCaloReco {

 public:
  RealisticCaloRecoSilicon(const std::string& name, ISvcLocator* svcLoc);

 protected:
  float reconstructEnergy(const edm4hep::CalorimeterHit* hit, int layer) const override;
} ;

#endif 
