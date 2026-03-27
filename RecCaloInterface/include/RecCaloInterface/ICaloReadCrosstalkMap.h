#ifndef RECINTERFACE_ICALOREADCROSSTALKMAP_H
#define RECINTERFACE_ICALOREADCROSSTALKMAP_H

// Gaudi
#include "GaudiKernel/IAlgTool.h"

/** @class ICaloReadCrosstalkMap k4Interface/include/RecCaloInterface/ICaloReadCrosstalkMap.h ICaloReadCrosstalkMap.h
 *
 *  Interface to the service reading the crosstalk map for ALLEGRO ECAL barrel.
 *
 *  @author Zhibo Wu
 */

class ICaloReadCrosstalkMap : virtual public IAlgTool {
public:
  DeclareInterfaceID(ICaloReadCrosstalkMap, 1, 0);

  virtual std::vector<uint64_t> const& getNeighbours(uint64_t cellID) const = 0;
  virtual std::vector<double> const& getCrosstalks(uint64_t cellID) const = 0;
};
#endif /* RECINTERFACE_ICALOREADCROSSTALKMAP_H */
