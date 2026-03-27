#ifndef RECINTERFACE_INOISECONSTTOOL_H
#define RECINTERFACE_INOISECONSTTOOL_H

// from Gaudi
#include "GaudiKernel/IAlgTool.h"

/** @class INoiseConstTool
 *
 *  Abstract interface for tools that return noise values (RMS and offset) per calorimeter cell
 *  @author Coralie Neubueser
 *  @date   2018-01
 */

class INoiseConstTool : virtual public IAlgTool {
public:
  DeclareInterfaceID(INoiseConstTool, 1, 0);

  /** Expected noise per cell in terms of sigma of Gaussian distibution.
   *   @param[in] aCellId of the cell of interest.
   *   return double.
   */
  virtual double getNoiseRMSPerCell(uint64_t aCellID) const = 0;

  /** Expected noise per cell in terms of mean of distibution.
   *   @param[in] aCellId of the cell of interest.
   *   return double.
   */
  virtual double getNoiseOffsetPerCell(uint64_t aCellID) const = 0;

  /** Expected noise per cell.
   *   @param[in] aCellId of the cell of interest.
   *   return [rms, offset]
   */
  virtual std::pair<double, double> getNoisePerCell(uint64_t aCellID) const = 0;
};

#endif /* RECINTERFACE_INOISECONSTTOOL_H */
