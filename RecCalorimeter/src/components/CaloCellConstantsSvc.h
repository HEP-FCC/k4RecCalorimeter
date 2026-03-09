// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCalorimeter/src/components/CaloCellConstantsSvc.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for data related to calorimeter cells.
 */


#ifndef RECCALORIMETER_CALOCELLCONSTANTSSVC_H
#define RECCALORIMETER_CALOCELLCONSTANTSSVC_H


#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "GaudiKernel/Service.h"
#include <mutex>


namespace k4::recCalo {


/**
 * @brief Holder for data related to calorimeter cells.
 *
 * Numerous tools deal with constants defined for each calorimeter cell.
 * For highly-granular calorimeters with many cells, such as the Allegro
 * ECal, this can consume a significant amount of memory.  But then if one
 * creates multiple instances of these tools --- perhaps to explore
 * the effects of different configurations --- then these data
 * may be duplicated, effectively wasting the memory allocated by the
 * duplicates.  This service provides a simple key/value registry allowing
 * the data to be stored outside of the tools themselves, and thus
 * providing the possibility of sharing.
 *
 * This can be thought of as supplying something akin to the ATLAS
 * detector store.  In principle, it would be nice to store such objects
 * the same manner as event data.  But podio is not so suitable, as it
 * does not allow for storing arbitrary types, and in general does not
 * really seem intended for the storage of data which is not meant
 * to be persistent.
 */
class CaloCellConstantsSvc
  : public extends<Service, ICaloCellConstantsSvc>
{
public:
  using base_class::base_class;


  /**
   * @brief Retrieve an object from the store, as a @c std::any.
   * @param key Key of the object.
   *
   * Retrieve an object of type @c T from the store that was recorded
   * with key @c key.  Returns a pointer to a @c std::any holding the
   * object, or @c nullptr if the key wasn't found.
   */
  virtual const std::any* getAnyObj (const std::string& key) const override;


  /**
   * @brief Record an object in the store, as a @c std::any.
   * @param key Key of the object.
   * @param obj The object to record, passed via move.
   *
   * Returns true if the object was successfully recorded, false if not
   * (an object was already recorded with the same key, for example).
   */
  virtual bool putAnyObj (const std::string& key, std::any&& obj) override;


private:
  /// The stored objects.
  std::map<std::string, std::any> m_objs;

  /// Guard access to the map.
  mutable std::mutex m_mutex;
};


} // namespace k4::recCalo


#endif // not RECCALORIMETER_CALOCELLCONSTANTSSVC_H
