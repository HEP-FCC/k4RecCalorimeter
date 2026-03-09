// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/ICaloCellConstantsSvc.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for data related to calorimeter cells.
 */


#ifndef RECCALOCOMMON_ICALOCELLCONSTANTSSVC_H
#define RECCALOCOMMON_ICALOCELLCONSTANTSSVC_H


#include "GaudiKernel/IInterface.h"
#include <typeinfo>
#include <string>
#include <memory>
#include <any>
#include <iostream>


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
class ICaloCellConstantsSvc
  : virtual public IInterface
{
public:
  DeclareInterfaceID (ICaloCellConstantsSvc, 1, 0);


  /**
   * @brief Retrieve an object from the store.
   * @param key Key of the object.
   *
   * Retrieve an object of type @c T from the store that was recorded
   * with key @c key.  Returns a pointer to the object, or @c nullptr
   * if either the key wasn't found or the type doesn't match
   * what was recorded.
   */
  template <class T>
  const T* getObj (const std::string& key) const
  {
    const std::any* a = getAnyObj (key);
    if (a) return std::any_cast<T> (a);
    return nullptr;
  }


  /**
   * @brief Record an object in the store.
   * @param key Key of the object.
   * @param obj The object to record, passed via move.
   *
   * Returns true if the object was successfully recorded, false if not
   * (an object was already recorded with the same key, for example).
   */
  template <class T>
  bool putObj (const std::string& key, T&& obj)
  {
    return putAnyObj (key, std::any (std::move (obj)));
  }


  /**
   * @brief Retrieve an object from the store, as a @c std::any.
   * @param key Key of the object.
   *
   * Retrieve an object of type @c T from the store that was recorded
   * with key @c key.  Returns a pointer to a @c std::any holding the
   * object, or @c nullptr if the key wasn't found.
   */
  virtual const std::any* getAnyObj (const std::string& key) const = 0;


  /**
   * @brief Record an object in the store, as a @c std::any.
   * @param key Key of the object.
   * @parma obj The object to record, passed via move.
   *
   * Returns true if the object was successfully recorded, false if not
   * (an object was already recorded with the same key, for example).
   */
  virtual bool putAnyObj (const std::string& key, std::any&& obj) = 0;
};


} // namespace k4::recCalo


#endif // not RECCALOCOMMON_ICALOCELLCONSTANTSSVC_H
