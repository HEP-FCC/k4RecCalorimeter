/**
 * @file RecCalorimeter/src/components/CaloCellConstantsSvc.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for data related to calorimeter cells.
 */

#include "CaloCellConstantsSvc.h"


DECLARE_COMPONENT(k4::recCalo::CaloCellConstantsSvc);


namespace k4::recCalo {


/**
 * @brief Retrieve an object from the store, as a @c std::any.
 * @param key Key of the object.
 */
const std::any* CaloCellConstantsSvc::getAnyObj (const std::string& key) const
{
  std::lock_guard lock (m_mutex);
  auto it = m_objs.find (key);
  if (it != m_objs.end()) return &it->second;
  return nullptr;
}


/**
 * @brief Record an object in the store, as a @c std::any.
 * @param key Key of the object.
 * @param obj The object to record, passed via move.
 */
bool CaloCellConstantsSvc::putAnyObj (const std::string& key, std::any&& obj)
{
  std::lock_guard lock (m_mutex);
  return m_objs.try_emplace (key, std::move(obj)).second;
}


} // namespace k4::recCalo


