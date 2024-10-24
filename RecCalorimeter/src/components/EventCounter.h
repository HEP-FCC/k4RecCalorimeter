#ifndef RECCALORIMETER_EVENTCOUNTER_H
#define RECCALORIMETER_EVENTCOUNTER_H

// Gaudi
#include "Gaudi/Algorithm.h"

class EventCounter : public Gaudi::Algorithm {
 public:
  EventCounter(const std::string& name, ISvcLocator* svcLoc);
  StatusCode initialize();
  StatusCode execute(const EventContext& ctx) const;
  StatusCode finalize();
 protected:
  Gaudi::Property<int> m_frequency{ this, "Frequency", 1, "How often to print the event number" };
  mutable int                  m_counter = 0;  // needed or not? to be checked
};

#endif
