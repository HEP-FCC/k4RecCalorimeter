#include "EventCounter.h"

DECLARE_COMPONENT(EventCounter)

EventCounter::EventCounter(const std::string& name, ISvcLocator* svcLoc) : Gaudi::Algorithm(name, svcLoc) {
}

StatusCode EventCounter::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure()) return sc;
  m_counter = 0;
  return StatusCode::SUCCESS;
}

StatusCode EventCounter::execute(const EventContext& ctx ) const {
  int evt = ctx.evt();
  // if ( ( m_frequency <= 1 ) || ( ( m_counter ) % m_frequency == 0 ) ) {
  if ( ( m_frequency <= 1 ) || ( evt % m_frequency == 0 ) ) {
    info() << "Processing event " << evt << endmsg;
    //info() << "Counter is " << m_counter << endmsg;
  }
  m_counter++;
  return StatusCode::SUCCESS;
}

StatusCode EventCounter::finalize() {
  info() << "Processed " << m_counter << " events" << endmsg;
  return Gaudi::Algorithm::finalize();
}
