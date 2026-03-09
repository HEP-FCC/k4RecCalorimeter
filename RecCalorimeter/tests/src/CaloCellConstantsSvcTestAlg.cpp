/**
 * @file RecCalorimeter/tests/src/CaloCellConstantsSvcTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Test for CaloCellConstantsSvc
 */


#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "k4FWCore/GaudiChecks.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"


namespace k4::recCalo {


class CaloCellConstantsSvcTestAlg
  : public Algorithm
{
public:
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  ServiceHandle<ICaloCellConstantsSvc> m_svc
  { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc", "" };
};


DECLARE_COMPONENT(k4::recCalo::CaloCellConstantsSvcTestAlg);


StatusCode CaloCellConstantsSvcTestAlg::initialize()
{
  K4_GAUDI_CHECK( m_svc.retrieve() );

  using payload_t = std::vector<int>;
  K4_GAUDI_CHECK (!m_svc->getObj<payload_t> ("test"));
  payload_t v1 {1, 2, 3};
  payload_t v2 (v1);
  K4_GAUDI_CHECK (m_svc->putObj ("test", std::move(v2)));
  K4_GAUDI_CHECK (v2.empty());
  const payload_t* v3 = nullptr;
  K4_GAUDI_CHECK ((v3 = m_svc->getObj<payload_t> ("test")));
  K4_GAUDI_CHECK (*v3 == v1);
  K4_GAUDI_CHECK (!m_svc->getObj<payload_t> ("test2"));
  K4_GAUDI_CHECK (!m_svc->getObj<std::vector<float> > ("test1"));
  K4_GAUDI_CHECK (!m_svc->putObj ("test", std::move(v2)));
  K4_GAUDI_CHECK ((v3 = m_svc->getObj<payload_t> ("test")));
  K4_GAUDI_CHECK (*v3 == v1);

  return StatusCode::SUCCESS;
}


StatusCode CaloCellConstantsSvcTestAlg::execute()
{
  return StatusCode::SUCCESS;
}


} // namespace k4::recCalo
