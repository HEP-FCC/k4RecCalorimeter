#include "CreateCaloClustersSlidingWindowFCCee.h"

// Gaudi
#include "GaudiKernel/PhysicalConstants.h"

// datamodel
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/Vector3f.h"

DECLARE_COMPONENT(CreateCaloClustersSlidingWindowFCCee)

CreateCaloClustersSlidingWindowFCCee::CreateCaloClustersSlidingWindowFCCee(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("clusters", m_clusters, "Handle for calo clusters (output collection)");
  declareProperty("clusterCells", m_clusterCells, "Handle for calo cluster cells (output collection)");
  declareProperty("towerTool", m_towerTool, "Handle for the tower building tool");
}

StatusCode CreateCaloClustersSlidingWindowFCCee::initialize() {
  if (Gaudi::Algorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }

  // retrieve tool that builds towers of calorimeter cells
  if (!m_towerTool.retrieve()) {
    error() << "Unable to retrieve the tower building tool." << endmsg;
    return StatusCode::FAILURE;
  }

  // get number of towers in theta and phi
  m_towerTool->towersNumber(m_nThetaTower, m_nPhiTower);
  debug() << "Number of calorimeter towers (theta x phi) : " << m_nThetaTower << " x " << m_nPhiTower << endmsg;

  // make sure that the number of towers in theta is larger than the seeding sliding window
  if (m_nThetaTower < m_nThetaWindow) {
    debug() << "Window size in theta too small!!! Window " << m_nThetaWindow << " # of theta towers " << m_nThetaTower
            << endmsg;
    m_nThetaTower = m_nThetaWindow;
  }

  info() << "CreateCaloClustersSlidingWindowFCCee initialized" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CreateCaloClustersSlidingWindowFCCee::execute(const EventContext&) const {

  // 1. Create calorimeter towers (calorimeter grid in theta phi, all layers merged)
  m_towers.assign(m_nThetaTower, std::vector<float>(m_nPhiTower, 0));
  // Create an output cluster collection
  auto clusters = m_clusters.createAndPut();
  edm4hep::CalorimeterHitCollection* clusterCells = nullptr;
  // If the user wants a new cell collection with only the clustered cells, then create it
  // If not, links will point from clusters to cells in their initial collections
  if (m_attachCells) {
    clusterCells = m_clusterCells.createAndPut();
  }
  // Build towers
  if (m_towerTool->buildTowers(m_towers, true) == 0) {
    debug() << "Empty cell collection." << endmsg;
    return StatusCode::SUCCESS;
  }

  // 2. Find local maxima with sliding window, build preclusters, calculate their barycentre position
  // calculate the sum of first m_nThetaWindow bins in theta, for each phi tower
  std::vector<float> sumOverTheta(m_nPhiTower, 0);
  for (int iTheta = 0; iTheta < m_nThetaWindow; iTheta++) {
    std::transform(sumOverTheta.begin(), sumOverTheta.end(), m_towers[iTheta].begin(), sumOverTheta.begin(),
                   std::plus<float>());
  }

  // preclusters with phi, theta weighted position and transverse energy
  m_preClusters.clear();
  int halfThetaPos = floor(m_nThetaPosition / 2.);
  int halfPhiPos = floor(m_nPhiPosition / 2.);
  float posX = 0;
  float posY = 0;
  float posZ = 0;
  float sumEnergyPos = 0;
  std::map<std::pair<uint, uint>, std::vector<edm4hep::CalorimeterHit>> cellsInTowersMap = m_towerTool->cellsInTowers();

  // final cluster window
  int halfThetaFin = floor(m_nThetaFinal / 2.);
  int halfPhiFin = floor(m_nPhiFinal / 2.);
  double sumEnergyFin = 0.;
  int idThetaFin = 0;
  int idPhiFin = 0;

  // loop over all Theta slices starting at the half of the first window
  int halfThetaWin = floor(m_nThetaWindow / 2.);
  int halfPhiWin = floor(m_nPhiWindow / 2.);
  float sumWindow = 0;
  float sumPhiSlicePrevThetaWin = 0;
  float sumPhiSliceNextThetaWin = 0;
  float sumFirstPhiSlice = 0;
  float sumLastPhiSlice = 0;
  bool toRemove = false;

  // one slice in theta = window, now only sum over window elements in phi
  for (int iTheta = halfThetaWin; iTheta < m_nThetaTower - halfThetaWin; iTheta++) {
    // sum the first window in phi
    for (int iPhiWindow = -halfPhiWin; iPhiWindow <= halfPhiWin; iPhiWindow++) {
      sumWindow += sumOverTheta[phiNeighbour(iPhiWindow)];
    }
    // loop over all the phi slices
    for (int iPhi = 0; iPhi < m_nPhiTower; iPhi++) {
      // if energy is above threshold, it may be a precluster
      if (sumWindow > m_energyThreshold) {
        // test local maximum in phi
        // check closest neighbour on the right
        if (sumOverTheta[phiNeighbour(iPhi - halfPhiWin)] < sumOverTheta[phiNeighbour(iPhi + halfPhiWin + 1)]) {
          toRemove = true;
        }
        // check closest neighbour on the left
        if (sumOverTheta[phiNeighbour(iPhi + halfPhiWin)] < sumOverTheta[phiNeighbour(iPhi - halfPhiWin - 1)]) {
          toRemove = true;
        }
        // test local maximum in theta
        // check closest neighbour on the right (if it is not the first window)
        if (iTheta > halfThetaWin) {
          for (int iPhiWindowLocalCheck = iPhi - halfPhiWin; iPhiWindowLocalCheck <= iPhi + halfPhiWin;
               iPhiWindowLocalCheck++) {
            sumPhiSlicePrevThetaWin += m_towers[iTheta - halfThetaWin - 1][phiNeighbour(iPhiWindowLocalCheck)];
            sumLastPhiSlice += m_towers[iTheta + halfThetaWin][phiNeighbour(iPhiWindowLocalCheck)];
          }
          if (sumPhiSlicePrevThetaWin > sumLastPhiSlice) {
            toRemove = true;
          }
        }
        // check closest neighbour on the left (if it is not the last window)
        if (iTheta < m_nThetaTower - halfThetaWin - 1) {
          for (int iPhiWindowLocalCheck = iPhi - halfPhiWin; iPhiWindowLocalCheck <= iPhi + halfPhiWin;
               iPhiWindowLocalCheck++) {
            sumPhiSliceNextThetaWin += m_towers[iTheta + halfThetaWin + 1][phiNeighbour(iPhiWindowLocalCheck)];
            sumFirstPhiSlice += m_towers[iTheta - halfThetaWin][phiNeighbour(iPhiWindowLocalCheck)];
          }
          if (sumPhiSliceNextThetaWin > sumFirstPhiSlice) {
            toRemove = true;
          }
        }
        sumFirstPhiSlice = 0;
        sumLastPhiSlice = 0;
        sumPhiSlicePrevThetaWin = 0;
        sumPhiSliceNextThetaWin = 0;
        // Build precluster
        if (!toRemove) {
          // Calculate barycentre position (usually smaller window used to reduce noise influence)
          posX = 0;
          posY = 0;
          posZ = 0;
          sumEnergyPos = 0;
          // weighted mean for position in theta and phi
          for (int ipTheta = iTheta - halfThetaPos; ipTheta <= iTheta + halfThetaPos; ipTheta++) {
            for (int ipPhi = iPhi - halfPhiPos; ipPhi <= iPhi + halfPhiPos; ipPhi++) {
              for (auto cell : cellsInTowersMap[std::make_pair(ipTheta, phiNeighbour(ipPhi))]) {
                posX += cell.getPosition().x * cell.getEnergy();
                posY += cell.getPosition().y * cell.getEnergy();
                posZ += cell.getPosition().z * cell.getEnergy();
                sumEnergyPos += cell.getEnergy();
              }
            }
          }
          // If too small energy in the position window, calculate the position in the whole sliding window
          // Assigns correct position for cases with maximum energy deposits close to the border in theta
          if (sumEnergyPos > m_energyThresholdFraction * m_energyThreshold) {
            posX /= sumEnergyPos;
            posY /= sumEnergyPos;
            posZ /= sumEnergyPos;
          } else {
            posX = 0;
            posY = 0;
            posZ = 0;
            sumEnergyPos = 0;
            for (int ipTheta = iTheta - halfThetaWin; ipTheta <= iTheta + halfThetaWin; ipTheta++) {
              for (int ipPhi = iPhi - halfPhiWin; ipPhi <= iPhi + halfPhiWin; ipPhi++) {
                for (auto cell : cellsInTowersMap[std::make_pair(ipTheta, phiNeighbour(ipPhi))]) {
                  posX += cell.getPosition().x * cell.getEnergy();
                  posY += cell.getPosition().y * cell.getEnergy();
                  posZ += cell.getPosition().z * cell.getEnergy();
                  sumEnergyPos += cell.getEnergy();
                }
              }
            }
            posX /= sumEnergyPos;
            posY /= sumEnergyPos;
            posZ /= sumEnergyPos;
          }
          // Calculate final cluster energy
          sumEnergyFin = 0;
          // Final cluster position
          idThetaFin = m_towerTool->idTheta(atan2(sqrt(posX * posX + posY * posY), posZ));
          idPhiFin = m_towerTool->idPhi(atan2(posY, posX));
          // Recalculating the energy within the final cluster size
          for (int ipTheta = idThetaFin - halfThetaFin; ipTheta <= idThetaFin + halfThetaFin; ipTheta++) {
            for (int ipPhi = idPhiFin - halfPhiFin; ipPhi <= idPhiFin + halfPhiFin; ipPhi++) {
              if (ipTheta >= 0 && ipTheta < m_nThetaTower) { // check if we are not outside of map in theta
                if (m_ellipseFinalCluster) {
                  if (pow((ipTheta - idThetaFin) / (m_nThetaFinal / 2.), 2) +
                          pow((ipPhi - idPhiFin) / (m_nPhiFinal / 2.), 2) <
                      1) {
                    sumEnergyFin += m_towers[ipTheta][phiNeighbour(ipPhi)];
                  }
                } else {
                  sumEnergyFin += m_towers[ipTheta][phiNeighbour(ipPhi)];
                }
              }
            }
          }
          // check if changing the barycentre did not decrease energy below threshold
          if (sumEnergyFin > m_energyThreshold) {
            precluster newPreCluster;
            newPreCluster.X = posX;
            newPreCluster.Y = posY;
            newPreCluster.Z = posZ;
            newPreCluster.theta = atan2(sqrt(posX * posX + posY * posY), posZ);
            newPreCluster.phi = atan2(posY, posX);
            newPreCluster.transEnergy = sumEnergyFin;
            m_preClusters.push_back(newPreCluster);
          }
        }
      }
      toRemove = false;
      // finish processing that window in phi, shift window to the next phi tower
      // substract first phi tower in current window
      sumWindow -= sumOverTheta[phiNeighbour(iPhi - halfPhiWin)];
      // add next phi tower to the window
      sumWindow += sumOverTheta[phiNeighbour(iPhi + halfPhiWin + 1)];
    }
    sumWindow = 0;
    // finish processing that slice, shift window to next theta tower
    if (iTheta < m_nThetaTower - halfThetaWin - 1) {
      // substract first theta slice in current window
      std::transform(sumOverTheta.begin(), sumOverTheta.end(), m_towers[iTheta - halfThetaWin].begin(),
                     sumOverTheta.begin(), std::minus<float>());
      // add next theta slice to the window
      std::transform(sumOverTheta.begin(), sumOverTheta.end(), m_towers[iTheta + halfThetaWin + 1].begin(),
                     sumOverTheta.begin(), std::plus<float>());
    }
  }

  debug() << "Pre-clusters size before duplicates removal: " << m_preClusters.size() << endmsg;

  // 4. Sort the preclusters according to the transverse energy (descending)
  std::sort(m_preClusters.begin(), m_preClusters.end(),
            [](precluster clu1, precluster clu2) { return clu1.transEnergy > clu2.transEnergy; });

  // 5. Remove duplicates
  for (auto it1 = m_preClusters.begin(); it1 != m_preClusters.end(); it1++) {
    // loop over all clusters with energy lower than it1 (sorting), erase if too close
    for (auto it2 = it1 + 1; it2 != m_preClusters.end();) {
      if ((abs(int(m_towerTool->idTheta((*it1).theta) - m_towerTool->idTheta((*it2).theta))) < m_nThetaDuplicates) &&
          ((abs(int(m_towerTool->idPhi((*it1).phi) - m_towerTool->idPhi((*it2).phi))) < m_nPhiDuplicates) ||
           (abs(int(m_towerTool->idPhi((*it1).phi) - m_towerTool->idPhi((*it2).phi))) >
            m_nPhiTower - m_nPhiDuplicates))) {
        m_preClusters.erase(it2);
      } else {
        it2++;
      }
    }
  }
  debug() << "Pre-clusters size after duplicates removal: " << m_preClusters.size() << endmsg;

  // 6. Create final clusters
  for (const auto clu : m_preClusters) {
    float clusterEnergy = clu.transEnergy / sin(clu.theta);
    // apply energy sharing correction (if flag set to true)
    if (m_energySharingCorrection) {
      int idThetaCl = m_towerTool->idTheta(clu.theta);
      int idPhiCl = m_towerTool->idPhi(clu.phi);
      // sum of energies in other clusters in each theta-phi tower of our current cluster (idThetaCl, idPhiCl)
      std::vector<std::vector<float>> sumEnergySharing;
      sumEnergySharing.assign(m_nThetaFinal, std::vector<float>(m_nPhiFinal, 0));
      // loop over all clusters and check if they have any tower in common with our current cluster
      for (const auto cluSharing : m_preClusters) {
        int idThetaClShare = m_towerTool->idTheta(cluSharing.theta);
        int idPhiClShare = m_towerTool->idPhi(cluSharing.phi);
        if (idThetaCl != idThetaClShare && idPhiCl != idPhiClShare) {
          // check for overlap between clusters
          if (abs(idThetaClShare - idThetaCl) < m_nThetaFinal &&
              ((abs(idPhiClShare - idPhiCl) < m_nPhiFinal) ||
               (abs(idPhiClShare - idPhiCl) > m_nPhiTower - m_nPhiFinal))) {
            // add energy in shared towers to sumEnergySharing[][]
            for (int iTheta = std::max(idThetaCl, idThetaClShare) - halfThetaFin;
                 iTheta <= std::min(idThetaCl, idThetaClShare) + halfThetaFin; iTheta++) {
              for (int iPhi = std::max(idPhiCl, idPhiClShare) - halfPhiFin;
                   iTheta <= std::min(idPhiCl, idPhiClShare) + halfPhiFin; iPhi++) {
                if (iTheta >= 0 && iTheta < m_nThetaTower) { // check if we are not outside of map in theta
                  sumEnergySharing[iTheta - idThetaCl + halfThetaFin][phiNeighbour(iPhi - idPhiCl + halfPhiFin)] +=
                      m_towers[iTheta][phiNeighbour(iPhi)] / sin(m_towerTool->theta(iTheta));
                }
              }
            }
          }
        }
      }
      // apply the actual correction: substract the weighted energy contributions in other clusters
      for (int iTheta = idThetaCl - halfThetaFin; iTheta <= idThetaCl + halfThetaFin; iTheta++) {
        for (int iPhi = idPhiCl - halfPhiFin; iPhi <= idPhiCl + halfPhiFin; iPhi++) {
          if (iTheta - idThetaCl + halfThetaFin >= 0)
            if (sumEnergySharing[iTheta - idThetaCl + halfThetaFin][phiNeighbour(iPhi - idPhiCl + halfPhiFin)] != 0) {
              float sumButOne =
                  sumEnergySharing[iTheta - idThetaCl + halfThetaFin][phiNeighbour(iPhi - idPhiCl + halfPhiFin)];
              float towerEnergy = m_towers[iTheta][phiNeighbour(iPhi)] / sin(m_towerTool->theta(iTheta));
              clusterEnergy -= towerEnergy * sumButOne / (sumButOne + towerEnergy);
            }
        }
      }
    }

    // save the clusters in our EDM
    // check ET threshold once more (ET could change with the energy sharing correction)
    if (clusterEnergy * sin(clu.theta) > m_energyThreshold) {
      auto cluster = clusters->create();
      cluster.setPosition(edm4hep::Vector3f(clu.X, clu.Y, clu.Z));
      cluster.setEnergy(clusterEnergy);
      debug() << "Attaching cells to the clusters." << endmsg;
      m_towerTool->attachCells(clu.theta, clu.phi, halfThetaFin, halfPhiFin, cluster, clusterCells,
                               m_ellipseFinalCluster);
      debug() << "Cluster theta: " << clu.theta << " phi: " << clu.phi << " x: " << cluster.getPosition().x
              << " y: " << cluster.getPosition().y << " z: " << cluster.getPosition().z
              << " energy: " << cluster.getEnergy() << " contains: " << cluster.hits_size() << " cells" << endmsg;
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode CreateCaloClustersSlidingWindowFCCee::finalize() { return Gaudi::Algorithm::finalize(); }

unsigned int CreateCaloClustersSlidingWindowFCCee::phiNeighbour(int aIPhi) const {
  if (aIPhi < 0) {
    return m_nPhiTower + aIPhi;
  } else if (aIPhi >= m_nPhiTower) {
    return aIPhi % m_nPhiTower;
  }
  return aIPhi;
}
