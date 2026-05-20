/*
 * Copyright (c) 2014-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef RECCALOCOMMON_ICELLPOSITIONSTOOL_H
#define RECCALOCOMMON_ICELLPOSITIONSTOOL_H

// Gaudi
#include "GaudiKernel/IAlgTool.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDSegmentation/BitFieldCoder.h" // CellID

namespace edm4hep {
class CalorimeterHit;
class CalorimeterHitCollection;
} // namespace edm4hep

/** @class ICellPositionsTool
 *
 *  Abstract interface to calorimeter cell positions.
 *
 *  @author Coralie Neubueser
 *  @date   2018-01
 */

namespace k4::recCalo {

class ICellPositionsTool : virtual public IAlgTool {
public:
  using CellID = dd4hep::DDSegmentation::CellID;

  DeclareInterfaceID(ICellPositionsTool, 1, 0);

  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection& outputColl) const = 0;

  virtual dd4hep::Position xyzPosition(const CellID aCellId) const = 0;
  virtual int layerId(const CellID aCellId) const = 0;
};

} // namespace k4::recCalo

#endif /* RECCALOCOMMON_ICELLPOSITIONSTOOL_H */
