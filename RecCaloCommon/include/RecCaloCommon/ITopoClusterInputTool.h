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
#ifndef RECCALOCOMMON_ITOPOCLUSTERINPUTTOOL_H
#define RECCALOCOMMON_ITOPOCLUSTERINPUTTOOL_H

#include "DDSegmentation/BitFieldCoder.h" // CellID

// Gaudi
#include "GaudiKernel/IAlgTool.h"

#include <unordered_map>

namespace edm4hep {
class CalorimeterHit;
}


namespace k4::recCalo {


/** @class ITopoClusterInputTool RecInterface/RecInterface/ITopoClusterInput.h ITopoClusterInputTool.h
 *
 *  Abstract interface to topo cluster input tool.
 *
 *  @author Coralie Neubueser
 */
class ITopoClusterInputTool : virtual public IAlgTool {
public:
  using CellID = dd4hep::DDSegmentation::CellID;

  DeclareInterfaceID(ITopoClusterInputTool, 1, 0);

  virtual StatusCode cellIDMap(std::unordered_map<CellID, double>& aCells) const = 0;
};


} // namespace k4::recCalo


#endif /* RECCALOCOMMON_ITOPOCLUSTERINPUTTOOL_H */
