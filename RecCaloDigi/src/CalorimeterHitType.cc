/*
 * Copyright (c) 2020-2024 Key4hep-Project.
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
#include "CalorimeterHitType.h"

#include <algorithm>
#include <string>

/** detailed string for calo type */
std::ostream& operator<<(std::ostream& os, const CHT& cht) {
  os << " calo hit type: ";

  switch (cht.caloType()) {
  case CHT::CaloType::em:
    os << " em,   ";
    break;
  case CHT::CaloType::had:
    os << " had,  ";
    break;
  case CHT::CaloType::muon:
    os << " muon, ";
    break;
  default:
    os << "  -  ,";
  }
  switch (cht.caloID()) {
  case CHT::CaloID::ecal:
    os << "ecal,  ";
    break;
  case CHT::CaloID::hcal:
    os << "hcal,  ";
    break;
  case CHT::CaloID::yoke:
    os << "yoke,  ";
    break;
  case CHT::CaloID::lcal:
    os << "lcal,  ";
    break;
  case CHT::CaloID::lhcal:
    os << "lhcal, ";
    break;
  case CHT::CaloID::bcal:
    os << "bcal,  ";
    break;
  default:
    os << "  -  ,";
  }
  switch (cht.layout()) {
  case CHT::Layout::any:
    os << "any,    ";
    break;
  case CHT::Layout::ring:
    os << "ring,   ";
    break;
  case CHT::Layout::endcap:
    os << "endcap, ";
    break;
  case CHT::Layout::barrel:
    os << "barrel, ";
    break;
  case CHT::Layout::plug:
    os << "plug,   ";
    break;
  default:
    os << "  -  ,";
  }
  os << " layer: " << cht.layer();

  return os;
}

/** Helper functions that should go to Marlinutil/CalorimeterHitTypes.hh */

CHT::Layout layoutFromString(const std::string& name) {
  std::string str(name);
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);

  if (str.find("ring") != std::string::npos)
    return CHT::Layout::ring;
  if (str.find("plug") != std::string::npos)
    return CHT::Layout::plug;
  if (str.find("endcap") != std::string::npos)
    return CHT::Layout::endcap;
  if (str.find("barrel") != std::string::npos)
    return CHT::Layout::barrel;

  return CHT::Layout::any;
}

CHT::CaloID caloIDFromString(const std::string& name) {
  std::string str(name);
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);

  if (str.find("ecal") != std::string::npos)
    return CHT::CaloID::ecal;
  if (str.find("hcal") != std::string::npos)
    return CHT::CaloID::hcal;
  if (str.find("yoke") != std::string::npos)
    return CHT::CaloID::yoke;
  if (str.find("lcal") != std::string::npos)
    return CHT::CaloID::lcal;
  if (str.find("lhcal") != std::string::npos)
    return CHT::CaloID::lhcal;
  if (str.find("bcal") != std::string::npos)
    return CHT::CaloID::bcal;

  return CHT::CaloID::unknown;
}

CHT::CaloType caloTypeFromString(const std::string& name) {
  std::string str(name);
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);

  if (str.find("em") != std::string::npos)
    return CHT::CaloType::em;
  if (str.find("had") != std::string::npos)
    return CHT::CaloType::had;
  if (str.find("muon") != std::string::npos)
    return CHT::CaloType::muon;

  return CHT::CaloType::unknown;
}
