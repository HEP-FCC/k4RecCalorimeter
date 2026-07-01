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

/** detailed string for calo type */
std::ostream& operator<<(std::ostream& os, const CHT& cht) {
  os << " calo hit type: ";

  switch (cht.caloType()) {
    case CHT::em:
      os << " em,   ";
      break;
    case CHT::had:
      os << " had,  ";
      break;
    case CHT::muon:
      os << " muon, ";
      break;
    default:
      os << "  -  ,";
  }
  switch (cht.caloID()) {
    case CHT::ecal:
      os << "ecal,  ";
      break;
    case CHT::hcal:
      os << "hcal,  ";
      break;
    case CHT::yoke:
      os << "yoke,  ";
      break;
    case CHT::lcal:
      os << "lcal,  ";
      break;
    case CHT::lhcal:
      os << "lhcal, ";
      break;
    case CHT::bcal:
      os << "bcal,  ";
      break;
    default:
      os << "  -  ,";
  }
  switch (cht.layout()) {
    case CHT::any:
      os << "any,    ";
      break;
    case CHT::ring:
      os << "ring,   ";
      break;
    case CHT::endcap:
      os << "endcap, ";
      break;
    case CHT::barrel:
      os << "barrel, ";
      break;
    case CHT::plug:
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
    return CHT::ring;
  if (str.find("plug") != std::string::npos)
    return CHT::plug;
  if (str.find("endcap") != std::string::npos)
    return CHT::endcap;
  if (str.find("barrel") != std::string::npos)
    return CHT::barrel;

  std::cout << " not found :" << str << " in " << name << std::endl;
  return CHT::any;
}

CHT::CaloID caloIDFromString(const std::string& name) {
  std::string str(name);
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);

  if (str.find("ecal") != std::string::npos)
    return CHT::ecal;
  if (str.find("hcal") != std::string::npos)
    return CHT::hcal;
  if (str.find("yoke") != std::string::npos)
    return CHT::yoke;
  if (str.find("lcal") != std::string::npos)
    return CHT::lcal;
  if (str.find("lhcal") != std::string::npos)
    return CHT::lhcal;
  if (str.find("bcal") != std::string::npos)
    return CHT::bcal;

  return CHT::unknown;
}

CHT::CaloType caloTypeFromString(const std::string& name) {
  std::string str(name);
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);

  if (str.find("em") != std::string::npos)
    return CHT::em;
  if (str.find("had") != std::string::npos)
    return CHT::had;
  if (str.find("muon") != std::string::npos)
    return CHT::muon;

  // jl: this should probably also have a separate "unknown" or "any" value?
  return CHT::em;
}
