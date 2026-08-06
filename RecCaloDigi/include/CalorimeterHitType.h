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
#ifndef CalorimeterHitType_h
#define CalorimeterHitType_h 1

#include <iostream>

/** Helper class for decoding/encoding lcio::CalorimeterHit types for the ILD
 *  detector. The encoding is: caloType + 10 * caloID + 1000 * layout + 10000 * layerNum <br>
 *  (see enums: CaloType, CaloID and Layout for possible values).<br>
 *  Example usage: <br>
 *  <pre>
 *     lcio::CalorimeterHit* cHit = .... ;
 *
 *     // set the type (e.g. in digitization )
 *     cHit->setType( CHT( CHT::em ,  CHT::ecal , CHT::plug , 12 ) ) ;
 *
 *     ...
 *
 *     CHT cht = cHit->getType() ;
 *
 *     //   sum energies for electromagentic, hadronic and tailcatcher:
 *     if( cht.is( CHT::em ) )
 *          e_em +=  cHit->getEnergy() ;
 *     else
 *       if ( cht.is(CHT::had ) )
 *          e_had += cHit->getEnergy() ;
 *       else
 *          e_muon += cHit->getEnergy() ;
 *
 *     // use only EcalPlug hits:
 *     if( cht.is( CHT::ecal) && cht.is( CHT::plug) )
 *
 *     // get the layer number (e.g. for calibration or clustering)
 *     unsigned l = cht.layer() ;
 *     // or directly :
 *     unsigned l = CHT(  cHit->getType() ).layer()  ;
 *
 *     // detailed print:
 *     std::cout <<  CHT(  cHit->getType() ) << std::endl ;
 *
 *  </pre>
 *
 *  F.Gaede, DESY, 12/2008
 */

class CHT {
public:
  /** calorimeter types */
  enum CaloType { em = 0, had = 1, muon = 2 };

  /** calo ids - specific to ILD */
  enum CaloID { unknown = 0, ecal = 1, hcal = 2, yoke = 3, lcal = 4, lhcal = 5, bcal = 6 };

  /** calo layout / subdetector  */
  enum Layout { any = 0, barrel = 1, endcap = 2, plug = 3, ring = 4 };

  /** C'tor for initialization from CalorimeterHit::getType()  */
  CHT(int type) : m_type(type) {}

  /** C'tor  for encoding the calo type inforamtion  */
  CHT(CaloType c, CaloID n, Layout l, unsigned lay)
      : m_type(c * fCaloType + n * fCaloID + l * fLayout + lay * fLayer) {}

  /** calorimeter type: CHT::em , CHT::had, CHT::muon */
  CaloType caloType() const { return (CaloType)(m_type % fCaloID); }

  /** calo ID - see enum CaloID for allowed values */
  CaloID caloID() const { return (CaloID)((m_type % fLayout) / fCaloID); }

  /** calo layout - see enum layout for allowed values */
  Layout layout() const { return (Layout)((m_type % fLayer) / fLayout); }

  /** calo layer of hit  */
  unsigned layer() const { return unsigned(m_type) / fLayer; }

  bool is(CaloType t) const { return caloType() == t; }

  bool is(CaloID n) const { return caloID() == n; }

  bool is(Layout l) const { return layout() == l; }

  /** automatic conversion to int */
  operator int() const { return m_type; }

  /** explicit conversion to int */
  int toInt() const { return m_type; }

protected:
  int m_type;

  static const int fCaloType = 1;
  static const int fCaloID   = 10;
  static const int fLayout   = 1000;
  static const int fLayer    = 10000;
};

/** detailed string for calo type */
std::ostream& operator<<(std::ostream& os, const CHT& cht);

/** Return Layout based on the collection name, e.g. if name contains tolower("endcap") CHT::endcap is returned. In case no known layout
    is found, CHT::any is returned.*/
CHT::Layout layoutFromString(const std::string& name);

/** Return caloID based on the collection name, e.g. if name contains tolower("HCal") CHT::hcal is returned. In case no known type
    is found, CHT::unknown is returned.*/
CHT::CaloID caloIDFromString(const std::string& name);

/** Return caloType from string, e.g. if name contains tolower("Had") CHT::had is returned. In case no known type
    is found, CHT::em is returned.*/
CHT::CaloType caloTypeFromString(const std::string& name);

#endif