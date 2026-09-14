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

/* NOTE ON DUPLICATION
 *
 * This file (and its .cc) originate in MarlinUtil and are also carried, in near-identical
 * form, by k4GaudiPandora. That package compiles them privately into its Gaudi plugin module
 * rather than exposing an installed target, so there is nothing for this package to link
 * against and the sources are duplicated here instead.
 *
 * This copy has since been made a little safer to use: the three enums are scoped
 * (enum class), and CaloType gained the "unknown" fallback that CaloID and Layout already
 * had, so a calorimeter type that matches nothing is no longer silently reported as
 * electromagnetic. The encoded integer is unchanged for all pre-existing values.
 *
 * The duplication should be resolved by promoting these to a single shared, installed
 * location that both packages depend on.
 */

#include <ostream>
#include <string>

/** Helper class for decoding/encoding lcio::CalorimeterHit types for the ILD
 *  detector. The encoding is: caloType + 10 * caloID + 1000 * layout + 10000 * layerNum <br>
 *  (see enums: CaloType, CaloID and Layout for possible values).<br>
 *  Example usage: <br>
 *  <pre>
 *     lcio::CalorimeterHit* cHit = .... ;
 *
 *     // set the type (e.g. in digitization )
 *     cHit->setType( CHT( CHT::CaloType::em, CHT::CaloID::ecal, CHT::Layout::plug, 12 ) ) ;
 *
 *     ...
 *
 *     CHT cht = cHit->getType() ;
 *
 *     //   sum energies for electromagentic, hadronic and tailcatcher:
 *     if( cht.is( CHT::CaloType::em ) )
 *          e_em +=  cHit->getEnergy() ;
 *     else
 *       if ( cht.is(CHT::CaloType::had ) )
 *          e_had += cHit->getEnergy() ;
 *       else
 *          e_muon += cHit->getEnergy() ;
 *
 *     // use only EcalPlug hits:
 *     if( cht.is( CHT::CaloID::ecal) && cht.is( CHT::Layout::plug) )
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
  /* The enumerator values below are a wire format, not free choices: they are the positional
   * digits of the integer stored in CalorimeterHit::type (see the encoding above), and they
   * must keep matching k4GaudiPandora's copy of this file, which reads back hits written here.
   * Renumbering any of them changes the type of every hit and silently desynchronises the two
   * packages.
   *
   * That is why "no match" sits in a different place in each enum. CaloID and Layout were
   * defined by ILD with theirs at 0 ("unknown" and "any" respectively). CaloType could not
   * follow, because 0 was already taken by "em", so its "unknown" is appended at the end
   * instead. The inconsistency is inherited, not deliberate.
   */

  /** calorimeter types */
  enum class CaloType { em = 0, had = 1, muon = 2, unknown = 3 };

  /** calo ids - specific to ILD */
  enum class CaloID { unknown = 0, ecal = 1, hcal = 2, yoke = 3, lcal = 4, lhcal = 5, bcal = 6 };

  /** calo layout / subdetector; "any" is this enum's equivalent of "unknown" */
  enum class Layout { any = 0, barrel = 1, endcap = 2, plug = 3, ring = 4 };

  /** C'tor for initialization from CalorimeterHit::getType()  */
  CHT(int type) : m_type(type) {}

  /** C'tor  for encoding the calo type inforamtion  */
  CHT(CaloType c, CaloID n, Layout l, unsigned lay)
      : m_type(static_cast<int>(c) * fCaloType + static_cast<int>(n) * fCaloID + static_cast<int>(l) * fLayout +
               lay * fLayer) {}

  /** calorimeter type: CHT::em , CHT::had, CHT::muon */
  CaloType caloType() const { return static_cast<CaloType>(m_type % fCaloID); }

  /** calo ID - see enum CaloID for allowed values */
  CaloID caloID() const { return static_cast<CaloID>((m_type % fLayout) / fCaloID); }

  /** calo layout - see enum layout for allowed values */
  Layout layout() const { return static_cast<Layout>((m_type % fLayer) / fLayout); }

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
  static const int fCaloID = 10;
  static const int fLayout = 1000;
  static const int fLayer = 10000;
};

/** detailed string for calo type */
std::ostream& operator<<(std::ostream& os, const CHT& cht);

/** Return Layout based on the collection name, e.g. if name contains tolower("endcap") CHT::Layout::endcap is
   returned. In case no known layout is found, CHT::Layout::any is returned.*/
CHT::Layout layoutFromString(const std::string& name);

/** Return caloID based on the collection name, e.g. if name contains tolower("HCal") CHT::CaloID::hcal is returned.
   In case no known type is found, CHT::CaloID::unknown is returned.*/
CHT::CaloID caloIDFromString(const std::string& name);

/** Return caloType from string, e.g. if name contains tolower("Had") CHT::CaloType::had is returned. In case no
    known type is found, CHT::CaloType::unknown is returned.*/
CHT::CaloType caloTypeFromString(const std::string& name);

#endif
