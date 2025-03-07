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

#include "Gaudi/Property.h"

// edm4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/Constants.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/TimeSeriesCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

#include "k4FWCore/Transformer.h"

#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
 
 struct CaloDigitizerFunc final
     : k4FWCore::Transformer<edm4hep::TimeSeriesCollection(const edm4hep::SimCalorimeterHitCollection&)> {
        CaloDigitizerFunc(const std::string& name, ISvcLocator* svcLoc)
       : Transformer(name, svcLoc, 
        {KeyValues("InputCollection", {"SimCaloHitsCollection"})},
        {KeyValues("OutputCollection", {"DigitsFloat"})}) {}
 
    StatusCode initialize() override{
        // Check if file exists
        if (m_SignalFileName.empty()) {
            error() << "Name of the file with the pulse shapes not provided!" << endmsg;
            return StatusCode::FAILURE;
        }
        if (gSystem->AccessPathName(m_SignalFileName.value().c_str())) {
            error() << "Provided file with the pulse shape values not found!" << endmsg;
            error() << "File path: " << m_SignalFileName.value() << endmsg;
            return StatusCode::FAILURE;
        }
        std::unique_ptr<TFile> pulseShapesFile(TFile::Open(m_SignalFileName.value().c_str(), "READ"));
        if (pulseShapesFile->IsZombie()) {
            error() << "Unable to open the file with the pulse shape values!" << endmsg;
            error() << "File path: " << m_SignalFileName.value() << endmsg;
            return StatusCode::FAILURE;
        } else {
            info() << "Using the following file with pulse shape values: "
                << m_SignalFileName.value() << endmsg;
        }

        // Read TTree and save the pulse shapes into map
        TTree* tree = nullptr;
        pulseShapesFile->GetObject(m_treename.value().c_str(), tree);
        ULong64_t readCellId;
        std::vector<float>* readPulse = nullptr;

        tree->Print("all");
        tree->Show(0);

        tree->SetBranchStatus("cellId", 1); tree->SetBranchAddress("cellId", &readCellId);
        tree->SetBranchStatus("PulseAmplitude", 1); tree->SetBranchAddress("PulseAmplitude", &readPulse);

        for (uint i = 0; i < tree->GetEntries(); i++) {
            tree->GetEntry(i);
            debug() << "cell ID: " << readCellId << endmsg;
            debug() << "PulseAmplitude size: " << readPulse->size() << endmsg;
            m_map.insert(std::pair<uint64_t, std::vector<float>>(readCellId, *readPulse));
        }

        info() << "Read " << m_map.size() << " cellIDs with pulse shapes from the file." << endmsg;
        info() << "Nentries in tree: " << tree->GetEntries() << endmsg;
        if (m_map.size() != static_cast<std::unordered_map<uint64_t, std::vector<float>>::size_type>(tree->GetEntries())) {
            error() << "Number of entries in the tree does not match the number of cellIDs read!" << endmsg;
            return StatusCode::FAILURE;
        }


        delete tree;
        pulseShapesFile->Close();
        return StatusCode::SUCCESS;
    }
   
    // This is the function that will be called to transform the data
   // Note that the function has to be const, as well as all pointers to collections
   // we get from the input
   edm4hep::TimeSeriesCollection operator()(const edm4hep::SimCalorimeterHitCollection& CaloHits) const override {
     info() << "calorimeter hits collection size: " << CaloHits.size() << endmsg;

     edm4hep::TimeSeriesCollection DigitsCollection;

     // Initialize a new map with the same keys as m_map but with values set to zero
    std::unordered_map<uint64_t, float> energyMap;
    for (const auto& entry : m_map) {
        energyMap[entry.first] = 0.0f;
    }

    if (m_map.size() != energyMap.size()) {
        error() << "Size of the map with pulse shapes does not match the size of the map with cellIDs!" << endmsg;
    }

    // Loop over CaloHits to accumulate energy
    for (const auto& hit : CaloHits) {
        uint64_t cellID = hit.getCellID();
        float energy = hit.getEnergy();
        
        if (energyMap.find(cellID) != energyMap.end()) {
            energyMap[cellID] += energy;
        } else {
            // If cellID is not in m_map, error!!
            error() << "CellID " << cellID << " not found in the map with pulse shapes!" << endmsg;
        }
    }


    // Loop over all hits and sum up the energy in the cells -> CaloHitContribution has timing info + energy from particle interactions + position of interaction -> but no cellID info.... DISCUSS w/ DENIS!!!!! No idea what the best way to include timing info is.... I would think that the timing info should be included in the digitized pulse (i.e. if within 25ns, it is in the first bin of the pulse but if within 50ns, it is in the second bin of the pulse, etc.). Not sure what the best way to code this would be.... maybe we calculate pulses on the fly? Or perhaps we need to create a vector of energy hits <E_t0, E_t1, E_t2, ...> where E_t0 is the summed energy within the t0 bin (first XX ns of sampling fraction), E_t1 is the summed energy within the t1 bin (second XX ns of sampling fraction), etc. Then, it is a matter of taking a dot product with a vector representing the pulse. However, one would have to take care to ensure that the pulse we create has the correct timing interval (or would it even matter since it is arbritary?)........

    // Loop over m_map to create DigitsCollection
    for (const auto& entry : energyMap) {
        uint64_t cellID = entry.first;
        float totalEnergy = entry.second;

        auto Digit = DigitsCollection.create();
        Digit.setCellID(cellID);
        Digit.setTime(0.0); // Placeholder for time info
        Digit.setInterval(m_samplingInterval); // Set the interval for the digitized pulse in ns

        const std::vector<float>& pulseAmplitudes = m_map.at(cellID);
        debug() << "cell ID: " << cellID << ", PulseAmplitude size: " << pulseAmplitudes.size() << endmsg;

        for (unsigned int i = 0; i < pulseAmplitudes.size(); i++) {
            Digit.addToAmplitude(totalEnergy * pulseAmplitudes[i]);
        }

       
    }

     return DigitsCollection;
   }

   StatusCode finalize() override{return StatusCode::SUCCESS;}


 
 private:
   /// Name of input root file that contains the TTree with cellID->vec<PulseAmplitude>
  Gaudi::Property<std::string> m_SignalFileName{this, "signalFileName", "PulseShape_map.root"};
  /// Name of the TTree in the input root file
  Gaudi::Property<std::string> m_treename{this, "treename", "Signal_shape"};
  /// Map to be used for the lookup of the pulse shapes
  Gaudi::Property<float> m_samplingInterval {this, "samplingInterval", 25.0, "Sampling interval in [ns]" };
  
  std::unordered_map<uint64_t, std::vector<float>> m_map;

 };
 
 DECLARE_COMPONENT(CaloDigitizerFunc)