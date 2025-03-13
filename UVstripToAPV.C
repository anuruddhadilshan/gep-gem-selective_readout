#include "DBread.h"

// #include "TVector2.h"
// #include "TMath.h"
#include <iostream>
#include <map>
#include <array>
#include <vector>
#include <set>
// #include <pair.h>


#include "GetchanMap.cxx"
#include "GetGemInfoMap.cxx"

#include "outStripsToMap.cpp"


// int StripToAPVpos(int strip, int axis, int modID){

//     // int nStrips = 

//     //total number of strips in the module
//     int nStrips = gemInfoMap[modID].nstripsuv[axis];

// }


// std::map<int, std::pair<std::map<int, apvInfoVals>, std::map<int, apvInfoVals>>> UVstripToAPV( std::map<int, std::pair<std::set<int>,std::set<int>>> modStrips){

std::map<int, std::pair<std::map<int, apvInfoVals>, std::map<int, apvInfoVals>>> makeUVstripToAPVmap( std::map<int, std::pair<std::set<int>,std::set<int>>> modStrips){

    std::cout << "Running makeUVstripToAPVmap()" << std::endl;

    // std::map<int, std::pair<std::set<int>,std::set<int>>> modStrips;// modID=>{Ustrips, Vstrips}

    // std::set<int>modIDs=modStripCounts~map;
    //  <modID < Ustrip, info> <
    std::map< int, //modID
        std::pair< 
                std::map<int, apvInfoVals>, //U strip number, APV info(vtpCrate, fiber, adc_ch)
                std::map<int, apvInfoVals> // Vstrip number, (vtpCrate, fiber, adc_ch)
    > > modUVstripAPVs;

    std::map<apvInfoKeys, apvInfoVals> apvInfoMap = GetchanMap("db_sbs.gemFT_TEST.txt");

    // std::map<int, gemInfo> gemInfoMap = GetGemInfoMap();

    // GetGemInfoMap();
    // printGEMinfoMap();

    std::cout << "size of modStrips " << modStrips.size() << std::endl;

    for (auto &[mod, uvStrips] : modStrips) {
        // for (mod)

        std::cout << "\nstarting for (auto &[mod, uvStrips] : modStrips) {" << std::endl;

        std::cout << "uStrips size " << uvStrips.first.size();
        std::cout << "vStrips size " << uvStrips.second.size();

        std::set<int>uStrips = uvStrips.first;
        std::set<int>vStrips = uvStrips.second;

        // int nStripsU = uStrips.size();
        // int nStripsV = vStrips.size();
        
        //strip Num => APV
        std::map<int, apvInfoVals> uStripInfo;
        std::map<int, apvInfoVals> vStripInfo;

        for (int stripID : uStrips){
            
            apvInfoKeys currStripKeys;
            currStripKeys.gemid = mod;
            currStripKeys.axis = 0;
            // apvInfoVals currStripInfo;

            // currStripKeys.pos = StripToAPVpos(int strip, int axis);
            currStripKeys.pos = stripID/128;//APV position
            // posInAPV = stripID%128;

            std::cout << "Checking modID: " << mod << " Strip: " << stripID 
          << " Axis: " << currStripKeys.axis << " Pos: " << currStripKeys.pos << std::endl;
if (apvInfoMap.find(currStripKeys) == apvInfoMap.end()) {
    std::cout << "  --> APV key NOT FOUND!" << std::endl;
}


            if (apvInfoMap.find(currStripKeys) != apvInfoMap.end()) {
                uStripInfo[stripID] = apvInfoMap[currStripKeys];  // Store unique strip -> APV mapping
            }
            std::cout<< "modID: " << mod << " UStrip: " << stripID << " | "
            << "VTP Crate: " << uStripInfo[stripID].vtpcrate
            << " Fiber: " << uStripInfo[stripID].fiber
            << " ADC Ch: " <<uStripInfo[stripID].adc_ch
            << " Invert: " << uStripInfo[stripID].invert << std::endl;
        }
       
        for (int stripID : vStrips){
            
            apvInfoKeys currStripKeys;
            currStripKeys.gemid = mod;
            currStripKeys.axis = 1;
            // apvInfoVals currStripInfo;

            // currStripKeys.pos = StripToAPVpos(int strip, int axis);
            currStripKeys.pos = stripID/128;//APV position
            // posInAPV = stripID%128;


            std::cout << "Checking modID: " << mod << " Strip: " << stripID 
          << " Axis: " << currStripKeys.axis << " Pos: " << currStripKeys.pos << std::endl;
if (apvInfoMap.find(currStripKeys) == apvInfoMap.end()) {
    std::cout << "  --> APV key NOT FOUND!" << std::endl;
}


            if (apvInfoMap.find(currStripKeys) != apvInfoMap.end()) {
                vStripInfo[stripID] = apvInfoMap[currStripKeys];  // Store unique strip -> APV mapping

                std::cout<< "modID: " << mod << " VStrip: " << stripID << " | "
                << "VTP Crate: " << vStripInfo[stripID].vtpcrate
                << " Fiber: " << vStripInfo[stripID].fiber
                << " ADC Ch: " <<vStripInfo[stripID].adc_ch
                << " Invert: " << vStripInfo[stripID].invert << std::endl;
            }
        }

        modUVstripAPVs[mod] = {uStripInfo, vStripInfo};

    }
    return modUVstripAPVs;
}


void OutputstripToAPVmap(const std::map<int, std::pair<std::map<int, apvInfoVals>, std::map<int, apvInfoVals>>>& modUVstripAPVs) {
    std::ofstream mapFile("StripAPVmap.txt");
    if (!mapFile.is_open()) {
        std::cerr << "Error: Could not create file StripAPVmap.txt\n";
        return;
    }

    // Print a descriptive header
    mapFile << "# Mapping of GEM strip numbers to APV electronics\n";
    mapFile << "# Format: Module | Axis (0=U, 1=V) | Strip | VTP Crate | Fiber (MPD ID) | ADC Channel | Invert\n\n";

    // Set column widths for better formatting
    int colWidth = 12;

    mapFile << std::left 
            << std::setw(colWidth) << "Module"
            << std::setw(colWidth) << "Axis"
            << std::setw(colWidth) << "Strip"
            << std::setw(colWidth) << "VTP Crate"
            << std::setw(colWidth) << "Fiber"
            << std::setw(colWidth) << "ADC Ch"
            << std::setw(colWidth) << "Invert"
            << "\n";

    mapFile << std::string(7 * colWidth, '-') << "\n"; // Line separator

    // Iterate over the module map
    for (const auto& [modID, stripPairs] : modUVstripAPVs) {
        const auto& uStrips = stripPairs.first;
        const auto& vStrips = stripPairs.second;

        // Output U-strips
        for (const auto& [strip, apvInfo] : uStrips) {
            mapFile << std::left 
                    << std::setw(colWidth) << modID
                    << std::setw(colWidth) << "0"  // Axis 0 = U
                    << std::setw(colWidth) << strip
                    << std::setw(colWidth) << apvInfo.vtpcrate
                    << std::setw(colWidth) << apvInfo.fiber
                    << std::setw(colWidth) << apvInfo.adc_ch
                    << std::setw(colWidth) << apvInfo.invert
                    << "\n";
        }

        // Output V-strips
        for (const auto& [strip, apvInfo] : vStrips) {
            mapFile << std::left 
                    << std::setw(colWidth) << modID
                    << std::setw(colWidth) << "1"  // Axis 1 = V
                    << std::setw(colWidth) << strip
                    << std::setw(colWidth) << apvInfo.vtpcrate
                    << std::setw(colWidth) << apvInfo.fiber
                    << std::setw(colWidth) << apvInfo.adc_ch
                    << std::setw(colWidth) << apvInfo.invert
                    << "\n";
        }
    }

    std::cout << "Parsed data written to StripAPVmap.txt\n";
}


void testUVstripToAPV() {
    std::cout << "testUVstripToAPV() is executing..." << std::endl;

    auto modStrips = outStripsToMap(); // This gives correct format
    auto modUVstripAPVs = makeUVstripToAPVmap(modStrips); // Pass correct format
    
    std::cout <<"\n\n Finished testUVstripToAPV" << std::endl;
}

void UVstripToAPV() {
    std::cout << "UVstripToAPV() is executing..." << std::endl;

    auto modStrips = outStripsToMap(); // Ensure it's retrieved before use
    std::cout << "Entering makeUVstripToAPVmap() with " << modStrips.size() << " modules." << std::endl;

    testUVstripToAPV();
    auto modUVstripAPVs = makeUVstripToAPVmap(modStrips); // Pass correct format

    OutputstripToAPVmap(modUVstripAPVs);

    
}