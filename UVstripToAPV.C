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


// int StripToAPVpos(int strip, int axis, int modID){

//     // int nStrips = 

//     //total number of strips in the module
//     int nStrips = gemInfoMap[modID].nstripsuv[axis];

// }


std::map<int, std::pair<std::map<int, apvInfoVals>, std::map<int, apvInfoVals>>> UVstripToAPV(){

    std::map<int, std::pair<std::set<int>,std::set<int>>> modStrips;// modID=>{Ustrips, Vstrips}

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

    for (auto &[mod, uvStrips] : modStrips) {
        // for (mod)

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

            if (apvInfoMap.find(currStripKeys) != apvInfoMap.end()) {
                uStripInfo[stripID] = apvInfoMap[currStripKeys];  // Store unique strip -> APV mapping
            }
        }
       
        for (int stripID : vStrips){
            
            apvInfoKeys currStripKeys;
            currStripKeys.gemid = mod;
            currStripKeys.axis = 1;
            // apvInfoVals currStripInfo;

            // currStripKeys.pos = StripToAPVpos(int strip, int axis);
            currStripKeys.pos = stripID/128;//APV position
            // posInAPV = stripID%128;

            if (apvInfoMap.find(currStripKeys) != apvInfoMap.end()) {
                uStripInfo[stripID] = apvInfoMap[currStripKeys];  // Store unique strip -> APV mapping
            }
        }

        modUVstripAPVs[mod] = {uStripInfo, vStripInfo};

    }
    return modUVstripAPVs;
}
