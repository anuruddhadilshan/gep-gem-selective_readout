#include "DBread.h"

<<<<<<< HEAD
#include "TVector2.h"
#include "TMath.h"
=======
// #include "TVector2.h"
// #include "TMath.h"
>>>>>>> origin/main
#include <iostream>
#include <map>
#include <array>
#include <vector>
#include <set>
<<<<<<< HEAD


#include "GetchanMap.cxx"


=======
// #include <pair.h>


#include "GetchanMap.cxx"
#include "GetGemInfoMap.cxx"


// int StripToAPVpos(int strip, int axis, int modID){

//     // int nStrips = 

//     //total number of strips in the module
//     int nStrips = gemInfoMap[modID].nstripsuv[axis];

// }


std::map<int, std::pair<std::pair<std::set<int>, std::set<apvInfoVals>>, std::pair<std::set<int>, std::set<apvInfoVals>>>> UVstripToAPV(){

    std::map< int, std::pair< std::set<int>, std::set<int> > > modStrips;

    // std::set<int>modIDs=modStripCounts~map;
    //  <modID < Ustrip, info> <
    std::map< int, //modID
        std::pair< 
                std::pair<std::set<int>, std::set<apvInfoVals>>, //U strip number, APV info(vtpCrate, fiber, adc_ch)
                std::pair <std::set<int>, std::set<apvInfoVals>> // Vstrip number, (vtpCrate, fiber, adc_ch)
    > > modUVstripAPVs;

    std::map<apvInfoKeys, apvInfoVals> apvInfoMap = GetchanMap("db_sbs.gemFT_TEST.txt");

    // std::map<int, gemInfo> gemInfoMap = GetGemInfoMap();

    // GetGemInfoMap();
    // printGEMinfoMap();

    for (auto [mod, uvStrips] : modStrips) {
        // for (mod)

        std::set<int>uStrips = uvStrips.first;
        std::set<int>vStrips = uvStrips.second;

        int nStripsU = uStrips.size();
        int nStripsV = vStrips.size();
        
        std::pair<std::set<int>, std::set<apvInfoVals>> uStripInfo;//<stripNum, info
        std::pair<std::set<int>, std::set<apvInfoVals>> vStripInfo;//<stripNum, info


        for (int stripID = 0; stripID < uStrips.size(); stripID++){

            apvInfoKeys currStripKeys;
            currStripKeys.gemid = mod;
            currStripKeys.axis = 0;
            // apvInfoVals currStripInfo;

            // currStripKeys.pos = StripToAPVpos(int strip, int axis);
            currStripKeys.pos = static_cast<int>(stripID/128);

            uStripInfo.first.emplace(stripID);
            uStripInfo.second.emplace(apvInfoMap[currStripKeys]);            
        }
        
        for (int stripID = 0; stripID < vStrips.size(); stripID++){ 
            apvInfoKeys currStripKeys;
            currStripKeys.gemid = mod;
            currStripKeys.axis = 1;
            // apvInfoVals currStripInfo;

            // currStripKeys.pos = StripToAPVpos(int strip, int axis);
            currStripKeys.pos = static_cast<int>(stripID/128);
            
            vStripInfo.first.emplace(stripID);
            vStripInfo.second.emplace(apvInfoMap[currStripKeys]);

        }
        modUVstripAPVs[mod].first.first.emplace(uStripInfo.first);
        modUVstripAPVs[mod].first.second.emplace(uStripInfo.first);
        modUVstripAPVs[mod].second.first.emplace(uStripInfo.second);
        modUVstripAPVs[mod].second.second.emplace(uStripInfo.second);

    }

    return modUVstripAPVs;


}
>>>>>>> origin/main
