// This script takes the local GEM database and ROI text file as inputs, and outputs APV lists, per ECal bin. The exact output format is TBD.

#include "DBread.h"
#include "ROIread.h"
#include "GEMModROItoStrips.h"
#include "GetGemInfoMap.cxx"
#include "TVector2.h"
#include "TH2.h"
#include "TCanvas.h"
#include <map>
#include <set>
#include <array>

void estimateROIreduction(const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt") {

    std::map<int, gemInfo> refModMap = GetGemInfoMap();
    DBread db{ db_local };
    ROIread roi{ roi_file };

    if (db.returnFileReadStatus() == -1 || roi.returnFileReadStatus() == -1) {
        std::cerr << "Exiting the program!!!" << std::endl << std::endl;
        return;
    }

    std::map<int, GEMLayer> map_GEMLayers = db.returnGEMLayerMap();
    std::map<int, std::map<int, ROI>> map_ROIsByBinsAndLayers = roi.return_ROIMap();

    int nECalBins = map_ROIsByBinsAndLayers.size();
    int nGEMLayers = map_GEMLayers.size();

    std::map<int, GEMLayerROItoStrips> map_GEMLayerROItoStrips;
    for (const auto& [layerNum, gemLayerInstance] : map_GEMLayers) {
        map_GEMLayerROItoStrips.emplace(layerNum, GEMLayerROItoStrips{gemLayerInstance});
    }

    std::map<int, std::map<int, std::pair<std::set<int>, std::set<int>>>> map_physicalUVStrips_byECalBin_byGEMMod;
    std::map<int, std::set<apvInfoVals>> ECalBinAPVinfo;
    std::set<apvInfoKeys> missingKeys;
    std::map<int, std::map<int, std::array<double, 2>>> ratioReducPerEcalbin;
    TH2D* uvReducHist = new TH2D("uvReducHist", "U vs V Strip Reduction", 100, 0, 1, 100, 0, 1);

    // NEW: Tracking for averages
    std::map<int, int> activeModulesPerBin;
    std::map<int, int> totalUstripsPerBin;
    std::map<int, int> totalVstripsPerBin;

    for (auto& [binNum, map_ROIbyLayer_forThisBin] : map_ROIsByBinsAndLayers) {
        std::map<int, std::pair<std::set<int>, std::set<int>>> map_allUVstripsSetsForAllModules_forThisECalBin;
        std::set<apvInfoVals> currECalBinAPVs;

        int activeModCount = 0;
        int totalUstrips = 0;
        int totalVstrips = 0;

        for (auto& [layerNum, gemLayerROItoStripsInstance] : map_GEMLayerROItoStrips) {
            std::map<int, std::pair<std::set<int>, std::set<int>>> map_thisLayer_allROIsForAllModules_forThisECalBin =
                gemLayerROItoStripsInstance.takeROI_givePhysicalUVStrips(map_ROIsByBinsAndLayers.at(binNum).at(layerNum));

            for (int modNum = 0; modNum < 14; ++modNum) {
                if (map_thisLayer_allROIsForAllModules_forThisECalBin.find(modNum) == map_thisLayer_allROIsForAllModules_forThisECalBin.end()) {
                    std::cout << "\nECalBin " << binNum << "\n!No strips on in mod: " << modNum << std::endl;
                    continue;
                }

                const auto& modUVstripPair = map_thisLayer_allROIsForAllModules_forThisECalBin.at(modNum);
                const std::set<int>& uStrips = modUVstripPair.first;
                const std::set<int>& vStrips = modUVstripPair.second;

                double uReduced = static_cast<double>(uStrips.size()) / refModMap[modNum].nstripsuv[0];
                double vReduced = static_cast<double>(vStrips.size()) / refModMap[modNum].nstripsuv[1];

                ratioReducPerEcalbin[binNum][modNum][0] = uReduced;
                ratioReducPerEcalbin[binNum][modNum][1] = vReduced;
                uvReducHist->Fill(uReduced, vReduced);

                if (!uStrips.empty() || !vStrips.empty()) {
                    ++activeModCount;
                    totalUstrips += uStrips.size();
                    totalVstrips += vStrips.size();
                }

                std::cout << "\nECalBin " << binNum << "\nmod: " << modNum << std::endl;
                std::cout << "nROI uStrips " << uStrips.size() << " nROI vStrips " << vStrips.size() << std::endl;
                std::cout << "nMod uStrips " << refModMap[modNum].nstripsuv[0]
                          << " nMod vStrips " << refModMap[modNum].nstripsuv[1] << std::endl;
                std::cout << "reduced U by " << uReduced << " reduced V by " << vReduced << std::endl;
            }
        }

        activeModulesPerBin[binNum] = activeModCount;
        totalUstripsPerBin[binNum] = totalUstrips;
        totalVstripsPerBin[binNum] = totalVstrips;
        map_physicalUVStrips_byECalBin_byGEMMod[binNum] = map_allUVstripsSetsForAllModules_forThisECalBin;
        ECalBinAPVinfo[binNum] = currECalBinAPVs;
    }

    TCanvas* ROIreducCanv = new TCanvas();
    uvReducHist->Draw();

    // Final summary output
    double sumActiveMods = 0;
    double sumUstrips = 0;
    double sumVstrips = 0;
    int numBins = activeModulesPerBin.size();

    for (const auto& [binNum, activeModCount] : activeModulesPerBin) {
        sumActiveMods += activeModCount;
        sumUstrips += totalUstripsPerBin[binNum];
        sumVstrips += totalVstripsPerBin[binNum];
    }

    double avgModsPerBin = sumActiveMods / numBins;
    double avgUstripsPerBin = sumUstrips / numBins;
    double avgVstripsPerBin = sumVstrips / numBins;

    // int totalModules = 14;
    // double totalPossibleStripsU = numBins * totalModules * refModMap[0].nstripsuv[0];
    // double totalPossibleStripsV = numBins * totalModules * refModMap[0].nstripsuv[1];
    double totalPossibleStripsU = 0;
double totalPossibleStripsV = 0;

for (int modNum = 0; modNum < 14; ++modNum) {
    totalPossibleStripsU += numBins * refModMap[modNum].nstripsuv[0];
    totalPossibleStripsV += numBins * refModMap[modNum].nstripsuv[1];
}

int totalStripsU_perBin = 0;
int totalStripsV_perBin = 0;

for (int modNum = 0; modNum < 14; ++modNum) {
    totalStripsU_perBin += refModMap[modNum].nstripsuv[0];
    totalStripsV_perBin += refModMap[modNum].nstripsuv[1];
}


    // std::cout << "\nTotal Ustrips per bin: " << totalPossibleStripsU<< std::endl;
    // std::cout << "\nTotal Vstrips per bin: " << totalPossibleStripsV<< std::endl;

    double avgFractionUsedU = sumUstrips / totalPossibleStripsU;
double avgFractionUsedV = sumVstrips / totalPossibleStripsV;

double avgStripsTurnedOffU = (1.0 - avgFractionUsedU) * (totalPossibleStripsU / numBins / 14);
    double avgStripsTurnedOffV = (1.0 - avgFractionUsedV) * (totalPossibleStripsV / numBins / 14);

    std::cout << "\n----- ROI Reduction Summary -----\n";
    std::cout << "Total U strips possible per bin: " << totalStripsU_perBin << std::endl;
    std::cout << "Total V strips possible per bin: " << totalStripsV_perBin << std::endl;
    std::cout << "Average active modules per ECal bin: " << avgModsPerBin << std::endl;
    std::cout << "Average U strips per ECal bin (totaled over all modules): " << avgUstripsPerBin << std::endl;
    std::cout << "Average V strips per ECal bin (totaled over all modules): " << avgVstripsPerBin << std::endl;
    std::cout << "Average fraction of U strips used: " << avgFractionUsedU << std::endl;
    std::cout << "Average fraction of V strips used: " << avgFractionUsedV << std::endl;
    std::cout << "Average U strips turned OFF per module per bin: " << avgStripsTurnedOffU << std::endl;
    std::cout << "Average V strips turned OFF per module per bin: " << avgStripsTurnedOffV << std::endl;
    std::cout << "Average U strips turned OFF for all modules per bin: " << avgStripsTurnedOffU*14 << std::endl;
    std::cout << "Average V strips turned OFF for all modules per bin: " << avgStripsTurnedOffV*14 << std::endl;

}
