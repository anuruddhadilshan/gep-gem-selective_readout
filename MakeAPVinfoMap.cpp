#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include <iomanip>

// #include "DBread.h"

// Define key-value structs
struct apvInfoKeys {
    int gemid, axis, pos;
    
    // Define operator< to allow this struct to be used as a map key
    bool operator<(const apvInfoKeys& other) const {
        if (gemid != other.gemid) return gemid < other.gemid;
        if (axis != other.axis) return axis < other.axis;
        return pos < other.pos;
    }
};

struct apvInfoVals {
    int vtpcrate, fiber, adc_ch;
};

// Global map to store APV information
std::map<apvInfoKeys, apvInfoVals> apvInfoMap;

// Function to read reference file and fill `apvInfoMap`
// void MakeRefMap(const char* refFile) {
void MakeRefMap(const char* refFile) {
    std::ifstream aRefFile(refFile);
    if (!aRefFile.is_open()) {
        std::cerr << "Error: Could not open file " << refFile << std::endl;
        return;
    }

    std::string currLine;
    int currentM = -1;  // Track current module ID

    while (std::getline(aRefFile, currLine)) {
        std::stringstream lineStream(currLine);
        std::string lineStarter;
        lineStream >> lineStarter;

        if (lineStarter == "##") {
            if (std::getline(aRefFile, currLine)) {
                size_t mPos = currLine.find(".m");
                size_t chanPos = currLine.find(".chanmap");

                if (mPos != std::string::npos && chanPos != std::string::npos) {
                    size_t start = mPos + 2;
                    size_t end = currLine.find('.', start);
                    if (end != std::string::npos) {
                        currentM = std::stoi(currLine.substr(start, end - start));
                    }
                }
            }
        } 
        else if (!currLine.empty() && currentM != -1) {
            apvInfoKeys currAPVkeys;
            apvInfoVals currAPVvals;
            currAPVkeys.gemid = currentM;

            std::string dumpStr;
            std::stringstream dataStream(currLine);
            dataStream >> currAPVvals.vtpcrate >> dumpStr >> currAPVvals.fiber >> dumpStr
                       >> currAPVvals.adc_ch >> dumpStr >> currAPVkeys.pos >> dumpStr >> currAPVkeys.axis;

            apvInfoMap[currAPVkeys] = currAPVvals;
        } else {
            currentM = -1;  // Reset if empty line encountered
        }
    }

    std::cout << "Finished filling apvInfoMap with " << apvInfoMap.size() << " entries.\n";
}

// Function to print `apvInfoMap`
void printAPVinfoMap() {
    std::cout << "\nPrinting APV Info Map:\n";
    for (const auto& [key, val] : apvInfoMap) {
        std::cout << "GEMId: " << key.gemid << ", Axis: " << key.axis
                  << ", Pos: " << key.pos << std::endl;
        std::cout << "  VTPcrate: " << val.vtpcrate << ", Fiber: " << val.fiber
                  << ", ADC_ch: " << val.adc_ch << std::endl;
        std::cout << "------------------------------------\n";
    }
}

// Function to write `apvInfoMap` to a file
void OutputRefMap() {
    std::ofstream mapFile("gemFT_Map_TEST.txt");
    if (!mapFile.is_open()) {
        std::cerr << "Error: Could not create file gemFT_Map_TEST.txt\n";
        return;
    }

    mapFile << "#apvInfoMap actually maps key and value structs to eachother\n" ;
    mapFile << "## the apvInfoKeys consists of: module id, axis(U/V depending on map), and  pos(along axis)\n";
    mapFile << "## the apvInfoVals consists of: VTPcrate, Fiber(MPD) ID, and  the ADC channel\n\n";

    // Set column widths for formatting
    int colWidth = 10; // Adjust as needed for alignment

    mapFile << std::left << std::setw(colWidth) << "GEMId"
            << std::setw(colWidth) << "Axis"
            << std::setw(colWidth) << "Pos"
            << std::setw(colWidth) << "VTPcrate"
            << std::setw(colWidth) << "Fiber"
            << std::setw(colWidth) << "ADC_ch"
            << "\n";
    
    mapFile << std::string(6 * colWidth, '-') << "\n"; // Create a line separator

    for (const auto& [key, val] : apvInfoMap) {
        mapFile << std::left << std::setw(colWidth) << key.gemid
                << std::setw(colWidth) << key.axis
                << std::setw(colWidth) << key.pos
                << std::setw(colWidth) << val.vtpcrate
                << std::setw(colWidth) << val.fiber
                << std::setw(colWidth) << val.adc_ch
                << "\n";
    }

    std::cout << "Parsed data written to gemFT_Map_TEST.txt\n";
}

const char* refFile = "db_sbs.gemFT_TEST.txt";


std::map <apvInfoKeys, apvInfoVals> GetAPVinfoMap(){
    MakeRefMap(refFile);

    // OutputRefMap();
    // printAPVinfoMap();

    // return 0;
    return apvInfoMap;
}


// void MakeAPVinfoMap() {
int MakeAPVinfoMap() {
    // const char* refFile = "db_sbs.gemFT_TEST.txt";
    MakeRefMap(refFile);

    OutputRefMap();
    printAPVinfoMap();

    return 0;
    // return apvInfoMap;
}
