#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>

// Function to read a reference file and return a map of module strips
std::map<int, std::pair<std::set<int>, std::set<int>>> makeStripsToMap(const char* refFile) {
    std::ifstream aRefFile(refFile);
    if (!aRefFile.is_open()) {
        std::cerr << "Error: Could not open file " << refFile << std::endl;
        return {}; // Return an empty map
    }

    std::map<int, std::pair<std::set<int>, std::set<int>>> modStrips; // Local storage for results
    std::string currLine;
    int currentM = -1; // Current module ID
    std::set<int> uStrips, vStrips;

    while (std::getline(aRefFile, currLine)) {
        std::stringstream lineStream(currLine);
        std::string lineStarter;
        lineStream >> lineStarter;

        if (lineStarter == "***") {
            if (currentM != -1) {
                // Save the previous module's strips
                modStrips[currentM] = {uStrips, vStrips};
                uStrips.clear();
                vStrips.clear();
            }
            std::string dumpStr;
            lineStream >> dumpStr //"Mod"
            >> dumpStr //"Num:"
            >> currentM; // Extract module ID

            // std::cout << "\ncurrent Module is " << currentM << std::endl;


        } else if (!currLine.empty() && currentM != -1) {

            
            std::string axis = lineStarter;
            int currStrip;
            std::string dumpStr;

            lineStream 
            // >> axis 
            >> dumpStr //Strip: 
            >> currStrip;

            if (axis == "U") {
                // std::cout << "Curr Strip: "<< currStrip << " on axis: " << axis << " of Mod " << currentM << std::endl; 
                uStrips.insert(currStrip);
            } else if (axis == "V") {
                // std::cout << "Curr Strip: "<< currStrip << " on axis: " << axis << "of Mod" << currentM<< std::endl;
                vStrips.insert(currStrip);
            }
        }
    }

    // Store the last module read
    if (currentM != -1) {
        modStrips[currentM] = {uStrips, vStrips};
    }

    std::cout << "Finished filling modStrips with " << modStrips.size() << " entries.\n";
    return modStrips;
}

// void outStripsToMap(){
    std::map<int, std::pair<std::set<int>, std::set<int>>>  outStripsToMap(){

        std::cout << "Running outStripsToMap() \n" << std::endl;

    const char* refFile = "outPhysicalStrips.txt";
    std::map<int, std::pair<std::set<int>, std::set<int>>> modStrips = makeStripsToMap(refFile);

    // // Output some results to verify correctness
    // for (const auto& [modID, stripPairs] : modStrips) {
    //     std::cout << "Module " << modID << " has " 
    //               << stripPairs.first.size() << " U-strips and "
    //               << stripPairs.second.size() << " V-strips.\n";
    // }
    return modStrips;
}