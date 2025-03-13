#include <iostream>

#include "FindAPV.h"
// #include "SBSGEMModule.h"
// #include "THaSubDetector.h"
// #include "ECalToHCal.h"

// #include "FindAPV.h"

// #include "MakeAPVinfoMap.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include <numbers>

// #include "TObject.h"
// #include "TString.h"


// #include "TVector2.h"
// #include "TVector3.h"

// #include "TRotation.h"
// #include "TH1D.h"
// #include "TH2D.h"
// #include "TF1.h"
// #include "TGraphErrors.h"
// #include "TClonesArray.h"
#include <algorithm>
#include <iomanip>
#include <math.h>


struct apvInfo {
  // could be arranged as multiDim arrays
  //vtpcrate[fiber][adc_ch]=gemID[axis][pos]
  //or TTrees?
  
  //use these 3 as the "keys"
  int gemid; int axis; int pos;
  
  //fill these values they id the APV
  int vtpcrate;
  int fiber; //NOTE: == mpd_id
  int adc_ch;
  int invert;
  
  int strip; //if available
  
  // Define operator< for sorting in std::set
  bool operator<(const apvInfo& other) const;
  
  void print() const {
    std::cout << "GEMId: " << gemid << std::endl;
    std::cout << "Axis: " << axis << std::endl;
    std::cout << "Pos: " << pos << std::endl;
    std::cout << "invert: " << invert << std::endl;
    std::cout << "VTPcrate: " << vtpcrate << std::endl;
    std::cout << "fiber: " << fiber << std::endl;
    std::cout << "adc_ch: " << adc_ch << std::endl;
  }
  
};

//module id, axis, position on axis
struct apvInfoKeys {
  int gemid, axis, pos;

  void print() const {
    std::cout 
    << "GEMId: " << gemid 
    << " Axis: " << axis 
    << " Pos: " << pos 
    << std::endl;
  };
  bool operator<(const apvInfoKeys& other) const {
    if (gemid != other.gemid) return gemid < other.gemid;
    if (axis != other.axis) return axis < other.axis;
    return pos < other.pos;
  }
};

//vtpcrate, fiber, adc_ch
struct apvInfoVals{
  //fill these values they id the APV
  int vtpcrate;
  int fiber; //NOTE: == mpd_id
  int adc_ch;
  int invert;
  void print() const {
    std::cout << "VTPcrate: " << vtpcrate 
    << " Fiber: " << fiber 
    << " ADC_ch: " << adc_ch 
    << " invert? " << invert 
    << std::endl; 
  }
};

std::map<apvInfoKeys, apvInfoVals> apvInfoMap; 


//changes invert to either 1 or -1 to be used as a factor
int SetInvert(int invert){
  if (invert == 0){return 1;}
  if (invert == 1){return -1;}
}



//Takes in txt files
void MakechanMap(const char* refFile) {
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
                     >> currAPVvals.adc_ch >> dumpStr >> currAPVkeys.pos 
                     >> currAPVvals.invert
                     >> currAPVkeys.axis;

          
          
          
                    //  currAPVvals.invert = SetInvert(currAPVvals.invert);



                     
          apvInfoMap[currAPVkeys] = currAPVvals;
      } else {
          currentM = -1;  // Reset if empty line encountered
      }
  }

  std::cout << "Finished filling apvInfoMap with " << apvInfoMap.size() << " entries.\n";
}

  // Function to print `apvInfoMap`
void printAPVinfoMap() {
  std::cout << "Printing APV Info Map:\n";
  for (const auto& [key, val] : apvInfoMap) {
      std::cout << "GEMId: " << key.gemid << ", Axis: " << key.axis
                << ", Pos: " << key.pos << std::endl;
      std::cout << "VTPcrate: " << val.vtpcrate << ", Fiber: " << val.fiber
                << ", invert: " << val.invert
                << ", ADC_ch: " << val.adc_ch << std::endl;
      std::cout << "------------------------------------\n";
  }
};


// Function to write `apvInfoMap` to a file
void OutputAPVinfoMap() {
  std::ofstream mapFile("APV_Map_TEST.txt");
  if (!mapFile.is_open()) {
      std::cerr << "Error: Could not create file APV_Map_TEST.txt\n";
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

  std::cout << "Parsed data written to APV_Map_TEST.txt\n";
}

void TestAPVInfoMap(const char* refFile){
  // APVFinder *anAPVFinder = new APVFinder();
  MakechanMap(refFile);
  printAPVinfoMap();
  OutputAPVinfoMap();
}

std::map<apvInfoKeys, apvInfoVals> GetchanMap(const char* refFile){
  // MakechanMap("db_sbs.gemFT_TEST.txt");
  MakechanMap(refFile);
  
  return apvInfoMap;
}