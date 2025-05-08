#include <iostream>

// #include "FindAPV.h"
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


      //Set default values for decode map parameters:
  int fN_APV25_CHAN = 128;
  int fMPDMAP_ROW_SIZE = 9;

  //arrays to hold raw data from one APV card:
  // fStripAPV.resize( MAXNSAMP_PER_APV );
  // fRawStripAPV.resize( MAXNSAMP_PER_APV );
  // fRawADC_APV.resize( MAXNSAMP_PER_APV );

  //default to 
  //fMAX2DHITS = 250000;
  // fMAX2DHITS = 10000;

  // fAPVmapping = SBSGEM::kUVA_XY; //default to UVA X/Y style APV mapping, but require this in the database::
  int fAPVmapping = 2; 

  // fPxU = cos(UAngle);
  // fPyU = sin(UAngle);
  // fPxV = cos(VAngle);
  // fPyV = sin(VAngle);

  int nMods = 14; //number of modules

  double stripOffset = .004;

  double fUStripOffset=.004; double fVStripOffset=.004;
  std::map<int, std::array<double, 3>> ModPositionmap = {
    {0, {0., 0., 0.}},
    {1, {0., 0., 0.075}},
    {2, {0., 0., 0.202}},
    {3, {0., 0., 0.341}},
    {4, {0., 0., 0.484}},
    {5, {0., 0., 0.624}},

    {6, {-0.766, 0., 0.769}},
    {7, {-0.256, 0., 0.737}},
    {8, {0.265, 0., 0.769}},
    {9, {0.765, 0., 0.737}},

    {10, {-0.766, 0., 0.897}},
    {11, {-0.256, 0., 0.865}},
    {12, {0.256, 0., 0.897}},
    {13, {0.765, 0., 0.865}}
};
  


  struct gemInfo{
    // int modID; made into key to gemInfoMap
    int layer;
    int apvmap;
    // TVector3 position;
    std::array<double, 3> position;
    // double xPos, yPos, zPos;
    // TVector2 angle; xyz ang
    // TVector3 size;
    std::array<double, 3> size;
    // double xSize, ySize, zSize;


    std::array<double, 2> uvangles;
    // double uangle;  double vangle;

    std::array<double, 2> uvoffsets;
    // double uoffset;  double voffset;

    std::array<int, 2> nstripsuv;
    // int nstripsu;  int nstripsv;

    std::array<int, 2> NuvAPVs;
    // int NuAPVs; int NvAPVs;
  

  
    void print() const{
        // std::cout << "ModID: " << modID << std::endl;
        std::cout << "Layer: " << layer << std::endl;
        std::cout << "APVmap: " << apvmap << std::endl;
        std::cout << "Position: (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
        // std::cout << "Angle: (" << angle.X() << ", " << angle.Y() << ","
        // << angle.Z() << ")" << std::endl;
        std::cout << "Size: (" << size[0] << ", " << size[1] << ", " << size[2] << ")" << std::endl;
        std::cout << "Uangle: " << uvangles[0] << std::endl;
        std::cout << "Vangle: " << uvangles[1] << std::endl;
        std::cout << "Uoffset: " << uvoffsets[0] << std::endl;
        std::cout << "Voffset: " << uvoffsets[1] << std::endl;
        std::cout << "NstripsU: " << nstripsuv[0] << std::endl;
        std::cout << "NstripsV: " << nstripsuv[1] << std::endl;
        std::cout << "UAPVs: " << NuvAPVs[0] << std::endl;
        std::cout << "VAPVs: " << NuvAPVs[1] << std::endl;
      };
  
    };

std::map<int, gemInfo> gemInfoMap;

  

int GetLayerOfMod(int modID){
  if (modID >= 0 && modID < 6) { return modID; }
    else if (modID >= 6 && modID < 10) { return 6; }
    else if (modID >= 10 && modID < 14) { return 7; }
    return -1;
}

int GetMod_apvmap(int modID){
  if(modID>=0 && modID<6){ return 2;}
  else{return 1;}
};



  
  // std::array<double, 3> GetModDimensions(int apvmap){
  std::array<double, 3> GetModDimensions(int modID){
    if (modID<6){
      //in meters
      return std::array<double, 3> {1.5, .4, .001};
    }

    else if (modID>=6){
      return std::array<double, 3> {.512, .6144, .001};
    }
  }

  // TVector3 GetModPostion(int modID){}

  //TODO: maybe make func to pull position of each from reference files in case of changes in measurements
  //see db_sbs.gemFT.dat  
  //Positions in internal GEM coordinate system?
  

  std::array<double, 2> GetUVang(int modID){
    if(modID ==0){ return std::array<double, 2>{180., 135.};}
    else if(modID ==1){ return std::array<double, 2>{180., -135.};}
    else if (modID >=2 && modID <= 5){return std::array<double, 2>{150., -150.};}
    else if (modID >=6 && modID < 14){return std::array<double, 2>{0., -90.};}
  }
  
  std::array<int, 2> GetNstrips(int modID){
    if(modID ==0 || modID == 1){ 
      return std::array<int, 2>{3968, 3456};}

    else if (modID >=2 && modID < 6){
      return std::array<int, 2>{3840, 3840};}
    
    else if (modID >=6 && modID < 14){
      return std::array<int, 2>{1280, 1536};}
  }

  std::array<int, 2> GetNAPVs(int modID){
      
      //128 strips per APV
      // return GetNstrips(modID)/128;
      auto currStripNums = GetNstrips(modID);
      return std::array<int, 2>{currStripNums[0]/128, currStripNums[1]/128};
      
      // NOTE: if(modID ==0 || modID == 1){ 
      //   return TVector2(31, 27);}
  
      // NOTE: else if (modID >=2 && modID < 5){
      //   return TVector2(30, 30);}
      
      // NOTE: else if (modID >=6 && modID < 14){
      //   return TVector2(10, 12);}
    }


void fillGEMInfoMap(){
    for(int i=0; i<nMods; i++){
  
      gemInfoMap[i].layer = GetLayerOfMod(i);
  
      gemInfoMap[i].apvmap = GetMod_apvmap(i);
  
      gemInfoMap[i].position = ModPositionmap[i];
  
      gemInfoMap[i].size = GetModDimensions(i);
  
      gemInfoMap[i].uvangles = GetUVang(i);
      
      gemInfoMap[i].uvoffsets= {stripOffset, stripOffset};
      
      gemInfoMap[i].nstripsuv=GetNstrips(i);
  
      gemInfoMap[i].NuvAPVs = GetNAPVs(i);
  
    }
  
  }
  
  
  void printGEMinfoMap() {
    std::cout << "\nPrinting APV Info Map:\n";
    for (const auto& [key, val] : gemInfoMap) {
        std::cout << "modID: " << key << ", Layer: " << val.layer
        << ", Pos: " << val.position[0] <<" "<<val.position[1] << " " << val.position[2] <<std::endl;
        std::cout << ", Size: " << val.size[0] <<" "<<val.size[1] << " " << val.size[2]<< std::endl;
  
        std::cout << " Uangle: " << val.uvangles[0] << ", Vangle: " << val.uvangles[1] <<std::endl;
  
        std::cout<< "NstripsU: " << val.nstripsuv[0]
        << ", NstripsV: " << val.nstripsuv[1] <<
        std::endl;
        
        std::cout<< "NuAPVs: " << val.NuvAPVs[0]
        << ", NuAPVs: " << val.NuvAPVs[1] << std::endl;
        
        std::cout << "------------------------------------\n" <<std::endl;
    }
  }
  
  // Function to write `apvInfoMap` to a file
  void OutputGEMinfoMap() {
    std::ofstream mapFile("gemFT_Mod_Info.txt");
    if (!mapFile.is_open()) {
        std::cerr << "Error: Could not create file gemFT_Mod_Info.txt\n";
        return;
    }

    mapFile << "Lists helpful values for each GEM module accessible in gemInfoMap struct\n";
    mapFile << "Generated by GetGemInfoMap.cxx\n";
  
    // Set column widths for formatting
    int colWidth = 12;  // Adjust to align properly
  
    mapFile << std::left 
      << std::setw(colWidth) << "modID"
      << std::setw(colWidth) << "layer"
      << std::setw(22) << "Pos"  // Extra width for coordinates
      << std::setw(26) << "Size"  // Extra width for size
      << std::setw(colWidth) << "Uangle"
      << std::setw(colWidth) << "Vangle"
      << std::setw(colWidth) << "NstripsU"
      << std::setw(colWidth) << "NstripsV"
      << std::setw(colWidth) << "NuAPVs"
      << std::setw(colWidth) << "NvAPVs"
      << "\n";
  
    mapFile << std::string(140, '-') << "\n"; 
  
    for (const auto& [key, val] : gemInfoMap) {
      std::ostringstream posStream, sizeStream;
      posStream << "(" << val.position[0] << ", " << val.position[1] << ", " << val.position[2] << ")";
      sizeStream << "(" << val.size[0] << ", " << val.size[1] << ", " << val.size[2] << ")";
  
      mapFile << std::left 
        << std::setw(colWidth) << key
        << std::setw(colWidth) << val.layer
        << std::setw(22) << posStream.str()   
        << std::setw(26) << sizeStream.str() 
        << std::setw(colWidth) << val.uvangles[0]
        << std::setw(colWidth) << val.uvangles[1]
        << std::setw(colWidth) << val.nstripsuv[0]
        << std::setw(colWidth) << val.nstripsuv[1]
        << std::setw(colWidth) << val.NuvAPVs[0]
        << std::setw(colWidth) << val.NuvAPVs[1]
        << "\n";
    }
  
    // std::cout << "Parsed data written to gemFT_Map_TEST.txt\n";
    std::cout << "Parsed data written to gemFT_Mod_Info.txt\n";
  }
  

  void fillGEMInfoMap_fromFile(const std::string& filename = "db_FT_local.dat") {
    std::ifstream infile(filename);
    if (!infile) {
      std::cerr << "Error: Could not open file " << filename << "\n";
      return;
    }
  
    std::string line;
    while (std::getline(infile, line)) {
      if (line.empty() || line[0] == '#') continue;
  
      std::istringstream iss(line);
      std::string key;
      iss >> key;
  
      if (key.substr(0, 1) != "m") continue;
  
      // Example: m0.uangle 180.0
      size_t dot = key.find('.');
      if (dot == std::string::npos) continue;
  
      std::string modStr = key.substr(1, dot - 1);
      std::string field = key.substr(dot + 1);
      int modID = std::stoi(modStr);
  
      auto& info = gemInfoMap[modID];
  
      if (field == "layer") {
        iss >> info.layer;
      } else if (field == "uangle") {
        iss >> info.uvangles[0];
      } else if (field == "vangle") {
        iss >> info.uvangles[1];
      } else if (field == "uoffset") {
        iss >> info.uvoffsets[0];
      } else if (field == "voffset") {
        iss >> info.uvoffsets[1];
      } else if (field == "nstripsu") {
        iss >> info.nstripsuv[0];
      } else if (field == "nstripsv") {
        iss >> info.nstripsuv[1];
      } else if (field == "size") {
        iss >> info.size[0] >> info.size[1] >> info.size[2];
      } else if (field == "position") {
        iss >> info.position[0] >> info.position[1] >> info.position[2];
      }
    }
  
    // Now compute number of APVs per axis
    for (auto& [modID, info] : gemInfoMap) {
      info.apvmap = GetMod_apvmap(modID); // Still using your existing logic
      info.NuvAPVs[0] = info.nstripsuv[0] / fN_APV25_CHAN;
      info.NuvAPVs[1] = info.nstripsuv[1] / fN_APV25_CHAN;
    }
  }
  
    
    
void TestGEMInfoMap(){
  // APVFinder *anAPVFinder = new APVFinder();

  // fillGEMInfoMap();
  fillGEMInfoMap();
  printGEMinfoMap();
  OutputGEMinfoMap();
}

// void GetGemInfoMap(){
std::map<int, gemInfo> GetGemInfoMap(){

    // TestGEMInfoMap();
    // fillGEMInfoMap();
    fillGEMInfoMap_fromFile();
    OutputGEMinfoMap();
    return gemInfoMap;
}

