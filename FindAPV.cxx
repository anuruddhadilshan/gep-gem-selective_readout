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

using namespace std;
// using namespace SBSGEM;
using namespace APVFindingSpace;


  APVFinder::APVFinder() 
  {
      //Set default values for decode map parameters:
  fN_APV25_CHAN = 128;
  fMPDMAP_ROW_SIZE = 9;

  //arrays to hold raw data from one APV card:
  // fStripAPV.resize( MAXNSAMP_PER_APV );
  // fRawStripAPV.resize( MAXNSAMP_PER_APV );
  // fRawADC_APV.resize( MAXNSAMP_PER_APV );

  //default to 
  //fMAX2DHITS = 250000;
  // fMAX2DHITS = 10000;

  // fAPVmapping = SBSGEM::kUVA_XY; //default to UVA X/Y style APV mapping, but require this in the database::
  fAPVmapping = 2; 

  // fPxU = cos(UAngle);
  // fPyU = sin(UAngle);
  // fPxV = cos(VAngle);
  // fPyV = sin(VAngle);

  nMods = 13; //number of modules

  stripOffset = .004;

  fUStripOffset=.004; fVStripOffset=.004;
  ModPositionmap = {
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
  };


  void APVFinder::XY_ROI::print() const{
    std::cout << "ECalBin: " << ECalBin << std::endl;
    std::cout << "GEMLayer: " << GEMLayer << std::endl;
    std::cout << "xMin: " << xMin << " yMin: " << yMin << std::endl;
    std::cout << "xMax: " << xMax << " yMax: " << yMax << std::endl;
    std::cout << "lineNumber: " << lineNumber << std::endl;
    std::cout << "hitNumber: " << hitNumber << std::endl;
  };

  void APVFinder::UV_ROI::print() const {
    std::cout << "GEMLayer: " << GEMLayer << std::endl;
    std::cout << "MinModID: " << MinModID << std::endl;
    std::cout << "MaxModID: " << MaxModID << std::endl;
    std::cout << "Umin: " << Umin << std::endl;
    std::cout << "Vmin: " << Vmin << std::endl;
    std::cout << "Umax: " << Umax << std::endl;
    std::cout << "Vmax: " << Vmax << std::endl;
    std::cout << "uMinAPVid: " << uMinAPVid << std::endl;
    std::cout << "vMinAPVid: " << vMinAPVid << std::endl;
    std::cout << "uMaxAPVid: " << uMaxAPVid << std::endl;
    std::cout << "vMaxAPVid: " << vMaxAPVid << std::endl;
  };

  void APVFinder::gemInfo::print() const{
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

  // Define operator< for sorting in std::set
  bool APVFinder::apvInfo::operator<(const apvInfo& other) const {
    if (gemid != other.gemid) return gemid < other.gemid;
    if (axis != other.axis) return axis < other.axis;
    return pos < other.pos;
  }

  void APVFinder::apvInfo::print() const {
    std::cout << "GEMId: " << gemid << std::endl;
    std::cout << "Axis: " << axis << std::endl;
    std::cout << "Pos: " << pos << std::endl;
    std::cout << "invert: " << invert << std::endl;
    std::cout << "VTPcrate: " << vtpcrate << std::endl;
    std::cout << "fiber: " << fiber << std::endl;
    std::cout << "adc_ch: " << adc_ch << std::endl;
  };

  bool APVFinder::apvInfoKeys::operator<(const apvInfoKeys& other) const {
    if (gemid != other.gemid) return gemid < other.gemid;
    if (axis != other.axis) return axis < other.axis;
    return pos < other.pos;
}

void APVFinder::apvInfoKeys::print() const {
  std::cout 
  << "GEMId: " << gemid 
  << " Axis: " << axis 
  << " Pos: " << pos 
  << std::endl;
};

void APVFinder::apvInfoVals::print() const {
  std::cout << "VTPcrate: " << vtpcrate 
  << " Fiber: " << fiber 
  << " ADC_ch: " << adc_ch 
  << " invert? " << invert 
  << std::endl; 
};

//   void APVFinder::printAPVinfoMap() {
//     std::cout << "Printing APV Info Map:\n";
//     for (const auto& [key, val] : apvInfoMap) {
//         std::cout << "GEMId: " << key.gemid << ", Axis: " << key.axis
//                   << ", Pos: " << key.pos << std::endl;
//         std::cout << "VTPcrate: " << val.vtpcrate << ", Fiber: " << val.fiber
//                   << ", ADC_ch: " << val.adc_ch << std::endl;
//         std::cout << "------------------------------------\n";
//     }
// };



  bool APVFinder::stripInfo::operator<(const stripInfo& other) const {
    if (stripNumber != other.stripNumber) return stripNumber < other.stripNumber;
    if (axis != other.axis) return axis < other.axis;
  };

  void APVFinder::stripInfo::print() const {
    std::cout << "\nStripNumber: " << stripNumber << std::endl;
    std::cout << "Axis: " << axis << std::endl;
    std::cout << "APV_ch: " << APV_ch << std::endl;
    std::cout << "MPD ID: " << mpd_id << std::endl;
    std::cout << "ADC ID: " << adcID << std::endl;
  };

  void APVFinder::roiAPVinfo::print() const {
    std::cout << "\nuMaxAPVinfoKeys: " << std::endl;
    uMaxAPVinfoKeys.print();
    std::cout << "\nuMaxAPVinfoVals: " << std::endl;
    uMaxAPVinfoVals.print();

    std::cout << "\nuMinAPVinfoKeys: " << std::endl;
    uMinAPVinfoKeys.print();
    std::cout << "\nuMinAPVinfoVals: " << std::endl;
    uMinAPVinfoVals.print();
    
    std::cout << "\nvMaxAPVinfoKeys: " << std::endl;
    vMaxAPVinfoKeys.print();
    std::cout << "\nvMaxAPVinfoVals: " << std::endl;
    vMaxAPVinfoVals.print();

    std::cout << "\nvMinAPVinfoKeys: " << std::endl;
    vMinAPVinfoKeys.print();
    std::cout << "\nvMinAPVinfoVals: " << std::endl;
    vMinAPVinfoVals.print();
  };
  //   std::cout << "uMaxAPVinfoKeys: " << std::endl;
  //   uMaxAPVinfoKeys.print();
  //   std::cout << "\nuMinAPVinfoKeys: " << std::endl;
  //   uMinAPVinfoKeys.print();
  //   std::cout << "\nvMaxAPVinfoKeys: " << std::endl;
  //   vMaxAPVinfoKeys.print();
  //   std::cout << "\nvMinAPVinfoKeys: " << std::endl;

  //   std::cout << "uMaxAPVinfoVals: " << std::endl;
  //   uMaxAPVinfoVals.print();
  //   std::cout << "\nuMinAPVinfoVals: " << std::endl;
  //   uMinAPVinfoVals.print();
  //   std::cout << "\nvMaxAPVinfoVals: " << std::endl;
  //   vMaxAPVinfoVals.print();
  //   std::cout << "\nvMinAPVinfoVals: " << std::endl;
  //   vMinAPVinfoVals.print();
  // };



int APVFinder::GetLayerOfMod(int modID){
  if (modID >= 0 && modID < 6) { return modID; }
    else if (modID >= 6 && modID < 10) { return 6; }
    else if (modID >= 10 && modID < 14) { return 7; }
    return -1;
}

int APVFinder::GetMod_apvmap(int modID){
  if(modID>=0 && modID<6){ return 2;}
  else{return 1;}
};



void APVFinder::fillGEMInfoMap(){
  for(int i=0; i<=nMods; i++){

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


void APVFinder::printGEMinfoMap() {
  std::cout << "\nPrinting APV Info Map:\n";
  for (const auto& [key, val] : APVFinder::gemInfoMap) {
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
void APVFinder::OutputGEMinfoMap() {
  std::ofstream mapFile("gemFT_Map_TEST.txt");
  if (!mapFile.is_open()) {
      std::cerr << "Error: Could not create file gemFT_Map_TEST.txt\n";
      return;
  }

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

  for (const auto& [key, val] : APVFinder::gemInfoMap) {
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

  std::cout << "Parsed data written to gemFT_Map_TEST.txt\n";
}

//changes invert to either 1 or -1 to be used as a factor
int APVFinder::SetInvert(int invert){
  if (invert == 0){return 1;}
  if (invert == 1){return -1;}
}



//Takes in txt files
void APVFinder::MakeRefMap(const char* refFile) {
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

          currAPVvals.invert = SetInvert(currAPVvals.invert);

          apvInfoMap[currAPVkeys] = currAPVvals;
      } else {
          currentM = -1;  // Reset if empty line encountered
      }
  }

  std::cout << "Finished filling apvInfoMap with " << apvInfoMap.size() << " entries.\n";
}

  // Function to print `apvInfoMap`
void APVFinder::printAPVinfoMap() {
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
void APVFinder::OutputAPVinfoMap() {
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


void APVFinder::printRoiAPVMap() {
  std::cout << "Printing ROI APV Map:\n";
  for (const auto& [key, val] : roiAPVmap) {

      std::cout << "LineNumber: " << key << std::endl;

      
      val.print();
  }
};

// Function to write `apvInfoMap` to a file
void APVFinder::OutputRoiAPVMap() {
  std::ofstream mapFile("ROI_APV_Map_TEST.txt");
  if (!mapFile.is_open()) {
      std::cerr << "Error: Could not create file ROI_APV_Map_TEST.txt\n";
      return;
  }

  mapFile << "# roiAPVmap contains mapping of ROIs to APV information\n";
  mapFile << "## Each line represents an ROI with its corresponding APV IDs and VTP crate, fiber, and ADC channel info\n\n";

  // Set column widths for formatting
  int colWidth = 12;

  mapFile << std::left << std::setw(colWidth) << "LineNum"
          // << std::setw(colWidth) << "GEMLayer"
          << std::setw(colWidth) << "MinModID"
          << std::setw(colWidth) << "MaxModID"
          << std::setw(colWidth) << "UminAPV"
          << std::setw(colWidth) << "VminAPV"
          << std::setw(colWidth) << "UmaxAPV"
          << std::setw(colWidth) << "VmaxAPV"
          << std::setw(colWidth) << "VTPcrate"
          << std::setw(colWidth) << "Fiber"
          << std::setw(colWidth) << "ADC_ch"
          << std::setw(colWidth) << "invert"
          << "\n";

  mapFile << std::string(11 * colWidth, '-') << "\n"; // Create a line separator

  for (const auto& [lineNum, roiInfo] : roiAPVmap) {
      mapFile << std::left
              << std::setw(colWidth) << lineNum
              << std::setw(colWidth) << roiInfo.uMinAPVinfoVals.vtpcrate
              << std::setw(colWidth) << roiInfo.uMinAPVinfoVals.fiber
              << std::setw(colWidth) << roiInfo.uMinAPVinfoVals.adc_ch
              << std::setw(colWidth) << roiInfo.uMinAPVinfoVals.invert
              << std::setw(colWidth) << roiInfo.uMaxAPVinfoVals.vtpcrate
              << std::setw(colWidth) << roiInfo.uMaxAPVinfoVals.fiber
              << std::setw(colWidth) << roiInfo.uMaxAPVinfoVals.adc_ch
              << std::setw(colWidth) << roiInfo.uMaxAPVinfoVals.invert
              << std::setw(colWidth) << roiInfo.vMinAPVinfoVals.vtpcrate
              << std::setw(colWidth) << roiInfo.vMinAPVinfoVals.fiber
              << std::setw(colWidth) << roiInfo.vMinAPVinfoVals.adc_ch
              << std::setw(colWidth) << roiInfo.vMinAPVinfoVals.invert
              << std::setw(colWidth) << roiInfo.vMaxAPVinfoVals.vtpcrate
              << std::setw(colWidth) << roiInfo.vMaxAPVinfoVals.fiber
              << std::setw(colWidth) << roiInfo.vMaxAPVinfoVals.adc_ch
              << std::setw(colWidth) << roiInfo.vMaxAPVinfoVals.invert
              << "\n";
  }

  std::cout << "Parsed data written to ROI_APV_Map_TEST.txt\n";
}

// UInt_t numLines = 0;

// void LoadLine(UInt_t currLine){
// void LoadLine(std::string LineStr){
// std::vector<HitInfo> LoadLine(std::string &LineStr){
  APVFinder::XY_ROI APVFinder::LoadLine(std::string &LineStr){

  std::stringstream lineStream(LineStr);

  APVFinder::XY_ROI hit;
  lineStream >> hit.ECalBin >> hit.GEMLayer >> hit.xMin >> hit.xMax >> hit.yMin >> hit.yMax;

  if (std::isnan(hit.xMin) || std::isnan(hit.xMax) || std::isnan(hit.yMin) || std::isnan(hit.yMax)) {
    // std::cerr << "Skipping line with NaN values: " << LineStr << std::endl;
    hit.lineNumber = -1;  // Mark as invalid (or handle it as needed)
  };

  return hit;
}


  APVFinder::apvInfoVals APVFinder::GetAPV(int gemid, int axis, int pos){ 
    apvInfoKeys thisKey;
    thisKey.gemid = gemid;
    thisKey.axis = axis;
    thisKey.pos = pos;
    return APVFinder::apvInfoMap[thisKey];
  }


  
  std::array<double, 3> APVFinder::GetModDimensions(int apvmap){
    if (apvmap==2){
      //in meters
      return std::array<double, 3> {1.5, .4, .001};
    }
    if (apvmap==1){
      return std::array<double, 3> {.512, .6144, .001};
    }
  }

  // TVector3 GetModPostion(int modID){}

  //TODO: maybe make func to pull position of each from reference files in case of changes in measurements
  //see db_sbs.gemFT.dat
  //Positions in internal GEM coordinate system?
  

  std::array<double, 2> APVFinder::GetUVang(int modID){
    if(modID ==0){ return std::array<double, 2>{180., 135.};}
    else if(modID ==1){ return std::array<double, 2>{180., -135.};}
    else if (modID >=2 && modID < 5){return std::array<double, 2>{150., -150.};}
    else if (modID >=6 && modID < 14){return std::array<double, 2>{0., -90.};}
  }
  

  std::array<double, 2> APVFinder::LayerUVang(int layer){
    if(layer ==0){ return std::array<double, 2>{180., 135.};}
    else if(layer ==1){ return std::array<double, 2>{180., -135.};}
    else if (layer >=2 && layer <6){return std::array<double, 2>{150., -150.};}
    else if (layer>=6){return std::array<double, 2>{0., -90.};}
  }



  //Get Projection operators given angles u and v
  void APVFinder::SetProjOps(std::array<double, 2> UVangles){
    auto fUAngle=UVangles[0]; auto fVAngle=UVangles[1];

    APVFinder::fPxU = cos( fUAngle * M_PI/180.);
    APVFinder::fPyU = sin( fUAngle * M_PI/180.);
    APVFinder::fPxV = cos( fVAngle * M_PI/180.);
    APVFinder::fPyV = sin( fVAngle * M_PI/180.);
  }

  //.0004m = 400 micrometers
  // double GetPitch(){return .0004;}

  // double APVFinder::GetUVoffset(int modID, char *axis){
  double APVFinder::GetUVoffset(int modID, int axis){
    if(modID ==0 || modID == 1){
      if (axis == 0){ return .0176;}
      else if (axis == 1){return 0.;}
    }
    else if (modID >=2 && modID < 5){return .0108;}
    else if (modID >=6 && modID < 14){return 0.;}
  }
  
  std::array<double, 2> APVFinder::GetOffset(int modID){
    if(modID ==0 || modID == 1){ return std::array<double, 2>{.0176, 0.};}
    else if (modID >=2 && modID < 5){return std::array<double, 2>{.0108, .0108};}
    else if (modID >=6 && modID < 14){return std::array<double, 2>{0., 0.};}
  }
 
  // int APVFinder::GetNstrips(int modID, const std::string& axis){
  int APVFinder::GetNstrips(int modID, int axis){
    if(modID ==0 || modID == 1){
      if (axis == 0){ return 3968;}
      else if (axis == 1){return 3456;}

    else if (modID >=2 && modID < 5){return 3840;}
    
    else if (modID >=6 && modID < 14){
      if (axis == 0){ return 1280;}
      else if (axis == 1){return 1536;}
    }
    }
  }

  // Int_t SBSGEMModule::GetStripNumber( UInt_t rawstrip, UInt_t pos, UInt_t invert ){
  int APVFinder::GetStripNumber( int rawstrip, int pos, int invert ){
    int RstripNb = APVMAP[fAPVmapping][rawstrip];
    RstripNb = RstripNb + (127-2*RstripNb)*invert;
    int RstripPos = RstripNb + 128*pos;
  
    // ?????????????????
    // if( fIsMC ){
    //   return rawstrip + 128*pos;
    // }
    
    return RstripPos;
  }
  
  std::array<int, 2> APVFinder::GetNstrips(int modID){
    if(modID ==0 || modID == 1){ 
      return std::array<int, 2>{3968, 3456};}

    else if (modID >=2 && modID < 6){
      return std::array<int, 2>{3840, 3840};}
    
    else if (modID >=6 && modID < 14){
      return std::array<int, 2>{1280, 1536};}
  }


  std::array<int, 2> APVFinder::GetNAPVs(int modID){
      
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
    
    int APVFinder::UorVtoAPVid(int modID, std::array<int, 2> localUVPos,
      int axis //U/X=0, V/Y=1
    ){
      double localPos;
      
      double APVsize = 128*GetPitch();

      int numAPVs;

      if (axis == 0){
        numAPVs = int(GetNAPVs(modID)[0]);
        localPos=localUVPos[0];}
      else if (axis == 1){
        numAPVs = int(GetNAPVs(modID)[1]);
        localPos=localUVPos[1];} 
    
      //numAPVs along this axis

      // TVector2 APVedges

      int layerAPVid;

      //NOTE: PER Layer APV id;
      if(modID >= 0 && modID < 6){
        for( int i = 0; i <= numAPVs; i++){
          // APVedges= TVector2(i, (i+1))*APVsize;
          if (localPos>=i*APVsize && 
            localPos<(i+1)*APVsize)
            { layerAPVid = i; break;}
         }
      }
      else if(modID >= 6 && modID < 10){
        modID -= 6;
        for( int i = 0; i <= numAPVs; i++)
          {if (localPos>=i*APVsize && 
            localPos<(i+1)*APVsize)
            { layerAPVid = i; break;}}
      }
      else if (modID >= 10 && modID < 14){
        modID-=10;
        for( int i = 0; i <= numAPVs; i++)
          {if (localPos>=i*APVsize && 
            localPos<(i+1)*APVsize)
            { layerAPVid = i; break;}}
       };
      

      int currMod = 0;
      int globalAPVid=layerAPVid;

      while(currMod < modID) {
        int currNAPVs;
        if (axis == 0){
          currNAPVs = int(GetNAPVs(currMod)[0]);
          localPos=localUVPos[0];}
        else if (axis == 1){
          currNAPVs = int(GetNAPVs(currMod)[1]);
          localPos=localUVPos[1];} 

        globalAPVid += currNAPVs;
        currMod++;
      }
      return globalAPVid;
      
    }


  // TODO: for now just copied from SBSGEMModule 
  std::array<double, 2> APVFinder::XYtoUV( std::array<double, 2> XY ){

    double Xtemp = XY[0];
    double Ytemp = XY[1];
  
    double Utemp = Xtemp*fPxU + Ytemp*fPyU;
    double Vtemp = Xtemp*fPxV + Ytemp*fPyV;
  
    return std::array<double, 2> {Utemp,Vtemp};
  }

  //transforms with respect to the center of the layer 
  //count down from top dividing up mods along each layer 


  // int APVFinder::GetModId(int layerID, std::array<int, 2> hitPos, std::array<int, 3> modDims){
  int APVFinder::GetModId(int layerID, double hitX){

    std::cout << "\n#Calling GetModId(layer, xGlobal): " << std::endl;

    double xSize=APVFinder::GetModDimensions(layerID)[0];

    std::cout << "xSize: " << xSize << std::endl;

    std::array<double, 3> modPos;
    int modID = -1;
    if(layerID > 5){
      for (int i = 6; i < 14; i++){
          std::array<double, 3> modPos = ModPositionmap[i];
          modPos[0] *= -1; //flip to match GEM coordinate system

          // std::cout << "modID: " << i << " modPos: " << modPos[0] << std::endl;
          // modPod[0] += 2*xSize; //make top 
          // give 2 strip leeway:
          double highEnd = modPos[0] + xSize/2;
          double lowEnd = modPos[0] - xSize/2;
          // std::cout << "lowEnd: " << lowEnd << " highEnd: " << highEnd << std::endl;
          // if (hitX >= modPos[0] - xSize/2 && hitX <= modPos[0] + xSize/2)
          if (hitX >= lowEnd && hitX <= highEnd)
            {modID = i; break;}

            // if (hitX >= modPos[0]*((i-6)%4) - 10*GetPitch() && hitX <= modPos[0]*((i-6)%4) + 10*GetPitch())

      }
    if (modID == -1){
        std::cerr << "Error: Can't find modID for hitX=" << hitX << std::endl;
        exit(-1);
    }
    } else {
      // For layers where a single module is expected, assign the layer ID (or reconsider this logic)
      modID = layerID;
    }
    std::cout << "xGlobal:" << hitX << " modID: " << modID << std::endl;
      return modID;
  }



  APVFinder::UV_ROI APVFinder::XYtoUV_ROI(APVFinder::XY_ROI hitXY){
 
    std::cout << "\n#Calling XYtoUV_ROI(): " << std::endl;


    UV_ROI thisUV_ROI;

    int layerID = hitXY.GEMLayer; 
    thisUV_ROI.GEMLayer = layerID;

    std::cout << "LayerID: " << layerID << std::endl;
    
    std::array<int, 2> hitUV;
    
    // int modID;  
    
    std::array<double , 3> modDims;


    int apvmap; //type of APV config/module type
    // if apvmap= 0; //INFN
    // if apvmap= 1; //UVA XY
    // if apvmap= 2; //UVA U/V(layers 2-5) or X/W(layers 0-1)
  
    //NOTE: Layers 0 and 1 are UVA X/W
    // if (layerID in range(0,6)){apvmap=2;}
    if (layerID >=0 && layerID<6 )
    // {apvmap=2; std::array<double, 3>modDims {1.5, .4, .001};}
    {apvmap=2; 
      modDims = {1.5, .4, .001};}
    // else if (layerID >=6 && layerID<14){apvmap=1; std::array<double, 3>modDims {.512, .6144, .001};};
    else if (layerID >=6 && layerID<14){apvmap=1; 
      modDims = {.512, .6144, .001};};

    std::cout << "APVMap(style): " << apvmap << std::endl;
    std::cout << "ModDims: " << modDims[0] << " " << modDims[1] << " " << modDims[2] << std::endl;


    // std::cout << "\nminX: " << hitXY.xMin //<< std::endl; std::cout 
    // << " minXmodID: " << GetModId(layerID, hitXY.xMin) << std::endl;

    // std::cout << "maxX: " << hitXY.xMax //<< std::endl; std::cout 
    // << " maxXmodID: " << GetModId(layerID, hitXY.xMax) << std::endl;
    
    // std::cout << "minY: " << hitXY.yMin //<< std::endl; std::cout 
    // << " minYmodID: " << GetModId(layerID, hitXY.yMin) << std::endl;

    // std::cout << "maxY: " << hitXY.yMax //<< std::endl; std::cout 
    // << " maxXmodID: " << GetModId(layerID, hitXY.yMax) << std::endl;


    
    // auto localX=hitPos[0]-mod[modID].position;
    //relative to center of module

    
    //Take 4 corners of the XY ROI
    std::array<double, 2>minXminY {hitXY.xMin, hitXY.yMin};
    std::array<double, 2>maxXmaxY {hitXY.xMax, hitXY.yMax};
    std::array<double, 2>minXmaxY {hitXY.xMin, hitXY.yMax};
    std::array<double, 2>maxXminY {hitXY.xMax, hitXY.yMin};

    
    // thisUV_ROI.MinModID = GetModId(layerID, minXY, modDims);
    // thisUV_ROI.MaxModID = GetModId(layerID, maxXY, modDims);
    thisUV_ROI.MinModID = GetModId(layerID, minXminY[0]);
    thisUV_ROI.MaxModID = GetModId(layerID, maxXmaxY[1]);

    std::cout<< "\nMinModID(corresponding with xCoord of minXminY): " << thisUV_ROI.MinModID << std::endl;
    std::cout<< "MaxModID(corresponding with yCoord of maxXmaxY): " << thisUV_ROI.MaxModID << std::endl;

    
    // std::array<int, 2> currUVangs(GetUVang(layerID));
    std::array<double, 2>currUVangs (APVFinder::GetUVang(thisUV_ROI.MinModID));

    std::cout << "\nUV Angles(For MinModID): " << currUVangs[0] << " " << currUVangs[1] << std::endl;
    
    SetProjOps(currUVangs);

    //TODO: make local XYtoUV
    // thisUV_ROI.Umin = XYtoUV(maxXminY)[0];
    // thisUV_ROI.Vmin = XYtoUV(maxXmaxY)[1];
    // thisUV_ROI.Umax = XYtoUV(minXmaxY)[0];
    // thisUV_ROI.Vmax = XYtoUV(minXminY)[1];

    // std::cout << "\nUmin: " << thisUV_ROI.Umin << " Vmin: " << thisUV_ROI.Vmin << std::endl;
    // std::cout << "Umax: " << thisUV_ROI.Umax << " Vmax: " << thisUV_ROI.Vmax << std::endl;

    //Turn Coordinates into APV positions along their axes
    // thisUV_ROI.UminPos = GetAPVpos(thisUV_ROI.Umin);
    // thisUV_ROI.VminPos = GetAPVpos(thisUV_ROI.Vmin);
    // thisUV_ROI.UmaxPos = GetAPVpos(thisUV_ROI.Umax);
    // thisUV_ROI.VmaxPos = GetAPVpos(thisUV_ROI.Vmax);
    
    thisUV_ROI.uMinAPVid = GetAPVpos(thisUV_ROI.MinModID, hitXY.xMax, hitXY.yMin, 0); 
    thisUV_ROI.vMinAPVid = GetAPVpos(thisUV_ROI.MinModID, hitXY.xMax, hitXY.yMax, 1);
    thisUV_ROI.uMaxAPVid = GetAPVpos(thisUV_ROI.MinModID, hitXY.xMin, hitXY.yMax, 0); 
    thisUV_ROI.vMaxAPVid = GetAPVpos(thisUV_ROI.MinModID, hitXY.xMin, hitXY.yMin, 1);

    std::cout << "\nuMinAPVid: " << thisUV_ROI.uMinAPVid << std::endl;
    std::cout << "vMinAPVid: " << thisUV_ROI.vMinAPVid << std::endl;
    std::cout << "uMaxAPVid: " << thisUV_ROI.uMaxAPVid << std::endl;
    std::cout << "vMaxAPVid: " << thisUV_ROI.vMaxAPVid << std::endl;
    
   
    // TVector2 thisUV_ROI=XYtoUV(hitPos);  

    //NOTE: UV_ROI should now have
    // layerid, modid, lineNumber(from earlier in FindAPV(), 
    // Umin, Vmin, Umax, Vmax
    return thisUV_ROI;
  }


  // int APVFinder::GetAPVpos(int layerID, int axis, double Coord){

  //local u/v position in module
  // int APVFinder::GetAPVpos(UV_ROI theUVROI){

  // std::array<double, 2> APVFinder::GetXYLocal(XY_ROI theXYROI){
  std::array<double, 2> APVFinder::GetXYLocal(int layerID, double xGlobal, double yGlobal){

    //NOTE: x
    
    std::cout << "\n#Calling GetXYLocal(): " << std::endl;
    std::cout << "GlobalXY: " << xGlobal << " " << yGlobal << std::endl;

    int thisModID = GetModId(layerID, xGlobal);

    // **Check for Invalid Module ID**
    if (thisModID < 0 || thisModID >= gemInfoMap.size()) {
        std::cerr << "Error: Invalid modID (" << thisModID << ") for xGlobal: " << xGlobal << "\n";
        return {NAN, NAN}; // Return NaN to indicate failure
    }

    // Get module properties
    double modX = gemInfoMap[thisModID].position[0]; // Module center X
    double modY = gemInfoMap[thisModID].position[1]; // Module center Y
    double modLength = gemInfoMap[thisModID].size[0]; // Module width

    // **Shift xGlobal so that right edge is at x = 0**
    double localX = xGlobal - (modX + modLength / 2); // Right edge of mod is x=0
    double localY = yGlobal - modY; // Keep local Y relative to module center

    std::cout << "LocalXY: " << localX << " " << localY << std::endl;

    return {localX, localY};

    // std::cout << "\n#Calling GetXYLocal(): " << std::endl;
    // std::cout << "GlobalXY: " << xGlobal << " " << yGlobal << std::endl;

    // //TODO: shift from 0(center of layer) to top of mod
    // int thisModID = GetModId(layerID, xGlobal);

    // double localX;
    // double localY = yGlobal - gemInfoMap[thisModID].position[1];

    // if (thisModID < 6){ 
    //   localX = xGlobal + gemInfoMap[thisModID].size[0]/2;}

    // else if (thisModID >= 6){
    //   //NOTE: 2 mods = 2 whole mods
    //   // return theXYROI.x - 2*(2*gemInfoMap[thisModID].size[0]);}
    //   localX = xGlobal + 2*(gemInfoMap[thisModID].size[0]);}

    //   std::cout << "LocalXY: " << localX << " " << localY << std::endl;

    // return std::array<double, 2>{localX, localY};

    }
  
    double APVFinder::GetXstripCenterLocal(int layerID, double xGlobal){

    //NOTE: x
    
    std::cout << "\n#Calling GetXYLocal(): " << std::endl;
    std::cout << "GlobalX: " << xGlobal  << std::endl;

    int thisModID = GetModId(layerID, xGlobal);

    // **Check for Invalid Module ID**
    if (thisModID < 0 || thisModID >= gemInfoMap.size()) {
        std::cerr << "Error: Invalid modID (" << thisModID << ") for xGlobal: " << xGlobal << "\n";
        // return {NAN, NAN}; // Return NaN to indicate failure
        return NAN;
    }

    // Get module properties
    double modX = gemInfoMap[thisModID].position[0]; // Module center X
    double modY = gemInfoMap[thisModID].position[1]; // Module center Y
    double modLength = gemInfoMap[thisModID].size[0]; // Module width

    double modTopPos = (modX)+(modLength/2); //top edge of module

     
    // **Shift xGlobal so that top edge is at x = 0 of module **
    double localX = xGlobal + modTopPos;  

    double stripCenter = localX+(GetUdx(GetModId(layerID, xGlobal)/2)); 
    
    
    std::cout << "LocalX: " << localX << std::endl;

    return localX;

    // std::cout << "\n#Calling GetXYLocal(): " << std::endl;
    // std::cout << "GlobalXY: " << xGlobal << " " << yGlobal << std::endl;

    // //TODO: shift from 0(center of layer) to top of mod
    // int thisModID = GetModId(layerID, xGlobal);

    // double localX;
    // double localY = yGlobal - gemInfoMap[thisModID].position[1];

    // if (thisModID < 6){ 
    //   localX = xGlobal + gemInfoMap[thisModID].size[0]/2;}

    // else if (thisModID >= 6){
    //   //NOTE: 2 mods = 2 whole mods
    //   // return theXYROI.x - 2*(2*gemInfoMap[thisModID].size[0]);}
    //   localX = xGlobal + 2*(gemInfoMap[thisModID].size[0]);}

    //   std::cout << "LocalXY: " << localX << " " << localY << std::endl;

    // return std::array<double, 2>{localX, localY};

    }

// bool APVFinder::aboveStripCenter(modID, xGlobal, yGlobal){
//   if (yGlobal<0) 
//   GetXstripCenterLocal(layerID, xGlobal)% GetUdx(modID)
// }

  // int APVFinder::GetAPVpos(modID, posXY){
  // int APVFinder::GetAPVpos(int layerID, double xGlobal, double yGlobal, int axis){
  int APVFinder::GetAPVpos(int modID, double xGlobal, double yGlobal, int axis){

    std::cout << "\n#Calling GetAPVpos(): " << std::endl;
    std::cout << "modID: " <<   modID << std::endl;

    // int layerID = GetLayerOfMod(modID);
    int layerID = gemInfoMap[modID].layer;
    std::cout << "LayerID: " << layerID << std::endl;

    SetProjOps(LayerUVang(modID));
    // std::cout << "SetProjOps(" << LayerUVang.print() << ")" << std::endl;

    // std::array<double, 2> localXY = GetXYLocal(layerID, xGlobal, yGlobal);
    std::array<double, 2> localXY = GetXYLocal(layerID, xGlobal, yGlobal);   
    
    
    
    // if (GetXstripCenterLocal(layerID, xGlobal)% GetUdx(modID) < 0){//if Xlocal is>strip center
    //   int localStrip;
    //   localStrip = localXY[0]/GetUdx(modID);
    // }
    // // else if (GetXstripCenterLocal(layerID, xGlobal)% GetUdx(modID) < 0){//if Xlocal is>strip center
    
    // int localStrip;

    // if (axis == 0){
    //   std::cout << "localStrip in UV coords: " << localXY[0]/GetUdx(modID) << std::endl;
    //   localStrip = int(localXY[0]/GetUdx(modID));}
    // else if (axis == 1){
    //   std::cout << "localStrip in UV coords: " << localXY[0]/GetVdx(modID) << std::endl;
    //   localStrip = int(localXY[1]/GetVdx(modID));
    //   }
    
    // if (hitX < hitCenter && hitY>0){ localStrip+=1;}
    
    // else if (hitX >= hitCenter && hitY<=0){localStrip-=1;}


    int localStrip;
    if (axis == 0){
      std::cout << "localStrip in UV coords: " << localXY[0]/GetUdx(modID) << std::endl;
      localStrip = int(localXY[0]/GetUdx(modID));}
    else if (axis == 1){
      std::cout << "localStrip in UV coords: " << localXY[0]/GetVdx(modID) << std::endl;
      localStrip = int(localXY[1]/GetVdx(modID));

      std::cout << "LocalStripid: " << localStrip << std::endl;

      int APVpos = localStrip%128;
      std::cout << "APVpos: " << APVpos << std::endl;
    
    return APVpos;
    // return localStrip%128;

    // if (axis == 0){ //U/X
    //   return int((localUV[0]+GetUVoffset(layerID, axis))/GetPitch());}
      
    // else if (axis == 1){ //V/Y
    //   return int((localUV[1]+GetUVoffset(layerID, axis))/GetPitch());}

    }
  }

    // std::map<double, std::array<int, 2>> xToUVmap;
    // std::map<double, std::array<int, 2>> yToUVmap;

    
    // int APVFinder::GetAPVpos(int layerID){
    //   return std::array<int, 2>{GetAPVpos(layerID, 0), GetAPVpos(layerID, 1)};
    // }

  // bool APVFinder::isAboveStripCenter(int modID, double xLocal){
  //   if (xLocal % GetUdx(modID) > 0) {return true;}
  //   else {return false;};
  // }


    double APVFinder::GetUdx(int modID){
      // SetProjOps(GetUVang(modID));
      std::cout<< "\n#Calling GetUdx(modID): " << modID << std::endl;
      Udy=gemInfoMap[modID].size[1];
      std::cout<< "Udy(size): " << Udy << std::endl;

      double Uang = gemInfoMap[modID].uvangles[0];
      Udx=(Udy*sin(Uang))/cos(Uang);
      std::cout<< "Udx: " << Udx << std::endl;
      return Udx;
    }

    double APVFinder::GetVdx(int modID){
      // SetProjOps(GetUVang(modID));
      std::cout<< "\n#Calling GetVdx(modID): " << modID << std::endl;
      Vdy=gemInfoMap[modID].size[1];
      std::cout<< "Vdy(size): " << Udy << std::endl;
      double Vang = gemInfoMap[modID].uvangles[1];
      Vdx=(Vdy*sin(Vang))/cos(Vang);
      std::cout<< "Vdx: " << Udx << std::endl;
      return Vdx;
    }


    // void APVFinder::MakeXtoUVstripMap(int layerID){
      
    //   SetProjOps(LayerUVang(layerID));

    //   double Udx = 

    //   // int numStrips = GetNstrips(layerID, 0);
    //   int numStrips = GetNstrips(layerID, 0);
    //   for (int i=0; i<numStrips; i++){
    //     XtoUVmap[i] = GetAPVpos(layerID, i, 0);
    //   }
    // }

    // void APVFinder::MakeYtoUVstripMap(int layerID){
    //   // int numStrips = GetNstrips(layerID, 1);
    //   int numStrips = GetNstrips(layerID, 1);
    //   for (int i=0; i<numStrips; i++){
    //     YtoUVmap[i] = GetAPVpos(layerID, i, 1);
    //   }
    // }


  //"Inverse" of hitpos from SBSGEMModule::find_clusters_1D()
  // int hitposToStripID( std::array<int, 2> hitpos, int &Nstrips ){
    //NOTE: dont need b/c find APV by spatial extent instead

    //TODO: UV or XY hitpos?

    // int istrip = int(
    // (hitpos-offset)/pitch 
    // + .5*fNstrips - .5);

    // return istrip;
  // }





void APVFinder::TestGEMInfoMap(){
  // APVFinder *anAPVFinder = new APVFinder();

  fillGEMInfoMap();
  printGEMinfoMap();
  OutputGEMinfoMap();
}

void APVFinder::TestAPVInfoMap(const char* refFile){
  // APVFinder *anAPVFinder = new APVFinder();

  MakeRefMap(refFile);
  printAPVinfoMap();
  OutputAPVinfoMap();
}

// ################################################################

// FindAPV general outline
 // 1) xyHit => uvHit 
      // i) LoadLine => ECalBin GEMLayer xMin xMax yMin yMax
      //ii)(layer, xMin, xMax, yMin, yMax) => (uMin, uMax, vMin, vMax)


    // 2) uvHit => APVinfo
      //need apvInfoKeys: module id, axis, position on axis
      // i) GetModId(layer, xMin, xMax)=>minMod,maxMod
      // ii) TODO: Get INTEGER position on axis assuming "-x" is positive for u,v
      // iii) GetAPV(min/max Mod, 0, min/max UAPVnum)=> minU APV info, max U APV info, minV APV info, maxV APV info
      // iv) pushback(minU APV info, max U APV info, minV APV info, maxV APV info)

// void FindAPV(const TDatime& date, const char *aHitFile){



  //GEM and APV print tests
void FindAPV(){
  
  APVFinder *anAPVFinder = new APVFinder();
  
  anAPVFinder->TestGEMInfoMap();
  anAPVFinder->TestAPVInfoMap("db_sbs.gemFT_TEST.txt");

}
  
  
// int FindAPV(ifstream aHitFile){

// int FindAPV(const char* aHitFile){
// void FindAPV(){
//   const char* aHitFile = "ECALtoGEMhitEx.txt";

//   APVFinder *anAPVFinder = new APVFinder();
  
//   anAPVFinder->fillGEMInfoMap(); anAPVFinder->MakeRefMap("db_sbs.gemFT_TEST.txt");

//   // apvInfoMap = anAPVFinder->apvInfoMap;

//   // int chanPerAPV=APVFinder->GetstripsPerAPV();  
//   int chanPerAPV=128;


//   std::ifstream hitFile(aHitFile);  

//   if (!hitFile.is_open()) {  // Check if the file opened successfully
//       std::cerr << "Unable to open file\n";
//       // return 1;
//   }
 

//   // ReadDB(date); 
  
//   // DBfile = 

//   // makeAPVid(date);

//   std::string currLine;

//   // std::vector<LocalROIInfo> localHitMap; 

//   // skip first line(Column headers)
//   std::getline(hitFile, currLine);

//   int hitCount = 0;//0
//   int lineCount = hitCount+2;//1
  
//   while (std::getline(hitFile, currLine)) {
    
//     APVFinder::UV_ROI thisUV_ROI;

//     // 1) xyHit => uvHit 
//       // i) LoadLine => ECalBin GEMLayer xMin xMax yMin yMax
      
//       APVFinder::XY_ROI currHit = anAPVFinder->LoadLine(currLine);
//       // HitInfo currHit = LoadLine(currLine);
      
//       if (currHit.lineNumber == -1) {continue;}//skip line if it has NaN values

//       std::cout << "\n\n--------------------------------------------------------------------"<<"\nLine: " << lineCount << "  Hit Entry number " << hitCount << std::endl;
      
//       currHit.lineNumber=lineCount;//Add the txt file line number just for sake of indexing/debugging
      
//       currHit.hitNumber=hitCount;

//       currHit.print();
      
//       //ii)(layer, xMin, xMax, yMin, yMax) => (uMin, uMax, vMin, vMax)
//       thisUV_ROI = anAPVFinder->XYtoUV_ROI(currHit);


//       APVFinder::apvInfoKeys UminAPVinfoKeys;
//       UminAPVinfoKeys.gemid = thisUV_ROI.MinModID;
//       UminAPVinfoKeys.axis = 0;
//       UminAPVinfoKeys.pos = thisUV_ROI.uMinAPVid;

//       APVFinder::apvInfoKeys UmaxAPVinfoKeys;
//       UmaxAPVinfoKeys.gemid = thisUV_ROI.MaxModID;
//       UmaxAPVinfoKeys.axis = 0;
//       UmaxAPVinfoKeys.pos = thisUV_ROI.uMaxAPVid;

//       APVFinder::apvInfoKeys VminAPVinfoKeys;
//       VminAPVinfoKeys.gemid = thisUV_ROI.MinModID;
//       VminAPVinfoKeys.axis = 1;
//       VminAPVinfoKeys.pos = thisUV_ROI.vMinAPVid;

//       APVFinder::apvInfoKeys VmaxAPVinfoKeys;
//       VmaxAPVinfoKeys.gemid = thisUV_ROI.MaxModID;
//       VmaxAPVinfoKeys.axis = 1;
//       VmaxAPVinfoKeys.pos = thisUV_ROI.vMaxAPVid;

     
//       APVFinder::roiAPVinfo thisroiAPVinfo;
//       thisroiAPVinfo.uMaxAPVinfoVals = anAPVFinder->apvInfoMap[UmaxAPVinfoKeys];
//       std::cout << "\nUmaxAPVinfoKeys: " << UmaxAPVinfoKeys.gemid << " " << UmaxAPVinfoKeys.axis << " " << UmaxAPVinfoKeys.pos << std::endl;
//       std::cout << "UmaxAPVinfoVals: " << thisroiAPVinfo.uMaxAPVinfoVals.vtpcrate<<" "<< thisroiAPVinfo.uMaxAPVinfoVals.fiber << " " << thisroiAPVinfo.uMaxAPVinfoVals.adc_ch << std::endl;

//       thisroiAPVinfo.uMinAPVinfoVals = anAPVFinder->apvInfoMap[UminAPVinfoKeys];
//       std::cout << "\nUminAPVinfoKeys: " << UminAPVinfoKeys.gemid << " " << UminAPVinfoKeys.axis << " " << UminAPVinfoKeys.pos << std::endl;
//       std::cout << "UminAPVinfoVals: " << thisroiAPVinfo.uMinAPVinfoVals.vtpcrate<<" "<< thisroiAPVinfo.uMinAPVinfoVals.fiber << " " << thisroiAPVinfo.uMinAPVinfoVals.adc_ch << std::endl;

//       thisroiAPVinfo.vMaxAPVinfoVals = anAPVFinder->apvInfoMap[VmaxAPVinfoKeys];
//       std::cout << "\nVmaxAPVinfoKeys: " << VmaxAPVinfoKeys.gemid << " " << VmaxAPVinfoKeys.axis << " " << VmaxAPVinfoKeys.pos << std::endl;
//       std::cout << "VmaxAPVinfoVals: " << thisroiAPVinfo.vMaxAPVinfoVals.vtpcrate<<" "<< thisroiAPVinfo.vMaxAPVinfoVals.fiber << " " << thisroiAPVinfo.vMaxAPVinfoVals.adc_ch << std::endl;

//       thisroiAPVinfo.vMinAPVinfoVals = anAPVFinder->apvInfoMap[VminAPVinfoKeys];
//       std::cout << "\nVminAPVinfoKeys: " << VminAPVinfoKeys.gemid << " " << VminAPVinfoKeys.axis << " " << VminAPVinfoKeys.pos << std::endl;
//       std::cout << "VminAPVinfoVals: " << thisroiAPVinfo.vMinAPVinfoVals.vtpcrate<<" "<< thisroiAPVinfo.vMinAPVinfoVals.fiber << " " << thisroiAPVinfo.vMinAPVinfoVals.adc_ch << std::endl;

      
//       //store the current hits line number and APV info
//       // anAPVFinder->roiAPVmap[lineCount] = thisroiAPVinfo;
//       anAPVFinder->roiAPVmap[hitCount] = thisroiAPVinfo;

      
//     // std::cout << "\nLine: " << lineCount << std::endl;

//     // UminAPVinfoKeys.print(); UmaxAPVinfoKeys.print(); VminAPVinfoKeys.print(); VmaxAPVinfoKeys.print(); 

//     std::cout << "\nROI:______________________________" << std::endl;
//     currHit.print();
//     thisroiAPVinfo.print();

//     lineCount++;
//     hitCount++;

//   }

//   // anAPVFinder->printRoiAPVMap();
//   // anAPVFinder->OutputRoiAPVMap();

//   hitFile.close();
//   // return 0;

  
// }


