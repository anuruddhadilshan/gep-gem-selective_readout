#include <iostream>

#include "SBSGEMModule.h"
// #include "THaSubDetector.h"
// #include "ECalToHCal.h"

// #include "FindGem.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>

// #include "TObject.h"
#include "TString.h"

// #include "THaAnalysisObject.h"
#include "TVector2.h"
#include "TVector3.h"

#include "TDatime.h"
#include "THaEvData.h"
#include "THaApparatus.h"
#include "THaRun.h"
#include "TRotation.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include <algorithm>
#include <iomanip>


using namespace std;
using namespace SBSGEM;

//This should not be hard-coded, I think, but read in from the database (or perhaps not, if it never changes? For now we keep it hard-coded)
// const int APVMAP[128] = {
// 1, 33, 65, 97, 9, 41, 73, 105, 17, 49, 81, 113, 25, 57, 89, 121, 3, 35, 67, 99, 11, 43, 75, 107, 19, 51, 83, 115, 27, 59, 91, 123, 5, 37, 69, 101, 13, 45, 77, 109, 21, 53, 85, 117, 29, 61, 93, 125, 7, 39, 71, 103, 15, 47, 79, 111, 23, 55, 87, 119, 31, 63, 95, 127, 
// 0, 32, 64, 96, 8, 40, 72, 104, 16, 48, 80, 112, 24, 56, 88, 120, 2, 34, 66, 98, 10, 42, 74, 106, 18, 50, 82, 114, 26, 58, 90, 122, 4, 36, 68, 100, 12, 44, 76, 108, 20, 52, 84, 116, 28, 60, 92, 124, 6, 38, 70, 102, 14, 46, 78, 110, 22, 54, 86, 118, 30, 62, 94, 126};


// class APVFinder {
namespace APVFinder {

  //Based on SBSGEMModule

  

  // struct inverseMap invMap_t {
  //   UInt_t posGlobal;
  //   UInt_t gem_id;
  //   // UInt_t crate;
  //   // UInt_t slot;
  //   UInt_t mpd_id;
  //   UInt_t adc_id;
  //   UInt_t i2c;
  //   UInt_t invert;
  //   UInt_t axis; //needed to add axis to the decode map
  //   UInt_t index;
  // };

// struct HitInfo{
struct XY_ROI{
  UInt_t ECalBin;
  UInt_t GEMLayer;
  Double_t xMin;
  Double_t xMax;
  Double_t yMin;
  Double_t yMax;
  int lineNumber;


  void print() const{
    std::cout << "ECalBin: " << ECalBin << std::endl;
    std::cout << "GEMLayer: " << GEMLayer << std::endl;
    std::cout << "xMin: " << xMin << std::endl;
    std::cout << "yMin: " << yMin << std::endl;
    std::cout << "xMax: " << xMax << std::endl;
    std::cout << "yMax: " << yMax << std::endl;
    std::cout << "lineNumber: " << lineNumber << std::endl;
  }
};

struct UV_ROI{
  UInt_t GEMLayer;
  int lineNumber;
  
  int MinModID;
  int MaxModID;

  Double_t Umin;
  Double_t Umax;
  Double_t Vmin;
  Double_t Vmax;

  UInt_t uMinAPVid;
  UInt_t vMinAPVid;
  
  UInt_t uMaxAPVid;
  UInt_t vMaxAPVid;

  //TODO: 
  int vtpcrate; 
  int fiber; // ids which MPD on a VTPcrate
  int adc_ch; //each APV HDMI connected within an MPD

  
  void print() const {
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

    std::cout << "VTPcrate: " << vtpcrate << std::endl;
    std::cout << "fiber" << fiber << std::endl;
    std::cout << "adc_ch" << fiber;
  }
  };


std::vector<UV_ROI> ROI_APV_map;




int fN_APV25_CHAN = 128;
int fMPDMAP_ROW_SIZE = 9;
int MAXNSAMP_PER_APV;

//Arrays to temporarily hold raw data from ONE APV card:
  std::vector<UInt_t> fStripAPV;
  std::vector<UInt_t> fRawStripAPV;
  std::vector<Int_t> fRawADC_APV;

  int fMAX2DHITS;
int fAPVmapping;

// void Init_Finder(){
//   //Set default values for decode map parameters:
//   fN_APV25_CHAN = 128;
//   fMPDMAP_ROW_SIZE = 9;

//   //arrays to hold raw data from one APV card:
//   fStripAPV.resize( MAXNSAMP_PER_APV );
//   fRawStripAPV.resize( MAXNSAMP_PER_APV );
//   fRawADC_APV.resize( MAXNSAMP_PER_APV );

//   //default to 
//   //fMAX2DHITS = 250000;
//   fMAX2DHITS = 10000;

//   fAPVmapping = SBSGEM::kUVA_XY; //default to UVA X/Y style APV mapping, but require this in the database::

//   InitAPVMAP();

// }

struct gemInfo{

  int modID;
  int layer;
  int apvmap;
  TVector3 position;
  // TVector2 angle; xyz ang
  TVector3 size;

  double uangle;  double vangle;
  double uoffset;  double voffset;
  int nstripsu;  int nstripsv;

  int NuAPVs; int NvAPVs;
  
  Double_t fPxU = cos(uangle);            //U Strip X projection = cos( UAngle );
  Double_t fPyU = sin(uangle);            //U Strip Y projection = sin( UAngle );
  Double_t fPxV = cos(vangle);            //V Strip X projection = cos( VAngle );
  Double_t fPyV = sin(vangle);            //V Strip Y projection = sin( VAngle );

  // TVector3 alignedPos;  TVector3 alignedAng;

  void print(){
    std::cout << "ModID: " << modID << std::endl;
    // std::cout << "Layer: " << layer << std::endl;
    std::cout << "APVmap: " << apvmap << std::endl;
    std::cout << "Position: (" << position.X() << ", " << position.Y() << ", " << position.Z() << ")" << std::endl;
    // std::cout << "Angle: (" << angle.X() << ", " << angle.Y() << ","
    // << angle.Z() << ")" << std::endl;
    std::cout << "Size: (" << size.X() << ", " << size.Y() << ", " << size.Z() << ")" << std::endl;
    std::cout << "Uangle: " << uangle << std::endl;
    std::cout << "Vangle: " << vangle << std::endl;
    std::cout << "Uoffset: " << uoffset << std::endl;
    std::cout << "Voffset: " << voffset << std::endl;
    std::cout << "NstripsU: " << nstripsu << std::endl;
    std::cout << "NstripsV: " << nstripsv << std::endl;
    std::cout << "UAPVs: " << NuAPVs << std::endl;
    std::cout << "VAPVs: " << NvAPVs << std::endl;
  }

  };

  //first term is module ID
  std::map<int, gemInfo> gemInfoMap;



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

  int strip; //if available

  // Define operator< for sorting in std::set
  bool operator<(const apvInfo& other) const {
    if (gemid != other.gemid) return gemid < other.gemid;
    if (axis != other.axis) return axis < other.axis;
    return pos < other.pos;
  }

  void print() const {
    std::cout << "GEMId: " << gemid << std::endl;
    std::cout << "Axis: " << axis << std::endl;
    std::cout << "Pos: " << pos << std::endl;
    std::cout << "VTPcrate: " << vtpcrate << std::endl;
    std::cout << "fiber: " << fiber << std::endl;
    std::cout << "adc_ch: " << adc_ch << std::endl;
  };

};

// Int_t effChan = it->mpd_id << 4 | it->adc_id; //left-shift mpd id by 4 bits and take the bitwise OR with ADC_id to uniquely identify the APV card.

// std::vector<apvIDmap> apvIDs; 
// std::set<apvInfo> apvInfoMap; 

//effChan is first val
// std::map<Int_t, apvInfo> apvInfoMap; 

struct apvInfoKeys {
  int gemid; int axis; int pos;
};

struct apvInfoVals{
  //fill these values they id the APV
  int vtpcrate;
  int fiber; //NOTE: == mpd_id
  int adc_ch;
};

std::map<apvInfoKeys, apvInfoVals> apvInfoMap; 



struct stripInfo{
  int stripNumber;
  int axis; //0=U, 1=V
  int APV_ch, mpd_id, adcID;
  int index; //line in mapFile

  bool operator<(const stripInfo& other) const {
    if (stripNumber != other.stripNumber) return stripNumber < other.stripNumber;
    if (axis != other.axis) return axis < other.axis;
  }

  void print() const {
    std::cout << "StripNumber: " << stripNumber << std::endl;
    std::cout << "Axis: " << axis << std::endl;
    std::cout << "APV_ch: " << APV_ch << std::endl;
    std::cout << "MPD ID: " << mpd_id << std::endl;
    std::cout << "ADC ID: " << adcID << std::endl;
  };
};

std::set<stripInfo> stripInfoMap;

  // SBSGEM::APVmap_t fAPVmapping; //choose APV channel --> strip mapping; there are only three possible values supported for now (see SBSGEM::APVmap_t)

    //Decode map information: 
    // std::vector<sbsgemroi_t>    fGEMroi; //this may need to be modified
    std::vector<mpdmap_t>    fGEMroi; //this may need to be modified
    std::vector<Int_t>       fChanMapData;

  std::array<std::vector<UInt_t>, 4 > APVMAP;

  //some convenience maps: are these actually used yet? 
  std::map<Int_t, Int_t> fAPVch_by_Ustrip;
  std::map<Int_t, Int_t> fAPVch_by_Vstrip;
  std::map<Int_t, Int_t> fMPDID_by_Ustrip;
  std::map<Int_t, Int_t> fMPDID_by_Vstrip;
  std::map<Int_t, Int_t> fADCch_by_Ustrip;
  std::map<Int_t, Int_t> fADCch_by_Vstrip;

  //Constant, module-specific parameters:
  UShort_t fModule; // Module index within a tracker. Should be unique! Since this is a GEM module class, this parameter should be unchanging
  UShort_t fLayer;  // Layer index of this module within a tracker. Since this is a GEM module class, this parameter should be unchanging

  UInt_t fNstripsU; // Total number of strips in this module along the generic "U" axis
  UInt_t fNstripsV; // Total number of strips in this module along the generic "V" axis

  //To be determined from channel map/strip count information:
  UInt_t fNAPVs_U; //Number of APV cards per module along "U" strip direction; this is typically 8, 10, or 12, but could be larger for U/V GEMs
  UInt_t fNAPVs_V; //Number of APV cards per module along "V" strip direction; 
  //UInt_t fNTimeSamplesADC; //Number of ADC time samples (this could be variable in principle, but should be the same for all strips within a module within a run) redundant with fN_MPD_TIME_SAMP

  std::vector<UInt_t> fStrip;  //Strip index of hit (these could be "U" or "V" generalized X and Y), assumed to run from 0..N-1
  std::vector<SBSGEM::GEMaxis_t>  fAxis;  //We just made our enumerated type that has two possible values, makes the code more readable (maybe)


  UInt_t fChan_TimeStamp_low;
  UInt_t fChan_TimeStamp_high;
  UInt_t fChan_MPD_EventCount;




Int_t extraReadDB(const TDatime& date ){
  //Based on SBSGEMModule::ReadDatabase
  std::cout << "[APVFinder::ReadDB]" << std::endl;

  // hitmap_t thisHit;

  // mpdmap_t thisHit; //from SBSGEMModule


  Int_t status;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  const DBRequest request[] =  {
    { "chanmap",        &fChanMapData,        kIntV, 0, 0, 0}, // mandatory: decode map info
    { "apvmap",         &fAPVmapping,    kUInt, 0, 1, 1}, //optional, allow search up the tree if all modules in a setup have the same APV mapping;

    { "layer",          &fLayer,         kUShort, 0, 0, 0}, // mandatory: logical tracking layer must be specified for every module:
    { "nstripsu",       &fNstripsU,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along U axis
    { "nstripsv",       &fNstripsV,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along V axis
    { "uangle",         &fUAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "U" strips wrt X axis
    { "vangle",         &fVAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "V" strips wrt X axis
    { "uoffset",        &fUStripOffset, kDouble, 0, 1, 1}, //optional: position of first U strip
    { "voffset",        &fVStripOffset, kDouble, 0, 1, 1}, //optional: position of first V strip
    { "upitch",         &fUStripPitch,  kDouble, 0, 0, 1}, //mandatory: Pitch of U strips
    { "vpitch",         &fVStripPitch,  kDouble, 0, 0, 1}, //mandatory: Pitch of V strips
  };

  status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree

  // SetProjOps(UVangs);

  fAPVch_by_Ustrip.clear();
  fAPVch_by_Vstrip.clear();
  fMPDID_by_Ustrip.clear();
  fMPDID_by_Vstrip.clear();
  fADCch_by_Ustrip.clear();
  fADCch_by_Vstrip.clear();
  

  Int_t nentry = fChanMapData.size()/fMPDMAP_ROW_SIZE;

  for( Int_t mapline = 0; mapline < nentry; mapline++ ){

    apvInfo currAPVinfo;
    stripInfo currStripInfo;


    currAPVinfo.vtpcrate  = fChanMapData[0+mapline*fMPDMAP_ROW_SIZE];

    // thisdata.slot   = fChanMapData[1+mapline*fMPDMAP_ROW_SIZE];//
    int currSlot = fChanMapData[1+mapline*fMPDMAP_ROW_SIZE];//


    //NOTE: fiber <=>mpd
    currAPVinfo.fiber = fChanMapData[2+mapline*fMPDMAP_ROW_SIZE];
    currStripInfo.mpd_id = fChanMapData[2+mapline*fMPDMAP_ROW_SIZE];
    
    currAPVinfo.gemid = fChanMapData[3+mapline*fMPDMAP_ROW_SIZE];
    
    currAPVinfo.adc_ch = fChanMapData[4+mapline*fMPDMAP_ROW_SIZE];

    currStripInfo.adcID = currAPVinfo.adc_ch;
    // thisdata.i2c    = fChanMapData[5+mapline*fMPDMAP_ROW_SIZE];
    currAPVinfo.pos    = fChanMapData[6+mapline*fMPDMAP_ROW_SIZE];
    // currAPVid.invert = fChanMapData[7+mapline*fMPDMAP_ROW_SIZE];
     int invert = fChanMapData[7+mapline*fMPDMAP_ROW_SIZE];

    currAPVinfo.axis   = fChanMapData[8+mapline*fMPDMAP_ROW_SIZE];

    currStripInfo.axis = currAPVinfo.axis;

    //TODO:?????????????????? whats this
    currStripInfo.index  = mapline;  


  //Populate relevant quantities mapped by strip index:
  for( int ich=0; ich<fN_APV25_CHAN; ich++ ){
    currStripInfo.stripNumber = GetStripNumber( ich, currAPVinfo.pos, invert );
  //   if( stripInfo.axis == 0 ){
  // fAPVch_by_Ustrip[strip] = ich;
  // fMPDID_by_Ustrip[strip] = currStripInfo.mpd_id;
  // fADCch_by_Ustrip[strip] = currStripInfo.adcID;
  //     } else {
  // fAPVch_by_Vstrip[strip] = ich;
  // fMPDID_by_Vstrip[strip] = currStripInfo.mpd_id;
  // fADCch_by_Vstrip[strip] = currStripInfo.adcID;
  //     }

    stripInfoMap.insert(currStripInfo);

    if (gemInfoMap.find(currAPVinfo.gemid) == gemInfoMap.end()){

      //Geometry info is required to be present in the database for each module:
      Int_t err = ReadGeometry( file, date, true );
      if( err ) {
        fclose(file);
        return err;
      }

        gemInfo currGeminfo;
        
        currGeminfo.modID=currAPVinfo.gemid;
        currGeminfo.layer = fLayer;
        currGeminfo.apvmap = 2; // I think all are UVa or XW so 2
        currGeminfo.position = TVector3(fXax,fYax,fZax);
        currGeminfo.angle = ;
        currGeminfo.size = fSize;
        currGeminfo.uangle=fUAngle;
        currGeminfo.vangle=fVAngle;
        currGeminfo.uoffset = fUStripOffset;
        currGeminfo.voffset = fVStripOffset;

        // nUVstrips=GetUVstrips(currGeminfo.modID)
        
        currGeminfo.nstripsu = fNstripsU;
        currGeminfo.nstripsv = fNstripsV;

        currGeminfo.NuAPVs = fNstripsU/fN_APV25_CHAN;
        currGeminfo.NvAPVs = fNstripsV/fN_APV25_CHAN;

        
        gemInfoMap[currAPVinfo.gemid]=currGeminfo; //
        
      }
    }
    //fiber=mpd_id
    Int_t effChan = currAPVinfo.fiber << 4 | currAPVinfo.adc_ch; //left-shift mpd id by 4 bits and take the bitwise OR with ADC_id to uniquely identify the APV card.

    // if apvInfoMap[effChan]

    // apvInfoMap.insert(currAPVinfo);
    // apvInfoMap[effChan] = currAPVinfo; 

    apvInfoKeys currAPVkeys;
    currAPVkeys.gemid = currAPVinfo.gemid;
    currAPVkeys.axis = currAPVinfo.axis;
    currAPVkeys.pos = currAPVinfo.pos;

    if (apvInfoMap.find(currAPVkeys)==apvInfoMap.end()){
      apvInfoVals currAPVvals;
      currAPVvals.vtpcrate = currAPVinfo.vtpcrate;
      currAPVvals.fiber = currAPVinfo.fiber;
      currAPVvals.adc_ch = currAPVinfo.adc_ch;
    }
  };

   //resize all the "decoded strip" arrays to their maximum possible values for this module:
  // UInt_t nstripsmax = fNstripsU + fNstripsV;
  
  // fStrip.resize( nstripsmax );
  // fAxis.resize( nstripsmax );
  
  }
}



Int_t ReadDB(const TDatime& date ){
  //Based on SBSGEMModule::ReadDatabase
  std::cout << "[APVFinder::ReadDB]" << std::endl;

  // hitmap_t thisHit;

  // mpdmap_t thisHit; //from SBSGEMModule


  Int_t status;

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  const DBRequest request[] =  {
    { "chanmap",        &fChanMapData,        kIntV, 0, 0, 0}, // mandatory: decode map info
    { "apvmap",         &fAPVmapping,    kUInt, 0, 1, 1}, //optional, allow search up the tree if all modules in a setup have the same APV mapping;

    { "layer",          &fLayer,         kUShort, 0, 0, 0}, // mandatory: logical tracking layer must be specified for every module:
    { "nstripsu",       &fNstripsU,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along U axis
    { "nstripsv",       &fNstripsV,     kUInt, 0, 0, 1}, //mandatory: number of strips in module along V axis
    { "uangle",         &fUAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "U" strips wrt X axis
    { "vangle",         &fVAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "V" strips wrt X axis
    { "uoffset",        &fUStripOffset, kDouble, 0, 1, 1}, //optional: position of first U strip
    { "voffset",        &fVStripOffset, kDouble, 0, 1, 1}, //optional: position of first V strip
    { "upitch",         &fUStripPitch,  kDouble, 0, 0, 1}, //mandatory: Pitch of U strips
    { "vpitch",         &fVStripPitch,  kDouble, 0, 0, 1}, //mandatory: Pitch of V strips
  };

  status = LoadDB( file, date, request, fPrefix, 1 ); //The "1" after fPrefix means search up the tree

  // SetProjOps(UVangs);

  fAPVch_by_Ustrip.clear();
  fAPVch_by_Vstrip.clear();
  fMPDID_by_Ustrip.clear();
  fMPDID_by_Vstrip.clear();
  fADCch_by_Ustrip.clear();
  fADCch_by_Vstrip.clear();
  

  Int_t nentry = fChanMapData.size()/fMPDMAP_ROW_SIZE;

  for( Int_t mapline = 0; mapline < nentry; mapline++ ){

    apvInfo currAPVinfo;
    stripInfo currStripInfo;


    currAPVinfo.vtpcrate  = fChanMapData[0+mapline*fMPDMAP_ROW_SIZE];

    // thisdata.slot   = fChanMapData[1+mapline*fMPDMAP_ROW_SIZE];//
    int currSlot = fChanMapData[1+mapline*fMPDMAP_ROW_SIZE];//


    //NOTE: fiber <=>mpd
    currAPVinfo.fiber = fChanMapData[2+mapline*fMPDMAP_ROW_SIZE];
    currStripInfo.mpd_id = fChanMapData[2+mapline*fMPDMAP_ROW_SIZE];
    
    currAPVinfo.gemid = fChanMapData[3+mapline*fMPDMAP_ROW_SIZE];
    
    currAPVinfo.adc_ch = fChanMapData[4+mapline*fMPDMAP_ROW_SIZE];

    currStripInfo.adcID = currAPVinfo.adc_ch;
    // thisdata.i2c    = fChanMapData[5+mapline*fMPDMAP_ROW_SIZE];
    currAPVinfo.pos    = fChanMapData[6+mapline*fMPDMAP_ROW_SIZE];
    // currAPVid.invert = fChanMapData[7+mapline*fMPDMAP_ROW_SIZE];
     int invert = fChanMapData[7+mapline*fMPDMAP_ROW_SIZE];

    currAPVinfo.axis   = fChanMapData[8+mapline*fMPDMAP_ROW_SIZE];

    currStripInfo.axis = currAPVinfo.axis;

    //TODO:?????????????????? whats this
    currStripInfo.index  = mapline;  


  //Populate relevant quantities mapped by strip index:
  for( int ich=0; ich<fN_APV25_CHAN; ich++ ){
    currStripInfo.stripNumber = GetStripNumber( ich, currAPVinfo.pos, invert );
  //   if( stripInfo.axis == 0 ){
  // fAPVch_by_Ustrip[strip] = ich;
  // fMPDID_by_Ustrip[strip] = currStripInfo.mpd_id;
  // fADCch_by_Ustrip[strip] = currStripInfo.adcID;
  //     } else {
  // fAPVch_by_Vstrip[strip] = ich;
  // fMPDID_by_Vstrip[strip] = currStripInfo.mpd_id;
  // fADCch_by_Vstrip[strip] = currStripInfo.adcID;
  //     }

    stripInfoMap.insert(currStripInfo);

    if (gemInfoMap.find(currAPVinfo.gemid) == gemInfoMap.end()){

      //Geometry info is required to be present in the database for each module:
      Int_t err = ReadGeometry( file, date, true );
      if( err ) {
        fclose(file);
        return err;
      }

        gemInfo currGeminfo;
        
        currGeminfo.modID=currAPVinfo.gemid;
        currGeminfo.layer = fLayer;
        currGeminfo.apvmap = 2; // I think all are UVa or XW so 2
        currGeminfo.position = TVector3(fXax,fYax,fZax);
        currGeminfo.angle = ;
        currGeminfo.size = fSize;
        currGeminfo.uangle=fUAngle;
        currGeminfo.vangle=fVAngle;
        currGeminfo.uoffset = fUStripOffset;
        currGeminfo.voffset = fVStripOffset;

        // nUVstrips=GetUVstrips(currGeminfo.modID)
        
        currGeminfo.nstripsu = fNstripsU;
        currGeminfo.nstripsv = fNstripsV;

        currGeminfo.NuAPVs = fNstripsU/fN_APV25_CHAN;
        currGeminfo.NvAPVs = fNstripsV/fN_APV25_CHAN;

        
        gemInfoMap[currAPVinfo.gemid]=currGeminfo; //
        
      }
    }
    //fiber=mpd_id
    Int_t effChan = currAPVinfo.fiber << 4 | currAPVinfo.adc_ch; //left-shift mpd id by 4 bits and take the bitwise OR with ADC_id to uniquely identify the APV card.

    // if apvInfoMap[effChan]

    // apvInfoMap.insert(currAPVinfo);
    // apvInfoMap[effChan] = currAPVinfo; 

    apvInfoKeys currAPVkeys;
    currAPVkeys.gemid = currAPVinfo.gemid;
    currAPVkeys.axis = currAPVinfo.axis;
    currAPVkeys.pos = currAPVinfo.pos;

    if (apvInfoMap.find(currAPVkeys)==apvInfoMap.end()){
      apvInfoVals currAPVvals;
      currAPVvals.vtpcrate = currAPVinfo.vtpcrate;
      currAPVvals.fiber = currAPVinfo.fiber;
      currAPVvals.adc_ch = currAPVinfo.adc_ch;
    }
  };

   //resize all the "decoded strip" arrays to their maximum possible values for this module:
  // UInt_t nstripsmax = fNstripsU + fNstripsV;
  
  // fStrip.resize( nstripsmax );
  // fAxis.resize( nstripsmax );
  
  }
}


//GEOMETRICAL PARAMETERS:
Double_t fUStripPitch;    //strip pitch along U, will virtually always be 0.4 mm
Double_t fVStripPitch;    //strip pitch along V, will virtually always be 0.4 mm
Double_t fUStripOffset;   //position of first U strip along the direction it measures:
Double_t fVStripOffset;   //position of first V sttrip alogn the direction it measures:
Double_t fUAngle;         //Angle between U strips and "X" axis of TRANSPORT coordinates;
Double_t fVAngle;         //Angle between V strips and "X" axis of TRANSPORT coordinates;
Double_t fPxU;            //U Strip X projection = cos( UAngle );
Double_t fPyU;            //U Strip Y projection = sin( UAngle );
Double_t fPxV;            //V Strip X projection = cos( VAngle );
Double_t fPyV;            //V Strip Y projection = sin( VAngle );


  //Aligned geometry
  // TVector3 alignPos; 
  // TVector3 alignAngs; 

  // Geometry from THaDetectorBase
  TVector3      fOrigin;    // Position of detector (m)
  Double_t      fSize[3];   // Detector size in x,y,z (m) - x,y are half-widths

  TVector3      fXax;       // X axis of the detector plane
  TVector3      fYax;       // Y axis of the detector plane
  TVector3      fZax;       // Normal to the detector plane



Int_t ReadGeometry( FILE *file, const TDatime &date, Bool_t required ){ //We start with a copy of THaDetectorBase::ReadGeometry and modify accordingly:
  // Read this detector's basic geometry information from the database.
  // Derived classes may override to read more advanced data.

  const char* const here = "ReadGeometry";

  vector<double> position, size, angles;
  Bool_t optional = !required;
  DBRequest request[] = {
    { "position", &position, kDoubleV, 0, optional, 0,
      "\"position\" (detector position [m])" },
    { "size",     &size,     kDoubleV, 0, optional, 1,
      "\"size\" (detector size [m])" },
    { "angle",    &angles,   kDoubleV, 0, true, 0,
      "\"angle\" (detector angles(s) [deg]" },
    { nullptr }
  };
  Int_t err = LoadDB( file, date, request );
  if( err )
    return kInitError;

    // const DBRequest alignRequest = DBRequest{

    //   { "angle",         &fUAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "U" strips wrt X axis
  
    //   { "angle",         &fVAngle,       kDouble, 0, 0, 1}, //mandatory: Angle of "V" strips wrt X axis
    // }
  
    // aligned = LoadDB( file, date, &alignRequest, fPrefix, 1 ); //The "1" after fPrefix means search up the tree

    if( !position.empty() ) {
      if( position.size() != 3 ) {
        Error( Here(here), "Incorrect number of values = %u for "
         "detector position. Must be exactly 3. Fix database.",
         static_cast<unsigned int>(position.size()) );
        return 1;
      }
      fOrigin.SetXYZ( position[0], position[1], position[2] );
    }
    else
      fOrigin.SetXYZ(0,0,0);
  
    if( !size.empty() ) {
      if( size.size() != 3 ) {
        Error( Here(here), "Incorrect number of values = %u for "
         "detector size. Must be exactly 3. Fix database.",
         static_cast<unsigned int>(size.size()) );
        return 2;
      }
      if( size[0] == 0 || size[1] == 0 || size[2] == 0 ) {
        Error( Here(here), "Illegal zero detector dimension. Fix database." );
        return 3;
      }
      if( size[0] < 0 || size[1] < 0 || size[2] < 0 ) {
        Warning( Here(here), "Illegal negative value for detector dimension. "
           "Taking absolute. Check database." );
      }
      fSize[0] = 0.5 * TMath::Abs(size[0]);
      fSize[1] = 0.5 * TMath::Abs(size[1]);
      fSize[2] = TMath::Abs(size[2]);
    }
    else
      fSize[0] = fSize[1] = fSize[2] = kBig;

      if( !angles.empty() ) {
        if( angles.size() != 1 && angles.size() != 3 ) {
          Error( Here(here), "Incorrect number of values = %u for "
           "detector angle(s). Must be either 1 or 3. Fix database.",
           static_cast<unsigned int>(angles.size()) );
          return 4;
        }
        // If one angle is given, it indicates a rotation about y, as before.
        // If three angles are given, they are interpreted as rotations about the X, Y, and Z axes, respectively:
        // 
        if( angles.size() == 1 ) {
          DefineAxes( angles[0] * TMath::DegToRad() );
        }
        else {
          TRotation RotTemp;
    
          // So let's review how to define the detector axes correctly.
    
          // THaDetectorBase::DetToTrackCoord(TVector3 p) returns returns p.X * fXax() + p.Y * fYax() + p.Z * fZax() + fOrigin
          // In the standalone code, we do Rot * (p) + fOrigin (essentially):
          // So in matrix form, when we do TRotation::RotateX(alpha), we get:
    
          // RotTemp * Point =  |  1    0            0         |    |  p.X()  |
          //                    |  0   cos(alpha) -sin(alpha)  | *  |  p.Y()  |
          //                    |  0   sin(alpha)  cos(alpha)  |    |  p.Z()  |
          // 
          // This definition ***appears**** to be consistent with the "sense" of the rotation as applied by the standalone code.
          // The detector axes are defined as the columns of the rotation matrix. We will have to test that it is working correctly, however:
          
          RotTemp.RotateX( angles[0] * TMath::DegToRad() );
          RotTemp.RotateY( angles[1] * TMath::DegToRad() );
          RotTemp.RotateZ( angles[2] * TMath::DegToRad() );
          
          fXax.SetXYZ( RotTemp.XX(), RotTemp.YX(), RotTemp.ZX() );
          fYax.SetXYZ( RotTemp.XY(), RotTemp.YY(), RotTemp.ZY() );
          fZax.SetXYZ( RotTemp.XZ(), RotTemp.YZ(), RotTemp.ZZ() );
        }
      } else
        DefineAxes(0);


      // gemInfo currGemInfo
    
      return 0;
  
}

// from THaDetectorBase::
void DefineAxes( Double_t rotation_angle )
{
  // Define detector orientation, assuming a tilt by rotation_angle (in rad)
  // around the y-axis

  fXax.SetXYZ( 1.0, 0.0, 0.0 );
  fXax.RotateY( rotation_angle );
  fYax.SetXYZ( 0.0, 1.0, 0.0 );
  fZax = fXax.Cross(fYax);

}


// UInt_t numLines = 0;

// void LoadLine(UInt_t currLine){
// void LoadLine(std::string LineStr){
// std::vector<HitInfo> LoadLine(std::string &LineStr){
XY_ROI LoadLine(std::string &LineStr){

  std::stringstream lineStream(LineStr);

  XY_ROI hit;
  lineStream >> hit.ECalBin >> hit.GEMLayer >> hit.xMin >> hit.xMax >> hit.yMin >> hit.yMax;

  if (std::isnan(hit.xMin) || std::isnan(hit.xMax) || std::isnan(hit.yMin) || std::isnan(hit.yMax)) {
    std::cerr << "Skipping line with NaN values: " << LineStr << std::endl;
    hit.lineNumber = -1;  // Mark as invalid (or handle it as needed)
  }

  return hit;
}



  // struct hitmap_t {
  //   UInt_t crate;
  //   UInt_t slot;
  //   UInt_t mpd_id;
  //   UInt_t gem_id;
  //   UInt_t adc_id;
  //   UInt_t i2c;
  //   UInt_t pos;
  //   UInt_t invert;
  //   UInt_t axis; //needed to add axis to the decode map
  //   UInt_t index;
  // };

  apvInfoVals GetAPV(int gemid, int axis, int pos){ 
    apvInfoKeys thisKey;
    thisKey.gem_id = gemid;
    thisKey.axis = axis;
    thisKey.pos = pos;

    return apvInfoMap[thisKey]
  }


  
  TVector3 GetModDimensions(int apvmap){
    if (apvmap==2){
      //in meters
      return TVector3(1.5, .4, .001);
    }
    if (apvmap==1){
      return TVector3(.512, .6144, .001);
    }
  }

  // TVector3 GetModPostion(int modID){}

  //TODO: maybe make func to pull position of each from reference files in case of changes in measurements
  //see db_sbs.gemFT.dat
  //Positions in internal GEM coordinate system.
  std::map<int, TVector3> ModPositionmap{
    {0, TVector3(0., 0., 0.)},
    {1, TVector3(0., 0., 0.075)},
    {2, TVector3(0., 0., 0.202)},
    {3, TVector3(0., 0., 0.341)},
    {4, TVector3(0., 0., 0.484)},
    {5, TVector3(0., 0., 0.624)},

    {6, TVector3(-.766, 0., 0.769)},
    // {6, TVector3(-.766, 0., 0.737)},
    {7, TVector3(-.256, 0., 0.737)},
    // {7, TVector3(-.256, 0., 0.769)},
    {8, TVector3(.265, 0., 0.769)},
    // {8, TVector3(.265, 0., 0.737)},
    {9, TVector3(.765, 0., 0.737)},
    // {9, TVector3(.765, 0., 0.769)},
    
    {10, TVector3(-0.766, 0, 0.897)},
    // {10, TVector3(-0.766, 0, 0.865)},
    {11, TVector3(-0.256, 0, 0.865)},
    // {11, TVector3(-0.256, 0, 0.897)},
    {12, TVector3(0.256, 0, 0.897)},
    // {12, TVector3(0.256, 0, 0.865)},
    {13, TVector3(0.765, 0, 0.865)}
    // {13, TVector3(0.765, 0, 0.897)}
  };
  

  TVector2 GetUVang(int modID){
    if(modID ==0){ return TVector2(180., 135.);}
    else if(modID ==1){ return TVector2(180., -135.);}
    else if (modID >=2 && modID < 5){return TVector2(150., -150.);}
    else if (modID >=6 && modID < 14){return TVector2(0., -90.);}
  }
  

  TVector2 LayerUVang(int layer){
    if(layer ==0){ return TVector2(180., 135.);}
    else if(layer ==1){ return TVector2(180., -135.);}
    else if (layer >=2 && layer <6){return TVector2(150., -150.);}
    else if (layer>=6){return TVector2(0., -90.);}
  }

  //from SBSGEMModule
  double fPxU;
  double fPyU;
  double fPxV;
  double fPyV;

  //Get Projection operators given angles u and v
  void SetProjOps(TVector2 UVangles){
    auto fUAngle=UVangles.X(); auto fVAngle=UVangles.X();

    fPxU = cos( fUAngle * TMath::DegToRad());
    fPyU = sin( fUAngle * TMath::DegToRad());
    fPxV = cos( fVAngle * TMath::DegToRad());
    fPyV = sin( fVAngle * TMath::DegToRad());
  }

  //.0004m = 400 micrometers
  double GetPitch(){return .0004;}

  double GetUVoffset(int modID, char *axis){
    if(modID ==0 || modID == 1){
      if (axis == "u"){ return .0176;}
      else if (axis == "v"){return 0.;}
    }
    else if (modID >=2 && modID < 5){return .0108;}
    else if (modID >=6 && modID < 14){return 0.;}
  }
  
  TVector2 GetOffset(int modID){
    if(modID ==0 || modID == 1){ return TVector2(.0176, 0.);}
    else if (modID >=2 && modID < 5){return TVector2(.0108, .0108);}
    else if (modID >=6 && modID < 14){return TVector2(0., 0.);}
  }
 
  int GetNstrips(int modID, char *axis){
    if(modID ==0 || modID == 1){
      if (axis == "u"){ return 3968;}
      else if (axis == "v"){return 3456;}

    else if (modID >=2 && modID < 5){return 3840;}
    
    else if (modID >=6 && modID < 14){
      if (axis == "u"){ return 1280;}
      else if (axis == "v"){return 1536;}
    }
    }
  }

  Int_t SBSGEMModule::GetStripNumber( UInt_t rawstrip, UInt_t pos, UInt_t invert ){
    Int_t RstripNb = APVMAP[fAPVmapping][rawstrip];
    RstripNb = RstripNb + (127-2*RstripNb)*invert;
    Int_t RstripPos = RstripNb + 128*pos;
  
    if( fIsMC ){
      return rawstrip + 128*pos;
    }
    
    return RstripPos;
  }
  
  TVector2 GetNstrips(int modID){
    if(modID ==0 || modID == 1){ 
      return TVector2(3968, 3456);}

    else if (modID >=2 && modID < 6){
      return TVector2(3840, 3840);}
    
    else if (modID >=6 && modID < 14){
      return TVector2(1280, 1536);}
  }

    //Get min and max of either u or v and returv the corresponding strip IDs
    TVector2 GetStripRange(){}

    //  int stripToAPVid(int modID, )
  
    TVector2 GetNAPVs(int modID){
      
      //128 strips per APV
      return GetNstrips(modID)/128;
      
      // NOTE: if(modID ==0 || modID == 1){ 
      //   return TVector2(31, 27);}
  
      // NOTE: else if (modID >=2 && modID < 5){
      //   return TVector2(30, 30);}
      
      // NOTE: else if (modID >=6 && modID < 14){
      //   return TVector2(10, 12);}
    }
    
    //TODO: implement invert later by total num of strips - LocalID
    //stripIDs in module along U,V counted from top to bottom
    TVector2 GetLocalStripIDs(int modID, TVector2 UVhit){ 
      
    //   //TODO: finish subtracting relative to TOP of module

    //   //TODO: THINK ABOUT THIS
    //   TVector2 LocalHitPos = 
    //   ModPositionmap[modID].XYvector() - GetModDimensions(modID).XYvector() // set pos to top of module
    //   + UVhit 
    //   ;


    //   TVector2 LocalStripIDs = TVector2(modID

    //   return GetLocalStripIDs;
    }

    //Get physical positions of max and min strips of the APV relative to top of layer
    TVector2 GetAPVposRange(int modID, TVector2){}

    
    int UorVtoAPVid(int modID, TVector2 localUVPos,
      int axis //U/X=0, V/Y=1
    ){
      double localPos;
      
      double APVsize = 128*GetPitch();

      int numAPVs;

      if (axis == 0){
        numAPVs = int(GetNAPVs(modID).X());
        localPos=localUVPos.X();}
      else if (axis == 1){
        numAPVs = int(GetNAPVs(modID).Y());
        localPos=localUVPos.Y();} 
    
      //numAPVs along this axis

      // TVector2 APVedges

      int layerAPVid;

      //NOTE: PER Layer APV id;
      if(modID >= 0 || modID < 6){
        for( int i = 0; i <= numAPVs; i++){
          // APVedges= TVector2(i, (i+1))*APVsize;
          if (localPos>=i*APVsize && 
            localPos<(i+1)*APVsize)
            { layerAPVid = i; break;}
         }
      }
      else if(modID >= 6 || modID < 10){
        modID -= 6;
        for( int i = 0; i <= numAPVs; i++)
          {if (localPos>=i*APVsize && 
            localPos<(i+1)*APVsize)
            { layerAPVid = i; break;}}
      }
      else if (modID >= 10 || modID < 14){
        modID-=10;
        for( int i = 0; i <= numAPVs; i++)
          {if (localPos>=i*APVsize && 
            localPos<(i+1)*APVsize)
            { layerAPVid = i; break;}}
       }
      

      int currMod = 0;
      int globalAPVid=layerAPVid;

      while(currMod < modID) {
        int currNAPVs;
        if (axis == 0){
          currNAPVs = int(GetNAPVs(currMod).X());
          localPos=localUVPos.X();}
        else if (axis == 1){
          currNAPVs = int(GetNAPVs(currMod).Y());
          localPos=localUVPos.Y();} 

        globalAPVid += currNAPVs;
        currMod++;
      }
      return globalAPVid;
      
    }

  int GetAPVid(int vtp, int fiber, int adc_channel){


  }



  // TODO: for now just copied from SBSGEMModule 
  TVector2 XYtoUV( TVector2 XY ){

    double Xtemp = XY.X();
    double Ytemp = XY.Y();
  
    double Utemp = Xtemp*fPxU + Ytemp*fPyU;
    double Vtemp = Xtemp*fPxV + Ytemp*fPyV;
  
    return TVector2(Utemp,Vtemp);
  }


  int GetModId(int layerID, TVector2 hitPos, TVector3 modDims){
    int modID = -1;
    if(layerID>5){
      for (int i=6; i<14; i++){
        TVector3 modPos =ModPositionmap[i];
        // modSize=modDims[0];
        
        //TODO: cconsider cases ON the bounds
        if (hitPos.X()> modPos[0]-modDims[0]
          &&
          hitPos.X() < modPos[0]+modDims[0])
          {modID = i; break;}
          
        }
      if (modID==-1){
        std::cerr <<"Error: Cant find modID"<<std::endl;
        exit(-1);}
      }
      return modID;
  }

  // TODO: use map file to associate adc etc directly with pos

  UV_ROI XYtoUV_ROI(XY_ROI hitXY){
 
    TVector2 minXY = TVector2(hitXY.xMin, hitXY.yMin);
    TVector2 maxXY = TVector2(hitXY.xMax, hitXY.yMax);

    UV_ROI thisUV_ROI;

    int layerID = hitXY.GEMLayer; thisUV_ROI.GEMLayer = layerID;
    
    TVector2 hitUV;
    
    // int modID;  
    
    TVector3 modDims;


    int apvmap; //type of APV config/module type
    // if apvmap= 0; //INFN
    // if apvmap= 1; //UVA XY
    // if apvmap= 2; //UVA U/V(layers 2-5) or X/W(layers 0-1)
  
    //NOTE: Layers 0 and 1 are UVA X/W
    // if (layerID in range(0,6)){apvmap=2;}
    if (layerID >=0 && layerID<6 )
    {apvmap=2; modDims=TVector3(1.5, .4, .001);}
    else{apvmap=1; modDims=TVector3(.512, .6144, .001);};


    thisUV_ROI.MinModID = GetModId(layerID, minXY, modDims);
    thisUV_ROI.MaxModID = GetModId(layerID, maxXY, modDims);

    // auto localX=hitPos[0]-mod[modID].position;
    //relative to center of module

    TVector2 currUVangs=GetUVang(layerID);

    SetProjOps(currUVangs);

    TVector2 minUV=XYtoUV(minXY); TVector2 maxUV=XYtoUV(maxXY);

    thisUV_ROI.Umin = minUV.X(); thisUV_ROI.Vmin = minUV.Y();
    thisUV_ROI.Umax = maxUV.X(); thisUV_ROI.Vmax = maxUV.Y();
    // TVector2 thisUV_ROI=XYtoUV(hitPos);  

    //NOTE: UV_ROI should now have
    // layerid, modid, lineNumber(from earlier in FindAPV(), 
    // Umin, Vmin, Umax, Vmax
    return thisUV_ROI;
  }




  // UInt_t UV_ROItoAPVid(int modID, TVector2 UVhit){
  //Stripwise APV finding
  UInt_t UV_ROItoAPVs(int modID, TVector2 UVhit){
    // RecodeAPV(adc)

    // if (modID == 0 || modID == 1){        }


    //TODO: either do by finding strip ID then which APV it relates to OR make a function to find the physical position range of APV 

    //TODO: IMPLEMENT THE APV position range

    // TVector2 modUVstripCount=GetNstrips(modID);

    // int stripsPerAPV=128;

    //   TVector2 LocalStripIds=



    // //NOTE: "Local" per module APVids
    // TVector2 LocalAPVids = LocalStripIds%GetNAPVs(modID);


    

    // //TODO: need to account for inverted?
  
    // UstripID=(UInt_t UrawStrip, hitUV[0], ThisMod.IsInverted()); 
    // VstripID=(UInt_t VrawStrip, hitUV[1], ThisMod.IsInverted()); 
    
    // UInt_t uAPVid, vAPVid;
    // //todo: Again theres a much better way to do this fix later
    // for (i in range(nAPVsU)){
    //   if (APV[i].minStrip <= UstripID || 
    //     UstripID <= APV[i].maxStrip)
    //     {uAPVid = i;break}
    //   else if(i>=nAPVsU)
    //   }
    
    // for (i in range(nAPVsV)){
    //   if (APV[i].minStrip <= VstripID || 
    //     UstripID <= APV[i].maxStrip)
    //     {vAPVid = i; break}
    //   else if(i>=nAPVsV)
    //   }        
      //GemPos measured top to muttom

  }


  // UV_ROI XYtoAPV(UInt_t layer, TVector2 hitPos){
  UV_ROI XYtoAPV(UInt_t layer, TVector2 hitPos){
//   // void XYtoAPV(TVector3 &hitPos){
//   // void XYtoAPV(XY_ROI HitInfo){

//     // TVector2(HitInfo.GEMLayer) hitPos;
//     UV_ROI thisUV_ROI;

//     thisUV_ROI.GEMLayer = layer;

//     TVector2 hitUV;
//     int modID;
  
//     TVector3 modDims;
    
//     // int layerID = GetLayer(hitPos[3]);
//     // int layerID = HitInfo.GEMLayer;
//     int layerID = layer;
    
//     int apvmap; //type of APV config/module type
//     // if apvmap= 0; //INFN
//     // if apvmap= 1; //UVA XY
//     // if apvmap= 2; //UVA U/V(layers 2-5) or X/W(layers 0-1)
//     //NOTE: Z coordinate offset between 2 consecutive GEM modules within a UVA-XY GEM layer is found to be 0.032 meters -> Anu on Dec 2nd 2021 after referring engineering drawings. Andrew has taken this value to be 0.036 meters in his initial guess values for SBS.
  
//     //NOTE: Layers 0 and 1 are UVA X/W
//     // if (layerID in range(0,6)){apvmap=2;}
//     if (layerID >=0 && layerID<6 )
//     {modID=layerID; apvmap=2; 
//       modDims=TVector3(1.5, .4, .001);
//     }
//     else{apvmap=1; modDims=TVector3(.512, .6144, .001);};
  
//     // if (layerID>=6 && layerID<10){}
  
//     //todo: choose different(more efficeint?) type 
//     // double modSize;//x size(up and down)
  
//     TVector2 modXY;
  
//     // if (apvmap == 1){
//     //   // auto modDims = TVector3(0.512, 0.6144, 0.001);
//     //   // modSize = modDims[0];
//     //   modXY = TVector2(0.512, 0.6144);
//     //   // modSize = modDims[0];      
//     // }
//     // else{
//     //   // auto modDims = TVector3(1.5, 0.4, 0.001);}
//     //   modXY = TVector2(1.5, 0.4);}
  
//     //Probably better way to do this
//     //Loop to search for modID given its position and surface area
//     // for (i in range(2, 26)){
  
//     //get module from position relative to center of the layer(remember -x)

//     TVector3 modPos;


//     //TODO: Need to make this arbitrary so can input min and max points seperately
    
//     for (int i=6; i<14; i++){
//       modPos=ModPositionmap[i];
//       // modSize=modDims[0];
      
//       //TODO: cconsider cases ON the bounds
//       if (hitPos[0]> modPos[0]-modDims[0]
//         ||
//         hitPos[0] < modPos[0]+modDims[0])
//         {modID = i; ThisMod=mod[i]; break}
        
//       else if (i == 14){
//         std::cerr <<"Error: Cant find modID"<<std::endl;
//         exit(-1);}
//       }
    
//     // auto localX=hitPos[0]-mod[modID].position;
//     //relative to center of module

//     TVector2 currUVangs=GetUVang(modID);

//     SetProjOps(currUVangs);
  
//     TVector2 minUV = XYtoUV(TVector2(HitInfo.xMin, HitInfo.yMin));
    
//     TVector2 maxUV = XYtoUV(TVector2(HitInfo.xMax, HitInfo.yMax));

  
//     // if (modID == 0 || modID == 1){        }
  
//     UstripID=(UInt_t UrawStrip, hitUV[0], ThisMod.IsInverted()); 
//     VstripID=(UInt_t VrawStrip, hitUV[1], ThisMod.IsInverted()); 
    
//     UInt_t uAPVid, vAPVid;
//     //todo: Again theres a much better way to do this fix later
//     for (i in range(nAPVsU)){
//       if (APV[i].minStrip <= UstripID || 
//         UstripID <= APV[i].maxStrip)
//         {uAPVid = i;break}
//       else if(i>=nAPVsU)
//       }
    
//     for (i in range(nAPVsV)){
//       if (APV[i].minStrip <= VstripID || 
//         UstripID <= APV[i].maxStrip)
//         {vAPVid = i; break}
//       else if(i>=nAPVsV)
//       }        
//       //GemPos measured top to muttom
//     
}




  //"Inverse" of hitpos from SBSGEMModule::find_clusters_1D()
  int hitposToStripID( TVector2 hitpos, int &Nstrips ){
    //NOTE: dont need b/c find APV by spatial extent instead

    //TODO: UV or XY hitpos?

    // int istrip = int(
    // (hitpos-offset)/pitch 
    // + .5*fNstrips - .5);

    // return istrip;
  }



//TODO: if using this method need to find mpd, adc
void Decode(){
  // Int_t effChan = it->mpd_id << 4 | it->adc_id; //left-shift mpd id by 4 bits and take the bitwise OR with ADC_id to uniquely identify the APV card.
  //mpd_id is not necessarily equal to slot, but that seems to be the convention in many cases

}


void Clear( Option_t* opt){ //we will want to clear out many more things too
  // Modify this a little bit so we only clear out the "hit counters", not necessarily the
  // arrays themselves, to make the decoding more efficient:

  THaSubDetector::Clear(opt);
  
  fNstrips_hit = 0;
  fNstrips_hitU = 0;
  fNstrips_hitV = 0;
  fNstrips_hitU_neg = 0;
  fNstrips_hitV_neg = 0;
  fNdecoded_ADCsamples = 0;
  fIsDecoded = false;

  fClustering1DIsDone = false;
  
  fTrackPassedThrough = 0;

  //numbers of strips passing basic zero suppression thresholds and timing cuts:
  fNstrips_keep = 0;
  fNstrips_keepU = 0;
  fNstrips_keepV = 0;
  //numbers of strips passing basic zero suppression thresholds, timing cuts, and higher max. sample and strip sum thresholds for
  // local max:
  fNstrips_keep_lmax = 0;
  fNstrips_keep_lmaxU = 0;
  fNstrips_keep_lmaxV = 0;
  
  
  fNclustU = 0;
  fNclustV = 0;
  fNclustU_pos = 0;
  fNclustV_pos = 0;
  fNclustU_neg = 0;
  fNclustV_neg = 0;
  fNclustU_total = 0;
  fNclustV_total = 0;
  //later we may need to check whether this is a performance bottleneck:
  fUclusters.clear();
  fVclusters.clear();
  fN2Dhits = 0;
  //similar here:
  fHits.clear();

  fTrigTime = 0.0;
  
  fCM_online.assign(fN_MPD_TIME_SAMP,0.0);

  fxcmin.clear();
  fxcmax.clear();
  fycmin.clear();
  fycmax.clear();
  
  //fStripAxis.clear();
  // fADCsamples1D.clear();
  // fStripTrackIndex.clear();
  // fRawADCsamples1D.clear();

  // fStripIsU.clear();
  // fStripIsV.clear();
  // fStripADCavg.clear();
  
  // fUstripIndex.clear();
  // fVstripIndex.clear();
  // fStrip.clear();
  // fAxis.clear();
  // fADCsamples.clear();
  // fRawADCsamples.clear();
  // fADCsums.clear();
  // fKeepStrip.clear();
  // fMaxSamp.clear();
  // fADCmax.clear();
  // fTmean.clear();
  // fTsigma.clear();
  // fTcorr.clear();

  //THaSubDetector::Clear(opt);
}

  TVector2 minXminY;
  TVector2 minXmaxY;
  TVector2 maxXminY;
  TVector2 maxXmaxY;
  
  std::vector<TVector2> roiSquare={minXminY, minXmaxY, maxXminY, maxXmaxY};

}


// ################################################################


// int FindAPV(ifstream aHitFile){

void FindAPV(const TDatime& date, const char *aHitFile){

  using namespace APVFinder;

  std::ifstream hitFile(aHitFile);  

  if (!hitFile.is_open()) {  // Check if the file opened successfully
      std::cerr << "Unable to open file\n";
      // return 1;
  }

  ReadDB(date); 
  
  // DBfile = 

  // makeAPVid(date);

  std::string currLine;

  // std::vector<LocalROIInfo> localHitMap; 

  // skip first line(Column headers)
  std::getline(hitFile, currLine);

  int lineCount = 1;
  
  while (std::getline(hitFile, currLine)) {

    // UV_ROI thisUV_ROI;
    // thisUV_ROI.lineNumber=lineCount;
    
    // std::stringstream lineStream(currLine);
    
    XY_ROI currHit = LoadLine(currLine);
    // HitInfo currHit = LoadLine(currLine);
    
    if (currHit.lineNumber == -1) {continue;}

    //Add the txt file line number just for sake of indexing/debugging
    currHit.lineNumber=lineCount;

    // thisUV_ROI.GEMLayer=currHit.GEMLayer;
    
    
    
    //Take 4 corners of the XY ROI 
    minXY = TVector2(currHit.xMin, currHit.yMin);
    maxXY = TVector2(currHit.xMax, currHit.yMax);
    minXmaxY = TVector2(currHit.xMin, currHit.yMax);
    maxXminY = TVector2(currHit.xMax, currHit.yMin);
    

    for (i=0; i< roiSquare.size(); i++) {
    }


    //Fill all but APVids
    UV_ROI thisUV_ROI = XYtoUV_ROI(currHit);

    GetModAPVnum(thisUV_ROI.
    );



    apvInfoVals minU_APVdef= GetAPV(minModID, 0, minUAPVnum))

    thisUV_ROI.lineNumber=lineCount;
    thisUV_ROI.GEMLayer=currHit.GEMLayer;


    TVector2 minUVpos = TVector2(thisUV_ROI.Umin, thisUV_ROI.Vmin);
    
    TVector2 maxUVpos = TVector2(thisUV_ROI.Umax, thisUV_ROI.Vmax);


    //TODO: Dont worry about inversion! only dictates direction of readout of APV strips 

    thisUV_ROI.uMinAPVid = UorVtoAPVid(thisUV_ROI.MinModID, minUVpos, 0);
    thisUV_ROI.vMinAPVid = UorVtoAPVid(thisUV_ROI.MinModID, minUVpos, 1);
    
    thisUV_ROI.uMaxAPVid = UorVtoAPVid(thisUV_ROI.MaxModID, maxUVpos, 0);
    thisUV_ROI.vMaxAPVid = UorVtoAPVid(thisUV_ROI.MaxModID, maxUVpos, 1);


    //TODO: make it so that all APVs btwn min max are explicitly active?

    thisUV_ROI.print();
    

    ROI_APV_map.push_back(thisUV_ROI);
    // ROI_APV_map.push(thisAPV);

    lineCount++;

    //NOTE: I dont think I need b/c redeclared every loop
    // currHit.clear();
    // thisUV_ROI.clear();

  }


  hitFile.close();
  // return 0;

  
}


