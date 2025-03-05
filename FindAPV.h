#ifndef FindAPV_H
#define FindAPV_H

#include "THaSubDetector.h"
#include <vector>
#include <set>
#include <map>
#include <array>
#include <deque>

#include "TVector2.h"
#include "TVector3.h"

//using namespace std;

class THaDetectorBase;
class THaEvData;
class THaRunBase;
class TH1D;
class TH2D;
class TF1;
class TClonesArray;


namespace APVFindingSpace {
// namespace APVFinder {

    //Based on SBSGEMModule

    // class APVFinder{
        // public:
  
    
  
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
      
        void print() const;
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
        void print() const;
    };

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
        
        // Double_t fPxU = cos(uangle);            //U Strip X projection = cos( UAngle );
        // Double_t fPyU = sin(uangle);            //U Strip Y projection = sin( UAngle );
        // Double_t fPxV = cos(vangle);            //V Strip X projection = cos( VAngle );
        // Double_t fPyV = sin(vangle);            //V Strip Y projection = sin( VAngle );
      
        // // TVector3 alignedPos;  TVector3 alignedAng;
      
        void print() const;
      
        };

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
        bool operator<(const apvInfo& other) const;
        
        void print() const ;
        
    };

    //module id, axis, position on axis
    struct apvInfoKeys {int gemid, axis, pos;};
      
    //vtpcrate, fiber, adc_ch
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

  bool operator<(const stripInfo& other) const ;

  void print() const;
};

std::set<stripInfo> stripInfoMap;


      



class APVFinder{
public:
    APVFinder();

    APVFinder(int nAPV25_CHAN, int MPDMAP_ROW_SIZE, int MAXNSAMP_PER_APV, int MAX2DHITS, int APVmapping);

    virtual ~APVFinder();

    void printAPVinfoMap();

    std::vector<UV_ROI> ROI_APV_map;

    int fN_APV25_CHAN;
    int fMPDMAP_ROW_SIZE;
    int fMAXNSAMP_PER_APV;

    //Arrays to temporarily hold raw data from ONE APV card:
    std::vector<UInt_t> fStripAPV;
    std::vector<UInt_t> fRawStripAPV;
    std::vector<Int_t> fRawADC_APV;

    int fMAX2DHITS;
    int fAPVmapping;

    //   class APVFinder{

    //first term is module ID
    std::map<int, gemInfo> gemInfoMap;

    //GEOMETRICAL PARAMETERS:
    Double_t fUStripPitch;    //strip pitch along U, will virtually always be 0.4 mm
    Double_t fVStripPitch;    //strip pitch along V, will virtually always be 0.4 mm
    Double_t fUStripOffset;   //position of first U strip along the direction it measures:
    Double_t fVStripOffset;   //position of first V sttrip alogn the direction it measures:
    Double_t fUAngle;         //Angle between U strips and "X" axis of TRANSPORT coordinates;
    Double_t fVAngle;         //Angle between V strips and "X" axis of TRANSPORT coordinates;

public: //Projection operators
    Double_t fPxU;            //U Strip X projection = cos( UAngle );
    Double_t fPyU;            //U Strip Y projection = sin( UAngle );
    Double_t fPxV;            //V Strip X projection = cos( VAngle );
    Double_t fPyV;            //V Strip Y projection = sin( VAngle );
    
    
    // TVector3 alignedPos;  TVector3 alignedAng;

    double UAngle, VAngle;


public:
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

    Int_t extraReadDB(const TDatime& date );
    Int_t ReadDB(const TDatime& date );


public: //GEOMETRICAL PARAMETERS:
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




    Int_t ReadGeometry( FILE *file, const TDatime &date, Bool_t required );

    void DefineAxes( Double_t rotation_angle );

    void MakeRefMap(const char* refFile);

    void writeOutRefMap();

    XY_ROI LoadLine(std::string &LineStr);
    
    apvInfoVals GetAPV(int gemid, int axis, int pos);

    TVector3 GetModDimensions(int apvmap);

    std::map<int, TVector3> ModPositionmap;

    public:
    TVector2 GetUVang(int modID);

    TVector2 LayerUVang(int layer);

    public: //from SBSGEMModule
        Double_t fPxU;
        Double_t fPyU;
        Double_t fPxV;
        Double_t fPyV;
    void SetProjOps(TVector2 UVangles);

    double GetPitch(){return .0004;}

    double GetUVoffset(int modID, char *axis);

    TVector2 GetOffset(int modID);

    int GetNstrips(int modID, char *axis);

    Int_t GetStripNumber( UInt_t rawstrip, UInt_t pos, UInt_t invert );

    TVector2 GetNstrips(int modID);

    TVector2 GetStripRange();

    TVector2 GetNAPVs(int modID);


    // TVector2 GetLocalStripIDs(int modID, TVector2 UVhit);

    int UorVtoAPVid(int modID, TVector2 localUVPos,
        int axis //U/X=0, V/Y=1
      );

    // TODO: for now just copied from SBSGEMModule 
    TVector2 XYtoUV( TVector2 XY );

    int GetModId(int layerID, TVector2 hitPos, TVector3 modDims);

    UV_ROI XYtoUV_ROI(XY_ROI hitXY);

    // UInt_t UV_ROItoAPVs(int modID, TVector2 UVhit);

    // UV_ROI XYtoAPV(UInt_t layer, TVector2 hitPos);

    int hitposToStripID( TVector2 hitpos, int &Nstrips );

    void Clear( Option_t* opt);

   //Reuse every hit so no global needed
    // TVector2 minXminY;
    // TVector2 minXmaxY;
    // TVector2 maxXminY;
    // TVector2 maxXmaxY;

    std::vector<TVector2> roiSquare;


    void FindAPV(const TDatime& date, const char *aHitFile);

};

    
//   }
  
  
}




#endif