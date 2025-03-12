#ifndef FindAPV_H
#define FindAPV_H

// #include "THaSubDetector.h"
#include <vector>
#include <set>
#include <map>
#include <array>
#include <deque>

// #include "TVector2.h"
// #include "TVector3.h"

//using namespace std;

class THaDetectorBase;
class THaEvData;
class THaRunBase;
class TH1D;
class TH2D;
class TF1;
class TClonesArray;


namespace APVFindingSpace{
      

class APVFinder{
public:
    APVFinder();

    // APVFinder(int nAPV25_CHAN, int MPDMAP_ROW_SIZE, int MAXNSAMP_PER_APV, int MAX2DHITS, int APVmapping);

    virtual ~APVFinder() {};

    struct XY_ROI{
        int ECalBin;
        int GEMLayer;
        double xMin;
        double xMax;
        double yMin;
        double yMax;
        int lineNumber;
        int hitNumber;
      
        void print() const;
    };
      
    struct UV_ROI{
        int GEMLayer;
        int lineNumber;
        
        int MinModID;
        int MaxModID;
        
        //UV hit position coordinates
        double Umin;  double Umax;
        double Vmin;  double Vmax;
        
        //positions along their axes
        int UminPos;  int UmaxPos;
        int VminPos;  int VmaxPos;
        
        int uMinAPVid;  int vMinAPVid;
        
        int uMaxAPVid;
        int vMaxAPVid;
        void print() const;
    };
    // std::vector<UV_ROI> ROI_APV_map;


    double stripOffset;

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
      

      
        void print() const;
      
        };

    // std::vector<int> modIDs;//simply store indices from 0 to any module

    std::map<int, gemInfo> gemInfoMap;


    int nMods; //number of Modules in FT

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
        
        void print() const ;
        
    };

    //module id, axis, position on axis
    struct apvInfoKeys {int gemid, axis, pos;
        void print() const;
        bool operator<(const apvInfoKeys& other) const;
    };
      
    //vtpcrate, fiber, adc_ch
    struct apvInfoVals{
    //fill these values they id the APV
    int vtpcrate;
    int fiber; //NOTE: == mpd_id
    int adc_ch;
    int invert;
    void print() const;
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
    
    struct roiAPVinfo{
        apvInfoKeys uMaxAPVinfoKeys;
        apvInfoKeys uMinAPVinfoKeys;
        apvInfoKeys vMaxAPVinfoKeys;
        apvInfoKeys vMinAPVinfoKeys;


        apvInfoVals uMaxAPVinfoVals;
        apvInfoVals uMinAPVinfoVals;
        apvInfoVals vMaxAPVinfoVals;
        apvInfoVals vMinAPVinfoVals;

        void print() const;
    };

          
    struct roiAPVinfoVals{
        int GEMLayer;
        int MinModID;
        int MaxModID;
        int UminAPVid;
        int VminAPVid;
        int UmaxAPVid;
        int VmaxAPVid;
        int lineNumber;
        int hitEntry;
      };
      
      struct roiAPVinfoKeys{
        int modID;
        int axis;
        int pos;
      };
  
  std::map<int, roiAPVinfo> roiAPVmap;

  void printRoiAPVMap();
  void OutputRoiAPVMap();

    // APV Configuration Constants
    int fN_APV25_CHAN;
    int fMPDMAP_ROW_SIZE;
    int fMAXNSAMP_PER_APV;
    int fMAX2DHITS;
    int fAPVmapping;
    
    
    //Arrays to temporarily hold raw data from ONE APV card:
    std::vector<int> fStripAPV;
    std::vector<int> fRawStripAPV;
    std::vector<int> fRawADC_APV;
    
    double UAngle, VAngle;
    
    // TVector3 alignedPos;  TVector3 alignedAng;


 //Decode map information: 
    std::vector<int>       fChanMapData;

    std::array<std::vector<int>, 4 > APVMAP;

    //some convenience maps: are these actually used yet? 
    std::map<int, int> fAPVch_by_Ustrip;
    std::map<int, int> fAPVch_by_Vstrip;
    std::map<int, int> fMPDID_by_Ustrip;
    std::map<int, int> fMPDID_by_Vstrip;
    std::map<int, int> fADCch_by_Ustrip;
    std::map<int, int> fADCch_by_Vstrip;

    //Constant, module-specific parameters:
    int fModule; // Module index within a tracker. Should be unique! Since this is a GEM module class, this parameter should be unchanging
    int fLayer;  // Layer index of this module within a tracker. Since this is a GEM module class, this parameter should be unchanging

    int fNstripsU; // Total number of strips in this module along the generic "U" axis
    int fNstripsV; // Total number of strips in this module along the generic "V" axis

    //To be determined from channel map/strip count information:
    int fNAPVs_U; //Number of APV cards per module along "U" strip direction; this is typically 8, 10, or 12, but could be larger for U/V GEMs
    int fNAPVs_V; //Number of APV cards per module along "V" strip direction; 
    //UInt_t fNTimeSamplesADC; //Number of ADC time samples (this could be variable in principle, but should be the same for all strips within a module within a run) redundant with fN_MPD_TIME_SAMP

    std::vector<int> fStrip;  //Strip index of hit (these could be "U" or "V" generalized X and Y), assumed to run from 0..N-1
    // std::vector<SBSGEM::GEMaxis_t>  fAxis;  //We just made our enumerated type that has two possible values, makes the code more readable (maybe)


    int fChan_TimeStamp_low;
    int fChan_TimeStamp_high;
    int fChan_MPD_EventCount;



    // int extraReadDB(const TDatime& date );
    // int ReadDB(const TDatime& date );


 //GEOMETRICAL PARAMETERS:
    double fUStripPitch;    //strip pitch along U, will virtually always be 0.4 mm
    double fVStripPitch;    //strip pitch along V, will virtually always be 0.4 mm
    double fUStripOffset;   //position of first U strip along the direction it measures:
    double fVStripOffset;   //position of first V sttrip alogn the direction it measures:
    double fUAngle;         //Angle between U strips and "X" axis of TRANSPORT coordinates;
    double fVAngle;         //Angle between V strips and "X" axis of TRANSPORT coordinates;
    double fPxU;            //U Strip X projection = cos( UAngle );
    double fPyU;            //U Strip Y projection = sin( UAngle );
    double fPxV;            //V Strip X projection = cos( VAngle );
    double fPyV;            //V Strip Y projection = sin( VAngle );


    //Aligned geometry
    // TVector3 alignPos; 
    // TVector3 alignAngs; 



 //Geometry related Funcs

    void fillGEMInfoMap();
    
    void printGEMinfoMap();
    void OutputGEMinfoMap();

    int SetInvert(int invert);
    
    void MakeModGeomRef(const char* refFile);


    int GetLayerOfMod(int modID);
    int GetMod_apvmap(int modID);

 //APV Reference Map functions
    void MakeRefMap(const char* refFile);

    void printAPVinfoMap();
    
    void OutputAPVinfoMap();
    


// Start hit file processing###################################

 
    XY_ROI LoadLine(std::string &LineStr);
    
    apvInfoVals GetAPV(int gemid, int axis, int pos);

    std::array<double, 3> GetModDimensions(int apvmap);


    // std::map<int, TVector3> ModPositionmap;
    std::map<int, std::array<double, 3>> ModPositionmap;


 //Get geometry relevant to hit 
    std::array<double, 2> GetUVang(int modID);

    std::array<double, 2> LayerUVang(int layer);


    void SetProjOps(std::array<double, 2> UVangles);

    double GetPitch(){return .0004;}

    // double GetUVoffset(int modID, char *axis);
    // double GetUVoffset(int modID, const std::string& axis);
    double GetUVoffset(int modID, int axis);

    std::array<double, 2> GetOffset(int modID);

    // int GetNstrips(int modID, char *axis);
    int GetNstrips(int modID, int);

    int GetStripNumber( int rawstrip, int pos, int invert );

    std::array<int, 2> GetNstrips(int modID);

    std::array<int, 2> GetStripRange();

    std::array<int, 2> GetNAPVs(int modID);


    // TVector2 GetLocalStripIDs(int modID, TVector2 UVhit);

    int UorVtoAPVid(int modID, std::array<int, 2> localUVPos,
        int axis //U/X=0, V/Y=1
      );


    // TODO: for now just copied from SBSGEMModule 
    std::array<double, 2> XYtoUV( std::array<double, 2> XY );

    // int GetModId(int layerID, std::array<int, 2> hitPos, TVector3 modDims);
    int GetModId(int layerID, double hitX);

    UV_ROI XYtoUV_ROI(XY_ROI hitXY);

    bool isAboveStripCenter(int modID, double xLocal);

      
    double GetXstripCenterLocal(int layerID, double xGlobal);
    double GetXLocal(int layerID, double xGlobal);
    // std::array<double, 2> GetXYLocal(XY_ROI theXYROI);
    std::array<double, 2> GetXYLocal(int layerID, double xGlobal, double yGlobal);

    // int uvCoordToPos(double Coord);//UVdouble positions to INTEGER position on axis assuming "-x" is positive for u,v
    
    // int GetAPVpos(int layerID, int axis, double Coord);//UVdouble positions to INTEGER position on axis assuming "-x" is positive for u,v

    // int GetAPVpos(UV_ROI theUVROI);//UVdouble positions to INTEGER position on axis assuming "-x" is positive for u,v

    // int GetAPVpos(int layerID, double xGlobal, double yGlobal, int axis);
    int GetAPVpos(int modID, double xGlobal, double yGlobal, int axis);

    // UInt_t UV_ROItoAPVs(int modID, TVector2 UVhit);

    // UV_ROI XYtoAPV(UInt_t layer, TVector2 hitPos);



    std::map<double, std::array<int, 2>> xToUVmap;
    std::map<double, std::array<int, 2>> yToUVmap;

    double Udx, Udy; //dx for 1 u strip, dy for 1 u strip
    double Vdx, Vdy; //dx for 1 v strip, dy for 1 v strip

    // double GetUdx(int layerID); //get length in x direction of a U strip
    double GetUdx(int modID); //get length in x direction of a U strip
    // double GetVdx(int layerID); //get length in x direction
    double GetVdx(int modID); //get length in x direction


    int hitposToStripID( std::array<int, 2> hitpos, int &Nstrips );

    //TODO: make a clear func
    // void Clear( Option_t* opt);

   //Reuse every hit so no global needed
    // TVector2 minXminY;
    // TVector2 minXmaxY;
    // TVector2 maxXminY;
    // TVector2 maxXmaxY;
    
    // std::vector<TVector2> roiSquare;

    void TestGEMInfoMap();
    void TestAPVInfoMap(const char* refFile);


    // void FindAPV(const TDatime& date, const char *aHitFile);
    void FindAPV();

};

    
//   }
  




}
#endif