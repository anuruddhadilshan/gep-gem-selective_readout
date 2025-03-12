#include "DBread.h"

double fPxU, fPyU, fPxV, fPyV;

double mod_max_x, mod_min_x, mod_max_y, mod_min_y;

double uangle, vangle;

double roi_max_x {0.23};
double roi_min_x {-0.15};
double roi_max_y {0.14};
double roi_min_y {-0.25};

TVector2 XYtoUV( TVector2 XY )
{
  double Xtemp = XY.X();
  double Ytemp = XY.Y();

  double Utemp = Xtemp*fPxU + Ytemp*fPyU;
  double Vtemp = Xtemp*fPxV + Ytemp*fPyV;

  return TVector2(Utemp, Vtemp);
}

int return_istrip( double hitpos, int Nstrips, double pitch, double offset )
{
  return static_cast<int>(std::round( (hitpos - offset)/pitch + 0.5*Nstrips - 0.5 ) );
}

bool isPointInModule( const double x, const double y  )
{
  return ( x <= mod_max_x && x >= mod_min_x && y <= mod_max_y && y >= mod_min_y );
}

  

void testdb( const double x, const double y)
{
  //DBread dbread {};

  uangle = 150.0;
  vangle = -150.0;

  fPxU = cos( uangle * TMath::DegToRad() );
  fPyU = sin( uangle * TMath::DegToRad() );
  fPxV = cos( vangle * TMath::DegToRad() );
  fPyV = sin( vangle * TMath::DegToRad() );

  TVector2 XY{x, y};
    
  TVector2 UV = XYtoUV( XY  );

  std::cout << "U: " << UV.X() << "; V: " << UV.Y() << endl;

  mod_max_x = 0 + 1.5/2;
  mod_min_x = 0 - 1.5/2;

  mod_max_y = 0 + 0.4/2;
  mod_min_y = 0 - 0.4/2;
  
  if ( isPointInModule(x, y) )
  {
    std::cout << "Point is inside Module" << endl;
  }
  else
  {
    std::cout << "Point is outside Module" << endl;
  }

  int  nstrips = 3840;
  double pitch = 0.0004;
  double offset = 0.0108;
  
  std::cout << "U strip no: " << return_istrip(UV.X(), nstrips, pitch, offset ) << endl;
  std::cout << "V strip no: " << return_istrip(UV.Y(), nstrips, pitch, offset ) << endl;

  TVector2 XY_ROI_1 { roi_max_x, roi_max_y };
  TVector2 XY_ROI_2 { roi_max_x, roi_min_y };
  TVector2 XY_ROI_3 { roi_min_x, roi_min_y };
  TVector2 XY_ROI_4 { roi_min_x, roi_max_y };

  TVector2 UV_ROI_1 = XYtoUV( XY_ROI_1 );
  TVector2 UV_ROI_2 = XYtoUV( XY_ROI_2 );
  TVector2 UV_ROI_3 = XYtoUV( XY_ROI_3 );
  TVector2 UV_ROI_4 = XYtoUV( XY_ROI_4 );

  std::cout << "UV ROI 1 " << UV_ROI_1.X() << " " << UV_ROI_1.Y() << endl;
  std::cout << "UV ROI 2 " << UV_ROI_2.X() << " " << UV_ROI_2.Y() << endl;
  std::cout << "UV ROI 3 " << UV_ROI_3.X() << " " << UV_ROI_3.Y() << endl;
  std::cout << "UV ROI 4 " << UV_ROI_4.X() << " " << UV_ROI_4.Y() << endl;

  std::map < int, std::pair< int, int > > roi_UV;

  roi_UV.emplace( 1, std::make_pair( return_istrip(UV_ROI_1.X(), nstrips, pitch, offset), return_istrip(UV_ROI_1.Y(), nstrips, pitch, offset) ) );
  roi_UV.emplace( 2, std::make_pair( return_istrip(UV_ROI_2.X(), nstrips, pitch, offset), return_istrip(UV_ROI_2.Y(), nstrips, pitch, offset) ) );
  roi_UV.emplace( 3, std::make_pair( return_istrip(UV_ROI_3.X(), nstrips, pitch, offset), return_istrip(UV_ROI_3.Y(), nstrips, pitch, offset) ) );
  roi_UV.emplace( 4, std::make_pair( return_istrip(UV_ROI_4.X(), nstrips, pitch, offset), return_istrip(UV_ROI_4.Y(), nstrips, pitch, offset) ) );

  int roi_U_strip_max {-1};
  int roi_U_strip_min {10000}; 
  int roi_V_strip_max {-1};
  int roi_V_strip_min {10000};

  for ( const auto& [key, pair] : roi_UV )
  {
    std::cout << "ROI_" << key << ": (" << pair.first << ", " << pair.second << ")" << endl;

    roi_U_strip_max = ( pair.first > roi_U_strip_max ) ? pair.first : roi_U_strip_max;
    roi_V_strip_max = ( pair.second > roi_V_strip_max ) ? pair.second : roi_V_strip_max;

    roi_U_strip_min = ( pair.first < roi_U_strip_min ) ? pair.first : roi_U_strip_min;
    roi_V_strip_min = ( pair.second < roi_V_strip_min ) ? pair.second : roi_V_strip_min;
  }

  std::cout << endl;
  std::cout << "Max U strip: " << roi_U_strip_max << endl;
  std::cout << "Min U strip: " << roi_U_strip_min << endl;
  std::cout << "Max V strip: " << roi_V_strip_max << endl;
  std::cout << "Min V strip: " << roi_V_strip_min << endl;

}