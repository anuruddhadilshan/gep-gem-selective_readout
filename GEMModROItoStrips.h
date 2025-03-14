#include "DBread.h"

#ifndef GEMMODROITOSTRIPS_H
#define GEMMODROITOSTRIPS_H



class GEMModROItoStrips
{

private:

const GEMModule& m_GEMModule;

const double m_fPxU;
const double m_fPyU;
const double m_fPxV;
const double m_fPyV;

std::set < int > m_uStrips; // Set to store the list of U strips within the ROI.
std::set < int > m_vStrips; // Set to store the list of V strips within the ROI.
std::pair < std::set<int>, std::set<int> > m_uvStrips_pair;

public:

	GEMModROItoStrips( const GEMModule& GEM_module ) 
	: m_GEMModule{GEM_module}, m_fPxU{cos(m_GEMModule.uAngle*TMath::DegToRad())}, m_fPyU{sin(m_GEMModule.uAngle*TMath::DegToRad())},
	m_fPxV{cos(m_GEMModule.vAngle*TMath::DegToRad())}, m_fPyV{sin(m_GEMModule.vAngle*TMath::DegToRad())}
	{	
		
	}

	bool isROIWithinModule( const ROI& roi )
	{
		return ( roi.xMin >= m_GEMModule.xMin && roi.xMax <= m_GEMModule.xMax && roi.yMin >= m_GEMModule.yMin && roi.yMax <= m_GEMModule.yMax );
	}

	TVector2 XYtoUV( TVector2 XY )
	{
	  double Xtemp = XY.X();
	  double Ytemp = XY.Y();

	  double Utemp = Xtemp*m_fPxU + Ytemp*m_fPyU;
	  double Vtemp = Xtemp*m_fPxV + Ytemp*m_fPyV;

	  return TVector2(Utemp, Vtemp);
	}

	int return_Ustrip( double hitpos )
	{
  		int ustrip = static_cast<int>(std::round( (hitpos - m_GEMModule.uOffset)/m_GEMModule.pitch + 0.5*m_GEMModule.nStripsu - 0.5 ) );

  		if ( ustrip >= 0 ) return ustrip;
  		else return 0;
	}

	int return_Vstrip( double hitpos )
	{
  		int vstrip = static_cast<int>(std::round( (hitpos - m_GEMModule.vOffset)/m_GEMModule.pitch + 0.5*m_GEMModule.nStripsv - 0.5 ) );

  		if ( vstrip >= 0 ) return vstrip;
  		else return 0;
	}

	std::pair< std::set<int>, std::set<int> > calcAndReturn_UandVstripSetPair_forROI( const ROI& roi_const )
	{
		// Let's make a copy 
		ROI roi = roi_const;

		// First we have to substract the GEM module position offsets in the X and Y directions from ROI boundary values.
		// So that they are defined w.r.t the GEM module coordinate system.


		roi.xMin = roi.xMin - m_GEMModule.xPos;
		roi.xMax = roi.xMax - m_GEMModule.xPos;
		roi.yMin = roi.yMin - m_GEMModule.yPos;
		roi.yMax = roi.yMax - m_GEMModule.yPos;


		// Now this for XY layers as ROI could be spread across many GEM modules in the layer.
		// IF both the xMin and xMax value of the ROI is OUTSIDE the module X value limits, that must mean the module is COMPLETELY
		// outside the ROI. So we will just return two empty sets.
		//if ( (roi.xMin < m_GEMModule.xMin || roi.xMin > m_GEMModule.xMax) && (roi.xMax < m_GEMModule.xMin || roi.xMax > m_GEMModule.xMax) )
		if ( (roi.xMax > m_GEMModule.xMax && roi.xMin > m_GEMModule.xMax) || (roi.xMax < m_GEMModule.xMin && roi.xMin < m_GEMModule.xMin) )
		{
			m_uStrips.clear();
			m_vStrips.clear();
			m_uvStrips_pair = std::make_pair(m_uStrips, m_vStrips);

			return m_uvStrips_pair; // Return an pair with empty sets.
		}

		if ( roi.xMin < m_GEMModule.xMin ) roi.xMin = m_GEMModule.xMin;
		if ( roi.xMax > m_GEMModule.xMax ) roi.xMax = m_GEMModule.xMax;
		if ( roi.yMin < m_GEMModule.yMin ) roi.yMin = m_GEMModule.yMin;
		if ( roi.yMax > m_GEMModule.yMax ) roi.yMax = m_GEMModule.yMax;
 		
		TVector2 XY_ROI_1 { roi.xMax, roi.yMax };
		TVector2 XY_ROI_2 { roi.xMax, roi.yMin };
		TVector2 XY_ROI_3 { roi.xMin, roi.yMin };
		TVector2 XY_ROI_4 { roi.xMin, roi.yMax };

		TVector2 UV_ROI_1 = XYtoUV( XY_ROI_1 );
		TVector2 UV_ROI_2 = XYtoUV( XY_ROI_2 );
		TVector2 UV_ROI_3 = XYtoUV( XY_ROI_3 );
		TVector2 UV_ROI_4 = XYtoUV( XY_ROI_4 );

		std::map < int, std::pair< int, int > > roi_UVstrip;

		roi_UVstrip.emplace( 1, std::make_pair( return_Ustrip(UV_ROI_1.X()), return_Vstrip(UV_ROI_1.Y()) ) );
		roi_UVstrip.emplace( 2, std::make_pair( return_Ustrip(UV_ROI_2.X()), return_Vstrip(UV_ROI_2.Y()) ) );
		roi_UVstrip.emplace( 3, std::make_pair( return_Ustrip(UV_ROI_3.X()), return_Vstrip(UV_ROI_3.Y()) ) );
		roi_UVstrip.emplace( 4, std::make_pair( return_Ustrip(UV_ROI_4.X()), return_Vstrip(UV_ROI_4.Y()) ) );

		int roi_U_strip_max {-1};
  		int roi_U_strip_min {100000}; 
  		int roi_V_strip_max {-1};
  		int roi_V_strip_min {100000};

  		for ( const auto& [key, pair] : roi_UVstrip )
  		{
		    //std::cout << "ROI_" << key << ": (" << pair.first << ", " << pair.second << ")" << endl;

		    roi_U_strip_max = ( pair.first > roi_U_strip_max ) ? pair.first : roi_U_strip_max;
		    roi_V_strip_max = ( pair.second > roi_V_strip_max ) ? pair.second : roi_V_strip_max;

		    roi_U_strip_min = ( pair.first < roi_U_strip_min ) ? pair.first : roi_U_strip_min;
		    roi_V_strip_min = ( pair.second < roi_V_strip_min ) ? pair.second : roi_V_strip_min;
  		}

  		// std::cout << "U strip " << roi_U_strip_min << " " << roi_U_strip_max << endl;
  		// std::cout << "V strip " << roi_V_strip_min << " " << roi_V_strip_max << endl;

  		m_uStrips.clear();
  		for ( int iUstrip = roi_U_strip_min; iUstrip </* = */ roi_U_strip_max; iUstrip++ ) m_uStrips.insert(iUstrip);
  		
  		m_vStrips.clear();
  		for ( int iVstrip = roi_V_strip_min; iVstrip </* = */ roi_V_strip_max; iVstrip++ ) m_vStrips.insert(iVstrip);
  		
  		m_uvStrips_pair = std::make_pair(m_uStrips, m_vStrips);
  		
  		return m_uvStrips_pair;
  		
	}

};

class GEMLayerROItoStrips
{
private:

	const GEMLayer& m_GEMLayer;

	const int m_LayerNum;

	const std::map< int, GEMModule >& m_GEMModulesInLayerMap;

	const int m_NmodulesInLayer;

	std::map< int, GEMModROItoStrips > m_GEMModulesROItoStripsMapLayer;

	std::map< int, std::pair< std::set<int>, std::set<int> > > m_LayerModule_UVstripPairs;

public:

	GEMLayerROItoStrips( const GEMLayer& GEM_layer ) : m_GEMLayer{ GEM_layer }, m_LayerNum{ m_GEMLayer.layerNum },
	m_GEMModulesInLayerMap{ m_GEMLayer.gemModulesLayer }, m_NmodulesInLayer{ static_cast<int>(m_GEMModulesInLayerMap.size()) }
	{
		for ( const auto& [mod_num, gemMod] : m_GEMModulesInLayerMap )
		{
			m_GEMModulesROItoStripsMapLayer.emplace( mod_num, GEMModROItoStrips{gemMod} );			
		}

		//std::cout << "Number of GEM modules: " << m_NmodulesInLayer << endl;

		//for ( const auto& [mod_num, class_1] : m_GEMModulesROItoStripsMapLayer ) std::cout << "Mod num: " << mod_num << endl;

	}

	std::map< int, std::pair< std::set<int>, std::set<int> > >& takeROI_givePhysicalUVStrips( const ROI& roi )
	{

		if ( m_NmodulesInLayer == 1 ) 
		{
			int modNum = m_GEMModulesROItoStripsMapLayer.begin()->first;

			std::pair< std::set<int>, std::set<int> > moduleUVstripPair = (m_GEMModulesROItoStripsMapLayer.at(modNum)).calcAndReturn_UandVstripSetPair_forROI( roi );

			m_LayerModule_UVstripPairs[modNum] = moduleUVstripPair;

			return m_LayerModule_UVstripPairs;

		}

		else if ( m_NmodulesInLayer == 4 )
		{
			for ( auto& [modNum, gemModROItoStrips] : m_GEMModulesROItoStripsMapLayer )
			{
				std::pair< std::set<int>, std::set<int> > moduleUVstripPair = gemModROItoStrips.calcAndReturn_UandVstripSetPair_forROI( roi );

				if ( !moduleUVstripPair.first.empty() && !moduleUVstripPair.second.empty()  ) m_LayerModule_UVstripPairs[modNum] = moduleUVstripPair;
				else m_LayerModule_UVstripPairs.erase(modNum);
			}

			return m_LayerModule_UVstripPairs;
		}

		else
		{
			std::cerr << "ERROR, GEMLayerROItoStrips - INVALID NUMBER OF MODULES IN THE GEM LAYER" << endl;

			return m_LayerModule_UVstripPairs;
		}

	}

};



#endif