// This script takes the local GEM database and ROI text file as inputs, and outputs APV lists, per ECal bin. The exact output format is TBD.

#include "DBread.h"
#include "ROIread.h"
#include "GEMModROItoStrips.h"



int generate_ROI_APVlist( const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt" )
{

	DBread db{ db_local };

	ROIread roi{ roi_file };

	if ( db.returnFileReadStatus() == -1 || roi.returnFileReadStatus() == -1 )
	{
		std::cerr << "Exiting the program!!!" << endl << endl;

		return -1;
	}

	std::map< int, GEMLayer > map_GEMLayers = db.returnGEMLayerMap();

	std::map< int, std::map< int, ROI > > map_ROIsByBinsAndLayers = roi.return_ROIMap();

	int nECalBins = map_ROIsByBinsAndLayers.size();
	std::cout << "N ECal bins: " << nECalBins << endl;

	int nGEMLayers = map_GEMLayers.size();
	std::cout << "N GEM layers: " << nGEMLayers << endl;

	std::map< int, GEMLayerROItoStrips > map_GEMLayerROItoStrips;

	for ( const auto& [layerNum, gemLayerInstance] : map_GEMLayers )
	{
		map_GEMLayerROItoStrips.emplace( layerNum, GEMLayerROItoStrips{gemLayerInstance} );
	}

	// FINAL OUTPUT FOR Rafael to convert to GEM VTP, FIBER, and APV info.
	std::map < int /*ECalBinNo*/,
	 std::map< int /*GEMModNumAsInDB*/,
	  std::pair< std::set<int> /*setOfUstripsForModule*/, std::set<int> /*SetofVstripsforMOdule*/> > > map_physicalUVStrips_byECalBin_byGEMMod;

	for ( auto& [binNum, map_ROIbyLayer_forThisBin] : map_ROIsByBinsAndLayers )
	{
		std::map <int, std::pair < std::set<int>, std::set<int> >> map_allUVstripsSetsForAllModules_forThisECalBin;

		//std::cout << "***Filling bin: " << binNum << endl;

		for ( auto& [layerNum, gemLayerROItoStripsInstance] : map_GEMLayerROItoStrips )
		{
			std::map <int, std::pair < std::set<int>, std::set<int> >> map_thisLayer_allROIsForAllModules_forThisECalBin = gemLayerROItoStripsInstance.takeROI_givePhysicalUVStrips( /*map_ROIbyLayer_forThisBin.at( layerNum )*/ (map_ROIsByBinsAndLayers.at(binNum).at(layerNum)) );
			
			//std::cout << "Layer Num: " << layerNum << "  Numer of ROI modules: " << map_thisLayer_allROIsForAllModules_forThisECalBin.size() << endl;
			
			for (const auto& [modNum, modUVstripPair] : map_thisLayer_allROIsForAllModules_forThisECalBin) 
			{
				// If key doesn't exist in finalMap, insert it
			    if (map_allUVstripsSetsForAllModules_forThisECalBin.find(modNum) == map_allUVstripsSetsForAllModules_forThisECalBin.end())
			    {
			    	map_allUVstripsSetsForAllModules_forThisECalBin[modNum] = modUVstripPair;
			    } 
			}
		}

		map_physicalUVStrips_byECalBin_byGEMMod[binNum] = map_allUVstripsSetsForAllModules_forThisECalBin;

	}

	std::ofstream outchanfile( "outPhysicalStrips.txt" );

	for ( const auto& [binNum, uvStripSetbyModule] : map_physicalUVStrips_byECalBin_byGEMMod )
	{
		outchanfile << endl << "### BIN NUMBER: " << binNum << " ###" << endl;

		for ( const auto& [modNUm, uvStripPair] : uvStripSetbyModule )
		{
			outchanfile << endl << "*** Mod Num: " << modNUm << " ***" << endl;

	 		for ( const auto& ustrip : uvStripPair.first ) outchanfile << "U strip: " << ustrip << endl;

			for ( const auto& vstrip : uvStripPair.second ) outchanfile << "V strip: " << vstrip << endl;
		}
	}

	outchanfile.close();

	
	////
	// std::map<int, GEMModule> layer_gem_map = (map_GEMLayers.at(2)).gemModulesLayer;

	// std::cout << "NGEM modules: " << layer_gem_map.size() << endl;

	// GEMLayerROItoStrips gemLayer { map_GEMLayers.at(7) }; //*

	// // GEMModROItoStrips gemToStrips { layer_gem_map.at(2)  }; //((map_GEMLayersAndModules[2]).gemModulesLayer)[0]

	// // if ( gemToStrips.isROIWithinModule( (map_ROIsByBinsAndLayers.at(110).at(2)) ) ) std::cout << "ROI is within layer" << endl;
	// // else std::cout << "ROI is outside layer" << endl;

	// // std::pair < std::set<int>, std::set<int> > stripSetMap = gemToStrips.calcAndReturn_UandVstripSetPair_forROI( (map_ROIsByBinsAndLayers.at(110).at(2)) );

	// // std::cout << "N U strips: " << stripSetMap.first.size() << endl;
	// // std::cout << "N V strips: " << stripSetMap.second.size() << endl;

	// std::map <int, std::pair < std::set<int>, std::set<int> >> stripSetMap = gemLayer.takeROI_givePhysicalUVStrips( (map_ROIsByBinsAndLayers.at(6).at(7)) ); //*

	// std::map <int, std::pair < std::set<int>, std::set<int> >> stripSetMap2 = gemLayer.takeROI_givePhysicalUVStrips( (map_ROIsByBinsAndLayers.at(5).at(7)) );

	// // for ( const auto& ustrip : stripSetMap.at(2).first )
	// // {
	// // 	std::cout << "U strip: " << ustrip << endl;
	// // }

	// // for ( const auto& vstrip : stripSetMap.at(2).second )
	// // {
	// // 	std::cout << "V strip: " << vstrip << endl;
	// // }

	// std::cout << "Number of mods in ROI bin 6: " << stripSetMap.size() << endl; //*
	// std::cout << "Number of mods in ROI bin 5: " << stripSetMap2.size() << endl; //*
	// for ( const auto& [modnum, pair] : stripSetMap )
	// {
	// 	std::cout << endl << "*** Mod Num: " << modnum << " ***" << endl;

	// 	for ( const auto& ustrip : pair.first ) std::cout << "U strip: " << ustrip << endl;

	// 	for ( const auto& vstrip : pair.second ) std::cout << "V strip: " << vstrip << endl;
	// }

	////

	


	return 0;
}