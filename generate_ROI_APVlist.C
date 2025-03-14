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
		std::cerr << "Exiting the program!!!" << std::endl << std::endl;

		return -1;
	}

	std::map< int, GEMLayer > map_GEMLayers = db.returnGEMLayerMap();

	std::map< int, std::map< int, ROI > > map_ROIsByBinsAndLayers = roi.return_ROIMap();

	int nECalBins = map_ROIsByBinsAndLayers.size();
	std::cout << "N ECal bins: " << nECalBins << std::endl;

	int nGEMLayers = map_GEMLayers.size();
	std::cout << "N GEM layers: " << nGEMLayers << std::endl;

	std::map< int, GEMLayerROItoStrips > map_GEMLayerROItoStrips;

	for ( const auto& [layerNum, gemLayerInstance] : map_GEMLayers )
	{
		map_GEMLayerROItoStrips.emplace( layerNum, GEMLayerROItoStrips{gemLayerInstance} );
	}











	
	// FINAL OUTPUT FOR Rafael to convert to GEM VTP, FIBER, and APV info.
	
	std::map < int /*ECalBinNo*/,
	std::map< int /*GEMModNumAsInDB*/,
	 std::pair< std::set<int> /*setOfUstripsForModule*/, std::set<int> /*SetofVstripsforMOdule*/> > > map_physicalUVStrips_byECalBin_byGEMMod;


	 //RREDITS:******
	 std::map < int /*ECalBinNo*/, std::set<apvInfoVals>> ECalBinAPVinfo;
	 std::set<apvInfoKeys> missingKeys;




	for ( auto& [binNum, map_ROIbyLayer_forThisBin] : map_ROIsByBinsAndLayers )
	{
		std::map <int, std::pair < std::set<int>, std::set<int> >> map_allUVstripsSetsForAllModules_forThisECalBin;

		//std::cout << "***Filling bin: " << binNum << std::endl;
		
		
		//RREDIT:Start
		std::set<apvInfoVals> currECalBinAPVs;

		//RREDIT:End



		for ( auto& [layerNum, gemLayerROItoStripsInstance] : map_GEMLayerROItoStrips )
		{	
			
			std::map <int, std::pair < std::set<int>, std::set<int> >> map_thisLayer_allROIsForAllModules_forThisECalBin = gemLayerROItoStripsInstance.takeROI_givePhysicalUVStrips( /*map_ROIbyLayer_forThisBin.at( layerNum )*/ (map_ROIsByBinsAndLayers.at(binNum).at(layerNum)) );
			
			//std::cout << "Layer Num: " << layerNum << "  Numer of ROI modules: " << map_thisLayer_allROIsForAllModules_forThisECalBin.size() << std::endl;
			
			for (const auto& [modNum, modUVstripPair] : map_thisLayer_allROIsForAllModules_forThisECalBin) 
			{
				// If key doesn't exist in finalMap, insert it
			    if (map_allUVstripsSetsForAllModules_forThisECalBin.find(modNum) == map_allUVstripsSetsForAllModules_forThisECalBin.end())
			    {
			    	map_allUVstripsSetsForAllModules_forThisECalBin[modNum] = modUVstripPair;
			    } 


				// RREDIT:Start
				std::set<int>uStrips = modUVstripPair.first;
				std::set<int>vStrips = modUVstripPair.second;
				//RREDIT:END

				for (int stripID : uStrips){
					apvInfoKeys currStripKeys;
					currStripKeys.gemid = modNum;
					currStripKeys.axis = 0;


					currStripKeys.pos = stripID/128;//APV position
					// posInAPV = stripID%128;

					std::cout << "ECalBin: " << binNum << " Checking modID: " << modNum << " Strip: " << stripID 
					<< " Axis: " << currStripKeys.axis << " Pos: " << currStripKeys.pos << std::endl;
					if (apvInfoMap.find(currStripKeys) == apvInfoMap.end()) {
						std::cout << "\nFor ECal Bin " << binNum <<" For strip " << stripID  << "  --> APV key: " << currStripKeys.gemid <<", "<< currStripKeys.axis <<", "<<
						currStripKeys.pos << " NOT FOUND!\n" << std::endl;
						missingKeys.emplace(currStripKeys);
					}
		
		
					if (apvInfoMap.find(currStripKeys) != apvInfoMap.end()) {
						currECalBinAPVs.emplace(apvInfoMap[currStripKeys]);  // Store unique strip -> APV mapping
					}

				}

				for (int stripID : vStrips){
					apvInfoKeys currStripKeys;
					currStripKeys.gemid = modNum;
					currStripKeys.axis = 1;


					currStripKeys.pos = (stripID/128);//APV position(-1 so index starts at 0)
					// posInAPV = stripID%128;

					std::cout << "ECalBin: " << binNum << " Checking modID: " << modNum << " Strip: " << stripID 
					<< " Axis: " << currStripKeys.axis << " Pos: " << currStripKeys.pos << std::endl;
					
					if (apvInfoMap.find(currStripKeys) == apvInfoMap.end()) {
						std::cout << "\nFor ECal Bin " << binNum <<" For strip " << stripID  << "  --> APV key: " << currStripKeys.gemid <<", "<< currStripKeys.axis <<", "<<
						currStripKeys.pos << " NOT FOUND!\n" << std::endl;
						missingKeys.emplace(currStripKeys);
					}
		
		
					if (apvInfoMap.find(currStripKeys) != apvInfoMap.end()) {
						currECalBinAPVs.emplace(apvInfoMap[currStripKeys]);  // Store unique strip -> APV mapping
					}

				}
				//RREDITS:End
		
			}

		map_physicalUVStrips_byECalBin_byGEMMod[binNum] = map_allUVstripsSetsForAllModules_forThisECalBin;

		}

		ECalBinAPVinfo[binNum] = currECalBinAPVs;//RREDIT
	}

	std::ofstream outchanfile( "outPhysicalStrips.txt" );

	for ( const auto& [binNum, uvStripSetbyModule] : map_physicalUVStrips_byECalBin_byGEMMod )
	{
		outchanfile << std::endl << "### BIN NUMBER: " << binNum << " ###" << std::endl;

		for ( const auto& [modNUm, uvStripPair] : uvStripSetbyModule )
		{
			outchanfile << std::endl << "*** Mod Num: " << modNUm << " ***" << std::endl;

	 		for ( const auto& ustrip : uvStripPair.first ) outchanfile << "U strip: " << ustrip << std::endl;

			for ( const auto& vstrip : uvStripPair.second ) outchanfile << "V strip: " << vstrip << std::endl;
		}
	}

	outchanfile.close();
	
	

	//RREDIT:START


	std::ofstream APVmapfile( "EcalBinToAPVmap.txt");

	APVmapfile << "# Lists APV information corresponding to each ECal Bin\n";
    APVmapfile << "# Format: ECal Bin | VTP Crate | Fiber (MPD ID) | ADC Channel\n\n";

	for ( auto &[binNum, APVinfo] : ECalBinAPVinfo){
		APVmapfile << std::endl << "### ECal BIN NUMBER: " << binNum << " ###" << std::endl;

		int colWidth = 12;

    	APVmapfile << std::left 
		<< std::setw(colWidth) << "ECalBin"
		<< std::setw(colWidth) << "VTP Crate"
		<< std::setw(colWidth) << "Fiber"
		<< std::setw(colWidth) << "ADC Ch"
		<< "\n";

    	APVmapfile << std::string(7 * colWidth, '-') << "\n";

		for (auto &APV:APVinfo){
			APVmapfile << std::left 
			<< std::setw(colWidth) << binNum
			<< std::setw(colWidth) << APV.vtpcrate
			<< std::setw(colWidth) << APV.fiber
			<< std::setw(colWidth) << APV.adc_ch
			<< "\n";
		}
	}

	APVmapfile.close();
	//RREDITS:End

	std::cout<< "\n\n\n!!!!!!!!!! missing keys!!!!!!!!!!\n" << std::endl;
	for (auto &key: missingKeys){ key.print();}

	


	return 0;
}