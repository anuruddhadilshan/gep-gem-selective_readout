#include <TMath.h>
#include <TLine.h>
#include <TString.h>
#include <TH2I.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TRandom.h>


#include <TLatex.h>

#include <map>
#include <set>
#include <string>
#include <iostream>

#include "GEMModROItoStrips.h"
#include "DBread.h"
#include "ROIread.h"
#include "showgemhit_for_ecalbin.C"
// #include "GetGemInfoMap.cxx"

// ---- Geometry settings ----
std::map<int, double> u_angles = {
	{0, 0}, {1, 0},
	{2, 60}, {3, 60}, {4, 60}, {5, 60},
	{6, 0}, {7, 0}
};

std::map<int, double> v_angles = {
	{0, 0}, {1, 0},
	{2, -60}, {3, -60}, {4, -60}, {5, -60},
	{6, 90}, {7, 90}
};

std::map<int, int> nstrips_u = {
	{0, 3968}, {1, 3968},
	{2, 3840}, {3, 3840}, {4, 3840}, {5, 3840},
	{6, 5120}, {7, 5120}
};

std::map<int, int> nstrips_v = {
	{0, 3456}, {1, 3456},
	{2, 3840}, {3, 3840}, {4, 3840}, {5, 3840},
	{6, 1536}, {7, 1536}
};


int uDrawn=0, vDrawn=0;


// ---- Strip drawing function ----
void drawGEMStrip(double mod_x, double mod_y, TString axis, double angle_deg,
                  int strip_num, double pitch = 0.004, int n_strips = 256, double length = 0.8)
{
	double angle = angle_deg * TMath::DegToRad();
	int center_strip = n_strips / 2;
	double offset = (strip_num - center_strip) * pitch;
	double dx = offset * TMath::Cos(angle + TMath::Pi()/2);
	double dy = offset * TMath::Sin(angle + TMath::Pi()/2);
	double x_center = mod_x + dx;
	double y_center = mod_y + dy;
	double x1 = x_center - 0.5 * length * TMath::Cos(angle);
	double y1 = y_center - 0.5 * length * TMath::Sin(angle);
	double x2 = x_center + 0.5 * length * TMath::Cos(angle);
	double y2 = y_center + 0.5 * length * TMath::Sin(angle);


	TMarker *thisHit = new TMarker();

	thisHit->SetX(x_center);
	thisHit->SetY(y_center);

	thisHit->SetMarkerSize(20);
	thisHit->Draw("same");

	TLine *stripLine = new TLine(x1, y1, x2, y2);


	// // stripLine->SetLineWidth(2);
	// stripLine->SetLineWidth(2);
	// stripLine->SetLineColor(axis.EqualTo("U", TString::kIgnoreCase) ? kCyan : kPink);

	if (axis.EqualTo("U", TString::kIgnoreCase)) {
		stripLine->SetLineColor(kCyan);
		stripLine->SetLineWidth(2);
		// stripLine->SetLineStyle(1);  // Solid
		stripLine->SetLineStyle(3);  

		uDrawn+=1;
	} else {
		stripLine->SetLineColor(kPink);
		stripLine->SetLineWidth(1);
		stripLine->SetLineStyle(2);  // Dashed
		vDrawn+=1;
		
	}


	stripLine->Draw("same");
}

// ---- GEM Strip ROI Visualization ----
void showgemstrips_for_ecalbin(const std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerLayer)
// void showgemstrips_for_ecalbin(std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerLayer)
{

	std::cout << "\n\nCalling showgemstrips_for_ecalbin()" <<std::endl;

	static int canvas_id = 0;
	TCanvas* C = new TCanvas(Form("canvas_%d", canvas_id++), "Active ROI Strips in GEM Planes", 1600, 1000);
	C->Divide(4, 2);


	// int uDrawn=0, vDrawn=0;
	

	for (int layer = 0; layer < 8; layer++) {

		std::cout<< "\n\nLayer " << layer << "\n-------------------------------------";

		TPad* pad = (TPad*) C->cd(layer + 1);  // layer+1 because ROOT pads are 1-indexed
		pad->SetName(Form("Layer_%d", layer));  // Properly sets a name for the pad
		C->cd(layer + 1);

		TH2I* dummy = new TH2I(Form("dummy_layer%d_%d", layer, canvas_id),
		                       Form("Layer %d;X (m);Y (m)", layer),
		                       100, -0.5, 0.5, 100, -1.0, 1.0);
		dummy->SetStats(0);
		dummy->Draw();

		auto it = StripsPerLayer.find(layer);
		if (it == StripsPerLayer.end()) continue;

		double mod_x = 0.0;
		double mod_y = 0.0;



		
		double u_angle = u_angles[layer];
		double v_angle = v_angles[layer];
		double pitch = 0.004;
		int n_u = nstrips_u[layer];
		int n_v = nstrips_v[layer];



		for (int u : it->second.first)
			drawGEMStrip(mod_x, mod_y, "U", u_angle, u, pitch, n_u);
			// uDrawn+=1;

		for (int v : it->second.second)
			drawGEMStrip(mod_x, mod_y, "V", v_angle, v, pitch, n_v);
			// vDrawn+=1;


		std::cout<< "\nActive Ustrips = " << uDrawn << " Active Ustrips = " << vDrawn << std::endl;

		uDrawn=0; vDrawn=0;


		TLatex* latex = new TLatex();
		latex->SetNDC(); // normalized coordinates (0–1)
		latex->SetTextSize(0.05);
		latex->DrawLatex(0.1, 0.9, Form("Layer %d", layer));
			

		delete dummy;
	}

	std::cout<< "\nActive Ustrips = " << uDrawn << " Active Ustrips = " << vDrawn << std::endl;

	C->Update();

}



std::map< int, GEMLayerROItoStrips > map_GEMLayerROItoStrips;


// std::map< int, GEMLayerROItoStrips > Get_BinLayerStrips_Map( const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt" ){
void Get_BinLayerStrips_Map( const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt" ){

	DBread db{ db_local };

	ROIread roi{ roi_file };

	if ( db.returnFileReadStatus() == -1 || roi.returnFileReadStatus() == -1 )
	{
		std::cerr << "Exiting the program!!!" << std::endl << std::endl;

		// return -1;
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


		
	std::map< int, std::map< int, GEMLayer > > map_ROIsStripsByBinsAndLayers;
	
	}

}


std::map< int, GEMLayerROItoStrips > Get_LayerStrips_Map( const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt" )
{

	DBread db{ db_local };

	ROIread roi{ roi_file };

	if ( db.returnFileReadStatus() == -1 || roi.returnFileReadStatus() == -1 )
	{
		std::cerr << "Exiting the program!!!" << std::endl << std::endl;

		// return -1;
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

	return map_GEMLayerROItoStrips;
}




//Main func
int showGEMstripsHit_for_ecalbin(const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt") {
	DBread db{ db_local };
	ROIread roi{ roi_file };

	if (db.returnFileReadStatus() == -1 || roi.returnFileReadStatus() == -1) {
		std::cerr << "Exiting the program!!!" << std::endl;
		return -1;
	}

	std::map<int, GEMLayer> map_GEMLayers = db.returnGEMLayerMap();
	std::map<int, std::map<int, ROI>> map_ROIsByBinsAndLayers = roi.return_ROIMap();

	std::string usrinput;
	std::cout << "Enter ECal bin number to visualize (Enter for all bins, or 'q' to quit): ";
	std::getline(std::cin, usrinput);

	if (usrinput == "q") {
		std::cout << "Quitting visualization.\n";
		return 0;
	}

	bool single_bin_mode = false;
	int chosen_bin = -1;
	if (!usrinput.empty()) {
		try {
			chosen_bin = std::stoi(usrinput);
			single_bin_mode = true;
		} catch (const std::invalid_argument&) {
			std::cerr << "Invalid input. Please enter a valid number or leave blank.\n";
			return -1;
		}
	}

	for (const auto& [bin, layerROIs] : map_ROIsByBinsAndLayers) {
    std::cout << "\n\nBin Number: " << bin << std::endl;

    std::map<int, std::pair<std::set<int>, std::set<int>>> currBinStripsPerLayer;

    for (const auto& [modNum, roi] : layerROIs) {
        if (!map_GEMLayers.count(modNum)) continue;

        // GEMLayer thisGEMLayer = map_GEMLayers[modNum];
        // GEMLayerROItoStrips roiTool(thisGEMLayer);

        // auto modUVstripMap = roiTool.takeROI_givePhysicalUVStrips(roi);

        // Determine GEM layer index from module number
        int layer = -1;
        if (modNum <= 5) layer = modNum;
        else if (modNum <= 9) layer = 6;
        else if (modNum <= 13) layer = 7;


		GEMLayer thisGEMLayer = map_GEMLayers[layer];
        GEMLayerROItoStrips roiTool(thisGEMLayer);

        auto modUVstripMap = roiTool.takeROI_givePhysicalUVStrips(roi);

        if (layer < 0) continue;

        for (const auto& [m, uvpair] : modUVstripMap) {
            currBinStripsPerLayer[layer].first.insert(uvpair.first.begin(), uvpair.first.end());
            currBinStripsPerLayer[layer].second.insert(uvpair.second.begin(), uvpair.second.end());
        }
    }

    for (const auto& [layer, uvpair] : currBinStripsPerLayer) {
        std::cout << "Layer " << layer << "\n" << std::string(40, '-') << "\n"
                  << "U: " << uvpair.first.size() << " | V: " << uvpair.second.size() << std::endl;
    }

    if (!single_bin_mode || bin == chosen_bin) {
        std::cout << "\nCalling showgemstrips_for_ecalbin()\n";
        showgemstrips_for_ecalbin(currBinStripsPerLayer);
        if (single_bin_mode) break;
    }
}


	return 0;
}




// 		std::string usrinput;
// 		std::cout << "Enter ECal bin number to visualize (Enter for all bins, or 'q' to quit): ";
// 		std::getline(std::cin, usrinput);

// 		if (usrinput == "q") {
// 			std::cout << "Quitting visualization.\n";
// 			return 0;
// 		}

// 		bool single_bin_mode = false;
// 		int chosen_bin = -1;
// 		if (!usrinput.empty()) {
// 			try {
// 				chosen_bin = std::stoi(usrinput);
// 				single_bin_mode = true;

// 			} catch (const std::invalid_argument&) {
// 				std::cerr << "Invalid input. Please enter a valid number or leave blank.\n";
// 				return -1;
// 			}
// 		}



// 	std::map<int, GEMLayerROItoStrips> map_GEMLayerROItoStrips= Get_LayerStrips_Map(db_local, roi_file);



// 	for (auto&[bin, layer] : map_GEMLayerROItoStrips){

// 		std::cout << "\n\nbin Number: " << bin << std::endl;

// 		if(single_bin_mode==true){ 
// 			if (bin==chosen_bin){
// 				showgemstrips_for_ecalbin(layer.GetUVpairs());
// 			}
// 		}
// 		else{
// 			showgemstrips_for_ecalbin(layer.GetUVpairs());
// 		}
// 	}

// 	return 0;
// }




// for (auto& [binNum, map_ROIbyLayer_forThisBin] : map_ROIsByBinsAndLayers) {
// 	std::map<int, std::pair<std::set<int>, std::set<int>>> thisLayer_ROI_Strips;

// 	std::cout << "Processing ECal bin: " << binNum << std::endl;

// 	for (auto& [modNum, gemLayerROItoStripsInstance] : map_GEMLayerROItoStrips) {
// 		int layer = -1;
// 		if (modNum >= 0 && modNum <= 5) layer = modNum;
// 		else if (modNum >= 6 && modNum <= 9) layer = 6;
// 		else if (modNum >= 10 && modNum <= 13) layer = 7;
// 		else continue;

// 		const ROI& roi = map_ROIbyLayer_forThisBin.at(modNum);
// 		auto strip_pair = gemLayerROItoStripsInstance.takeROI_givePhysicalUVStrips(roi);

// 		thisLayer_ROI_Strips[layer].first.insert(strip_pair[layer].first.begin(), strip_pair[layer].first.end());
// 		thisLayer_ROI_Strips[layer].second.insert(strip_pair[layer].second.begin(), strip_pair[layer].second.end());
// 	}
	
// 	std::cout << "  Layers with strips for ECal bin " << binNum << ":\n";
// 	for (const auto& [layer, pair] : thisLayer_ROI_Strips) {
// 		std::cout << "    Layer " << layer << " — U: " << pair.first.size()
// 				<< " | V: " << pair.second.size() << "\n";
// 	}

// 	showgemstrips_for_ecalbin(thisLayer_ROI_Strips);  // <---- Add this here

// 	// // LIMIT: For testing only show certain bin
// 	// if (loopNum==1){
// 	// 	break;
// 	// 	showgemstrips_for_ecalbin(thisLayer_ROI_Strips);
// 	// }
// 	// loopNum +=1;

// 	showgemstrips_for_ecalbin(thisLayer_ROI_Strips);
	
// 	// break;

// }
	

// // 	for (auto& [layerNum, gemLayerInstance] : map_GEMLayers) {
// // 		map_GEMLayerROItoStrips.emplace(layerNum, GEMLayerROItoStrips{gemLayerInstance});
// // 	}

// // 	std::cout << "N ECal bins: " << map_ROIsByBinsAndLayers.size() << std::endl;
// // 	std::cout << "N GEM layers: " << map_GEMLayers.size() << std::endl;


// // 	std::string usrinput;
// // 	std::cout << "Enter ECal bin number to visualize (Enter for all bins, or 'q' to quit): ";
// // 	std::getline(std::cin, usrinput);

// // 	if (usrinput == "q") {
// // 		std::cout << "Quitting visualization.\n";
// // 		return 0;
// // 	}

// // 	bool single_bin_mode = false;
// // 	int chosen_bin = -1;
// // 	if (!usrinput.empty()) {
// // 		try {
// // 			chosen_bin = std::stoi(usrinput);
// // 			single_bin_mode = true;
// // 		} catch (const std::invalid_argument&) {
// // 			std::cerr << "Invalid input. Please enter a valid number or leave blank.\n";
// // 			return -1;
// // 		}
// // 	}



// // 	// for (const auto& [binNum, map_ROIbyLayer_forThisBin] : map_ROIsByBinsAndLayers) {
// // 	int bin_count = 0;
// // 	int loopNum = 0;



// // 	for (auto& [binNum, map_ROIbyLayer_forThisBin] : map_ROIsByBinsAndLayers) {
// // 		std::map<int, std::pair<std::set<int>, std::set<int>>> thisLayer_ROI_Strips;
	
// // 		std::cout << "Processing ECal bin: " << binNum << std::endl;
	
// // 		for (auto& [modNum, gemLayerROItoStripsInstance] : map_GEMLayerROItoStrips) {
// // 			int layer = -1;
// // 			if (modNum >= 0 && modNum <= 5) layer = modNum;
// // 			else if (modNum >= 6 && modNum <= 9) layer = 6;
// // 			else if (modNum >= 10 && modNum <= 13) layer = 7;
// // 			else continue;

// // 			const ROI& roi = map_ROIbyLayer_forThisBin.at(modNum);
// // 			auto strip_pair = gemLayerROItoStripsInstance.takeROI_givePhysicalUVStrips(roi);

// // 			thisLayer_ROI_Strips[layer].first.insert(strip_pair[layer].first.begin(), strip_pair[layer].first.end());
// // 			thisLayer_ROI_Strips[layer].second.insert(strip_pair[layer].second.begin(), strip_pair[layer].second.end());
// // 		}
		
// // 		std::cout << "  Layers with strips for ECal bin " << binNum << ":\n";
// // 		for (const auto& [layer, pair] : thisLayer_ROI_Strips) {
// // 			std::cout << "    Layer " << layer << " — U: " << pair.first.size()
// // 					<< " | V: " << pair.second.size() << "\n";
// // 		}

// // 		showgemstrips_for_ecalbin(thisLayer_ROI_Strips);  // <---- Add this here

// // 		// // LIMIT: For testing only show certain bin
// // 		// if (loopNum==1){
// // 		// 	break;
// // 		// 	showgemstrips_for_ecalbin(thisLayer_ROI_Strips);
// // 		// }
// // 		// loopNum +=1;

// // 		showgemstrips_for_ecalbin(thisLayer_ROI_Strips);
		
// // 		// break;
	
// // 	}

// // 	return 0;
// // }
