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


#include "GetGemInfoMap.cxx"




std::map<int, gemInfo> refModMap = GetGemInfoMap();

// // ---- Geometry settings ----
// std::map<int, double> u_angles = {
// 	{0, 0}, {1, 0},
// 	{2, 60}, {3, 60}, {4, 60}, {5, 60},
// 	{6, 0}, {7, 0}
// };

// std::map<int, double> v_angles = {
// 	{0, 0}, {1, 0},
// 	{2, -60}, {3, -60}, {4, -60}, {5, -60},
// 	{6, 90}, {7, 90}
// };

// std::map<int, int> nstrips_u = {
// 	{0, 3968}, {1, 3968},
// 	{2, 3840}, {3, 3840}, {4, 3840}, {5, 3840},
// 	{6, 5120}, {7, 5120}
// };

// std::map<int, int> nstrips_v = {
// 	{0, 3456}, {1, 3456},
// 	{2, 3840}, {3, 3840}, {4, 3840}, {5, 3840},
// 	{6, 1536}, {7, 1536}
// };


int uDrawn=0, vDrawn=0;


// ---- Strip drawing function ----
// void drawGEMStrip(double mod_x, double mod_y, TString axis, double angle_deg,int strip_num, double pitch = 0.004, int n_strips = 256, double length = 0.8)
// void drawGEMStrip(double mod_x, double mod_y, TString axis, double angle_deg, int strip_num, double pitch = 0.004, int n_strips = 1, double mod_size_x = 0.8, double mod_size_y = 0.4) // or no defaults
void drawGEMStrip(double mod_x, double mod_y, TString axis, double angle_deg,
	int strip_num, double pitch, int n_strips,
	double mod_size_x, double mod_size_y)
{
double angle = angle_deg * TMath::DegToRad();
int center_strip = n_strips / 2;
double offset = (strip_num - center_strip) * pitch;

// Full-length endpoints before clipping
double dx = offset * TMath::Cos(angle + TMath::Pi()/2);
double dy = offset * TMath::Sin(angle + TMath::Pi()/2);
double x_center = mod_x + dx;
double y_center = mod_y + dy;
double half_len = 0.5 * mod_size_y * 1.2; // generous initial length

double x1 = x_center - half_len * TMath::Cos(angle);
double y1 = y_center - half_len * TMath::Sin(angle);
double x2 = x_center + half_len * TMath::Cos(angle);
double y2 = y_center + half_len * TMath::Sin(angle);

// Box edges
double xmin = mod_x - mod_size_x / 2.0;
double xmax = mod_x + mod_size_x / 2.0;
double ymin = mod_y - mod_size_y / 2.0;
double ymax = mod_y + mod_size_y / 2.0;

// Lambda to clip a point to the rectangle
auto inside = [&](double x, double y) {
return x >= xmin && x <= xmax && y >= ymin && y <= ymax;
};

// Liang-Barsky like clipping
auto clipLine = [&](double& x1, double& y1, double& x2, double& y2) -> bool {
double dx = x2 - x1;
double dy = y2 - y1;

double p[4] = {-dx, dx, -dy, dy};
double q[4] = {x1 - xmin, xmax - x1, y1 - ymin, ymax - y1};

double u1 = 0, u2 = 1;

for (int i = 0; i < 4; i++) {
if (p[i] == 0 && q[i] < 0) return false;
if (p[i] != 0) {
  double u = q[i] / p[i];
  if (p[i] < 0) u1 = std::max(u1, u);
  else         u2 = std::min(u2, u);
}
}

if (u1 > u2) return false;

double new_x1 = x1 + u1 * dx;
double new_y1 = y1 + u1 * dy;
double new_x2 = x1 + u2 * dx;
double new_y2 = y1 + u2 * dy;

x1 = new_x1; y1 = new_y1;
x2 = new_x2; y2 = new_y2;

return true;
};

if (!clipLine(x1, y1, x2, y2)) return; // Skip drawing if line is totally outside

// Draw center marker
TMarker *thisHit = new TMarker(x_center, y_center, 20);
thisHit->SetMarkerSize(0.3);
thisHit->SetMarkerColor(kGray + 2);
thisHit->Draw("same");

// Draw clipped line
TLine *stripLine = new TLine(x1, y1, x2, y2);
if (axis.EqualTo("U", TString::kIgnoreCase)) {
stripLine->SetLineColor(kCyan);
stripLine->SetLineWidth(2);
stripLine->SetLineStyle(3);
uDrawn++;
} else {
stripLine->SetLineColor(kPink);
stripLine->SetLineWidth(1);
stripLine->SetLineStyle(2);
vDrawn++;
}
stripLine->Draw("same");
}


// ---- GEM Strip ROI Visualization ----
// void showgemstrips_for_ecalbin(const std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerLayer)
// void showgemstrips_for_ecalbin(std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerLayer)
void showgemstrips_for_ecalbin(const std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerMod)
{

	std::cout << "\n\nStarting showgemstrips_for_ecalbin()" <<std::endl;

	std::cout <<"number of keys = "<< StripsPerMod.size() << std::endl;

	std::cout << "Keys include: "<<std::endl;
	for (auto& [key, val] : StripsPerMod){ std::cout << key << ", "<<std::endl;}

	static int canvas_id = 0;
	TCanvas* C = new TCanvas(Form("canvas_%d", canvas_id++), "Active ROI Strips in GEM Planes", 1600, 1000);
	C->Divide(4, 2);


	// int uDrawn=0, vDrawn=0;
	

	for (int layer = 0; layer < 8; layer++) {

		std::cout << "\n making subPad for layer " << layer <<std::endl;

		// TPad* pad = (TPad*) C->cd(layer + 1);  // layer+1 because ROOT pads are 1-indexed
		// pad->SetName(Form("Layer_%d", layer));  // Properly sets a name for the pad
		C->cd(layer + 1);


		// std::array<double, 3> layerCenter {0., 0., 0.};

		std::array<double,3> modSize = gemInfoMap[layer].size;
		
		std::cout << "\nMod Size: " << modSize[0] << ", " << modSize[1] << ", " << modSize[2] << std::endl;
		
		double layerSizeX;

		double xMin, xMax;//really the -y direction for module coords 
		double yMin, yMax;

		if(layer<6){ layerSizeX = modSize[0];}
		if (layer>5){layerSizeX = modSize[0] *4;}

		std::cout << "layerSizeX = " << layerSizeX << std::endl;

		
		xMin = -modSize[1]/2.; xMax = modSize[1]/2.;
		yMin = -layerSizeX/2.; yMax = layerSizeX/2. ;// first size Dim so it plots the x vals down corresponding to mod coords 
		
		std::cout << "\nxMin, yMin: " << xMin << ", " << yMin 
		<< "\nxMax, yMax: " << xMax << ", " << yMax << std::endl;



		// TODO:: !!!!! IF LAYER = .... DR%AW ... BIG!


		// TODO:: !!! x and y switch b/c coordinates?
		// 	do they match UV strips



		TH2I* dummy = new TH2I(
			Form("dummy_layer%d_%d", layer, canvas_id),
			Form("Layer %d;X (m);Y (m)", layer),
			// 100, -0.5, 0.5, 100, -1.0, 1.0
			// 100, -0.5*100, 0.5*100, 100, -1.0*100, 1.0*100
			100, xMin, xMax, 100, yMin, yMax
		);
		dummy->SetStats(0);
		dummy->Draw();

		// TLatex* latex = new TLatex();
		// latex->SetNDC(); // normalized coordinates (0–1)
		// latex->SetTextSize(0.05);
		// latex->DrawLatex(0.1, 0.9, Form("Layer %d", layer));

	}

	for (auto&[mod, uvStripPair] : StripsPerMod){
		// int mod = 0; mod < gemInfoMap.size(); mod++) {
	// for (int mod = 0; mod < gemInfoMap.size(); mod++) {

		std::cout<< "\nmod " << mod << "\n-------------------------------------";

		// TPad* pad = (TPad*) C->cd(mod + 1);  // layer+1 because ROOT pads are 1-indexed
		// pad->SetName(Form("Layer_%d", mod));  // Properly sets a name for the pad
		// C->cd(mod + 1);

		// TH2I* dummy = new TH2I(
		// 	Form("dummy_mod%d_%d", mod, canvas_id),
		// 	Form("Mod %d;X (m);Y (m)", mod),
		// 	100, -0.5, 0.5, 100, -1.0, 1.0
		// );
		// dummy->SetStats(0);
		// dummy->Draw();


		int layer = -1;
		if (mod <= 5) layer = mod;
		else if (mod <= 9) layer = 6;
		else if (mod <= 13) layer = 7;
		if (layer < 0) continue;

		std::cout << "\n Layer is " << layer << std::endl;


		TPad* pad = (TPad*) C->cd(layer + 1);//so can name layer 1, 2...
		std::cout << "\nChanging to subpad(layer plot) Named " << pad->GetName() <<std::endl;
		std::cout << "\nChanging to subpad(layer plot) Numbererd " << pad->GetNumber() <<std::endl;
		std::cout << "\nChanging to Canvas " << pad->GetCanvasID() <<std::endl;

		pad->cd();  



		// auto it = StripsPerMod.find(mod);
		
		// if (it == StripsPerMod.end()) {
		// 	// std::cout << " --> Skipping mod " << mod << ", not in StripsPerMod\n";
		// 	std::cout << "\n !!!!--> StripsPerMod.end() at mod " << mod << std::endl;

		// 	continue;
		// }

		// if (it->second.first.empty() && it->second.second.empty()) {
		// 	std::cout << "\n !!!!--> Skipping mod " << mod << ", no active U or V strips\n" << std::endl;

		// 	continue;
		// }

		// std::cout << "\nShould not see this message if mod " << mod << " was skipped";

		// // double mod_x = 0.0;
		// // double mod_y = 0.0;

		auto thisModInfo = gemInfoMap[mod];
		
		double mod_x = thisModInfo.position[0];		
		
		double mod_y = thisModInfo.position[1];		


		
		double u_angle = thisModInfo.uvangles[0];
		double v_angle = thisModInfo.uvangles[1];
		double pitch = 0.004;
		int n_u = thisModInfo.nstripsuv[0];
		int n_v = thisModInfo.nstripsuv[1];



		double mod_size_x = thisModInfo.size[0];
		double mod_size_y = thisModInfo.size[1];

		for (auto& uStrip : uvStripPair.first){
			drawGEMStrip(mod_x, mod_y, "U", u_angle, uStrip, pitch, n_u, mod_size_x, mod_size_y);
			// drawGEMStrip(mod_x, mod_y, "U", u_angle, uStrip, pitch, n_u);

		// for (int u : it->second.first){
			// if (!it->second.first.empty()){
			// {	std::cout << "\nDrawing U striups for mod " << mod <<std::endl;
				// drawGEMStrip(mod_x, mod_y, "U", u_angle, u, pitch, n_u);
				// uDrawn+=1;
			// }
		}


		for (auto& vStrip : uvStripPair.second){
			drawGEMStrip(mod_x, mod_y, "V", v_angle, vStrip, pitch, n_v, mod_size_x, mod_size_y);
// drawGEMStrip(mod_x, mod_y, "V", v_angle, vStrip, pitch, n_v);
		// for (int v : it->second.second){
			// if (!it->second.second.empty()){
			// // {std::cout << "\nDrawing V striups for mod " << mod <<std::endl;
			// 	drawGEMStrip(mod_x, mod_y, "V", v_angle, v, pitch, n_v);
			// // vDrawn+=1;
			// }
		}


		std::cout<< "\nActive Ustrips = " << uDrawn << " Active Vstrips = " << vDrawn << std::endl;

		uDrawn=0; vDrawn=0;


		// TLatex* latex = new TLatex();
		// latex->SetNDC(); // normalized coordinates (0–1)
		// latex->SetTextSize(0.05);
		// latex->DrawLatex(0.1, 0.9, Form("Layer %d", layer));
	}

	// std::cout<< "\nActive Ustrips = " << uDrawn << " Active Ustrips = " << vDrawn << std::endl;

	C->Update();

}





//Main func
int showGEMstripsHit_for_ecalbin(const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt") {
	
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


	for ( auto& [binNum, map_ROIbyLayer_forThisBin] : map_ROIsByBinsAndLayers )
	{
		std::map <int, std::pair < std::set<int>, std::set<int> >> map_allUVstripsSetsForAllModules_forThisECalBin;

		//std::cout << "***Filling bin: " << binNum << std::endl;
		
		
		

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
		
			}

		map_physicalUVStrips_byECalBin_byGEMMod[binNum] = map_allUVstripsSetsForAllModules_forThisECalBin;

		}
	}

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

	

	for ( const auto& [binNum, uvStripSetbyModule] : map_physicalUVStrips_byECalBin_byGEMMod )
	{
		std::cout << "\n### BIN NUMBER: " << binNum << " ###" << std::endl;

		if (!single_bin_mode || binNum == chosen_bin) {
			std::cout << "\nCalling showgemstrips_for_ecalbin() from main script\n";
			showgemstrips_for_ecalbin(uvStripSetbyModule);
		}
		if (single_bin_mode) { break;}


		// showgemstrips_for_ecalbin(uvStripSetbyModule);


		// for ( const auto& [modNUm, uvStripPair] : uvStripSetbyModule )

		
		// {
		// 	std::cout << "\n*** Mod Num: " << modNUm << " ***" << std::endl;

	 	// 	for ( const auto& ustrip : uvStripPair.first ) {

		// 	}

		// 	for ( const auto& vstrip : uvStripPair.second ) {

		// 	}
		// }
	}

	


	return 0;

}