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

#include <fstream>

#include "GEMModROItoStrips.h"
#include "DBread.h"
#include "ROIread.h"
#include "showgemhit_for_ecalbin.C"


#include "GetGemInfoMap.cxx"

// namespace GEMstripsHit_for_ecalbin{

std::map<int, gemInfo> refModMap = GetGemInfoMap();

int global_canvas_id = 0;


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


void drawGEMStrip(int modNum, int axis, int strip_num, int center_strip)
{
	if(strip_num%10==0){//Testing with greater	plotted spacing for visibility
	auto& modInfo = refModMap[modNum];

	double mod_x = modInfo.position[0];
	double mod_y = modInfo.position[1];
	double mod_size_x = modInfo.size[0];
	double mod_size_y = modInfo.size[1];
	double angle = modInfo.uvangles[axis] * TMath::DegToRad();
	int n_strips = modInfo.nstripsuv[axis];
	double pitch = 0.0004;


	double layerSizeX=0;
	double layerSizeY=0;

	// if(modNum<6){ layerSizeX = mod_size_x;}
	// if (modNum>=6){layerSizeX = mod_size_x *4;}



	// std::cout << "\nnumber of strips in mod " << n_strips << std::endl;
	// std::cout << "\nThis strip number " << strip_num << std::endl;


	// Offset from center strip
	// int center_strip = n_strips / 2;




	
	// double offset = (strip_num - center_strip) * pitch;//from center
	double offset = (center_strip-strip_num) * pitch;//b/c strips counted downwards?



	//TODO: consider specifics of y case


	if(modNum<6){ layerSizeX = mod_size_x;
		// offset-=.0108;//see .dat file
		}
		if (modNum>=6){layerSizeX = mod_size_x *4;}

		double dx_offset, dy_offset;

	// Offset is applied perpendicular to the strip direction
	// double dx_offset = std::sin(angle);
	// double dy_offset = -std::cos(angle);

	// double dx_offset = std::cos(angle - TMath::Pi()/2);
	// double dy_offset = std::sin(angle - TMath::Pi()/2);
	// double dx_offset = std::cos(angle + TMath::Pi()/2);
	// double dy_offset = std::sin(angle + TMath::Pi()/2);

	if (axis == 0) {
		// from 0 to pi/2 (U strips)
		dx_offset = std::cos(angle + TMath::Pi()/2);
		dy_offset = std::sin(angle + TMath::Pi()/2);
	} else {
		//from 0 to 3pi/2
		dx_offset = std::cos(angle - TMath::Pi()/2);
		dy_offset = std::sin(angle - TMath::Pi()/2);
	}


// 	Double_t fPxU = cos(uangle);            //U Strip X projection = cos( UAngle );
//   Double_t fPyU = sin(uangle);            //U Strip Y projection = sin( UAngle );
//   Double_t fPxV = cos(vangle);            //V Strip X projection = cos( VAngle );
//   Double_t fPyV = sin(vangle);            //V Strip Y projection = sin( VAngle );



	// double zero_strip_offset = -center_strip * pitch;
	// double strip_relative_offset = strip_num * pitch;
	
	// double total_offset = zero_strip_offset + strip_relative_offset;
	
	// double x_center = mod_x + total_offset * dx_offset;
	// double y_center = mod_y + total_offset * dy_offset;
	double x_center = mod_x + offset * dx_offset ;
	double y_center = mod_y + offset * dy_offset;
	// Compute starting edge of strip 0 (first strip)
	
	

	
	// double half_len = 0.5 * mod_size_y;
	// Give a generous length for the line, then clip
	// double half_len = 0.5 * mod_size_y * 1.2;
	double half_length = 0.5 * (
		fabs(mod_size_x * cos(angle)) +
		fabs(mod_size_y * sin(angle))
	);
	

	double x1 = x_center - mod_size_x/2 * std::cos(angle);
	// double x1 = x_center - half_length * std::cos(angle);
	// double x1 = x_center - layerSizeX/2 * std::cos(angle);
	double y1 = y_center - mod_size_y/2 * std::sin(angle);
	// double y1 = y_center - half_length * std::sin(angle);
	double x2 = x_center + mod_size_x/2 * std::cos(angle);
	// double x2 = x_center + half_length * std::cos(angle);
	// double x2 = x_center + layerSizeX/2 * std::cos(angle);
	double y2 = y_center + mod_size_y/2 * std::sin(angle);
	// double y2 = y_center + half_length * std::sin(angle);
	


	if (modNum < 3){//5-7 are kinda fine focus on way off first

		std::cout << "\n------------\nMod Num " << modNum << std::endl;
		std::cout << "axis: " << axis << std::endl;
		std::cout << "angle: " << angle*TMath::RadToDeg() << std::endl;

		std::cout << "\nnumber of strips in mod " << n_strips << std::endl;
		std::cout << "Center strip: " << center_strip << std::endl;
		std::cout << "This strip number " << strip_num << std::endl;
		std::cout << " Strip offset: " << strip_num-center_strip << std::endl;

		std::cout << "\noffset: " << offset << std::endl;
		std::cout << "dx offset: " << dx_offset << std::endl;
		std::cout << "dy offset: " << dy_offset << std::endl;
		std::cout << "offset * dx_offset: " << offset * dx_offset << std::endl;
		std::cout << "offset * dy_offset: " << offset * dy_offset << std::endl;
	
		std::cout << "\nmod Pos size(modCoords) x y: " << mod_x << " " << mod_y<<std::endl;
		std::cout << "layer size(modCoords) x y: " << layerSizeX << " " << mod_size_y<<std::endl;
	
		std::cout << "\nxCenter yCenter: " << x_center << " " << y_center <<std::endl;
	};





		// Optional: mark center point (useful for debug)
		TMarker* centerMark = new TMarker(y_center, x_center, 20);
		centerMark->SetMarkerSize(0.3);
		// centerMark->SetMarkerColor(kGray + 2);
		centerMark->Draw("same");


		
// 	// 	std::cout << "x1 y1: " << x1 << " " << y1 <<std::endl;
// 	// 	std::cout << "x2 y2: " << x2 << " " << y2 <<std::endl;
// 	// 	}

// 		double layer_center_x, layer_center_y;

// if (modNum < 6) {
//     layer_center_x = mod_x;
//     layer_center_y = mod_y;
// } else {
//     std::vector<int> layerMods = (modNum <= 9) ? std::vector<int>{6,7,8,9} : std::vector<int>{10,11,12,13};
//     double xsum = 0.0, ysum = 0.0;
//     for (int m : layerMods) {
//         xsum += refModMap[m].position[0];
//         ysum += refModMap[m].position[1];
//     }
//     layer_center_x = xsum / layerMods.size();
//     layer_center_y = ysum / layerMods.size();
// }

// 		x_center -= layer_center_x;
// 		y_center -= layer_center_y;
// 		x1 -= layer_center_x;
// 		y1 -= layer_center_y;
// 		x2 -= layer_center_x;
// 		y2 -= layer_center_y;




// 		// Optional: mark center point (useful for debug)
// 		TMarker* centerMark = new TMarker(y_center, x_center, 20);
// 		centerMark->SetMarkerSize(0.3);
// 		// centerMark->SetMarkerColor(kGray + 2);
// 		centerMark->Draw("same");


	
	// TLine* xCenterLine = new TLine(0, -layerSizeX/2, 0, +layerSizeX/2);
	TLine* xCenterLine = new TLine(0, -mod_size_x/2, 0, +mod_size_x/2);
	TLine* yCenterLine = new TLine(mod_y-mod_size_y/2, 0, mod_y+mod_size_y/2, 0);
	// TLine* xCenterLine = new TLine(0, mod_y-mod_size_y/2, 0, mod_y+mod_size_y/2);
	// TLine* yCenterLine = new TLine(-layerSizeX/2, 0, +layerSizeX/2, 0);


	xCenterLine->Draw("SAME");
	yCenterLine->Draw("SAME");




	// Module boundaries in global coordinates

	

	double xmin = - layerSizeX / 2.0;
	double xmax = + layerSizeX / 2.0;

	double ymin = -mod_size_y / 2.0;
	double ymax = mod_size_y / 2.0;

	// std::cout << "\nclipping to xmin ymin xmax ymax" << "\n"<< xmin << " "<< ymin << " "<< xmax << " "<< ymax <<std::endl;

	// Liang-Barsky clipping
	auto clipLine = [&](double& x1, double& y1, double& x2, double& y2) -> bool {
		double dy = y2 - y1;
		double dx = x2 - x1;

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

	if (!clipLine(x1, y1, x2, y2)) return; // Strip is entirely outside module


	
	// Optional: mark center point (useful for debug)
	// TMarker* centerMark = new TMarker(x_center, y_center, 20);
	// centerMark->SetMarkerSize(0.3);
	// centerMark->SetMarkerColor(kGray + 2);
	// centerMark->Draw("same");

	// Draw the clipped strip
	TLine* stripLine = new TLine(y1, x1, y2, x2); //Drawn 
	if (axis == 0) { // U
		stripLine->SetLineColor(kCyan);
		// stripLine->SetLineStyle(2);
		// stripLine->SetLineStyle(9);
		stripLine->SetLineWidth(1);
		// stripLine->SetLineWidth(.3);
		uDrawn++;
	} else {         // V
		stripLine->SetLineColor(kPink);
		// stripLine->SetLineStyle(2);
		// stripLine->SetLineStyle(9);
		stripLine->SetLineWidth(1);
		// stripLine->SetLineWidth(.3);
		vDrawn++;
	}
	stripLine->Draw("same");
}
}



// ---- GEM Strip ROI Visualization ----
// void showgemstrips_for_ecalbin(const std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerLayer)
// void showgemstrips_for_ecalbin(std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerLayer)
// void showgemstrips_for_ecalbin(int ecalBinNum, const std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerMod)
TCanvas* showgemstrips_for_ecalbin(int ecalBinNum, const std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerMod)
{

	std::cout << "\n\nStarting showgemstrips_for_ecalbin()" <<std::endl;

	std::cout <<"number of keys = "<< StripsPerMod.size() << std::endl;

	std::cout << "Keys include: "<<std::endl;
	for (auto& [key, val] : StripsPerMod){ 
		std::cout << key << ": "
		 << "\nnumber of U strips " << val.first.size() << ", number of V strips " << val.second.size()<<std::endl;
	}



	// static int canvas_id = 0;
	int canvas_id = 0;
	// TCanvas* C = new TCanvas(Form("canvas_%d", canvas_id++), "Active ROI Strips in GEM Planes", 1600, 1000);
	TCanvas* C = new TCanvas(Form("canvas_%d", global_canvas_id++), "Active ROI Strips in GEM Planes", 1600, 1000);

	
	
	C->SetTitle("Bin ");

	C->Divide(4, 2);

	// Draw ECal Bin number at the top of the canvas
TLatex latex;
latex.SetNDC(); // Use normalized device coordinates
latex.SetTextSize(0.04);
latex.SetTextAlign(22); // Center alignment
latex.DrawLatex(0.5, 0.97, Form("ECal Bin %d", ecalBinNum));


	// int uDrawn=0, vDrawn=0;
	

	for (int layer = 0; layer < 8; layer++) {

		std::cout << "\n making subPad for layer " << layer 
		<< "\n--------------------------------------"
		<<std::endl;

		// TPad* pad = (TPad*) C->cd(layer + 1);  // layer+1 because ROOT pads are 1-indexed
		// pad->SetName(Form("Layer_%d", layer));  // Properly sets a name for the pad
		C->cd(layer + 1);


		// std::array<double, 3> layerCenter {0., 0., 0.};

		std::array<double,3> modSize = refModMap[layer].size;
		std::array<double,3> modPos = refModMap[layer].position;
		
		
		std::cout << "\nMod Pos: " << modPos[0] << ", " << modPos[1] << ", " << modPos[2] << std::endl;

		std::cout << "\nMod Size: " << modSize[0] << ", " << modSize[1] << ", " << modSize[2] << std::endl;
		
		
		double layerSizeX;

		double xMin, xMax;//really the -y direction for module coords 
		double yMin, yMax;

		if(layer<6){ layerSizeX = modSize[0];}
		if (layer>5){layerSizeX = modSize[0] *4;}

		std::cout << "layerSizeX = " << layerSizeX << std::endl;

		
		// xMin = -modSize[1]/2.; xMax = modSize[1]/2.;
		// yMin = -layerSizeX/2.; yMax = layerSizeX/2. ;// first size Dim so it plots the x vals down corresponding to mod coords 

		double mod_size_x = modSize[0]; // actual horizontal extent
		double mod_size_y = modSize[1]; // actual vertical extent

		xMin = -layerSizeX / 2.;// *3;
		xMax = +layerSizeX / 2.;// *3;
		yMin = -mod_size_y / 2.;// *3;
		yMax = +mod_size_y / 2.;// *3;

		
		std::cout << "\nxMin, yMin: " << xMin << ", " << yMin 
		<< "\nxMax, yMax: " << xMax << ", " << yMax << std::endl;



		// TODO:: !!!!! IF LAYER = .... DR%AW ... BIG!


		// TODO:: !!! x and y switch b/c coordinates?
		// 	do they match UV strips



		TH2I* dummy = new TH2I(
			Form("dummy_layer%d_%d", layer, canvas_id),
			Form("Layer %d;Y (m);X (m)", layer),
			// 100, -0.5, 0.5, 100, -1.0, 1.0
			// 100, -0.5*100, 0.5*100, 100, -1.0*100, 1.0*100
			100, yMin, yMax, 100, xMin, xMax
		);
		dummy->SetStats(0);
		dummy->Draw();

		// TLatex* latex = new TLatex();
		// latex->SetNDC(); // normalized coordinates (0–1)
		// latex->SetTextSize(0.05);
		// latex->DrawLatex(0.1, 0.9, Form("Layer %d", layer));

	}

	for (auto&[mod, uvStripPair] : StripsPerMod){
		// int mod = 0; mod < refModMap.size(); mod++) {
	// for (int mod = 0; mod < refModMap.size(); mod++) {

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




		// DEAL WITH IF NAN



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

		// // double mod_y = 0.0;

		auto thisModInfo = refModMap[mod];
		
		double mod_x = thisModInfo.position[0];		
		
		double mod_y = thisModInfo.position[1];		


		
		double u_angle = thisModInfo.uvangles[0];
		double v_angle = thisModInfo.uvangles[1];
		double pitch = 0.004;
		int n_u = thisModInfo.nstripsuv[0];
		int n_v = thisModInfo.nstripsuv[1];



		double mod_size_x = thisModInfo.size[0];
		double mod_size_y = thisModInfo.size[1];


		//NOTE: Plot V first as cyan(v) under pink(u) is less visible/readable

		std::cout << "\nU,V angles are " << u_angle <<", "<< v_angle <<std::endl;

		std::cout << "\n N uStrips = " << uvStripPair.first.size() << std::endl;
		std::cout << "N vStrips = " << uvStripPair.second.size() << std::endl;

		int vStripLoopCalls=0;
		for (auto& vStrip : uvStripPair.second){
		// for (auto& vStrip : std::array{10, 20 , 30 , 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000}){
			// std::cout << "\n number of V strips " << uvStripPair.second.size()<<std::endl;
			// std::cout << "\nDrawing vStrip " << vStrip << std::endl;

			vStripLoopCalls+=1;
			int local_center_v = (*uvStripPair.second.begin() + *uvStripPair.second.rbegin()) / 2;
			drawGEMStrip(mod, 1, vStrip, local_center_v);
			// drawGEMStrip(mod_x, mod_y, "V", v_angle, vStrip, pitch, n_v, mod_size_x, mod_size_y);
// drawGEMStrip(mod_x, mod_y, "V", v_angle, vStrip, pitch, n_v);
		// for (int v : it->second.second){
			// if (!it->second.second.empty()){
			// // {std::cout << "\nDrawing V striups for mod " << mod <<std::endl;
			// 	drawGEMStrip(mod_x, mod_y, "V", v_angle, v, pitch, n_v);
			// // vDrawn+=1;
			// }.q

		}

		int uStripLoopCalls=0;
		
		// for (int uStrip=0 , uStrip=uStripuvStripPair.first.size; uStrip++){¸
		for (auto& uStrip : uvStripPair.first){
		// for (auto& uStrip : std::array{10, 20 , 30 , 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000}){

			// std::cout << "\nnumber of U strips " << uvStripPair.first.size() <<std::endl;
			// std::cout << "\nDrawing uStrip" << uStrip << std::endl;

			uStripLoopCalls+=1;
			// std::cout << "\ncalled uStrip loop " << uStripLoopCalls << std::endl;
			int local_center_u = (*uvStripPair.first.begin() + *uvStripPair.first.rbegin()) / 2;
			drawGEMStrip(mod, 0, uStrip, local_center_u);
			// drawGEMStrip(mod, 0, uStripuvStripPair.first[0]);
			// drawGEMStrip(mod_x, mod_y, "U", u_angle, uStrip, pitch, n_u, mod_size_x, mod_size_y);
			// drawGEMStrip(mod_x, mod_y, "U", u_angle, uStrip, pitch, n_u);

		// for (int u : it->second.first){
			// if (!it->second.first.empty()){
			// {	std::cout << "\nDrawing U striups for mod " << mod <<std::endl;
				// drawGEMStrip(mod_x, mod_y, "U", u_angle, u, pitch, n_u);
				// uDrawn+=1;
			// }
		}



		std::cout<< "\nActive Ustrips = " << uDrawn << " Active Vstrips = " << vDrawn << std::endl;

		std::cout << "Number times called uLoop vLoop " << uStripLoopCalls << " " << vStripLoopCalls <<std::endl;

		uDrawn=0; vDrawn=0;

		uStripLoopCalls=0;
		vStripLoopCalls=0;


		// TLatex* latex = new TLatex();
		// latex->SetNDC(); // normalized coordinates (0–1)
		// latex->SetTextSize(0.05);
		// latex->DrawLatex(0.1, 0.9, Form("Layer %d", layer));
	}

	for (auto& [mod, uvStripPair] : StripsPerMod) {
		std::cout << "Module " << mod 
				  << ": U strips = " << uvStripPair.first.size() 
				  << ", V strips = " << uvStripPair.second.size() << std::endl;
	}

	// std::cout<< "\nActive Ustrips = " << uDrawn << " Active Ustrips = " << vDrawn << std::endl;

	// C->SetBit(kMustCleanup); // ROOT will clean safely
	C->Update();
	// gSystem->ProcessEvents();
	// std::cout << "Press enter to exit...\n";
	// std::cin.get();
	return C;
}


int exportPDF(TCanvas * C){
	C->SaveAs("GEM_ECalBin.pdf");
	return 0;
}


//Main func
int showGEMstripsHit_for_ecalbin(const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt") {

	// using namespace GEMstripsHit_for_ecalbin;
	
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
	std::cout << "Enter ECal bin numbers to visualize separated by spaces (or 'q' to quit): ";
	std::getline(std::cin, usrinput);

	// Quit option
	if (usrinput == "q" || usrinput == "Q") {
		std::cout << "Quitting visualization.\n";
		return 0;
	}

	// Prepare list of selected bins
	std::vector<int> selected_bins;

	// If user input is empty, we will process ALL bins
	bool all_bins_mode = usrinput.empty();

	if (!all_bins_mode) {
		std::istringstream iss(usrinput);
		int bin;
		while (iss >> bin) {
			selected_bins.push_back(bin);
		}
		if (selected_bins.empty()) {
			std::cerr << "Invalid input. No valid bin numbers detected.\n";
			return -1;
		}
	}

	
	// // --- NEW CLEAN INPUT: Ask about PDF saving ---
	// bool save_as_pdf = false;
	// std::string usr_save_as_pdf;
	// std::cout << "Save as PDF and quit? (y/n): ";
	// std::getline(std::cin, usr_save_as_pdf);
	// if (usr_save_as_pdf == "y" || usr_save_as_pdf == "Y") {
	// 	save_as_pdf = true;
	// }
	
	
	// bool interactive_view = false;
	// std::string usr_interactive_view;
	// std::cout << "Visualize with ROOT interactive window? (y/n): ";
	// std::getline(std::cin, usr_interactive_view);
	// if (usr_interactive_view == "y" || usr_interactive_view == "Y") {
	// 	interactive_view = true;
	// }






	// -- Printing info --
std::cout << "\nNumber of ECalBins is " << map_physicalUVStrips_byECalBin_byGEMMod.size() << std::endl;

	if (all_bins_mode) {
		std::cout << "\n### ALL BINS MODE ###" << std::endl;

		for (const auto& [binNum, uvStripSetbyModule] : map_physicalUVStrips_byECalBin_byGEMMod) {
			std::cout << "\n### BIN NUMBER: " << binNum << " ###" << std::endl;
			std::cout << "Number of Modules " << uvStripSetbyModule.size() << std::endl;

			TCanvas* canvas = showgemstrips_for_ecalbin(binNum, uvStripSetbyModule);

			if (!canvas) {
				std::cerr << "Error: canvas is null for bin " << binNum << "!" << std::endl;
				continue;
			}

			canvas->Update();
			canvas->Draw();
			gSystem->ProcessEvents();
		}
	} else {
		std::cout << "\n### SELECTED BINS MODE ###" << std::endl;

		for (const int binNum : selected_bins) {
			if (map_physicalUVStrips_byECalBin_byGEMMod.find(binNum) == map_physicalUVStrips_byECalBin_byGEMMod.end()) {
				std::cerr << "Warning: Bin number " << binNum << " not found in data.\n";
				continue;
			}

			const auto& uvStripSetbyModule = map_physicalUVStrips_byECalBin_byGEMMod[binNum];

			std::cout << "\n### BIN NUMBER: " << binNum << " ###" << std::endl;
			std::cout << "Number of Modules " << uvStripSetbyModule.size() << std::endl;

			TCanvas* canvas = showgemstrips_for_ecalbin(binNum, uvStripSetbyModule);

			if (!canvas) {
				std::cerr << "Error: canvas is null for bin " << binNum << "!" << std::endl;
				continue;
			}

			canvas->Update();
			canvas->Draw();
			gSystem->ProcessEvents();
		}
	}

	return 0;
}