/*
 * Visualizes the active GEM U/V strips for one or more ECal bins by plotting
 * the corresponding Region of Interest (ROI) strips per GEM module layer.
 * 
 * The code reads in a static local GEM database and an ROI definition file,
 * then maps the physical strip numbers per GEM module. Each selected ECal bin
 * is rendered as a multipanel ROOT canvas with strip overlays, module geometry,
 * and ROI boxes for each GEM layer.
 * 
 * Use interactively or export the visualizations to PDF.
 * 
 * --------------------
 * To Run:
 * 
 * In ROOT, run the following command:
 * 
 *   analyzer .x showGEMstripsHit_for_ecalbin.cpp
 * 
 * You will be prompted to input ECal bin numbers (e.g., `10 11 12`) or press Enter to visualize all.
 * Then select whether to export to PDF or display the ROOT canvas interactively.
 * 
 * Ensure that the files `db_FT_local.dat` and your ROI file (default: `ROI_GEP3_FT_1.txt`)
 * are present in the same directory or specify them by name.
 *
 *
 * Rafael Ruiz  
 * Jefferson Lab – SULI Program  
 * May 2025
 */

#include <TMath.h>
#include <TLine.h>
#include <TH2I.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TLatex.h>

#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

#include "GEMModROItoStrips.h"
#include "DBread.h"
#include "ROIread.h"

// #include <TString.h>
// #include <TRandom.h>
// #include "showgemhit_for_ecalbin.C"


#include "GetGemInfoMap.cxx"

// namespace GEMstripsHit_for_ecalbin{

std::map<int, gemInfo> refModMap = GetGemInfoMap();

int global_canvas_id = 0;



double getOffset(int mod, int axis) {
    // Axis: 0 = U, 1 = V

    // Define the offsets per module
    const double u_offsets[14] = {
        0.0176, // m0
        0.0176, // m1
        0.0108, // m2
        0.0108, // m3
        0.0108, // m4
        0.0108, // m5
        0.0,    // m6
        0.0,    // m7
        0.0,    // m8
        0.0,    // m9
        0.0,    // m10
        0.0,    // m11
        0.0,    // m12
        0.0     // m13
    };

    const double v_offsets[14] = {
        0.0,    // m0
        0.0,    // m1
        0.0108, // m2
        0.0108, // m3
        0.0108, // m4
        0.0108, // m5
        0.0,    // m6
        0.0,    // m7
        0.0,    // m8
        0.0,    // m9
        0.0,    // m10
        0.0,    // m11
        0.0,    // m12
        0.0     // m13
    };

    if (mod < 0 || mod >= 14) {
        std::cerr << "Invalid module number: " << mod << std::endl;
        return 0.0;
    }

    return (axis == 0) ? u_offsets[mod] : v_offsets[mod];
}


int uDrawn=0, vDrawn=0;



void drawGEMStrip(int modNum, int axis, int strip_num, int ROIcenter_strip)
{

	//NOTE:: calculates x and y values on canvas coordinates (+x to the right, +y up) then plots them in gem coordinates (+x down, +y left)
	

	if(strip_num%20==0){//Testing with greater	plotted spacing for visibility
	// if(strip_num){//Testing with greater	plotted spacing for visibility
	auto& modInfo = refModMap[modNum];

	double mod_x = modInfo.position[0];
	double mod_y = modInfo.position[1];
	double mod_size_x = modInfo.size[0];
	double mod_size_y = modInfo.size[1];
	double angle = modInfo.uvangles[axis] * TMath::DegToRad();

	
	// angle = -angle; // clockwise 

	int n_strips = modInfo.nstripsuv[axis];
	double pitch = 0.0004;

	// std::cout << "strip_num " << strip_num << std::endl;


	double layerSizeX=0;
	double layerSizeY=0;

	if(modNum<6){ layerSizeX = mod_size_x;}
	if (modNum>=6){layerSizeX = mod_size_x *4;}

	int center_strip = (n_strips - 1) / 2.0 ;

	// if (modNum){std::cout << "\nCenter strip " << center_strip << std::endl;}



	//NOTE: indexes strip number so that the center strip(assumed to be approx at 0 on the module) is 0 and converts to physical distance with pitch

	
	double strip_num_offset = strip_num - center_strip;


	if (modNum<6){std::cout << "\n strip_num " << strip_num << " on axis " << axis << std::endl;
	std::cout << "\nCenter strip " << center_strip << std::endl;
	}

	
	double distFromCenter = (strip_num_offset) * pitch;//

	double fPx, fPy;
	double dx_offset, dy_offset;

	//NOTE:clockwise rotation is positive
	fPx = std::cos(angle);
	fPy = std::sin(angle);



	// NOTE: OFFSET PERPENDICULAR TO STRIP DIRECTION

	dx_offset = std::cos(angle + (TMath::Pi()/2)); // Offset direction
	dy_offset = std::sin(angle + (TMath::Pi()/2));



	//NOTE: xy coords of the center of the current strip(shifted perpendicular to the strip)
	double x_center = (distFromCenter * fPx  
	+ (mod_x
		-getOffset(modNum, axis)
	));

	
	double y_center = (mod_y
		-getOffset(modNum, axis) 
		+ distFromCenter * fPy) ;


	TMarker* OffsetMark = new TMarker(distFromCenter * fPy, distFromCenter * fPx, 20);
	OffsetMark->SetMarkerSize(0.3);
	OffsetMark->SetMarkerColor(kBlue);
	// OffsetMark->Draw("same");

	TMarker* centerMark = new TMarker(y_center, x_center, 20);
	centerMark->SetMarkerSize(0.3);
	// centerMark->Draw("same");

	// double half_length = 0.5 * (TMath::Sqrt(
	// 	TMath::Sq(fabs(mod_size_x * cos(angle))) +
	// 	TMath::Sq(fabs(mod_size_y * sin(angle))))
	// );

	double half_length = mod_size_x / 2.0;
	

	double x1 = x_center - half_length * dx_offset;
	double x2 = x_center + half_length * dx_offset;

	double y1, y2;
	
	if (modNum>=6 && axis==0){y1 = y_center - mod_size_y/2;
	y2 = y_center + mod_size_y/2;
	}
	else {y1 = y_center - half_length * dy_offset;	
		y2 = y_center + half_length * dy_offset;
		}

	// X center line (vertical): fix "x" value at 0 (your flipped), span full Y range
	TLine* xCenterLine = new TLine(
		-mod_size_y / 2.0, 0,
		+mod_size_y / 2.0, 0
	);

	// Y center line (horizontal): fix "y" value at 0 (your flipped), span full X range
	TLine* yCenterLine = new TLine(
		0, -mod_size_x / 2.0,
		0, +mod_size_x / 2.0
	);

	// xCenterLine->Draw("SAME");
	// yCenterLine->Draw("SAME");




	// Module boundaries in global coordinates
	double xmin = - layerSizeX / 2.0;
	double xmax = + layerSizeX / 2.0;

	double ymin = -mod_size_y / 2.0;
	double ymax = mod_size_y / 2.0;

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


	// Draw the clipped strip
	// TLine* stripLine = new TLine(x1, y1, x2, y2); //Drawn 
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

TCanvas* showgemstrips_for_ecalbin(int ecalBinNum, std::map<int, ROI> binROI, const std::map<int, std::pair<std::set<int>, std::set<int>>>& StripsPerMod)
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

		
		xMin = -layerSizeX / 2. -.1;// *3;
		xMax = +layerSizeX / 2. +.1;// *3;
		yMin = -mod_size_y / 2. -.1;// *3;
		yMax = +mod_size_y / 2. +.1;// *3;
		// xMin = std::floor(-layerSizeX / 2. *10)/10;// *3;
		// xMax = std::ceil(+layerSizeX / 2. *10)/10;// *3;
		// yMin = std::floor(-mod_size_y / 2. *10)/10;// *3;
		// yMax = std::ceil(+mod_size_y / 2. *10)/10;// *3;
		
		
		
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

		// dummy->Fill(binROI.xMin, binROI.yMin, 0);
		// dummy->Fill(binROI.xMax, binROI.yMax, 0);

		dummy->SetStats(0);
		dummy->Draw();


		double Layer_yMin = -layerSizeX / 2.;// *3;
		double Layer_yMax = +layerSizeX / 2.;// *3;
		double Layer_xMin = -mod_size_y / 2.;// *3;
		double Layer_xMax = +mod_size_y / 2.;// *3;
		TLine* Layer_left = new TLine(Layer_xMin, Layer_yMin, Layer_xMin, Layer_yMax);
			TLine* Layer_right = new TLine(Layer_xMax, Layer_yMin, Layer_xMax, Layer_yMax);
			TLine* Layer_top = new TLine(Layer_xMin, Layer_yMax, Layer_xMax, Layer_yMax);
			TLine* Layer_bottom = new TLine(Layer_xMin, Layer_yMin, Layer_xMax, Layer_yMin);
		
			// roi_left->SetLineColor(kGreen);
			// roi_right->SetLineColor(kGreen);
			// roi_top->SetLineColor(kGreen);
			// roi_bottom->SetLineColor(kGreen);
		
			Layer_left->SetLineWidth(2);
			Layer_right->SetLineWidth(2);
			Layer_top->SetLineWidth(2);
			Layer_bottom->SetLineWidth(2);
		
			Layer_left->Draw("same");
			Layer_right->Draw("same");
			Layer_top->Draw("same");
			Layer_bottom->Draw("same");



		// TLatex* latex = new TLatex();
		// latex->SetNDC(); // normalized coordinates (0–1)
		// latex->SetTextSize(0.05);
		// latex->DrawLatex(0.1, 0.9, Form("Layer %d", layer));

		// if (binROI.find(layer) != binROI.end()) {
		// 	const ROI& roi = binROI.at(layer);
		
		// 	double xMinROI = roi.yMin; // swap X and Y because your dummy hist Y is geom X
		// 	double xMaxROI = roi.yMax;
		// 	double yMinROI = roi.xMin;
		// 	double yMaxROI = roi.xMax;

		// 	// //xy already in GEM coords?
		// 	// double xMinROI = roi.xMin; // swap X and Y because your dummy hist Y is geom X
		// 	// double xMaxROI = roi.xMax;
		// 	// double yMinROI = roi.yMin;
		// 	// double yMaxROI = roi.yMax;

		// 	std::cout << "\nROI xMin, yMin: " << xMinROI << ", " << yMinROI
		// 	<< "\nROI xMax, yMax: " << xMaxROI << ", " << yMaxROI << std::endl;
		
		// 	TLine* roi_left = new TLine(xMinROI, yMinROI, xMinROI, yMaxROI);
		// 	TLine* roi_right = new TLine(xMaxROI, yMinROI, xMaxROI, yMaxROI);
		// 	TLine* roi_top = new TLine(xMinROI, yMaxROI, xMaxROI, yMaxROI);
		// 	TLine* roi_bottom = new TLine(xMinROI, yMinROI, xMaxROI, yMinROI);
		
		// 	roi_left->SetLineColor(kGreen);
		// 	roi_right->SetLineColor(kGreen);
		// 	roi_top->SetLineColor(kGreen);
		// 	roi_bottom->SetLineColor(kGreen);
		
		// 	roi_left->SetLineWidth(2);
		// 	roi_right->SetLineWidth(2);
		// 	roi_top->SetLineWidth(2);
		// 	roi_bottom->SetLineWidth(2);
		
		// 	roi_left->Draw("same");
		// 	roi_right->Draw("same");
		// 	roi_top->Draw("same");
		// 	roi_bottom->Draw("same");
		// }

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

		pad->SetGridx();
		pad->SetGridy();

		pad->cd();  

		// // double mod_y = 0.0;

		auto thisModInfo = refModMap[mod];
		
		double mod_x = thisModInfo.position[0];		
		
		double mod_y = thisModInfo.position[1];		


		
		double u_angle = thisModInfo.uvangles[0];
		double v_angle = thisModInfo.uvangles[1];
		double pitch = 0.0004;
		int n_u = thisModInfo.nstripsuv[0];
		int n_v = thisModInfo.nstripsuv[1];



		double mod_size_x = thisModInfo.size[0];
		double mod_size_y = thisModInfo.size[1];


		//NOTE: Plot V first as cyan(v) under pink(u) is less visible/readable

		std::cout << "\nU,V angles are " << u_angle <<", "<< v_angle <<std::endl;

		std::cout << "\n N uStrips = " << uvStripPair.first.size() << std::endl;
		std::cout << "N vStrips = " << uvStripPair.second.size() << std::endl;

		// std::cout << "\n N uStrips in Mod " << thisModInfo.nstripsuv[0] << std::endl;
		std::cout << "N vStrips in Mod " << thisModInfo.nstripsuv[1] << std::endl;	

		int local_center_v = (*uvStripPair.second.begin() + *uvStripPair.second.rbegin()) / 2;

		std::cout << "min v strip " << *uvStripPair.second.begin() << std::endl;
		std::cout << "max v strip " << *uvStripPair.second.rbegin() << std::endl;
		std::cout << "\nlocal center v strip " << local_center_v << std::endl;

		int vStripLoopCalls=0;
		for (auto& vStrip : uvStripPair.second){
		// for (auto& vStrip : std::array{10, 20 , 30 , 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000}){
			// std::cout << "\n number of V strips " << uvStripPair.second.size()<<std::endl;
			// std::cout << "\nDrawing vStrip " << vStrip << std::endl;

			vStripLoopCalls+=1;
			// int local_center_v = (*uvStripPair.second.begin() + *uvStripPair.second.rbegin()) / 2;
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


		std::cout << "\n N uStrips in Mod " << thisModInfo.nstripsuv[0] << std::endl;
		

		int local_center_u = (*uvStripPair.first.begin() + *uvStripPair.first.rbegin()) / 2;

		std::cout << "min v strip " << *uvStripPair.second.begin() << std::endl;
		std::cout << "max v strip " << *uvStripPair.second.rbegin() << std::endl;
		std::cout << "\nlocal center v strip " << local_center_v << std::endl;


		int uStripLoopCalls=0;
		
		// for (int uStrip=0 , uStrip=uStripuvStripPair.first.size; uStrip++){¸
		for (auto& uStrip : uvStripPair.first){
		// for (auto& uStrip : std::array{10, 20 , 30 , 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000}){

			// std::cout << "\nnumber of U strips " << uvStripPair.first.size() <<std::endl;
			// std::cout << "\nDrawing uStrip" << uStrip << std::endl;

			uStripLoopCalls+=1;
			// std::cout << "\ncalled uStrip loop " << uStripLoopCalls << std::endl;
			// int local_center_u = (*uvStripPair.first.begin() + *uvStripPair.first.rbegin()) / 2;
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

	
		if (binROI.find(layer) != binROI.end()) {
			const ROI& roi = binROI.at(layer);
		
			double xMinROI = roi.yMin; // swap X and Y because your dummy hist Y is geom X
			double xMaxROI = roi.yMax;
			double yMinROI = roi.xMin;
			double yMaxROI = roi.xMax;

			// //xy already in GEM coords?
			// double xMinROI = roi.xMin; // swap X and Y because your dummy hist Y is geom X
			// double xMaxROI = roi.xMax;
			// double yMinROI = roi.yMin;
			// double yMaxROI = roi.yMax;

			std::cout << "\nROI xMin, yMin: " << xMinROI << ", " << yMinROI
			<< "\nROI xMax, yMax: " << xMaxROI << ", " << yMaxROI << std::endl;
		
			TLine* roi_left = new TLine(xMinROI, yMinROI, xMinROI, yMaxROI);
			TLine* roi_right = new TLine(xMaxROI, yMinROI, xMaxROI, yMaxROI);
			TLine* roi_top = new TLine(xMinROI, yMaxROI, xMaxROI, yMaxROI);
			TLine* roi_bottom = new TLine(xMinROI, yMinROI, xMaxROI, yMinROI);
		
			roi_left->SetLineColor(1);
			roi_right->SetLineColor(1);
			roi_top->SetLineColor(1);
			roi_bottom->SetLineColor(1);
			// roi_left->SetLineColor(kP6Violet);
			// roi_right->SetLineColor(kP6Violet);
			// roi_top->SetLineColor(kP6Violet);
			// roi_bottom->SetLineColor(kP6Violet);
		
			roi_left->SetLineWidth(2);
			roi_right->SetLineWidth(2);
			roi_top->SetLineWidth(2);
			roi_bottom->SetLineWidth(2);
		
			roi_left->Draw("same");
			roi_right->Draw("same");
			roi_top->Draw("same");
			roi_bottom->Draw("same");
		}


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

	// C->SetBit(kMustCleanup); // ROOT will clean safely
	C->Update();
	// gSystem->ProcessEvents();
	// std::cout << "Press enter to exit...\n";
	// std::cin.get();
	return C;
}


int exportPDF(TCanvas *C, int binNum) {
    C->SaveAs(Form("GEM_ECalBin%d.pdf", binNum));
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


	//NOTE: start "cataloging"

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




	//NOTE: start EcalBin to strip mapping
	
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

// 		std::cout << "For bin " << binNum << " collected modules:\n";
// for (const auto& [modNum, modUVstripPair] : map_allUVstripsSetsForAllModules_forThisECalBin) {
//     std::cout << "Module " << modNum 
//               << " has U strips: " << modUVstripPair.first.size()
//               << ", V strips: " << modUVstripPair.second.size() << std::endl;
// }


		}
	}
	
	
	std::string usrinput;
	std::cout << "Enter ECal bin numbers to visualize separated by spaces, 'Return' to visualize all ECalBins or 'q' to quit: ";
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

	
	// --- NEW CLEAN INPUT: Ask about PDF saving ---
	bool save_as_pdf = false;
	std::string usr_save_as_pdf;
	std::cout << "Save as PDF(y) or show Active Canvas(n)? (y/n): ";
	std::getline(std::cin, usr_save_as_pdf);
	if (usr_save_as_pdf == "y" || usr_save_as_pdf == "Y") {
		save_as_pdf = true;
	}
	
	
	// bool interactive_view = false;
	// std::string usr_interactive_view;
	// std::cout << "Visualize with ROOT interactive window? (y/n): ";
	// std::getline(std::cin, usr_interactive_view);
	// if (usr_interactive_view == "y" || usr_interactive_view == "Y") {
	// 	interactive_view = true;
	// }






	// -- Printing info --
std::cout << "\nNumber of ECalBins is " << map_physicalUVStrips_byECalBin_byGEMMod.size() << std::endl;



	// NOTE:Start using Cataloged info

	if (all_bins_mode) {
		std::cout << "\n### ALL BINS MODE ###" << std::endl;

		for (const auto& [binNum, uvStripSetbyModule] : map_physicalUVStrips_byECalBin_byGEMMod) {
			std::cout << "\n### BIN NUMBER: " << binNum << " ###" << std::endl;
			std::cout << "Number of Modules " << uvStripSetbyModule.size() << std::endl;

			TCanvas* canvas = showgemstrips_for_ecalbin(binNum, map_ROIsByBinsAndLayers[binNum], uvStripSetbyModule);

			if (!canvas) {
				std::cerr << "Error: canvas is null for bin " << binNum << "!" << std::endl;
				continue;
			}

			canvas->Update();

			if (save_as_pdf==true){
				exportPDF(canvas, binNum);
				// return 0;
			}

			// canvas->Draw();
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

			TCanvas* canvas = showgemstrips_for_ecalbin(binNum, map_ROIsByBinsAndLayers[binNum], uvStripSetbyModule);

			if (!canvas) {
				std::cerr << "Error: canvas is null for bin " << binNum << "!" << std::endl;
				continue;
			}

			canvas->Update();

			if (save_as_pdf==true){
				exportPDF(canvas, binNum);
				// return 0;
			}

			// canvas->Draw();
			gSystem->ProcessEvents();
		}
	}

	return 0;
}