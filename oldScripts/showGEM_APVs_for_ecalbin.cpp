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

// std::map<int, gemInfo> refModMap = GetGemInfoMap();
std::map<int, gemInfo> refModMap = GetGemInfoMap();

// auto refAPVkeys = apvInfoKeys{};
// auto apvInfoKeys = returnAPVInfoKeys();


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
int uAPVs=0, vAPVs=0;



void draw_yx_GEMStrip(int modNum, int axis, int strip_num, int ROIcenter_strip)
{

	//NOTE:: calculates x and y values on canvas coordinates (+x to the right, +y up) then plots them in gem coordinates (+x down, +y left)
	

	// if(strip_num%20==0){//Testing with greater	plotted spacing for visibility
	if(strip_num){//Testing with greater	plotted spacing for visibility
	auto& modInfo = gemInfoMap[modNum];

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

	if (modNum){std::cout << "\nCenter strip " << center_strip << std::endl;}



	//NOTE: indexes strip number so that the center strip(assumed to be approx at 0 on the module) is 0 and converts to physical distance with pitch

	if (modNum){std::cout << "\n strip_num " << strip_num << " on axis " << axis << std::endl;}

	double strip_num_offset = strip_num - center_strip;

	
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
	double x_center = (distFromCenter * fPx  + (mod_x+getOffset(modNum, axis)));

	
	double y_center = (mod_y+getOffset(modNum, axis) + distFromCenter * fPy) ;


	TMarker* OffsetMark = new TMarker(distFromCenter * fPy, distFromCenter * fPx, 20);
	OffsetMark->SetMarkerSize(0.3);
	OffsetMark->SetMarkerColor(kBlue);
	// OffsetMark->Draw("same");

	TMarker* centerMark = new TMarker(y_center, x_center, 20);
	centerMark->SetMarkerSize(0.3);
	// centerMark->Draw("same");

	double half_length = 0.5 * (TMath::Sqrt(
		TMath::Sq(fabs(mod_size_x * cos(angle))) +
		TMath::Sq(fabs(mod_size_y * sin(angle))))
	);
	

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
        if (strip_num%(n_strips/2) -1){
        // TMarker* APVMark = new TMarker(y2, distFromCenter * fPx, 20);
        TMarker* APVMark = new TMarker(y1, x1, 21);
        APVMark->SetMarkerSize(0.3);
        APVMark->SetMarkerColor(kBlack);
        // APVMark->SetMarkerStyle(21);
        APVMark->Draw("same");
        uAPVs++;
        }


	} else {         // V
		stripLine->SetLineColor(kPink);
		// stripLine->SetLineStyle(2);
		// stripLine->SetLineStyle(9);
		stripLine->SetLineWidth(1);
		// stripLine->SetLineWidth(.3);
		vDrawn++;
        }
        if (strip_num%(n_strips/2) -1){
        // TMarker* APVMark = new TMarker(y1, distFromCenter * fPx, 20);
        TMarker* APVMark = new TMarker(y2, x2, 21);
        APVMark->SetMarkerSize(0.3);
        APVMark->SetMarkerColor(kBlack);
        // APVMark->SetMarkerStyle(21);
        APVMark->Draw("same");
        vAPVs++;
        }
	stripLine->Draw("same");
    }

}
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

	if (modNum){std::cout << "\nCenter strip " << center_strip << std::endl;}



	//NOTE: indexes strip number so that the center strip(assumed to be approx at 0 on the module) is 0 and converts to physical distance with pitch

	if (modNum){std::cout << "\n strip_num " << strip_num << " on axis " << axis << std::endl;}

	double strip_num_offset = strip_num - center_strip;

	
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
		-getOffset(modNum, axis)//shifted upward toward -x b/c negative offset
	));

	
	double y_center = (mod_y
		-getOffset(modNum, axis)//shifted upward toward -x b/c negative offset 
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
		
			roi_left->SetLineColor(kGreen);
			roi_right->SetLineColor(kGreen);
			roi_top->SetLineColor(kGreen);
			roi_bottom->SetLineColor(kGreen);
		
			roi_left->SetLineWidth(2);
			roi_right->SetLineWidth(2);
			roi_top->SetLineWidth(2);
			roi_bottom->SetLineWidth(2);
		
			roi_left->Draw("same");
			roi_right->Draw("same");
			roi_top->Draw("same");
			roi_bottom->Draw("same");
		}

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


int exportPDF(TCanvas * C){
	C->SaveAs("GEM_ECalBin.pdf");
	return 0;
}


void drawAPV_TOPDOWN(apvInfoKeys APVinfo){

	// if(strip_num%20==0){//Testing with greater	plotted spacing for visibility
		// if(strip_num){//Testing with greater	plotted spacing for visibility
	auto modNum = APVinfo.gemid;
	auto& modInfo = gemInfoMap[modNum];
	auto axis = APVinfo.axis;

	// int nAPVs = GetNAPVs(modNum)[axis];
	int nAPVs = modInfo.nstripsuv[axis]/128;
	
	// std::cout << "\nModNum: " << APVinfo.gemid
	// << ", axis " << axis
	// << ", pos " << APVinfo.pos
	// << ", nAPVs " << nAPVs
	// << std::endl;

	double mod_x = modInfo.position[0];
	double mod_y = modInfo.position[1];
	double mod_size_x = modInfo.size[0];
	double mod_size_y = modInfo.size[1];
	// double angle = modInfo.uvangles[axis] * TMath::DegToRad();


	double apvPos;

	// if (axis == 0) {
	// 	apvPos = (mod_x-mod_size_x / 2.0) //top of mod
	// 	- getOffset(APVinfo.gemid, axis)
	// 	+(APVinfo.pos)*128*0.0004;
	
	// }
	// else if (axis == 1){
	// apvPos = (mod_x+mod_size_x / 2.0) //top of mod
	// -getOffset(APVinfo.gemid, axis)
	// -APVinfo.pos*128*0.0004;
	// }

	int topAPV=nAPVs-1;

	int centerAPV = topAPV/2.0;
	

	int APVchanPos =APVinfo.pos;

	

	// double APVoffsetPos = //offset form top of module
	
	



	// double APVoffsetPos = //offset form top of module
	// topAPV - APVchanPos //counting from 0
	// -nAPVs/2.0; //translate to center of module
	// APVchanPos //counting from 0
	// -nAPVs/2.0; //translate to center of module

	// int APVoffsetPos = ((nAPVs/2-1))-APVchanPos;
	// int APVoffsetTop = ((nAPVs/2-1))-APVchanPos; //delta from top of module
	// // APVchanPos = APVoffsetPos;

	// ;

	double APVposOffset = 0.0;
	;
	


	double edgePos;
	
	if (modNum<6){

		// edgePos = -edgePos;//count top to bottom???

		if (modNum==0){// so v is on right 
			if (axis == 1) {edgePos = mod_size_y / 2.0; 
			}
			else if (axis == 0) {edgePos = -mod_size_y / 2.0;
		}
		else if(modNum>0 && modNum<6){
			if (axis == 0) {edgePos = mod_size_y / 2.0; 
			}
			else if (axis == 1) {edgePos = -mod_size_y / 2.0;
			}}

		//if need can say 1/2 on left and 1/2 on right
		


		APVposOffset = mod_size_x / 2.0 //physical dist
		- (APVchanPos*128*0.0004) //counting from 0
		;

		apvPos = mod_x + APVposOffset //top of mod
		- getOffset(APVinfo.gemid, axis)
		;

		std::cout << "\nAPVchanPos " << APVchanPos << std::endl;
		std::cout << "nAPVs " << nAPVs << std::endl;
		std::cout << "\nAPVpos " << apvPos << std::endl;
		std::cout << "axis " << axis << std::endl;
		std::cout << "\nAPVposOffset " << APVposOffset << std::endl;
		std::cout << "\nedgePos " << edgePos << std::endl;
		std::cout << "Mod num" << modNum << std::endl;


	}
	

	else if (modNum>=6) { 

		if (axis==1){//just flipping names for ease of input
		// if (axis==0){  //top of mod
			
			apvPos = (mod_y - mod_size_y / 2.0)
			+((APVchanPos)*128*0.0004);

			edgePos = (mod_x+mod_size_x / 2.0) //right edge of mod
			;

		}

		else if (axis==0){

			// (mod_x+mod_size_x / 2.0) //top of mod
			// -(APVchanPos)*128*0.0004;
			//top of mod
			edgePos = mod_x
			+(APVchanPos*128*0.0004)
			-(mod_size_x / 2.0);

			apvPos = mod_y+mod_size_y/ 2.0;
			// apvPos = mod_x+mod_size_x/ 2.0;
		}
		
	// apvPos = -apvPos;//count top to bottom
	}
	
	//TODO: check if this is correct 
	
	
	
	if (modNum < 6) {
		std::cout << "\nAngle " << modInfo.uvangles[axis] << std::endl;
	

	std::cout  << "\nModNum: " << APVinfo.gemid<< std::endl;
	std::cout << "\nAPV pos " << apvPos << std::endl;
	std::cout << "edgePos " << edgePos << std::endl;
	std::cout << "mod_x " << mod_x << std::endl;
	std::cout << "mod_y " << mod_y << std::endl;

	std::cout << "\nModNum: " << APVinfo.gemid
	<< ", axis " << axis
	<< ", pos " << APVinfo.pos
	<< ", nStrips " << modInfo.nstripsuv[axis]
	<< ", nAPVs " << nAPVs
	<< ", offset" << getOffset(APVinfo.gemid, axis)
	<< std::endl;

	}


	// apvPos = //count up in + of axes proj on x,y
	// (mod_x-mod_size_x / 2.0) //bottom of mod
	// -getOffset(APVinfo.gemid, axis)//b/c given positive in db
	// +(APVinfo.pos)*128*0.0004 //count up in + of axes proj on x,y
	// ;


	// apvPos += 128*0.0004/2.0; //center of APV

	TMarker* APVmark;


	//-!!! b/c flips!


	if (modNum>5){
	APVmark = new TMarker(-apvPos, edgePos, 21);
	}
	else if (modNum<6){

	APVmark = new TMarker(-edgePos, apvPos, 21);
		// -edgePos, -apvPos, 21);
	}
	APVmark->SetMarkerStyle(36);

	if (axis == 0) {
		APVmark->SetMarkerColor(kCyan);
	} else {
		// V
		APVmark->SetMarkerColor(kPink);
	}
	
	// APVmark->SetMarkerSize(0.3);
	APVmark->SetMarkerSize(1.3);
	APVmark->Draw("same");
	// }
	
}

}




void drawAPV(apvInfoKeys APVinfo){

	auto modNum = APVinfo.gemid;
	auto& modInfo = gemInfoMap[modNum];
	auto axis = APVinfo.axis;

	// int nAPVs = GetNAPVs(modNum)[axis];
	int nAPVs = modInfo.nstripsuv[axis]/128;
	double mod_x = modInfo.position[0];
	double mod_y = modInfo.position[1];
	double mod_size_x = modInfo.size[0];
	double mod_size_y = modInfo.size[1];


	double apvPos;


	int topAPV=nAPVs-1;

	int centerAPV = topAPV/2.0;
	

	int APVchanPos =APVinfo.pos;


	// int APVchanOffset = -nAPVs + APVchanPos;
	// int APVchanOffset = nAPVs - APVchanPos;
	int APVchanOffset = topAPV - APVchanPos;
	// int APVchanOffset = -topAPV + APVchanPos;

	//now can go from bottom up of mod up
	// APVchanOffset = - APVchanOffset ;

	

	double APVposOffset = 0.0;
	

	double edgePos;
	
	if (modNum<6){

		// edgePos = -edgePos;//count top to bottom???

		if (modNum==0){// so v is on right 
			if (axis == 1) {edgePos = mod_size_y / 2.0; 
			}
			else if (axis == 0) {edgePos = -mod_size_y / 2.0;
			}
		}
		else if(modNum>0 && modNum<6){
			if (axis == 0) {edgePos = mod_size_y / 2.0; 
			}
			else if (axis == 1) {edgePos = -mod_size_y / 2.0;
			}}

		//if need can say 1/2 on left and 1/2 on right


		// edgePos= - edgePos;//try this then flip xy
		


		APVposOffset = 
		// -(APVchanPos*128*0.0004) //counting from 0
		// -(APVchanOffset*128*.0004)	

		+(
			(APVchanPos)*128
			// +64
		)
		*0.0004
		-mod_size_x / 2.0
		// +(APVchanOffset*128*.0004)	
		;

		apvPos = + mod_x- getOffset(APVinfo.gemid, axis)
		// mod_size_x / 2.0 //physical dist
		
		// -mod_size_x / 2.0 //physical dist
		-APVposOffset //top of mod
		
		;

		// apvPos=-apvPos;

		
		std::cout << "\nAPVchanPos " << APVchanPos << std::endl;
		std::cout << "APV offset(index from top) " << APVchanOffset << std::endl;
		std::cout << "nAPVs " << nAPVs << std::endl;
		std::cout << "\nAPVpos " << apvPos << std::endl;
		std::cout << "centerAPV " << centerAPV << std::endl;
		std::cout << "axis " << axis << std::endl;
		std::cout << "\nAPVposOffset " << APVposOffset << std::endl;
		std::cout << "\nedgePos " << edgePos << std::endl;
		std::cout << "Mod num" << modNum << std::endl;

		std::cout << "(center index)*128*0.0004=" << (centerAPV)*128*0.0004 << std::endl;
		std::cout << "(APVchanPos)*128*0.0004=" << (APVchanPos)*128*0.0004 << std::endl;
		std::cout << "APVchanOffset*128*.0004=" << APVchanOffset*128*.0004 << std::endl;


		std::cout << "mod size x = " << mod_size_x << std::endl;
		std::cout << "nAPVs*128*pitch = " << nAPVs*(128)*.0004 << std::endl;


	}
	

	else if (modNum>=6) { 

		if (axis==1){//just flipping names for ease of input
			apvPos = (mod_y - mod_size_y / 2.0)
			+((APVchanPos)*128*0.0004);

			edgePos = (mod_x+mod_size_x / 2.0); //right edge of mod
			

		}
		else if (axis==0){
			edgePos = mod_x
			+(APVchanPos*128*0.0004)
			-(mod_size_x / 2.0);

			apvPos = mod_y+mod_size_y/ 2.0;
		}
		
	// apvPos = -apvPos;//count top to bottom
	}
	
	//TODO: check if this is correct 
	
	
	
	if (modNum < 6) {//Print info
		std::cout << "\nAngle " << modInfo.uvangles[axis] << std::endl;
	
	std::cout  << "\nModNum: " << APVinfo.gemid<< std::endl;
	std::cout << "\nAPV pos " << apvPos << std::endl;
	std::cout << "edgePos " << edgePos << std::endl;
	std::cout << "mod_x " << mod_x << std::endl;
	std::cout << "mod_y " << mod_y << std::endl;

	std::cout << "\nModNum: " << APVinfo.gemid
	<< ", axis " << axis
	<< ", pos " << APVinfo.pos
	<< ", nStrips " << modInfo.nstripsuv[axis]
	<< ", nAPVs " << nAPVs
	<< ", offset" << getOffset(APVinfo.gemid, axis)
	<< std::endl;

	}





	// apvPos += 128*0.0004/2.0; //center of APV

	TMarker* APVmark;

	if (modNum>5){
	APVmark = new TMarker(-apvPos, edgePos, 21);
	}
	else if (modNum<6){

	APVmark = new TMarker(edgePos, apvPos, 21);
		// -edgePos, -apvPos, 21);
	}
	APVmark->SetMarkerStyle(36);

	if (axis == 0) {
		APVmark->SetMarkerColor(kCyan);
	} else {
		// V
		APVmark->SetMarkerColor(kPink);
	}
	
	// APVmark->SetMarkerSize(0.3);
	APVmark->SetMarkerSize(1.3);
	APVmark->Draw("same");
	// }
	
}




void showAPVs_for_ecalbin(
	TCanvas * C,
	int ecalBinNum, 
	std::set<apvInfoKeys> binAPVkeysByModule
	) {
    // showAPVs_for_ecalbin(ecalBinNum, StripsPerMod);
    std::cout << "showAPVs_for_ecalbin() called" << std::endl;

	
	
	
	for (auto & key : binAPVkeysByModule) {
		int mod = key.gemid;
		std::cout << "\nModule " << key.gemid
		<< ", axis " << key.axis
		<< ", pos " << key.pos << std::endl;
		std::cout << std::endl;

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

		drawAPV(key);

		// if (key== binAPVkeysByModule.){}

	}
	
	C->Update();

}



//NOTE: want n/2 to have edges of APVs at 1st and last strip
// if (stripNum%(numStrips/2))




//Main func
int showGEM_APVs_for_ecalbin(const std::string& db_local = "db_FT_local.dat", const std::string& roi_file = "ROI_GEP3_FT_1.txt") {

	// using namespace GEMstripsHit_for_ecalbin;
	
	DBread db{ db_local };

	ROIread roi{ roi_file };

	if ( db.returnFileReadStatus() == -1 || roi.returnFileReadStatus() == -1 )
	{
		std::cerr << "Exiting the program!!!" << std::endl << std::endl;

		return -1;
	}

	// auto refAPVMap = db.returnAPVInfoMap();
	auto refAPVMap = apvInfoMap;


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


	 
	 //RREDITS:******
	 std::map < int /*ECalBinNo*/, std::set<apvInfoVals>> ECalBinAPVvals;
	 std::set<struct apvInfoKeys> missingKeys;

	//NOTE:Dont need Mod Num b/c it is in the key of the map
	 std::map<int /*ECalBinNo*/, std::set<apvInfoKeys>> ECalBinAPVkeys;


	for ( auto& [binNum, map_ROIbyLayer_forThisBin] : map_ROIsByBinsAndLayers )
	{
		std::map <int, 
			std::pair < std::set<int>, /*setOfUstripsForModule*/
			std::set<int> /*setOfVstripsForModule*/
			>
		> map_allUVstripsSetsForAllModules_forThisECalBin;

		//std::cout << "***Filling bin: " << binNum << std::endl;

		//RREDIT:Start
		std::set<apvInfoKeys> currECalBinAPV_Keys;

		std::set<apvInfoVals>currECalbinAPV_Vals; 
		
		//NOTE: dont need Mod Num b/c it is in the key of the map
		// std::set<apvInfoVals> currECalbinAPV_Info; 

		//RREDIT:End
		
		

		for ( auto& [layerNum, gemLayerROItoStripsInstance] : map_GEMLayerROItoStrips )
		{	
			std::map <int, /*ModNum*/
				std::pair < 
					std::set<int>, /*setOfUstripsForModule*/
					std::set<int> /*setOfUstripsForModule*/
				>
			> CurrLayer_ROI_UVstrips_perMod_forCurrECalBin 

			= gemLayerROItoStripsInstance.takeROI_givePhysicalUVStrips( /*map_ROIbyLayer_forThisBin.at( layerNum )*/ (map_ROIsByBinsAndLayers.at(binNum).at(layerNum)) );
			
			//std::cout << "Layer Num: " << layerNum << "  Numer of ROI modules: " << map_thisLayer_allROIsForAllModules_forThisECalBin.size() << std::endl;
			
			for (const auto& [modNum, modUVstripPair] : CurrLayer_ROI_UVstrips_perMod_forCurrECalBin) 
			{
				// If key doesn't exist in finalMap, insert it
			    if (map_allUVstripsSetsForAllModules_forThisECalBin.find(modNum) == map_allUVstripsSetsForAllModules_forThisECalBin.end())
			    {
			    	map_allUVstripsSetsForAllModules_forThisECalBin[modNum] = modUVstripPair;
			    } 
			
				{map_physicalUVStrips_byECalBin_byGEMMod[binNum] = map_allUVstripsSetsForAllModules_forThisECalBin;}


				// RREDIT:Start
				std::set<int>uStrips = modUVstripPair.first;
				std::set<int>vStrips = modUVstripPair.second;
				//RREDIT:END

				// std::set <int, std::set<apvInfoVals> > currModAPV_Keys;

				for (int stripID : uStrips){
					apvInfoKeys currStripKeys;
					currStripKeys.gemid = modNum;
					currStripKeys.axis = 0;


					currStripKeys.pos = stripID/128;//APV position

					// posInAPV = stripID%128;

					// std::cout << "ECalBin: " << binNum << " Checking modID: " << modNum << " Strip: " << stripID 
					// << " Axis: " << currStripKeys.axis << " Pos: " << currStripKeys.pos << std::endl;

					if (refAPVMap.find(currStripKeys) == refAPVMap.end()) {
						std::cout << "\nFor ECal Bin " << binNum <<" For strip " << stripID  << "  --> APV key: " << currStripKeys.gemid <<", "<< currStripKeys.axis <<", "<<
						currStripKeys.pos << " NOT FOUND!\n" << std::endl;
						missingKeys.emplace(currStripKeys);
					}


					if (refAPVMap.find(currStripKeys) != refAPVMap.end()) {
						currECalbinAPV_Vals.emplace(refAPVMap[currStripKeys]);  // Store unique strip -> APV mapping
						currECalBinAPV_Keys.emplace(currStripKeys);
						// currModAPV_Keys.emplace(currStripKeys);

					}

				}

				; 
				for (int stripID : vStrips){
					struct apvInfoKeys currStripKeys;
					currStripKeys.gemid = modNum;
					currStripKeys.axis = 1;
					


					currStripKeys.pos = (stripID/128);//APV position(-1 so index starts at 0)
					// posInAPV = stripID%128;

					// std::cout << "ECalBin: " << binNum << " Checking modID: " << modNum << " Strip: " << stripID 
					// << " Axis: " << currStripKeys.axis << " Pos: " << currStripKeys.pos << std::endl;
					
					if (refAPVMap.find(currStripKeys) == refAPVMap.end()) {
						std::cout << "\nFor ECal Bin " << binNum <<" For strip " << stripID  << "  --> APV key: " << currStripKeys.gemid <<", "<< currStripKeys.axis <<", "<<
						currStripKeys.pos << " NOT FOUND!\n" << std::endl;
						missingKeys.emplace(currStripKeys);
					}	


				if (refAPVMap.find(currStripKeys) != refAPVMap.end()) {
					currECalbinAPV_Vals.emplace(refAPVMap[currStripKeys]);  // Store unique strip -> APV mapping
					
					currECalBinAPV_Keys.emplace(currStripKeys);
					// currModAPV_Keys.emplace(currStripKeys);

					}

				}

				// currECalBinAPVs.emplace(modNum, currModAPVs);
				// currECalBinAPV_Keys.emplace(currModAPV_Keys);
				
				//RREDITS:End
		

			}	


		map_physicalUVStrips_byECalBin_byGEMMod[binNum] = map_allUVstripsSetsForAllModules_forThisECalBin;

		}

		ECalBinAPVvals[binNum] = currECalbinAPV_Vals;//RREDIT

		ECalBinAPVkeys[binNum] = currECalBinAPV_Keys;

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



	// NOTE:Start using Cataloged info

	if (all_bins_mode) {
		std::cout << "\n### ALL BINS MODE ###" << std::endl;

		TCanvas* canvas;

		for (const auto& [binNum, uvStripSetbyModule] : map_physicalUVStrips_byECalBin_byGEMMod) {
			std::cout << "\n### BIN NUMBER: " << binNum << " ###" << std::endl;
			std::cout << "Number of Modules " << uvStripSetbyModule.size() << std::endl;

			canvas = showgemstrips_for_ecalbin(binNum, map_ROIsByBinsAndLayers[binNum], uvStripSetbyModule);

			if (!canvas) {
				std::cerr << "Error: canvas is null for bin " << binNum << "!" << std::endl;
				continue;
			}

			canvas->Update();
			canvas->Draw();
			gSystem->ProcessEvents();
		}

		
		for (const auto& [binNum, keys] : ECalBinAPVkeys){

			showAPVs_for_ecalbin(canvas, binNum, keys);

			showAPVs_for_ecalbin(canvas, binNum, ECalBinAPVkeys[binNum]);


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
			
			showAPVs_for_ecalbin(canvas, binNum, ECalBinAPVkeys[binNum]); 


			
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
