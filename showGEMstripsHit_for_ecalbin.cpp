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

#include <TString.h>
#include <TRandom.h>
#include "showgemhit_for_ecalbin.C"


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


	// if (modNum<6){std::cout << "\n strip_num " << strip_num << " on axis " << axis << std::endl;
	// std::cout << "\nCenter strip " << center_strip << std::endl;
	// }

	
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
	// double xmin = - layerSizeX / 2.0;
	// double xmax = + layerSizeX / 2.0;

	// double ymin = -mod_size_y / 2.0;
	// double ymax = mod_size_y / 2.0;
	double xmin = mod_x - layerSizeX / 2.0;
	double xmax = mod_x + layerSizeX / 2.0;
	double ymin = mod_y - mod_size_y / 2.0;
	double ymax = mod_y + mod_size_y / 2.0;


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
    std::cout << "\n\nStarting showgemstrips_for_ecalbin()" << std::endl;
    std::cout << "number of keys = " << StripsPerMod.size() << std::endl;
    std::cout << "Keys include: " << std::endl;
    for (auto& [key, val] : StripsPerMod) {
        std::cout << key << ": "
                  << "\nnumber of U strips " << val.first.size()
                  << ", number of V strips " << val.second.size() << std::endl;
    }

    int canvas_id = 0;
    TCanvas* C = new TCanvas(Form("canvas_%d", global_canvas_id++), "Active ROI Strips in GEM Planes", 1600, 1000);
    C->SetTitle(Form("Active ROI Strips for ECal Bin %d", ecalBinNum));
    C->Divide(4, 2);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(22);
    latex.DrawLatex(0.5, 0.97, Form("ECal Bin %d", ecalBinNum));

    for (int layer = 0; layer < 8; layer++) {
        std::cout << "\nMaking subPad for layer " << layer << "\n--------------------------------------" << std::endl;
        C->cd(layer + 1);
        gPad->SetGridx();
        gPad->SetGridy();

        std::vector<int> mods_in_layer;
        for (const auto& [modNum, modInfo] : refModMap) {
            int modLayer = (modNum <= 5) ? modNum : (modNum <= 9) ? 6 : (modNum <= 13) ? 7 : -1;
            if (modLayer == layer) mods_in_layer.push_back(modNum);
        }
        if (mods_in_layer.empty()) continue;

        double xMin = 1e9, xMax = -1e9;
        double yMin = 1e9, yMax = -1e9;

        for (int mod : mods_in_layer) {
            auto& info = refModMap[mod];
            double mx = info.position[0];
            double my = info.position[1];
            double sx = info.size[0];
            double sy = info.size[1];
            double lx = (mod <= 5) ? sx : sx * 4;

            xMin = std::min(xMin, mx - lx / 2.0);
            xMax = std::max(xMax, mx + lx / 2.0);
            yMin = std::min(yMin, my - sy / 2.0);
            yMax = std::max(yMax, my + sy / 2.0);
        }

        double layer_xMin = xMin;
        double layer_xMax = xMax;
        double layer_yMin = yMin;
        double layer_yMax = yMax;

        double pad_x = 0.1 * (xMax - xMin);
        double pad_y = 0.1 * (yMax - yMin)+0.1;
        xMin -= pad_x; xMax += pad_x;
        yMin -= pad_y; yMax += pad_y;

        std::cout << "\nxMin, yMin: " << xMin << ", " << yMin << "\nxMax, yMax: " << xMax << ", " << yMax << std::endl;

        TH2I* dummy = new TH2I(
            Form("dummy_layer%d_%d", layer, canvas_id),
            Form("Layer %d;Y (m);X (m)", layer),
            100, yMin, yMax, 100, xMin, xMax
        );

		dummy->SetDirectory(0);  // Prevent memory leak warning by detaching from ROOT directory
        dummy->SetStats(0);
        dummy->Draw();

        TLine* box_left   = new TLine(layer_yMin, layer_xMin, layer_yMin, layer_xMax);
        TLine* box_right  = new TLine(layer_yMax, layer_xMin, layer_yMax, layer_xMax);
        TLine* box_bottom = new TLine(layer_yMin, layer_xMin, layer_yMax, layer_xMin);
        TLine* box_top    = new TLine(layer_yMin, layer_xMax, layer_yMax, layer_xMax);

        for (auto* line : {box_left, box_right, box_top, box_bottom}) {
            line->SetLineColor(kBlack);
            line->SetLineWidth(2);
            line->Draw("same");
        }

        // if (binROI.find(layer) != binROI.end()) {
        //     const ROI& roi = binROI.at(layer);
        //     double xMinROI = roi.yMin;
        //     double xMaxROI = roi.yMax;
        //     double yMinROI = roi.xMin;
        //     double yMaxROI = roi.xMax;

        //     TLine* roi_left   = new TLine(xMinROI, yMinROI, xMinROI, yMaxROI);
        //     TLine* roi_right  = new TLine(xMaxROI, yMinROI, xMaxROI, yMaxROI);
        //     TLine* roi_top    = new TLine(xMinROI, yMaxROI, xMaxROI, yMaxROI);
        //     TLine* roi_bottom = new TLine(xMinROI, yMinROI, xMaxROI, yMinROI);

        //     for (auto* line : {roi_left, roi_right, roi_top, roi_bottom}) {
        //         line->SetLineColor(kBlack);
        //         line->SetLineStyle(2);
        //         line->SetLineWidth(2);
        //         line->Draw("same");
        //     }
        // }
    }

    for (const auto& [mod, uvStrips] : StripsPerMod) {
        int layer = (mod <= 5) ? mod : (mod <= 9) ? 6 : (mod <= 13) ? 7 : -1;
        if (layer < 0) continue;

        C->cd(layer + 1);

        if (uvStrips.first.empty() && uvStrips.second.empty()) {
            std::cerr << "[INFO] Skipping module " << mod << " in layer " << layer << " due to empty U and V strips." << std::endl;
            continue;
        }

        int centerU = uvStrips.first.empty() ? 0 : (*uvStrips.first.begin() + *uvStrips.first.rbegin()) / 2;
        int centerV = uvStrips.second.empty() ? 0 : (*uvStrips.second.begin() + *uvStrips.second.rbegin()) / 2;

        int rangeU = uvStrips.first.empty() ? -1 : (*uvStrips.first.rbegin() - *uvStrips.first.begin());
        int rangeV = uvStrips.second.empty() ? -1 : (*uvStrips.second.rbegin() - *uvStrips.second.begin());

        if (!uvStrips.second.empty() && rangeV < 2000 && rangeV >= 0) {
            for (int vStrip : uvStrips.second)
                drawGEMStrip(mod, 1, vStrip, centerV);
        } else if (!uvStrips.second.empty()) {
            std::cerr << "[WARNING] Skipping V strip draw for module " << mod << " due to large or invalid range: " << rangeV << std::endl;
        }

        if (!uvStrips.first.empty() && rangeU < 2000 && rangeU >= 0) {
            for (int uStrip : uvStrips.first)
                drawGEMStrip(mod, 0, uStrip, centerU);
        } else if (!uvStrips.first.empty()) {
            std::cerr << "[WARNING] Skipping U strip draw for module " << mod << " due to large or invalid range: " << rangeU << std::endl;
        }


		if (binROI.find(layer) != binROI.end()) {
            const ROI& roi = binROI.at(layer);
            double xMinROI = roi.yMin;
            double xMaxROI = roi.yMax;
            double yMinROI = roi.xMin;
            double yMaxROI = roi.xMax;

            TLine* roi_left   = new TLine(xMinROI, yMinROI, xMinROI, yMaxROI);
            TLine* roi_right  = new TLine(xMaxROI, yMinROI, xMaxROI, yMaxROI);
            TLine* roi_top    = new TLine(xMinROI, yMaxROI, xMaxROI, yMaxROI);
            TLine* roi_bottom = new TLine(xMinROI, yMinROI, xMaxROI, yMinROI);

            for (auto* line : {roi_left, roi_right, roi_top, roi_bottom}) {
                line->SetLineColor(kBlack);
                // line->SetLineStyle(2);
                // line->SetLineWidth(2);
                line->Draw("same");
            }
        }
    }

    C->Update();
    return C;
}



int exportPDF(TCanvas *C, int binNum) {
    C->SaveAs(Form("GEM_ECalBin%d.pdf", binNum));
    return 0;
}



//Main func
int showGEMstripsHit_for_ecalbin(const std::string& db_local = "db_FT_local.dat", 
	// const std::string& roi_file = "ROI_GEP3_FT_1.txt") {
	const std::string& roi_file = "ROI_GEP1_FT_realdat.txt") {

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
			
			
			
			for (const auto& [modNum, modUVstripPair] : map_thisLayer_allROIsForAllModules_forThisECalBin) 
			{
				

				// If key doesn't exist in finalMap, insert it
			    if (map_allUVstripsSetsForAllModules_forThisECalBin.find(modNum) == map_allUVstripsSetsForAllModules_forThisECalBin.end())
			    {
			    	map_allUVstripsSetsForAllModules_forThisECalBin[modNum] = modUVstripPair;
			    } 
		
			}

		map_physicalUVStrips_byECalBin_byGEMMod[binNum] = map_allUVstripsSetsForAllModules_forThisECalBin;

		// std::cout << "For bin " << binNum << " collected modules:\n";
		// for (const auto& [modNum, modUVstripPair] : map_allUVstripsSetsForAllModules_forThisECalBin) {
		// 	std::cout << "Module " << modNum 
		// 			<< " has U strips: " << modUVstripPair.first.size()
		// 			<< ", V strips: " << modUVstripPair.second.size() << std::endl;
		// }


		}
	}


	std::cout << "\nFilled all ROI, strip info etc" << std::endl;
	
	
	std::string usrinput;
	std::cout << "\nEnter ECal bin numbers to visualize separated by spaces, 'Return' to visualize all ECalBins or 'q' to quit: ";
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
	std::cout << "\nSave as PDF(y) or show Active Canvas(n)? (y/n): ";
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

			// if (save_as_pdf==true){
			// 	exportPDF(canvas, binNum);
			// 	// return 0;
			// }

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

			// if (save_as_pdf==true){
			// 	exportPDF(canvas, binNum);
			// 	// return 0;
			// }

			// canvas->Draw();
			gSystem->ProcessEvents();
		}
	}

	return 0;
}