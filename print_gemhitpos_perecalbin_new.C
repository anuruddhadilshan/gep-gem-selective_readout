#include "ecalbin2.h"
#include <vector>
#include <algorithm>
#include <fstream>

// Define the struct to hold the Region Of Interest (ROI) rectangle boundaries.
struct ROI
{
	double xMin, xMax, yMin, yMax;

	// Constructor requiring parameters
    ROI(double x_min, double x_max, double y_min, double y_max) 
        : xMin(x_min), xMax(x_max), yMin(y_min), yMax(y_max) {}
};

// Function to compute the X% ROI from a TH2F
ROI computeROI( TH2F* h2d, double percentage = 0.95 )
{
	// If there are no entries in the histogram, return an un-physical ROI.
	if ( h2d->GetEntries() == 0 )
	{
		return {std::numeric_limits<double>::quiet_NaN(), 
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
	}

	int nBinsX = h2d->GetNbinsX();
	int nBinsY = h2d->GetNbinsY();

	// Store (content, binX, binY)
	std::vector<std::tuple<double, int, int>> histBinData;
	double totalHitsInHist = 0.0;

	// Fill the vector with bin content and compute total number of hits
	for ( int i = 1; i <= nBinsX; i++ )
	{
		for ( int j = 1; j <= nBinsY; j++ )
		{
			double binContent = h2d->GetBinContent(i, j);

			if ( binContent > 0 )
			{
				histBinData.push_back({binContent, i, j});

				totalHitsInHist += binContent;
			}
		}
	}

	// Sort bins by content in the descending order.
	std::sort(histBinData.begin(), histBinData.end(), [](const auto& a, const auto& b) {
    	return std::get<0>(a) > std::get<0>(b); // Sort by the first element in descending order
	});

	// Accumulate bins until reaching X% of total hits
	double accumulated = 0.0;
	int cutoffIndex = 0;

	for (; cutoffIndex < histBinData.size(); cutoffIndex++)
	{
		accumulated += std::get<0>(histBinData[cutoffIndex]);
		if ( accumulated / totalHitsInHist >= percentage ) break;
	}

	// Find min/max x and y bins for the selected ROI
	double xMin = 1e9;
	double xMax = -1e9;
	double yMin = 1e9;
	double yMax = -1e9;

	for ( int k = 0; k <= cutoffIndex; k++ )
	{
		int binX = std::get<1>(histBinData[k]);
		int binY = std::get<2>(histBinData[k]);

		double xCenter = h2d->GetXaxis()->GetBinCenter(binX);
		double yCenter = h2d->GetYaxis()->GetBinCenter(binY);

		xMin = std::min(xMin, xCenter);
		xMax = std::max(xMax, xCenter);
		yMin = std::min(yMin, yCenter);
		yMax = std::max(yMax, yCenter);
	}
	// Return ROI boundary rectangle
	return {xMin, xMax, yMin, yMax};
}

// Function to find ROI using max bin content expansion method
ROI computeROI_MaxBinExpand(TH2F* h2d, double thresholdFactor = 0.01) 
{
    int nBinsX = h2d->GetNbinsX();
    int nBinsY = h2d->GetNbinsY();

    // Find the bin with the highest content
    int maxBinX = -1, maxBinY = -1;
    double maxContent = 0.0;

    for (int i = 1; i <= nBinsX; i++) {
        for (int j = 1; j <= nBinsY; j++) {
            double content = h2d->GetBinContent(i, j);
            if (content > maxContent) {
                maxContent = content;
                maxBinX = i;
                maxBinY = j;
            }
        }
    }

    if (maxBinX == -1 || maxBinY == -1) {
        std::cerr << "No valid max bin found!" << std::endl;
        return {std::numeric_limits<double>::quiet_NaN(), 
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()}; // Return empty ROI
    }

    double threshold = maxContent * thresholdFactor;

    // Expand in four directions
    int leftX = maxBinX, rightX = maxBinX;
    int downY = maxBinY, upY = maxBinY;

    while (leftX > 1 && h2d->GetBinContent(leftX - 1, maxBinY) > threshold) leftX--;
    while (rightX < nBinsX && h2d->GetBinContent(rightX + 1, maxBinY) > threshold) rightX++;
    while (downY > 1 && h2d->GetBinContent(maxBinX, downY - 1) > threshold) downY--;
    while (upY < nBinsY && h2d->GetBinContent(maxBinX, upY + 1) > threshold) upY++;

    // Convert bin indices to coordinates
    leftX = std::max(leftX - 1, 1);
	rightX = std::min(rightX + 1, nBinsX);
	downY = std::max(downY - 1, 1);
	upY = std::min(upY + 1, nBinsY);

    double xMin = h2d->GetXaxis()->GetBinLowEdge(leftX);
    double xMax = h2d->GetXaxis()->GetBinUpEdge(rightX);
    double yMin = h2d->GetYaxis()->GetBinLowEdge(downY);
    double yMax = h2d->GetYaxis()->GetBinUpEdge(upY);

    return {xMin, xMax, yMin, yMax};
}

// Function to find ROI using max bin content expansion method
ROI computeROI_MaxBinExpand_2(TH2F* h2d, double maxbin_thresholdFactor = 0.1, double thresholdFactor = 0.01) 
{
    // If there are no entries in the histogram, return an un-physical ROI.
	if ( h2d->GetEntries() == 0 )
	{
		return {std::numeric_limits<double>::quiet_NaN(), 
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
	}

	int nBinsX = h2d->GetNbinsX();
	int nBinsY = h2d->GetNbinsY();

	// Store (content, binX, binY)
	std::vector<std::tuple<double, int, int>> histBinData;
	
	// Fill the vector with bin content and compute total number of hits
	for ( int i = 1; i <= nBinsX; i++ )
	{
		for ( int j = 1; j <= nBinsY; j++ )
		{
			double binContent = h2d->GetBinContent(i, j);

			if ( binContent > 0 )
			{
				histBinData.push_back({binContent, i, j});				
			}
		}
	}

	// Sort bins by content in the descending order.
	std::sort(histBinData.begin(), histBinData.end(), [](const auto& a, const auto& b) {
    	return std::get<0>(a) > std::get<0>(b); // Sort by the first element in descending order
	});

	// See if there is at least ~z% max bin's hits in one of the neighbouring bins to make sure we are not selecting one of random hot channles/regoins.
	// If not move, on to the next max bin.
	int maxBinX = -1, maxBinY = -1;
	double threshold;
	
	for (const auto& element : histBinData) 
	{
        double binContent;
        int ix, iy;
        std::tie(binContent, ix, iy) = element; // Unpacking tuple

        if ( h2d->GetBinContent(ix + 1, iy) / binContent >= maxbin_thresholdFactor || h2d->GetBinContent(ix - 1, iy) / binContent >= maxbin_thresholdFactor || h2d->GetBinContent(ix, iy + 1) / binContent >= maxbin_thresholdFactor || h2d->GetBinContent(ix, iy - 1) / binContent >= maxbin_thresholdFactor ) 
        {
        	maxBinX = ix;
        	maxBinY = iy;
        	threshold = binContent * thresholdFactor;
        	break;
        }
    }

    if (maxBinX == -1 || maxBinY == -1) 
    {
        std::cerr << "No valid max bin found!" << std::endl;
        return {std::numeric_limits<double>::quiet_NaN(), 
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()}; // Return empty ROI
    }

    // Expand in four directions
    int leftX = maxBinX, rightX = maxBinX;
    int downY = maxBinY, upY = maxBinY;

    while (leftX > 1 && h2d->GetBinContent(leftX - 1, maxBinY) > threshold) leftX--;
    while (rightX < nBinsX && h2d->GetBinContent(rightX + 1, maxBinY) > threshold) rightX++;
    while (downY > 1 && h2d->GetBinContent(maxBinX, downY - 1) > threshold) downY--;
    while (upY < nBinsY && h2d->GetBinContent(maxBinX, upY + 1) > threshold) upY++;

    // Convert bin indices to coordinates with adding one bin (~APV) safety margins.
    leftX = std::max(leftX - 1, 1);
	rightX = std::min(rightX + 1, nBinsX);
	downY = std::max(downY - 1, 1);
	upY = std::min(upY + 1, nBinsY);

    double xMin = h2d->GetXaxis()->GetBinLowEdge(leftX);
    double xMax = h2d->GetXaxis()->GetBinUpEdge(rightX);
    double yMin = h2d->GetYaxis()->GetBinLowEdge(downY);
    double yMax = h2d->GetYaxis()->GetBinUpEdge(upY);

    return {xMin, xMax, yMin, yMax};
}

// Function to check if a point is inside the ROI.
bool isInsideROI(const double x, const double y, const ROI& roi)
{
	return ( x >= roi.xMin && x <= roi.xMax && y >= roi.yMin && y <= roi.yMax );
}

// Function to write `vvROI` to a text file
void writeROIToFile(const std::vector<std::vector<ROI>>& vvROI, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Write header
    outFile << "# ECalBin GEMLayer xMin xMax yMin yMax\n";

    // Loop over the ROI vector and write data
    for (size_t i = 0; i < vvROI.size(); i++)
    {  // Loop over ECal bins
        for (size_t j = 0; j < vvROI[i].size(); j++)
        {  // Loop over GEM layers
            const ROI& roi = vvROI[i][j];

            // Write data in a structured format
            // IMPORTANT: To reflect the actual physical setup, we have drawn the Transport coordinate system X axisi on the histogram vertical axis and the Y axisi on histogram horizontal axis.
            // But in order to compare with the DB GEM positions we should put these back to transport coordiates.
            outFile << i << " " << j << " " << roi.yMin << " " << roi.yMax << " " << roi.xMin << " " << roi.xMax << "\n";
        }
    }

    outFile.close();
    std::cout << "ROI data successfully written to " << filename << std::endl;
}



void print_gemhitpos_perecalbin_new( const char* simfilename = "dummy" )
{
	
	TChain* T = new TChain("T");

	// T->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep_sim/g4sbs_out/GEP3_elastic_targzoff_9cm_preinit_006mmlead_job_1*.root");
	// T->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep_sim/g4sbs_out/GEP3_elastic_targzoff_9cm_preinit_006mmlead_job_2*.root");
	// T->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep_sim/g4sbs_out/GEP3_elastic_targzoff_9cm_preinit_006mmlead_job_4*.root");
	// T->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep_sim/g4sbs_out/GEP3_elastic_targzoff_9cm_preinit_006mmlead_job_5*.root");

	T->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep_sim/g4sbs_out/GEP3_elastic_targzoff_9cm_preinit_006mmlead_job_*.root");

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("Earm.ECalTF1.hit.nhits", 1);
	T->SetBranchStatus("Earm.ECalTF1.hit.row", 1);
	T->SetBranchStatus("Earm.ECalTF1.hit.col", 1);
	T->SetBranchStatus("Earm.ECalTF1.hit.cell", 1);
	T->SetBranchStatus("Earm.ECalTF1.hit.sumedep", 1);
	T->SetBranchStatus("Earm.ECalTF1.hit.xcell", 1);
	T->SetBranchStatus("Earm.ECalTF1.hit.ycell", 1);
	T->SetBranchStatus("Harm.FT.hit.nhits", 1);
	T->SetBranchStatus("Harm.FT.hit.plane", 1);
	T->SetBranchStatus("Harm.FT.hit.x", 1);
	T->SetBranchStatus("Harm.FT.hit.y", 1);
	T->SetBranchStatus("Harm.FT.Track.ntracks", 1);
	T->SetBranchStatus("Harm.FT.Track.MID", 1);
	T->SetBranchStatus("Harm.FT.Track.NumHits", 1);
	//T->SetBranchStatus("Harm.HCalScint.hit.sumedep", 1);

	int n_ecalhits {0};
	int n_gemhits {0};
	std::vector<int>* ecalhit_row = nullptr;
	std::vector<int>* ecalhit_col = nullptr;
	std::vector<int>* ecalhit_cell = nullptr;
	std::vector<double>* ecalhit_sumedep = nullptr;
	std::vector<double>* ecalhit_xcell = nullptr;
	std::vector<double>* ecalhit_ycell = nullptr;
	std::vector<int>* gemhit_plane = nullptr;
	std::vector<double>* gemhit_x = nullptr;
	const double g4sbsFT_xoffset { 0.165 }; // This need to be added to the Harm.FT.hit.x value to get the GEM hit pos in X with the origin of the coordinate system centered at the GEM module origin.
	// This is because G4SBS TRASPORT coordinate origin is at the heigh of the beam heigh and NOT at the height of the FT first GEM layer's active area center.	
	std::vector<double>* gemhit_y = nullptr;
	int n_gemtracks {0};
	std::vector<int>* gemtrack_mid = nullptr;
	std::vector<int>* gemtrack_numhits = nullptr;
	//std::vector<double>* hcalScintHit_sumedep = nullptr;

	T->SetBranchAddress("Earm.ECalTF1.hit.nhits", &n_ecalhits);
	T->SetBranchAddress("Earm.ECalTF1.hit.row", &ecalhit_row);
	T->SetBranchAddress("Earm.ECalTF1.hit.col", &ecalhit_col);
	T->SetBranchAddress("Earm.ECalTF1.hit.cell", &ecalhit_cell);
	T->SetBranchAddress("Earm.ECalTF1.hit.sumedep", &ecalhit_sumedep);
	T->SetBranchAddress("Earm.ECalTF1.hit.xcell", &ecalhit_xcell);
	T->SetBranchAddress("Earm.ECalTF1.hit.ycell", &ecalhit_ycell);
	T->SetBranchAddress("Harm.FT.hit.nhits", &n_gemhits);
	T->SetBranchAddress("Harm.FT.hit.plane", &gemhit_plane);
	T->SetBranchAddress("Harm.FT.hit.x", &gemhit_x);
	T->SetBranchAddress("Harm.FT.hit.y", &gemhit_y);
	T->SetBranchAddress("Harm.FT.Track.ntracks", &n_gemtracks);
	T->SetBranchAddress("Harm.FT.Track.MID", &gemtrack_mid);
	T->SetBranchAddress("Harm.FT.Track.NumHits", &gemtrack_numhits);
	//T->SetBranchStatus("Harm.HCalScint.hit.sumedep", &hcalScintHit_sumedep);

	// Initialize the ECalbin class.
	Ecalbin ecalbin{};

	const int n_gemlayers = 8;
	const int n_ecalbins = ecalbin.numOfECalBinsNew();

	// Creating a vector of vectors.
	std::vector<std::vector<TH2F*>> gemhithists_perecalbin( n_ecalbins, std::vector<TH2F*>(n_gemlayers, nullptr) );

	// Creating Canvas
	std::vector<TCanvas*> canvas_perecalbin( n_ecalbins, nullptr );	

	// Initilize the histograms and Canvases
	for ( int i = 0; i < n_ecalbins; i++ )
	{
		for ( int j = 0; j < n_gemlayers; j++ )
		{
			std::string histName = "h2dhit_ecalbin_" + std::to_string(i) + "_gemlayer_" + std::to_string(j);

			if ( j <= 5 )
			{
				gemhithists_perecalbin[i][j] = new TH2F(histName.c_str(), histName.c_str(), 12,-0.3, 0.3, 34, -0.85, 0.85); 
			}

			else
			{
				gemhithists_perecalbin[i][j] = new TH2F(histName.c_str(), histName.c_str(), 16,-0.4, 0.4, 44, -1.1, 1.1); 
			}				
		}

		std::string canvasName = "ecalbin_" + std::to_string(i);

		canvas_perecalbin[i] = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 1000);
		//canvas_perecalbin[i]->Divide(4,2);
	}

	
	long event_no {0};

	// First loop to fill the histograms.
	while ( T->GetEntry(event_no++) )
	{
		// Some Track quality cuts on the simulation data.
		if ( n_gemtracks == 0 || gemtrack_mid->at(0) != 0 || gemtrack_numhits->at(0) < 3 ) continue;

		// Check if the ECal arrays have data
		if ( ecalhit_row->empty() || ecalhit_col->empty() || ecalhit_cell->empty() || ecalhit_sumedep->empty() || ecalhit_xcell->empty() || ecalhit_ycell->empty() ) continue;
			
		// Check if the array sizes match
		if ( ecalhit_row->size() != n_ecalhits || ecalhit_col->size() != n_ecalhits || ecalhit_cell->size() != n_ecalhits || ecalhit_sumedep->size() != n_ecalhits || ecalhit_xcell->size() != n_ecalhits || ecalhit_ycell->size() != n_ecalhits ) continue;

		// See if the GEM array sizes match.
		if ( gemhit_x->size() != n_gemhits || gemhit_y->size() != n_gemhits || gemhit_plane->size() != n_gemhits ) continue;

		// Do ECal clustering.
		ecalbin.clearECalClusVec();
		ecalbin.ECalClustering(n_ecalhits, ecalhit_row, ecalhit_col, ecalhit_cell, ecalhit_sumedep, ecalhit_xcell, ecalhit_ycell);

		// Skip this event if no ECal clusters exists.
		if ( !ecalbin.does_ecalclus_exist() ) continue;

		int ecal_bin_no = ecalbin.highestE_ecalbinno();

		//std::cout << "Ev no: " << event_no << " : ECal bin: " << ecal_bin_no << '\n';

		// Fill the GEM 2D hit histograms.
		for ( int ihit = 0; ihit < n_gemhits; ihit++ )
		{
			//gemhitlayer_hists[ gemhit_plane->at(ihit) - 1 ]->Fill( gemhit_y->at(ihit), gemhit_x->at(ihit) );
			gemhithists_perecalbin[ecal_bin_no][gemhit_plane->at(ihit) - 1]->Fill( gemhit_y->at(ihit), gemhit_x->at(ihit) + g4sbsFT_xoffset );
		}

	}

	// Make a vector of vectors with ROI.
	std::vector<std::vector<ROI>> vvROI(n_ecalbins, std::vector<ROI>(n_gemlayers, ROI(0,0,0,0)));

	// Loop through the histograms and make ROIs for each of them.
	for ( int i = 0; i < n_ecalbins; i++ )
	{
		for ( int j = 0; j < n_gemlayers; j++ )
		{
			//vvROI[i][j] = computeROI( gemhithists_perecalbin[i][j], 0.95 );
			vvROI[i][j] = computeROI_MaxBinExpand_2( gemhithists_perecalbin[i][j], 0.1, 0.01 );
			std::cout << "ECal bin " << i << " GEM layer " << j << " ROI, x_min, x_max, y_min, y_max:" << vvROI[i][j].xMin << vvROI[i][j].xMax << vvROI[i][j].yMin << vvROI[i][j].yMax << '\n';
		}
	}

	for ( int i_can = 0; i_can < n_ecalbins; i_can++ )
	{
		canvas_perecalbin[i_can]->cd();
		// canvas_perecalbin[i_can]->SetGrid();
		// canvas_perecalbin[i_can]->DrawFrame(-0.6, -1.5, 0.7, 1.5, "ECal Detector;X Position (m);Y Position (m)");

		// Create a large left pad manually
    	TPad *leftPad = new TPad("leftPad", "Left Pad", 0.0, 0.0, 0.33, 1.0); // Merging first column
    	//leftPad->cd();
    	leftPad->Draw();
    	// // // leftPad->cd();
    	// leftPad->SetGrid();
    	//gemhithists_perecalbin[i_can][1]->Draw("COLZ");
    	leftPad->cd();

    	int iECalBin = i_can;
	    leftPad->DrawFrame(-0.6, -1.5, 0.7, 1.5, Form("ECal Bin: %i;X Position (m);Y Position (m)", iECalBin));

	    TList* ecalListPerBin = ecalbin.returnECalDiagramPerBinTListPtr( iECalBin );

	 	ecalListPerBin->Draw("same");
   		// //leftPad->cd();
   		// leftPad->Update();

   		 // Create the remaining 8 pads in a 4×2 grid (excluding the leftmost column)
	    canvas_perecalbin[i_can]->cd();
    	TPad* pad[8];
    	int padIndex = 0;
	    for (int i = 0; i < 2; i++) 
	    { // Two rows
	        for (int j = 1; j < 5; j++) 
	        { // Skipping the first column
	            double xlow = 0.33 + 0.16 * (j - 1); // Start from 0.2 since leftPad occupies 0-0.2
	            double ylow = 0.5 * (1 - i);
	            double xup = xlow + 0.16;
	            double yup = ylow + 0.5;
	            pad[padIndex] = new TPad(Form("pad%d", padIndex + 1), Form("Pad %d", padIndex + 1), xlow, ylow, xup, yup);
	            pad[padIndex]->Draw();
	            padIndex++;
	        }
	    }

		for ( int i_gemlayer = 0; i_gemlayer < n_gemlayers; i_gemlayer++ )
		{
			// canvas_perecalbin[i_can]->cd(i_gemlayer+1);
			pad[i_gemlayer]->cd();
			pad[i_gemlayer]->SetGrid();

			gemhithists_perecalbin[i_can][i_gemlayer]->Draw("COLZ");

			TBox* roiBox = new TBox(vvROI[i_can][i_gemlayer].xMin, vvROI[i_can][i_gemlayer].yMin, vvROI[i_can][i_gemlayer].xMax, vvROI[i_can][i_gemlayer].yMax);
			roiBox->SetLineColor(kRed);
			roiBox->SetLineWidth(2);
			roiBox->SetFillStyle(0);
			roiBox->Draw("same");
		}

		canvas_perecalbin[i_can]->cd();
		canvas_perecalbin[i_can]->Update();

		if ( i_can == 0 ) canvas_perecalbin[i_can]->Print("output_gemhits_perecalbin_wROI_1perc_wTrackCuts_NEWECalBinning.pdf(");
		else if ( i_can == n_ecalbins - 1 ) canvas_perecalbin[i_can]->Print("output_gemhits_perecalbin_wROI_1perc_wTrackCuts_NEWECalBinning.pdf)");
		else canvas_perecalbin[i_can]->Print("output_gemhits_perecalbin_wROI_1perc_wTrackCuts_NEWECalBinning.pdf");
	}

	writeROIToFile( vvROI, "ROI_GEP3_FT_1.txt" );

}