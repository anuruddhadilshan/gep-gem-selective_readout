// This script will generate GEM FT ROIs by analyzing replayed real/simulation data.

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

//vtpcrate, fiber, adc_ch
struct apvInfoVals{
	//fill these values they id the APV
	int vtpcrate;
	int fiber; //NOTE: == mpd_id
	int adc_ch;
	// int invert;

	bool operator<(const apvInfoVals& other) const {
        if (vtpcrate != other.vtpcrate) return vtpcrate < other.vtpcrate;
        if (fiber != other.fiber) return fiber < other.fiber;
        return adc_ch < other.adc_ch;
    }

    void print() const {
        std::cout << "VTPcrate: " << vtpcrate 
                  << " Fiber: " << fiber 
                  << " ADC_ch: " << adc_ch 
                  << std::endl; 
	}
};



constexpr int max_NECalClus {100};
constexpr int max_NGEMhits {10000};


void FT_ROI_fromRealData()
{

	TChain* C = new TChain("T");

	C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30.root");
	C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_1.root");
	C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_2.root");
	C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_3.root");
	C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_4.root");
	C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_5.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_1.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_2.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_3.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_4.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_5.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_6.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_7.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_8.root");
	// C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse_9.root");
	//C->Add("/lustre24/expphy/volatile/halla/sbs/adr/gep-parsed/GEP1/LH2/bigprod_Apr30/GEP1_LH2_bigbrod_Apr30_May6parse.root");
	//C->Add("/lustre24/expphy/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/bigprod_Apr28/rootfiles/gep5_fullreplay_3320_stream2_2_seg9_9.root");
	// std::cout << "Begin making the gloabl cut..." << std::endl;
	// C->Draw(">>elist", 
	// 	"(((((abs(heep.dxECAL-0.03)<0.06&&abs(heep.dyECAL-0.023)<0.08)&&(abs(heep.dt_ADC-115)<20.0&&abs(heep.dpp+0.12)<0.1))&&(sbs.gemFT.track.chi2ndf[0]<30))&&((sbs.gemFT.track.nhits[0]>4||sbs.gemFT.track.chi2ndf[0]<15)))&&(abs(sbs.tr.x[0]+sbs.tr.th[0]*sbs.z_bcp[30]-sbs.x_bcp[30]+0.076)<0.21))&&(abs(sbs.tr.y[0]+sbs.tr.ph[0]*sbs.z_bcp[30]-sbs.y_bcp[30])<0.21)"
	// 	, "entrylist");
	// TEntryList *elist = (TEntryList*)gDirectory->Get("elist");

	// check how many entries passed the cut
	std::cout << "Events for analysis: " << C->GetEntries() << std::endl;

	// C->SetEntryList(elist);

	TH2F* h2_FT_layer [8];

	for ( int i = 0; i < 8; i++ ) h2_FT_layer[i] = new TH2F( Form("h2_FT_layer_%i", i), Form("h2_FT_layer_%i; Y (m); X (m)", i), 100,-0.4, 0.4, 100, -1.1, 1.1 );

	C->SetBranchStatus("*", 0);
	C->SetBranchStatus("Ndata.earm.ecal.clus_blk.id", 1);
	C->SetBranchStatus("earm.ecal.clus_blk.id", 1);
	//C->SetBranchStatus("Ndata.earm.ecal.clus.e", 1);
	C->SetBranchStatus("earm.ecal.clus_blk.e", 1);
	C->SetBranchStatus("Ndata.sbs.gemFT.hit.layer", 1);
	C->SetBranchStatus("sbs.gemFT.hit.layer", 1);
	//C->SetBranchStatus("Ndata.sbs.gemFT.module", 1);
	//C->SetBranchStatus("sbs.gemFT.module", 1);
	//C->SetBranchStatus("Ndata.sbs.gemFT.hit.xglobal", 1);
	C->SetBranchStatus("sbs.gemFT.hit.xglobal", 1);
	//C->SetBranchStatus("Ndata.sbs.gemFT.hit.yglobal", 1);
	C->SetBranchStatus("sbs.gemFT.hit.yglobal", 1);
	C->SetBranchStatus("sbs.gemFT.hit.trackindex", 1);

	int nECalClusBlks;
	double ecalclusblk_id [max_NECalClus];
	double ecalclusblk_e [max_NECalClus];
	int nGemHits;
	double gemhit_layer [max_NGEMhits];
	//double gemhit_module [max_NGEMhits];
	double gemhit_xglobal [max_NGEMhits];
	double gemhit_yglobal [max_NGEMhits];
	double gemhit_trackindex [max_NGEMhits];

	C->SetBranchAddress("Ndata.earm.ecal.clus_blk.id", &nECalClusBlks);
	C->SetBranchAddress("earm.ecal.clus_blk.id", ecalclusblk_id);
	C->SetBranchAddress("earm.ecal.clus_blk.e", ecalclusblk_e);
	C->SetBranchAddress("Ndata.sbs.gemFT.hit.layer", &nGemHits);
	C->SetBranchAddress("sbs.gemFT.hit.layer", gemhit_layer);
	//C->SetBranchAddress("sbs.gemFT.hit.module", gemhit_module);
	C->SetBranchAddress("sbs.gemFT.hit.xglobal", gemhit_xglobal);
	C->SetBranchAddress("sbs.gemFT.hit.yglobal", gemhit_yglobal); 
	C->SetBranchAddress("sbs.gemFT.hit.trackindex", gemhit_trackindex);

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

		canvas_perecalbin[i] = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
	}

	long nevent = 0;

	while ( C->GetEntry(nevent++) )
	{
		int maxCellOfMaxCluster = int( ecalclusblk_id[0] ); 
		int ecalBinNo_thisEvent = ecalbin.ecalBinByCellNo_New(maxCellOfMaxCluster);

		for ( int ihit = 0; ihit < nGemHits; ihit++ )
		{
			if ( gemhit_trackindex[ihit] != 0 ) continue; // Only plot hits in the best track.

			gemhithists_perecalbin[ecalBinNo_thisEvent][gemhit_layer[ihit]]->Fill( gemhit_yglobal[ihit], gemhit_xglobal[ihit] );		

		}

	}

	// TCanvas* canvas = new TCanvas("hists","hists",1200,800);
	// canvas->Divide(4,2);

	// for ( int i = 0; i < 8; i++ )
	// {
	// 	canvas->cd(i+1);

	// 	h2_FT_layer[i]->Draw("COLZ");
	// }
	
	// canvas->SaveAs("FT_GEP1_allhisforlayers.pdf");

	// Make a vector of vectors with ROI.
	std::vector<std::vector<ROI>> vvROI(n_ecalbins, std::vector<ROI>(n_gemlayers, ROI(0,0,0,0)));

	// Loop through the histograms and make ROIs for each of them.
	for ( int i = 0; i < n_ecalbins; i++ )
	{
		for ( int j = 0; j < n_gemlayers; j++ )
		{
			//vvROI[i][j] = computeROI( gemhithists_perecalbin[i][j], 0.95 );
			vvROI[i][j] = computeROI_MaxBinExpand_2( gemhithists_perecalbin[i][j], 0.1, 0.0001 );
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

		if ( i_can == 0 ) canvas_perecalbin[i_can]->Print("realdat_GEP1.pdf(");
		else if ( i_can == n_ecalbins - 1 ) canvas_perecalbin[i_can]->Print("realdat_GEP1.pdf)");
		else canvas_perecalbin[i_can]->Print("realdat_GEP1.pdf");
	}

	writeROIToFile( vvROI, "ROI_GEP1_FT_realdat.txt" );


	// Now let's see what VTP, MPD, and ADC are getting hits within these ROIs for each ECal bin.
	std::map < int /*ECalBinNo*/, std::set<apvInfoVals>> ECalBinAPVinfo;

	C->SetBranchStatus("*", 0);
	C->SetBranchStatus("sbs.gemFT.hit.trackindex", 1);
	C->SetBranchStatus("sbs.gemFT.hit.crate_U", 1);
	C->SetBranchStatus("sbs.gemFT.hit.mod_U", 1);
	C->SetBranchStatus("sbs.gemFT.hit.adc_id_U", 1);
	

}