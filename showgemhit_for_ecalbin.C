// Show the GEM hits on layers for a given ECal bin.

#include "ecalbin.h"




void showgemhit_for_ecalbin( const char* simfilename )
{

	
	//for ( auto bin_no : list_of_bins_ecal ) std::cout << bin_no << '\n';

	// Continue with reading the TTree and making the plots.
	// TFile* simrootfile = new TFile(simfilename, "READ");

	// TTree* T = (TTree*)simrootfile->Get("T");

	// Initialize the ECalbin class.
	Ecalbin ecalbin{};

	TChain* T = new TChain("T");

	T->Add("/lustre24/expphy/volatile/halla/sbs/nhunt/gep_elas_crystal/GEP3_crystal_align_ecal_job_5*.root");

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
	std::vector<double>* gemhit_y = nullptr;

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

	std::vector<TH2F*> gemhitlayer_hists;

	for (int i = 0; i < 8; i++)
	{
		gemhitlayer_hists.push_back( new TH2F( Form("h2Dhit_layer_%d", i+1), Form("2D hits in GEM layer: %d", i+1), 100,-0.5, 0.5, 100, -1.5, 1.5) );
	}

	TCanvas* C = new TCanvas("canvas", "Hit Locations in Detector Planes", 1200, 800);
	C->Divide(4,2);

	std::string usrinput;

	while ( true )
	{
		std::cout << "Enter the ECal bin number (or 'q' to quit): ";
		std::cin >> usrinput;

		// Check for quit condition
		if (usrinput == "q")
		{
			std::cout << "'q' entered. Exiting the program. \n";
			break;
		}

		// Parse ECal bin number
		int ecalbin_no;

		try
		{
			ecalbin_no = std::stoi(usrinput);
		}
		catch ( const std::invalid_argument& )
		{
			std::cerr << "Invalid input. Please enter a valid bin number or 'q' to quit.\n";
			continue;
		}

		// See if the bin# entered is within the list of bins in ECal.
		if ( ecalbin.doesECalBinNotExists(ecalbin_no) )
		{
			std::cerr << "The ECal bin number enered is not within the list of bins." << '\n';
			continue;
		}

		
		long event_no {0};
		
		long no_of_ecalclus_inbin{0};

		// Reset GEM histos.
		for ( int i = 0; i < 8; i++ ) gemhitlayer_hists[i]->Reset();	

		// Event loop
		while ( T->GetEntry(event_no++) )
		{
			// Check if the ECal arrays have data
			if ( ecalhit_row->empty() || ecalhit_col->empty() || ecalhit_cell->empty() || ecalhit_sumedep->empty() || ecalhit_xcell->empty() || ecalhit_ycell->empty() ) continue;
			
			// Check if the array sizes match
			if ( ecalhit_row->size() != n_ecalhits || ecalhit_col->size() != n_ecalhits || ecalhit_cell->size() != n_ecalhits || ecalhit_sumedep->size() != n_ecalhits || ecalhit_xcell->size() != n_ecalhits || ecalhit_ycell->size() != n_ecalhits ) continue;
			
			// Do ECal clustering.
			ecalbin.clearECalClusVec();
			ecalbin.ECalClustering(n_ecalhits, ecalhit_row, ecalhit_col, ecalhit_cell, ecalhit_sumedep, ecalhit_xcell, ecalhit_ycell);

    	// Check if there is an ECal cluster.
    	if ( !ecalbin.does_ecalclus_exist() ) continue;
    		
    	// Check if the current ECal cluster belongs to the bin number entered by the user.
    	if ( ecalbin.highestE_ecalbinno() != ecalbin_no ) continue;

			++no_of_ecalclus_inbin;

			// See if all the GEM arrays have data.
			if ( gemhit_x->empty() || gemhit_x->empty() || gemhit_plane->empty() ) continue;

			// See if the GEM array sizes match.
			if ( gemhit_x->size() != n_gemhits || gemhit_y->size() != n_gemhits || gemhit_plane->size() != n_gemhits ) continue;

			// Fill the GEM 2D hit histograms.
			for ( int ihit = 0; ihit < n_gemhits; ihit++ )
			{
				gemhitlayer_hists[ gemhit_plane->at(ihit) - 1 ]->Fill( gemhit_y->at(ihit), gemhit_x->at(ihit) );
			}

			std::cout << "ECal clus no: " << no_of_ecalclus_inbin << '\n';
		}

		std::cout << "Total number of ECal clusters in the bin: " << no_of_ecalclus_inbin << '\n';

		for ( int i = 0; i < 8; i++ )
		{
			C->cd(i+1);
			gemhitlayer_hists[i]->Draw("COLZ");
		}
	
		C->Update();

	}
	

}