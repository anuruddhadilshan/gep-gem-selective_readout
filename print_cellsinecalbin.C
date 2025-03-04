#include "ecalbin2.h"

void print_cellsinecalbin()
{
	// Initialize the ECalbin class.
	Ecalbin ecalbin{};

	std::map< int, std::set<int> > cells_by_binglobal_ecal_NEW = ecalbin.returnCellsByBinGlobalNewMap();
	std::map< int, double > xcell_ecal = ecalbin.returnXcellsEcalMap();
	std::map< int, double > ycell_ecal = ecalbin.returnYcellsEcalMap();
	std::set< int > list_of_bins_ecal_NEW = ecalbin.returnListOfBinsECalNewMap();

	int nECalBins = cells_by_binglobal_ecal_NEW.size(); //ecalbin.numOfECalBins();

	std::cout << "Number of ECal bins: " << nECalBins << endl;

	TCanvas* c[nECalBins];

	//TList* ecalList = ecalbin.returnECalDiagramTListPtr();
	TList* ecalListPerBin;

	int i = 0;


	// for ( const auto & binAndcells : cells_by_binglobal_ecal )
	// {
	// 	c[i] = new TCanvas( Form("c%i", i), Form("ECal Bin %i", binAndcells.first), 500, 800 );
	// 	c[i]->cd();
	// 	c[i]->DrawFrame(-0.6, -1.5, 0.7, 1.5, Form("ECal Bin: %i;X Position (m);Y Position (m)", binAndcells.first));

	// 	ecalList->Draw("same");

	// 	for ( const auto & cell : cells_by_binglobal_ecal[binAndcells.first] )
	//     {
	//     	TString text = TString::Format("%i", binAndcells.first);
	// 		TText *t = new TText(xcell_ecal[cell], ycell_ecal[cell], text);
	// 		t->SetTextSize(0.01);
	// 		t->SetTextColor(kBlue);
	// 		t->Draw("same");
	//     }

	// 	c[i]->Update();

	// 	if ( i == 0 ) c[i]->Print("ecalbin_display.pdf(");
	// 	else if ( i == nECalBins - 1 ) c[i]->Print("ecalbin_display.pdf)");
	// 	else c[i]->Print("ecalbin_display.pdf");

	// 	i++;
	
	// }

	for ( const auto& binNum : list_of_bins_ecal_NEW )
	{
		c[i] = new TCanvas( Form("c%i", i), Form("ECal Bin %i", binNum), 500, 800 );
	 	c[i]->cd();
	 	c[i]->DrawFrame(-0.6, -1.5, 0.7, 1.5, Form("ECal Bin: %i;X Position (m);Y Position (m)", binNum));

	 	ecalListPerBin = ecalbin.returnECalDiagramPerBinTListPtr( binNum );

	 	ecalListPerBin->Draw("same");

	 	c[i]->Update();

	 	if ( i == 0 ) c[i]->Print("ecalbin_display2_newbinning.pdf(");
	 	else if ( i == nECalBins - 1 ) c[i]->Print("ecalbin_display2_newbinning.pdf)");
	 	else c[i]->Print("ecalbin_display2_newbinning.pdf");
	 	i++;
	}


}