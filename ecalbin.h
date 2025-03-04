#include <fstream>
#include <cmath>
#include "TVector3.h"

#ifndef ECALBIN_h
#define ECALBIN_h

//using namespace std;

struct calblock
{
	int cell;
	int row;
	int col;
	double x;
	double y;
	double E; //smeared energy
	double Etrue; //true energy
	double nphe; //estimated nphe
};

struct calcluster
{
	int node; // "node" index
	int cell; // cell index of center cell
	int bin; //bin number of cluster = binx + biny * nbinsx
	int binx, biny;
	double Esum_true;
	double nphesum;
	double Esum; //cluster energy sum
	double Esum_norm; //cluster energy sum / elastic peak energy
};

class Ecalbin
{

private:

	//const double Rcal = 4.9; //m, to back surface
	const double Mp = 0.938272; //GeV
	const double PI = TMath::Pi();

	const double xspace = 0.04292; // ECal cell width in m.
	const double yspace = 0.04293; // ECal cell height in m.

	double Ebeam_default = 10.688;
 	double ECALtheta_default = 29.75*TMath::DegToRad();
 	double ECALdist_default = 4.7;

 	double normfac_ECAL=0.88;

	TVector3 ecal_zaxis{std::sin(ECALtheta_default),0,std::cos(ECALtheta_default)};
	TVector3 ecal_yaxis{0,1,0};
	TVector3 ecal_xaxis = ecal_yaxis.Cross(ecal_zaxis).Unit();
	TVector3 ecal_pos = ECALdist_default * ecal_zaxis;

	//Make every cell a cluster center:
  
	std::set<int> list_of_nodes_ecal;
	std::map<int, int> cluster_centers_ecal_by_node; //Record the center cell of each cluster!
	std::map<int, int> nodes_ecal_by_cluster_center; //key = cell, mapped value = node
	std::map<int, set<int> > cells_logic_sums_ecal; //mapping between node numbers and cell numbers
	std::map<int, double> logic_mean_ecal; //mean peak positions by node number
	std::map<int, double> logic_sigma_ecal; //peak width by node number
	std::map<int, double> threshold_ecal; //threshold by node number
	std::map<std::pair<int,int>, int > cell_rowcol_ecal; //cell numbers mapped by unique row and column pairs
	std::map<int,std::set<int> > nodes_cells_ecal; //mapping of nodes by cell number:

	std::map<int,double> Eprime_expect_by_node_ecal;

	int ncells_ecal;
	std::set<int> cells_ecal;      //If cells are somehow not in order: 
	std::map<int,int> rows_cells_ecal;
	std::map<int,int> cols_cells_ecal;
	  
	std::map<int,double> xcell_ecal;
	std::map<int,double> ycell_ecal;
	std::map<int,int> nearestcol_up_ecal;
	std::map<int,int> nearestcol_down_ecal;
	//  map<int,int> bins_cells_ecal; //This is for later
	std::set<int> list_of_bins_ecal; //list of bins containing at least one cluster center
	  
	std::map<int,int> binx_by_node_ecal;
	std::map<int,int> biny_by_node_ecal;
	std::map<int,int> binglobal_by_node_ecal;

	std::map<int,int> binx_by_cell_ecal;
	std::map<int,int> biny_by_cell_ecal;
	std::map<int,int> binglobal_by_cell_ecal;

	std::map<int,std::set<int> > cells_by_binglobal_ecal;
	std::map<int,std::set<int> > nodes_by_binglobal_ecal;

	//keep track of min and max x by row number:
	double ycellmin = 1.0e9;
	double ycellmax = -1.0e9;
	double xcellmin = 1.0e9;
	double xcellmax = -1.0e9;
	std::map<int,double> ycell_rows;
	//map<int,double> cellsize_rows;
	std::map<int,double> xcellmin_rows;
	std::map<int,double> xcellmax_rows;
	  
	int minrow=1000, maxrow=-1;
	 
	std::set<int> rows_ecal;
	std::map<int,std::set<int> > columns_rows_ecal;
	  
	std::map<int,double> elastic_peak_new_ecal;
	std::map<int,double> sigma_new_ecal;
	std::map<int,double> threshold_new_ecal;

	std::ifstream mapfile_ecal;
	  
	//Cell runs from 0 to N-1

	TRandom3 num{0};

	vector<calcluster> ECALclusters;

	TList* ecalList = new TList();

	TList* ecalListPerBin = new TList();


public:

	Ecalbin( const std::string& ecal_mapfile = "/w/halla-scshelf2102/sbs/adr/g4sbs/g4sbs/database/ecal_gep_blockmap.txt" ) : mapfile_ecal(ecal_mapfile)
	{
		TString currentline;

		//Get cell mapping info:
		while( currentline.ReadLine( mapfile_ecal ) )
		{
			//    cout << currentline << endl;
		    if( !currentline.BeginsWith("#") )
		    {
		      TObjArray *tokens = currentline.Tokenize(",");
		      int ntokens = tokens->GetEntries();

		      if( ntokens == 5 )
		      {
				TString scell = ( (TObjString*) (*tokens)[0] )->GetString();
				TString srow  = ( (TObjString*) (*tokens)[1] )->GetString();
				TString scol  = ( (TObjString*) (*tokens)[2] )->GetString();
				TString sxcell = ( (TObjString*) (*tokens)[3] )->GetString();
				TString sycell = ( (TObjString*) (*tokens)[4] )->GetString();

				int cellnum = scell.Atoi();
				int row = srow.Atoi();
				int col = scol.Atoi();
				double xcell = sxcell.Atof();
				double ycell = sycell.Atof();

				std::pair<int, int> rowcoltemp(row,col);

				rows_ecal.insert(row);
				columns_rows_ecal[row].insert(col); //if we're lazy we can ASSUME that column within a row runs from 0 to ncolumns - 1.
						
				cells_ecal.insert( cellnum );
				cell_rowcol_ecal[rowcoltemp] = cellnum;
				rows_cells_ecal[cellnum] = row;
				cols_cells_ecal[cellnum] = col;
				xcell_ecal[cellnum] = xcell/100.0;
				ycell_ecal[cellnum] = ycell/100.0;

				if( ycell_rows.empty() || ycell/100.0 < ycellmin ) ycellmin = ycell/100.0;
				if( ycell_rows.empty() || ycell/100.0 > ycellmax ) ycellmax = ycell/100.0;
						  
				ycell_rows[row] = ycell/100.0;

				maxrow = (row > maxrow) ? row : maxrow;
				minrow = (row < minrow) ? row : minrow;
				
				// if( row < 30 ) {
				//   cellsize_rows[row] = 4.0; //cm
				// } else {
				//   cellsize_rows[row] = 4.2; //cm
				// }
						
				if( xcellmin_rows.empty() || xcell/100.0 < xcellmin_rows[row] )
				{
				  xcellmin_rows[row] = xcell/100.0;
				}

				if( xcellmax_rows.empty() || xcell/100.0 > xcellmax_rows[row] )
				{
				  xcellmax_rows[row] = xcell/100.0;
				}

				xcellmin = xcell/100.0 < xcellmin ? xcell/100.0 : xcellmin;
				xcellmax = xcell/100.0 > xcellmax ? xcell/100.0 : xcellmax;		
		      }
		    }
		}

		//we want to superimpose a rectangular binning on the ECAL clusters:
		ncells_ecal = cells_ecal.size();
		//The SAFEST approach to forming the cluster groupings is to go by rows and columns:
		for( auto row : rows_ecal )
		{
			for( auto col : columns_rows_ecal[row] )
		    {
		      std::pair<int,int> rowcoltemp( row, col );
		      int cell = cell_rowcol_ecal[rowcoltemp];

		      //grab position, although I think we are mainly interested in x (horizontal) here:
		      double xcell = xcell_ecal[cell];
		      double ycell = ycell_ecal[cell];
		      
		      if( row+1 < rows_ecal.size() )
		      { //if the row above exists, find nearest column:
				int nearestcol=-1;
				double mindiff = 1000.0;
				for( auto colup : columns_rows_ecal[row+1] )
				{
				  std::pair<int,int> rowcolup(row+1,colup);
				  int cellup = cell_rowcol_ecal[rowcolup];
				  double xup = xcell_ecal[cellup];

				  bool goodx = fabs(xup-xcell)<0.5*0.043 && fabs(xup-xcell)<mindiff;
			  
				  if( goodx )
				  {
				    nearestcol = colup;
				    mindiff = fabs(xup-xcell);
				  }
				}

				nearestcol_up_ecal[cell] = nearestcol;
		      } 
		      else 
		      {
			  	nearestcol_up_ecal[cell] = -1;
		      }

		      if( row > 0 )
		      { //if the row below exists, find nearest column below:
			  	int nearestcol = -1;
				double mindiff = 1000.0;
				for( auto coldown : columns_rows_ecal[row-1] )
				{
				  std::pair<int,int> rowcoldown(row-1,coldown);
				  int celldown = cell_rowcol_ecal[rowcoldown];
				  double xdown = xcell_ecal[celldown];
				  bool goodx = fabs(xdown-xcell)<0.5*0.043 && fabs(xdown-xcell) < mindiff;
			 	  if( goodx )
			 	  {
				  	nearestcol = coldown;
				    mindiff = abs(xdown-xcell);
			  	  }
				}
				nearestcol_down_ecal[cell] = nearestcol;
		      } 
		      else 
		      {
				nearestcol_down_ecal[cell] = -1;
		      }	
		    }
		}

		int node = 0;

		
		//We impose a regular rectangular binning on ECAL to account for its non-rectangular shape
		//since the minimum x and y cell positions occur for EDGE blocks; let's set the low edge of the grid at xmin + block spacing/2 etc.
		double binxmin = xcellmin + 0.5*xspace;
		double binymin = ycellmin + 0.5*yspace;

		int nbinsx = 0, nbinsy = 0;
		while( binxmin + 3.0*xspace*nbinsx < xcellmax ) nbinsx++;
		while( binymin + 3.0*yspace*nbinsy < ycellmax ) nbinsy++; 

		double binxmax = binxmin + 3.0*xspace*nbinsx;
		double binymax = binymin + 3.0*yspace*nbinsy;
		  
		//  double binxmax = xcellmax + 0.5*xspace;
		//double binymax = ycellmax + 0.5*yspace;

		//this rounds to the nearest integer number of bins and adds 1
		//int nbinsx = int( (binxmax-binxmin)/(2.0*xspace) + 0.5) + 1;
		//int nbinsy = int( (binymax-binymin)/(2.0*yspace) + 0.5) + 1;

		int nbinstot = nbinsx*nbinsy;
		  
		//now implement 3x5 groupings and binning:
		for( auto cell : cells_ecal )
		{
			int row = rows_cells_ecal[cell];
		    int col = cols_cells_ecal[cell];

		    double xtemp = xcell_ecal[cell];
		    double ytemp = ycell_ecal[cell];

		    //int binxtemp = int( (xtemp - xcellmin + 0.5*xspace)/(2.0*xspace) );
		    
		    if( row > 0 && row + 1 < rows_ecal.size() && 	col > 0 && col + 1 < columns_rows_ecal[row].size() )
		    { //not an edge block; make a cluster!
		      list_of_nodes_ecal.insert(node); //"node" starts at ZERO and goes to Nnodes - 1
		      cluster_centers_ecal_by_node[node] = cell;
		      nodes_ecal_by_cluster_center[cell] = node;

		      //here we are calculating the global position of the cell center coordinate in the g4sbs coordinate system:
		      TVector3 clusterpos_global = ecal_pos + xtemp * ecal_xaxis + ytemp * ecal_yaxis;

		      //Calculate the polar scattering angle for this cluster center assuming a point target at the origin
		      double clustertheta = acos( clusterpos_global.Z() / clusterpos_global.Mag() );

		      //Calculate the expected scattered electron energy for this scattering angle:
		      double clusterEprime = Ebeam_default/(1.0 + Ebeam_default/Mp * (1.0-cos(clustertheta)));

		      //store expected scattered electron energy for this cluster center:
		      Eprime_expect_by_node_ecal[node] = clusterEprime; 

		      //figure out the bins:
		      int binxtemp = int( (xtemp - binxmin)/(3.0*xspace) );
		      int binytemp = int( (ytemp - binymin)/(3.0*yspace) );
		      int binglobaltemp = binxtemp + binytemp * nbinsx;
		      if( binxtemp >= nbinsx ) cout << "warning: cluster center with bin x >= nbinsx, this will screw up global bin counting" << endl;

		      list_of_bins_ecal.insert( binglobaltemp );

		      binx_by_node_ecal[node] = binxtemp;
		      biny_by_node_ecal[node] = binytemp;
		      binglobal_by_node_ecal[node] = binglobaltemp;
		      binx_by_cell_ecal[cell] = binxtemp;
		      biny_by_cell_ecal[cell] = binytemp;
		      binglobal_by_cell_ecal[cell] = binglobaltemp;

		      logic_mean_ecal[node] = Eprime_expect_by_node_ecal[node]*normfac_ECAL;
		      logic_sigma_ecal[node] = 0.06*logic_mean_ecal[node];
		      
		      cells_by_binglobal_ecal[binglobaltemp].insert( cell );
		      nodes_by_binglobal_ecal[binglobaltemp].insert( node );
		      
		      //cells_logic_sums_ecal[node].insert( cell );
		      
		      for( int rowj=row-2; rowj<=row+2; rowj++ )
		      {
				int colmin=col-1, colmax=col+1;
				if( rowj < row )
				{
				  colmin = nearestcol_down_ecal[cell]-1;
				  colmax = colmin + 2;
				} 
				else if( rowj > row )
				{
				  colmin = nearestcol_up_ecal[cell]-1;
				  colmax = colmin + 2;
				}

				//check valid row and column indices:
				if( rowj >= 0 && rowj < rows_ecal.size() )
				{ 
			 		for( int colj = colmin; colj <= colmax; colj++ )
			 		{
			    		if( colj >= 0 && colj < columns_rows_ecal[rowj].size() )
			    		{
					    	std::pair<int,int> rowcolj(rowj,colj);
					    	auto testcell = cell_rowcol_ecal.find(rowcolj);
							if( testcell != cell_rowcol_ecal.end() )
					      	{
								cells_logic_sums_ecal[node].insert( testcell->second );
								nodes_cells_ecal[cell].insert(node);
					      	}
			    		}
			  		}
				}
		       }
		      
		      node++;
			}	
		}
	}

	static bool CompareCalBlocks( const calblock &b1, const calblock &b2 )
	{
  		return b1.E > b2.E;
	}

	static bool CompareCalClusters( const calcluster &c1, const calcluster &c2 )
	{
  		return c1.Esum > c2.Esum; 
	}


	void ECalClustering( int n_ecalhits, std::vector<int>* ecalhit_row, std::vector<int>* ecalhit_col, std::vector<int>* ecalhit_cell, std::vector<double>* ecalhit_sumedep, std::vector<double>* ecalhit_xcell, std::vector<double>* ecalhit_ycell )
	{
		// Loop to fill the ECal block vectors
		double nphe {0.};
		vector<calblock> ECALblocks_all;
		vector<calblock> ECALblocks_unused;

		for ( int ihit = 0; ihit < n_ecalhits; ihit++ )
		{
			int rowhit = ecalhit_row->at(ihit);
			int colhit = ecalhit_col->at(ihit);
			std::pair<int,int> rowcolhit( rowhit, colhit );

			int cellhit = ecalhit_cell->at(ihit);

			double edep = ecalhit_sumedep->at(ihit);

			double mean = 528.*edep; //from C16_GEANT_report.pdf, 528 phe/GeV

			double nphe = num.Gaus(mean, sqrt(mean));

			double esmear = nphe/528.;

			calblock thisblock;
			thisblock.cell = cellhit;
			thisblock.row = rowhit;
			thisblock.col = colhit;
			thisblock.x = ecalhit_xcell->at(ihit);
			thisblock.y = ecalhit_ycell->at(ihit);
			thisblock.E = esmear;
			thisblock.Etrue = edep;
			thisblock.nphe = nphe;

			ECALblocks_all.push_back( thisblock );
			ECALblocks_unused.push_back( thisblock );

			//g_ECalClus->SetPoint( g_ECalClus->GetN(), thisblock.y, thisblock.x );
			//std::cout << thisblock.y << " " << thisblock.x << '\n';
		}

		while ( !ECALblocks_unused.empty() )
		{
			// Sort remaining ECal blocks by energy
			std::sort( ECALblocks_unused.begin(), ECALblocks_unused.end(), CompareCalBlocks ); //This should sort in descending order by energy

			auto maxblk = ECALblocks_unused.begin();

			//Only form a cluster if the current local maximum is in the list of possible cluster centers; i.e., if it is not a edge block.
			if( nodes_ecal_by_cluster_center.find((*maxblk).cell) != nodes_ecal_by_cluster_center.end() )
			{
				calcluster thiscluster;
				thiscluster.cell = (*maxblk).cell;
				thiscluster.node = nodes_ecal_by_cluster_center[(*maxblk).cell];
				thiscluster.bin = binglobal_by_cell_ecal[(*maxblk).cell];
				thiscluster.binx = binx_by_cell_ecal[(*maxblk).cell];
				thiscluster.biny = biny_by_cell_ecal[(*maxblk).cell];
				thiscluster.Esum = 0.0;
				thiscluster.Esum_true = 0.0;
				thiscluster.nphesum = 0.0;
					  
				//Next we need to loop on all blocks, add up the energies from any blocks in this cluster, and erase unused blocks in this cluster
				// so they won't also form clusters:
				ECALblocks_unused.erase( maxblk );

				set<int> blockstoerase;
				for( int iblk=0; iblk<ECALblocks_all.size(); iblk++ )
				{
			    	calblock blk_i = ECALblocks_all[iblk];
			    	if( cells_logic_sums_ecal[thiscluster.node].find( blk_i.cell ) != cells_logic_sums_ecal[thiscluster.node].end() )
			    	{
				      blockstoerase.insert(blk_i.cell);
				      //don't double-count the energy of the seed:
				      thiscluster.Esum += blk_i.E;
				      thiscluster.Esum_true += blk_i.Etrue;
				      thiscluster.nphesum += blk_i.nphe;
			    	}
			  	}

			  	auto itblk = ECALblocks_unused.begin();

			  	while( itblk != ECALblocks_unused.end() )
			  	{
			    	if( blockstoerase.find((*itblk).cell) != blockstoerase.end() )
			    	{
			     		itblk = ECALblocks_unused.erase( itblk );
			    	}
			    	else 
			    	{
			      		++itblk;
			    	}
			  	}
				thiscluster.Esum_norm = thiscluster.Esum / logic_mean_ecal[thiscluster.node];
				ECALclusters.push_back( thiscluster );
			} 
			else 
			{ 
				//std::cout << "Removing blk from list...\n";
				// current maximum is an edge block; erase without forming a cluster:
		  		ECALblocks_unused.erase( maxblk );
			}

		} // done with ECal clustering loop.

		//Sort ECAL clusters by energy (existing comparison operator is based on smeared energy sum, NOT normalized):
    	std::sort( ECALclusters.begin(), ECALclusters.end(), CompareCalClusters );
	}

	bool does_ecalclus_exist()
	{
		return ECALclusters.size() > 0;
	}

	int highestE_ecalbinno()
	{
		return (ECALclusters.at(0)).bin;
	}

	bool doesECalBinNotExists( const int input_ecalbin_no )
	{
		return list_of_bins_ecal.find(input_ecalbin_no) == list_of_bins_ecal.end();
	}

	void clearECalClusVec()
	{
		ECALclusters.clear();
	}

	int numOfECalBins()
	{
		return cells_by_binglobal_ecal.size();
	}

	std::map< int, std::set<int> >& returnCellsByBinGlobalMap()
	{
		return cells_by_binglobal_ecal;
	}

	void drawECal2DCellDiagram()
	{
		// Create a ROOT Canvas
    	TCanvas *c1 = new TCanvas("c1", "ECal Detector Front Face", 800, 1600);
	    c1->SetGrid();
	    c1->DrawFrame(-0.6, -1.5, 0.7, 1.5, "ECal Detector;X Position (m);Y Position (m)");

	    // Loop through the cells and draw them as TBoxes
	    for (const auto &cell : xcell_ecal) {
	        int cellnum = cell.first;
	        double x_center = cell.second;
	        double y_center = ycell_ecal[cellnum];

	        // Calculate boundaries
	        double x_min = (x_center - xspace / 2.0);
	        double x_max = (x_center + xspace / 2.0);
	        double y_min = y_center - yspace / 2.0;
	        double y_max = y_center + yspace / 2.0;

	        // Create a TBox to represent the ECal cell
	        TBox *cellBox = new TBox(x_min, y_min, x_max, y_max);
	        cellBox->SetLineColor(kBlack);
	        cellBox->SetLineWidth(1);
	        cellBox->SetFillStyle(0); // Transparent fill
	        cellBox->Draw("same");
	    }

	    // Loop through the nodes and draw a point where the node is at.

	  //   for ( const auto & cellbynode : cluster_centers_ecal_by_node )
	  //   {
	  //   	//int nodenum = cellbynode.first;
	  //   	int cellnum = cellbynode.second;
	  //   	double x_node = xcell_ecal[cellnum];
	  //   	double y_node = ycell_ecal[cellnum];

	  //   	TMarker *marker = new TMarker(x_node, y_node, 20); // 20 is a filled circle
			// marker->SetMarkerSize(1);
			// marker->SetMarkerColor(kBlue);
			// marker->Draw("same");

	  //   }

	    for ( const auto & cellsinbin : cells_by_binglobal_ecal )
	    {
	    	for ( const auto & cell : cells_by_binglobal_ecal[cellsinbin.first] )
	    	{
	    		TString text = TString::Format("%i", cellsinbin.first);
			    TText *t = new TText(xcell_ecal[cell], ycell_ecal[cell], text);
			    t->SetTextSize(0.01);
			    t->SetTextColor(kBlue);
			    t->Draw("same");
	    	}
	    }

	    // Update canvas
	    c1->Update();
	}

	TList* returnECalDiagramTListPtr()
	{
		// Clear list if it was previously used.
		ecalList->Clear();

		// Loop through the cells and draw them as TBoxes
	    for (const auto &cell : xcell_ecal) 
	    {
	        int cellnum = cell.first;
	        double x_center = cell.second;
	        double y_center = ycell_ecal[cellnum];

	        // Calculate boundaries
	        double x_min = (x_center - xspace / 2.0);
	        double x_max = (x_center + xspace / 2.0);
	        double y_min = y_center - yspace / 2.0;
	        double y_max = y_center + yspace / 2.0;

	        // Create a TBox to represent the ECal cell
	        TBox *cellBox = new TBox(x_min, y_min, x_max, y_max);
	        cellBox->SetLineColor(kBlack);
	        cellBox->SetLineWidth(1);
	        cellBox->SetFillStyle(0); // Transparent fill
	        
	        ecalList->Add(cellBox); // Store in the list.
	    }

	    return ecalList;
	}

	TList* returnECalDiagramPerBinTListPtr( const int binNum )
	{
		// Clear list if it was previously used.
		ecalListPerBin->Clear();

		if ( list_of_bins_ecal.find(binNum) == list_of_bins_ecal.end() ) // See whether the input ECal bin number exists.
	    {
	    	std::cerr << "The ECal bin number: " << binNum << " is not a valid ECal bin number.\n";

	    	return ecalListPerBin;   	
	    }

		// Loop through the cells and draw them as TBoxes
	    for (const auto &cell : xcell_ecal) 
	    {
	        int cellnum = cell.first;
	        double x_center = cell.second;
	        double y_center = ycell_ecal[cellnum];

	        // Calculate boundaries
	        double x_min = (x_center - xspace / 2.0);
	        double x_max = (x_center + xspace / 2.0);
	        double y_min = y_center - yspace / 2.0;
	        double y_max = y_center + yspace / 2.0;

	        // Create a TBox to represent the ECal cell
	        TBox *cellBox = new TBox(x_min, y_min, x_max, y_max);
	        cellBox->SetLineColor(kBlack);
	        cellBox->SetLineWidth(1);
	        cellBox->SetFillStyle(0); // Transparent fill

	        if ( binglobal_by_cell_ecal[cellnum] == binNum ) // See if the current cell belongs to the input ECal bin number.
	        {
	        	// cellBox->SetLineColor(kRed);
	        	// cellBox->SetLineWidth(3);
	        	// cellBox->SetFillStyle(4100);
	        	// cellBox->SetFillColor(kBlue);

	        	TString text = TString::Format("%i", cellnum);
	        	const double text_Xpos = xcell_ecal[cellnum] - xspace/2 ;
	        	const double text_Ypos =  ycell_ecal[cellnum] - yspace/4;
				TText *t = new TText(text_Xpos, text_Ypos, text);
				t->SetTextSize(0.01);
				t->SetTextColor(kBlue);

				ecalListPerBin->Add(t);
	        }
	        
	        ecalListPerBin->Add(cellBox); // Store in the list.
	    }

	    return ecalListPerBin;

	}

	std::map< int, double >& returnXcellsEcalMap()
	{
		return xcell_ecal;
	}

	std::map< int, double >& returnYcellsEcalMap()
	{
		return ycell_ecal;
	}

	std::set < int >& returnListOfBinsECalMap()
	{
		return list_of_bins_ecal;
	}

};

#endif