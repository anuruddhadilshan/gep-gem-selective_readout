#ifndef ROIREAD_H
#define ROIREAD_H


// Define the struct to hold the Region Of Interest (ROI) rectangle boundaries.
struct ROI
{
	double xMin, xMax, yMin, yMax;

	// Constructor requiring parameters
    ROI(double x_min, double x_max, double y_min, double y_max) 
        : xMin(x_min), xMax(x_max), yMin(y_min), yMax(y_max) {}
};


class ROIread
{

private:

	const std::string m_roi_filename;
	std::ifstream m_roi_file;

	int m_roifile_readstatus;

	std::map< int, std::map< int, ROI > > m_roi;


public:

	ROIread( const std::string& roi_file )
	: m_roi_filename{roi_file}, m_roi_file{roi_file}
	{
		std::cout << endl;
		std::cout << "*** Reading and copying ROIs for each GEM layers, for each ECal bin, from file: " << m_roi_filename << " ***" << endl;
		
		if ( readROIfile() == -1 ) 
		{
			m_roifile_readstatus = -1;

			std::cerr << "ROI file not read properly. STOP and debug!!!" << endl;
		}
		else
		{
			m_roifile_readstatus = 0;

			std::cout << "*** Finished reading the ROI file ***" << endl;
		}

		std::cout << endl;
		
	}

	int readROIfile()
	{
		if ( !m_roi_file.is_open() ) 
		{
			std::cerr << "ERROR: Could not open the ROI file " << m_roi_filename << endl;

			return -1;
		}

		TString currentline;

		while ( currentline.ReadLine(m_roi_file) )
		{
			// Ignore comment lines or empty lines
		 	if ( !currentline.BeginsWith("#") )
		 	{
		 		TObjArray *tokens = currentline.Tokenize(" ");

		 		int ntokens = tokens->GetEntries();

		 		if ( ntokens == 6 )
		 		{
		 			TString sbin = ( (TObjString*)(*tokens)[0] )->GetString();
		 			TString slayer = ( (TObjString*)(*tokens)[1] )->GetString();
		 			TString sxmin = ( (TObjString*)(*tokens)[2] )->GetString();
		 			TString sxmax = ( (TObjString*)(*tokens)[3] )->GetString();
		 			TString symin = ( (TObjString*)(*tokens)[4] )->GetString();
		 			TString symax = ( (TObjString*)(*tokens)[5] )->GetString();

		 			int bin = sbin.Atoi();
		 			int layer = slayer.Atoi();
		 			double xmin = sxmin.Atof();
		 			double xmax = sxmax.Atof();
		 			double ymin = symin.Atof();
		 			double ymax = symax.Atof();

		 			ROI thisROI{ xmin, xmax, ymin, ymax };

		 			m_roi[bin].emplace( layer, thisROI );		 			
		 		}
		 	}

		}

		std::cout << "Bin Layer xMin xMax yMin yMax" << endl;

		for ( const auto& [binNum, layer_roi] : m_roi )
		{
			for ( const auto& [layerNum, roi] : layer_roi )
			{
				std::cout << binNum << " " << layerNum << " " << roi.xMin << " " << roi.xMax << " " << roi.yMin << " " << roi.yMax << endl;
			}
		}

		return 0;
	}

	int returnFileReadStatus()
	{
		return m_roifile_readstatus;
	}

	std::map< int, std::map< int, ROI > >& return_ROIMap() { return m_roi; }

};
	
#endif