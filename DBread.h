// This class will read the local DB file and construct all the arrays and maps to be used in the main/final analysis.

#ifndef DBREAD_H
#define DBREAD_H

struct GEMModule
{
	int modNum;
	int modLayer;
	double xPos, yPos, zPos;
	double xSize, ySize, zSize;
	double xMin, xMax, yMin, yMax;
	double uAngle, vAngle;
	double uOffset, vOffset;
	int nStripsu, nStripsv;
	double pitch = 0.0004; // We know this will not change for ANY of the SBS GEMs so let's hard-code it.
};

struct GEMLayer
{
	int layerNum;
	std::map < int, GEMModule > gemModulesLayer;
};


class DBread
{

private: 

int m_nlayers = 0;
int m_nmodules = 0;

int m_dbfile_readstatus;

double m_m0_uangle, m_m0_vangle, m_m1_uangle, m_m1_vangle, m_m2_uangle, m_m2_vangle, m_m3_uangle, m_m3_vangle, m_m4_uangle, m_m4_vangle, m_m5_uangle, m_m5_vangle, m_m6_uangle, m_m6_vangle;
// double m_m2_uoffset, m_m2_voffset;
// double m_m2_nstripsu, m_m2_nstripsv;
// double m_m2_pos_x, m_m2_pos_y, m_m2_pos_z;
// double m_m2_size_x, m_m2_size_y, m_m2_size_z;
double m_upithc, m_vpithc;

std::map < int, GEMModule > m_GEMModules;
std::map < int, GEMLayer > m_GEMLayers; 

std::string m_dbFileName;
std::ifstream m_dbFile;

public:

	DBread( const std::string& db_local = "db_FT_local.dat" )
	: m_dbFileName { db_local }, m_dbFile{ db_local }
	{
		std::cout << endl;
		std::cout << "*** Reading local GEM database: " << db_local << " ***" << endl;

		if ( readDBfile() == -1) 
		{
			m_dbfile_readstatus = -1;

			std::cerr << "Database not read properly. STOP and debug!!!\n";
		}
		else
		{
			m_dbfile_readstatus = 0;

			std::cout << "*** Finished reading local GEM database ***" << endl;
		}

		
		std::cout << endl;
	}

	int readDBfile()
	{
		if ( !m_dbFile.is_open() ) 
		{
			std::cerr << "ERROR: Could not open the DB file " << m_dbFileName << endl;

			return -1;
		}

		TString currentline1;

		while ( currentline1.ReadLine( m_dbFile ) )
		{
			// Ignore comment lines or empty lines
			if ( !currentline1.BeginsWith("#") )
			{
				TObjArray *tokens = currentline1.Tokenize(" ");

				int ntokens = tokens->GetEntries();

				if ( ntokens > 1 )
				{
					TString skey = ( (TObjString*)(*tokens)[0] )->GetString();

					if ( skey == "nlayers" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_nlayers = sval.Atoi();
						std::cout << "Number of Tracking Layers: " << m_nlayers << endl;
					}

					if ( skey == "nmodules" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_nmodules = sval.Atoi();
						std::cout << "Number of GEM modules " << m_nmodules << endl;
					}
				}

				delete tokens;
			}			
		}

		// Resetting the file pointer to go to the beginning of the file.
		m_dbFile.clear();                 // Clear EOF flag
		m_dbFile.seekg(0, std::ios::beg); // Seek back to beginning

		if ( m_nlayers == 0 || m_nmodules == 0 )
		{
			std::cerr << "EROOR: Number of tracking layers OR GEM modules not defined in the DB or set to 0 in the DB!!!!" << endl;
			return -1;
		}

		
		for ( int iMod = 0; iMod < m_nmodules; iMod++ )
		{
			TString currentline2;

			GEMModule thisMod;
			thisMod.modNum = iMod;

			while ( currentline2.ReadLine( m_dbFile ) )
			{	
				// Ignore comment lines or empty lines
				if ( !currentline2.BeginsWith("#") )
				{
					TObjArray *tokens = currentline2.Tokenize(" ");

					int ntokens = tokens->GetEntries();

					if ( ntokens > 1 )
					{
						TString skey = ( (TObjString*)(*tokens)[0] )->GetString();

						if ( skey == Form("m%i.layer", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.modLayer = sval.Atoi();
							std::cout << Form("m%i.layer: ", iMod) << thisMod.modLayer << endl;
						}

						if ( skey == Form("m%i.uangle", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.uAngle = sval.Atof();
							std::cout << Form("m%i.uangle: ", iMod) << thisMod.uAngle << endl;
						}

						if ( skey == Form("m%i.vangle", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.vAngle = sval.Atof();
							std::cout << Form("m%i.vangle: ", iMod) << thisMod.vAngle << endl;
						}

						if ( skey == Form("m%i.uoffset", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.uOffset = sval.Atof();
							std::cout << Form("m%i.uoffset: ", iMod) << thisMod.uOffset << endl;
						}

						if ( skey == Form("m%i.voffset", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.vOffset = sval.Atof();
							std::cout << Form("m%i.voffset: ", iMod) << thisMod.vOffset << endl;
						}

						if ( skey == Form("m%i.nstripsu", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.nStripsu = sval.Atoi();
							std::cout << Form("m%i.nstripsu: ", iMod) << thisMod.nStripsu << endl;
						}

						if ( skey == Form("m%i.nstripsv", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.nStripsv = sval.Atoi();
							std::cout << Form("m%i.nstripsv: ", iMod) << thisMod.nStripsv << endl;
						}

						if ( skey == Form("m%i.position", iMod) )
						{
							TString sval1 = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.xPos = sval1.Atof();
							TString sval2 = ( (TObjString*)(*tokens)[2] )->GetString();
							thisMod.yPos = sval2.Atof();
							TString sval3 = ( (TObjString*)(*tokens)[3] )->GetString();
							thisMod.zPos = sval3.Atof();

							std::cout << Form("m%i.position: ", iMod) << thisMod.xPos << " " << thisMod.yPos << " " << thisMod.zPos << endl;
						}

						if ( skey == Form("m%i.size", iMod) )
						{
							TString sval1 = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.xSize = sval1.Atof();
							TString sval2 = ( (TObjString*)(*tokens)[2] )->GetString();
							thisMod.ySize = sval2.Atof();
							TString sval3 = ( (TObjString*)(*tokens)[3] )->GetString();
							thisMod.zSize = sval3.Atof();

							std::cout << Form("m%i.size: ", iMod) << thisMod.xSize << " " << thisMod.ySize << " " << thisMod.zSize << endl;
						}

					}

					delete tokens;

				}
			}

			// Resetting the file pointer to go to the beginning of the file.
			m_dbFile.clear();                 // Clear EOF flag
			m_dbFile.seekg(0, std::ios::beg); // Seek back to beginnings

			thisMod.xMin = - thisMod.xSize/2; //thisMod.xPos - thisMod.xSize/2;
			thisMod.xMax = thisMod.xSize/2; //thisMod.xPos + thisMod.xSize/2;
			thisMod.yMin = - thisMod.ySize/2; //thisMod.yPos - thisMod.ySize/2;
			thisMod.yMax = thisMod.ySize/2; //thisMod.yPos + thisMod.ySize/2;

			m_GEMModules.emplace( iMod, thisMod );

		}
		
		for ( int iLayer = 0; iLayer < m_nlayers; iLayer++ )
		{
			GEMLayer thisLayer;
			thisLayer.layerNum = iLayer;

			for ( const auto& [modNum, mod] : m_GEMModules )
			{
				if ( mod.modLayer == iLayer )
				{
					(thisLayer.gemModulesLayer).emplace( mod.modNum, mod );
				}
			}

			m_GEMLayers.emplace( iLayer, thisLayer );

			if ( (thisLayer.gemModulesLayer).size() > 0 ) std::cout << "* Layer: " << iLayer << " constructed." << endl;
		}

		// TString currentline;

		// while ( currentline.ReadLine( m_dbFile ) )
		// {
		// 	// Ignore comment lines or empty lines
		// 	if ( !currentline.BeginsWith("#") )
		// 	{
		// 		TObjArray *tokens = currentline.Tokenize(" ");

		// 		int ntokens = tokens->GetEntries();

		// 		if ( ntokens > 1 )
		// 		{
		// 			TString skey = ( (TObjString*)(*tokens)[0] )->GetString();

		// 			if ( skey == "m2.uangle" )
		// 			{
		// 				TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_m2_uangle = sval.Atoi();
		// 				std::cout << "sbs.gemFT.m2.uangle " << m_m2_uangle << endl;
		// 			}

		// 			if ( skey == "m2.vangle" )
		// 			{
		// 				TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_m2_vangle = sval.Atoi();
		// 				std::cout << "sbs.gemFT.m2.vangle " << m_m2_vangle << endl;
		// 			}

		// 			if ( skey == "m2.uoffset" )
		// 			{
		// 				TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_m2_uoffset = sval.Atof();
		// 				std::cout << "sbs.gemFT.m2.uoffset " << m_m2_uoffset << endl;
		// 			}

		// 			if ( skey == "m2.voffset" )
		// 			{
		// 				TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_m2_voffset = sval.Atof();
		// 				std::cout << "sbs.gemFT.m2.voffset " << m_m2_voffset << endl;
		// 			}

		// 			if ( skey == "m2.nstripsu" )
		// 			{
		// 				TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_m2_nstripsu = sval.Atoi();
		// 				std::cout << "sbs.gemFT.m2.nstripsu " << m_m2_nstripsu << endl;
		// 			}

		// 			if ( skey == "m2.nstripsv" )
		// 			{
		// 				TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_m2_nstripsv = sval.Atoi();
		// 				std::cout << "sbs.gemFT.m2.nstripsv " << m_m2_nstripsv << endl;
		// 			}

		// 			if ( skey == "m2.position" )
		// 			{
		// 				TString sval1 = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_m2_pos_x = sval1.Atof();

		// 				TString sval2 = ( (TObjString*)(*tokens)[2] )->GetString();
		// 				m_m2_pos_y = sval2.Atof();

		// 				TString sval3 = ( (TObjString*)(*tokens)[3] )->GetString();
		// 				m_m2_pos_z = sval3.Atof();

		// 				std::cout << "sbs.gemFT.m2.position " << m_m2_pos_x << " " << m_m2_pos_y << " " << m_m2_pos_z << endl;
		// 			}

		// 			if ( skey == "m2.size" )
		// 			{
		// 				TString sval1 = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_m2_size_x = sval1.Atof();

		// 				TString sval2 = ( (TObjString*)(*tokens)[2] )->GetString();
		// 				m_m2_size_y = sval2.Atof();

		// 				TString sval3 = ( (TObjString*)(*tokens)[3] )->GetString();
		// 				m_m2_size_z = sval3.Atof();

		// 				std::cout << "sbs.gemFT.m2.size " << m_m2_size_x << " " << m_m2_size_y << " " << m_m2_size_z << endl;
		// 			}

		// 			if ( skey == "upithc" )
		// 			{
		// 				TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_upithc = sval.Atoi();
		// 				std::cout << "sbs.gemFT.upithc " << m_upithc << endl;
		// 			}

		// 			if ( skey == "vpithc" )
		// 			{
		// 				TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
		// 				m_vpithc = sval.Atoi();
		// 				std::cout << "sbs.gemFT.vpithc " << m_vpithc << endl;
		// 			}

		// 		}
		// 	}

		// }

		return 0;
	}

	
	int return_uang()
	{
		return m_m2_uangle;
	}	

	int return_vang()
	{
		return m_m2_vangle;
	}

	int returnFileReadStatus()
	{
		return m_dbfile_readstatus;
	}

	std::map< int, GEMLayer >& returnGEMLayerMap() { return m_GEMLayers; }

};

#endif