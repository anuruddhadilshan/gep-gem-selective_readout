// This class will read the local DB file and construct all the arrays and maps to be used in the main/final analysis.

#ifndef DBREAD_H
#define DBREAD_H



class DBread
{

private: 


double m_m0_uangle, m_m0_vangle, m_m1_uangle, m_m1_vangle, m_m2_uangle, m_m2_vangle, m_m3_uangle, m_m3_vangle, m_m4_uangle, m_m4_vangle, m_m5_uangle, m_m5_vangle, m_m6_uangle, m_m6_vangle;
double m_m2_uoffset, m_m2_voffset;
double m_m2_nstripsu, m_m2_nstripsv;
double m_m2_pos_x, m_m2_pos_y, m_m2_pos_z;
double m_m2_size_x, m_m2_size_y, m_m2_size_z;
double m_upithc, m_vpithc;

std::string m_dbFileName;
std::ifstream m_dbFile;

public:

	DBread( const std::string& db_local = "db_FT_local.dat" )
	: m_dbFileName { db_local }, m_dbFile{ db_local }
	{

		if ( readDBfile() == -1) std::cerr << "Database not read properly. STOP and debug!!!\n";

	}

	int readDBfile()
	{
		if ( !m_dbFile.is_open() ) 
		{
			std::cerr << "ERROR: Could not open the DB file " << m_dbFileName << endl;

			return -1;
		}

		TString currentline;

		while ( currentline.ReadLine( m_dbFile ) )
		{
			// Ignore comment lines or empty lines
			if ( !currentline.BeginsWith("#") )
			{
				TObjArray *tokens = currentline.Tokenize(" ");

				int ntokens = tokens->GetEntries();

				if ( ntokens >1 )
				{
					TString skey = ( (TObjString*)(*tokens)[0] )->GetString();

					if ( skey == "m2.uangle" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_m2_uangle = sval.Atoi();
						std::cout << "sbs.gemFT.m2.uangle " << m_m2_uangle << endl;
					}

					if ( skey == "m2.vangle" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_m2_vangle = sval.Atoi();
						std::cout << "sbs.gemFT.m2.vangle " << m_m2_vangle << endl;
					}

					if ( skey == "m2.uoffset" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_m2_uoffset = sval.Atof();
						std::cout << "sbs.gemFT.m2.uoffset " << m_m2_uoffset << endl;
					}

					if ( skey == "m2.voffset" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_m2_voffset = sval.Atof();
						std::cout << "sbs.gemFT.m2.voffset " << m_m2_voffset << endl;
					}

					if ( skey == "m2.nstripsu" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_m2_nstripsu = sval.Atoi();
						std::cout << "sbs.gemFT.m2.nstripsu " << m_m2_nstripsu << endl;
					}

					if ( skey == "m2.nstripsv" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_m2_nstripsv = sval.Atoi();
						std::cout << "sbs.gemFT.m2.nstripsv " << m_m2_nstripsv << endl;
					}

					if ( skey == "m2.position" )
					{
						TString sval1 = ( (TObjString*)(*tokens)[1] )->GetString();
						m_m2_pos_x = sval1.Atof();

						TString sval2 = ( (TObjString*)(*tokens)[2] )->GetString();
						m_m2_pos_y = sval2.Atof();

						TString sval3 = ( (TObjString*)(*tokens)[3] )->GetString();
						m_m2_pos_z = sval3.Atof();

						std::cout << "sbs.gemFT.m2.position " << m_m2_pos_x << " " << m_m2_pos_y << " " << m_m2_pos_z << endl;
					}

					if ( skey == "m2.size" )
					{
						TString sval1 = ( (TObjString*)(*tokens)[1] )->GetString();
						m_m2_size_x = sval1.Atof();

						TString sval2 = ( (TObjString*)(*tokens)[2] )->GetString();
						m_m2_size_y = sval2.Atof();

						TString sval3 = ( (TObjString*)(*tokens)[3] )->GetString();
						m_m2_size_z = sval3.Atof();

						std::cout << "sbs.gemFT.m2.size " << m_m2_size_x << " " << m_m2_size_y << " " << m_m2_size_z << endl;
					}

					if ( skey == "upithc" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_upithc = sval.Atoi();
						std::cout << "sbs.gemFT.upithc " << m_upithc << endl;
					}

					if ( skey == "vpithc" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_vpithc = sval.Atoi();
						std::cout << "sbs.gemFT.vpithc " << m_vpithc << endl;
					}

				}
			}

		}

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

};

#endif