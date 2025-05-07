// This class will read the local DB file and construct all the arrays and maps to be used in the main/final analysis.

#ifndef DBREAD_H
#define DBREAD_H


//just adding these for sake of getting rid of the correspondoning error highlighting in vscode
//can comment out or delete later
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>

#include "TString.h"
#include "TObjString.h"


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


//module id, axis, position on axis
struct apvInfoKeys {
	int gemid, axis, pos;

	void print() const {
	std::cout 
	<< "GEMId: " << gemid 
	<< " Axis: " << axis 
	<< " Pos: " << pos 
	<< std::endl;
	};
	bool operator<(const apvInfoKeys& other) const {
	if (gemid != other.gemid) return gemid < other.gemid;
	if (axis != other.axis) return axis < other.axis;
	return pos < other.pos;
	}
};

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

std::map<apvInfoKeys, apvInfoVals> apvInfoMap; 

// Function to write `apvInfoMap` to a file
void OutputAPVinfoMap() {
std::ofstream mapFile("APV_Map_TEST.txt");
if (!mapFile.is_open()) {
	std::cerr << "Error: Could not create file APV_Map_TEST.txt\n";
	return;
}

mapFile << "#apvInfoMap actually maps key and value structs to eachother\n" ;
mapFile << "## the apvInfoKeys consists of: module id, axis(U/V depending on map), and  pos(along axis)\n";
mapFile << "## the apvInfoVals consists of: VTPcrate, Fiber(MPD) ID, and  the ADC channel\n\n";

// Set column widths for formatting
int colWidth = 10; // Adjust as needed for alignment

mapFile << std::left << std::setw(colWidth) << "GEMId"
		<< std::setw(colWidth) << "Axis"
		<< std::setw(colWidth) << "Pos"
		<< std::setw(colWidth) << "VTPcrate"
		<< std::setw(colWidth) << "Fiber"
		<< std::setw(colWidth) << "ADC_ch"
		<< "\n";

mapFile << std::string(6 * colWidth, '-') << "\n"; // Create a line separator

for (const auto& [key, val] : apvInfoMap) {
	mapFile << std::left << std::setw(colWidth) << key.gemid
			<< std::setw(colWidth) << key.axis
			<< std::setw(colWidth) << key.pos
			<< std::setw(colWidth) << val.vtpcrate
			<< std::setw(colWidth) << val.fiber
			<< std::setw(colWidth) << val.adc_ch
			<< "\n";
}

std::cout << "Parsed data written to APV_Map_TEST.txt\n";
}	





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
		std::cout << std::endl;
		std::cout << "*** Reading local GEM database: " << db_local << " ***" << std::endl;

		if ( readDBfile() == -1) 
		{
			m_dbfile_readstatus = -1;

			std::cerr << "Database not read properly. STOP and debug!!!\n";
		}
		else
		{
			m_dbfile_readstatus = 0;

			std::cout << "*** Finished reading local GEM database ***" << std::endl;
		}

		
		std::cout << std::endl;
	}

	int readDBfile()
	{
		if ( !m_dbFile.is_open() ) 
		{
			std::cerr << "ERROR: Could not open the DB file " << m_dbFileName << std::endl;

			return -1;
		}

		std::cout << "DEBUG: First 5 lines of DB file:\n";
TString debugLine;
int debugCounter = 0;
while (debugLine.ReadLine(m_dbFile) && debugCounter < 5) {
    std::cout << debugLine << std::endl;
    debugCounter++;
}
m_dbFile.clear();
m_dbFile.seekg(0, std::ios::beg);



		TString currentline1;
		int currentM = -1;  // Current module ID used for getting chanMap info

bool inAPVSection = false;  // Flag to detect APV map section

		while ( currentline1.ReadLine( m_dbFile ) )
		{
			
			
			// std::cout << "\n\n currentline1 is: " << 
			// currentline1 << std::endl; //RREDIT



			// Ignore comment lines or empty lines
			if ( !currentline1.BeginsWith("#") 
			// ||currentline1.BeginsWith("##")//marks start of each module's chanMap data
			)
			{
				TObjArray *tokens = currentline1.Tokenize(" ");

				int ntokens = tokens->GetEntries();

				if ( ntokens > 1 )
				{
					TString skey = ( (TObjString*)(*tokens)[0] )->GetString();

					// std::cout << "\n\nskey is " << skey << std::endl;

					if ( skey == "nlayers" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_nlayers = sval.Atoi();
						std::cout << "Number of Tracking Layers: " << m_nlayers << std::endl;
					}

					if ( skey == "nmodules" )
					{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_nmodules = sval.Atoi();
						std::cout << "Number of GEM modules " << m_nmodules << std::endl;
					}

					// Detect start of APV mapping section
					if (skey.Contains("GEMId") && ntokens >= 6) {
						inAPVSection = true;  // We are now inside the APV mapping section
						delete tokens;
						continue;  // Skip this line (header row)
					}
			
					if (inAPVSection) {
						if (skey.Contains("-")) {  
							delete tokens;
							continue;  // Skip separator line (dashes)
						}
			
						// Process APV mappings (Only numeric values should be here)
						if (ntokens >= 6) {
							apvInfoKeys currAPVkeys;
							apvInfoVals currAPVvals;
			
							currAPVkeys.gemid = ((TObjString*)(*tokens)[0])->GetString().Atoi();
							currAPVkeys.axis = ((TObjString*)(*tokens)[1])->GetString().Atoi();
							currAPVkeys.pos = ((TObjString*)(*tokens)[2])->GetString().Atoi();
							currAPVvals.vtpcrate = ((TObjString*)(*tokens)[3])->GetString().Atoi();
							currAPVvals.fiber = ((TObjString*)(*tokens)[4])->GetString().Atoi();
							currAPVvals.adc_ch = ((TObjString*)(*tokens)[5])->GetString().Atoi();
			
							// Store in map
							apvInfoMap[currAPVkeys] = currAPVvals;
			
							// Debug Output
							// std::cout << "APV Mapping: Module " << currAPVkeys.gemid
							// 		  << ", Axis " << currAPVkeys.axis
							// 		  << ", Pos " << currAPVkeys.pos
							// 		  << ", VTPcrate " << currAPVvals.vtpcrate
							// 		  << ", Fiber " << currAPVvals.fiber
							// 		  << ", ADC_ch " << currAPVvals.adc_ch
							// 		  << std::endl;
						}
				}

				delete tokens;
			}			
		}
	}

		// Resetting the file pointer to go to the beginning of the file.
		m_dbFile.clear();                 // Clear EOF flag
		m_dbFile.seekg(0, std::ios::beg); // Seek back to beginning

		if ( m_nlayers == 0 || m_nmodules == 0 )
		{
			std::cerr << "EROOR: Number of tracking layers OR GEM modules not defined in the DB or set to 0 in the DB!!!!" << std::endl;
			return -1;
		}

		// **Processes GEM Modules**
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
							std::cout << Form("m%i.layer: ", iMod) << thisMod.modLayer << std::endl;
						}

						if ( skey == Form("m%i.uangle", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.uAngle = sval.Atof();
							std::cout << Form("m%i.uangle: ", iMod) << thisMod.uAngle << std::endl;
						}

						if ( skey == Form("m%i.vangle", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.vAngle = sval.Atof();
							std::cout << Form("m%i.vangle: ", iMod) << thisMod.vAngle << std::endl;
						}

						if ( skey == Form("m%i.uoffset", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.uOffset = sval.Atof();
							std::cout << Form("m%i.uoffset: ", iMod) << thisMod.uOffset << std::endl;
						}

						if ( skey == Form("m%i.voffset", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.vOffset = sval.Atof();
							std::cout << Form("m%i.voffset: ", iMod) << thisMod.vOffset << std::endl;
						}

						if ( skey == Form("m%i.nstripsu", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.nStripsu = sval.Atoi();
							std::cout << Form("m%i.nstripsu: ", iMod) << thisMod.nStripsu << std::endl;
						}

						if ( skey == Form("m%i.nstripsv", iMod) )
						{
							TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.nStripsv = sval.Atoi();
							std::cout << Form("m%i.nstripsv: ", iMod) << thisMod.nStripsv << std::endl;
						}

						if ( skey == Form("m%i.position", iMod) )
						{
							TString sval1 = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.xPos = sval1.Atof();
							TString sval2 = ( (TObjString*)(*tokens)[2] )->GetString();
							thisMod.yPos = sval2.Atof();
							TString sval3 = ( (TObjString*)(*tokens)[3] )->GetString();
							thisMod.zPos = sval3.Atof();

							std::cout << Form("m%i.position: ", iMod) << thisMod.xPos << " " << thisMod.yPos << " " << thisMod.zPos << std::endl;
						}

						if ( skey == Form("m%i.size", iMod) )
						{
							TString sval1 = ( (TObjString*)(*tokens)[1] )->GetString();
							thisMod.xSize = sval1.Atof();
							TString sval2 = ( (TObjString*)(*tokens)[2] )->GetString();
							thisMod.ySize = sval2.Atof();
							TString sval3 = ( (TObjString*)(*tokens)[3] )->GetString();
							thisMod.zSize = sval3.Atof();

							std::cout << Form("m%i.size: ", iMod) << thisMod.xSize << " " << thisMod.ySize << " " << thisMod.zSize << std::endl;
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

			if ( (thisLayer.gemModulesLayer).size() > 0 ) std::cout << "* Layer: " << iLayer << " constructed." << std::endl;
		}

		std::cout << "Finished filling apvInfoMap with " << apvInfoMap.size() << " entries.\n";

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


	std::map <apvInfoKeys, apvInfoVals> GetAPVinfoMap(){
		
		return apvInfoMap;
	}
	
};
#endif