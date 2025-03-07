#include "DBread.h"

double fPxU, fPyU, fPxV, fPyV;

std::map < int, std::vector<double> > layer_mod_bounds;

// std::vector<double> roi {0.1, 0.45, -0.2, 0.2};
std::vector<double> roi {0.25, 0.65, -0.25, 0.15};

void xylayer()
{
	layer_mod_bounds.emplace(0, std::vector<double>{ 
    (-0.766 - 0.512 / 2), 
    (-0.766 + 0.512 / 2), 
    (-0.6144 / 2), 
    (0.6144 / 2) 
	});

	layer_mod_bounds.emplace(1, std::vector<double>{ 
    (-0.256 - 0.512 / 2), 
    (-0.256 + 0.512 / 2), 
    (-0.6144 / 2), 
    (0.6144 / 2) 
	});

	layer_mod_bounds.emplace(2, std::vector<double>{ 
    (0.256 - 0.512 / 2), 
    (0.256 + 0.512 / 2), 
    (-0.6144 / 2), 
    (0.6144 / 2) 
	});

	layer_mod_bounds.emplace(3, std::vector<double>{ 
    (0.766 - 0.512 / 2), 
    (0.766 + 0.512 / 2), 
    (-0.6144 / 2), 
    (0.6144 / 2) 
	});

	TVector2 XY_ROI_1 { roi.at(1), roi.at(3) };
  	TVector2 XY_ROI_2 { roi.at(1), roi.at(2) };
  	TVector2 XY_ROI_3 { roi.at(0), roi.at(2) };
  	TVector2 XY_ROI_4 { roi.at(0), roi.at(3) };

  	std::set <int> roi_modules;

	// First we should check within what module numbers within the layer, our ROI is spanning across.
	for ( const auto& [modNum, posVec] : layer_mod_bounds )
	{
		// if ( ( roi.at(0) >= posVec.at(0) || roi.at(1) <= posVec.at(1) ) ) // && roi.at(2) >= posVec.at(2) && roi.at(3) <= posVec.at(3) )
		// {
		// 	std::cout << "ROI is WITHIN module " << modNum << endl; 
		// }
		// else
		// {
		// 	std::cout << "ROI is NOT within module " << modNum << endl;
		// }

		if ( roi.at(0) >= posVec.at(0) && roi.at(0) <= posVec.at(1) && roi.at(2) >= posVec.at(2) && roi.at(3) <= posVec.at(3) )
		{
			std::cout << "ROI X min is within module " << modNum << endl;
			roi_modules.insert(modNum);
		}

		if ( roi.at(1) >= posVec.at(0) && roi.at(1) <= posVec.at(1) && roi.at(2) >= posVec.at(2) && roi.at(3) <= posVec.at(3) )
		{
			std::cout << "ROI X max is within module " << modNum << endl;
			roi_modules.insert(modNum);
		}


	}

	for ( const auto& modnum : roi_modules ) std::cout << "ROI is in mod: " << modnum << endl;

	
}