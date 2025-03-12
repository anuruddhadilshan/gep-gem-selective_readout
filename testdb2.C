#include "DBread.h"


void testdb2()
{
	DBread db{"db_FT_local.dat"};

	std::map < int, GEMLayer > gemLayerMap = db.returnGEMLayerMap();

	int nLayers = gemLayerMap.size();

	std::cout << "Num of GEM Layers: " << nLayers << endl;

	std::cout << "GEM Layer 2 mod pos z: " << (((gemLayerMap[2] ).gemModulesLayer)[2]).zPos << endl;

	std::cout << "GEM Layer 7 mod 12 pos x: " << (((gemLayerMap[7] ).gemModulesLayer)[12]).xPos << endl;
	std::cout << "GEM Layer 7 mod 13 pos z: " << (((gemLayerMap[7] ).gemModulesLayer)[13]).zPos << endl;

	

}