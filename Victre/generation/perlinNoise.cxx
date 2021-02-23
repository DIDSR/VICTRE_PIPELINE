/*
 * perlinNoise.cxx
 *
 *  Created on: Apr 24, 2015
 *      Author: cgg
 */
 
#include "perlinNoise.hxx"

double randPerm[256][3] = {
	{0.7321, 0.3079, 0.6076},
	{-0.7627, -0.6454, 0.0418},
	{-0.6154, 0.7791, 0.1197},
	{0.1981, -0.0670, -0.9779},
	{-0.4692, -0.8294, -0.3032},
	{0.9552, -0.0901, 0.2820},
	{-0.7339, 0.6183, -0.2814},
	{-0.5067, 0.1114, -0.8549},
	{0.9405, -0.2913, -0.1750},
	{-0.4266, -0.5980, -0.6785},
	{-0.6682, 0.0107, 0.7439},
	{0.4419, -0.4260, -0.7895},
	{-0.1338, 0.0352, 0.9904},
	{-0.2090, 0.4668, 0.8593},
	{-0.4557, 0.8243, -0.3359},
	{-0.9385, 0.0787, 0.3363},
	{0.1633, 0.1045, 0.9810},
	{-0.2892, 0.7112, 0.6407},
	{-0.8552, -0.4496, -0.2578},
	{0.2314, -0.0248, 0.9726},
	{-0.6555, 0.4092, -0.6347},
	{-0.9242, 0.2578, 0.2819},
	{0.6887, 0.2845, 0.6669},
	{0.5184, 0.7331, -0.4402},
	{-0.0003, 0.0628, -0.9980},
	{-0.1868, -0.5701, 0.8000},
	{0.8902, -0.0258, -0.4549},
	{0.4641, -0.0894, -0.8813},
	{-0.8292, 0.5297, -0.1783},
	{0.8903, -0.4453, -0.0955},
	{-0.6731, 0.7060, -0.2205},
	{-0.2349, 0.7489, -0.6197},
	{0.1263, -0.3793, -0.9166},
	{0.4793, 0.7522, -0.4523},
	{-0.6545, 0.4021, -0.6403},
	{-0.6249, -0.6880, 0.3691},
	{-0.6274, -0.7774, 0.0445},
	{-0.5314, -0.4337, 0.7277},
	{-0.4342, 0.0934, -0.8960},
	{-0.8424, -0.1205, 0.5252},
	{-0.6197, 0.6002, -0.5057},
	{0.6002, 0.4990, 0.6251},
	{-0.1960, 0.3326, 0.9224},
	{-0.4486, 0.7984, -0.4017},
	{-0.3491, -0.6902, -0.6338},
	{0.3503, -0.6902, 0.6332},
	{0.0054, 0.6057, 0.7957},
	{-0.3844, -0.6576, 0.6479},
	{0.2750, 0.7956, -0.5398},
	{0.0821, 0.9815, 0.1727},
	{0.1032, 0.6621, 0.7422},
	{-0.9786, -0.0871, -0.1866},
	{-0.1505, -0.8425, -0.5173},
	{-0.9588, 0.0868, 0.2704},
	{-0.1923, 0.5672, 0.8008},
	{-0.7584, -0.1827, -0.6256},
	{-0.5757, 0.6563, 0.4878},
	{0.5126, 0.3276, -0.7937},
	{-0.2219, 0.8150, 0.5353},
	{-0.6929, 0.2327, 0.6825},
	{0.8388, -0.3683, -0.4009},
	{-0.8311, 0.5492, -0.0871},
	{0.4851, 0.1866, -0.8543},
	{0.6251, 0.3983, -0.6712},
	{-0.4376, 0.8912, -0.1198},
	{0.0918, 0.9372, 0.3364},
	{0.1645, -0.5193, 0.8386},
	{0.6724, -0.6643, 0.3264},
	{-0.5813, 0.4588, 0.6720},
	{0.5544, 0.0660, 0.8296},
	{-0.0322, 0.9961, -0.0825},
	{0.5919, -0.2674, -0.7604},
	{-0.7931, -0.0770, 0.6043},
	{-0.2930, -0.5563, -0.7776},
	{-0.5586, -0.4115, -0.7201},
	{-0.6794, -0.6753, 0.2871},
	{0.7238, -0.2761, -0.6324},
	{0.9019, 0.3714, 0.2206},
	{0.0955, -0.9812, -0.1678},
	{-0.5817, -0.1072, -0.8063},
	{-0.2312, 0.9726, -0.0263},
	{0.5445, 0.7721, 0.3276},
	{0.7307, -0.4999, 0.4649},
	{0.8150, -0.5755, -0.0673},
	{-0.1454, -0.0326, -0.9888},
	{0.0044, 0.8159, 0.5781},
	{0.6092, 0.5649, -0.5566},
	{0.8413, 0.5303, -0.1046},
	{-0.0908, 0.1134, 0.9894},
	{-0.9169, 0.1078, -0.3843},
	{-0.6039, 0.2603, -0.7534},
	{-0.6366, 0.6036, 0.4799},
	{-0.2517, -0.9380, -0.2383},
	{-0.7743, 0.3626, 0.5186},
	{0.6240, -0.6339, -0.4569},
	{-0.5114, -0.4063, -0.7572},
	{0.0686, -0.9293, -0.3630},
	{-0.6761, -0.3906, -0.6248},
	{0.2617, 0.9569, 0.1258},
	{0.3088, -0.8581, 0.4102},
	{-0.7167, -0.6472, -0.2597},
	{0.2905, 0.5797, -0.7613},
	{-0.1163, 0.9880, 0.1013},
	{0.1774, -0.8937, -0.4121},
	{0.8012, 0.1367, 0.5826},
	{0.7638, -0.2467, -0.5965},
	{0.8406, 0.2085, -0.5000},
	{-0.4977, 0.0109, 0.8673},
	{0.1136, -0.3604, 0.9259},
	{-0.0440, -0.6937, -0.7189},
	{0.0848, -0.6461, 0.7586},
	{0.4835, 0.7675, -0.4210},
	{-0.5589, -0.7499, -0.3539},
	{0.3558, 0.1458, -0.9231},
	{-0.4700, -0.7746, -0.4232},
	{-0.7533, -0.6374, 0.1620},
	{0.9568, 0.2172, 0.1931},
	{0.2123, 0.6267, 0.7498},
	{-0.2250, 0.7109, 0.6663},
	{0.3891, 0.6807, 0.6207},
	{-0.6638, -0.6985, 0.2675},
	{-0.6995, -0.0115, -0.7146},
	{-0.5760, -0.5446, -0.6097},
	{0.0210, 0.9670, 0.2541},
	{0.4853, -0.7614, 0.4297},
	{-0.3825, 0.2933, 0.8761},
	{-0.1845, -0.5498, -0.8146},
	{-0.0640, -0.5714, -0.8182},
	{-0.9255, 0.1666, 0.3401},
	{0.5173, 0.7951, -0.3165},
	{-0.7908, -0.3113, -0.5270},
	{0.7585, 0.0755, -0.6473},
	{-0.6738, -0.0001, -0.7389},
	{-0.2761, -0.5959, 0.7541},
	{0.8625, -0.4371, -0.2550},
	{0.9803, -0.1719, 0.0972},
	{0.1651, 0.8711, 0.4625},
	{-0.7632, 0.2654, 0.5892},
	{-0.9250, 0.3507, 0.1462},
	{0.3641, 0.7392, -0.5666},
	{0.8223, -0.4537, 0.3434},
	{0.9014, -0.1115, 0.4184},
	{0.9918, 0.1271, -0.0110},
	{0.1672, 0.8429, -0.5115},
	{-0.7638, 0.6284, -0.1471},
	{0.0135, 0.6228, -0.7822},
	{0.1823, 0.9755, 0.1233},
	{0.5311, 0.4948, 0.6878},
	{-0.4295, 0.8889, 0.1595},
	{0.2574, -0.9654, -0.0428},
	{-0.3177, -0.8373, -0.4450},
	{-0.8000, 0.5946, 0.0808},
	{0.4582, -0.0408, -0.8879},
	{-0.3659, 0.2641, 0.8924},
	{-0.8128, 0.5693, -0.1239},
	{-0.3630, -0.9288, -0.0744},
	{-0.1302, 0.5226, 0.8426},
	{-0.9377, 0.2284, -0.2617},
	{-0.8814, 0.2604, 0.3941},
	{-0.5333, 0.2432, -0.8102},
	{0.2519, 0.9312, -0.2635},
	{0.5826, 0.4235, 0.6937},
	{-0.5009, -0.2137, -0.8387},
	{0.0946, 0.2682, 0.9587},
	{-0.0795, -0.8703, 0.4861},
	{0.3588, -0.7481, 0.5583},
	{-0.7646, 0.4739, 0.4368},
	{-0.3202, -0.7916, -0.5204},
	{0.3130, -0.8557, -0.4120},
	{-0.8588, -0.3781, -0.3458},
	{-0.7770, -0.6251, -0.0738},
	{-0.4739, 0.1896, 0.8599},
	{-0.6780, -0.5672, -0.4676},
	{0.1215, 0.4195, 0.8996},
	{-0.1547, -0.6755, -0.7210},
	{0.5227, 0.1780, -0.8337},
	{-0.9390, 0.3289, 0.1000},
	{0.2292, -0.8200, -0.5245},
	{-0.0877, -0.6893, -0.7192},
	{-0.6912, 0.4608, -0.5566},
	{-0.8712, 0.4851, -0.0753},
	{0.2609, -0.7249, 0.6376},
	{-0.4690, 0.8826, 0.0336},
	{0.9610, -0.0747, -0.2663},
	{-0.3423, -0.8989, -0.2734},
	{-0.4653, -0.1613, -0.8703},
	{0.5476, -0.1433, -0.8244},
	{0.2540, 0.6240, 0.7390},
	{0.4097, -0.8361, 0.3650},
	{0.7497, -0.1641, -0.6411},
	{-0.8612, -0.2650, 0.4337},
	{0.1382, 0.3693, 0.9190},
	{-0.1560, 0.9805, -0.1194},
	{0.5390, 0.8262, 0.1637},
	{-0.8426, -0.0430, 0.5368},
	{-0.9044, 0.4214, -0.0669},
	{-0.3714, 0.9242, 0.0889},
	{-0.9331, 0.1464, -0.3284},
	{0.0086, 0.8716, 0.4902},
	{-0.7630, 0.5044, 0.4042},
	{-0.6998, -0.4462, -0.5578},
	{0.2384, -0.1265, -0.9629},
	{0.2794, 0.2446, 0.9285},
	{-0.6599, -0.7477, 0.0747},
	{0.7670, -0.4827, 0.4226},
	{0.6568, 0.0349, -0.7533},
	{0.3527, -0.3232, 0.8782},
	{0.4858, 0.8656, -0.1216},
	{-0.8790, 0.3887, 0.2762},
	{-0.8942, 0.4476, -0.0018},
	{-0.9945, -0.0543, 0.0894},
	{-0.7792, -0.3553, 0.5163},
	{0.0909, 0.9899, -0.1089},
	{0.6303, -0.6459, 0.4308},
	{-0.9626, -0.2425, -0.1206},
	{-0.4803, 0.2245, 0.8479},
	{-0.9055, -0.4231, 0.0312},
	{0.2403, 0.5935, -0.7681},
	{0.0634, 0.2091, -0.9758},
	{0.4437, -0.1539, -0.8829},
	{0.8685, -0.0428, 0.4938},
	{0.2929, -0.6249, -0.7236},
	{0.5312, -0.7727, 0.3476},
	{0.5005, 0.4707, -0.7266},
	{-0.8272, -0.2612, 0.4975},
	{0.5923, -0.4378, 0.6763},
	{-0.0769, 0.9183, 0.3884},
	{0.1037, -0.4437, 0.8901},
	{0.7844, 0.0410, -0.6189},
	{0.1837, 0.7913, -0.5832},
	{-0.6373, 0.7662, -0.0823},
	{-0.1453, -0.6667, 0.7310},
	{0.8440, -0.3984, 0.3591},
	{0.2133, 0.9769, -0.0115},
	{0.8048, -0.4878, 0.3381},
	{0.5298, 0.1217, 0.8393},
	{-0.2886, 0.2390, 0.9271},
	{-0.5014, 0.7109, 0.4932},
	{-0.6136, -0.6999, -0.3656},
	{-0.7680, -0.3627, -0.5279},
	{-0.3133, 0.3102, 0.8976},
	{-0.6481, -0.1708, -0.7422},
	{0.0870, 0.9850, 0.1491},
	{0.9637, 0.2520, -0.0885},
	{0.0111, 0.5720, -0.8202},
	{0.0320, -0.9797, -0.1981},
	{0.3254, 0.5055, -0.7991},
	{0.6722, -0.0630, -0.7377},
	{-0.4864, -0.4289, -0.7612},
	{0.3559, 0.8705, -0.3398},
	{0.8758, 0.4350, 0.2093},
	{-0.8788, -0.4454, 0.1713},
	{-0.5031, -0.7045, 0.5006},
	{0.2308, -0.0383, -0.9723},
	{-0.3010, -0.9137, 0.2732},
	{-0.2747, 0.9577, 0.0858}
};

perlinNoise::perlinNoise(boost::program_options::variables_map vm, int32_t inSeed, const char* type){
	 
	 // read values from options
	 char varName[80];
	 strcpy(varName, type);
	 strcat(varName, ".frequency");
	 frequency = vm[varName].as<double>();
	 strcpy(varName, type);
	 strcat(varName, ".lacunarity");
	 lacunarity = vm[varName].as<double>();
	 strcpy(varName, type);
	 strcat(varName, ".persistence");
	 persistence = vm[varName].as<double>();
	 
	 numOctaves = vm["perlin.numOctaves"].as<int>();
	 xNoiseGen = (int32_t)vm["perlin.xNoiseGen"].as<int>();
	 yNoiseGen = (int32_t)vm["perlin.yNoiseGen"].as<int>();
	 zNoiseGen = (int32_t)vm["perlin.zNoiseGen"].as<int>();
	 seedNoiseGen = (int32_t)vm["perlin.seedNoiseGen"].as<int>();
	 shiftNoiseGen = (int32_t)vm["perlin.shiftNoiseGen"].as<int>();
	 seed = inSeed;
};

perlinNoise::perlinNoise(boost::program_options::variables_map vm, const char* type){
	 
	 // read values from options
	 char varName[80];
	 strcpy(varName, type);
	 strcat(varName, ".frequency");
	 frequency = vm[varName].as<double>();
	 strcpy(varName, type);
	 strcat(varName, ".lacunarity");
	 lacunarity = vm[varName].as<double>();
	 strcpy(varName, type);
	 strcat(varName, ".persistence");
	 persistence = vm[varName].as<double>();
	 
	 numOctaves = vm["perlin.numOctaves"].as<int>();
	 xNoiseGen = (int32_t)vm["perlin.xNoiseGen"].as<int>();
	 yNoiseGen = (int32_t)vm["perlin.yNoiseGen"].as<int>();
	 zNoiseGen = (int32_t)vm["perlin.zNoiseGen"].as<int>();
	 seedNoiseGen = (int32_t)vm["perlin.seedNoiseGen"].as<int>();
	 shiftNoiseGen = (int32_t)vm["perlin.shiftNoiseGen"].as<int>();
	 seed = 334;	// default seed
};

perlinNoise::perlinNoise(boost::program_options::variables_map vm, int32_t inSeed, double freq, double lac, double pers, int oct){
	frequency = freq;
	lacunarity = lac;
	persistence = pers;
	numOctaves = oct;
	seed = inSeed;
	xNoiseGen = (int32_t)vm["perlin.xNoiseGen"].as<int>();
	yNoiseGen = (int32_t)vm["perlin.yNoiseGen"].as<int>();
	zNoiseGen = (int32_t)vm["perlin.zNoiseGen"].as<int>();
	seedNoiseGen = (int32_t)vm["perlin.seedNoiseGen"].as<int>();
	shiftNoiseGen = (int32_t)vm["perlin.shiftNoiseGen"].as<int>();
}
	
inline double perlinNoise::makeInt32Range(double x){
	if(x >= 1073741824.0){
		return (2.0*fmod(x, 1073741824.0)) - 1073741824.0;
	} else if (x <= -1073741824.0){
		return (2.0*fmod(x, 1073741824.0)) + 1073741824.0;
	} else {
		return x;
	}
}

inline double perlinNoise::pInterp(double x){
	double x3 = x*x*x;
	return (6.0*x3*x*x - 15.0*x3*x + 10.0*x3);
}

inline double perlinNoise::linInterp(double lbound, double rbound, double x){
	return ((1.0 - x)*lbound + x*rbound);
}

double perlinNoise::gradientNoise(double x, double y, double z, 
	int32_t ia, int32_t ib, int32_t ic, int32_t mySeed){
		
	int32_t vectorInd = (xNoiseGen*ia + yNoiseGen*ib + zNoiseGen*ic
		+ seedNoiseGen*mySeed) & 0xffffffff;
		
	vectorInd ^= (vectorInd >> shiftNoiseGen);
	vectorInd &= 0xff;
	
	double xGrad = randPerm[vectorInd][0];
	double yGrad = randPerm[vectorInd][1];
	double zGrad = randPerm[vectorInd][2];
	
	double dx = x - (double)ia;
	double dy = y - (double)ib;
	double dz = z - (double)ic;
	
	return (xGrad*dx + yGrad*dy + zGrad*dz)*2.12;
};

double perlinNoise::coherentNoise(double x, double y, double z, int32_t mySeed){
	
	// make integer coordinate surrounding cube
	int32_t xlb = (x > 0.0 ? (int32_t)x : (int32_t)x - 1);
	int32_t xub = xlb + 1;
	int32_t ylb = (y > 0.0 ? (int32_t)y : (int32_t)y - 1);
	int32_t yub = ylb + 1;
	int32_t zlb = (z > 0.0 ? (int32_t)z : (int32_t)z - 1);
	int32_t zub = zlb + 1;
	
	// interpolation based on location in cube
	double xInterp = pInterp(x - (double)xlb);
	double yInterp = pInterp(y - (double)ylb);
	double zInterp = pInterp(z - (double)zlb);
	
	// calculate gradient noise at cube points and interpolate to target location
	double na,nb;
	double ix0,ix1,iy0,iy1;
	
	na = gradientNoise(x, y, z, xlb, ylb, zlb, mySeed);
	nb = gradientNoise(x, y, z, xub, ylb, zlb, mySeed);
	ix0 = linInterp(na, nb, xInterp);
	na = gradientNoise(x, y, z, xlb, yub, zlb, mySeed);
	nb = gradientNoise(x, y, z, xub, yub, zlb, mySeed);
	ix1 = linInterp(na, nb, xInterp);
	iy0 = linInterp(ix0, ix1, yInterp);
	na = gradientNoise(x, y, z, xlb, ylb, zub, mySeed);
	nb = gradientNoise(x, y, z, xub, ylb, zub, mySeed);
	ix0 = linInterp(na, nb, xInterp);
	na = gradientNoise(x, y, z, xlb, yub, zub, mySeed);
	nb = gradientNoise(x, y, z, xub, yub, zub, mySeed);
	ix1 = linInterp(na, nb, xInterp);
	iy1 = linInterp(ix0, ix1, yInterp);
	
	return linInterp(iy0, iy1, zInterp);
}


void perlinNoise::setSeed(int32_t inSeed){
	seed = inSeed;
}

double perlinNoise::getNoise(double* r){
	
	double nval = 0.0;
	double signal = 0.0;
	double myPersistence = 1.0;
	
	int32_t mySeed;
	
	double xv = r[0]*frequency;
	double yv = r[1]*frequency;
	double zv = r[2]*frequency;
	
	for(int32_t i=0; i<numOctaves; i++){
		
		xv = makeInt32Range(xv);
		yv = makeInt32Range(yv);
		zv = makeInt32Range(zv);
		
		mySeed = (seed + i) & 0xffffffff;
		signal = coherentNoise(xv, yv, zv, mySeed);
		nval += signal*myPersistence;
		
		xv *= lacunarity;
		yv *= lacunarity;
		zv *= lacunarity;
		myPersistence *= persistence;
	}
	
	return nval;
}

	
	
	 
