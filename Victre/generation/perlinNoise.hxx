/*
 * perlinNoise.hxx
 *
 *  Created on: Apr 24, 2015
 *      Author: cgg
 */

#ifndef PERLINNOISE_HXX_
#define PERLINNOISE_HXX_

#ifndef __BOOST__
	#define __BOOST__
	#include <boost/random.hpp>
	#include <boost/math/distributions/beta.hpp>
	#include <boost/program_options.hpp>
#endif

class perlinNoise{
	 
	private:
		double frequency;
		double lacunarity;
		double persistence;
		int numOctaves;
		int32_t seed;
		int32_t xNoiseGen,yNoiseGen,zNoiseGen,seedNoiseGen,shiftNoiseGen;
		double makeInt32Range(double x);
		double pInterp(double x);
		double linInterp(double lbound, double rbound, double x);
		double coherentNoise(double x, double y, double z, int32_t mySeed);
		double gradientNoise(double x, double y, double z, 
			int32_t ia, int32_t ib, int32_t ic, int32_t mySeed);
	
	public:
		double getNoise(double* r);
		void setSeed(int32_t inSeed);
		perlinNoise(boost::program_options::variables_map vm, int32_t inSeed, const char* type);
		perlinNoise(boost::program_options::variables_map vm, int32_t inSeed, double freq, double lac, double pers, int oct);
		perlinNoise(boost::program_options::variables_map vm, const char* type);
		
};

#endif /* PERLINNOISE_HXX_ */
		
	
		
