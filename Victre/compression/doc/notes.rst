Notes
=====

- Both breastCompress and FeBiO are multi-threaded applications which respect the environment variable $OMP_NUM_THREADS
- Currently for the purposes of compression modeling all tissue is treated as either fat, glandular, or muscle
- Several parameters that users may want to adjust cannot be specified as a command line option and must be changed in the source code file breastCompress.cxx.  Some of these parameters include

  - Elasticity parameters

    - *fatModulus* - Young's Modulus for fat tissue, default value = 5 kPa
    - *fatPoisson* - Poisson's Ratio for fat tissue, default value = 0.49
    - *fatDensity* - density of fat tissue, default value = 0.000001 kg/mm^3
    - *glandModulus* - Young's Modulus for glandular tissue, default value = 15.0 kPa
    - *glandPoisson* - Poisson's Ratio for glandular tissue, default value = 0.49
    - *glandDensity* - density of glandular tissue, default value = 0.000001 kg/mm^3

  - Paddle geometry

    - *paddleHeight* - thickness of paddle, default value = 1 mm

  - Tetgen parameters - the tetrahedral mesh is generated using the Tetgen library.  Mesh parameters are specified in a tetgenbehavior struct.  See the TetGet documentatio for more details http://wias-berlin.de/software/tetgen/1.5/doc/manual/index.html

  - FeBiO parameters - all finite element parameters are written to the FeBiO XML input file automatically by the code.  To change these parameters consult the FeBiO documentation available at https://febio.org/support/ and the relevant sections of breastCompress.cxx.


