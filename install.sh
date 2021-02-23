



CUDA_INC=/projects01/sysadmin/stuartb/cuda/cuda-10.2/cuda-samples/common/inc/





echo "
============================================
    __   __ ___  ___  _____  ___  ___       
    \ \ / /|_ _|/ __||_   _|| _ \| __|      
     \ V /  | || (__   | |  |   /| _|       
      \_/  |___|\___|  |_|  |_|_\|___|      
  ___  ___  ___  ___  _     ___  _  _  ___ 
 | _ \|_ _|| _ \| __|| |   |_ _|| \| || __|
 |  _/ | | |  _/| _| | |__  | | | .\` || _| 
 |_|  |___||_|  |___||____||___||_|\_||___|
============================================

WELCOME TO THE VICTRE PIPELINE INSTALLATION SCRIPT.

We are going to check the requirements.
"

check_command () {
    if ! command -v $1 &> /dev/null
    then
        exists=0
        return
    fi
    exists=1
}

RED='\033[1;31m'
GREEN='\033[1;32m'
NC='\033[0m' # No Color

printf "Checking nvcc: \t\t\t"
check_command nvcc
NVCC=$exists
if [ $exists = 1 ]
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You won't be able to compile MCGPU on this computer.\n"
fi

printf "Checking gcc: \t\t\t"
check_command gcc
GCC=$exists
if [ $exists = 1 ]
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You won't be able to compile many things on this computer."
fi

printf "Checking cmake: \t\t"
check_command cmake
CMAKE=$exists
if [ $exists = 1 ]
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You won't be able to compile breast generation, compression and mass generation on this computer."
fi

printf "Checking python: \t\t"
check_command python
PYTHON=$exists
if [ $exists = 1 ]
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You won't be able to run the pipeline in this computer."
fi

printf "Checking febio2.lnx64: \t\t"
check_command febio2.lnx64
FEBIO=$exists
if [ $exists = 1 ]
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You won't be able to run the breast compression in this computer."
fi

printf "Checking VTK_DIR \t\t"
if [ ! -z "$VTK_DIR" ];
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You might not be able to compile breast generation, compression and mass generation on this computer."
fi

printf "Checking boost \t\t\t"
if [[ "$C_INCLUDE_PATH" == *"boost"* && "$CPLUS_INCLUDE_PATH" == *"boost"* ]]; then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You might not be able to compile breast generation, compression and mass generation on this computer."
fi


echo "
Before you proceed, you need to compile some parts
of the pipeline."

while [ true ] ; do
    echo "
What do you want to compile?
    1. Breast generation
    2. Breast compression
    3. Breast mass
    4. MCGPU projection
    5. Reconstruction software

Press any other key to exit...
    "

    read -n 1 -t 15 a
    printf "\n"
    case $a in
    1* )
        echo "Compiling breast generation"
        rm -rf ./Victre/generation/build 
        mkdir ./Victre/generation/build
        (cd ./Victre/generation/build && cmake ..)
        (cd ./Victre/generation/build && make)
    ;;
    
    2* )     
        echo "Compiling breast compression"
        rm -rf ./Victre/compression/build 
        mkdir ./Victre/compression/build
        (cd ./Victre/compression/build && cmake ..)
        (cd ./Victre/compression/build && make)
    ;;

    3* )     
        echo "Compiling breast mass"
        rm -rf ./Victre/breastMass/build 
        mkdir ./Victre/breastMass/build
        (cd ./Victre/breastMass/build && cmake ..)
        (cd ./Victre/breastMass/build && make)
    ;;

    4* )
        ( cd ./Victre/projection && nvcc MC-GPU_v1.5b.cu -o MC-GPU_v1.5b.x -m64 -O3 -use_fast_math -DUSING_MPI -I. -I$CUDA_INC -I $MPI_INCLUDE -L $MPI_LIB -lmpi -lz --ptxas-options=-v -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70)
    ;;

    5* )     
        (cd ./Victre/reconstruction && make)
        (cd ./Victre/reconstruction && gcc extract_projections.c -O -o extract_projections_RAW.x)
    ;;

    * )     break;;
    esac
done

# read -p "Press enter to continue"
