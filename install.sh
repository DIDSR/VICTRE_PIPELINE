##########################################################################
# UNCOMMENT THESE VARIABLES TO INDICATE CUSTOM PATHS FOR THESE LIBRARIES #
##########################################################################

# CUDA_INC=/usr/local/cuda/include 
# CUDA_SAMPLES=/usr/local/cuda/samples/common/inc
# MPI_INCLUDE=/usr/include/openmpi
# BOOST_OPTIONS_LIB=/usr/lib64 # route to libboost_program_options
# VTK_DIR=/usr/local/VTK-build


# ------------------ DO NOT MODIFY --------------------------

RED='\033[1;31m'
GREEN='\033[1;32m'
BLUE='\033[1;36m'
NC='\033[0m' # No Color

echo -e "$GREEN
==================================================
        __   __ ___  ___  _____  ___  ___       
        \ \ / /|_ _|/ __||_   _|| _ \| __|      
         \ V /  | || (__   | |  |   /| _|       
          \_/  |___|\___|  |_|  |_|_\|___|      
     ___  ___  ___  ___  _     ___  _  _  ___ 
    | _ \|_ _|| _ \| __|| |   |_ _|| \| || __|
    |  _/ | | |  _/| _| | |__  | | | .\` || _| 
    |_|  |___||_|  |___||____||___||_|\_||___|
    
==================================================
$NC
${BLUE}Citation:${NC} \"Evaluation of Digital Breast Tomosynthesis as Replacement of Full-Field Digital Mammography Using an In Silico Imaging Trial.\" Aldo Badano, Ph. D., Christian G. Graff, Ph. D., Andreu Badal, Ph. D., Diksha Sharma, M. Sc., Rongping Zeng, Ph. D., Frank W. Samuelson, Ph. D., Stephen Glick, Ph. D., and Kyle J. Myers, Ph. D. JAMA Network Open. 2018;1(7):e185474; doi:10.1001/jamanetworkopen.2018.5474.

${BLUE}VICTRE team:${NC} Aldo Badano, Ph. D., Christian G. Graff, Ph. D., Andreu Badal, Ph. D., Diksha Sharma, M. Sc., Rongping Zeng, Ph. D., Aunnasha Sengupta, Miguel A. Lago, Ph. D., Frank W. Samuelson, Ph. D., Stephen Glick, Ph. D., and Kyle J. Myers, Ph. D.
"
read -n 1 -r -s -p $'Press enter to continue...\n'
check_command () {
    if ! command -v $1 &> /dev/null
    then
        exists=0
        return
    fi
    exists=1
}


export C_INCLUDE_PATH=$BOOST_INCLUDE:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${C_INCLUDE_PATH}:${CPLUS_INCLUDE_PATH}

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
    printf "\t*** You won't be able to compile many things on this computer.\n"
fi

printf "Checking cmake: \t\t"
check_command cmake
CMAKE=$exists
if [ $exists = 1 ]
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You won't be able to compile breast generation, compression and mass generation on this computer.\n"
fi

printf "Checking python: \t\t"
check_command python
PYTHON=$exists
if [ $exists = 1 ]
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You won't be able to run the pipeline in this computer.\n"
fi

printf "Checking febio2: \t\t"
check_command febio2
FEBIO=$exists
if [ $exists = 1 ]
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You won't be able to run the breast compression in this computer.\n"
fi

printf "Checking VTK_DIR \t\t"
if [ ! -z "$VTK_DIR" ];
then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You might not be able to compile breast generation, compression and mass generation on this computer.\n"
fi

printf "Checking boost \t\t\t"
if [[ "$C_INCLUDE_PATH" == *"boost"* && "$CPLUS_INCLUDE_PATH" == *"boost"* ]]; then
    printf "${GREEN}FOUND${NC}\n"
else
    printf "${RED}NOT FOUND${NC}"
    printf "\t*** You might not be able to compile breast generation, compression and mass generation on this computer.\n"
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
        printf "
What is your GPU architecture?
    1. Maxwell: GTX 745-980, Quadro Kxxx
    2. Pascal: GTX 1070, Tesla P100
    3. Volta: Tesla V100
    4. Turing: GTX 16X0, RTX 20X0, Quadro RTX, Tesla T4
    5. Ampere: RTX 30X0, RTX A6000, Tesla A100
"
        read -n 1 -t 15 a
        case $a in
            1* ) GPU_OPTIONS="-gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_53,code=sm_53" ;;
            2* ) GPU_OPTIONS="-gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_62,code=sm_62" ;;
            3* ) GPU_OPTIONS="-gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_72,code=sm_72" ;;
            4* ) GPU_OPTIONS="-gencode=arch=compute_75,code=sm_75" ;;
            5* ) GPU_OPTIONS="-gencode=arch=compute_80,code=sm_80" ;;
        esac
        printf "\n\n"
        read -n 1 -p "Compile with MPI support? (y/n): " yn
        case $yn in
            [Yy]* ) ( cd ./Victre/projection && nvcc MC-GPU_v1.5b.cu -o MC-GPU_v1.5b.x -m64 -O2 -DUSING_MPI -I. -I$CUDA_INC -I$CUDA_SAMPLES -I $MPI_INCLUDE -lmpi -lz --ptxas-options=-v $GPU_OPTIONS );;
            [Nn]* ) ( cd ./Victre/projection && nvcc MC-GPU_v1.5b.cu -o MC-GPU_v1.5b.x -m64 -O2 -I. -I$CUDA_INC -I$CUDA_SAMPLES -lz --ptxas-options=-v $GPU_OPTIONS );;
            * ) echo "No answer";;
        esac
    ;;

    5* )     
        (cd ./Victre/reconstruction && make)
        (cd ./Victre/reconstruction && gcc extract_projections.c -O -o extract_projections_RAW.x)
    ;;

    * )     break;;
    esac
done

# read -p "Press enter to continue"
