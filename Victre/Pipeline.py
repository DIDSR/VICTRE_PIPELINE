"""!
Documentation for the Victre pipeline class.
"""


import numpy as np
import os
from termcolor import colored, cprint
import shutil
from os.path import isfile, join
from os import walk
import contextlib
import pathlib
import glob
import progressbar
import h5py
import subprocess
import asyncio
from string import Template
import random
import time
from . import Constants
import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
import copy
import datetime
from pydicom.encaps import encapsulate
import re
import gzip
from scipy import interpolate


class Pipeline:
    """!
    Victre Pipeline class. **This is all you need!**
    """

    def __init__(self,
                 ips={"cpu": "localhost", "gpu": "localhost"},
                 seed=None,
                 results_folder="./results",
                 phantom_file=None,
                 spectrum_file="./Victre/projection/spectrum/W28kVp_Rh50um_Be1mm.spc",
                 lesion_file=None,
                 materials=None,
                 roi_sizes=None,
                 arguments_generation=dict(),
                 arguments_mcgpu=dict(),
                 arguments_recon=dict(),
                 flatfield_DBT=None,
                 flatfield_DM=None,
                 density=None):
        """!
        Object constructor for the Victre pipeline class

        @param ips Dictionary with two IP addresses to run the pipeline: "gpu" for the projection process. "cpu" for the reconstruction.
        @param seed Random seed used to generate or read the phantom
        @param results_folder Path to folder to be used when saving the results
        @param phantom_file Path to file containing the phantom to be loaded
        @param spectrum_file Path to file containing the spectrum used to project in MCGPU
        @param lesion_file Path to file containing the lesion to be inserted (in HDF5 format)
        @param materials Dictionary including the materials to be used during projection
        @param roi_sizes Dictionary with the ROI sizes for the extraction
        @param arguments_generation Arguments to be overriden for the breast phantom generation
        @param arguments_mcgpu Arguments to be overridden for the projection in MCGPU
        @param arguments_recon Arguments to be overridden for the reconstruction algorithm
        @param flatfield_DBT Path to the flatfield file for the DBT reconstruction
        @param flatfield_DM Path to the flatfield file for the digital mammography
        @param density [EXPERIMENTAL] Percentage of dense tissue of the phantom to be generated, this will adjust the compression thickness too
        @return None
        """

        if seed is None:
            self.seed = int(time.time())
        else:
            self.seed = seed

        self.ips = ips
        self.lesion_file = lesion_file
        self.lesions = []
        self.lesion_locations = {"dbt": [], "dm": []}
        self.results_folder = results_folder
        self.roi_sizes = roi_sizes
        self.candidate_locations = None

        random.seed(self.seed)

        self.arguments_mcgpu = Constants.VICTRE_DEFAULT_MCGPU
        self.arguments_mcgpu["spectrum_file"] = spectrum_file
        self.arguments_mcgpu["phantom_file"] = phantom_file
        self.arguments_mcgpu["output_file"] = "{:s}/{:d}/projection".format(
            self.results_folder, self.seed)

        locations = None

        if phantom_file is None:
            if os.path.exists("{:s}/{:d}/pcl_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found phantom with lesions information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pcl_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                locations = np.loadtxt(
                    "{:s}/{:d}/pcl_{:d}.loc".format(self.results_folder, self.seed, self.seed)).tolist()

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pcl_{:d}.raw.gz".format(
                    self.results_folder, seed, seed)

            elif os.path.exists("{:s}/{:d}/pc_{:d}_crop.mhd".format(self.results_folder, seed, seed)):
                cprint("Found cropped phantom information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pc_{:d}_crop.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                self.candidate_locations = np.loadtxt(
                    "{:s}/{:d}/pc_{:d}_crop.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}_crop.raw.gz".format(
                    self.results_folder, seed, seed)
            elif os.path.exists("{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found compressed phantom information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                self.candidate_locations = np.loadtxt(
                    "{:s}/{:d}/pc_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}.raw.gz".format(
                    self.results_folder, seed, seed)

            elif os.path.exists("{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found phantom generation information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                self.candidate_locations = np.loadtxt(
                    "{:s}/{:d}/p_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/p_{:d}.raw.gz".format(
                    self.results_folder, seed, seed)

        self.arguments_mcgpu.update(arguments_mcgpu)

        self.materials = materials
        if self.materials is None:
            self.materials = Constants.VICTRE_DEFAULT_MATERIALS

        self.arguments_recon = dict(
            number_projections=self.arguments_mcgpu["number_projections"],
            detector_elements=self.arguments_mcgpu["image_pixels"][0],
            detector_elements_perpendicular=self.arguments_mcgpu["image_pixels"][1],
            pixel_size=self.arguments_mcgpu["image_size"][0] /
            self.arguments_mcgpu["image_pixels"][0],
            distance_source=self.arguments_mcgpu["distance_source"],
            rotation_axis_distance=self.arguments_mcgpu["rotation_axis_distance"],
            detector_offset=0.000,
            orbit_projection=50.0,
            voxels_x=self.arguments_mcgpu["number_voxels"][1],
            voxels_y=self.arguments_mcgpu["number_voxels"][0],
            voxels_z=self.arguments_mcgpu["number_voxels"][2],
            voxel_size=self.arguments_mcgpu["voxel_size"][0],
            recon_pixel_size=self.arguments_mcgpu["image_size"][0] /
            self.arguments_mcgpu["image_pixels"][0],
            recon_thickness=0.1,
            volume_center_offset_x=0,
            angular_rotation_first=self.arguments_mcgpu["angular_rotation_first"],
            projections_angle=self.arguments_mcgpu["projections_angle"],
            flatfield_file=flatfield_DBT,
            projection_file="{:s}/{:d}/projection_{:s}pixels_{:d}proj.raw".format(
                self.results_folder,
                self.seed,
                'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                self.arguments_mcgpu["number_projections"]),
            one=1,
            reconstruction_file="{:s}/{:d}/reconstruction{:d}.raw".format(
                self.results_folder,
                self.seed,
                self.seed)
        )

        self.flatfield_DBT = flatfield_DBT
        self.flatfield_DM = flatfield_DM

        if self.flatfield_DBT is None and os.path.exists("{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                self.results_folder,
                self.seed,
                'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                self.arguments_mcgpu["number_projections"])):
            self.arguments_recon["flatfield_file"] = "{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                self.results_folder,
                self.seed,
                'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                self.arguments_mcgpu["number_projections"])
            self.flatfield_DBT = "{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                self.results_folder,
                self.seed,
                'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                self.arguments_mcgpu["number_projections"])

        if self.flatfield_DM is None and os.path.exists("{:s}/{:d}/flatfield_DM{:d}.raw".format(
                self.results_folder,
                self.seed, self.seed)):
            self.flatfield_DM = "{:s}/{:d}/flatfield_DM{:d}.raw".format(
                self.results_folder,
                self.seed, self.seed)

        self.arguments_recon.update(arguments_recon)

        self.arguments_generation = Constants.VICTRE_DENSE  # dense by default

        self.arguments_generation["outputDir"] = os.path.abspath("{:s}/{:d}/".format(
            self.results_folder, self.seed))
        self.arguments_generation["seed"] = self.seed

        if density is not None:
            fat = np.max([0.4, np.min([0.95, 1 - density])])
            ranges = {}
            for key in Constants.DENSITY_RANGES.keys():
                interp = interpolate.interp1d(
                    Constants.DENSITY_RANGES["targetFatFrac"], Constants.DENSITY_RANGES[key])
                ranges[key] = np.round(float(interp(fat)), 2)

            ranges["numBackSeeds"] = int(
                ranges["numBackSeeds"])  # this should be integer

            self.arguments_generation.update(ranges)
            if fat >= 0.75:  # increase the kVp when breast has low density
                # this is hardcoded here, careful
                self.arguments_mcgpu["spectrum_file"] = "./Victre/projection/spectrum/W30kVp_Rh50um_Be1mm.spc"
                self.arguments_mcgpu["fam_beam_aperture"][1] = 11.2

            self.arguments_mcgpu["number_histories"] = ranges["number_histories"]

        self.arguments_generation.update(arguments_generation)

        self.recon_size = dict(
            x=np.ceil(self.arguments_recon["voxels_x"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            y=np.ceil(self.arguments_recon["voxels_y"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            z=np.ceil(self.arguments_recon["voxels_z"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_thickness"]).astype(int)
        )

        if phantom_file is not None:
            splitted = phantom_file.split('/')
            path = '/'.join(splitted[:-1])
            filename = splitted[-1].split('.')[0]

            if os.path.exists("{:s}/{:s}.mhd".format(path, filename)):
                cprint("Found phantom information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:s}.mhd".format(path, filename))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]
            if os.path.exists("{:s}/{:s}.loc".format(path, filename)):
                locations = list(np.loadtxt(
                    "{:s}/{:s}.loc".format(path, filename)))

        if locations is not None:
            self.insert_lesions(locations=locations,
                                save_phantom=False)

        self.arguments_mcgpu["source_position"][1] = self.arguments_mcgpu["number_voxels"][1] * \
            self.arguments_mcgpu["voxel_size"][1] / 2

        # self.arguments_mcgpu["number_voxels"]

        os.makedirs("{:s}".format(self.results_folder), exist_ok=True)
        os.makedirs("{:s}/{:d}".format(self.results_folder,
                                       self.seed), exist_ok=True)

    def project(self, clean=True, do_flatfield=0):
        """!
            Method that runs MCGPU to project the phantom.

            @param clean If True, it will delete the contents of the output folder before projecting.
            @param do_flatfield If > 0, it will generate an empty flat field projection
        """

        def get_gpu_memory():
            def _output_to_list(x): return x.decode('ascii').split('\n')[:-1]

            ACCEPTABLE_AVAILABLE_MEMORY = 1024
            COMMAND = "nvidia-smi --query-gpu=memory.free --format=csv"
            memory_free_info = _output_to_list(
                subprocess.check_output(COMMAND.split()))[1:]
            memory_free_values = [int(x.split()[0])
                                  for i, x in enumerate(memory_free_info)]

            return memory_free_values

        if do_flatfield > 0:
            filename = "flatfield"
            empty_phantom = np.zeros(
                self.arguments_mcgpu["number_voxels"], np.uint8)
            with gzip.open("{:s}/{:d}/empty_phantom.raw.gz".format(
                    self.results_folder, self.seed), "wb") as gz:
                gz.write(empty_phantom)

            prev_flatfield_DBT, prev_flatfield_DM = None, None

            if os.path.exists("{:s}/{:d}/{:s}_DM{:d}.raw".format(self.results_folder, self.seed, filename, self.seed)):
                prev_flatfield_DM = np.fromfile("{:s}/{:d}/{:}_DM.raw".format(self.results_folder, self.seed, filename),
                                                dtype="float32").reshape(2,
                                                                         self.arguments_recon["detector_elements_perpendicular"],
                                                                         self.arguments_recon["detector_elements"])
            if os.path.exists("{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                self.results_folder,
                self.seed,
                'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                self.arguments_mcgpu["number_projections"])
            ):
                prev_flatfield_DBT = np.fromfile("{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                    self.results_folder,
                    self.seed,
                    'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                    self.arguments_mcgpu["number_projections"]),
                    dtype="float32").reshape(self.arguments_mcgpu["number_projections"],
                                             self.arguments_mcgpu["image_pixels"][0],
                                             self.arguments_mcgpu["image_pixels"][1])
        else:
            filename = "projection"

        # %% PROJECTION
        if clean:
            shutil.rmtree("{:s}/{:d}/{:s}_*".format(self.results_folder,
                                                    self.seed,
                                                    filename), ignore_errors=True)

        # check for MPI-compiled MCGPU
        command = "ldd ./Victre/projection/MC-GPU_v1.5b.x"

        mpi = False
        if "mpi" in str(subprocess.run(command.split(), stdout=subprocess.PIPE).stdout):
            mpi = True

        phantom_config = "{:s}/{:d}/input_{:s}.in".format(
            self.results_folder, self.seed, filename)

        with open("./Victre/projection/configs/template_mcgpu.tpl", "r") as f:
            src = Template(f.read())
            template_arguments = copy.deepcopy(self.arguments_mcgpu)
            if do_flatfield > 0:
                template_arguments["phantom_file"] = "{:s}/{:d}/empty_phantom.raw.gz".format(
                    self.results_folder,
                    self.seed)
                template_arguments["number_histories"] *= Constants.FLATFIELD_DOSE_MULTIPLIER
            template_arguments["output_file"] = "{:s}/{:d}/{:s}".format(
                self.results_folder, self.seed, filename)

            # from MBytes to Bytes and reduce 500MB for extra room
            # this would be for the first GPU
            gpu_ram = (get_gpu_memory()[0] - 500) * 1024 * 1024

            if gpu_ram < template_arguments["number_voxels"][0] * template_arguments["number_voxels"][1] * template_arguments["number_voxels"][2]:
                template_arguments["low_resolution_voxel_size"] = [1, 1, 0]

            for key in template_arguments.keys():
                if type(template_arguments[key]) is list:
                    template_arguments[key] = ' '.join(
                        map(str, template_arguments[key]))
            result = src.substitute(template_arguments)

        materials_write = []
        for mat in self.materials:
            materials_write.append("{:s} density={:f} voxelId={:s}".format(mat["material"],
                                                                           mat["density"],
                                                                           ','.join(map(str, mat["voxel_id"]))))

        with open(phantom_config, "w") as f:
            f.write(result)
            f.writelines(s + '\n' for s in materials_write)

        mpistr = " "
        if mpi:
            mpistr = " mpirun -v -n {:d} ".format(
                self.arguments_mcgpu["number_gpus"])
        command = "cd {:s} && time{:s}./Victre/projection/MC-GPU_v1.5b.x {:s}".format(
            os.getcwd(),
            mpistr,
            phantom_config
        )

        if self.ips["gpu"] == "localhost":
            ssh_command = command
        else:
            ssh_command = "ssh -Y {:s} \"{:s}\"".format(
                self.ips["gpu"], command)
        # print(os.system("cd {:s} && ls -la ./Victre/projection/MC-GPU_v1.5b.x > log.txt".format(os.getcwd())))
        # os.system(ssh_command, sh)
        cprint("Initializing MCGPU for {:s}...".format(filename), 'cyan')

        completed = 0

        process = subprocess.Popen(ssh_command, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)

        bar = None
        with open("{:s}/{:d}/output_{:s}.out".format(self.results_folder, self.seed, filename), "wb") as f:
            while True:
                output = process.stdout.readline().decode("utf-8")
                if output == "" and process.poll() is not None:
                    break
                elif "!!DBT!! Simulating first" in output.strip():
                    cprint(
                        "Starting DM projection, this may take a few minutes...", 'cyan')
                elif "Simulating tomographic projection" in output.strip():
                    if completed == 0:
                        cprint("Starting DBT projection...", 'cyan')
                        bar = progressbar.ProgressBar(
                            max_value=self.arguments_mcgpu["number_projections"])
                        bar.update(0)
                    completed += 1
                    bar.update(completed)
                # rc = process.poll()
                f.write(output.encode('utf-8'))
                f.flush()

        if completed != self.arguments_mcgpu["number_projections"]:
            cprint("\nError while projecting, check the output_{:s}.out file in the results folder".format(filename),
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Projection error")

        bar.finish()
        cprint("Projection finished!", 'green', attrs=['bold'])

        command = "cd {:s} && ./Victre/reconstruction/extract_projections_RAW.x {:s} {:d} 0001 {:s}/{:d}/{:s}".format(
            os.getcwd(),
            ' '.join(map(str, self.arguments_mcgpu["image_pixels"])),
            self.arguments_mcgpu["number_projections"],
            self.results_folder,
            self.seed,
            filename)

        stream = os.popen(command)

        with open("{:s}/{:d}/output_{:s}.out".format(self.results_folder, self.seed, filename), "ab+") as f:
            f.write(stream.read().encode('utf-8'))

        with contextlib.suppress(FileNotFoundError):
            os.remove(
                "{:s}/{:d}/{:s}_0000".format(self.results_folder, self.seed, filename))
            os.remove(
                "{:s}/{:d}/{:s}_DM{:d}.raw".format(self.results_folder, self.seed, filename, self.seed))

        os.rename("{:s}/{:d}/{:s}_0000.raw".format(self.results_folder, self.seed, filename),
                  "{:s}/{:d}/{:s}_DM{:d}.raw".format(self.results_folder, self.seed, filename, self.seed))

        with open("{:s}/{:d}/{:s}_DM{:d}.mhd".format(self.results_folder, self.seed, filename, self.seed), "w") as f:
            src = Template(Constants.MHD_FILE)
            template_arguments = copy.deepcopy(self.mhd)
            template_arguments["ElementSpacing"] = [self.arguments_mcgpu["image_size"][0] / self.arguments_mcgpu["image_pixels"][0] * 10,  # cm to mm
                                                    self.arguments_mcgpu["image_size"][1] / self.arguments_mcgpu["image_pixels"][1] * 10]
            template_arguments["DimSize"] = self.arguments_mcgpu["image_pixels"]
            template_arguments["ElementType"] = "MET_FLOAT"
            template_arguments["NDims"] = 2
            template_arguments["ElementDataFile"] = "{:s}_DM{:d}.raw".format(
                filename, self.seed)
            template_arguments["Offset"] = [0, 0, 0]

            for key in template_arguments.keys():
                if type(template_arguments[key]) is list:
                    template_arguments[key] = ' '.join(
                        map(str, template_arguments[key]))
            result = src.substitute(template_arguments)
            f.write(result)

        for i in range(self.arguments_mcgpu["number_projections"]):
            with contextlib.suppress(FileNotFoundError):
                os.remove(
                    "{:s}/{:d}/{:s}_{:04d}.raw".format(self.results_folder, self.seed, filename, i + 1))
                os.remove(
                    "{:s}/{:d}/{:s}_{:04d}".format(self.results_folder, self.seed, filename, i + 1))

        if do_flatfield > 0:
            os.remove("{:s}/{:d}/empty_phantom.raw.gz".format(
                self.results_folder, self.seed))

            if prev_flatfield_DM is not None:
                curr_flatfield_DM = np.fromfile("{:s}/{:d}/flatfield_DM{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                                dtype="float32").reshape(2,
                                                                         self.arguments_recon["detector_elements_perpendicular"],
                                                                         self.arguments_recon["detector_elements"])

                prev_flatfield_DM += curr_flatfield_DM / \
                    do_flatfield / Constants.FLATFIELD_DOSE_MULTIPLIER

                prev_flatfield_DM.tofile(
                    "{:s}/{:d}/flatfield_DM{:d}.raw".format(self.results_folder, self.seed, self.seed))

            if prev_flatfield_DBT is not None:
                curr_flatfield_DBT = np.fromfile("{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                    self.results_folder,
                    self.seed,
                    'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                    self.arguments_mcgpu["number_projections"]),
                    dtype="float32").reshape(self.arguments_mcgpu["number_projections"],
                                             self.arguments_mcgpu["image_pixels"][0],
                                             self.arguments_mcgpu["image_pixels"][1])

                prev_flatfield_DBT += curr_flatfield_DBT / \
                    do_flatfield / Constants.FLATFIELD_DOSE_MULTIPLIER

                prev_flatfield_DBT.tofile("{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                    self.results_folder,
                    self.seed,
                    'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                    self.arguments_mcgpu["number_projections"]))

        elif self.arguments_recon["flatfield_file"] is None:
            with contextlib.suppress(FileNotFoundError):
                os.remove(
                    "{:s}/{:d}/flatfield_DM{:d}.raw".format(self.results_folder, self.seed, self.seed).format(
                        self.results_folder,
                        self.seed,
                        'x'.join(
                            map(str, self.arguments_mcgpu["image_pixels"])),
                        self.arguments_mcgpu["number_projections"]))
                os.remove(
                    "{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                        self.results_folder,
                        self.seed,
                        'x'.join(
                            map(str, self.arguments_mcgpu["image_pixels"])),
                        self.arguments_mcgpu["number_projections"]))

            # number of iterations to average the flatfield
            for n in range(Constants.FLATFIELD_REPETITIONS):
                cprint("Flatfield file not specified, projecting {:d}/{:d}...".format(
                    n + 1, Constants.FLATFIELD_REPETITIONS), 'cyan')
                self.project(do_flatfield=Constants.FLATFIELD_REPETITIONS)

            self.arguments_recon["flatfield_file"] = "{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                self.results_folder,
                self.seed,
                'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                self.arguments_mcgpu["number_projections"])

            self.flatfield_DM = "{:s}/{:d}/flatfield_DM{:d}.raw".format(
                self.results_folder, self.seed, self.seed)

        if do_flatfield == 0:
            # normalize with flatfield
            curr_flatfield_DM = np.fromfile(self.flatfield_DM,
                                            dtype="float32").reshape(2,
                                                                     self.arguments_recon["detector_elements_perpendicular"],
                                                                     self.arguments_recon["detector_elements"])
            projection_DM = np.fromfile("{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                        dtype="float32").reshape(2,
                                                                 self.arguments_recon["detector_elements_perpendicular"],
                                                                 self.arguments_recon["detector_elements"])

            projection_DM = np.divide(curr_flatfield_DM, projection_DM)

            projection_DM.tofile(
                "{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed))

    def reconstruct(self):
        """!
            Method that runs the reconstruction code for the DBT volume
        """

        # %% RECONSTRUCTION
        with open("./Victre/reconstruction/configs/parameters.tpl", "r") as f:
            src = Template(f.read())
            template_arguments = copy.deepcopy(self.arguments_recon)
            result = src.substitute(template_arguments)

        with open("{:s}/{:d}/input_recon.in".format(self.results_folder, self.seed), "w") as f:
            f.write(result)

        command = "cd {:s} && \
            ./Victre/reconstruction/FBP {:s}/{:d}/input_recon.in".format(
            os.getcwd(),
            self.results_folder,
            self.seed
        )

        if self.ips["cpu"] == "localhost":
            ssh_command = command
        else:
            ssh_command = "ssh -Y {:s} \"{:s}\"".format(
                self.ips["cpu"], command)
        # print(ssh_command)
        # res = os.popen(ssh_command).read()

        self.recon_size = dict(
            x=np.ceil(self.arguments_recon["voxels_x"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            y=np.ceil(self.arguments_recon["voxels_y"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            z=np.ceil(self.arguments_recon["voxels_z"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_thickness"]).astype(int)
        )

        cprint("Initializing reconstruction, this may take a few minutes...", 'cyan')

        completed = 0

        process = subprocess.Popen(ssh_command, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)

        bar = None
        with open("{:s}/{:d}/output_recon.out".format(self.results_folder, self.seed), "wb") as f:
            while True:
                output = process.stdout.readline().decode("utf-8")
                if output == "" and process.poll() is not None:
                    break
                elif "Image slice" in output.strip():
                    if completed == 0:
                        cprint("Starting reconstruction...", 'cyan')
                        bar = progressbar.ProgressBar(
                            max_value=self.recon_size["y"])
                        bar.update(0)
                    completed += 1
                    bar.update(completed)
                    progressbar.streams.flush()
                # rc = process.poll()
                f.write(output.encode('utf-8'))
                f.flush()

        if completed != self.recon_size["y"]:
            cprint("\nError while reconstructing, check the output_recon.out file",
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Reconstruction error")

        bar.finish()

        self.mhd["ElementDataFile"] = "reconstruction{:d}.raw".format(
            self.seed)
        self.mhd["Offset"] = [0, 0, 0]
        self.mhd["DimSize"] = [self.recon_size["x"],
                               self.recon_size["y"],
                               self.recon_size["z"]]
        self.mhd["ElementType"] = "MET_DOUBLE"
        self.mhd["ElementSpacing"] = [self.arguments_recon["recon_pixel_size"] * 10,  # cm to mm
                                      self.arguments_recon["recon_pixel_size"] * 10,
                                      self.arguments_recon["recon_thickness"] * 10]

        with open("{:s}/{:d}/reconstruction{:d}.mhd".format(
                self.results_folder,
                self.seed,
                self.seed), "w") as f:
            src = Template(Constants.MHD_FILE)
            template_arguments = copy.deepcopy(self.mhd)
            for key in template_arguments.keys():
                if type(template_arguments[key]) is list:
                    template_arguments[key] = ' '.join(
                        map(str, template_arguments[key]))
            result = src.substitute(template_arguments)
            f.write(result)

        cprint("Reconstruction finished!", 'green', attrs=['bold'])

    def get_coordinates_dbt(self, vx_location):
        """!
            Method to get the corresponding coordinates in the DBT volume from the voxelized coordinates

            @param vx_location Coordinates in the voxel/phantom space
            @return Coordinates in the DBT space
        """
        location = vx_location.copy()

        location[2] = location[2] - self.arguments_recon["detector_offset"]

        location[0] = location[0] * self.arguments_recon["voxel_size"] / \
            self.arguments_recon["pixel_size"]
        location[1] = location[1] * self.arguments_recon["voxel_size"] / \
            self.arguments_recon["pixel_size"]
        location[2] = location[2] * self.arguments_recon["voxel_size"] / \
            self.arguments_recon["recon_thickness"]  # in mm

        # mirror Y axis
        location[1] = self.recon_size["x"] - location[1]

        # interchange X and Y
        location[0], location[1] = location[1], location[0]
        # location = [location[1], location[0], location[2], location[3]]

        return location

    def get_coordinates_dm(self, vx_location):
        """!
            Method to get the corresponding coordinates in the DM volume from the voxelized coordinates

            @param vx_location Coordinates in the voxel/phantom space
            @return Coordinates in the DM space
        """
        location = vx_location.copy()

        detector_z = - (self.arguments_mcgpu["distance_source"] -
                        self.arguments_mcgpu["source_position"][2])  # -20 / 10

        location[2] = location[2] - self.arguments_recon["detector_offset"]

        location[0] = location[0] * \
            self.arguments_recon["voxel_size"]
        location[1] = location[1] * \
            self.arguments_recon["voxel_size"]
        location[2] = location[2] * \
            self.arguments_recon["voxel_size"]

        # cropped phantom length in Y dimension (mm)
        crop_phan_lenY_mm = self.arguments_recon["voxels_x"] * \
            self.arguments_recon["voxel_size"]
        # detector length in Y dimension (mm)
        det_lenY_mm = self.arguments_recon["detector_elements"] * \
            self.arguments_recon["pixel_size"]

        alpha = (detector_z -
                 self.arguments_mcgpu["source_position"][2]) / \
            (location[2] - self.arguments_mcgpu["source_position"][2])

        location[0] = self.arguments_mcgpu["source_position"][0] + \
            alpha * (location[0] -
                     self.arguments_mcgpu["source_position"][0])
        location[1] = self.arguments_mcgpu["source_position"][1] + \
            alpha * (location[1] -
                     self.arguments_mcgpu["source_position"][1])

        det_origin = [0,
                      ((det_lenY_mm - crop_phan_lenY_mm) * 0.5) /
                      self.arguments_recon["pixel_size"]]

        location[0] = int(location[0] / self.arguments_recon["pixel_size"])
        location[1] = int(location[1] / self.arguments_recon["pixel_size"])

        location[0] += det_origin[0]
        location[1] += det_origin[1]

        # we figured out by looking at the voxels and pixels that Y
        location[1] = self.arguments_recon["detector_elements"] - location[1]

        return location[0], location[1]

    def save_DICOM(self):
        """!
            Saves the DBT generated reconstruction (if available) in DICOM format
        """
        def save_DICOM_one(data, count):
            # Populate required values for file meta information
            file_meta = FileMetaDataset()
            file_meta.MediaStorageSOPClassUID = pydicom.uid.generate_uid()
            # file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
            file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()
            file_meta.ImplementationClassUID = pydicom.uid.PYDICOM_IMPLEMENTATION_UID

            # file_meta.ClinicalTrialProtocolName = 'VICTRE'
            # file_meta.ClinicalTrialSiteName = 'FDA'

            # file_meta.Manufacturer = 'VICTRE'
            # file_meta.OrganExposed = 'BREAST'
            # file_meta.Modality = "DBT"

            # Create the FileDataset instance (initially no data elements, but file_meta
            # supplied)
            ds = FileDataset("{:s}/{:d}/DICOM/{:03d}.dcm".format(self.results_folder, self.seed, count), {},
                             file_meta=file_meta, preamble=b"\0" * 128)

            # Add the data elements -- not trying to set all required here. Check DICOM
            # standard
            # ds.SamplesPerPixel = 1
            # ds.PhotometricInterpretation = "MONOCHROME2"
            ds.PixelRepresentation = 0
            ds.HighBit = 15
            ds.BitsStored = 16
            ds.BitsAllocated = 16
            ds.SmallestImagePixelValue = 0
            ds[0x00280106].VR = 'US'
            ds.LargestImagePixelValue = 65535
            ds[0x00280107].VR = 'US'

            ds.PatientName = "VICTRE 2.0"
            ds.PatientID = str(self.seed)

            ds.Columns = data.shape[0]
            ds.Rows = data.shape[1]

            # ds.ImagesInAcquisition = self.recon_size["z"]

            ds.PixelData = data.tobytes()

            # Set the transfer syntax
            # ds.is_little_endian = True
            # ds.is_implicit_VR = True

            # Set creation date/time
            dt = datetime.datetime.now()
            ds.ContentDate = dt.strftime('%Y%m%d')
            # long format with micro seconds
            timeStr = dt.strftime('%H%M%S.%f')
            ds.ContentTime = timeStr

            # print("Writing test file",
            #       "./results/{:d}/DICOM/{:03d}.dcm".format(self.seed, count))
            # ds.save_as("./results/{:d}/DICOM/{:03d}.dcm".format(self.seed, count))

            pydicom.filewriter.dcmwrite(
                "{:s}/{:d}/DICOM/{:03d}.dcm".format(self.results_folder, self.seed, count), ds)

        os.makedirs("{:s}/{:d}/DICOM/".format(self.results_folder,
                                              self.seed), exist_ok=True)

        pixel_array = np.fromfile("{:s}/{:d}/reconstruction{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                  dtype="float64").reshape(self.recon_size["z"], self.recon_size["x"], self.recon_size["y"])

        pixel_array = (2**16 * (pixel_array - pixel_array.min()) /
                       (pixel_array.max() - pixel_array.min())).astype(np.uint16)

        for s in progressbar.progressbar(range(pixel_array.shape[0])):
            save_DICOM_one(np.squeeze(pixel_array[s, :, :]), s)

    def save_ROIs(self, roi_sizes=None, clean=True):
        """!
            Saves the generated ROIs (absent and present) in RAW and HDF5 formats

            @param roi_sizes Size of the ROIs for the defined lesion types
            @param clean If True, the existing ROI folder will be deleted
        """

        if len(self.lesion_locations["dbt"]) == 0:
            cprint("There are no ROIs!", 'red')
            return

        if roi_sizes is not None:
            self.roi_sizes = roi_sizes

        # SAVE DBT ROIs
        if clean:
            shutil.rmtree(
                "{:s}/{:d}/ROIs".format(self.results_folder, self.seed), ignore_errors=True)

        os.makedirs("{:s}/{:d}/ROIs/".format(self.results_folder,
                                             self.seed), exist_ok=True)

        hf = h5py.File(
            "{:s}/{:d}/ROIs.h5".format(self.results_folder, self.seed), 'w')
        hfdbt = hf.create_group("dbt")

        pixel_array = np.fromfile("{:s}/{:d}/reconstruction{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                  dtype="float64").reshape(self.recon_size["z"], self.recon_size["y"], self.recon_size["x"])

        for idx, lesion in enumerate(self.lesion_locations["dbt"]):
            lesion_type = np.abs(lesion[3])
            roi = pixel_array[1 + lesion[2] - int(np.ceil(self.roi_sizes[lesion_type][2] / 2)):1 + lesion[2] + int(np.floor(self.roi_sizes[lesion_type][2] / 2)),
                              lesion[1] - int(np.ceil(self.roi_sizes[lesion_type][1] / 2)):lesion[1] + int(np.floor(self.roi_sizes[lesion_type][1] / 2)),
                              lesion[0] - int(np.ceil(self.roi_sizes[lesion_type][0] / 2)):lesion[0] + int(np.floor(self.roi_sizes[lesion_type][0] / 2))]
            # with open("./results/{:d}/ROIs/ROI_{:03d}_type{:d}.raw".format(self.seed, idx, lesion_type), 'wb') as f:
            roi.astype(np.dtype('<f8')).tofile(
                "{:s}/{:d}/ROIs/ROI_DBT_{:02d}_type{:d}.raw".format(self.results_folder, self.seed, idx, lesion[3]))
            hfdbt.create_dataset("{:d}".format(idx),
                                 data=roi.astype(np.float32), compression="gzip", compression_opts=9)
        hfdbt.create_dataset("lesion_type",
                             data=np.array(self.lesion_locations["dbt"])[:, 3])

        # SAVE DM ROIs
        pixel_array = np.fromfile("{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                  dtype="float32").reshape(2,
                                                           self.arguments_recon["detector_elements_perpendicular"],
                                                           self.arguments_recon["detector_elements"])

        hfdm = hf.create_group("dm")

        for idx, lesion in enumerate(self.lesion_locations["dm"]):
            lesion_type = np.abs(lesion[2])
            roi = pixel_array[0,
                              lesion[0] - int(np.ceil(self.roi_sizes[lesion_type][0] / 2)):lesion[0] + int(np.floor(self.roi_sizes[lesion_type][0] / 2)),
                              lesion[1] - int(np.ceil(self.roi_sizes[lesion_type][1] / 2)):lesion[1] + int(np.floor(self.roi_sizes[lesion_type][1] / 2))]
            # with open("./results/{:d}/ROIs/ROI_{:03d}_type{:d}.raw".format(self.seed, idx, lesion_type), 'wb') as f:
            roi.tofile(
                "{:s}/{:d}/ROIs/ROI_DM_{:02d}_type{:d}.raw".format(self.results_folder, self.seed, idx, lesion[2]))
            # dm_rois.append(roi)

            hfdm.create_dataset("{:d}".format(idx),
                                data=roi, compression="gzip", compression_opts=9)

        hfdm.create_dataset("lesion_type",
                            data=np.array(self.lesion_locations["dm"])[:, 2])

        hf.close()

        cprint("ROIs saved!", 'green', attrs=['bold'])

    def generate_spiculated(self, seed, size):
        """!
            Generates a spiculated mass using the breastMass software

            @param seed Seed to be used when generating the mass
            @param size Size of the mass to be used in the breastMass config file
            @return None. The result is saved in the `lesions` subfolder
        """
        with open("./Victre/breastMass/configs/spiculated.tpl", "r") as f:
            src = Template(f.read())
            result = src.substitute({"alpha": size, "seed": seed})

        os.makedirs("./lesions/", exist_ok=True)
        os.makedirs("./lesions/spiculated/", exist_ok=True)

        with open("./lesions/spiculated/input_breastMass_{:d}.in".format(seed), "w") as f:
            f.write(result)

        command = "cd ./lesions/spiculated/ && ../../Victre/breastMass/breastMass -c input_breastMass_{:d}.in".format(
            seed)

        # print(ssh_command)
        cprint("Generating mass (seed={:d}, size={:.2f})...".format(
            seed, size), 'cyan')
        os.system(command)

        generated_files = self.get_folder_contents("./lesions/spiculated/")

        side = None
        for name in generated_files:
            s = re.search(".*\/mass_{:d}_([0-9]*)\.raw".format(seed), name)
            if s is not None:
                side = int(s[1])

        lesion_raw = np.fromfile("./lesions/spiculated/mass_{:d}_{:d}.raw".format(seed, side),
                                 dtype=np.int8).reshape(side, side, side)

        # clean files
        for name in generated_files:
            if ".raw" in name or ".cfg" in name or ".vti" in name or ".in" in name or "core" in name:
                os.remove(name)

        # save in HDF
        with h5py.File("./lesions/spiculated/mass_{:d}_size{:.2f}.h5".format(seed, size), "w") as hf:
            hf.create_dataset("volume", data=lesion_raw, compression="gzip")
            hf.create_dataset("seed", data=seed)
            hf.create_dataset("size", data=size)

        cprint("Generation finished!", 'green', attrs=['bold'])

    def insert_lesions(self, lesion_type=None, n=-1, lesion_file=None, lesion_size=None, locations=None, roi_sizes=None, save_phantom=True):
        """!
            Inserts the specified number of lesions in the phantom.

            @param lesion_type Constant with the desired lesion type. Check available lesion types and materials in the Constants file.
            @param n Number of lesions to be added
            @param lesion_file Path to file including the lesion to be inserted (in HDF5 format). If specified, it will overrite the lesion file specified in the constructor.
            @param lesion_size If lesion_file is a raw file, lesion_size indicates the size of this file
            @param locations List of coordinates in the voxel/phantom space where the lesions will be inserted. If not specified, random locations will be generated.
            @param roi_sizes Size of the region of interest to be calculated to avoid overlapping with other tissues and check out of bounds locations

            @return None. A phantom file will be saved inside the results folder with the corresponding raw phantom. Three files will be generated: `pcl_SEED.raw.gz` with the raw data, `pcl_SEED.mhd` with the information about the raw data, and `pcl_SEED.loc` with the voxel coordinates of the lesion centers.

        """
        if self.lesion_file is None and lesion_file is None and save_phantom is True:
            cprint(
                "There is no lesion to insert, just adding lesion locations...", color="cyan")
            # raise Exceptions.VictreError("No lesion file has been specified")

        lesion, phantom = None, None

        if lesion_file is not None:
            self.lesion_file = lesion_file

        # read self.arguments_mcgpu compressed
        with gzip.open(self.arguments_mcgpu["phantom_file"], 'rb') as f:
            phantom = f.read()

        phantom = np.fromstring(phantom, dtype=np.uint8).reshape(
            self.arguments_mcgpu["number_voxels"][2],
            self.arguments_mcgpu["number_voxels"][1],
            self.arguments_mcgpu["number_voxels"][0])

        if self.lesion_file is not None:
            if n == -1:
                if locations is not None:
                    n = len(locations)
                else:
                    n = 1

            if save_phantom:
                cprint("Inserting {:d} non-overlapping lesions...".format(n),
                       'cyan')
            else:
                cprint("Retrieving {:d} lesion locations...".format(n), 'cyan')

            if "h5" in self.lesion_file:
                with h5py.File(self.lesion_file, "r") as hf:
                    lesion = hf["volume"][()]
            else:  # raw
                with open(self.lesion_file, "rb") as f:
                    lesion = f.read()
                lesion = np.fromstring(
                    lesion, dtype=np.uint8).reshape(lesion_size)

        if roi_sizes is None and lesion is not None:
            roi_shape = lesion.shape
        elif roi_sizes is not None:
            self.roi_sizes = roi_sizes

        if locations is not None:
            for cand in locations:
                cand_type = lesion_type
                if cand_type is None:
                    cand_type = cand[3]

                if lesion is None:
                    lesion_shape = self.roi_sizes[np.abs(cand_type)]

                if lesion is not None:
                    lesion_shape = lesion.shape
                    roi = phantom[int(cand[0] - lesion_shape[0] / 2):int(cand[0] + lesion_shape[0] / 2),
                                  int(cand[2] - lesion_shape[2] / 2):int(cand[2] + lesion_shape[2] / 2),
                                  int(cand[1] - lesion_shape[1] / 2):int(cand[1] + lesion_shape[1] / 2)]
                    phantom[int(cand[0] - lesion_shape[0] / 2):int(cand[0] + lesion_shape[0] / 2),
                            int(cand[2] - lesion_shape[2] / 2):int(cand[2] + lesion_shape[2] / 2),
                            int(cand[1] - lesion_shape[1] / 2):int(cand[1] + lesion_shape[1] / 2)][lesion == 1] = Constants.LESION_MATERIALS[np.abs(cand_type)]

                self.lesions.append(np.array([cand[0],
                                              cand[1],
                                              cand[2],
                                              cand_type
                                              ]))
        else:
            if self.candidate_locations is not None:
                # from mm to voxels
                for idx, cand in enumerate(self.candidate_locations):
                    self.candidate_locations[idx] = [int(np.round((cand[0] - self.mhd["Offset"][0]) /
                                                                  self.mhd["ElementSpacing"][0])),
                                                     int(np.round((cand[2] - self.mhd["Offset"][2]) /
                                                                  self.mhd["ElementSpacing"][2])),
                                                     int(np.round((cand[1] - self.mhd["Offset"][1]) /
                                                                  self.mhd["ElementSpacing"][1]))]
                Constants.INSERTION_MAX_TRIES = len(self.candidate_locations)
                Constants.INSERTION_MAX_TOTAL_ATTEMPTS = 1000
                # current_candidate = 0

            roi_shape = self.roi_sizes[lesion_type]
            c = 0
            current_seed = self.seed
            random.seed(current_seed)
            max_attempts = Constants.INSERTION_MAX_TOTAL_ATTEMPTS
            while c < n and max_attempts >= 0:
                found = False
                roi = None
                cand = None
                loc = None
                attempts = 0
                bar = progressbar.ProgressBar(
                    max_value=Constants.INSERTION_MAX_TRIES)
                while not found and max_attempts > 0:
                    attempts += 1
                    bar.update(attempts)
                    if attempts == bar.max_value:  # if too many attempts
                        bar.finish()
                        attempts = 0
                        max_attempts -= 1

                        cprint(
                            "Too many attempts at inserting, restarting the insertion! ({:d} remaining)".format(max_attempts), 'red')
                        with gzip.open(self.arguments_mcgpu["phantom_file"], 'rb') as f:
                            phantom = f.read()

                        phantom = np.fromstring(phantom, dtype=np.uint8).reshape(
                            self.arguments_mcgpu["number_voxels"][2],
                            self.arguments_mcgpu["number_voxels"][1],
                            self.arguments_mcgpu["number_voxels"][0])
                        current_seed += 1
                        random.seed(current_seed)  # try with a different seed

                        # rollback
                        self.lesions = self.lesions[:-c]
                        c = 0

                        if self.candidate_locations is not None:
                            np.random.shuffle(self.candidate_locations)

                        if max_attempts == 0:
                            raise Exceptions.VictreError(
                                "Insertion attempts exceeded")

                        bar = progressbar.ProgressBar(max_value=bar.max_value)
                        continue

                    if self.candidate_locations is not None:
                        cand = (
                            self.candidate_locations[attempts] - np.array(lesion.shape) / 2).astype(int)
                    else:
                        cand = [
                            random.randint(0, phantom.shape[0] - roi_shape[0]),
                            random.randint(0, phantom.shape[2] - roi_shape[2]),
                            random.randint(0, phantom.shape[1] - roi_shape[1])]

                    loc = {"dm": self.get_coordinates_dm([
                        cand[1] + lesion.shape[1] / 2,
                        cand[2] + lesion.shape[2] / 2,
                        cand[0] + lesion.shape[0] / 2]),
                        "dbt": self.get_coordinates_dbt([
                            cand[1] + lesion.shape[1] / 2,
                            cand[2] + lesion.shape[2] / 2,
                            cand[0] + lesion.shape[0] / 2])}

                    if np.any(np.array(loc["dm"]) < np.array(roi_shape[:2])) or \
                       np.any(np.array(loc["dbt"]) < np.array(roi_shape)):
                        continue

                    roi = phantom[cand[0]:cand[0] + lesion.shape[0],
                                  cand[2]:cand[2] + lesion.shape[2],
                                  cand[1]:cand[1] + lesion.shape[1]]

                    # check if lesion volume is too close to air, skin, nipple and muscle
                    if not np.any(np.array(roi.shape) < lesion.shape) and \
                        not (np.any([np.any(roi == x) for x in np.append(Constants.FORBIDDEN_OVERLAP,
                                                                         list(Constants.LESION_MATERIALS.values()))])):
                        found = True

                phantom[cand[0]:cand[0] + lesion.shape[0],
                        cand[2]:cand[2] + lesion.shape[2],
                        cand[1]:cand[1] + lesion.shape[1]][lesion == 1] = Constants.LESION_MATERIALS[lesion_type]

                self.lesions.append(np.array([int(cand[0] + lesion.shape[0] / 2),
                                              int(cand[1] +
                                                  lesion.shape[1] / 2),
                                              int(cand[2] +
                                                  lesion.shape[2] / 2),
                                              lesion_type
                                              ]))

                c += 1

                bar.finish()

        for cand in self.lesions:
            loc = {"dm": self.get_coordinates_dm([
                cand[1],
                cand[2],
                cand[0]]),
                "dbt": self.get_coordinates_dbt([
                    cand[1],
                    cand[2],
                    cand[0]])}
            self.lesion_locations["dm"].append(
                list(np.round([loc["dm"][0], loc["dm"][1], cand[3]]).astype(int)))

            self.lesion_locations["dbt"].append(
                list(np.round([loc["dbt"][0], loc["dbt"][1], loc["dbt"][2], cand[3]]).astype(int)))

        if lesion is not None:
            np.savetxt("{:s}/{:d}/pcl_{:d}.loc".format(self.results_folder, self.seed, self.seed),
                       np.asarray(self.lesions), fmt="%d")

            # with h5py.File("phantom/pcl_{:d}_crop.h5".format(self.seed), "w") as hf:
            #     hf.create_dataset("phantom", data=phantom.astype(
            #         np.uint8), compression="gzip")

            # save new phantom file
            if save_phantom:
                cprint("Saving new phantom...", 'cyan')

                # We save the phantom in gzip to reduce needed disk space
                with gzip.open("{:s}/{:d}/pcl_{:d}.raw.gz".format(self.results_folder, self.seed, self.seed), "wb") as gz:
                    gz.write(phantom)

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pcl_{:d}.raw.gz".format(
                    self.results_folder, self.seed, self.seed)

                with open("{:s}/{:d}/pcl_{:d}.mhd".format(self.results_folder, self.seed, self.seed), "w") as f:
                    src = Template(Constants.MHD_FILE)
                    template_arguments = copy.deepcopy(self.mhd)
                    template_arguments["ElementDataFile"] = "pcl_{:d}.raw.gz".format(
                        self.seed)
                    for key in template_arguments.keys():
                        if type(template_arguments[key]) is list:
                            template_arguments[key] = ' '.join(
                                map(str, template_arguments[key]))
                    result = src.substitute(template_arguments)
                    f.write(result)

                cprint("Insertion finished!", 'green', attrs=['bold'])

    def add_absent_ROIs(self, lesion_type, n=1, locations=None, roi_sizes=None):
        """!
            Adds the specified number of absent regions of interest.

            @param lesion_type Constant with the desired lesion type. Check available lesion types and materials in the Constants file.
            @param n Number of lesions to be added
            @param locations List of coordinates in the voxel/phantom space where the lesions will be inserted. If not specified, random locations will be generated.
            @param roi_sizes Size of the region of interest to be calculated to avoid overlapping with other tissues and check out of bounds locations
            @return None. A location file will be saved inside the `phantom` folder with the corresponding seed. Negative lesion type means absent ROI.
        """
        with gzip.open(self.arguments_mcgpu["phantom_file"], 'rb') as f:
            phantom = f.read()

        if roi_sizes is not None:
            self.roi_sizes = roi_sizes

        roi_shape = self.roi_sizes[lesion_type]

        phantom = np.fromstring(phantom, dtype=np.uint8).reshape(
            self.arguments_mcgpu["number_voxels"][2],
            self.arguments_mcgpu["number_voxels"][1],
            self.arguments_mcgpu["number_voxels"][0])

        if locations is not None:
            for cand in locations:
                roi = phantom[int(cand[0] - roi_shape[0] / 2):int(cand[0] + roi_shape[0] / 2),
                              int(cand[2] - roi_shape[2] / 2):int(cand[2] + roi_shape[2] / 2),
                              int(cand[1] - roi_shape[1] / 2):int(cand[1] + roi_shape[1] / 2)]

                self.lesions.append(np.array([cand[0],
                                              cand[1],
                                              cand[2],
                                              -lesion_type  # -1 for absent
                                              ]))
        else:
            c = 0
            while c < n:
                found = False
                roi = None
                cand = None
                loc = None
                while not found:
                    cand = [
                        random.randint(0, phantom.shape[0] - roi_shape[0]),
                        random.randint(0, phantom.shape[2] - roi_shape[2]),
                        random.randint(0, phantom.shape[1] - roi_shape[1])]

                    loc = {"dm": self.get_coordinates_dm([cand[1] + roi_shape[1] / 2,
                                                          cand[2] +
                                                          roi_shape[2] / 2,
                                                          cand[0] + roi_shape[0] / 2]),
                           "dbt": self.get_coordinates_dbt([cand[1] + roi_shape[1] / 2,
                                                            cand[2] +
                                                            roi_shape[2] / 2,
                                                            cand[0] + roi_shape[0] / 2])}

                    if np.any(np.array(loc["dm"]) < np.array(roi_shape[:2])) or \
                       np.any(np.array(loc["dbt"]) < np.array(roi_shape)):
                        continue

                    roi = phantom[cand[0]:cand[0] + roi_shape[0],
                                  cand[2]:cand[2] + roi_shape[2],
                                  cand[1]:cand[1] + roi_shape[1]]

                    # check if lesion volume is too close to air, skin, nipple, muscle or lesion
                    if not (np.any([np.any(roi == x)
                                    for x in np.append(Constants.FORBIDDEN_OVERLAP,
                                                       Constants.LESION_MATERIALS)])):
                        found = True

                self.lesions.append(np.array([int(cand[0] + roi_shape[0] / 2),
                                              int(cand[1] + roi_shape[1] / 2),
                                              int(cand[2] + roi_shape[2] / 2),
                                              -lesion_type  # -1 for absent
                                              ]))

                self.lesion_locations["dm"].append(
                    list(np.round([loc["dm"][0], loc["dm"][1], -lesion_type]).astype(int)))

                self.lesion_locations["dbt"].append(
                    list(np.round([loc["dbt"][0], loc["dbt"][1], loc["dbt"][2], -lesion_type]).astype(int)))

                c += 1

        np.savetxt("{:s}/{:d}/pcl_{:d}.loc".format(self.results_folder, self.seed, self.seed),
                   self.lesions, fmt="%d")

    def generate_phantom(self):
        """!
            Runs breast phantom generation.

            @return None. A phantom file will be saved inside the results folder with the corresponding raw phantom. Two files will be generated: `p_SEED.raw.gz` with the raw data, and `p_SEED.mhd` with the information about the raw data.
        """
        generation_config = "{:s}/{:d}/input_generation.in".format(
            self.results_folder, self.seed)

        with open("./Victre/generation/configs/template_generation.tpl", "r") as f:
            src = Template(f.read())
            template_arguments = copy.deepcopy(self.arguments_generation)
            result = src.substitute(template_arguments)

        with open(generation_config, "w") as f:
            f.write(result)

        full_path = os.path.abspath(generation_config)

        command = "cd {:s} && ./Victre/generation/breastPhantomMain -c {:s}".format(
            os.getcwd(),
            full_path
        )

        if self.ips["cpu"] == "localhost":
            ssh_command = command
        else:
            ssh_command = "ssh -Y {:s} \"{:s}\"".format(
                self.ips["cpu"], command)

        cprint("Starting phantom generation (seed = {:d}), this will take some time...".format(
            self.seed), 'cyan')

        completed = 0

        process = subprocess.Popen(ssh_command, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)

        with open("{:s}/{:d}/output_generation.out".format(self.results_folder, self.seed), "wb") as f:
            while True:
                output = process.stdout.readline().decode("utf-8")
                f.write(output.encode('utf-8'))
                f.flush()
                if output == "" and process.poll() is not None:
                    break
                elif "Error extracting eigenfunctions" in output:
                    break
                completed += 1

        if not os.path.exists("{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, self.seed, self.seed)):
            cprint("\nError while generating, check the output_generation.out file in the results folder",
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Generation error")

        cprint("Generation finished!", 'green', attrs=['bold'])
        self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/p_{:d}.raw.gz".format(
            self.results_folder, self.seed, self.seed)

        self.mhd = self._read_mhd(
            "{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
        self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
        self.arguments_mcgpu["voxel_size"] = [
            x / 10 for x in self.mhd["ElementSpacing"]]
        self.arguments_recon["voxels_x"] = self.arguments_mcgpu["number_voxels"][1]
        self.arguments_recon["voxels_y"] = self.arguments_mcgpu["number_voxels"][0]
        self.arguments_recon["voxels_z"] = self.arguments_mcgpu["number_voxels"][2]

        self.arguments_recon["voxel_size"] = self.arguments_mcgpu["voxel_size"][0]

        self.recon_size = dict(
            x=np.ceil(self.arguments_recon["voxels_x"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            y=np.ceil(self.arguments_recon["voxels_y"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            z=np.ceil(self.arguments_recon["voxels_z"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_thickness"]).astype(int)
        )

        self.arguments_mcgpu["source_position"][1] = self.arguments_mcgpu["number_voxels"][1] * \
            self.arguments_mcgpu["voxel_size"][1] / 2

        self.lesions = []

        self.candidate_locations = np.loadtxt(
            "{:s}/{:d}/p_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

    def compress_phantom(self, thickness=None):
        """!
            Runs the FEBio compression.

            @param thickness Specifies the objective thickness for the phantom to be compressed (in cm)
            @return None. A phantom file will be saved inside the results folder with the corresponding raw phantom. Two files will be generated: `pc_SEED.raw.gz` with the raw data, and `pc_SEED.mhd` with the information about the raw data.
        """
        if thickness is None:
            thickness = int(self.arguments_generation["compressionThickness"])

        command = "cd {:s} && ./Victre/compression/build/breastCompressMain -s {:d} -t {:f} -d {:s}/{:d}".format(
            os.getcwd(),
            self.seed,
            thickness,
            self.results_folder,
            self.seed
        )

        if self.ips["cpu"] == "localhost":
            ssh_command = command
        else:
            ssh_command = "ssh -Y {:s} \"{:s}\"".format(
                self.ips["cpu"], command)

        cprint("Starting phantom compression, this will take some time...", 'cyan')

        completed = 0

        process = subprocess.Popen(ssh_command, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)

        with open("{:s}/{:d}/output_compression.out".format(self.results_folder, self.seed), "wb") as f:
            while True:
                output = process.stdout.readline().decode("utf-8")
                if output == "" and process.poll() is not None:
                    break

                completed += 1

                f.write(output.encode('utf-8'))
                f.flush()

        if completed == 0:
            cprint("\nError while compressing, check the output_compression.out file in the results folder",
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Generation error")

        cprint("Compression finished!", 'green', attrs=['bold'])
        self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}.raw.gz".format(
            self.results_folder, self.seed, self.seed)

        self.mhd = self._read_mhd(
            "{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
        self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]

        self.arguments_recon["voxels_x"] = self.arguments_mcgpu["number_voxels"][1]
        self.arguments_recon["voxels_y"] = self.arguments_mcgpu["number_voxels"][0]
        self.arguments_recon["voxels_z"] = self.arguments_mcgpu["number_voxels"][2]

        self.candidate_locations = np.loadtxt(
            "{:s}/{:d}/pc_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

    def crop(self, size=None):
        """!
            Runs breast phantom cropping.

            @return None. A phantom file will be saved inside the results folder with the corresponding raw phantom. Two files will be generated: `pc_SEED_crop.raw.gz` with the raw data, and `pc_SEED_crop.mhd` with the information about the raw data.
        """

        cprint("Cropping phantom...", 'cyan')

        with gzip.open(self.arguments_mcgpu["phantom_file"], 'rb') as f:
            phantom = f.read()

        phantom = np.fromstring(phantom, dtype=np.uint8).reshape(
            self.arguments_mcgpu["number_voxels"][2],
            self.arguments_mcgpu["number_voxels"][1],
            self.arguments_mcgpu["number_voxels"][0])

        # crop from top to bottom (and bottom to top) when the plates start/end
        crop = {"from": [0, 0, 0], "to": list(phantom.shape)}
        for x in range(phantom.shape[0]):
            if(np.any(phantom[x, :, -1] == 50)):
                crop["from"][0] = x
                break
        for x in range(phantom.shape[0] - 1, 0, -1):
            if(np.any(phantom[x, :, -1] == 50)):
                crop["to"][0] = x
                break

        # crop from pectoral muscle towards nipple when the plates start
        for z in range(phantom.shape[1]):
            if(np.any(phantom[crop["to"][0], :, z] == 50)):
                crop["from"][2] = z
                break
            if(np.any(phantom[crop["from"][0], :, z] == 50)):
                crop["from"][2] = z
                break

        phantom = phantom[crop["from"][0]:crop["to"][0],
                          crop["from"][1]:crop["to"][1],
                          crop["from"][2]:crop["to"][2]]

        self.arguments_mcgpu["number_voxels"] = [phantom.shape[2],
                                                 phantom.shape[1],
                                                 phantom.shape[0]]

        self.mhd["DimSize"] = self.arguments_mcgpu["number_voxels"]

        self.arguments_recon["voxels_x"] = self.arguments_mcgpu["number_voxels"][1]
        self.arguments_recon["voxels_y"] = self.arguments_mcgpu["number_voxels"][0]
        self.arguments_recon["voxels_z"] = self.arguments_mcgpu["number_voxels"][2]

        with gzip.open("{:s}/{:d}/pc_{:d}_crop.raw.gz".format(self.results_folder, self.seed, self.seed), 'wb') as f:
            f.write(np.ascontiguousarray(phantom))

        self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}_crop.raw.gz".format(
            self.results_folder, self.seed, self.seed)

        prevOffset = copy.deepcopy(self.mhd["Offset"])

        self.mhd["ElementDataFile"] = "pc_{:d}_crop.raw.gz".format(
            self.seed)
        self.mhd["Offset"][0] = self.mhd["Offset"][0] + \
            crop["from"][0] * self.mhd["ElementSpacing"][0]
        self.mhd["Offset"][1] = self.mhd["Offset"][1] + \
            crop["from"][1] * self.mhd["ElementSpacing"][1]
        self.mhd["Offset"][2] = self.mhd["Offset"][2] + \
            crop["from"][2] * self.mhd["ElementSpacing"][2]

        with open("{:s}/{:d}/pc_{:d}_crop.mhd".format(self.results_folder, self.seed, self.seed), "w") as f:
            src = Template(Constants.MHD_FILE)
            template_arguments = copy.deepcopy(self.mhd)
            for key in template_arguments.keys():
                if type(template_arguments[key]) is list:
                    template_arguments[key] = ' '.join(
                        map(str, template_arguments[key]))
            result = src.substitute(template_arguments)
            f.write(result)

        if self.candidate_locations is not None:
            for cand in self.candidate_locations:
                cand[0] = ((cand[0] - prevOffset[0]) / self.mhd["ElementSpacing"][0] -
                           crop["from"][0]) * self.mhd["ElementSpacing"][0] + self.mhd["Offset"][0]
                cand[1] = ((cand[1] - prevOffset[1]) / self.mhd["ElementSpacing"][1] -
                           crop["from"][1]) * self.mhd["ElementSpacing"][1] + self.mhd["Offset"][1]
                cand[2] = ((cand[2] - prevOffset[2]) / self.mhd["ElementSpacing"][2] -
                           crop["from"][2]) * self.mhd["ElementSpacing"][2] + self.mhd["Offset"][2]
        np.savetxt("{:s}/{:d}/pc_{:d}_crop.loc".format(self.results_folder,
                                                       self.seed,
                                                       self.seed),
                   self.candidate_locations,
                   delimiter=',')

    @staticmethod
    def get_folder_contents(folder):
        """!
            Gets a list of files in the given folder

            @param folder Path to the folder to be processed
            @return List with files inside the given folder
        """
        dir_folder = pathlib.Path(folder)
        files = []
        for currentFile in dir_folder.iterdir():
            files.append(join(folder, currentFile.name))

        return files

    @staticmethod
    def _read_mhd(filename):
        data = {}
        with open(filename, "r") as f:
            for line in f:
                s = re.search(
                    "([a-zA-Z]*) = (.*)", line)
                data[s[1]] = s[2]

                if " " in data[s[1]]:
                    data[s[1]] = data[s[1]].split(' ')
                    for i in range(len(data[s[1]])):
                        if data[s[1]][i].replace(".", "").replace("-", "").isnumeric():
                            if "." in data[s[1]][i]:
                                data[s[1]][i] = float(data[s[1]][i])
                            else:
                                data[s[1]][i] = int(data[s[1]][i])
                else:
                    if data[s[1]].replace(".", "").replace("-", "").isnumeric():
                        if "." in data[s[1]]:
                            data[s[1]] = float(data[s[1]])
                        else:
                            data[s[1]] = int(data[s[1]])
        return data
