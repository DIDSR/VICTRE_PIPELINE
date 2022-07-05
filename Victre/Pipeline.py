"""
==================================================
                VICTRE PIPELINE
==================================================

 Author: Miguel A. Lago
		 miguel.lago@fda.hhs.gov


				    DISCLAIMER

 This software and documentation (the "Software") were
 developed at the Food and Drug Administration (FDA) by
 employees of the Federal Government in the course of
 their official duties. Pursuant to Title 17, Section
 105 of the United States Code, this work is not subject
 to copyright protection and is in the public domain.
 Permission is hereby granted, free of charge, to any
 person obtaining a copy of the Software, to deal in the
 Software without restriction, including without
 limitation the rights to use, copy, modify, merge,
 publish, distribute, sublicense, or sell copies of the
 Software or derivatives, and to permit persons to whom
 the Software is furnished to do so. FDA assumes no
 responsibility whatsoever for use by other parties of
 the Software, its source code, documentation or compiled
 executables, and makes no guarantees, expressed or
 implied, about its quality, reliability, or any other
 characteristic. Further, use of this code in no way
 implies endorsement by the FDA or confers any advantage
 in regulatory decisions. Although this software can be
 redistributed and/or modified freely, we ask that any
 derivative works bear some notice that they are derived
 from it, and any modified versions bear some notice that
 they have been modified.

 More information: https://github.com/DIDSR/VICTRE_PIPELINE

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
from datetime import date
from string import Template
import random
import time
from . import Constants, Exceptions
import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
import copy
import datetime
from pydicom.encaps import encapsulate
import re
import gzip
from scipy import interpolate


class Pipeline:
    """
        Object constructor for the Victre pipeline class

        :param ips: Dictionary with two IP addresses to run the pipeline: "gpu" for the projection process. "cpu" for the reconstruction.
        :param seed: Random seed used to generate or read the phantom
        :param results_folder: Path to folder to be used when saving the results
        :param phantom_file: Path to file containing the phantom to be loaded
        :param spectrum_file: Path to file containing the spectrum used to project in MCGPU
        :param lesion_file: Path to file containing the lesion to be inserted (in HDF5 format)
        :param materials: Dictionary including the materials to be used during projection
        :param roi_sizes: Dictionary with the ROI sizes for the extraction
        :param arguments_generation: Arguments to be overriden for the breast phantom generation
        :param arguments_mcgpu: Arguments to be overridden for the projection in MCGPU
        :param arguments_recon: Arguments to be overridden for the reconstruction algorithm
        :param flatfield_DBT: Path to the flatfield file for the DBT reconstruction
        :param flatfield_DM: Path to the flatfield file for the digital mammography
        :param density: [EXPERIMENTAL] Percentage of dense tissue of the phantom to be generated, this will adjust the compression thickness too
        :param verbosity: True will output the progress of each process and steps
        :returns: None
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
                 arguments_spiculated=dict(),
                 arguments_mcgpu=dict(),
                 arguments_recon=dict(),
                 flatfield_DBT=None,
                 flatfield_DM=None,
                 density=None,
                 verbosity=True):

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
        self.verbosity = verbosity

        random.seed(self.seed)

        self.arguments_mcgpu = Constants.VICTRE_DEFAULT_MCGPU
        self.arguments_mcgpu["spectrum_file"] = spectrum_file
        self.arguments_mcgpu["phantom_file"] = phantom_file
        self.arguments_mcgpu["output_file"] = "{:s}/{:d}/projection".format(
            self.results_folder, self.seed)
        self.arguments_mcgpu["random_seed"] = self.seed

        self.arguments_spiculated = Constants.VICTRE_DEFAULT_SPICULATED_MASS
        self.arguments_spiculated["seed"] = self.seed

        locations = None
        self.mhd = {
            "ObjectType": "Image",
            "NDims": 2,
            "BinaryData": "True",
            "BinaryDataByteOrderMSB": "False",
            "CompressedData": "False",
            "TransformMatrix": "1 0 0 0 1 0 0 0 1",
            "Offset": "0 0 0",
            "CenterOfRotation": "0 0 0",
            "ElementSpacing": "0.085 0.085",
            "DimSize": "3000 1500",
            "AnatomicalOrientation": "???",
            "ElementType": "MET_FLOAT",
            "ObjectType": "Image",
            "ElementDataFile": ""
        }

        if phantom_file is None:
            if os.path.exists("{:s}/{:d}/pcl_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found phantom with lesions information!",
                       'cyan') if self.verbosity else None
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pcl_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                locations = np.loadtxt(
                    "{:s}/{:d}/pcl_{:d}.loc".format(self.results_folder, self.seed, self.seed)).tolist()

                if os.path.exists("{:s}/{:d}/pc_{:d}_crop.loc".format(self.results_folder, seed, seed)):
                    self.candidate_locations = np.loadtxt(
                        "{:s}/{:d}/pc_{:d}_crop.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

                # from mm to voxels
                self.candidate_locations = self._mm_to_voxels(
                    self.candidate_locations)

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pcl_{:d}.raw.gz".format(
                    self.results_folder, seed, seed)

            elif os.path.exists("{:s}/{:d}/pc_{:d}_crop.mhd".format(self.results_folder, seed, seed)):
                cprint("Found cropped phantom information!",
                       'cyan') if self.verbosity else None
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pc_{:d}_crop.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                if os.path.exists("{:s}/{:d}/pc_{:d}_crop.loc".format(self.results_folder, seed, seed)):
                    self.candidate_locations = np.loadtxt(
                        "{:s}/{:d}/pc_{:d}_crop.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

                # from mm to voxels
                self.candidate_locations = self._mm_to_voxels(
                    self.candidate_locations)

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}_crop.raw.gz".format(
                    self.results_folder, seed, seed)
            elif os.path.exists("{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found compressed phantom information!",
                       'cyan') if self.verbosity else None
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                if os.path.exists("{:s}/{:d}/pc_{:d}.loc".format(self.results_folder, seed, seed)):
                    self.candidate_locations = np.loadtxt(
                        "{:s}/{:d}/pc_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

                # from mm to voxels
                self.candidate_locations = self._mm_to_voxels(
                    self.candidate_locations)

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}.raw.gz".format(
                    self.results_folder, seed, seed)

            elif os.path.exists("{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found phantom generation information!",
                       'cyan') if self.verbosity else None
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                self.candidate_locations = np.loadtxt(
                    "{:s}/{:d}/p_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

                # from mm to voxels
                self.candidate_locations = self._mm_to_voxels(
                    self.candidate_locations)

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/p_{:d}.raw.gz".format(
                    self.results_folder, seed, seed)

        self.arguments_mcgpu.update(arguments_mcgpu)

        # cm to mm
        self.arguments_spiculated["imgRes"] = self.arguments_mcgpu["voxel_size"][0] * 10

        self.arguments_spiculated.update(arguments_spiculated)

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
        self.arguments_recon["flatfield_file"] = flatfield_DBT

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

        if density is not None:
            fat = np.max([0.4, np.min([0.95, 1 - density])])
            ranges = {}
            for key in Constants.DENSITY_RANGES.keys():
                interp = interpolate.interp1d(
                    Constants.DENSITY_RANGES["targetFatFrac"], Constants.DENSITY_RANGES[key])
                ranges[key] = np.round(float(interp(fat)), 2)

            ranges["compartment_numBackSeeds"] = int(
                ranges["compartment_numBackSeeds"])  # this should be integer
            ranges["compartment_maxSkinScale"] = int(
                ranges["compartment_maxSkinScale"])  # this should be integer

            self.arguments_generation.update(ranges)
            if fat >= 0.75:  # increase the kVp when breast has low density
                # this is hardcoded here, careful
                self.arguments_mcgpu["spectrum_file"] = "./Victre/projection/spectrum/W30kVp_Rh50um_Be1mm.spc"
                self.arguments_mcgpu["fam_beam_aperture"][1] = 11.2

            self.arguments_mcgpu["number_histories"] = ranges["number_histories"]

        self.arguments_generation.update(arguments_generation)
        self.arguments_generation["seed"] = self.seed
        self.arguments_generation["outputDir"] = os.path.abspath("{:s}/{:d}/".format(
            self.results_folder, self.seed))

        self.recon_size = dict(
            x=np.ceil(self.arguments_recon["voxels_x"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            y=np.ceil(self.arguments_recon["voxels_y"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            z=np.ceil(self.arguments_recon["voxels_z"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_thickness"]).astype(int)
        )

        os.makedirs("{:s}".format(self.results_folder), exist_ok=True)
        os.makedirs("{:s}/{:d}".format(self.results_folder,
                                       self.seed), exist_ok=True)

        if phantom_file is not None:
            splitted = phantom_file.split('/')
            path = '/'.join(splitted[:-1])
            filename = splitted[-1].split('.')[0]

            # shutil.copy(phantom_file,
            #             "{:s}/{:d}".format(self.results_folder, self.seed))
            # os.chmod(
            #     "{:s}/{:d}/{:s}.raw.gz".format(self.results_folder, self.seed, filename), 0o664)

            if os.path.exists("{:s}/{:s}.mhd".format(path, filename)):
                cprint("Found phantom information!",
                       'cyan') if self.verbosity else None
                self.mhd = self._read_mhd(
                    "{:s}/{:s}.mhd".format(path, filename))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                # shutil.copy("{:s}/{:s}.mhd".format(path, filename),
                #             "{:s}/{:d}".format(self.results_folder, self.seed))
            if os.path.exists("{:s}/{:s}.loc".format(path, filename)):
                try:
                    locations = np.loadtxt(
                        "{:s}/{:s}.loc".format(path, filename))
                except:
                    pass
                # shutil.copy("{:s}/{:s}.loc".format(path, filename),
                #             "{:s}/{:d}".format(self.results_folder, self.seed))
                # os.chmod(
                #     "{:s}/{:d}/{:s}.loc".format(self.results_folder, self.seed, filename), 0o664)

        self.arguments_mcgpu["source_position"][1] = self.arguments_mcgpu["number_voxels"][1] * \
            self.arguments_mcgpu["voxel_size"][1] / 2

        if locations is not None:
            if not (type(locations[0]) is list or type(locations[0]) is np.ndarray):
                locations = [locations]
            self.insert_lesions(locations=locations,
                                save_phantom=False)

        # self.arguments_mcgpu["number_voxels"]

    def project(self, flatfield_correction=True, clean=True, do_flatfield=0, for_presentation=False):
        """
            Method that runs MCGPU to project the phantom.

            :param flatfield_correction: If True, the projections will be corrected using a given flatfield.
                                   It will be generated if not found and not given.
            :param clean: If True, it will delete the contents of the output folder before projecting.
            :param do_flatfield: If > 0, it will generate an empty flat field projection.
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
            del empty_phantom

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
        elif for_presentation:
            with gzip.open(self.arguments_mcgpu["phantom_file"], 'rb') as f:
                phantom = f.read()
            phantom = np.fromstring(phantom, dtype=np.uint8).reshape(
                self.arguments_mcgpu["number_voxels"][2],
                self.arguments_mcgpu["number_voxels"][1],
                self.arguments_mcgpu["number_voxels"][0])
            phantom[phantom != 0] = Constants.PHANTOM_MATERIALS["adipose"]
            with gzip.open("{:s}/{:d}/presentation.raw.gz".format(
                    self.results_folder, self.seed), "wb") as gz:
                gz.write(phantom)
            del phantom
            filename = "presentation"
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
            elif for_presentation:
                template_arguments["phantom_file"] = "{:s}/{:d}/presentation.raw.gz".format(
                    self.results_folder,
                    self.seed)
                # template_arguments["number_projections"] = 1
                # if self.arguments_mcgpu["number_projections"] > 1:
                #     template_arguments["number_histories"] = self.arguments_mcgpu["number_histories"] * \
                #         self.arguments_mcgpu["number_projections"] * 2 / 3
            template_arguments["output_file"] = "{:s}/{:d}/{:s}".format(
                self.results_folder, self.seed, filename)

            # from MBytes to Bytes and reduce 500MB for extra room
            # this would be for the first GPU
            gpu_ram = (get_gpu_memory()[0] - 500) * 1024 * 1024

            # if the binary tree has not been set
            if template_arguments["low_resolution_voxel_size"] == [0, 0, 0]:
                if gpu_ram < template_arguments["number_voxels"][0] * template_arguments["number_voxels"][1] * template_arguments["number_voxels"][2] or \
                        template_arguments["number_voxels"][0] * template_arguments["number_voxels"][1] * template_arguments["number_voxels"][2] > 2**32:
                    template_arguments["low_resolution_voxel_size"] = [1, 1, 1]

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

        cprint("Initializing MCGPU for {:s}...".format(
            filename), 'cyan') if self.verbosity else None

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
                        "[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Starting DM projection, this may take a few minutes...", 'cyan') if self.verbosity else None
                elif "Simulating tomographic projection" in output.strip():
                    if completed == 0:
                        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Starting DBT projection...",
                               'cyan') if self.verbosity else None
                        bar = progressbar.ProgressBar(
                            max_value=template_arguments["number_projections"]) if self.verbosity else None
                        bar.update(0) if self.verbosity else None
                    completed += 1
                    bar.update(completed) if self.verbosity else None
                # rc = process.poll()
                f.write(output.encode('utf-8'))
                f.flush()

        if template_arguments["number_projections"] > 1 and completed != template_arguments["number_projections"]:
            cprint("\nError while projecting, check the output_{:s}.out file in the results folder (seed = {:d})".format(filename, self.seed),
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Projection error")

        bar.finish() if bar is not None and self.verbosity else None
        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Projection finished!", 'green', attrs=[
               'bold']) if self.verbosity else None

        if self.arguments_mcgpu["number_projections"] > 1:
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

        if template_arguments["number_projections"] > 1:
            os.rename("{:s}/{:d}/{:s}_0000.raw".format(self.results_folder, self.seed, filename),
                      "{:s}/{:d}/{:s}_DM{:d}.raw".format(self.results_folder, self.seed, filename, self.seed))
        else:
            os.rename("{:s}/{:d}/{:s}.raw".format(self.results_folder, self.seed, filename),
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
        with contextlib.suppress(FileNotFoundError):
            os.remove(
                "{:s}/{:d}/{:s}".format(self.results_folder, self.seed, filename))

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

            if prev_flatfield_DBT is not None and self.arguments_mcgpu["number_projections"] > 1:
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

        elif flatfield_correction and (self.arguments_recon["flatfield_file"] is None or self.flatfield_DM is None):
            with contextlib.suppress(FileNotFoundError):
                os.remove(
                    "{:s}/{:d}/flatfield_DM{:d}.raw".format(self.results_folder, self.seed, self.seed).format(
                        self.results_folder,
                        self.seed,
                        'x'.join(
                            map(str, self.arguments_mcgpu["image_pixels"])),
                        self.arguments_mcgpu["number_projections"]))
                if self.arguments_mcgpu["number_projections"] > 1:
                    os.remove(
                        "{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                            self.results_folder,
                            self.seed,
                            'x'.join(
                                map(str, self.arguments_mcgpu["image_pixels"])),
                            self.arguments_mcgpu["number_projections"]))

            # number of iterations to average the flatfield
            for n in range(Constants.FLATFIELD_REPETITIONS):
                cprint("Flatfield files not specified, projecting {:d}/{:d}...".format(
                    n + 1, Constants.FLATFIELD_REPETITIONS), 'cyan') if self.verbosity else None
                self.project(do_flatfield=Constants.FLATFIELD_REPETITIONS)

            if self.arguments_mcgpu["number_projections"] > 1:
                self.arguments_recon["flatfield_file"] = "{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                    self.results_folder,
                    self.seed,
                    'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                    self.arguments_mcgpu["number_projections"])

            self.flatfield_DM = "{:s}/{:d}/flatfield_DM{:d}.raw".format(
                self.results_folder, self.seed, self.seed)

        if for_presentation:
            os.remove("{:s}/{:d}/presentation.raw.gz".format(
                self.results_folder, self.seed))
            self.project(flatfield_correction=False, clean=clean)
            projection_DM = np.fromfile("{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                        dtype="float32").reshape(2,
                                                                 self.arguments_recon["detector_elements_perpendicular"],
                                                                 self.arguments_recon["detector_elements"])
            presentation_tmp = np.fromfile("{:s}/{:d}/presentation_DM{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                           dtype="float32").reshape(2,
                                                                    self.arguments_recon["detector_elements_perpendicular"],
                                                                    self.arguments_recon["detector_elements"])
            os.remove(
                "{:s}/{:d}/presentation_DM{:d}.raw".format(self.results_folder, self.seed, self.seed))
            # os.rename(
            #     "{:s}/{:d}/presentation_DM{:d}.mhd".format(
            #         self.results_folder, self.seed, self.seed),
            #     "{:s}/{:d}/projection_DM{:d}.mhd".format(self.results_folder, self.seed, self.seed))
            np.seterr(divide='ignore', invalid='ignore')
            projection_DM = 1 / projection_DM / presentation_tmp
            projection_DM.tofile(
                "{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed))

        if do_flatfield == 0:
            # normalize with flatfield
            projection_DM = np.fromfile("{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                        dtype="float32").reshape(2,
                                                                 self.arguments_recon["detector_elements_perpendicular"],
                                                                 self.arguments_recon["detector_elements"])
            np.seterr(divide='ignore', invalid='ignore')
            if flatfield_correction and self.flatfield_DM is not None:
                curr_flatfield_DM = np.fromfile(self.flatfield_DM,
                                                dtype="float32").reshape(2,
                                                                         self.arguments_recon["detector_elements_perpendicular"],
                                                                         self.arguments_recon["detector_elements"])

                projection_DM = np.true_divide(
                    curr_flatfield_DM, projection_DM)
            else:
                projection_DM = np.true_divide(1, projection_DM)

            projection_DM[projection_DM == np.inf] = 0
            projection_DM[np.isnan(projection_DM)] = 0
            projection_DM.tofile(
                "{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed))
            if len(self.lesion_locations["dm"]) > 0:
                np.savetxt("{:s}/{:d}/projection_DM{:d}.loc".format(self.results_folder, self.seed, self.seed),
                           np.asarray(self.lesion_locations["dm"]), fmt="%d")

    def reconstruct(self):
        """
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

        self.recon_size = dict(
            x=np.ceil(self.arguments_recon["voxels_x"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            y=np.ceil(self.arguments_recon["voxels_y"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_pixel_size"]).astype(int),
            z=np.ceil(self.arguments_recon["voxels_z"] * self.arguments_recon["voxel_size"] /
                      self.arguments_recon["recon_thickness"]).astype(int)
        )

        cprint("Initializing reconstruction, this may take a few minutes...",
               'cyan') if self.verbosity else None

        completed = 0

        process = subprocess.Popen(ssh_command, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)

        bar = None
        finished = False
        with open("{:s}/{:d}/output_recon.out".format(self.results_folder, self.seed), "wb") as f:
            while True:
                output = process.stdout.readline().decode("utf-8")
                if output == "" and process.poll() is not None:
                    break
                elif "Image slice" in output.strip():
                    if completed == 0:
                        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Starting reconstruction...",
                               'cyan') if self.verbosity else None
                        bar = progressbar.ProgressBar(
                            max_value=self.recon_size["y"]) if self.verbosity else None
                        bar.update(0) if self.verbosity else None
                    completed += 1
                    bar.update(completed) if self.verbosity else None
                    progressbar.streams.flush() if self.verbosity else None
                elif "Total execution time elapsed" in output.strip():
                    finished = True
                # rc = process.poll()
                f.write(output.encode('utf-8'))
                f.flush()

        if not finished or completed != self.recon_size["y"]:
            cprint("\nError while reconstructing, check the output_recon.out file (seed = {:d})".format(self.seed),
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Reconstruction error")

        bar.finish() if self.verbosity else None

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

        if len(self.lesion_locations["dbt"]) > 0:
            np.savetxt("{:s}/{:d}/reconstruction{:d}.loc".format(self.results_folder, self.seed, self.seed),
                       np.asarray(self.lesion_locations["dbt"]), fmt="%d")

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Reconstruction finished!", 'green',
               attrs=['bold']) if self.verbosity else None

    def reverse_dm_coordinates(self, dm_location):
        location = dm_location.copy()

        location[1] = self.arguments_mcgpu["image_pixels"][0] - location[1]

        # cropped phantom length in Y dimension (mm)
        crop_phan_lenY_mm = self.arguments_recon["voxels_x"] * \
            self.arguments_recon["voxel_size"]
        # detector length in Y dimension (mm)
        det_lenY_mm = self.arguments_recon["detector_elements"] * \
            self.arguments_recon["pixel_size"]

        det_origin = [0,
                      ((det_lenY_mm - crop_phan_lenY_mm) * 0.5) /
                      self.arguments_recon["pixel_size"]]

        detector_z = - (self.arguments_mcgpu["distance_source"] -
                        self.arguments_mcgpu["source_position"][2])  # -20 / 10

        location[0] -= det_origin[0]
        location[1] -= det_origin[1]

        location[0] = location[0] * self.arguments_recon["pixel_size"]
        location[1] = location[1] * self.arguments_recon["pixel_size"]

        orig_location = location

        locations = []

        for s in range(self.arguments_mcgpu['number_voxels'][0]):
            location = [orig_location[0], orig_location[1],
                        s * self.arguments_recon["voxel_size"]]
            alpha = (detector_z -
                     self.arguments_mcgpu["source_position"][2]) / \
                (location[2] - self.arguments_mcgpu["source_position"][2])
            location[0] = (location[0] - self.arguments_mcgpu["source_position"][0]) / \
                alpha + self.arguments_mcgpu["source_position"][0]
            location[1] = (location[1] - self.arguments_mcgpu["source_position"][1]) / \
                alpha + self.arguments_mcgpu["source_position"][1]

            location[0] = location[0] / \
                self.arguments_recon["voxel_size"]
            location[1] = location[1] / \
                self.arguments_recon["voxel_size"]
            location[2] = location[2] / \
                self.arguments_recon["voxel_size"]

            location[2] = location[2] - self.arguments_recon["detector_offset"]
            location = [int(n) for n in location]
            # locations.append([int(n) for n in location])
            if location[2] < self.arguments_mcgpu['number_voxels'][2] \
                    and location[1] < self.arguments_mcgpu['number_voxels'][1]\
                    and location[0] < self.arguments_mcgpu['number_voxels'][0]\
                    and location[0] >= 0 and location[1] >= 0 and location[2] >= 0:
                locations.append(location)

        return locations

    def reverse_dbt_coordinates(self, dbt_location):
        location = dbt_location.copy()

        # interchange X and Y
        location[0], location[1] = location[1], location[0]

        # mirror Y axis
        location[1] = self.recon_size["x"] - location[1]

        location[0] = location[0] / self.arguments_recon["voxel_size"] * \
            self.arguments_recon["pixel_size"]
        location[1] = location[1] / self.arguments_recon["voxel_size"] * \
            self.arguments_recon["pixel_size"]
        location[2] = location[2] / self.arguments_recon["voxel_size"] * \
            self.arguments_recon["recon_thickness"]  # in mm

        location[2] = location[2] + self.arguments_recon["detector_offset"]

        return location

    def get_coordinates_dbt(self, vx_location):
        """
            Method to get the corresponding coordinates in the DBT volume from the voxelized coordinates

            :param vx_location: Coordinates in the voxel/phantom space
            :returns: Coordinates in the DBT space
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

        return [int(np.round(x)) for x in location]

    def get_coordinates_dm(self, vx_location):
        """
            Method to get the corresponding coordinates in the DM volume from the voxelized coordinates

            :param vx_location: Coordinates in the voxel/phantom space
            :returns: Coordinates in the DM space
        """
        location = vx_location.copy()

        detector_z = - (self.arguments_mcgpu["distance_source"] -
                        self.arguments_mcgpu["source_position"][2])  # -20 / 10

        location[2] = location[2] - self.arguments_recon["detector_offset"]

        pixel_size = self.arguments_mcgpu["image_size"][0] / \
            self.arguments_mcgpu["image_pixels"][0]

        location[0] = location[0] * self.arguments_mcgpu["voxel_size"][0]
        location[1] = location[1] * self.arguments_mcgpu["voxel_size"][1]
        location[2] = location[2] * self.arguments_mcgpu["voxel_size"][2]

        # cropped phantom length in Y dimension
        crop_phan_lenY = self.arguments_mcgpu["number_voxels"][1] * \
            self.arguments_mcgpu["voxel_size"][1]

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
                      ((self.arguments_mcgpu["image_size"][0] - crop_phan_lenY) * 0.5) /
                      pixel_size]

        location[0] = location[0] / pixel_size
        location[1] = location[1] / pixel_size

        location[0] += det_origin[0]
        location[1] += det_origin[1]

        # we figured out by looking at the voxels and pixels that Y
        location[1] = self.arguments_mcgpu["image_pixels"][0] - location[1]

        return [int(np.round(location[0])), int(np.round(location[1]))]

    def save_DICOM(self, modality="dbt"):
        """
            Saves the DM or DBT images in DICOM format. If present, lesion location will be
            stored in a custom tag 0x009900XX where XX is the lesion number.

            :param modality: Modality to save: dbt or dm
        """
        def save_DICOM_one(data, count):

            # Populate required values for file meta information
            file_meta = FileMetaDataset()
            file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
            # file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
            file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid(
                "1.3.6.1.4.1.9590.100.1.1.")
            file_meta.ImplementationClassUID = "1.3.6.1.4.1.9590.100.1.0.100.4.0"
            # file_meta.FileMetaInformationGroupLength = 208

            # pydicom.uid.PYDICOM_IMPLEMENTATION_UID

            # Create the FileDataset instance (initially no data elements, but file_meta
            # supplied)
            ds = FileDataset("{:s}/{:d}/DICOM/{:03d}.dcm".format(self.results_folder, self.seed, count), {},
                             file_meta=file_meta, preamble=b"\0" * 128)

            ds.SOPClassUID = file_meta.MediaStorageSOPClassUID
            ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID

            # Add the data elements -- not trying to set all required here. Check DICOM
            # standard
            ds.SamplesPerPixel = 1
            ds.PhotometricInterpretation = "MONOCHROME2"
            ds.PixelRepresentation = 0
            ds.HighBit = 15
            ds.BitsStored = 16
            ds.BitsAllocated = 16
            ds.SmallestImagePixelValue = 0

            ds[0x00280106].VR = 'US'
            ds.LargestImagePixelValue = 65535
            ds[0x00280107].VR = 'US'

            ds.ImageRotation = 90

            ds.Manufacturer = 'VICTRE'
            ds.OrganExposed = 'BREAST'
            ds.Modality = "MG"

            ds.PatientName = "VICTRE/FDA"
            ds.PatientID = str(self.seed)
            ds.PatientComments = 'Density: {:.2f}%'.format(
                (1 - self.arguments_generation["targetFatFrac"]) * 100)
            ds.PatientState = "No lesions" if len(
                self.lesion_locations[modality]) == 0 else "With lesions"

            ds.ClinicalTrialProtocolName = "VICTRE"
            ds.ClinicalTrialSiteName = "FDA"

            ds.AccessionNumber = ' '
            ds.AcquisitionContextSequence = ''
            ds.AnatomicRegionSequence = ''
            ds.BurnedInAnnotation = 'NO'
            ds.ClinicalTrialProtocolID = ' '
            ds.ClinicalTrialSiteID = ' '
            ds.ClinicalTrialSponsorName = ' '
            ds.ClinicalTrialSubjectID = ' '
            ds.ClinicalTrialSubjectReadingID = ' '
            ds.ImageLaterality = 'R'
            ds.ImagerPixelSpacing = "{:f}\{:f}".format(
                self.arguments_recon["recon_pixel_size"] * 10, self.arguments_recon["recon_pixel_size"] * 10)
            ds.InstanceNumber = ''
            ds.PatientBirthDate = ''
            ds.PatientOrientation = 'P\H'
            ds.PatientSex = 'F'
            ds.PresentationIntentType = 'FOR PROCESSING'
            ds.ReferringPhysicianName = 'Virtual'
            ds.RescaleIntercept = 0
            ds.RescaleSlope = 1
            ds.RescaleType = 'US'
            ds.SeriesNumber = []
            ds.StudyID = ' '

            ds.ViewCodeSequence = ''

            ds.InstitutionName = 'FDA'
            ds.InstitutionalDepartmentName = 'DIDSR'
            ds.SoftwareVersions = 'MC-GPU_1.5b'

            ds.ImageType = 'ORIGINAL\PRIMARY'
            ds.ImageComments = "SA" if len(
                self.lesion_locations[modality]) == 0 else "SP"
            ds.LossyImageCompression = '00'
            ds.ConversionType = 'SYN'

            ds.DetectorType = 'DIRECT'
            ds.DetectorConfiguration = 'AREA'
            ds.DetectorDescription = 'a-Se, {:.2f} micron'.format(
                self.arguments_mcgpu["detector_thickness"] * 10000)  # cm to um
            ds.DetectorActiveShape = 'RECTANGLE'

            # 28 kVp for dense and hetero; 30 kVp for scattered  and fatty
            ds.KVP = '28' if self.arguments_generation["targetFatFrac"] < 0.75 else "30"
            ds.ExposureInmAs = 3.5  # ??
            ds.AnodeTargetMaterial = 'TUNGSTEN'
            ds.FilterType = 'FLAT'
            ds.FilterMaterial = 'RHODIUM'
            # cm to mm
            ds.FilterThicknessMinimum = self.arguments_mcgpu["antiscatter_grid_ratio"][0] * 10

            # cm to mm from source to detector center
            ds.DistanceSourceToDetector = self.arguments_mcgpu["distance_source"] * 10
            # cm to mm from source to the breast support side
            ds.DistanceSourceToPatient = self.arguments_mcgpu["source_position"][2] * 10
            ds.PositionerType = 'MAMMOGRAPHIC'
            ds.DerivationDescription = 'float64 to uint16 bit conversion'

            ds.Columns = data.shape[0]
            ds.Rows = data.shape[1]

            ds.SeriesDescription = modality.upper()
            ds.BodyPartExamined = 'BREAST'
            ds.AcquisitionNumber = count
            ds.InstanceNumber = count

            ds.ImagesInAcquisition = self.recon_size["z"]

            block = ds.private_block(
                0x0099, 'VICTRE/Lesion Information', create=True)

            for idx, lesion in enumerate(self.lesion_locations[modality]):
                if lesion[-1] > 0:
                    block.add_new(idx + 1, 'ST', ' '.join(str(item)
                                                          for item in lesion))
            ds.PixelData = data.tobytes()

            # Set the transfer syntax
            # ds.is_little_endian = True
            # ds.is_implicit_VR = True

            # Set creation date/time
            dt = datetime.datetime.now()
            ds.StudyDate = dt.strftime("%Y%m%d")

            ds.StudyTime = dt.strftime("%H%M")
            ds.ContentDate = dt.strftime('%Y%m%d')
            # long format with micro seconds
            timeStr = dt.strftime('%H%M%S.%f')
            ds.ContentTime = timeStr

            ds.fix_meta_info()

            pydicom.filewriter.dcmwrite(
                "{:s}/{:d}/DICOM_{:s}/{:03d}.dcm".format(
                    self.results_folder, self.seed, modality, count), ds,
                write_like_original=False)

        os.makedirs("{:s}/{:d}/DICOM_{:s}/".format(self.results_folder,
                                                   self.seed,
                                                   modality), exist_ok=True)

        if modality == "dbt":
            pixel_array = np.fromfile("{:s}/{:d}/reconstruction{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                      dtype="float64").reshape(self.recon_size["z"], self.recon_size["x"], self.recon_size["y"])
            pixel_array = np.clip(((2**16 - 1) * pixel_array),
                                  0, 2**16 - 1).astype(np.uint16)
        else:
            pixel_array = np.fromfile("{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                      dtype="float32").reshape(2, self.arguments_mcgpu["image_pixels"][0], self.arguments_mcgpu["image_pixels"][1])
            # pixel_array = ((2**16 - 1) * (pixel_array - np.nanmin(pixel_array)) /
            #                (np.nanmax(pixel_array) - np.nanmin(pixel_array))).astype(np.uint16)
            # pixel_array = (scaling["toUInt16"] * (scaling["offset"] + (pixel_array -
            #                                                            scaling["meanAdditiveNoise"]) * scaling["conversionFactorDM"])).astype(np.uint16)
        pixel_array = np.iinfo(np.uint16).max * (pixel_array - np.nanmin(
            pixel_array)) / (np.nanmax(pixel_array) - np.nanmin(pixel_array))
        bar = progressbar.ProgressBar(
            max_value=pixel_array.shape[0]) if self.verbosity else None
        for s in range(pixel_array.shape[0]):
            bar.update(s) if self.verbosity else None
            save_DICOM_one(np.squeeze(
                pixel_array[s, :, :]).astype(np.uint16), s)
        bar.finish() if self.verbosity else None

    def save_ROIs(self, roi_sizes=None, clean=True, save_folder=None):
        """
            Saves the generated ROIs (absent and present) in RAW and HDF5 formats

            :param roi_sizes: Size of the ROIs for the defined lesion types
            :param clean: If True, the existing ROI folder will be deleted
        """

        if len(self.lesion_locations["dbt"]) == 0:
            cprint("There are no ROIs!", 'red') if self.verbosity else None
            return

        if roi_sizes is not None:
            self.roi_sizes = roi_sizes

        if save_folder is None:
            save_folder = self.results_folder

        # SAVE DBT ROIs
        if clean:
            shutil.rmtree(
                "{:s}/{:d}/ROIs".format(save_folder, self.seed), ignore_errors=True)

        os.makedirs("{:s}/{:d}/ROIs/".format(save_folder,
                                             self.seed), exist_ok=True)

        hf = h5py.File(
            "{:s}/{:d}/ROIs.h5".format(save_folder, self.seed), 'w')

        if os.path.exists("{:s}/{:d}/reconstruction{:d}.raw".format(self.results_folder, self.seed, self.seed)):
            mhd = self._read_mhd(
                "{:s}/{:d}/reconstruction{:d}.mhd".format(self.results_folder, self.seed, self.seed))
            hfdbt = hf.create_group("dbt")
            hfdbt_loc = hf.create_group("dbt_locations")

            pixel_array = np.fromfile("{:s}/{:d}/reconstruction{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                      dtype="float64").reshape(int(mhd["DimSize"][2]), int(mhd["DimSize"][1]), int(mhd["DimSize"][0]))

            for idx, lesion in enumerate(self.lesion_locations["dbt"]):
                lesion_type = np.abs(lesion[3])
                roi = pixel_array[1 + lesion[2] - int(np.ceil(self.roi_sizes[lesion_type][2] / 2)):1 + lesion[2] + int(np.floor(self.roi_sizes[lesion_type][2] / 2)),
                                  lesion[1] - int(np.ceil(self.roi_sizes[lesion_type][1] / 2)): lesion[1] + int(np.floor(self.roi_sizes[lesion_type][1] / 2)),
                                  lesion[0] - int(np.ceil(self.roi_sizes[lesion_type][0] / 2)): lesion[0] + int(np.floor(self.roi_sizes[lesion_type][0] / 2))]
                # with open("./results/{:d}/ROIs/ROI_{:03d}_type{:d}.raw".format(self.seed, idx, lesion_type), 'wb') as f:
                roi.astype(np.dtype('<f8')).tofile(
                    "{:s}/{:d}/ROIs/ROI_DBT_{:02d}_type{:d}.raw".format(save_folder, self.seed, idx, lesion[3]))
                hfdbt.create_dataset("{:d}".format(idx),
                                     data=roi.astype(np.float32), compression="gzip", compression_opts=9)
                hfdbt_loc.create_dataset("{:d}".format(idx), data=lesion)

            hfdbt.create_dataset("lesion_type",
                                 data=np.array(self.lesion_locations["dbt"])[:, 3])

        # SAVE DM ROIs
        pixel_array = np.fromfile("{:s}/{:d}/projection_DM{:d}.raw".format(self.results_folder, self.seed, self.seed),
                                  dtype="float32").reshape(2,
                                                           self.arguments_recon["detector_elements_perpendicular"],
                                                           self.arguments_recon["detector_elements"])

        hfdm = hf.create_group("dm")
        hfdm_loc = hf.create_group("dm_locations")

        for idx, lesion in enumerate(self.lesion_locations["dm"]):
            lesion_type = np.abs(lesion[2])
            roi = pixel_array[0,
                              lesion[0] - int(np.ceil(self.roi_sizes[lesion_type][0] / 2)):lesion[0] + int(np.floor(self.roi_sizes[lesion_type][0] / 2)),
                              lesion[1] - int(np.ceil(self.roi_sizes[lesion_type][1] / 2)):lesion[1] + int(np.floor(self.roi_sizes[lesion_type][1] / 2))]
            # with open("./results/{:d}/ROIs/ROI_{:03d}_type{:d}.raw".format(self.seed, idx, lesion_type), 'wb') as f:
            roi.tofile(
                "{:s}/{:d}/ROIs/ROI_DM_{:02d}_type{:d}.raw".format(save_folder, self.seed, idx, lesion[2]))
            # dm_rois.append(roi)

            hfdm.create_dataset("{:d}".format(idx),
                                data=roi, compression="gzip", compression_opts=9)
            hfdm_loc.create_dataset("{:d}".format(idx), data=lesion)

        hfdm.create_dataset("lesion_type",
                            data=np.array(self.lesion_locations["dm"])[:, 2])

        hf.close()

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] ROIs saved!", 'green', attrs=[
               'bold']) if self.verbosity else None

    def generate_spiculated(self, seed=None, size=None):
        """
            Generates a spiculated mass using the breastMass software

            :param seed: Seed to be used when generating the mass
            :param size: Size of the mass to be used in the breastMass config file
            :returns: None. The result is saved in the `lesions` subfolder
        """

        if seed is not None:
            self.arguments_spiculated["seed"] = seed
        if size is not None:
            self.arguments_spiculated["alpha"] = size

        with open("./Victre/breastMass/configs/spiculated.tpl", "r") as f:
            src = Template(f.read())
            result = src.substitute(self.arguments_spiculated)

        os.makedirs("{:s}/lesions/".format(self.results_folder), exist_ok=True)
        os.makedirs(
            "{:s}/lesions/spiculated/".format(self.results_folder), exist_ok=True)

        with open("{:s}/lesions/spiculated/input_breastMass_{:d}.in".format(self.results_folder, seed), "w") as f:
            f.write(result)

        command = "cd {:s}/lesions/spiculated/ && {:s}/Victre/breastMass/build/breastMass -c input_breastMass_{:d}.in".format(
            self.results_folder,
            os.getcwd(),
            self.arguments_spiculated["seed"])

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Generating mass (seed={:d}, size={:.2f})...".format(
            self.arguments_spiculated["seed"], self.arguments_spiculated["alpha"]), 'cyan') if self.verbosity else None
        os.system(command)

        generated_files = self.get_folder_contents(
            "{:s}/lesions/spiculated/".format(self.results_folder))

        side = None
        for name in generated_files:
            s = re.search(
                ".*\/mass_{:d}_([0-9]*)\.raw".format(self.arguments_spiculated["seed"]), name)
            if s is not None:
                side = int(s[1])

        lesion_raw = np.fromfile("{:s}/lesions/spiculated/mass_{:d}_{:d}.raw".format(self.results_folder, self.arguments_spiculated["seed"], side),
                                 dtype=np.uint8).reshape(side, side, side)

        # clean files
        for name in generated_files:
            if ".raw" in name or ".cfg" in name or ".vti" in name or ".in" in name or "core" in name:
                os.remove(name)

        # save in HDF
        with h5py.File("{:s}/lesions/spiculated/mass_{:d}_size{:.2f}.h5".format(
                self.results_folder,
                self.arguments_spiculated["seed"],
                self.arguments_spiculated["alpha"]), "w") as hf:
            hf.create_dataset("volume", data=lesion_raw,
                              compression="gzip")
            hf.create_dataset("seed", data=self.arguments_spiculated["seed"])
            hf.create_dataset("size", data=self.arguments_spiculated["alpha"])

        self.lesion_file = "{:s}/lesions/spiculated/mass_{:d}_size{:.2f}.h5".format(
            self.results_folder,
            self.arguments_spiculated["seed"],
            self.arguments_spiculated["alpha"])

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Generation finished!", 'green', attrs=[
               'bold']) if self.verbosity else None

    def insert_lesions(self, lesion_type=None, n=-1, lesion_file=None, lesion_size=None, locations=None, roi_sizes=None, save_phantom=True):
        """
            Inserts the specified number of lesions in the phantom.

            :param lesion_type: Constant with the desired lesion type. Check available lesion types and materials in the Constants file.
            :param n: Number of lesions to be added
            :param lesion_file: Path to file including the lesion to be inserted (in HDF5 format). If specified, it will overrite the lesion file specified in the constructor.
            :param lesion_size: If lesion_file is a raw file, lesion_size indicates the size of this file
            :param locations: List of coordinates in the voxel/phantom space where the lesions will be inserted. If not specified, random locations will be generated.
            :param roi_sizes: Size of the region of interest to be calculated to avoid overlapping with other tissues and check out of bounds locations

            :returns: None. A phantom file will be saved inside the results folder with the corresponding raw phantom. Three files will be generated: `pcl_SEED.raw.gz` with the raw data, `pcl_SEED.mhd` with the information about the raw data, and `pcl_SEED.loc` with the voxel coordinates of the lesion centers.

        """
        if self.lesion_file is None and lesion_file is None and save_phantom is True:
            cprint(
                "There is no lesion to insert, just adding lesion locations...", color="cyan") if self.verbosity else None
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
                cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Inserting {:d} non-overlapping lesions...".format(n),
                       'cyan') if self.verbosity else None
            else:
                cprint("Retrieving {:d} lesion locations...".format(
                    n), 'cyan') if self.verbosity else None

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
                    if cand_type == 0:
                        cand_type = 2

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

                loc = {"dm": self.get_coordinates_dm([
                    cand[1],
                    cand[2],
                    cand[0]]),
                    "dbt": self.get_coordinates_dbt([
                        cand[1],
                        cand[2],
                        cand[0]])}

                self.lesion_locations["dm"].append(
                    list(np.round([loc["dm"][0], loc["dm"][1], cand_type]).astype(int)))

                self.lesion_locations["dbt"].append(
                    list(np.round([loc["dbt"][0], loc["dbt"][1], loc["dbt"][2], cand_type]).astype(int)))
        else:
            current_seed = self.seed
            np.random.seed(current_seed)

            if self.candidate_locations is not None:
                Constants.INSERTION_MAX_TRIES = len(self.candidate_locations)
                Constants.INSERTION_MAX_TOTAL_ATTEMPTS = 1000
                np.random.shuffle(self.candidate_locations)
                # current_candidate = 0

            roi_shape = self.roi_sizes[lesion_type]
            c = 0

            max_attempts = Constants.INSERTION_MAX_TOTAL_ATTEMPTS
            while c < n and max_attempts >= 0:
                found = False
                roi = None
                cand = None
                loc = None
                attempts = 0
                bar = progressbar.ProgressBar(
                    max_value=Constants.INSERTION_MAX_TRIES) if self.verbosity else None
                while not found and max_attempts > 0:
                    attempts += 1
                    bar.update(attempts) if self.verbosity else None
                    if attempts == Constants.INSERTION_MAX_TRIES:  # if too many attempts
                        bar.finish() if self.verbosity else None
                        attempts = 0
                        max_attempts -= 1

                        cprint(
                            "Too many attempts at inserting, restarting the insertion! ({:d} remaining)".format(max_attempts), 'red') if self.verbosity else None
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
                        self.lesion_locations["dm"] = self.lesion_locations["dm"][:-c]
                        self.lesion_locations["dbt"] = self.lesion_locations["dbt"][:-c]
                        c = 0

                        if self.candidate_locations is not None:
                            np.random.shuffle(self.candidate_locations)

                        if max_attempts == 0:
                            raise Exceptions.VictreError(
                                "Insertion attempts exceeded")

                        bar = progressbar.ProgressBar(
                            max_value=bar.max_value) if self.verbosity else None
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

                    # check if the locations in DM and DBT are inside the ROI
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

                self.lesion_locations["dm"].append(
                    list(np.round([loc["dm"][0], loc["dm"][1], lesion_type]).astype(int)))

                self.lesion_locations["dbt"].append(
                    list(np.round([loc["dbt"][0], loc["dbt"][1], loc["dbt"][2], lesion_type]).astype(int)))

                c += 1

                bar.finish() if self.verbosity else None

        if lesion is not None:
            np.savetxt("{:s}/{:d}/pcl_{:d}.loc".format(self.results_folder, self.seed, self.seed),
                       np.asarray(self.lesions), fmt="%d")

            # save new phantom file
            if save_phantom:
                cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Saving new phantom...",
                       'cyan') if self.verbosity else None

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

                cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Insertion finished!", 'green', attrs=[
                       'bold']) if self.verbosity else None

    def add_absent_ROIs(self, lesion_type, n=1, locations=None, roi_sizes=None, save_locations=True):
        """
            Adds the specified number of lesion-absent regions of interest.

            :param lesion_type: Constant with the desired lesion type. Check available lesion types and materials in the Constants file.
            :param n: Number of lesions to be added
            :param locations: List of coordinates in the voxel/phantom space where the lesions will be inserted. If not specified, random locations will be generated.
            :param roi_sizes: Size of the region of interest to be calculated to avoid overlapping with other tissues and check out of bounds locations
            :returns: None. A location file will be saved inside the `phantom` folder with the corresponding seed. Negative lesion type means absent ROI.
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
                loc = {"dm": self.get_coordinates_dm([cand[1] + roi_shape[1] / 2,
                                                      cand[2] +
                                                      roi_shape[2] / 2,
                                                      cand[0] + roi_shape[0] / 2]),
                       "dbt": self.get_coordinates_dbt([cand[1] + roi_shape[1] / 2,
                                                        cand[2] +
                                                        roi_shape[2] / 2,
                                                        cand[0] + roi_shape[0] / 2])}

                self.lesion_locations["dm"].append(
                    list(np.round([loc["dm"][0], loc["dm"][1], -lesion_type]).astype(int)))

                self.lesion_locations["dbt"].append(
                    list(np.round([loc["dbt"][0], loc["dbt"][1], loc["dbt"][2], -lesion_type]).astype(int)))
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

                    if np.any(np.array(loc["dm"][:2]) < np.array(roi_shape[:2])) or \
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
        if save_locations:
            np.savetxt("{:s}/{:d}/pcl_{:d}.loc".format(self.results_folder, self.seed, self.seed),
                       self.lesions, fmt="%d")

    def generate_phantom(self):
        """
            Runs breast phantom generation.

            :returns: None. A phantom file will be saved inside the results folder with the corresponding raw phantom. Two files will be generated: `p_SEED.raw.gz` with the raw data, and `p_SEED.mhd` with the information about the raw data.
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

        command = "cd {:s} && ./Victre/generation/build/breastPhantomMain -c {:s}".format(
            os.getcwd(),
            full_path
        )

        if self.ips["cpu"] == "localhost":
            ssh_command = command
        else:
            ssh_command = "ssh -Y {:s} \"{:s}\"".format(
                self.ips["cpu"], command)

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Starting phantom generation (seed = {:d}), this will take some time...".format(
            self.seed), 'cyan') if self.verbosity else None

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

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Generation finished!", 'green', attrs=[
               'bold']) if self.verbosity else None
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

        self.candidate_locations = self._mm_to_voxels(np.loadtxt(
            "{:s}/{:d}/p_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist())

    def compress_phantom(self, thickness=None):
        """
            Runs the FEBio compression.

            :param thickness: Specifies the objective thickness for the phantom to be compressed (in cm)
            :returns: None. A phantom file will be saved inside the results folder with the corresponding raw phantom. Two files will be generated: `pc_SEED.raw.gz` with the raw data, and `pc_SEED.mhd` with the information about the raw data.
        """
        if thickness is None:
            # thickness = int(self.arguments_generation["compressionThickness"])
            interp = interpolate.interp1d(
                Constants.DENSITY_RANGES["breastHeight"],
                Constants.DENSITY_RANGES["compressionThickness"],
                fill_value="extrapolate")
            thickness = np.round(float(interp(
                self.arguments_mcgpu["number_voxels"][2] * self.arguments_mcgpu["voxel_size"][2] * 10)), 2)

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

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Starting phantom compression, this will take some time...",
               'cyan') if self.verbosity else None

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

                if "fault" in output:
                    completed = 0
                    break

        if completed == 0 or not os.path.exists("{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, self.seed, self.seed)):
            cprint("\nError while compressing, check the output_compression.out file in the results folder",
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Compression error")

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "] Compression finished!", 'green', attrs=[
               'bold']) if self.verbosity else None
        self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}.raw.gz".format(
            self.results_folder, self.seed, self.seed)

        self.mhd = self._read_mhd(
            "{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
        self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]

        self.arguments_recon["voxels_x"] = self.arguments_mcgpu["number_voxels"][1]
        self.arguments_recon["voxels_y"] = self.arguments_mcgpu["number_voxels"][0]
        self.arguments_recon["voxels_z"] = self.arguments_mcgpu["number_voxels"][2]

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

        self.candidate_locations = self._mm_to_voxels(np.loadtxt(
            "{:s}/{:d}/pc_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist())

    def crop(self, size=None):
        """
            Runs breast phantom cropping.

            :returns: None. A phantom file will be saved inside the results folder with the corresponding raw phantom. Two files will be generated: `pc_SEED_crop.raw.gz` with the raw data, and `pc_SEED_crop.mhd` with the information about the raw data.
        """

        cprint("[" + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") +
               "] Cropping phantom...", 'cyan') if self.verbosity else None

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

        self.candidate_locations = np.loadtxt(
            "{:s}/{:d}/pc_{:d}.loc".format(self.results_folder, self.seed, self.seed), delimiter=',').tolist()

        if self.candidate_locations is not None:
            for cand in self.candidate_locations:
                cand[0] = ((cand[0] - prevOffset[0]) / self.mhd["ElementSpacing"][0] -
                           crop["from"][0]) * self.mhd["ElementSpacing"][0] + self.mhd["Offset"][0]
                cand[1] = ((cand[1] - prevOffset[1]) / self.mhd["ElementSpacing"][1] -
                           crop["from"][1]) * self.mhd["ElementSpacing"][1] + self.mhd["Offset"][1]
                cand[2] = ((cand[2] - prevOffset[2]) / self.mhd["ElementSpacing"][2] -
                           crop["from"][2]) * self.mhd["ElementSpacing"][2] + self.mhd["Offset"][2]
            #saving in mm
            np.savetxt("{:s}/{:d}/pc_{:d}_crop.loc".format(self.results_folder,
                                                           self.seed,
                                                           self.seed),
                       self.candidate_locations,
                       delimiter=',')

            self.candidate_locations = self._mm_to_voxels(
                self.candidate_locations)

    def get_dm_segmentation(self, roi=None, selected_materials=[]):
        with gzip.open(self.arguments_mcgpu["phantom_file"], 'rb') as f:
            phantom = f.read()
        phantom = np.fromstring(phantom, dtype=np.uint8).reshape(
            self.arguments_mcgpu["number_voxels"][2],
            self.arguments_mcgpu["number_voxels"][1],
            self.arguments_mcgpu["number_voxels"][0])

        if roi is None:
            roi = [[0, 0], self.arguments_mcgpu["image_pixels"][::-1]]

        mask = [[[] for _ in range(roi[1][1] - roi[0][1])]
                for _ in range(roi[1][0] - roi[0][0])]

        for row in progressbar.progressbar(range(roi[0][0], roi[1][0])):
            for col in range(roi[0][1], roi[1][1]):

                vx_location2 = self.reverse_dm_coordinates([row, col])
                for loc in vx_location2:
                    mat = phantom[loc[2],
                                  loc[1],
                                  loc[0]]
                    if mat != 0 and (len(selected_materials) == 0 or mat in selected_materials):
                        mask[row - roi[0][0]][col - roi[0][1]].append(mat)

        return mask

    def get_DBT_segmentation(self):
        with gzip.open(self.arguments_mcgpu["phantom_file"], 'rb') as f:
            phantom = f.read()
        phantom = np.fromstring(phantom, dtype=np.uint8).reshape(
            self.arguments_mcgpu["number_voxels"][2],
            self.arguments_mcgpu["number_voxels"][1],
            self.arguments_mcgpu["number_voxels"][0])

        recon_mhd = self._read_mhd("{:s}/{:d}/reconstruction{:d}.mhd".format(self.results_folder,
                                                                             self.seed,
                                                                             self.seed))

        mask = np.zeros(
            [recon_mhd["DimSize"][2], recon_mhd["DimSize"][1], recon_mhd["DimSize"][0]], dtype=np.uint8)

        for x in progressbar.progressbar(range(mask.shape[0])):
            for y in range(mask.shape[1]):
                for z in range(mask.shape[2]):
                    try:
                        vx_location = [int(x)
                                       for x in self.reverse_dbt_coordinates([z, y, x])]
                        mask[x, y, z] = phantom[vx_location[2],
                                                vx_location[1],
                                                vx_location[0]]
                    except:
                        pass
        return mask

    @staticmethod
    def get_folder_contents(folder):
        """
            Gets a list of files in the given folder

            :param folder: Path to the folder to be processed
            :returns: List with files inside the given folder
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

    def _mm_to_voxels(self, locations):
        if locations is not None:
            for idx, cand in enumerate(locations):
                locations[idx] = [int(np.round((cand[0] - self.mhd["Offset"][0]) /
                                               self.mhd["ElementSpacing"][0])),
                                  int(np.round((cand[2] - self.mhd["Offset"][2]) /
                                               self.mhd["ElementSpacing"][2])),
                                  int(np.round((cand[1] - self.mhd["Offset"][1]) /
                                               self.mhd["ElementSpacing"][1]))]
        return locations
