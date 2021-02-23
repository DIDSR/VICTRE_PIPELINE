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
import matplotlib.pyplot as plt
import cv2
import cc3d
import progressbar
import h5py
import subprocess
import asyncio
from string import Template
import random
import time
from . import Constants, Exceptions
import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
import tempfile
import datetime
from pydicom.encaps import encapsulate
import re
import gzip


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
                 arguments_mcgpu=dict(),
                 arguments_recon=dict(),
                 arguments_generation=dict(),
                 flatfield_DBT=None,
                 flatfield_DM=None):
        """!
        Object constructor for the Victre pipeline class

        @param ips Dictionary with two IP addresses to run the pipeline: "gpu" for the projection process. "cpu" for the reconstruction.
        @param seed Random seed used to generate or read the phantom
        @param results_folder Path to folder to be used when saving the results
        @param phantom_file Path to file containing the phantom to be loaded
        @param spectrum_file Path to file containing the spectrum used to project in MCGPU
        @param lesion_file Path to file containing the lesion to be inserted (in HDF5 format)
        @param materials Dictionary including the materials to be used during projection
        @param arguments_mcgpu Arguments to be overridden for the projection in MCGPU
        @param arguments_recon Arguments to be overridden for the reconstruction algorithm
        @param flatfield_file Path to the flatfield file for projection
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
        self.roi_sizes=roi_sizes

        random.seed(self.seed)

        self.arguments_mcgpu = dict(
            number_histories=1.02e10,
            random_seed=31415990,
            selected_gpu=0,
            number_gpus=1,
            gpu_threads=128,
            histories_per_thread=5000,
            spectrum_file=spectrum_file,
            source_position=[0.00001, 4.825, 63.0],
            source_direction=[0.0, 0.0, -1.0],
            fam_beam_aperture=[15.0, 7.4686667],
            euler_angles=[90.0, -90.0, 180.0],
            focal_spot=0.0300,
            angular_blur=0.18,
            collimate_beam="YES",
            output_file="{:s}/{:d}/projection".format(
                self.results_folder, self.seed),
            image_pixels=[3000, 1500],
            image_size=[25.50, 12.75],
            distance_source=65.00,
            image_offset=[0, 0],
            detector_thickness=0.02,
            mean_free_path=0.004027,
            k_edge_energy=[12658.0, 11223.0, 0.596, 0.00593],
            detector_gain=[50.0, 0.99],
            additive_noise=5200.0,
            cover_thickness=[0.10, 1.9616],
            antiscatter_grid_ratio=[5.0, 31.0, 0.0065],
            antiscatter_strips=[0.00089945, 1.9616],
            antiscatter_grid_lines=0,
            number_projections=25,
            rotation_axis_distance=60.0,
            projections_angle=2.083333333333,
            angular_rotation_first=-25.0,
            rotation_axis=[1.0, 0.0, 0.0],
            axis_translation=0,
            detector_fixed="YES",
            simulate_both="YES",
            tally_material_dose="YES",
            tally_voxel_dose="NO",
            output_dose_filename="mc-gpu_dose.dat",
            roi_voxel_dose_x=[1, 751],
            roi_voxel_dose_y=[1, 1301],
            roi_voxel_dose_z=[250, 250],
            phantom_file=phantom_file,
            voxel_geometry_offset=[0, 0, 0],
            number_voxels=[810, 1920, 745],
            voxel_size=[0.005, 0.005, 0.005],
            low_resolution_voxel_size=[1, 1, 1]
        )
        if phantom_file is None:
            if os.path.exists("{:s}/{:d}/pcl_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found phantom with lesions information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pcl_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                self.lesions = np.loadtxt(
                    "{:s}/{:d}/pcl_{:d}.loc".format(self.results_folder, self.seed, self.seed))

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pcl_{:d}.raw.gz".format(
                    self.results_folder, seed, seed)
            elif os.path.exists("{:s}/{:d}/pc_{:d}_crop.mhd".format(self.results_folder, seed, seed)):
                cprint("Found cropped phantom information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pc_{:d}_crop.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}_crop.raw.gz".format(
                    self.results_folder, seed, seed)
            elif os.path.exists("{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found compressed phantom information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}.raw.gz".format(
                    self.results_folder, seed, seed)

            elif os.path.exists("{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, seed, seed)):
                cprint("Found phantom generation information!", 'cyan')
                self.mhd = self._read_mhd(
                    "{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
                self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
                self.arguments_mcgpu["voxel_size"] = [
                    x / 10 for x in self.mhd["ElementSpacing"]]

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
            reconstruction_file="{:s}/{:d}/reconstruction.raw".format(
                self.results_folder,
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

        if self.flatfield_DM is None and os.path.exists("{:s}/{:d}/flatfield_DM.raw".format(
                self.results_folder,
                self.seed)):
            self.flatfield_DM = "{:s}/{:d}/flatfield_DM.raw".format(
                self.results_folder,
                self.seed)

        self.arguments_recon.update(arguments_recon)

        self.arguments_generation = Constants.VICTRE_DENSE

        self.arguments_generation["outputDir"] = os.path.abspath("{:s}/{:d}/".format(
            self.results_folder, self.seed))
        self.arguments_generation["seed"] = self.seed
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
                locations = list(np.loadtxt("{:s}/{:s}.loc".format(path, filename)))
                self.insert_lesions(locations=locations,
                                    save_phantom=False)

        # self.arguments_mcgpu["number_voxels"]

        os.makedirs("{:s}".format(self.results_folder), exist_ok=True)
        os.makedirs("{:s}/{:d}".format(self.results_folder,
                                       self.seed), exist_ok=True)

    def project(self, clean=True, do_flatfield=0):
        """!
            Method that runs MCGPU to project the phantom.

            @param clean If True, it will delete the contents of the output folder before projecting.
        """

        if do_flatfield > 0:
            filename = "flatfield"
            empty_phantom = np.zeros(
                self.arguments_mcgpu["number_voxels"], np.uint8)
            with gzip.open("{:s}/{:d}/empty_phantom.raw.gz".format(
                    self.results_folder, self.seed), "wb") as gz:
                gz.write(empty_phantom)

            prev_flatfield_DBT, prev_flatfield_DM = None, None

            if os.path.exists("{:s}/{:d}/{:s}_DM.raw".format(self.results_folder, self.seed, filename)):
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

        phantom_config = "{:s}/{:d}/input_{:s}.in".format(
            self.results_folder, self.seed, filename)

        with open("./Victre/projection/configs/template_mcgpu.tpl", "r") as f:
            src = Template(f.read())
            template_arguments = self.arguments_mcgpu.copy()
            if do_flatfield > 0:
                template_arguments["phantom_file"] = "{:s}/{:d}/empty_phantom.raw.gz".format(
                    self.results_folder,
                    self.seed)
            template_arguments["output_file"] = "{:s}/{:d}/{:s}".format(
                self.results_folder, self.seed, filename)
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

        command = "cd {:s} && time mpirun -v -n {:d} ./Victre/projection/MC-GPU_v1.5b.x {:s}".format(
            os.getcwd(),
            self.arguments_mcgpu["number_gpus"],
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
                                   stderr=subprocess.PIPE)

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
            output = process.stderr.readline().decode("utf-8")
            with open("{:s}/{:d}/output_{:s}.out".format(self.results_folder, self.seed, filename), "ab+") as f:
                f.write(output.encode('utf-8'))
            cprint("\nError while projecting, check the output file in the results folder",
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
                "{:s}/{:d}/{:s}_DM.raw".format(self.results_folder, self.seed, filename))

        os.rename("{:s}/{:d}/{:s}_0000.raw".format(self.results_folder, self.seed, filename),
                  "{:s}/{:d}/{:s}_DM.raw".format(self.results_folder, self.seed, filename))

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
                curr_flatfield_DM = np.fromfile("{:s}/{:d}/flatfield_DM.raw".format(self.results_folder, self.seed),
                                                dtype="float32").reshape(2,
                                                                         self.arguments_recon["detector_elements_perpendicular"],
                                                                         self.arguments_recon["detector_elements"])

                prev_flatfield_DM += curr_flatfield_DM / do_flatfield

                prev_flatfield_DM.tofile(
                    "{:s}/{:d}/flatfield_DM.raw".format(self.results_folder, self.seed))

            if prev_flatfield_DBT is not None:
                curr_flatfield_DBT = np.fromfile("{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                    self.results_folder,
                    self.seed,
                    'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                    self.arguments_mcgpu["number_projections"]),
                    dtype="float32").reshape(self.arguments_mcgpu["number_projections"],
                                             self.arguments_mcgpu["image_pixels"][0],
                                             self.arguments_mcgpu["image_pixels"][1])

                prev_flatfield_DBT += curr_flatfield_DBT / do_flatfield

                prev_flatfield_DBT.tofile("{:s}/{:d}/flatfield_{:s}pixels_{:d}proj.raw".format(
                    self.results_folder,
                    self.seed,
                    'x'.join(map(str, self.arguments_mcgpu["image_pixels"])),
                    self.arguments_mcgpu["number_projections"]))

        elif self.arguments_recon["flatfield_file"] is None:
            with contextlib.suppress(FileNotFoundError):
                os.remove(
                    "{:s}/{:d}/flatfield_DM.raw".format(self.results_folder, self.seed).format(
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

            self.flatfield_DM = "{:s}/{:d}/flatfield_DM.raw".format(
                self.results_folder, self.seed)

        if do_flatfield == 0:
            # normalize with flatfield
            curr_flatfield_DM = np.fromfile(self.flatfield_DM,
                                            dtype="float32").reshape(2,
                                                                     self.arguments_recon["detector_elements_perpendicular"],
                                                                     self.arguments_recon["detector_elements"])
            projection_DM = np.fromfile("{:s}/{:d}/projection_DM.raw".format(self.results_folder, self.seed),
                                        dtype="float32").reshape(2,
                                                                 self.arguments_recon["detector_elements_perpendicular"],
                                                                 self.arguments_recon["detector_elements"])

            projection_DM = np.divide(curr_flatfield_DM, projection_DM)

            projection_DM.tofile(
                "{:s}/{:d}/projection_DM.raw".format(self.results_folder, self.seed))

    def reconstruct(self):
        """!
            Method that runs the reconstruction code for the DBT volume
        """

        # %% RECONSTRUCTION
        with open("./Victre/reconstruction/configs/parameters.tpl", "r") as f:
            src = Template(f.read())
            template_arguments = self.arguments_recon.copy()
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
                                   stderr=subprocess.PIPE)

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
            output = process.stderr.readline().decode("utf-8")
            with open("{:s}/{:d}/output_recon.out".format(self.results_folder, self.seed), "ab+") as f:
                f.write(output.encode('utf-8'))
            cprint("\nError while reconstructing, check the input file",
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Reconstruction error")

        bar.finish()
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

        pixel_array = np.fromfile("{:s}/{:d}/reconstruction.raw".format(self.results_folder, self.seed),
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
            self.roi_sizes=roi_sizes

        # SAVE DBT ROIs
        if clean:
            shutil.rmtree(
                "{:s}/{:d}/ROIs".format(self.results_folder, self.seed), ignore_errors=True)

        os.makedirs("{:s}/{:d}/ROIs/".format(self.results_folder,
                                             self.seed), exist_ok=True)

        hf = h5py.File(
            "{:s}/{:d}/ROIs.h5".format(self.results_folder, self.seed), 'w')
        hfdbt = hf.create_group("dbt")

        pixel_array = np.fromfile("{:s}/{:d}/reconstruction.raw".format(self.results_folder, self.seed),
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
        pixel_array = np.fromfile("{:s}/{:d}/projection_DM.raw".format(self.results_folder, self.seed),
                                  dtype="float32").reshape(2,
                                                           self.arguments_recon["detector_elements_perpendicular"],
                                                           self.arguments_recon["detector_elements"])

        # invert LUT
        pixel_array = np.max(pixel_array) - pixel_array

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

    def get_folder_contents(self, folder):
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

    def insert_lesions(self, lesion_type=None, n=1, lesion_file=None, lesion_size=None, locations=None, roi_sizes=None, save_phantom=True):
        """!
            Inserts the specified number of lesions in the phantom.

            @param lesion_type Constant with the desired lesion type. Check available lesion types and materials in the Constants file.
            @param n Number of lesions to be added
            @param lesion_file Path to file including the lesion to be inserted (in HDF5 format). If specified, it will overrite the lesion file specified in the constructor.
            @param lesion_size If lesion_file is a raw file, lesion_size indicates the size of this file
            @param locations List of coordinates in the voxel/phantom space where the lesions will be inserted. If not specified, random locations will be generated.
            @param roi_sizes Size of the region of interest to be calculated to avoid overlapping with other tissues and check out of bounds locations
            @return None. A location file will be saved inside the `phantom` folder with the corresponding seed. Positive lesion type means absent ROI.
        """
        if self.lesion_file is None and lesion_file is None:
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
            if locations is not None:
                n = len(locations)
            cprint(
                "Inserting {:d} non-overlapping lesions...".format(n), 'cyan')

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
                loc = {"dm": self.get_coordinates_dm([
                    cand[1] + lesion_shape[1] / 2,
                    cand[2] + lesion_shape[2] / 2,
                    cand[0] + lesion_shape[0] / 2]),
                    "dbt": self.get_coordinates_dbt([
                        cand[1] + lesion_shape[1] / 2,
                        cand[2] + lesion_shape[2] / 2,
                        cand[0] + lesion_shape[0] / 2])}

                self.lesion_locations["dm"].append(
                    list(np.round([loc["dm"][0], loc["dm"][1], cand_type]).astype(int)))

                self.lesion_locations["dbt"].append(
                    list(np.round([loc["dbt"][0], loc["dbt"][1], loc["dbt"][2], cand_type]).astype(int)))
        else:
            roi_shape = self.roi_sizes[lesion_type]
            c = 0
            current_seed = self.seed
            random.seed(current_seed)
            max_attempts = 10
            while c < n and max_attempts >= 0:
                found = False
                roi = None
                cand = None
                loc = None
                attempts = 0
                bar = progressbar.ProgressBar(max_value=200)
                while not found and max_attempts >= 0:
                    attempts += 1
                    bar.update(attempts)
                    if attempts == bar.max_value:  # if too many attempts
                        attempts = 1
                        max_attempts -= 1
                        bar.finish()
                        bar = progressbar.ProgressBar(max_value=bar.max_value)
                        cprint(
                            "\nToo many attempts at inserting, restarting the insertion! ({:d} remaining)".format(max_attempts), 'red')
                        with gzip.open(self.arguments_mcgpu["phantom_file"], 'rb') as f:
                            phantom = f.read()

                        phantom = np.fromstring(phantom, dtype=np.uint8).reshape(
                            self.arguments_mcgpu["number_voxels"][2],
                            self.arguments_mcgpu["number_voxels"][1],
                            self.arguments_mcgpu["number_voxels"][0])
                        current_seed += 1
                        random.seed(current_seed)  # try with a different seed

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
                    if not (np.any([np.any(roi == x)
                                    for x in np.append(Constants.FORBIDDEN_OVERLAP,
                                                       Constants.LESION_MATERIALS)])):
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

                bar.finish()

        if lesion is not None:
            np.savetxt("{:s}/{:d}/pcl_{:d}.loc".format(self.results_folder, self.seed, self.seed),
                       self.lesions, fmt="%d")

            # with h5py.File("phantom/pcl_{:d}_crop.h5".format(self.seed), "w") as hf:
            #     hf.create_dataset("phantom", data=phantom.astype(
            #         np.uint8), compression="gzip")

            # save new phantom file
            if save_phantom:
                cprint("Saving new phantom...", 'cyan')

                # UNCOMMENT TO USE SYSTEM'S GZIP
                # phantom.tofile(
                #     "{:s}/{:d}/phantom.raw".format(self.results_folder, self.seed))
                # os.system("gzip -f {:s}/{:d}/phantom.raw".format(
                #     self.results_folder, self.seed))

                # We save the phantom in gzip to reduce needed disk space
                with gzip.open("{:s}/{:d}/pcl_{:d}.raw.gz".format(self.results_folder, self.seed, self.seed), "wb") as gz:
                    gz.write(phantom)

                self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pcl_{:d}.raw.gz".format(
                    self.results_folder, self.seed, self.seed)

                with open("{:s}/{:d}/pcl_{:d}.mhd".format(self.results_folder, self.seed, self.seed), "w") as f:
                    src = Template(Constants.MHD_FILE)
                    template_arguments = self.mhd.copy()
                    template_arguments["ElementDataFile"] = "pcl_{:d}.raw"
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
                   np.asarray(self.lesions), fmt="%d")

    def generate_phantom(self):
        generation_config = "{:s}/{:d}/input_generation.in".format(
            self.results_folder, self.seed)

        with open("./Victre/generation/configs/template_generation.tpl", "r") as f:
            src = Template(f.read())
            template_arguments = self.arguments_generation.copy()
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

        cprint("Starting phantom generation, this will take some time...", 'cyan')

        completed = 0

        process = subprocess.Popen(ssh_command, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        with open("{:s}/{:d}/output_generation.out".format(self.results_folder, self.seed), "wb") as f:
            while True:
                output = process.stdout.readline().decode("utf-8")
                if output == "" and process.poll() is not None:
                    break

                completed += 1
                f.write(output.encode('utf-8'))
                f.flush()

        if completed == 0:
            output = process.stderr.readline().decode("utf-8")
            with open("{:s}/{:d}/output_generation.out".format(self.results_folder, self.seed), "ab+") as f:
                f.write(output.encode('utf-8'))
            cprint("\nError while generating, check the output file in the results folder",
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Generation error")

        cprint("Generation finished!", 'green', attrs=['bold'])
        self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/p_{:d}.raw.gz".format(
            self.results_folder, self.seed, self.seed)

        self.mhd = self._read_mhd(
            "{:s}/{:d}/p_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
        self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
        self.arguments_mcgpu["voxel_size"] = self.mhd["ElementSpacing"]
        self.arguments_recon["voxel_size"] = self.arguments_mcgpu["voxel_size"][0]

    def compress_phantom(self, thickness):
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
                                   stderr=subprocess.PIPE)

        with open("{:s}/{:d}/output_compression.out".format(self.results_folder, self.seed), "wb") as f:
            while True:
                output = process.stdout.readline().decode("utf-8")
                if output == "" and process.poll() is not None:
                    break

                completed += 1

                f.write(output.encode('utf-8'))
                f.flush()

        if completed == 0:
            output = process.stderr.readline().decode("utf-8")
            with open("{:s}/{:d}/output_compression.out".format(self.results_folder, self.seed), "ab+") as f:
                f.write(output.encode('utf-8'))
            cprint("\nError while compressing, check the output file in the results folder",
                   'red', attrs=['bold'])
            raise Exceptions.VictreError("Generation error")

        cprint("Compression finished!", 'green', attrs=['bold'])
        self.arguments_mcgpu["phantom_file"] = "{:s}/{:d}/pc_{:d}.raw.gz".format(
            self.results_folder, self.seed, self.seed)

        self.mhd = self._read_mhd(
            "{:s}/{:d}/pc_{:d}.mhd".format(self.results_folder, self.seed, self.seed))
        self.arguments_mcgpu["number_voxels"] = self.mhd["DimSize"]
        self.arguments_mcgpu["voxel_size"] = self.mhd["ElementSpacing"]
        self.arguments_recon["voxel_size"] = self.arguments_mcgpu["voxel_size"][0]

        self.arguments_recon["voxels_x"] = self.arguments_mcgpu["number_voxels"][1]
        self.arguments_recon["voxels_y"] = self.arguments_mcgpu["number_voxels"][0]
        self.arguments_recon["voxels_z"] = self.arguments_mcgpu["number_voxels"][2]

    def crop(self, size=None):
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

        with open("{:s}/{:d}/pc_{:d}_crop.mhd".format(self.results_folder, self.seed, self.seed), "w") as f:
            src = Template(Constants.MHD_FILE)
            template_arguments = self.mhd.copy()
            for key in template_arguments.keys():
                if type(template_arguments[key]) is list:
                    template_arguments[key] = ' '.join(
                        map(str, template_arguments[key]))
            result = src.substitute(template_arguments)
            f.write(result)

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
