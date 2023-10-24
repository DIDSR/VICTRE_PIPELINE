## M-SYNTH GENERATION EXAMPLE

The file `msynth_VICTRE.py` contains the projection script for a given phantom an certain input parameters:

* `--results`: results folder where the projection will be stored
* `--seed`: phantom seed to be used
* `--density`: breast density to get default parameters (default: dense)
* `--dose`: number of histories for the MCGPU dose (optional)
* `--simulate_lesions`: flag to add 1 lesion to the phantom (optional)
* `--lesion_file`: path to the lesion file to be inserted 
* `--lesion_density`: factor to multiply by the default lesion density (optional)
* `--phantom_file`: path to the phantom file to be projected
* `--flatfield_file`: path to the flatfield file to be used (optional)
* `--lesion_only`: flag for a lesion only/segmentation projection (optional)

The file `msynth_project.sh` provides an example to run a single case. Change the file paths accordingly and run the file with `sh msynth_project.sh`. The results folder will be created with a subfolder with the `seed` you used. Inside you will find the `projection_DM` files with the in silico digital mammography. The corresponding `.raw` and `.mhd` files can be opened with image software such as `ImageJ`.

More information about the M-SYNTH dataset:
> https://github.com/DIDSR/msynth-release/

## Citation:

```
@article{sizikova2023knowledge,
  title={Knowledge-based in silico models and dataset for the regulatory evaluation of mammography AI for a range of breast characteristics, lesion conspicuities and doses},
  author={Sizikova, Elena and Saharkhiz, Niloufar and Sharma, Diksha and Lago, Miguel and Sahiner, Berkman and Delfino, Jana G. and Badano, Aldo},
  journal={Advances in Neural Information Processing Systems},
  volume={},
  pages={16764--16778},
  year={2023}
}
```
