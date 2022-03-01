import glob
import os
from pathlib import Path

import numpy as np
import skimage.io
import yaml


def read_parameters(parameter_file):
    """Reads in default parameters and replaces user defined parameters."""
    current_path = Path(os.path.dirname(os.path.realpath(__file__)))

    param_base_file = Path(current_path).joinpath("base", "parameters.yml")

    with open(param_base_file, 'r') as yml_f:
        parameters = yaml.safe_load(yml_f)

    with open(parameter_file) as file:
        parameters_local = yaml.safe_load(file)

    # overwrite global parameters with local setting
    for key in parameters_local:
        parameters[key] = parameters_local[key]

    return parameters


def read_image(filename):
    """Reads an image and reshapes to channel last."""
    img_ = skimage.io.imread(filename)

    if len(img_.shape) <= 2:
        img_ = np.array([img_, img_])

    if img_.shape[0] < min(img_.shape[1], img_.shape[2]):
        print("Warning: channel is on the first dimension of the image.")
        img = img_.reshape(img_.shape[1], img_.shape[2], img_.shape[0])
    else:
        img = img_

    return img


def write_dict_to_yml(yml_file, d):
    """Writes a dictionary to a file in yml format."""
    yml_file = Path(yml_file)
    p = Path(yml_file.parent)
    p.mkdir(parents=True, exist_ok=True)

    with open(yml_file, 'w+') as yml_f:
        yml_f.write(yaml.dump(d, Dumper=yaml.Dumper))

    return True


def create_path_recursively(path):
    """Creates a path. Creates missing parent folders."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)

    return True


def get_tif_list(path):
    path = str(path)

    return glob.glob(path + "*.tif")