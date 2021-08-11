import os
import yaml


def read_parameters(parameter_file):

    parameters = dict()

    if os.path.exists("../feature_extraction/base/parameters.yml"):
        with open("../feature_extraction/base/parameters.yml") as file:
            parameters = yaml.load(file, Loader=yaml.FullLoader)


    #with open("base/parameters.yml") as file:
    #    parameters = yaml.load(file, Loader=yaml.FullLoader)
    with open(parameter_file) as file:
        parameters_local = yaml.load(file, Loader=yaml.FullLoader)

    # overwrite global parameters with local setting
    for key in parameters_local:
        parameters[key] = parameters_local[key]
    
    return parameters


def get_cellpose_segmentation(filename):

    return 0


def get_features_from_cellpose_seg():

    return 0

