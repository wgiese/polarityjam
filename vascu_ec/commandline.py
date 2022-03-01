import json
import os
from pathlib import Path

import numpy as np

from vascu_ec.feature_extraction import get_image_for_segmentation, get_features_from_cellpose_seg_multi_channel
from vascu_ec.logging import get_logger
from vascu_ec.utils.io import read_parameters, read_image, create_path_recursively, get_tif_list
from vascu_ec.utils.plot import plot_seg_channels
from vascu_ec.utils.seg import load_or_get_cellpose_segmentation


def run(args, print_args=True):
    # read args
    param_file = args.param
    filepath = args.in_file
    filename = args.filename_prefix
    output_path = args.out_path

    # print info
    if print_args:
        get_logger().info(
            "Arguments provided: %s" % json.dumps(
                {"param": param_file, "in_file": filepath, "out_path": output_path, "filename": filename},
                sort_keys=True, indent=4
            )
        )

    # process args
    if not filename:
        filename, _ = os.path.splitext(os.path.basename(filepath))

    create_path_recursively(output_path)

    parameters = read_parameters(param_file)
    get_logger().info("Parameters: %s" % json.dumps(parameters, sort_keys=True, indent=4))

    # start routine
    img = read_image(filepath)
    img_seg = get_image_for_segmentation(parameters, img)

    plot_seg_channels(img_seg, output_path, filename)

    cellpose_mask = load_or_get_cellpose_segmentation(parameters, img_seg, filepath)

    get_logger().info("Number of cell labels: %s" % np.max(cellpose_mask))

    properties_df = get_features_from_cellpose_seg_multi_channel(parameters, img, cellpose_mask, filename, output_path)

    get_logger().info("Head of created dataset: \n %s" % properties_df.head())

    properties_df.to_csv(output_path + filename.split(".")[0] + ".csv")


def run_stack(args):
    # read args
    param_file = args.param
    inpath = args.in_path
    output_path = args.out_path

    # print info
    get_logger().info(
        "Arguments provided: %s" % json.dumps(
            {"param": param_file, "in_path": inpath, "out_path": output_path},
            sort_keys=True, indent=4
        )
    )

    file_list = get_tif_list(inpath)

    for ind, filepath in enumerate(file_list):
        filename = Path(filepath).parts[-1]
        stem, file_extension = os.path.splitext(filepath)

        if not ((file_extension != ".tif") or (file_extension != ".tiff")):
            continue

        get_logger().info("Processing fil with: file stem  %s and file extension: %s" % (stem, file_extension))

        args.infile = filename

        # call
        run(args, print_args=False)


def run_key(args):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    mpl.rcParams['figure.dpi'] = 300
    import sys
    import os
    import glob

    sys.path.append("../vascu_ec/")
    import functions

    if len(sys.argv) > 1:
        parameter_file = sys.argv[1]
    else:
        parameter_file = "base/parameters.yml"

    print("Read parameters from: %s" % parameter_file)

    parameters = functions.read_parameters(parameter_file)
    key_file = pd.read_csv(parameters['key_file'])
    output_path_base = parameters['output_path']
    summary_df = pd.DataFrame()
    summary_ind = 0

    for index, row in key_file.iterrows():

        input_path = parameters['input_path'] + str(row['folder_name'])
        # subfolder = "/results_%sdyn_%sh_%s/" % (row["shear_stress"], row['flow_exposure_time'], row['treatment'])
        subfolder = str(row["short_name"]) + "/"
        output_path = output_path_base + subfolder
        if not os.path.exists(output_path):
            os.mkdir(output_path)

        print("input path: %s" % input_path)
        print("output path: %s" % output_path)

        file_list = glob.glob(input_path + "*.tif")
        merged_properties_df = pd.DataFrame()

        for filepath in file_list:
            print(filepath)

        for ind, filepath in enumerate(file_list):

            # output_path = parameters['output_path']

            print("processing", filepath.split("/")[-1].split(".")[0])
            filename = filepath.split("/")[-1]

            img = functions.read_image(parameters, filepath)
            img_seg = functions.get_image_for_segmentation(parameters, img)
            if parameters["channel_nucleus"] >= 0:
                fig, ax = plt.subplots(1, 2)
                ax[0].imshow(img_seg[0, :, :])
                ax[1].imshow(img_seg[1, :, :])
                plt.savefig(output_path + filename + "-seg.png")
            else:
                fig, ax = plt.subplots()
                ax.imshow(img_seg[:, :])
                plt.savefig(output_path + filename + "-seg.png")

            cellpose_mask = functions.get_cellpose_segmentation(parameters, img_seg)
            print("Number of cell labels: ")
            print(len(np.unique(cellpose_mask)))
            # nuclei_mask = functions.get_nuclei_mask(parameters, img, cellpose_mask)
            # golgi_mask = functions.get_golgi_mask(parameters, img, cellpose_mask)
            properties_df = functions.get_features_from_cellpose_seg(parameters, img, cellpose_mask, filename,
                                                                     output_path)

            print(properties_df.head())
            # properties_df.to_csv("image_%s.csv" % ind)
            properties_df.to_csv(output_path + filename.split(".")[0] + ".csv")

            if len(merged_properties_df.index) < 10:
                merged_properties_df = properties_df.copy()
            else:
                merged_properties_df = merged_properties_df.append(properties_df, ignore_index=True)

            summary_df.at[summary_ind, "folder_name"] = row["folder_name"]
            summary_df.at[summary_ind, "short_name"] = row["short_name"]
            summary_df.at[summary_ind, "filepath"] = filepath
            summary_df.at[summary_ind, "cell_number"] = len(np.unique(cellpose_mask))

            summary_df.to_csv(output_path_base + "summary_table" + ".csv")

            summary_ind += 1

            # merged_properties_df.to_csv(output_path + "merged.csv")
