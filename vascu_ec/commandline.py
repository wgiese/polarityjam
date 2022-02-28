from vascu_ec.feature_extraction import get_image_for_segmentation, get_features_from_cellpose_seg
from vascu_ec.logging import get_logger
from vascu_ec.utils.io import read_parameters, read_image
import matplotlib.pyplot as plt
import numpy as np

from vascu_ec.utils.seg import get_cellpose_segmentation


def run(args):
    #mpl.rcParams['figure.dpi'] = 300


    get_logger().info("Read parameters from: %s" % args.param)
    parameters = read_parameters(args.param)

    get_logger().info("Parameters: %s" % parameters)

    filename = parameters['image_filename']
    input_path = parameters['input_path']
    output_path = parameters['output_path']

    filepath = input_path + filename
    img = read_image(filepath)
    img_seg = get_image_for_segmentation(parameters, img)

    if parameters["channel_nucleus"] >= 0:
        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(img_seg[0, :, :])
        ax[1].imshow(img_seg[1, :, :])
        plt.savefig(output_path + filename + "-seg.png")
    else:
        fig, ax = plt.subplots()
        ax.imshow(img_seg[:, :])
        plt.savefig(output_path + filename + "-seg.png")

    cellpose_mask = get_cellpose_segmentation(parameters, img_seg)
    get_logger().info("Number of cell labels: %s" % np.max(cellpose_mask))

    properties_df = get_features_from_cellpose_seg(parameters, img, cellpose_mask, filename, output_path)

    get_logger().info(properties_df.head())

    properties_df.to_csv(output_path + filename.split(".")[0] + ".csv")


def run_stack(args):
    import numpy as np
    import matplotlib.pyplot as plt
    import numpy as np
    import skimage.io

    mpl.rcParams['figure.dpi'] = 300
    import sys
    import os
    import glob

    sys.path.append("../vascu_ec/")
    import functions

    if len(sys.argv) > 1:
        parameter_file = sys.argv[1]
    else:
        if sys.platform.startswith("win"):
            parameter_file = "..\\vascu_ec\\local\\parameters.yml"
        else:
            parameter_file = "../vascu_ec/local/parameters.yml"
    print("Read parameters from: %s" % parameter_file)

    parameters = functions.read_parameters(parameter_file)
    input_path = parameters["input_path"]

    file_list = glob.glob(input_path + "*.tif")
    # merged_properties_df = pd.DataFrame()

    for ind, filepath in enumerate(file_list):

        output_path = parameters['output_path']

        if sys.platform.startswith("win"):
            print("processing", filepath.split("\\")[-1].split(".")[0])
            filename = filepath.split("\\")[-1]
        else:
            print("processing", filepath.split("/")[-1].split(".")[0])
            filename = filepath.split("/")[-1]

        filepath_, file_extension = os.path.splitext(filepath)
        print("filename_: %s" % filepath_)
        print("file extension: %s" % file_extension)

        if not ((file_extension != ".tif") or (file_extension != ".tiff")):
            continue

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

        if os.path.exists(filepath_ + "_seg.npy"):
            # in case an annotated mask is available
            cellpose_seg = np.load(filepath_ + "_seg.npy", allow_pickle=True)
            cellpose_mask = cellpose_seg.item()['masks']
            if parameters["clear_border"]:
                cellpose_mask = skimage.segmentation.clear_border(cellpose_mask)
        else:
            cellpose_mask = functions.get_cellpose_segmentation(parameters, img_seg)
        print("Number of cell labels: ")
        print(np.max(cellpose_mask))
        # nuclei_mask = functions.get_nuclei_mask(parameters, img, cellpose_mask)
        # golgi_mask = functions.get_golgi_mask(parameters, img, cellpose_mask)
        properties_df = functions.get_features_from_cellpose_seg(parameters, img, cellpose_mask, filename, output_path)

        print(properties_df.head())
        # properties_df.to_csv("image_%s.csv" % ind)
        properties_df.to_csv(output_path + filename.split(".")[0] + ".csv")

        # if len(merged_properties_df.index) < 10:
        #    merged_properties_df = properties_df.copy()
        # else:
        #    merged_properties_df = merged_properties_df.append(properties_df, ignore_index=True)

        # merged_properties_df.to_csv(output_path + "merged.csv")


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
