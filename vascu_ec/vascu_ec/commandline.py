import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from vascu_ec.feature_extraction import get_image_for_segmentation, get_features_from_cellpose_seg_multi_channel
from vascu_ec.utils.io import read_parameters, read_image, create_path_recursively, get_tif_list, read_key_file
from vascu_ec.utils.plot import plot_seg_channels, plot_cellpose_masks, set_figure_dpi
from vascu_ec.utils.seg import load_or_get_cellpose_segmentation
from vascu_ec.vascu_ec_logging import get_logger


# todo: copy parameter file to output folder.
# todo: add a logger handler to log into the file in the output folder.

# todo: Rename. Postponed for now due to missing junction-features.
# todo: run the app on a server? Ask for resources! Missing features?

def run(args):
    set_figure_dpi()

    # read args
    param_file = args.param
    filepath = args.in_file
    filename = args.filename_prefix
    output_path = args.out_path

    # print info
    get_logger().info(
        "Arguments provided: %s" % json.dumps(
            {"param": param_file, "in_file": filepath, "out_path": output_path, "filename": filename},
            sort_keys=True, indent=4
        )
    )

    # process args
    if not filename:
        filename, _ = os.path.splitext(os.path.basename(filepath))

    parameters = read_parameters(param_file)
    get_logger().info("Parameters: %s" % json.dumps(parameters, sort_keys=True, indent=4))

    # start routine
    _run(filepath, parameters, output_path, filename)


def _run(infile, parameters, output_path, fileout_name):
    create_path_recursively(output_path)

    # read input
    img = read_image(infile)
    img_seg = get_image_for_segmentation(parameters, img)

    # plot input
    plot_seg_channels(img_seg, output_path, fileout_name)

    # basic segmentation
    cellpose_mask = load_or_get_cellpose_segmentation(parameters, img_seg, infile)

    # plot cellpose mask
    plot_cellpose_masks(img_seg, cellpose_mask, output_path, fileout_name)

    # feature extraction
    properties_df = get_features_from_cellpose_seg_multi_channel(parameters, img, cellpose_mask, fileout_name, output_path)

    get_logger().info("Head of created dataset: \n %s" % properties_df.head())

    # write output
    fileout_base, _ = os.path.splitext(fileout_name)
    properties_df.to_csv(str(Path(output_path).joinpath(fileout_base + ".csv")))

    return properties_df, cellpose_mask


def run_stack(args):
    set_figure_dpi()

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

    # process
    parameters = read_parameters(param_file)
    get_logger().info("Parameters: %s" % json.dumps(parameters, sort_keys=True, indent=4))

    file_list = get_tif_list(inpath)

    for filepath in file_list:
        filepath = Path(filepath)
        filename = filepath.stem + filepath.suffix

        if not ((filepath.suffix != ".tif") or (filepath.suffix != ".tiff")):
            continue

        get_logger().info(
            "Processing file with: file stem  \"%s\" and file extension: \"%s\"" % (filepath.stem, filepath.suffix)
        )

        # start routine
        _run(filepath, parameters, output_path, filename)


def run_key(args):
    set_figure_dpi()

    # read args
    param_file = args.param
    in_path = args.in_path
    inkey = args.in_key
    output_path_base = args.out_path

    # print info
    get_logger().info(
        "Arguments provided: %s" % json.dumps(
            {"param": param_file, "in_key": inkey, "out_path": output_path_base},
            sort_keys=True, indent=4
        )
    )

    # convert
    output_path_base = Path(output_path_base)
    in_path = Path(in_path)

    # process
    parameters = read_parameters(param_file)
    get_logger().info("Parameters: %s" % json.dumps(parameters, sort_keys=True, indent=4))
    key_file = read_key_file(inkey)

    # empty DF summarizing overall results
    summary_df = pd.DataFrame()

    for index, row in key_file.iterrows():
        # current stack input sub folder
        cur_sub_path = str(row['folder_name'])
        if cur_sub_path.startswith(os.path.sep):
            cur_sub_path = cur_sub_path[1:-1]
        input_path = in_path.joinpath(cur_sub_path)

        # current stack output sub-folder
        cur_sub_out_path = str(row["short_name"])
        output_path = output_path_base.joinpath(cur_sub_out_path)

        # empty results dataset
        merged_properties_df = pd.DataFrame()

        file_list = get_tif_list(input_path)
        for filepath in file_list:
            filepath = Path(filepath)
            filename = filepath.stem + filepath.suffix

            if not ((filepath.suffix != ".tif") or (filepath.suffix != ".tiff")):
                continue

            get_logger().info(
                "Processing fil with: file stem  %s and file extension: %s" % (filepath.stem, filepath.suffix)
            )

            # single run
            properties_df, cellpose_mask = _run(filepath, parameters, output_path, filename)

            # append condition
            properties_df["condition"] = row["short_name"]

            if merged_properties_df.empty:
                merged_properties_df = properties_df.copy()
            else:
                merged_properties_df = pd.concat([merged_properties_df, properties_df], ignore_index=True)

            summary_df.at[index, "folder_name"] = row["folder_name"]
            summary_df.at[index, "short_name"] = row["short_name"]
            summary_df.at[index, "filepath"] = filepath
            summary_df.at[index, "cell_number"] = np.max(cellpose_mask)

            summary_df.to_csv(str(output_path_base.joinpath("summary_table" + ".csv")))
