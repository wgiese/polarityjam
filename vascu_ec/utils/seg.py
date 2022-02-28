import cellpose.models
import skimage.segmentation


def get_cellpose_model(use_gpu):
    """Gets the cellpose default model"""
    return cellpose.models.Cellpose(gpu=use_gpu, model_type='cyto')


def get_cellpose_segmentation(parameters, im_seg):
    """Gets the cellpose segmentation"""
    model = get_cellpose_model(parameters["use_gpu"])
    if parameters["channel_nucleus"] >= 0:
        channels = [0, 1]
    else:
        channels = [0, 0]

    masks, flows, styles, diams = model.eval(im_seg, channels=channels)

    if parameters["clear_border"]:
        masks = skimage.segmentation.clear_border(masks)

    return masks