import numpy as np
from skimage.future.graph import RAG

from polarityjam.polarityjam_logging import get_logger


def orientation_graph_nf(img):
    """Gets the RegionAdjacencyGraph for an instance segmentation """
    rag = RAG(img.astype("int"))
    rag.remove_node(0)
    return rag


def remove_islands(frame_graph, mask):
    """Remove unconnected cells (Cells without neighbours)."""
    # Get list of islands - nodes with no neighbours and remove them
    list_of_islands = []
    for nodes in frame_graph.nodes:
        if len(list(frame_graph.neighbors(nodes))) == 0:
            list_of_islands.append(nodes)

    get_logger().info("Removed number of islands: %s" % len(list_of_islands))

    # remove islands from image and graph
    for elemet in np.unique(list_of_islands):
        frame_graph.remove_node(elemet)
        mask[:, :][mask[:, :] == elemet] = 0

    return frame_graph, mask
