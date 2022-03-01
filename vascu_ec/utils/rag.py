from skimage.future.graph import RAG
from skimage.measure import regionprops


def orientation_graph_nf(img):
    """Gets the RegionAdjacencyGraph for an Image """
    rag = RAG(img.astype("int"))
    rag.remove_node(0)
    return (rag)


def orientation_graph(img):
    """DescribeMe"""
    rag = RAG(img.astype("int"))
    rag.remove_node(0)

    regions = regionprops(img.astype("int"))
    for region in regions:
        rag.nodes[region['label']]['orientation'] = region['orientation']
        rag.nodes[region['label']]['area'] = region['area']
        rag.nodes[region['label']]['polarity'] = region['major_axis_length'] / region['minor_axis_length']
        rag.nodes[region['label']]['aspect_ratio'] = (region["bbox"][2] - region["bbox"][0]) / (
                region["bbox"][3] - region["bbox"][1])
    return (rag)
