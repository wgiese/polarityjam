from skimage.future.graph import RAG


def orientation_graph_nf(img):
    """Gets the RegionAdjacencyGraph for an instance segmentation """
    rag = RAG(img.astype("int"))
    rag.remove_node(0)
    return rag
