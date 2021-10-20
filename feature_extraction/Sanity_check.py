#Small python script to test functions included in repo, give an example of how they might be used
#Note that nearest neigbourhood function takes only one focal cell and a region adjcanceny grpah per run
#Think this is the way to support what ever sort of indexing/iteration over cells you want
#read in image
frame100 = vasc_masks[100,:,:]
# remove edges from image
frame100_no_edges = remove_edges(frame100)
#get graph of image with useful attributes assigned to node
frame100_graph = orientation_graph(frame100_no_edges)
#remove nodes and cell masks that are islands - have no neighbours!
no_singleton_graph,masks_no_islands = remove_islands(frame100_graph,frame100_no_edges)
#access a node with its corresponding mask label - the pixel value of the label
#do some measurements idk
nearest_neighbour_comparisions(no_singleton_graph,31)
#Carry out morans I on area
features_list, weights = morans_data_prep(no_singleton_graph,'area')
global_morans,local_morans = run_morans(features_list,weights)
