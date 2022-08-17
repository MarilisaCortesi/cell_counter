# cell_counter
Script that can be used to segment and count the cells in images from transwell invasion assays. 

It takes as input the name of the folder in which the images are saved (line 57) and the kind of experiment (2D or 3D at line 61) and produces:
1. a PDF file in which the original DAPI and segmented images are compared (path and name to be set by user at line 58)
2. an excel spreadsheet (path and name to be set by user at line 59) containing the number of cells retrieved for each image
3. a pickle file  (path of folder to be set by user at line 60) for each image containing the different stages of the elaboration process.

In most cases the experiment kind needs to be set to 2D. The 3D mode was added to analyse images acquired using an organotypic model in which different kinds of cells are co-cultured and only one specific type (tagged with GFP) needs to be counted. When using the 3D mode the DAPI and GFP images need to be saved separately except for the name of the dye.
