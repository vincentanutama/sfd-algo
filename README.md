# sfd-algo
Finds neighbor sequence for Spatial First Difference estimator (Druckenmiller and Hsiang 2018). The algorithm starts with the leftmost polygon in a shapefile and samples neighboring polygons in the west-to-east direction. In the presence of multiple neighbors, the algorithm is greedy in that it extracts the longest sequence of contiguous neighbors. This greedy algorithm improves efficiency, defined as percentage of non-singleton polygons across all sampling channels.

MWE is provided using continental US counties from Database of Global Administrative Areas (GADM). To replicate, download the repo and run `example.R`. The code will call the algorithm in `f_sfdalgo.R`. The figure `sfdusa.png` plots the sequence of neighboring counties under different rotation schemes.
