# phylogenetic_tree
Python code to run UPGMA algorithm on WNT3 gene distance after alignment.

Usage: python3.7 phylotree.py <Distance Matrix file path> <optional: Debug message>
## Two agruments:
## 1. File path for DistanceMatrix data file
## 2. (Optional) - Debug message enabled
  
Distance Matric File Format:
<Species_0> <space separated distance matrix row values>
<Species_1> <space separated distance matrix row values>
...
<Species_n> <space separated distance matrix row values>

Code Description:
The Python code comprises of three main sections –
1.	Getting input arguments for the text file containing full (NxN) Distance Matrix. This section then reads the file and loads the Species list and Distance matrix to ‘spec_list[]’ and ‘dist_matrix[]’ arrays. The size of the matrix is checked for consistency. Finally, the full matrix is converted to a diagonal matrix is reduce iterations and duplicate work.
2.	Iterative reduction of the Distance Matrix using UPGMA Algorithm. This section will run n iterations (n = num_of_species - 1), where each iteration will find the minimum distance and combine two species. During combination, the new species name is generated as parent_species1+’-‘+ parent_species2. Distance Matrix is reduced 1 row and 1 column in each iteration by following two steps – Update Row and Column of one of the Parent_species (with lower index) to reflect the distances of the new species, and Delete Row and Column of the other Parent_species. Each time a combination of two species occur, a new entry is appended into the ‘ptree_list’ where the entry is a tuple containing (New_species, Parent_species1, Parent_species2, Distance)
3.	Visualization of the Phylogenetic Tree by using Matplotlib package. The diagram is plotted on a X-range of 0-1 and Y-range of 0-1, and therefore the scaling ratio for x-axis (based on Max Distance in Distance Matrix), y-axis (based on the number of species in the original Distance_Matrix) and font size (based on y-axis scaling). The ptree_list is traversed to get an ordered species list for starting x,y locations of the plot. Then for each entry in ptree_list, 3-line segments are added to the figure. Then figure is saved to a file named “Phylotree.png”.
