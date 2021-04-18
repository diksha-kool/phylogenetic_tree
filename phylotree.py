#!/usr/bin/python3.7
import re
import string
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.lines as lines

#####################################################################################
####                                  Functions                                  ####
#####################################################################################
## Function to return matrix value from the diagonal section 
def get_diag_data(l_matrix, r_idx, c_idx):
    #flip row<-->Col if Col is greater than the Row
    #confine within the diagonal matrix
    if(r_idx > c_idx):
        return (l_matrix[r_idx][c_idx])
    else:
        return (l_matrix[c_idx][r_idx])

## Function to return index of node in the P-tree
def search_spec(l_ptree, l_node):
    for idx in range(len(l_ptree)):
        if(l_ptree[idx][0] == l_node):
            return idx
    return -1

## Function traverse P-tree and find ordered list of primary set of species
def order_spec(l_ptree, l_node, l_ord_list):
    spec_idx = search_spec(l_ptree, l_node)
    (spec, parent0, parent1, dist) = l_ptree[spec_idx]

    #Check for parent0 if it is a leaf node (primary species)
    if(re.match(".+-.+",parent0)):
        #traverse sub-tree
        order_spec(l_ptree, parent0, l_ord_list)
    else:
        #append leaf node to ordered list
        l_ord_list.append(parent0)
        
    #Check for parent1 if it is a leaf node (primary species)
    if(re.match(".+-.+",parent1)):
        #traverse sub-tree
        order_spec(l_ptree, parent1, l_ord_list)
    else:
        #append leaf node to ordered list
        l_ord_list.append(parent1)
        
#####################################################################################
####                          Get input arguments                                ####
#####################################################################################
## Two agruments:
## 1. File path for DistanceMatrix data file
## 2. (Optional) - Debug message enabled
debug_en = 0
if (len(sys.argv) == 1):
    print("Please enter path for Input file!")
    print("Input file format:!")
    print("<species_0 <space separated dist values>!")
    print("<species_1 <space separated dist values>!")
    print("...")
    print("<species_n <space separated dist values>!")
    sys.exit()
if (len(sys.argv) >= 2):
    filepath = sys.argv[1]
if(len(sys.argv) >= 3):
    debug_en = int(sys.argv[2])
if (len(sys.argv) >= 4):
    print("Too many arguments!")
    sys.exit()

#####################################################################################
####                Read and load Species list, Distance matrix                  ####
#####################################################################################
spec_list = []
dist_matrix = []
d_matrix = []
#Open Input file and extract species list and distance matrix
print("Loading Distance Matrix...")
inputfilehnd = open(filepath)
for inputfileline in inputfilehnd:
    word_list = []
    word_find = re.findall("(\S+)", inputfileline)
    if(len(word_find) > 0):
        for word in word_find:
            word_list.append(word)
    spec_list.append(word_list[0])
    tmp = []
    for dist in word_list[1:]:
        tmp.append(float(dist))
    d_matrix.append(tmp)
inputfilehnd.close

#Check Species list and Distance Matrix dimentions
d_matrix_good = 0
if(len(spec_list) > 0) and (len(spec_list) == len(d_matrix)):
    d_matrix_good = 1
    for mat_row in d_matrix:
        if(len(mat_row) != len(spec_list)):
            d_matrix_good = 0
            break
if (d_matrix_good == 0):
    print(spec_list,"\n", d_matrix)
    print("Distance Matrix is not valid!")
    sys.exit()    

print("Valid Distance Matrix loaded!\n")

#Convert full matrix to diagonal matrix (reduce iterations)
for r_idx in range(len(d_matrix)):
    del d_matrix[r_idx][r_idx+1:]

#####################################################################################
####                      Reduce Matrix using UPGMA Algo                         ####
#####################################################################################
ptree_list = []
num_spec = len(spec_list)
#Max distance will be used later for visualization (image x-scaling)
max_dist = 0

print("Applying UPGMA Algorithm on the Distance Matrix!\n")
#number of iterations = Num Species-1
for loopidx in range(num_spec-1):
    #find min dist
    #assigning row=1 col=0 value as the min value
    min_dist = d_matrix[1][0]
    min_dist_sp1_id = 1 #Row index of min value
    min_dist_sp2_id = 0 #Col index of min value
    #Continue searching from row=2
    for r_idx in range(2, len(d_matrix)):
        #Searching in the diagonal matrix
        if(min(d_matrix[r_idx][:r_idx]) < min_dist):
            min_dist = min(d_matrix[r_idx][:r_idx])
            min_dist_sp1_id = r_idx
            min_dist_sp2_id = d_matrix[r_idx].index(min_dist)
    min_dist = round(min_dist/2, 10)
    #Capture max_distance in the last iteration
    if(loopidx == num_spec-2):
        max_dist = min_dist

    #Make entry in phylo-tree
    #Each entry is a tuple (new_species, parent_sp1, parent_sp2, distance)
    new_spec = spec_list[min_dist_sp1_id] + '-' + spec_list[min_dist_sp2_id]
    ptree_list.append((new_spec, spec_list[min_dist_sp1_id], spec_list[min_dist_sp2_id], min_dist))
    
    #Update the row and column of the diagonal matrix for parent_species with lower index
    #Delete the row and column of the diagonal matrix for parent_species with lower index    
    #|0
    #|  0
    #|----0 <- update row and column
    #|    | 0
    #|    |   0
    #|----------0 <- delete row and column
    #|    |     | 0
    #|____|_____|___0
    
    #Find row/col index to update, and row/col index to delete
    if(min_dist_sp1_id > min_dist_sp2_id):
        spec_idx_upd = min_dist_sp2_id
        spec_idx_del = min_dist_sp1_id
    else:
        spec_idx_upd = min_dist_sp1_id
        spec_idx_del = min_dist_sp2_id

    #Modify the spec_list at the update index
    spec_list[spec_idx_upd] = new_spec
    
    #Modify the d_matrix for the row and col at the update_index of the diagonal matrix
    #Don't touch any other values inside d_matrix
    for r_idx in range(len(d_matrix)):
        for c_idx in range(r_idx):
            if(r_idx == spec_idx_upd):
                #Update with new_distance = average(parent_sp1 x current_sp, parent_sp2 x current_sp)
                #Current_sp index = c_idx
                d_matrix[r_idx][c_idx] = (get_diag_data(d_matrix, min_dist_sp1_id, c_idx)+
                                          get_diag_data(d_matrix, min_dist_sp2_id, c_idx))/2
            elif(c_idx == spec_idx_upd):
                #Update with new_distance = average(parent_sp1 x current_sp, parent_sp2 x current_sp)
                #Current_sp index = r_idx
                d_matrix[r_idx][c_idx] = (get_diag_data(d_matrix, r_idx, min_dist_sp1_id)+
                                          get_diag_data(d_matrix, r_idx, min_dist_sp2_id))/2
            else:
                pass
            
    #Delete the redundant entry in spec_list, at the delete_index
    del spec_list[spec_idx_del]    
    #Delete the redundant column of the diagonal d_matrix, at the delete_index
    for r_idx in range(len(d_matrix)):
        for c_idx in range(r_idx):
            if (c_idx == spec_idx_del):
                del d_matrix[r_idx][c_idx]
    #Delete the redundant row of the diagonal d_matrix, at the delete_index
    del d_matrix[spec_idx_del]            

    #Print debug messages per iteration, if enabled
    if(debug_en == 1):
        print("********Interation******** # "+ str(loopidx + 1))
        print(spec_list)
        print(d_matrix)
        print("")

#print final Phylogenetic tree
print("**************************Final P-Tree**************************")
for node in ptree_list:
    print(node[1],"+",node[2],"->",node[0],"  Dist =", node[3])
print("")

#####################################################################################
####                Prepare Species Order list for Visualization                 ####
#####################################################################################
print("Preparing for visualization...")
spec_ord_list = []
#Traverse the ptree list, starting from the last entry
#Last entry is where all species combine
order_spec(ptree_list, ptree_list[-1][0], spec_ord_list)

#####################################################################################
####                         Visualization of Phylo Tree                         ####
#####################################################################################
#Visulizing P-Tree (figure x-range 0-1; y-range 0-1)
#setting figure margins and scaling ratios
x_margin = 0.05
y_margin = 0.05
#Scale figure x-axis to accomodate max_distance within the range of 0-1
#while leaving margin
x_scale = round((1-(2*x_margin))/max_dist, 13)
#Scale figure y-axis to accomodate all species within the range of 0-1
#while leaving margin
y_scale = round((1-(2*y_margin))/len(spec_ord_list), 13)
#Scale figure figure text font size to range between 4 to 18,
#for confortable viewing
font_scale = round((y_scale - 0.015)/0.00375, 0) + 4
font_scale = max(font_scale, 4)
font_scale = min(font_scale, 18)

#Create a figure object using matplotlib.pyplot
fig = plt.figure()
fig.dpi = 300

#Use visual node start dictionary for generating x1,y1 and x2,y2 of the line segments
#initialize visual node start dictionary for primary species ordered after traversal
vis_start_dict = {}
vis_y = y_margin
for spec in spec_ord_list:
    #Add primary species name text along with the staring points of line segments
    fig.text(0, vis_y+0.005, spec, fontsize=font_scale)
    vis_start_dict[spec] = (0, vis_y)
    vis_y = vis_y + y_scale

#For each node entry in ptree, add 3 line segments
#   ---------  <- seg1
#            | 
#            | <- seg3
#            |
#   ---------  <- seg2
for (node, parent0, parent1, dist) in ptree_list:
    #get the two parent species's x,y points in the graph
    (parent0_x, parent0_y)  = vis_start_dict[parent0]
    (parent1_x, parent1_y)  = vis_start_dict[parent1]
    #calculate current node species's x,y points based on parents
    node_x = round(dist*x_scale, 6)
    node_y = round((parent0_y + parent1_y)/2, 6)
    #make entry of the current node in the visual node start dictionary
    vis_start_dict[node] = (node_x, node_y)
    #fig.text(node_x, node_y+0.005, dist, fontsize=font_scale)
    #add the 3 line segments in the figure
    fig.add_artist(lines.Line2D([parent0_x, node_x], [parent0_y, parent0_y]))
    fig.add_artist(lines.Line2D([parent1_x, node_x], [parent1_y, parent1_y]))
    fig.add_artist(lines.Line2D([node_x, node_x], [parent0_y, parent1_y]))
    
plt.ioff()
plt.savefig("Phylotree.png")
print("Output file Phylotree.png is generated!")
