library(tidyr)
library(circlize)
library(reshape2)
suppressPackageStartupMessages(library(dendextend))

#input file is the output of RSAT's compare matrices
infile_read <- readLines(file.choose())

########### NCOR CALC AND PARSING

NCOR_clusters <- function(infile, na_thresh = NULL) {
  #### PARSING
  
  #remove all the information at the top of the file, all of it begins with ';'
  infile_edit <- infile[-(grep("^[;].*", infile))]
  #take out tab and '#' characters
  infile_edit <- gsub("\t", ",", infile_edit)
  infile_edit <- gsub("#", "", infile_edit)
  
  #column names
  infile_df_names <- unlist(strsplit(infile_edit[1], split = ","))
  
  #remove those column names from the file itself so it can actually just be data
  infile_edit <- infile_edit[-1]
  frame_rows <- rep(0, length(infile_edit))
  rows_parsed <- rep(0, length(infile_edit))
  
  ####
  
  
  #To put each of the infile rows as their own strings so a data frame can be made later
  for(i in 1:length(infile_edit)){
    frame_rows[[i]] <- infile_edit[[i]]
  }
  
  infile_df <- data.frame(frame_rows)
  
  #new_dat contains the information in a data frame
  #ordered_dat is ordered by Ncor value
  #uses tidyr::separate to separate out the row values (previously all in a string)
  new_dat <- infile_df %>% separate(frame_rows, infile_df_names, sep = ",", extra = "merge")
  ordered_dat <- new_dat[order(new_dat$Ncor, decreasing = TRUE), ]
  
  #sparse_dat contains names and ncor
  
  map_dat <- ordered_dat[,c(1:4)]
  sparse_dat <- ordered_dat[,c(3:4, 6) ]
  
  print(map_dat)
  
  #grabbing only the names of TFs that we inputted
  #we want to know relatedness of the things we put in, not to other ones like
  #STAT1_STAT2 in the case of measuring IRFs
  #we can include it if we want and then remove anything with a substantial amount of NAs
  lg_TF_list <- unique(c(sparse_dat[,1], sparse_dat[,2]))
  max_size <- length(lg_TF_list)
  lg_TF_list
  
  write_map <- unique(map_dat)
  print(write_map)
  #This is going to hold all of our NCORs
  
  mat_calc <- matrix(data = NA, nrow = max_size, ncol = max_size)
  rownames(mat_calc) <- lg_TF_list
  colnames(mat_calc) <- lg_TF_list
  
  #go through all of the rows in our sparse_matrix
  #populate with NCOR values
  
  for(i in 1:nrow(sparse_dat)){
    df_name1 <- sparse_dat[i,1]
    df_name2 <- sparse_dat[i,2]
  
    loc1 <- which(rownames(mat_calc) == df_name1)
    loc2 <- which(rownames(mat_calc) == df_name2)
    
    if(length(loc2) == 0 | length(loc1) == 0){next}else{
      if(is.na(mat_calc[loc1, loc2])){mat_calc[loc1, loc2] <- as.numeric(sparse_dat[i,3])}}
    
  
  } #end of the for loop

  diag(mat_calc) <- 1
  
  #replace any NAs with 0s for comparisons that did not meet the NCOR cutoff
  #use colMeans for both to keep it square
  
  if(is.null(na_thresh)){
    mat_calc[is.na(mat_calc)] <- 0
  }else{
    mat_calc <- mat_calc[-(which(colMeans(is.na(mat_calc)) > na_thresh)), -(which(colMeans(is.na(mat_calc)) > na_thresh))]
    mat_calc[is.na(mat_calc)] <- 0
  }

  #sim_measure is 1 - distance. So the smaller the number the smaller the distance, therefore the greater similarity
  sim_measure <- as.dist(1 - mat_calc)

  
  #hclust functions with different methods of clustering
  return(sim_measure)
  
}

##### BRANCH COLOR FUNCTION
## Function written by Alex Poole
color.dendro <- function(n, clusMember, labelColors, default.color="black") {
  a <- attributes(n)
  if( is.leaf(n) ) {		
    ## clusMember - a numeric vector designating leaf grouping
    ## labelColors - a vector of colors for the above grouping
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    if(a$label %in% names(clusMember)) {
      attr(n, "nodePar") <- c(a$nodePar, list(lab.col = labCol, pch=NA))
      attr(n, "edgePar") <- c(a$nodePar, list(col = labCol))
    }
    else {
      attr(n, "nodePar") <- c(a$nodePar, list(col = default.color, pch=NA))
      attr(n, "edgePar") <- c(a$nodePar, list(col = default.color))
    }
  }
  else {
    n[[1]] = color.dendro(n[[1]], clusMember, labelColors,
                          default.color = default.color)
    n[[2]] = color.dendro(n[[2]], clusMember, labelColors,
                          default.color = default.color)
    a1 <- attributes(n[[1]])
    a2 <- attributes(n[[2]])		
    left.branch.color <- as.character(a1$edgePar["col"])
    right.branch.color <- as.character(a2$edgePar["col"])
    #possible colors
    possibleColors = c(left.branch.color, right.branch.color)
    
    #this has been edited
    if (any(possibleColors != default.color)) {
      attr(n, "edgePar") <- c(a$nodePar, list(col = setdiff(possibleColors
                                                            ,default.color)))
    }
    else {
      attr(n, "edgePar") <- c(a$nodePar, list(col = default.color))
    }
    
  }
  return(n)
}

########### VISUALIZING INFILE DATA USING CIRCLIZE

#by default no flag to discard NAs for this example
sim <- NCOR_clusters(infile_read)

#create the dendrogram
dend <- sim %>% hclust(method = "average") %>% as.dendrogram

####based on coloring function example
labels <- attr(sim, "Labels")
populated1 = sample(c(1,2), 19, replace=TRUE)
names(populated1) = labels
#
# ## Color dendrogram branches based on majority class underneath branch.
cc = color.dendro(dend, populated1, c("black", "black"))


n <- length(labels)
#This is to cut the tree up into the same number of sections as labels
ct = cutree(cc, n)
dend_height = attr(cc, "height")
dend_labels <- labels[order.dendrogram(cc)]


#Drawing the p value rectangles (mock)
#help from http://www.r-graph-gallery.com/122-a-circular-plot-with-the-circlize-package/
x1 = c(0.00003, 0.00000323, 0.00043443, 0.000000005, 0.2323, 0.000000009, 0.008,
       0.0000004, 0.00004, 0.0000000000019, 0.00009, 0.0008, 0.005, 0.000566, 0.000001,
       0.00008, 0.00999, 0.00002, 0.00007)

xd <- -log10(x1 / max(x1))


circos.par("canvas.xlim" = c(-1.5, 1.5), "canvas.ylim" = c(-1.5, 1.5))
circos.par(track.margin=c(0,0)) 

#superimposition if needed for another track
circos.initialize(factors = "a", xlim = c(0, n))


color_rect <- c("white", "lightblue")[populated1]

circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.045,
             panel.fun = function(x, y) {
               xlim = get.cell.meta.data("xlim")
               ylim = get.cell.meta.data("ylim")
               i = get.cell.meta.data("sector.numeric.index")
               circos.rect(1:n - 0.8, rep(0,n), 1:n - 0.2, xd, col = color_rect)
             })


#This is for the labels


circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.4, 
             panel.fun = function(x, y) {
               for(i in seq_len(n)) {
                 circos.text(i-0.5, 0, dend_labels[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                            cex = 0.5)
               }
             })



#This is for the actual tree

circos.track(ylim = c(0, dend_height), bg.border = NA, 
             track.height = 0.5, panel.fun = function(x, y) {
               circos.dendrogram(cc, max_height = dend_height)
             })

circos.clear()



########### FULL JASPAR SET EXAMPLE NO COLORS

jaspar_full_dataset <- read.table(file = "Input.tab", sep = "\t", header=T, row.names = 1)
write(row.names(jaspar_full_dataset), file = "full_rownames.txt", sep = " ")
full_ex <- as.dist(1 - jaspar_full_dataset)
full_h <- hclust(full_ex, method = "complete")

full_dend <- as.dendrogram(full_h)

#trimming the labels
full_labels <- attr(full_ex, "Labels")
full_labels <- gsub("^.*_", "",full_labels)
full_labels <- gsub("\\.", "_",full_labels)
full_labels <- full_labels[order.dendrogram(full_dend)]

full_len <- as.numeric(length(full_labels))

CHD4_targets <- rep(0, full_len)

for(i in 1:dim(padj_calc)[1]){
  for (j in 1:full_len) {
    if(padj_calc[i,8] == full_labels[j]){
      CHD4_targets[j] <- 1
    }
  }
}
CHD4_targets <-  CHD4_targets + 1
names(CHD4_targets) <- full_labels
CC_full_dendro <- color.dendro(full_dend, CHD4_targets, c("black", "indianred2"))
full_height <- attr(full_dend, "height")

#begin circular visualization 
circos.par("canvas.xlim" = c(-0.7, 0.7), "canvas.ylim" = c(-0.7, 0.7))
circos.initialize(factors = "a", xlim = c(0, full_len))


#labels

circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
             panel.fun = function(x, y) {
               for(i in seq_len(full_len)) {
                 circos.text(i-0.5, 0, full_labels[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             cex = 0.2)
               }
             })

#dendro
circos.track(ylim = c(0, full_height), bg.border = NA, 
             panel.fun = function(x, y) {
               circos.dendrogram(CC_full_dendro, max_height = full_height)
             })

circos.clear()

########### HAGMAN DREME CHD4 DATA SET WITH RECTANGLES
DREME_file <- read.table(file = "Input.txt", 
                         sep = "\t", header = T, stringsAsFactors = F )


Mixed_data_parser <- function(input_data) {
  
  edit_data <- input_data
  edit_data <- edit_data[,c(1:4, 7:11, 35)]
}

padj_parser <- function(input_data) {
  
  edit_data <- input_data
  edit_data <- edit_data[,c(1:6, 8:11, 35, 45, 51)]
}

pos_neg_vals <- Mixed_data_parser(DREME_file)
padj_calc <- padj_parser(DREME_file)
padj_calc[is.na(padj_calc)] <- 0


#This is to get the unique DREME motifs
unique(pos_neg_vals[,1:5])

#This is to use the dreme files NCOR ids to get part of the full size dist matrix
numseq <- unique(as.numeric(pos_neg_vals[,5]))
NCOR_ids <- unique(as.vector(pos_neg_vals[,10]))
names_JP <- unique(as.vector(pos_neg_vals[,9]))
NCOR_ids<- NCOR_ids[-1]
names_JP <- names_JP[-1]

#using match to grab the first index to avoid duplicates
pos_vals <- (as.vector(pos_neg_vals[(match(NCOR_ids, pos_neg_vals[,10])),3]))
neg_vals <- (as.vector(pos_neg_vals[(match(NCOR_ids, pos_neg_vals[,10])),4]))

new_ex <- as.matrix(full_ex)

loc <- which(rownames(new_ex) %in% NCOR_ids)
pruned_mat <- new_ex[loc,loc]
CHD4_dist <- as.dist(1 - pruned_mat)
CHD4_hc <- hclust(CHD4_dist, method = "average")

TF_name_labels <- rep(0, length(names_JP))
p_track <- rep(0, length(pos_vals))
n_track <- rep(0, length(neg_vals))

for(j in 1:dim(pos_neg_vals[1])){
  for (i in 1:length(names_JP)) {
    if(rownames(pruned_mat)[i] %in% pos_neg_vals[j,10]){
      TF_name_labels[i] <- pos_neg_vals[j,9]
    }
  }
}

for(j in 1:dim(pos_neg_vals[1])){
  for (i in 1:length(names_JP)) {
    if(rownames(pruned_mat)[i] %in% pos_neg_vals[j,10]){
      p_track[i] <- ((pos_neg_vals[j,3]) / numseq)
      n_track[i] <- ((pos_neg_vals[j,4]) / numseq)
    }
  }
}

padj_vals <- rep(0, length(NCOR_ids))
rlog_vals <- rep(0, length(NCOR_ids))

for(j in 1:dim(padj_calc[1])){
  for (i in 1:length(padj_vals)) {
    if(rownames(pruned_mat)[i] %in% padj_calc[j,11]){
      if(padj_calc[j,12] <= 0.1){
        padj_vals[i] <- padj_calc[j,12]
        rlog_vals[i] <- padj_calc[j,13]
      }

    }
  }
}

padj_vals <- formatC(padj_vals, format = "e", digits = 2)
rlog_vals <- formatC(rlog_vals, format = "e", digits = 2)


d4 <- as.dendrogram(CHD4_hc)

#To maintain label order and orientation

name_order <- TF_name_labels[order.dendrogram(d4)]
padj_vals <- padj_vals[order.dendrogram(d4)]
rlog_vals <- rlog_vals[order.dendrogram(d4)]

labels(d4) <- name_order
CHD4_labels <- name_order
total_factors <- length(CHD4_labels)


##example coloring using same method as above coloring
populated = as.numeric(padj_vals[which(padj_vals > 0)])
for(i in 1:length(populated)){
  if(populated[i] > 0){
    populated[i] <- 2
  }else{
    populated[i] <- 1
  }
  
}

for(i in 1:length(padj_vals)){
  if(padj_vals[i] == "0.00e+00"){
    padj_vals[i] <- ""
  }
}

names(populated) = CHD4_labels
CHD4_cc = color.dendro(d4, populated, c("black", "indianred2"))
CHD4_height = attr(CHD4_cc, "height")



dev.new()
circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
circos.initialize(factors = "a", xlim = c(0, total_factors))

color_pos <- c("white", "indianred3")[populated]
color_neg <- c("white", "dodgerblue3")[populated]

circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
             panel.fun = function(x, y) {
               for(i in seq_len(total_factors)) {
                 circos.text(i-0.5, 0, padj_vals[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             cex = 0.5)
               }
             })


#pos and neg
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               circos.rect(1:total_factors - 0.8, rep(0,total_factors), 1:total_factors, p_track + 1, col = color_pos)
               circos.rect(1:total_factors - 0.8, rep(0,total_factors), 1:total_factors, -(n_track) - 1, col = color_neg)
               
             })

#labels
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
             panel.fun = function(x, y) {
               for(i in seq_len(total_factors)) {
                 circos.text(i-0.5, 0, CHD4_labels[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             cex = 0.5)
               }
             })

#dendro
circos.track(ylim = c(0, CHD4_height), bg.border = NA, 
             panel.fun = function(x, y) {
               circos.dendrogram(CHD4_cc, max_height = CHD4_height)
             })

circos.clear()
