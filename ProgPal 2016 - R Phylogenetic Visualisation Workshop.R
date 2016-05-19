# =======================================================
# = R Phylogenetic Visualisation Workshop, ProgPal 2016 =
# =                  Roger Close                        =
# =======================================================

# https://github.com/rclose/ProgPal-R-Phylogenetics-Visualisation-Workshop

# This workshop borrows heavily from tutorials by Liam Revell (and his phytools.org blog) and Lars Schmitz:
# http://www.phytools.org/eqg/Exercise_3.2/
# https://ecomorph.wordpress.com/2014/10/09/phylogenetic-trees-in-r-4/
# http://lukejharmon.github.io/ilhabela/instruction/2015/07/05/plotting-methods/ #some very fancy visualisation methods


# ==========================================================
# = Install these packages if necessary (uncomment first!) =
# ==========================================================
# install.packages("ape")
# install.packages("geiger")
# install.packages("phytools")
# install.packages("paleotree")

#load packages used in this workshop
library(ape) #the most important package for manipulating trees in R. Many core functions!
library(phytools) #contains numerous useful tools, especially for visualising trees and comparative data
library(geiger) #tools for modelling trait evolution
library(paleotree) #for timescaling trees

#clear the workspace
rm(list = ls()) #ls() tells you all the objects that are held in memory, and rm removes them. Confusingly, 'list =' specifies a character vector.

# ==============================================================================
# = Before we start, remember to cite the R packages you use. It's the law!    =
# ==============================================================================
#you can usually get the citation info like this (includes BibTeX entry):
citation("phytools")


# =======================================================
# = Let's kick things off by plotting a simulated tree! =
# =======================================================
tree <- rtree(n = 30) #rtree is a simple function for generating random trees. This generates a random tree with 30 tips.
#you can get help on R functions by prepending a question mark to the function name (uncomment the next line or just type it in manually):
? rtree
plot(tree, edge.width = 2, label.offset = 0.1) #R knows to use plot.phylo() on an object of class 'phylo' when you type plot()
plot(tree, edge.width = 2, label.offset = 0.1, direction = "up", srt = -45, cex = 0.7)
plot(tree, type = "cladogram", edge.width = 2, label.offset = 0.1)
plot(tree, type = "fan", edge.width = 2, label.offset = 0.1)
plot(tree, type = "unrooted", edge.width = 2, cex = 0.6, font = 1)


# ======================================
# = Learn about tree structures in ape =
# ======================================
#objects of class 'phylo' are lists. We can get a preview of the structure like this:
tree
class(tree)
str(tree)

#now let's specify our own simple tree, so we can understand the internal structure of a phylo object
tshirt <- read.tree(text = "(((((Glyptodon,Chelonia),Stenopterygius),(Papilio,Octopus)),(Whats_this_lol,((Brontosaurus,Tyto),(Dimetrodon,Thunnus)))),Anomalocaris);")

#this uses the parenthetical 'Newick' notation to encode the relationships
plot(tshirt, type = "cladogram", edge.width = 2, label.offset = 0.1) #label.offset creates a gap between the tips and the labels

#within the list masquerading as a 'phylo' object we have a vector containing the tip labels:
tshirt$tip.label

#there are m = n - 1 nodes for a fully bifurcating tree:
tshirt$Nnode

#and a matrix defining the edges (branches); numbers represent nodes:
tshirt$edge

#you can see how they relate to the tree structure by calling nodelabels() and tiplabels() without specifying any options
plot(tshirt, type = "cladogram", edge.width = 2, label.offset = 0.5) #label.offset creates a gap between the tips and the labels
tiplabels(frame = "circle") #the tips have node-numbers, too!
nodelabels(frame = "circle") #internal nodes
edgelabels(frame = "circle") #row-numbers for edges


# ======================================
# = Reading and writing trees to files =
# ======================================

#we can now save this tree for posterity:
write.tree(tshirt, "progpaltopology.tre")

#and we can look at the structure of the file:
cat(readLines("progpaltopology.tre"))

#and read it back in
tshirt <- read.tree("progpaltopology.tre") #see also read.nexus if you're reading in a tree in Nexus (.nex) format.

#it's also possible to read and write trees in Nexus format (requires phytools):
writeNexus(tree, "progpaltopology.nex")
cat(readLines("progpaltopology.nex"), sep = "\n") #look at the structure


# ===============================
# = Binary and non-binary trees =
# ===============================
#we can make a simple tree with a polytomy:
t1 <- read.tree(text = "((A,B,C),D);")

#and plot it:
plot(t1, type = "cladogram")

#is it fully dichotomous?
is.binary.tree(t1)

#multi2di randomly resolves polytomies in a multichotomous tree to produce a binary tree:
t2 <- multi2di(t1)
plot(t2, type = "cladogram")

#now it will be binary:
is.binary.tree(t2)

#(note that randomly resolving polytomies [e.g. in your consensus tree] is not necessarily a good idea!)


# ============================================
# = Ultrametric versus non-ultrametric trees =
# ============================================
t3 <- pbtree(b = 1, d = 0.5, n = 40)
plot(t3, show.tip.label = F)
is.ultrametric(t3)

#we can arbitrarily ultrametricise this tree using the Grafen (1989) method (see ?compute.brlen for more details)
t4 <- compute.brlen(t3, method = "Grafen", power = 1)
plot(ladderize(t4), show.tip.label = FALSE)
is.ultrametric(t4)

#this is a crude branch-length scaling method that's really only suitable for visualisation
par(mfrow = c(1,3))
plot(ladderize(t3), show.tip.label = F); plot(ladderize(t3), use.edge.length = F, show.tip.label = F); plot(ladderize(t4), show.tip.label = F)
graphics.off()

#compute.brlen also allows you to specify branch lengths using a numeric vector (in this case, all 1's)
t4 <- compute.brlen(t3, method = rep(1, nrow(t3$edge)), power = 1)
plot(t4, show.tip.label = F); axisPhylo()
#as branch lengths are often considered to represent the amount of change occurring along a branch (morphological, sequence substitutions, etc.) rather than time, equal branch lengths are sometimes said to represent purely speciational change



# ==========================================
# = Manipulating trees and phenotypic data =
# ==========================================

#rotating clades in plotting (useful for producing figures)
plot(tshirt, edge.width = 2, label.offset = 0.5); nodelabels()
rt.13 <- rotate(tshirt, 13) #rotate node 13
plot(rt.13, edge.width = 2, label.offset = 0.5); nodelabels()

#reroot a tree using a different outgroup taxon
rr.tshirt <- root(tshirt, outgroup = "Stenopterygius")
plot(rr.tshirt); nodelabels()

#binding trees together
plot(bind.tree(t1,t2))

#we could also drop the extinct tips, leaving only extant taxa
#first, we need to calculate the $root.time (especially if none of your tips are at the present), like so:
t3$root.time <- max(diag(vcv(t3))) #this computes the variance-covariance matrix of the branch lengths and finds the maximum value on the diagonal

#drop extinct tips (or 'leaves')
t5 <- dropExtinct(t3, tol = 0.01)
plot(t5, show.tip.label = F); axisPhylo()

#or drop the extant tips only
t6 <- dropExtant(t3, tol = 0.01)
plot(t6, show.tip.label = F); axisPhylo()

#we can also drop random or specific tips:
plot(drop.tip(phy = tshirt, tip = sample(x = tshirt$tip.label, size = 5, replace = FALSE))) #drop 5 random tips
plot(drop.tip(tshirt, c("Tyto","Octopus"))) #drop specific tips

#conversely, we could extract clades using a similar procedure
subtree <- extract.clade(phy = tshirt, node = getMRCA(phy = tshirt, tip = c("Octopus","Glyptodon")))
plot(subtree)


#matching phenotypic data to trees
x <- fastBM(t3) #simulate the evolution of a continuous trait on the tree (returns values at tips)
x #look at the trait values
all.equal(names(x), t3$tip.label) #are there any mismatches between the names of your trait and the tips of the tree?
#you can also use the function name.check in geiger:
name.check(t3, x)

#what if your tree contains fewer taxa than your trait data?
t3a <- drop.tip(t3, sample(t3$tip, size = 4)) #randomly drop some taxa from the tree
all.equal(names(x), t3a$tip.label) #do they match now?
setdiff(names(x), t3a$tip.label) #set operation: which elements are present in names(x) but not in the tree tips?
#we can fix that
xa <- x[t3a$tip.label]
all.equal(names(xa), t3a$tip.label) #they're now matched

#what if your tree contains more taxa than your trait data?
xb <- x[!(names(x) %in% names(sample(x, 4)))] #randomly drop some taxa from trait vector
all.equal(names(xb), t3$tip.label) #do they match now?
dropme <- setdiff(t3$tip.label, names(xb)) #taxa present in tree but not in trait data
t3b <- drop.tip(t3, dropme)
all.equal(names(xb), t3b$tip.label) #do they match?

#what if your data isn't in the same order as your tree's tips?
x <- sample(x)
x
all.equal(names(x), t3$tip.label) #do they match?
x <- x[t3$tip.label] #reorder the data
all.equal(names(x), t3$tip.label) #do they match?

#a fun example:
graphics.off()
par(mfrow = c(1,3))
phenogram(t3, x, col = "cornflowerblue"); title(main = "Extant and Extinct")
phenogram(t6, x[t6$tip.label], col = "darkred", ylim = range(x)); title(main = "Extinct")
phenogram(t5, x[t5$tip.label], col = "darkgreen", ylim = range(x)); title(main = "Extant")

#turn off the graphics device
dev.off()

#colouring branches by trait colour using phytools functions
fancyTree(t3,type = "phenogram95", x = x, spread.cost = c(1,0)) #phenogram with 95% confidence intervals on maximum likelihood ancestral-state values
plotBranchbyTrait(ladderize(t3), x, mode = "tips") #colour edges of tree according to ancestral states reconstructed from tip values
contMap(ladderize(t3), x, type = "fan") #similar method for plotting evolution of continuous character onto tree

# ==================================================
# = Timescale a tree using a posteriori algorithms =
# ==================================================
#there are many ways to timescale trees a posteriori. Functions in paleotree are quite useful.
#most functions require first- and last-occurrence dates:
timeData <- matrix(data = c(3,0,130,0,0,500,70,0,250,30,510,2,0,115,0,0,490,65,0,240,10,500), nrow = 11, ncol = 2, dimnames = list(tshirt$tip.label, c("FAD","LAD"))) #this is simply a matrix of first and last appearance dates with rownames set to tip label names
timeData

#this function call timescales our dubious tree using dubious tip-ages and the Minimum Branch Length (MBL) algorithm
timescaled.tree <- timePaleoPhy(tree = tshirt, timeData = timeData, type = "mbl", vartime = 10)

#we can now plot our timescaled tree
plot(timescaled.tree, edge.width = 2, label.offset = 0.5)

#we can also add a basic timescale using this function from ape:
axisPhylo()

#there are now more elements listed in the tree object
timescaled.tree$root.time

#these are the edge lengths (durations of each branch)
timescaled.tree$edge.length

#we can see how they map into the tree like so:
edgelabels(text = timescaled.tree$edge.length, frame = "circle", cex = 0.6)

#i.e., like this:
cbind(timescaled.tree$edge, timescaled.tree$edge.length)

#in the 'strap' and 'geoscale' packages there are handy functions for affixing geological timescales to your trees
geoscalePhylo(timescaled.tree, timeData, units = c("Period","Epoch"), tick.scale = "no", boxes = "Period", cex.tip = 0.7, quat.rm = T, cex.ts = 0.6)



# ===========================
# = Making your tree pretty =
# ===========================

#this is how to ladderise your tree --- very useful for publications
plot(ladderize(tshirt, right = F), direction = "up", label.offset = 0.1)
plot(ladderize(tshirt, right = T), direction = "up", label.offset = 0.1)

#this is how to colour a clade
totally.a.real.clade <- c("Stenopterygius","Chelonia","Glyptodon") #vector of tip-labels defining your clade of interest
clade.edges <- which.edge(tshirt, totally.a.real.clade) #get the rownames of the edges (i.e. branches) in the tree pertaining to your clade
clade.cols <- rep("black", nrow(tshirt$edge)) #a vector of colours - one for each branch in the tree
clade.cols
clade.cols[clade.edges] <- "cornflowerblue" #replace values for branches associated with our clade
clade.cols
plot(tshirt, edge.width = 5, edge.color = clade.cols, label.offset = 0.2) #plot it
cbind(tshirt$edge, clade.cols)
nodelabels(); tiplabels()

#you could also change the line type (or both):
clade.lty <- rep(1, nrow(tshirt$edge)) #a vector of colours - one for each branch in the tree
clade.lty[clade.edges] <- 3 #replace values for branches associated with our clade
plot(tshirt, edge.width = 5, edge.lty = clade.lty, edge.color = clade.cols, label.offset = 0.2) #plot it

#this is how to add (coloured) shapes to the tips
plot(tshirt, label.offset = 0.2)
rainbow.tip.cols <- rainbow(length(tshirt$tip.label))
names(rainbow.tip.cols) <- tshirt$tip.label
rainbow.tip.cols
tiplabels(pch = 23, cex = 2, bg = rainbow.tip.cols)
rainbow.tip.cols["Papilio"] <- "cornflowerblue"
plot(tshirt, label.offset = 0.2); tiplabels(pch = 23, cex = 2, bg = rainbow.tip.cols)


#this is how to add a box around a clade


#this is how to map a trait onto the branches of a tree


#this is a neat little function that makes a vector of gradient colours reflecting your trait values (in this case, branch duration)
color.gradient <- function(x, colors=c(lo.col,mid.col,hi.col), colsteps=20) {
	return(colorRampPalette(colors)(colsteps)[findInterval(x, seq(min(x),max(x), length.out = colsteps))])
}

edge.cols <- color.gradient(timescaled.tree$edge.length, colors = c("red","pink","blue"), colsteps = 20)

#you can then plot them
plot(timescaled.tree, edge.width = 4, edge.color = edge.cols)

#this is how you annotate a node with a label
plot(timescaled.tree, x.lim = c(0,1.25*max(nodeHeights(timescaled.tree))))
nodelabels(c(":(",":'("), c(19,14), frame = "none", adj = c(2,-1.5))
nodelabels("=)", 12, frame = "circle", cex = 2, bg = "cornflowerblue", col = "yellow")
cladelabels(timescaled.tree, node = 19, "\nDefinitely A Real Clade", offset = -2) #the '\n' inserts a new line to improve spacing
cladelabels(timescaled.tree, node = 14, "\nOMG!")
cladelabels(timescaled.tree, node = 12, "LOL!", offset = 5, orientation = "horizontal", cex = 2)

#this is how you can map images onto your plots (e.g. silhouettes)



# =======================================================
# = This is how to compare two trees (a 'cophylo' plot) =
# =======================================================
#useful if you want to examine difference between tree topologies
tree1 <- read.tree(text = "(((((Glyptodon,Chelonia),Stenopterygius),(Papilio,Octopus)),(Whats_this_lol,((Brontosaurus,Tyto),(Dimetrodon,Thunnus)))),Anomalocaris);")
tree2 <- read.tree(text = "(((((Brontosaurus,Tyto),(Chelonia,Stenopterygius)),(Glyptodon,Dimetrodon)),Thunnus),((Whats_this_lol,(Papilio,Anomalocaris)),Octopus));"); plot(tree2)
plot(cophylo(tree1,tree2,rotate = TRUE))


# ========================
# = 'multiPhylo' objects =
# ========================
#you can also store multiple trees in a single object. This is useful if you're performing analyses on a large posterior sample of trees from a Bayesian analysis, for example.
trees <- pbtree(n = 20, nsim = 36, scale = 1) #simulate 10 stochastic trees, each with 6 tips using a birth-death process
trees
str(trees) #look at the internal structure
trees <- roundBranches(trees, 1) #round branch-length values for ease of viewing structure

write.tree(trees, file = "example.trees")
cat(readLines("example.trees"), sep = "\n")

par(mfrow = c(6,6), mar = c(0,0,0,0)); lapply(trees, function(x) plot(ladderize(x), direction = "up", label.offset = 0.1, show.tip.label = FALSE, edge.color = sample(rainbow(30))[1]))
graphics.off()

#you can also combine several phylo objects using c()
c(t1,t2,t3,t4)


# ==========
# = ggtree =
# ==========
#http://cobra20.fhcrc.org/packages/release/bioc/vignettes/ggtree/inst/doc/ggtree.html
#to install:
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")
# biocLite("EBImage")
library(ggtree); library(ggplot2)
vignette("ggtree", package = "ggtree")

data(chiroptera)
gzoom(chiroptera, grep("Plecotus", chiroptera$tip.label))

ggtree(tree, aes(color = branch.length)) +
	scale_color_continuous(low = "blue", high = "green") +
	theme(legend.position = "bottom")

pp <- ggtree(tree) %>% phylopic("79ad5f09-cf21-4c89-8e7d-0c82a00ce728", color = "steelblue", alpha = .3)
pp

pp %>% phylopic("79ad5f09-cf21-4c89-8e7d-0c82a00ce728", color = "#86B875", alpha = .8, node = 4, width = 0.1) %>%
	phylopic("79ad5f09-cf21-4c89-8e7d-0c82a00ce728", color = "darkcyan", alpha = .8, node = 17, width = 0.1)

beast_file <- system.file("examples/MCC_FluA_H3.tree", package = "ggtree")
beast_tree <- read.beast(beast_file)

fasta <- system.file("examples/FluA_H3_AA.fas", package = "ggtree")
msaplot(ggtree(beast_tree), fasta)


