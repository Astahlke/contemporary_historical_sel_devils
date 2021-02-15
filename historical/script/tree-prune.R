library("ape")

Mitchell_2014_Supplemental_Trees <- read.csv("msu176_Supplementary_Data/14_0013_S1.txt", sep = "\t")
timecalibrated <- Mitchell_2014_Supplemental_Trees[7,'Tree..Newick.format.']

timecal_tree <- read.tree(text = timecalibrated)

allfour_species <- c("Monodelphis_domestica", "Sarcophilus_harrisii", "Macropus_eugenii", "Phascolarctos_cinereus")
wallkoala_species <- c("Sarcophilus_harrisii", "Macropus_eugenii", "Phascolarctos_cinereus") #173 12 107 90
oppkoala_species <- c("Monodelphis_domestica", "Sarcophilus_harrisii", "Phascolarctos_cinereus") #173 12 107 90

pruned.tree <- drop.tip(timecal_tree, timecal_tree$tip.label[-match(allfour_species, timecal_tree$tip.label)])

plot(pruned.tree, show.tip.label=TRUE, use.edge.length = T, edge.width = 3)
edgelabels(pruned.tree$edge.length, bg="black", col="white", font = 2, height =3)
edgelabels(pruned.tree$edge.length, adj = c(0.5, -0.25), frame = "n")

