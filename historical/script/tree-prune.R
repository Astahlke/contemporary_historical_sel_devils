# install.packages("ape")
library("ape")

setwd("~/Code/devils/")


plot(maycollpaml, type = "phylogram")

Mitchell_2014_Supplemental_Trees <- read.csv("msu176_Supplementary_Data/14_0013_S1.txt", sep = "\t")
timecalibrated <- Mitchell_2014_Supplemental_Trees[7,'Tree..Newick.format.']

timecal_tree <- read.tree(text = timecalibrated)

allfour_species <- c("Monodelphis_domestica", "Sarcophilus_harrisii", "Macropus_eugenii", "Phascolarctos_cinereus")


pruned.tree <- drop.tip(timecal_tree, timecal_tree$tip.label[-match(allfour_species, timecal_tree$tip.label)])

plot(pruned.tree, show.tip.label=TRUE, use.edge.length = T, edge.width = 3)
edgelabels(pruned.tree$edge.length, bg="black", col="white", font = 2, height =3)
edgelabels(pruned.tree$edge.length, adj = c(0.5, -0.25), frame = "n")

sort(tree2_new$tip.label) #  "Dasyurus_maculatus" "Dasyurus_viverrinus" "Monodelphis_domestica" "Sarcophilus_harrisii"  "Phascolarctos_cinereus" ""Macropus_eugenii" "Thylacinus_cynocephalus"
kevins_species <- c("Dasyurus_maculatus", "Dasyurus_viverrinus", "Monodelphis_domestica", "Sarcophilus_harrisii", "Macropus_eugenii", "Phascolarctos_cinereus", "Thylacinus_cynocephalus")
wallkoala_species <- c("Sarcophilus_harrisii", "Macropus_eugenii", "Phascolarctos_cinereus") #173 12 107 90
oppkoala_species <- c("Monodelphis_domestica", "Sarcophilus_harrisii", "Phascolarctos_cinereus") #173 12 107 90

# which(tree2_new$tip.label=="Macropus_eugenii")

write.tree(pruned.tree, file = "~/Code/devils/add_koala/maycollado.oppkoaladevil.paml.nwk")
cophenetic(pruned.tree)
write.tree(pruned.tree)
plot(pruned.tree)
paml <- read.tree(text ="(MOD:82.48354,(SHA #1:66.69734,(PCI:52.77884, MEU:52.77882)));")
paml_3 <- read.tree(text="(Monodelphis_domestica:82.48354,(Sarcophilus_harrisii:66.69734,Macropus_eugenii:66.69732));")

plot(pruned.tree)
"(Monodelphis_domestica:82.48354,(Sarcophilus_harrisii:66.69734,(Phascolarctos_cinereus:52.77884,Macropus_eugenii:52.77882):13.9185):15.78622);"

### Graveyard
# MayColl <- read.nexus("PublishedPeerJall90concat.nex.con.tre")
# maycollpaml <- read.tree("~/Downloads/maycollpaml.txt")
# ens <- read.nexus("ensmbltree.nex")
# tree <- MayColl$con_all_compat
# species <-c("Monodelphis_domestica","Macropus_eugenii","Sarcophilus_harrisii", 
#            "Phascolarctos_cinereus_")
# 
# pruned.tree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
# write.tree(pruned.tree, file = "~/Code/devils/add_koala/maycollado.allfour.paml.nwk")
# 
# pruned.tree<-drop.tip(tree,tree$tip.label[-match(species[1:3], tree$tip.label)])
# write.tree(pruned.tree) #, file = "~/Code/devils/add_koala/maycollado.oppwalldevil.paml.nwk")
# 
# pruned.tree<-drop.tip(tree,tree$tip.label[-match(species[2:4], tree$tip.label)])
# write.tree(pruned.tree, file = "~/Code/devils/add_koala/maycollado.wallkoaladevil.paml.nwk")
# 
# pruned.tree<-drop.tip(tree,tree$tip.label[-match(species[c(1,3,4)], tree$tip.label)])
# write.tree(pruned.tree, file = "~/Code/devils/add_koala/maycollado.oppkoaladevil.paml.nwk")
# write.tree(maycollpaml)