library(usethis)
library(homologene)
library(sinew)

load("./data-raw/HUMAN_Background_data_raw.RData", verbose = T)

# Subset human ligand-receptor interactions to those with secreted ligands
uniprot_hsa_tab <- read.table("./data-raw/uniprot-filtered-organism-HSA.tab", sep = '\t', header = T, stringsAsFactors = FALSE, quote = "", comment.char = "")
uniprot_hsa_tab <- uniprot_hsa_tab[uniprot_hsa_tab$Subcellular.location..CC. != "",]
secreted <- uniprot_hsa_tab[base::grepl("Secreted",uniprot_hsa_tab$Subcellular.location..CC.),]
secreted <- do.call("c",lapply(secreted$Gene.names,function(x){
  out <- strsplit(x," ")[[1]]
  return(out)
}))
secreted <- secreted[!duplicated(secreted)]
secreted <- sort(secreted)

LR <- LR[which(LR$Ligand %in% secreted),]
Ligands <- intersect(Ligands,secreted)

Background_signaling_interactome$Source <- as.character(Background_signaling_interactome$Source)
Background_signaling_interactome$Target <- as.character(Background_signaling_interactome$Target)
intermediates$Gene <- as.character(intermediates$Gene)
non_interface_TFs$Gene <- as.character(non_interface_TFs$Gene)
TF_TF_interactions$Source <- as.character(TF_TF_interactions$Source)
TF_TF_interactions$Target <- as.character(TF_TF_interactions$Target)

# Rename human data objects
HUMAN_Background_signaling_interactome <- Background_signaling_interactome
HUMAN_dummy.var <- dummy.var
HUMAN_intermediates <- intermediates
HUMAN_Ligands <- Ligands
HUMAN_LR <- LR
HUMAN_non_interface_TFs <- non_interface_TFs
HUMAN_Receptors <- Receptors
HUMAN_TF_TF_interactions <- TF_TF_interactions
HUMAN_tf.db <- tf.db
HUMAN_tfs <- tfs

# Export data to "data" folder
usethis::use_data(HUMAN_Background_signaling_interactome, HUMAN_dummy.var, HUMAN_intermediates,
         HUMAN_Ligands, HUMAN_LR, HUMAN_non_interface_TFs, HUMAN_Receptors,
         HUMAN_TF_TF_interactions, HUMAN_tf.db, HUMAN_tfs, version = 3, overwrite = T)

# Get mouse orthologs of human genes in Background signaling interactome
backSig_Source <- homologene(HUMAN_Background_signaling_interactome$Source,inTax = 9606, outTax = 10090)
backSig_Target <- homologene(HUMAN_Background_signaling_interactome$Target,inTax = 9606, outTax = 10090)

tmp <- HUMAN_Background_signaling_interactome
tmp <- tmp[tmp$Source %in% backSig_Source$`9606` & tmp$Target %in% backSig_Target$`9606`,]

tmp <- merge(tmp,backSig_Source, by.x = "Source", by.y = "9606")
tmp <- merge(tmp,backSig_Target, by.x = "Target", by.y = "9606")

tmp <- tmp[,c(4,7,3)]
colnames(tmp) <- colnames(HUMAN_Background_signaling_interactome)
tmp <- tmp[order(tmp$Source,tmp$Target),]

MOUSE_Background_signaling_interactome <- tmp
MOUSE_Background_signaling_interactome <- MOUSE_Background_signaling_interactome[!duplicated(MOUSE_Background_signaling_interactome),]

# Get mouse orthologs of human genes in intermediates
inter_orth <- homologene(intermediates$Gene,inTax = 9606, outTax = 10090)

MOUSE_intermediates <- inter_orth[,"10090",drop = F]
colnames(MOUSE_intermediates) <- colnames(HUMAN_intermediates)
MOUSE_intermediates <- MOUSE_intermediates[!duplicated(MOUSE_intermediates),,drop = F]

# Get mouse orthologs of human genes in LR
LR_Ligand_orth <- homologene(HUMAN_LR$Ligand,inTax = 9606, outTax = 10090)
LR_Receptor_orth <- homologene(HUMAN_LR$Receptor,inTax = 9606, outTax = 10090)

tmp <- HUMAN_LR
tmp <- tmp[tmp$Ligand %in% LR_Ligand_orth$`9606` & tmp$Receptor %in% LR_Receptor_orth$`9606`,]

tmp <- merge(tmp,LR_Ligand_orth, by.x = "Ligand", by.y = "9606")
tmp <- merge(tmp,LR_Receptor_orth, by.x = "Receptor", by.y = "9606")

tmp <- tmp[,c(3,4,7)]
colnames(tmp) <- colnames(HUMAN_LR)
tmp$L_R <- paste(tmp$Ligand,tmp$Receptor,sep = "_")
tmp <- tmp[order(tmp$Ligand,tmp$Receptor),]

MOUSE_LR <- tmp
MOUSE_LR <- MOUSE_LR[!duplicated(MOUSE_LR),]

# Generate mouse ligands and receptor vectors

MOUSE_Ligands <- sort(unique(MOUSE_LR$Ligand))
MOUSE_Receptors <- sort(unique(MOUSE_LR$Receptor))

# New dummy var

MOUSE_dummy.var <- "Dummy"

# Get Mouse TF db and create TF-related objects

MOUSE_tf.db <- read.table("./data-raw/TFs_Mouse_AnimalTFDB.txt", header = T, sep = "\t", stringsAsFactors = F)

MOUSE_tfs <- sort(MOUSE_tf.db$Symbol)

## Get mouse orthologs of human genes in non_interface_TFs
non_interface_TFs_orth <- homologene(HUMAN_non_interface_TFs$Gene,inTax = 9606, outTax = 10090)

MOUSE_non_interface_TFs <- non_interface_TFs_orth[,"10090",drop = F]
colnames(MOUSE_non_interface_TFs) <- colnames(HUMAN_non_interface_TFs)
MOUSE_non_interface_TFs <- MOUSE_non_interface_TFs[!duplicated(MOUSE_non_interface_TFs),,drop = F]
MOUSE_non_interface_TFs <- MOUSE_non_interface_TFs[order(MOUSE_non_interface_TFs$Gene),,drop = F]
MOUSE_non_interface_TFs <- MOUSE_non_interface_TFs[MOUSE_non_interface_TFs$Gene %in% MOUSE_tfs,,drop = F]

## Get mouse orthologs of human genes in TF_TF_interactions
TF_TF_interactions_Source_orth <- homologene(HUMAN_TF_TF_interactions$Source,inTax = 9606, outTax = 10090)
TF_TF_interactions_Target_orth <- homologene(HUMAN_TF_TF_interactions$Target,inTax = 9606, outTax = 10090)

tmp <- HUMAN_TF_TF_interactions
tmp <- tmp[tmp$Source %in% TF_TF_interactions_Source_orth$`9606` & tmp$Target %in% TF_TF_interactions_Target_orth$`9606`,]

tmp <- merge(tmp,TF_TF_interactions_Source_orth, by.x = "Source", by.y = "9606")
tmp <- merge(tmp,TF_TF_interactions_Target_orth, by.x = "Target", by.y = "9606")

tmp <- tmp[,c(4,7,3)]
colnames(tmp) <- colnames(HUMAN_TF_TF_interactions)
tmp <- tmp[order(tmp$Source,tmp$Target),]

MOUSE_TF_TF_interactions <- tmp
MOUSE_TF_TF_interactions <- MOUSE_TF_TF_interactions[!duplicated(MOUSE_TF_TF_interactions),]
MOUSE_TF_TF_interactions <- MOUSE_TF_TF_interactions[MOUSE_TF_TF_interactions$Source %in% MOUSE_tfs & MOUSE_TF_TF_interactions$Target %in% MOUSE_tfs,,drop = F]

# Export data to "data" folder
usethis::use_data(MOUSE_Background_signaling_interactome, MOUSE_dummy.var, MOUSE_intermediates,
         MOUSE_Ligands, MOUSE_LR, MOUSE_non_interface_TFs, MOUSE_Receptors,
         MOUSE_TF_TF_interactions, MOUSE_tf.db, MOUSE_tfs, version = 3, overwrite = T)





