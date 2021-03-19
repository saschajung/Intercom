#' Background signaling interactome
#'
#' @name Background_signaling_interactome
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' Ligand receptor interactions
#'
#' @name LR
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' Ligands
#'
#' @name Ligands
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' Receptors
#'
#' @name Receptors
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' TF-TF interactions
#'
#' @name TF_TF_interactions
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' Dummy variable
#'
#' @name dummy.var
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' Intermediate signaling molecules
#'
#' @name intermediates
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' Non-interface TFs
#'
#' @name non_interface_TFs
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' TF database
#'
#' @name tf.db
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' TFs
#'
#' @name tfs
#' @docType data
#' @author Sascha Jung \email{sjung@@cicbiogune.es}
#' @keywords data
NULL

#' Reconstructs functional cell-cell communication networks
#' 
#' @description Reconstructs functional cell-cell communication networks for a given dataset.
#' @param data Matrix or data frame of expression values; rows have to correspond to genes and columns 
#' to cells. Row names **must** be gene symbols while column names can be any identifier for a cell.
#' Raw counts or normalized data is permitted.
#' @param anno.tbl Data frame annotating each cell with a cluster/cell type. Cell names have to be provided in
#' the first and cluster/cell type information in the second column.
#' @param species The organism from which the data was obtained. Currently, only "MOUSE" and "HUMAN" are supported
#' @param sighot.cutoff Cutoff parameter for SigHotSpotter. Can be between 0 and 1.
#' @param sighot.percentile Percentile parameter for SigHotSpotter. Can be between 0 and 100.
#' @param consv.thrs Fraction of cells per cluster/cell type that must express a ligand, receptor or TF to be considered.
#' @param ncores Number of cores to use. Depending on the machine, can be any integer value greater or equal to 1.
#' @param sig.cutoff Significance cutoff between 0 (weakest) and 1 (strictest). Default: 0.9
#' @param z.score.cutoff Cutoff parameter to determine significant associations between receptors and interface TFs.
#' Default: 2.
#' @param tissue.name A name of the dataset
#' @param out.path Path to a folder where the output should be stored.
#' @param temp.folder.name Name of the temporary folder to be created within the output path. Default: "temp"
#' @return List object containing the final interactome and auxiliary information
#' @export
InterCom <- function(data,anno.tbl,species,sighot.cutoff=0.1,sighot.percentile=70,consv.thrs=0.05,ncores=4,sig.cutoff=0.9,z.score.cutoff=2,tissue.name,temp.folder.name = "temp",out.path){
  
  #source("intercom_library_coexp.R")
  #source("for_plotting_network_functions.R")
  
  # suppressPackageStartupMessages({
  #   require(textshape)
  #   require(ggplot2)
  #   require(stutils)
  #   require(dplyr)
  #   require(doParallel)
  #   require(stringr)
  #   require(plyr)
  #   require(igraph)
  #   require(Matrix)
  #   require(reshape2)
  #   require(RSpectra)
  #   require(snow)
  #   require(taRifx)
  #   require(stats)
  #   require(gtools)
  #   require(data.table)
  #   require(rlist)
  # })
  
  system(base::paste0("mkdir ",out.path))
  system(base::paste0("mkdir ",out.path,"/",temp.folder.name))
  
  
  base::cat("Creating input parameters file\n\n")
  
  parms <- base::c("tissue.name","out.path","temp.folder.name","species","sighot.cutoff","sighot.percentile","consv.thrs","ncores")
  
  parms.file <- base::do.call(base::rbind,base::lapply(X = parms, function(x){
    return(base::paste0(x," = ",base::get(x)))
  }))
  
  utils::write.table(x = parms.file,file = base::paste0(out.path,"/","input_parameters_",base::Sys.Date(),".txt"), sep = "\t", row.names = F, quote = F,col.names = F)
  
  base::cat("Tissue : ",tissue.name,"\n")
  base::cat(" Preparing data\n")
  
  # Prepare data
  
  base::colnames(anno.tbl) <- base::c("cell.name","cell.type")
  
  celltype.freq <- base::as.data.frame(base::table(anno.tbl$cell.type))
  rm.celltype <- base::as.character(celltype.freq[base::which(celltype.freq$Freq <= 10),][,1])
  
  
  new.colnames <- base::sapply(base::colnames(data),function(x) {
    base::make.names(base::as.character(anno.tbl[["cell.type"]][x == base::as.character(anno.tbl[["cell.name"]])]),
               unique=FALSE)
  })
  base::colnames(data) <- new.colnames
  
  # Dynamic loading of maxsub function shared object required for running the get.max.cluster function.

  if(species == "MOUSE"){
    utils::data("MOUSE_Background_signaling_interactome","MOUSE_dummy.var","MOUSE_intermediates",
         "MOUSE_Ligands","MOUSE_LR","MOUSE_non_interface_TFs","MOUSE_Receptors","MOUSE_TF_TF_interactions",
         "MOUSE_tf.db","MOUSE_tfs",package = "InterCom")
  } else {
    if(species == "HUMAN"){
      utils::data("HUMAN_Background_signaling_interactome","HUMAN_dummy.var","HUMAN_intermediates",
           "HUMAN_Ligands","HUMAN_LR","HUMAN_non_interface_TFs","HUMAN_Receptors","HUMAN_TF_TF_interactions",
           "HUMAN_tf.db","HUMAN_tfs",package = "InterCom")
    } else {
      base::stop("Only the following species are supported: 'MOUSE', 'HUMAN'")
    }
  }
  utils::globalVariables(c("Background_signaling_interactome", "LR", "Ligands", "Receptors",
                            "TF_TF_interactions", "dummy.var", "intermediates",
                             "non_interface_TFs"))
  
  base::row.names(Background_signaling_interactome) <- 1:base::nrow(Background_signaling_interactome)
  non.rec.id <- base::which(Background_signaling_interactome$Source == dummy.var & !(Background_signaling_interactome$Target %in% Receptors))
  keep.id <- base::setdiff(base::row.names(Background_signaling_interactome),non.rec.id)
  Background_signaling_interactome <- Background_signaling_interactome[keep.id,]
  
  TF_TF_interactions <- taRifx::remove.factors(TF_TF_interactions)
  
  all.pops <- base::setdiff(base::unique(anno.tbl$cell.type),rm.celltype)
  anno.tbl <- anno.tbl[base::which(anno.tbl$cell.type %in% all.pops),]
  
  data.lig.exp <- get.gene.expr(exp.tbl = data,genes = base::intersect(Ligands,base::row.names(data)),cell.type = all.pops)
  base::colnames(data.lig.exp) <- base::paste0("Ligand.",base::colnames(data.lig.exp))
  
  L.frame <- dplyr::inner_join(x = data.lig.exp,y = LR[,-1], by = base::c("Ligand.gene" = "Ligand"))
  L.frame <- L.frame[base::which(L.frame$Ligand.exp.perc > consv.thrs),]
  
  base::save(list = base::c("data","anno.tbl","data.lig.exp","L.frame"),file = base::paste0(out.path,"/",temp.folder.name,"/temp_data.RData"))
  
  base::rm(data.lig.exp)
  
  all.pops <- base::setdiff(all.pops,base::c("Unknown","Uknown"))
  
  base::invisible(base::lapply(all.pops, function(celltype){
    cell.exp.tbl <- data[,base::which(base::colnames(data) == celltype)]
    base::saveRDS(object = cell.exp.tbl,file = base::paste0(out.path,"/",temp.folder.name,"/temp_",celltype,".Rds"))
    base::rm(cell.exp.tbl)
  }))
  
  base::rm(data)
  
  base::invisible(base::gc())
  
  base::lapply(X = all.pops,FUN = function(celltype){
    
    if(base::length(base::list.files(path = base::paste0(out.path,"/",temp.folder.name),pattern = base::paste0(celltype,"_results.RData"))) == 0){
      
      base::cat(base::paste0(" Celltype : ",celltype,"\n"))
      
      cell.exp.tbl <- base::readRDS(file = base::paste0(out.path,"/",temp.folder.name,"/temp_",celltype,".Rds"))
      
      base::cat("   Finding maximum sum subcluster in expression space\n")
      
      max.cluster.info <- get.cons.tfs(exp.tbl = cell.exp.tbl)
      if(base::class(max.cluster.info) == "list"){
        
        sig.input <- cell.exp.tbl[,max.cluster.info$tf.max.mat.cell]
        sig.input <- base::cbind.data.frame(base::row.names(sig.input),sig.input,stringsAsFactors = F)
        
        cons.tfs <- base::as.data.frame(base::unique(max.cluster.info$tf.count$Gene),stringsAsFactors = F)
        
        base::colnames(cons.tfs)[1] <- "Gene"
        
        cons.tfs$bool <- base::as.numeric(1)
        
        base::cat("   Starting SigHotSpotter analysis for the sub-cluster identified\n")
        
        hotspot.out <- SigHotSpotter_pipeline(idata = sig.input,species = species,cutoff = sighot.cutoff,DE_Genes = cons.tfs,percentile = sighot.percentile,ncores = ncores)
        
        base::save(list = base::c("max.cluster.info","sig.input","cons.tfs","hotspot.out"),file = base::paste0(out.path,"/",temp.folder.name,"/temp_",celltype,"_results.RData"))
        base::cat("   Saving results\n")
        
        base::rm(list = base::c("cell.exp.tbl","max.cluster.info","sig.input","cons.tfs","hotspot.out"))
        
      }else{
        base::cat(max.cluster.info,"\n")
        base::save(list = "max.cluster.info",file = base::paste0(out.path,"/",temp.folder.name,"/temp_",celltype,"_results.RData"))
        base::rm(max.cluster.info)
      }
    }
    

  })
  
  base::cat(" Collating Results\n")
  
  result.files <- base::list.files(path = base::paste0(out.path,"/",temp.folder.name),pattern = "_results.RData")
  
  all.pops <- base::as.character(base::sapply(result.files, function(file){
    tmp <- base::unlist(stringr::str_split(string = file,pattern = "_"))[-1]
    tmp <- tmp[-base::length(tmp)]
    return(base::paste(tmp,collapse="_"))
  }))
  
  collate <- base::lapply(X = all.pops,FUN = function(celltype){
    
    base::load(base::paste0(out.path,"/",temp.folder.name,"/temp_",celltype,"_results.RData"))
    
    if(base::class(max.cluster.info) == "list"){
      base::cat(" Collating data for ... ",celltype,"\n")

      out <- base::list()
      
      hotspot.recs <- base::unique(base::c(hotspot.out$active,hotspot.out$inactive))
      
      path.sums <- hotspot.out$path.sums
      
      path.sums <- path.sums[base::which(path.sums$z.score > 0),]
      
      iTF.targets <- hotspot.out$iTF.targets
      
      if(!base::is.null(path.sums)){
        R.TF.info <- dplyr::inner_join(x = path.sums, y = iTF.targets,by = "Gene")
        
        R.TF.info <- R.TF.info[base::which(R.TF.info$Receptor %in% hotspot.recs),]
        
        if(base::dim(R.TF.info)[1] > 0){
          
          cell.exp.tbl <- base::readRDS(file = base::paste0(out.path,"/",temp.folder.name,"/temp_",celltype,".Rds"))
          
          coexp.info <- gene.coexp(exp.tbl = sig.input[-1],gene.frame = R.TF.info[,base::c(1,2,5)],ncores = ncores)
          coexp.info$submat.coexp.perc <- coexp.info$coexp.count/(base::dim(sig.input)[2]-1)
          coexp.info$coexp.perc <- coexp.info$coexp.count/base::dim(cell.exp.tbl)[2]
          
          R.TF.info <- dplyr::inner_join(x = R.TF.info, y = coexp.info, by = base::c("Receptor","Gene","Target"))
          R.TF.info <- R.TF.info[base::which(R.TF.info$coexp.perc > consv.thrs),]
          
          out$maxmat <- max.cluster.info
          out$hotspot <- hotspot.out
          
          if(base::dim(R.TF.info)[1] > 0){
            
            out$R.TF.info <- R.TF.info
            
            rec.expr <- get.gene.expr(exp.tbl = cell.exp.tbl,genes = base::unique(R.TF.info$Receptor),cell.type = celltype)
            
            base::colnames(rec.expr) <- base::paste0("Receptor.",base::colnames(rec.expr))
            
            LR.frame <- dplyr::inner_join(x = L.frame,y = rec.expr, by = base::c("Receptor" = "Receptor.gene"))
            
            LR.frame <- dplyr::inner_join(x = LR.frame,y = R.TF.info,by = "Receptor")
            LR.frame <- base::unique(LR.frame[,base::c(1,2,4,6,7,9)])
            base::colnames(LR.frame) <- base::c("Lig.pop","Ligand","Lig.exp.perc","Receptor","Rec.pop","Rec.exp.perc")
            out$LR.frame <- LR.frame
            return(out)
            
          } else {
            return(NULL)
          }
        }else {
          return(NULL)
        }
      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
  })
  collate <- stats::setNames(object = collate,nm = all.pops)
  collate <- rlist::list.clean(.data = collate,fun = base::is.null)
  
  tissue.LR <- base::data.frame()
  tissue.R.TF <- base::data.frame()
  
  all.pops <- base::names(collate)
  
  for (Celltype in all.pops) {
    if(base::exists(x = "LR.frame",where = collate[[Celltype]])){
      tissue.LR <- base::rbind.data.frame(tissue.LR,collate[[Celltype]][["LR.frame"]],stringsAsFactors = F)
    }
    if(base::exists(x = "R.TF.info",where = collate[[Celltype]])){
      tissue.R.TF <- base::rbind.data.frame(tissue.R.TF,base::cbind.data.frame(Celltype,collate[[Celltype]][["R.TF.info"]],stringsAsFactors = F),stringsAsFactors = F)
    }
  }
  
  collate$tissue.LR <- tissue.LR
  collate$tissue.R.TF <- tissue.R.TF
  
  
  base::saveRDS(object = tissue.LR,file = base::paste0(out.path,"/tissue_LR_no_bootstrap.Rds"))
  utils::write.table(x = tissue.LR,file = base::paste0(out.path,"/tissue_LR_no_bootstrap.txt"), sep = "\t", row.names = F, quote = F)
  
  base::cat(" Find significant LR interactions\n")

  base::load(base::paste0(out.path,"/",temp.folder.name,"/temp_data.RData"))
  
  
  scored.LR <- scoringFun(data = data,tissue.R.TF = tissue.R.TF,tissue.LR = tissue.LR,LR = LR,sig.cutoff = sig.cutoff,z.score.cutoff = z.score.cutoff)

  collate$scored.LR <- scored.LR
  
  if(base::all(base::c("Lig.pop","Ligand","Receptor","Rec.pop") %in% base::colnames(tissue.LR))){
    
    final <- dplyr::inner_join(x = tissue.LR,y = scored.LR,by = base::c("Lig.pop","Ligand","Receptor","Rec.pop"))
    
    base::saveRDS(object = final,file = base::paste0(out.path,"/tissue_LR_scored.Rds"))
    
    utils::write.table(x = final,file = base::paste0(out.path,"/tissue_LR_scored.txt"), sep = "\t", row.names = F, quote = F)
    
    collate$final <- final
    
    output <- collate
    
    base::save(list = base::c("output"),file = base::paste0(out.path,"/output_",tissue.name,".RData"))
    
  } else{
    collate$final <- base::c()
  }
}
