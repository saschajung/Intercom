base::set.seed(100)

############################
#  Intercom Functions
############################

#' Get Gene Expression Subset
#'
#' The function derives a subset of a gene expression data frame for a given list of genes and optionally
#' cell types/clusters
#'
#' @param exp.tbl A data frame of gene expression values. Row names must be genes, column names must be cell types/clusters
#' @param genes A character vector of genes
#' @param cell.type A character vector of cell types/clusters
#' @return Data frame of gene expression values
#' @export
get.gene.expr <- function(exp.tbl,genes,cell.type=NULL){
  gene.exp.tbl <- exp.tbl[genes,,drop=FALSE]
  
  all.pops <- cell.type
  
  if(base::length(all.pops) == 1){
    cell.gene.exp <- gene.exp.tbl[,base::which(base::grepl(x = base::colnames(gene.exp.tbl),pattern = base::paste0("^",cell.type,"[\\.0-9]*$"),ignore.case = F)),drop=FALSE]
    cell.gene.abs.exp <- base::rowSums(cell.gene.exp)
    cell.gene.abs.exp <- base::cbind.data.frame(base::row.names(cell.gene.exp),cell.gene.abs.exp,stringsAsFactors = F)
    base::colnames(cell.gene.abs.exp) <- c("gene","abs.expr")
    
    cell.gene.bool <- base::as.data.frame(bool.data(exp.tbl = cell.gene.exp))
    
    cell.gene.cons <- base::rowSums(cell.gene.bool)
    cell.gene.cons <- base::cbind.data.frame(base::row.names(cell.gene.bool),cell.gene.cons,stringsAsFactors = F)
    base::colnames(cell.gene.cons) <- c("gene","cell.count")
    cell.gene.cons$exp.perc <- cell.gene.cons$cell.count/base::dim(cell.gene.bool)[2]
    cell.gene.out <- dplyr::inner_join(x = cell.gene.cons,y = cell.gene.abs.exp, by = "gene")
    cell.gene.out$celltype <- cell.type
    out <- cell.gene.out[,base::c(5,1:4)]
  }else{
    
    out <- base::do.call(base::rbind,base::lapply(X = all.pops,FUN = function(celltype1){
      cell.gene.exp <- gene.exp.tbl[,base::which(base::grepl(x = base::colnames(gene.exp.tbl),pattern = base::paste0("^",celltype1,"[\\.0-9]*$"),ignore.case = F))]
      
      cell.gene.abs.exp <- base::rowSums(cell.gene.exp) 
      cell.gene.abs.exp <- base::cbind.data.frame(base::row.names(cell.gene.exp),cell.gene.abs.exp,stringsAsFactors = F)
      base::colnames(cell.gene.abs.exp) <- c("gene","abs.expr")
      
      cell.gene.bool <- bool.data(exp.tbl = cell.gene.exp)
      
      cell.gene.cons <- base::rowSums(cell.gene.bool)
      cell.gene.cons <- base::cbind.data.frame(base::row.names(cell.gene.bool),cell.gene.cons,stringsAsFactors = F)
      base::colnames(cell.gene.cons) <- c("gene","cell.count")
      cell.gene.cons$exp.perc <- cell.gene.cons$cell.count/base::dim(cell.gene.bool)[2]
      
      cell.gene.out <- dplyr::inner_join(x = cell.gene.cons,y = cell.gene.abs.exp, by = "gene")
      cell.gene.out$celltype <- celltype1
      return(cell.gene.out[,base::c(5,1:4)])
      
    }))
  }
  base::invisible(base::gc())
  return(out)
}


#' Booleanize gene expression data
#' 
#' @description The function booleanizes the expression matrix based on expression threshold. 
#' It also removes genes not expression in any cell.
#' @param exp.tbl Gene expression dataframe.
#' @param expr.thrs Expression threshold.
#' @return Gene expression data frame with booleanized expression.
#' @export
bool.data <- function(exp.tbl,expr.thrs=0){
  bool.tbl <- base::matrix(data = base::as.numeric(exp.tbl > expr.thrs),nrow = base::nrow(exp.tbl),ncol = base::ncol(exp.tbl))
  base::colnames(bool.tbl) <- base::colnames(exp.tbl)
  base::rownames(bool.tbl) <- base::rownames(exp.tbl)
  bool.tbl <- bool.tbl[base::which(base::rowSums(bool.tbl) > 0),,drop=FALSE]
  return(bool.tbl)
}

#' Get conserved TFs
#' 
#' @description The function computes TFs expressed in a given percentage of cells.
#' @param exp.tbl Gene expression dataframe.
#' @param quantile Minimum fraction of cells expressing the Tf (default: 0.95)
#' @return List containing all TFs and data frame including only TFs expressed in the required fraction of cells
#' @export
get.cons.tfs <- function(exp.tbl,quantile = 0.95){
  tf.tbl <- exp.tbl[base::which(base::rownames(exp.tbl) %in% .pck_env$tfs),]
  tf.bool <- bool.data(exp.tbl = tf.tbl)
  freq_df <- base::apply(tf.bool,1,sum)
  freq_df <- base::data.frame(Gene = base::names(freq_df), Freq = freq_df,stringsAsFactors = FALSE)
  cutoff <- base::unname(stats::quantile(freq_df$Freq[freq_df$Freq > 0],quantile))
  return(base::list(tf.max.mat.cell = base::colnames(exp.tbl), tf.count = freq_df[freq_df$Freq >= cutoff,]))
}


#' Compute Co-expression Between Genes
#' 
#' @description The function calculates the coexpression of genes in each row of a given data frame using the gene expression table.
#' @param exp.tbl Gene expression matrix/data frame to calculate coexpression of genes.
#' @param gene.frame A data frame with genes. The coexpression is calculated for genes in each row.
#' @param ncores Number of cores to be used (Default : 4) 
#' @return gene.frame data frame with an additional column corresponding to the number of cells coexpressing the genes. 
#' @export
gene.coexp <- function(exp.tbl,gene.frame,ncores=4){
  gene.frame <- taRifx::remove.factors(gene.frame)
  base::row.names(gene.frame) <- 1:base::nrow(gene.frame)
  bool.exp.tbl <- bool.data(exp.tbl = exp.tbl)
  out <- base::do.call(base::rbind,base::lapply(X = 1:base::nrow(gene.frame), FUN = function(i){
    genes <- base::as.character(gene.frame[i,])
    len <- base::length(dplyr::intersect(genes,base::row.names(bool.exp.tbl)))
    if (len == base::length(genes)) {
      bool.gene.tbl <- bool.exp.tbl[genes,]
      bool.sum <- base::colSums(bool.gene.tbl)
      coexp.count <- base::length(bool.sum[bool.sum == base::ncol(gene.frame)])
      out.row <- base::cbind.data.frame(gene.frame[i,],coexp.count,stringsAsFactors = F)
      return(out.row)
    }else{
      return(NULL)
    }
  }))
  return(out)
}


############################
#  SighotSpotter Functions
############################


#' General Pipeline Call to SigHotSpotter
#'
#' The function computes compatibility scores for signaling intermediates
#'
#' @param species Currently supported species: "HUMAN", "MOUSE"
#' @param idata Data frame of input data. Rows correspond to genes, columns to cells
#' @param cutoff Maximum number of zero-value genes, above this cutoff the genes are excluded
#' @param DE_Genes (Differential) expression dataset as a one column data frame (1 for up-regulated, -1 for down-regulated genes)
#' @param percentile Predicted intermediates are taken into account above this threshold
#' @param invert_DE If the differential expression should be inverted (default: FALSE)
#' @param showprogress shows progress bar in shiny app if set to TRUE, set it to FALSE in batch mode without GUI (default: TRUE)
#' @param ncores Number of cores to use (default: 4)
#' @return Compatibility scores
#' @export
SigHotSpotter_pipeline <- function(species, idata, cutoff, DE_Genes, percentile, invert_DE = FALSE, showprogress = TRUE,ncores=4){
  
  subg=Data_preprocessing(input_data = idata,cutoff = cutoff,species = species)
  
  ## Calculate stationary distribution of the MC
  Steady_state_true=Markov_chain_stationary_distribution(subg)
  
  prob.matrix <- Matrix::summary(Steady_state_true$prob.matrix)
  
  edge.id.name <- base::cbind.data.frame(igraph::as_edgelist(graph = subg,names = T), igraph::as_edgelist(graph = subg,names = F),stringsAsFactors = F)
  base::colnames(edge.id.name) <- c("a","b","c","d")
  
  prob.matrix <- dplyr::inner_join(x = edge.id.name,y = prob.matrix, by = c("c" = "i" , "d" = "j"))
  prob.matrix <- prob.matrix[,base::c(1,2,5)]
  
  Steady_state_true <- Steady_state_true$SD
  
  ## Retrieves high probability intermediates
  int=high_probability_intermediates(x = Steady_state_true, intermediates = .pck_env$intermediates,percentile =  percentile)
  gintg=integrate_sig_TF(g = subg,x = Steady_state_true,deg = DE_Genes, non_interface_TFs = .pck_env$non_interface_TFs,TF_TF_interactions = .pck_env$TF_TF_interactions )
  
  # nTF=nonterminal_DE_TFs(g = gintg,deg = DE_Genes,non_interface_TFs = non_interface_TFs)
  
  target.TFs <- dplyr::inner_join(x = DE_Genes, y = .pck_env$TF_TF_interactions, by = c("Gene" = "Source"))
  target.TFs <- target.TFs[base::which(target.TFs$Target %in% igraph::V(gintg)$name),]
  if(base::nrow(target.TFs) == 0){
    return(NULL)
  }
  target.coexp <- gene.coexp(exp.tbl = idata,gene.frame = target.TFs[,base::c(1,3)],ncores = ncores)
  target.coexp <- dplyr::inner_join(x = target.coexp, y = target.TFs[,-2], by = c("Gene", "Target"))
  target.coexp$perc <- target.coexp$coexp.count/base::dim(idata)[2]
  target.coexp <- target.coexp[base::which(target.coexp$perc > 0.05 & target.coexp$Effect == 1),]
  iTF.target.info <- target.coexp[,base::c(1,2,4)]
  nTF <- iTF.target.info$Target
  
  nTF_scoring <- base::unique(nTF)
  base::names(nTF_scoring) <- nTF_scoring
  
  if (base::length(int) == 0){
    base::cat("No intermediates found. You may decrease the percentile in order to find intermediates. .\n")
    return(NULL)
  }else if (base::length(nTF) == 0){
    base::cat("No non-terminal conserved TFs found at lower conservation also. .\n")
    nTF <- DE_Genes
    return(NULL)
  }
  ## Computing compatibility scores
  score <- base::lapply(nTF_scoring,function(x){return(comp_score_tf(x,int,gintg))})
  score <- base::lapply(nTF,function(x){score[[x]]})
  if(base::is.null(score)){
    base::cat("No shortest path found. You may decrease the cutoff in order to find shortest path. .\n")
    return(NULL)
  }else{
    #converting the nested list into a matrix whose row sum will give probability of each intermediate
    score_m=(base::matrix(base::unlist(score), ncol=base::length(score), byrow=F))
    score_m_means=base::as.list(base::rowMeans(score_m))
    final_score=compatability_score(score_m_means,Steady_state_true,int)
    
    toiintA <- base::as.character(final_score[base::which(final_score$Activation_probability > 0.5),]$Gene)
    toiintI <- base::as.character(final_score[base::which(final_score$Activation_probability < 0.5),]$Gene)
    
    
    #pruning the integrated networks
    gintg.p=prun.int.g(gintg)
    
    #building networks for all intermediates for active signaling hotspots
    sp_int_A <- base::lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes,.pck_env$non_interface_TFs)
    
    #building networks for inactive signaling hotspots
    sp_int_I <- base::lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes,.pck_env$non_interface_TFs)
    
    #retrieve receptors for active intermediates & inactive:
    u.gr <- base::Reduce(igraph::graph.union,sp_int_A)
    if(base::class(u.gr) == "igraph"){
      aa <- igraph::incident(u.gr,"NICHE",mode="out")
      bb <- igraph::ends(u.gr,aa)
      active_receptors <- bb[,2]
      active_receptors <- dplyr::intersect(x = active_receptors,y = .pck_env$LR$Receptor)
    }else{
      active_receptors <- NULL
    }
    
    u.gr <- base::Reduce(igraph::graph.union,sp_int_I)
    if(base::class(u.gr) == "igraph"){
      aa <- igraph::incident(u.gr,"NICHE",mode="out")
      bb <- igraph::ends(u.gr,aa)
      inactive_receptors <- bb[,2]
      inactive_receptors <- dplyr::intersect(x = inactive_receptors, y = .pck_env$LR$Receptor)
    }else{
      inactive_receptors <- NULL
    }
    
    
    ### Link signaling molecules to receptors
    del <- igraph::incident(gintg.p,"NICHE",mode = "out")
    gintg.p.noNiche <- igraph::delete.edges(gintg.p,del)
    dists <- igraph::distances(gintg.p.noNiche,to = base::as.character(final_score$Gene), v = base::as.character(base::unique(base::c(active_receptors,inactive_receptors))), mode = "out", weights = NA, algorithm = "johnson")
    dists <- reshape2::melt(dists)
    dists <- dists[base::which(!base::is.infinite(dists$value)),]
    edge_df <- igraph::ends(gintg.p.noNiche,igraph::E(gintg.p.noNiche))
    weights <- base::as.numeric(base::apply(idata[edge_df[,1],-1],1,base::mean))*base::as.numeric(base::apply(idata[edge_df[,2],-1],1,base::mean))
    weights[base::is.na(weights)] <- 0
    g_test <- igraph::set_edge_attr(gintg.p.noNiche,"weight", value = weights)
    dists_weighted <- igraph::distances(g_test,to = base::as.character(final_score$Gene), v = base::as.character(base::unique(base::c(active_receptors,inactive_receptors))), mode = "out", weights = NULL, algorithm = "johnson")
    dists_weighted <- reshape2::melt(dists_weighted)
    dists_weighted <- dists_weighted[base::which(!base::is.infinite(dists_weighted$value)),]
    
    recs.to.perturb <- base::unique(base::c(active_receptors,inactive_receptors))
    
    base::cat("   Calculating shortest path weights .\n")
    
    path.sums <- path.prob.sum(subg = subg,iTF = base::unique(iTF.target.info$Gene), Receptors = .pck_env$Receptors,prob.matrix = prob.matrix,ncores = ncores)
    
    OUTPUT <- base::list("active"= active_receptors, "inactive"= inactive_receptors,"iTF.targets" = iTF.target.info,"final.score" = base::as.data.frame(final_score,stringsAsFactors = F),"path.sums" = path.sums, "rec.hotspot" = dists, "rec.hotspot.weighted" = dists_weighted)
    return(OUTPUT)
  }
}

#' Input Data to Markov Chain
#'
#' @description The function computes a Markov chain of the signaling interactome.
#' @param input_data Data frame of input data. Rows correspond to genes, columns to cells
#' @param cutoff Maximum number of zero-value genes, above this cutoff the genes are excluded
#' @param species The species from which the input data has been obtained (UNUSED)
#' @return Markov chain of the signaling interactome
#' @export
Data_preprocessing <- function(input_data,cutoff,species = NULL){
  
  b=input_data
  ##COMVERT CELLS WITH LESS THAT 1 FPKM TO 0
  b[b < 1] <- 0
  ##Add a new gene Dummy with expression value 1, Only works if Dummy is present in the initial network
  b[base::nrow(b)+1, ] <- base::c(.pck_env$dummy.var, base::rep(1,(base::ncol(b)-1)))
  #This is to convert chr into numericversion
  b[2:base::ncol(b)]<-base::as.data.frame(base::lapply(b[2:base::ncol(b)],as.numeric,b[2:base::ncol(b)]))
  ##Renaming the first column for finding the union
  base::colnames(b)[1] <- "Source"
  a=.pck_env$Background_signaling_interactome
  
  ## Removing vertices connected to Dummy which are not receptors
  non.recs <- base::which(a$Source == .pck_env$dummy.var & !(a$Target %in% .pck_env$Receptors))
  recs <- dplyr::setdiff(base::c(1:base::nrow(a)),non.recs)
  a <- a[recs,]
  
  
  ab=plyr::join(a,b,by=c("Source"),type="left",match="first")
  base::colnames(ab)[3:base::ncol(ab)]="source"
  base::colnames(b)[1] <- "Target"
  ab1=plyr::join(a,b,by=c("Target"),type="left",match="first")
  base::names(ab1) <- NULL
  base::names(ab) <- NULL
  ab=ab[,4:base::ncol(ab)]
  ab1=ab1[,4:base::ncol(ab1)]
  ab=base::as.matrix(ab)
  ab1=base::as.matrix(ab1)
  ########Elementwise product
  g=ab * ab1
  ########Sum of elementwise product
  sum_product_expression=base::rowSums(g, na.rm = FALSE, dims = 1)
  g3=base::cbind(a,sum_product_expression)
  ########Calculation of precentage of cells expressed
  h=base::rowSums(g != 0)
  percent_expressed=(h*100)/base::ncol(ab)
  g3=base::cbind(g3,percent_expressed)
  g3[base::is.na(g3)]<-0

  ######NETWORK preprocessing
  g <- igraph::graph.data.frame(base::as.data.frame(g3))
  del=igraph::E(g)[sum_product_expression==0|percent_expressed<base::as.numeric(cutoff)]
  g <- igraph::delete.edges(g,del)
  #SINCE THE TFs AND RECEPTORS ARE ALREADY CONNECTED TO DUMMY, REMOVE ANY NODE THAT HAS ZERO in degree or zero out degree
  #To ensure reachability for the Markov chain
  igraph::V(g)$degree=igraph::degree(g, v=igraph::V(g), mode = c("in"))
  #Select Nodes to be deleted
  del=igraph::V(g)[igraph::V(g)$degree==0]
  #delete vertices from graph
  while(base::length(del)!=0)
  {
    g <- igraph::delete.vertices(g,del)
    igraph::V(g)$degree=igraph::degree(g, v=igraph::V(g), mode = c("in"))
    del=igraph::V(g)[igraph::V(g)$degree==0]
  }
  #Same as above but remove nodes with with zero out degree
  igraph::V(g)$degree=igraph::degree(g, v=igraph::V(g), mode = c("out"))
  #Select Nodes to be deleted
  del=igraph::V(g)[igraph::V(g)$degree==0]
  while(base::length(del)!=0)
  {
    g <- igraph::delete.vertices(g,del)
    igraph::V(g)$degree=igraph::degree(g, v=igraph::V(g), mode = c("out"))
    del=igraph::V(g)[igraph::V(g)$degree==0]
  }
  #####TO EXTRACT THE LARGEST STRONGLY CONNECTED COMPONENT
  members <- igraph::membership(igraph::clusters(g, mode="strong"))
  SCC <- igraph::clusters(g, mode="strong")
  subg <- igraph::induced.subgraph(g, base::which(igraph::membership(SCC) == base::which.max(igraph::sizes(SCC))))
  subg=igraph::simplify(subg,edge.attr.comb=base::list("first"))
  subg
}

#' Obtain Stationary Distribution of Markov Chain
#'
#' @description The function computes the stationary distribution of a given Markov chain
#' @param subg An igraph object representing a Markov Chain. Must be ergodic.
#' @return Stationary distribution of Markov chain
#' @export
Markov_chain_stationary_distribution <- function(subg){ #The function takes the edgelist with probabilitys and computes the SS probability
  ####Write this subgraph as edgelist to make normal graph with ids ie.e names=F
  out <- base::list()
  transition_probability=base::as.data.frame(igraph::as_edgelist(subg,names = F),stringsAsFactors = FALSE)
  transition_probability$probability=base::paste(igraph::E(subg)$sum_product_expression)
  transition_probability[3]=base::as.numeric(transition_probability[[3]])
  myMatrix = Matrix::sparseMatrix(i = transition_probability[1:base::nrow(transition_probability),1], j = transition_probability[1:base::nrow(transition_probability),2],x = transition_probability[1:base::nrow(transition_probability),3])
  #Making a stochastic matrix
  myMatrix = (myMatrix)/Matrix::rowSums((myMatrix))
  el=RSpectra::eigs(Matrix::t(myMatrix),1,which="LR")
  SD=(base::abs(el$vectors))/base::sum(base::abs(el$vectors))
  SD=base::as.data.frame(SD,stringsAsFactors=FALSE)
  SD
  SD=base::cbind((base::as.data.frame(igraph::V(subg)$name,stringsAsFactors = FALSE)),SD)
  SD=base::as.data.frame(SD,stringsAsFactors=FALSE)
  base::colnames(SD)[1] <- "Gene"
  out_SS=base::paste("Steady_state",sep="")
  base::colnames(SD)[2] <- out_SS
  out$SD <- SD
  out$prob.matrix <- myMatrix
  return(out)
}

#----------------------------------------------------------------------------------
# TO CALCULATE THE COMPATABILITY SCORE FOR THE HIGH PROBABILITY INTERMEDIATES
#----------------------------------------------------------------------------------
#' Obtain high probability intermediates
#'
#' @description The function computes high probability intermediates
#' @param x A data frame of gene expression values. Has to contain a column called Gene
#' @param intermediates A one-column data frame of intermediate signaling molecules.
#' @param percentile The percentile above which an intermediate molecule is considered significant (range: 0-100)
#' @return A vector of high probability intermediates
#' @export
high_probability_intermediates <- function(x, intermediates, percentile)
{
  intermediates=base::unique(plyr::join(intermediates,x,by=c("Gene"),type="inner"))
  if(base::nrow(intermediates) == 0){
    return(c())
  }
  #Selecting top 90 percentile of intermediates with steady state probabilities
  percentile=base::as.numeric(percentile)/100
  SS_90percentile=base::as.vector(stats::quantile(intermediates[,2], base::as.numeric(percentile)))
  ##Shortlisting high-probability intermediates > 90 percentile of SS probability
  int=as.vector((base::subset(intermediates, intermediates[2] > SS_90percentile , select=c("Gene")))$Gene)
  int
}


#' Integrate Signaling and Transcriptonal Regulation
#'
#' @description The function integrates downstream TF targets in the signaling network
#' @param g Markov chain of a signaling interactome
#' @param x Steady state distribution of the Markov chain underlying graph g
#' @param deg Data frame of differentially expressed genes (one column only)
#' @param non_interface_TFs TFs that are not regulated by any active signaling pathway
#' @param TF_TF_interactions Regulatory interactions between transcription factors
#' @return Markov chain of the signaling interactome including downstream TF targets
#' @export
integrate_sig_TF <- function(g,x,deg, non_interface_TFs, TF_TF_interactions ){
  
  el=igraph::as_edgelist(g)
  graph_markov=base::as.data.frame(base::cbind(el,igraph::E(g)$Effect))
  base::colnames(graph_markov)=c("Source","Target","Effect")
  base::colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=plyr::join(deg,non_interface_TFs,by=c("Gene"),type="inner")
  #Get the phenotype from the non_interface_DE_TFs and map it to the TF-TF interaction network
  base::colnames(non_interface_DE_TFs)[1]="Target"
  DE_TF_TF_interactions_target=plyr::join(TF_TF_interactions,non_interface_DE_TFs,by=c("Target"),type="left")
  DE_TF_TF_interactions_target=stats::na.omit(DE_TF_TF_interactions_target)
  ab=DE_TF_TF_interactions_target
  base::names(ab)<-NULL
  ab=base::as.matrix(ab)
  ab[,3]=base::as.numeric(ab[,3])*base::as.numeric(ab[,4])
  ab=base::as.data.frame(ab)
  base::names(ab)=c("Source","Target","Effect","DEG")
  graph_markov$Effect=base::as.numeric(base::as.character(graph_markov$Effect))
  graph_markov=base::rbind(graph_markov,ab[1:3]) #merging the nTF interaction with appropriate sign Effect with the original graph
  base::colnames(x)[1] <- "Source"
  ab=plyr::join(graph_markov,x,by=c("Source"),type="left",match="first")
  base::colnames(ab)[3:base::ncol(ab)]="source"
  base::colnames(x)[1] <- "Target"
  ab1=plyr::join(graph_markov,x,by=c("Target"),type="left",match="first")
  base::names(ab1) <- NULL
  base::names(ab) <- NULL
  ab=ab[,4:base::ncol(ab)]
  ab1=ab1[,4:base::ncol(ab1)]
  #creating node SS as the edge property
  weight=base::as.numeric(base::as.matrix(ab))
  graph_markov=(base::cbind(graph_markov,weight))
  graph_markov[base::is.na(graph_markov)] <- 1  #Making TF-TF interactions dependent only on the expression status
  graph_markov$Effect=base::as.numeric(base::as.matrix((graph_markov$Effect)))
  g3 <- igraph::graph.data.frame(base::as.data.frame(graph_markov))
  #updating the graph attribute for the adjacency matrix i.e. product SS (weight) and effect
  igraph::E(g3)$weight=igraph::E(g3)$weight*igraph::E(g3)$Effect
  #deleting TF nodes with no indegree
  igraph::V(g3)$degree=igraph::degree(g3, v=igraph::V(g3), mode = base::c("in"))
  #Select Nodes to be deleted
  del=igraph::V(g3)[igraph::V(g3)$degree==0]
  #delete vertices from graph
  while(base::length(del)!=0)
  {
    g3 <- igraph::delete.vertices(g3,del)
    igraph::V(g3)$degree=igraph::degree(g3, v=igraph::V(g3), mode = base::c("in"))
    del=igraph::V(g3)[igraph::V(g3)$degree==0]
  }
  return(g3)
}

#Function for classifying nTF as up or down regulated
#' Classify non-interface TFs as up or down regulated
#'
#' @description Returns all non-interface TFs as up or down regulated
#' @param g Markov chain of a signaling interactome
#' @param deg Data frame of differentially expressed genes (one column only)
#' @param non_interface_TFs TFs that are not regulated by any active signaling pathway
#' @return Vector of differentially expressed non-interface TFs
#' @export
up_down_tfs <- function(g,deg,non_interface_TFs)
{
  base::colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=plyr::join(deg,non_interface_TFs,by=c("Gene"),type="inner") #IF THIS IS ZERO NEED TO ABORT
  nTF=non_interface_DE_TFs[1]
  base::names(nTF)=NULL
  nTF=base::as.vector(base::t(nTF))
  nTF<-dplyr::intersect(nTF,igraph::V(g)$name) #Some TFs must still be missing in the final g3
  nTF
}

#' Calculating path weights
#'
#' @description Returns a weight for a given path calculated as the product of edge weigths
#' @param path A path as defined by igraph
#' @param graph An igraph graph object
#' @return Numerical path weight
product_path_weight <- function(path, graph){
  edge_weights=igraph::E(graph, path=path)$weight
  v_Edge_weight=base::sum(edge_weights)
  if (v_Edge_weight == 0) {
    return(0)
  } else{
    return(base::prod(edge_weights[edge_weights!=0]))
  }
}

##Function for Compatability score
##This function takes s="source" (one source gene), t="target" (a vector of target gene),g="graph", l="adjacency matrix" as input and finds the shortest paths and passes the argument to spsplit.
##Then it gets the product of intermediates of each path for a source and all its targets and returs its product as final output.
#' Calculating compatibility scores
#'
#' @description Calculates the compatibility between a source node and a list of target nodes depending
#' on their expression values
#' @param s A source node name present in graph g
#' @param t A character vector of target node names present in graph g
#' @param g An igraph graph object
#' @return Compatibility score
#' @export
spcal_path_weight <- function(s,t,g){
  paths=(igraph::get.all.shortest.paths(g, s, t, mode = c("out"), weights=NA)$res)
  if (base::length(paths) == 0){
    return(0)
  } 
  
  weight=base::lapply(paths,product_path_weight,g)
  #s=skewness(unlist(weight))
  s=weight_probability(base::unlist(weight))
  #s=sum(unlist(weight))
  return(s)
}

#' Calculating the weight probability
#'
#' @description Calculates the actual compatibility score given a vector of weights
#' @param x Vector of weights across a path in the graph
#' @return Compatibility score
#' @export
weight_probability <- function(x)
{
  x=base::unlist(x)
  x_pos=x[x>0]
  x_neg=x[x<0]
  x_tot=base::sum(base::abs(x_pos),base::abs(x_neg))
  #probability of the intermediate to be compatible: closer to 1 more compatible it is and closer to zero more incompatible it is
  p_pos=base::sum(x_pos)/x_tot
}

#' Calculating compatibility scores
#'
#' @description Calculates the compatibility between a list of source nodes and a list of target nodes depending
#' on their expression values
#' @param s A vector of source node names present in graph g
#' @param t A character vector of target node names present in graph g
#' @param g An igraph graph object
#' @return Vector of compatibility scores for each source node s
#' @export
comp_score_tf <- function(t,s,g)
{
  comp_score=base::lapply(s,spcal_path_weight,t,g)
}

#' Calculating activation probabilities
#'
#' @description Calculates the probabilities that signaling molecules are active
#' @param x A list of compatibility scores
#' @param y The steady state distribution of a Markov chain
#' @param int Character vector of intermediate signaling molecules to consider
#' @return Data frame of activation probabilities per gene
#' @export
compatability_score <- function(x,y,int)
{
  x=base::unlist(x)
  x=base::cbind(base::as.data.frame(int),base::as.data.frame(base::unlist(x)))
  base::colnames(x)=c("Gene", "Activation_probability")
  x=plyr::join(x,y,by=c("Gene"),type="inner")
  x=x[base::order(x$Activation_probability,decreasing = TRUE),]
}

#' Calculating significance of receptor-TF paths
#'
#' @description Calculates the z-score for each path between receptors and downstream TFs. Uses
#' all paths as background
#' @param subg An igraph graph object
#' @param iTF A vector of interface TF names
#' @param Receptors A vector of receptor names
#' @param prob.matrix A matrix of probabilities for the interaction between genes i (rows) to j (columns)
#' @param ncores The number of cores used (default: 4)
#' @return Data frame of z-scores per path
#' @export
path.prob.sum <- function(subg,iTF,Receptors,prob.matrix,ncores=4){

  Receptors <- dplyr::intersect(Receptors,igraph::V(subg)$name)
  iTF <- dplyr::intersect(iTF,igraph::V(subg)$name)
  
  if(base::length(iTF) == 0){
    base::cat("  No interface TFs in graph.\n")
    return(NULL)
  }else{
    
    path.sum.frame <- base::do.call(base::rbind,parallel::mclapply(X = Receptors,FUN = function(rec){
      paths <- igraph::all_shortest_paths(graph = subg,from = rec, to = iTF,mode = "out",weights = NA)$res
      path.weights <- base::lapply(X = paths,FUN = function(path){
        path <- base::names(path)
        prod = 1
        if(base::length(path) > 1){
          for (j in 1:(base::length(path)-1)) {
            weight <- base::as.numeric(prob.matrix[base::which(prob.matrix$a == path[j] & prob.matrix$b == path[j+1]),][,3])
            prod = prod * weight
          }
        }
        base::names(prod) <- path[base::length(path)]
        return(prod)
      })
      base::names(path.weights) <- base::sapply(path.weights,function(x){names(x)})
      path.weights.sum <- stats::aggregate(base::unlist(path.weights),by = base::list(Gene = base::names(path.weights)), FUN = function(x){base::sum(x,na.rm = TRUE)}, simplify = FALSE, drop = FALSE)
      path.weights.sum$Receptor <- base::rep(rec,base::nrow(path.weights.sum))
      path.weights.sum <- path.weights.sum[,c("Receptor","Gene","x")]
      base::colnames(path.weights.sum)[3] <- "path.weight.sum"
      path.weights.sum$path.weight.sum <- base::as.numeric(path.weights.sum$path.weight.sum)
      return(path.weights.sum)
    },mc.cores = ncores))
    vals  <- path.sum.frame$path.weight.sum
    path.sum.frame$z.score <- (vals-base::mean(vals))/stats::sd(vals)
    # path.sum.frame <- path.sum.frame[which(path.sum.frame$z.score > 0),]
    return(path.sum.frame)
  }
}

#' Scoring ligand-receptor interactions
#'
#' @description Calculates the probability of interaction for each potential ligand-receptor association.
#' @param data A data frame of gene expression values
#' @param tissue.R.TF Data frame of significant receptor-TF associations.
#' @param tissue.LR Data frame of potential ligand-receptor interactions
#' @param LR Data frame of ALL ligand-receptor interactions
#' @param sig.cutoff Significance cutoff as a percentile (ranges from 0 (most permissive) to 1 (most strict)) (default: 0.9)
#' @param z.score.cutoff An optional z-score cutoff from bootstrapping (DEPRECATED; DO NOT CHANGE)
#' @return Data frame significant ligand-receptor interactions
#' @export
scoringFun <- function(data,tissue.R.TF,tissue.LR,LR,sig.cutoff = 0.9,z.score.cutoff = 0){
  pops <- dplyr::union(tissue.LR$Lig.pop,tissue.LR$Rec.pop)
  out.tissue.lr.scored <- base::list()
  for(rec.pop in pops){
    totalscores <- base::c()
    lr <- base::list()
    for(lig.pop in pops){
      base::print(base::paste0(lig.pop,"_",rec.pop))
      lrsub <- tissue.LR[tissue.LR$Lig.pop == lig.pop & tissue.LR$Rec.pop == rec.pop,]
      if(base::nrow(lrsub) == 0){
        next
      }
      dat_lig <- data[base::unique(lrsub$Ligand),base::which(base::colnames(data) == lig.pop)]
      dat_rec <- data[base::unique(lrsub$Receptor),base::which(base::colnames(data) == rec.pop)]
      lig_means <- base::apply(dat_lig,1,function(x){base::mean(x[x>0])})
      rec_means <- base::apply(dat_rec,1,function(x){base::mean(x[x>0])})
      lrsub$score <- base::apply(lrsub,1,function(x){base::as.numeric(lig_means[x[2]])*base::as.numeric(rec_means[x[4]])})
      lr <- base::append(lr,base::list(lrsub))
      
      test_lig <- data[LR$Ligand,base::which(base::colnames(data) == lig.pop)]
      test_lig_means <- base::apply(test_lig,1,function(x){base::mean(base::as.numeric(x[x>0]))})
      test_rec <- data[LR$Receptor,base::which(base::colnames(data) == rec.pop)]
      test_rec_means <- base::apply(test_rec,1,function(x){base::mean(base::as.numeric(x[x>0]))})
      scores <- base::unname(test_lig_means)*base::unname(test_rec_means)
      scores[base::is.na(scores)] <- 0
      totalscores <- base::c(totalscores,scores)
    }
    e <- stats::ecdf(totalscores)
    lr1 <- base::do.call(base::rbind,lr)
    lr1$significance <- e(lr1$score)
    out.tissue.lr.scored[[base::paste0(rec.pop)]] <- lr1
  }
  out.tissue.lr.scored.joined <- base::do.call(base::rbind,out.tissue.lr.scored)
  
  m <- tissue.R.TF[tissue.R.TF$z.score >= z.score.cutoff,base::c(1,2)]
  m <- m[!base::duplicated(m),]
  mm <- base::merge(out.tissue.lr.scored.joined,m,by.x = base::c("Rec.pop","Receptor"), by.y = base::c("Celltype","Receptor"))
  mm <- mm[mm$significance >= sig.cutoff,]
  mm <- mm[,base::c(3,4,5,2,1,6,7,8)]
  
  return(mm)
}






