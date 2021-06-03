

shortest_path_edges <- function(s,t,g){
  if (base::length(s) == 0){
    #cat("No intermediates found. You may decrease the percentile in order to find intermediates.") # decrease percentile cutoff
    return(NULL)
  } 
  if (base::length(t) == 0){
    #cat("No non-terminal differentially expressed TFs found. You may decrease the cutoff.") # decrease normal cutoff
    return(NULL)
  }
  paths=suppressWarnings((igraph::get.all.shortest.paths(g, s, t, mode = base::c("out"))$res))
  if (base::length(paths) == 0){
    #cat("No shortest path found. You may decrease the cutoff in order to find shortest path.") # decrease normal cutoff
    return(NULL)
  }
  edges=base::lapply(paths,edge_shortest,g)
  return(edges)
}

#function for retaining edges
edge_shortest <- function(path, graph)
{
  edges=igraph::E(graph, path=path)
}

#making the edge weight non-negative and taking the absolute
prun.int.g <- function(g){
  igraph::E(g)$weight=1-(base::abs(igraph::E(g)$weight))
  #function to delete edges from tf to dummy in the network
  del=igraph::incident(g, .pck_env$dummy.var, mode = base::c("in"))
  g <- igraph::delete.edges(g,del)
  g <- igraph::set.vertex.attribute(g,"name",.pck_env$dummy.var,"NICHE")
  return(g)
}

# function for converting the input graph and the ints to a shortestpat network
# @param a,i active and inactive source nodes or int
# @param g input graph on which shortest path network must be inferred i.e. gintg
# @param t terminal nodes
to_sp_net_int <- function(s,g,t,deg,non_interface_TFs){
  #changing the edge attributes of the integrated network
  #g=non_neg_weight(g)
  #removing the edges from dummy to TFs as it affectes the shortest paths
  #del=as_adj_edge_list(g, mode = c("in"))$Dummy
  #g <- delete.edges(g,del)
  #Shortest path edges
  edges_a=base::lapply(s,shortest_path_edges,t,g)
  edges_d=base::lapply("NICHE",shortest_path_edges,base::c(s),g)
  #classifying up and downregulated TFs
  up_t=up_down_tfs(g,deg[deg[2]==1,],non_interface_TFs)
  down_t=up_down_tfs(g,deg[deg[2]==-1,],non_interface_TFs)
  #subnetwork from SP edges
  sp_sub_net=igraph::subgraph.edges(g, base::c(base::unlist(edges_a),base::unlist(edges_d)), delete.vertices = T)
  #sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(s), value="int")
  #sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(up_t), value="upregulated")
  #sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(down_t), value="downregulated")
  return(sp_sub_net)
}
