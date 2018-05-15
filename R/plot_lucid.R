#' Plot Sankey diagram for integrative clustering
#'
#'\code{plot_lucid} generates a Sankey diagram for the results of integrative clustering based on an \code{IntClust} object.
#' @param x An \code{IntClust} class object
#' @param ... Additional graphics parameters
#' @export
#' @import networkD3
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Peng, C., Yang, Z., Conti, D.V.

plot_lucid <- function(x, ...) {

  name <- c(x$Gnames, paste0("IntClust", 1:x$K), x$Znames, "Outcome")
  nodegroup <- c(rep(2, length(x$Gnames)), rep(3, x$K), rep(4, length(x$Znames)), 5)
  Nodes <- cbind(name, nodegroup)

  source <- c(rep(0:(length(x$Gnames)-1), x$K), rep(length(x$Gnames):(length(x$Gnames)+x$K-1), each=length(x$Znames)), length(x$Gnames):(length(x$Gnames)+x$K-1))
  target <- c(rep(length(x$Gnames):(length(x$Gnames)+x$K-1), each=length(x$Gnames)), rep((length(x$Gnames)+x$K):(length(x$Gnames)+x$K+length(x$Znames)-1), x$K), rep(length(x$Gnames)+x$K+length(x$Znames), x$K))
  orgvalue <- c(x$beta[1,-1], x$beta[2,-1], x$mu[1, ], x$mu[2, ], x$gamma)
  linkgroup <- numeric()
  value <- numeric()
  for (i in 1:length(orgvalue)) {
    linkgroup[i] <- ifelse(orgvalue[i]<0, 0, 1)
  }
  for (i in 1:length(orgvalue)) {
    value[i] <- ifelse(orgvalue[i]<0, abs(orgvalue[i]), orgvalue[i])
  }
  Links <- cbind(source, target, value, linkgroup)

  nodes <- as.data.frame(Nodes)
  nodes$nodegroup <- as.factor(nodes$nodegroup)
  links <- as.data.frame(Links)
  links$linkgroup <- as.factor(links$linkgroup)

  my_color <- 'd3.scaleOrdinal() .domain(["0", "1", "2", "3", "4", "5"]) .range(["lightsteelblue", "royalblue", "mediumorchid", "gold", "cyan", "tan"])'

  SankeyFOCM <- sankeyNetwork(Links = links, Nodes = nodes,
                              Source = "source", Target = "target", Value = "value",
                              NodeID = "name", NodeGroup = "nodegroup",LinkGroup = "linkgroup",
                              fontSize = 10, nodeWidth = 20, colourScale = my_color)

  return(SankeyFOCM)
}