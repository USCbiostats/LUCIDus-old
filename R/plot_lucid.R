#' Plot Sankey diagram for integrative clustering
#'
#'\code{plot_lucid} generates a Sankey diagram for the results of integrative clustering based on an \code{IntClust} object.
#' @param x An \code{IntClust} class object
#' @param switch An indicator to do label switching with a descending order in gamma or not, the default is FALSE
#' @param colorScale D3 color scheme for the Sankey diagram
#' @export
#' @import networkD3
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Peng, C., Conti, D.V., Integrative latent cluster assignment using multi-omics data with phenotypic traits (under preparation).
#' @examples
#' # Run the model with covariates in the G->X path
#' IntClusCoFit1 <- est_lucid(G=G1,CoG=CoG,Z=Z1,Y=Y1,K=2,family="binary",Pred=TRUE)
#'
#' # Visualize the results of integrative clustering
#' plot_lucid(IntClusCoFit1)

plot_lucid <- function(x, switch = FALSE, colorScale=default) {

  default <- 'd3.scaleOrdinal() .domain(["0", "1", "2", "3", "4", "5", "6", "7", "8"]) .range(["lightsteelblue", "royalblue", "mediumorchid", "gold", "cyan", "tan", "darkgray", "red", "green"])'

  family <- x$family
  K <- x$K

  name <- c(x$Gnames, paste0("IntClust", 1:x$K), x$Znames, "Outcome")
  nodegroup <- c(rep(2, length(x$Gnames)), rep(3, x$K), rep(4, length(x$Znames)), 5)
  Nodes <- cbind(name, nodegroup)

  source <- c(rep(0:(length(x$Gnames)-1), x$K), rep(length(x$Gnames):(length(x$Gnames)+x$K-1), each=length(x$Znames)), length(x$Gnames):(length(x$Gnames)+x$K-1))
  target <- c(rep(length(x$Gnames):(length(x$Gnames)+x$K-1), each=length(x$Gnames)), rep((length(x$Gnames)+x$K):(length(x$Gnames)+x$K+length(x$Znames)-1), x$K), rep(length(x$Gnames)+x$K+length(x$Znames), x$K))
  orgvalue <- c(c(t(summary_lucid(x, switch = switch)$Beta[,-1])), c(t(summary_lucid(x, switch = switch)$Mu)), summary_lucid(x, switch = switch)[[3]])
  linkgroup <- numeric()
  value <- numeric()
  for (i in 1:length(orgvalue)) {
    linkgroup[i] <- ifelse(orgvalue[i]<0, 0, 1)
  }
  if(family == "binary") {
    for (i in (length(orgvalue)-K+1):length(orgvalue)) {
      if (orgvalue[i]==1) {
        linkgroup[i]=6
      }
      if (orgvalue[i]>1) {
        linkgroup[i]=7
      }
      if (orgvalue[i]<1) {
        linkgroup[i]=8
      }
    }
  }
  for (i in 1:length(orgvalue)) {
    value[i] <- ifelse(orgvalue[i]<0, abs(orgvalue[i]), orgvalue[i])
  }
  Links <- cbind(source, target, value, linkgroup)

  nodes <- as.data.frame(Nodes)
  nodes$nodegroup <- as.factor(nodes$nodegroup)
  links <- as.data.frame(Links)
  links$linkgroup <- as.factor(links$linkgroup)

  SankeyFOCM <- sankeyNetwork(Links = links, Nodes = nodes,
                              Source = "source", Target = "target", Value = "value",
                              NodeID = "name", NodeGroup = "nodegroup",LinkGroup = "linkgroup",
                              fontSize = 10, nodeWidth = 20, colourScale = colorScale)

  return(SankeyFOCM)
}
