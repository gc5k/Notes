SingleLocusTree <- function(n = 10, N = 1000, debug = FALSE) {
  
  # Length of each branch in the tree
  lengths <- rep(0, n)
  # Number of descendants for each node
  descendants <- rep(1, n)
  # Parent node for each node
  parent <- rep(0, n)
  # Descendant nodes for each node
  activeNodes <- 1:n
  nextNode = n + 1
  
  for (i in n:2)
  {
    # sample time to next coalescent event
    t <- rexp(1, rate = (i*(i-1)/2)/N)
    # all current nodes increase in length by time t
    lengths[activeNodes] <- lengths[activeNodes] + t
    # next we select two nodes to coalesce at random
    coalescence <- sample(activeNodes, 2)
    
    if (debug)
      cat("After ", t, "generations: Nodes ", coalescence[1], " and ", coalescence[2], "coalesce\n")
    
    # assign their parental nodes
    parent[coalescence] = nextNode
    
    # remove these from the list of active nodes
    activeNodes <- setdiff(activeNodes, coalescence)
    # create a new node
    parent[nextNode] <- 0
    lengths[nextNode] <- 0
    descendants[nextNode] <- descendants[coalescence[1]] + descendants[coalescence[2]]
    # And update list of active nodes
    activeNodes <- c(activeNodes, nextNode)
    nextNode <- nextNode + 1
  }
  list(descendants = descendants, parent = parent, lengths = lengths)
}
