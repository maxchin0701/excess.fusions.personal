describe.simmap2 <- function (tree, ...) {
                      if (hasArg(plot)) 
                        plot <- list(...)$plot
                      else plot <- FALSE
                      if (hasArg(check.equal)) 
                        check.equal <- list(...)$check.equal
                      else check.equal <- FALSE
                      if (hasArg(message)) 
                        message <- list(...)$message
                      else message <- FALSE
                      if (hasArg(ref.tree)) 
                        ref.tree <- list(...)$ref.tree
                      else ref.tree <- NULL
                      if (inherits(tree, "multiPhylo")) {
                        if (check.equal) {
                          TT <- sapply(tree, function(x, y) sapply(y, all.equal.phylo, 
                                                                   x), y = tree)
                          check <- all(TT)
                          if (!check) 
                            cat("Note: Some trees are not equal.\nA \"reference\" tree will be computed if none was provided.\n\n")
                        }
                        else check <- TRUE
                        YY <- getStates(tree, "both")
                        states <- c()
                        for(i in 1:length(tree)){
                          states <- c(states,colnames(tree[[i]]$mapped.edge))
                        }
                        states <- sort(unique(as.vector(states)))
                        if (is.null(ref.tree) && check) 
                          ZZ <- t(apply(YY, 1, function(x, levels, Nsim) summary(factor(x, 
                                                                                        levels))/Nsim, levels = states, Nsim = length(tree)))
                        else {
                          if (is.null(ref.tree)) {
                            cat("No reference tree provided & some trees are unequal.\nComputing majority-rule consensus tree.\n")
                            ref.tree <- consensus(tree, p = 0.5)
                          }
                          YYp <- matrix(NA, ref.tree$Nnode, length(tree), 
                                        dimnames = list(1:ref.tree$Nnode + Ntip(ref.tree), 
                                                        NULL))
                          for (i in 1:length(tree)) {
                            M <- matchNodes(ref.tree, tree[[i]])
                            jj <- sapply(M[, 2], function(x, y) if (x %in% 
                                                                    y) 
                              which(as.numeric(y) == x)
                              else NA, y = as.numeric(rownames(YY)))
                            YYp[, i] <- YY[jj, i]
                          }
                          ZZ <- t(apply(YYp, 1, function(x, levels) summary(factor(x[!is.na(x)], 
                                                                                   levels))/sum(!is.na(x)), levels = states))
                        }
                        XX <- countSimmap(tree, states, FALSE)
                        XX <- XX[, -(which(as.vector(diag(-1, length(states))) == 
                                             -1) + 1)]
                        AA <- t(sapply(unclass(tree), function(x) c(colSums(x$mapped.edge), 
                                                                    sum(x$edge.length))))
                        colnames(AA)[ncol(AA)] <- "total"
                        BB <- getStates(tree, type = "tips")
                        CC <- t(apply(BB, 1, function(x, levels, Nsim) summary(factor(x, 
                                                                                      levels))/Nsim, levels = states, Nsim = length(tree)))
                        x <- list(count = XX, times = AA, ace = ZZ, tips = CC, 
                                  tree = tree, ref.tree = if (!is.null(ref.tree)) ref.tree else NULL)
                        class(x) <- "describe.simmap"
                      }
                      else if (inherits(tree, "phylo")) {
                        XX <- countSimmap(tree, message = FALSE)
                        YY <- getStates(tree)
                        states <- sort(unique(YY))
                        AA <- setNames(c(colSums(tree$mapped.edge), sum(tree$edge.length)), 
                                       c(colnames(tree$mapped.edge), "total"))
                        AA <- rbind(AA, AA/AA[length(AA)])
                        rownames(AA) <- c("raw", "prop")
                        x <- list(N = XX$N, Tr = XX$Tr, times = AA, states = YY, 
                                  tree = tree)
                        class(x) <- "describe.simmap"
                      }
                      if (message) 
                        print(x)
                      if (plot) 
                        plot(x)
                      x
                    }
