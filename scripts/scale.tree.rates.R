#### LOAD PACKAGES ####
library(phytools)
library(ape)


#### FUNCTIONs ####
scaleTreeRates <- function(tree,tip.states,
                           max.ratio,
                           nbins,
                           max.transition = 1,
                           model = "ARD",
                           var.start = F,
                           ...){
  
  # Arguments
  # tree: phylogenetic tree for analysis
  # tip.states: vector of tip states associated with tree (can be unordered)
  # max.ratio: maximum ratio of scalars to one, integer
  # nbins: number of bins above and below one, integer
  # model: model for likelihood calculation, character (any argument that fitMk takes) or matrix
  # var.start: whether or not to increment scalar values at root of tree. If true, all scalar
  #            values will be tested and best tree returned. If false, root scalar is set to one.
  # optional arguments: fixedQ: provide Qmatrix with pre-estimated rates
  
  #Get tip states
  x <- tip.states[order(factor(names(tip.states), levels=tree$tip.label))]
  
  if(hasArg(fixedQ)){
    
    print("using supplied Q-matrix")
    
    Q <- list(...)$fixedQ
  } else {
    
    print("estimating rates")
    
    #fit model and extract rates + index matrix
    fit <- fitMkNew(tree, x, model=model)
    rates <- fit$rates
    
    print("building Q-matrix")
    
    #Build matrix
    Q <- fit$index.matrix
    Q[is.na(Q)] <- 0
    
    for(i in 1:length(rates)){
      
      Q[Q == i] <- rates[i]
        
    }
    
    diag(Q) <- -rowSums(Q)
  }
  
  #Get max deviation from parent
  set.range <- max.transition
  
  #Set bins
  bins <- c(1/seq(from = max.ratio, to = 1,length.out = nbins + 1),
            seq(from = 1, to = max.ratio, length.out = nbins + 1))
  bins <- unique(bins)
  
  #Set scalar
  tree$scalar <- rep(1,times = length(tree$edge.length))
  
  trees <- list()
  
  #Determine whether to increment or start from 1
  if(var.start == F){
    
    start.bins <- 1
    
  } else {
    
    start.bins <- bins
    
  }
  
  #Loop through starting scalars and edges
  for(i in 1:length(start.bins)){
    
    print(paste0("starting scalar: ",start.bins[i]))
    #Get start state
    start.scalar <- start.bins[i]
    
    for(j in 1:length(tree$edge.length)){
      
      print(paste0("starting scalar ", start.bins[i] ,", edge ",j))
      
      #Check if edge is starting edge
      if(tree$edge[j,1] == length(tree$tip.label)+1){
        
        #Set previous state to start state
        prev.scalar <- start.scalar
        
      } else {
        
        #Get parent edge
        prev.edge <- which(tree$edge[,2] == tree$edge[j,1])
        
        #Get parent state
        prev.scalar <- tree$scalar[prev.edge]
        
      }
      
      #Determine range of scalars possible
      if(which(bins == prev.scalar) - set.range <= 0){
        
        #Get number of bins lower
        low.range <- which(bins == prev.scalar) - 1
        
        #Determine possible values
        poss.scalars <- bins[(which(bins == prev.scalar) - low.range):
                               (which(bins == prev.scalar)+set.range)]
        
      } else if(which(bins == prev.scalar) + set.range >= length(bins)){
        
        #Get number of bins higher
        high.range <- length(bins) - which(bins == prev.scalar)
        
        #Determine possible values
        poss.scalars <- bins[(which(bins == prev.scalar) - set.range):
                               (which(bins == prev.scalar)+high.range)]
        
      } else {
        
        poss.scalars <- bins[(which(bins == prev.scalar) - set.range):
                               (which(bins == prev.scalar)+set.range)]
        
      }
      
      #Get temporary tree
      temp.tree <- tree
      
      #Assign previous scalar to temporary tree
      temp.tree$scalar[j] <- prev.scalar
      
      #vector for likelihoods
      liks <- c()
      
      for(k in 1:length(poss.scalars)){
        
        #Set current scalar
        cur.scalar <- poss.scalars[k]
        
        #Reassign temp.tree
        temp.tree <- tree
        
        #Assign current scalar to temporary tree
        temp.tree$scalar[j] <- cur.scalar
        
        temp.tree$edge.length <- temp.tree$edge.length * temp.tree$scalar
        
        liks[k] <- fitMkNew(temp.tree, 
                         x, 
                         model=model, 
                         fixedQ = Q)$logLik
        
      }
      
      #Get best index
      index.best <- which(liks == max(liks))
      
      #check for multiple scalars with max likelihood
      if(length(index.best) > 1){

        index.best = which(poss.scalars == prev.scalar)

      }
      
      #Extract best scalar
      best.scalar <- poss.scalars[index.best]
      
      #Assign best scalar to scalar object in tree
      tree$scalar[j] <- best.scalar
      
    }
    
    #Assign to list of trees
    trees[[i]] <- tree
  }
  
  if(var.start == F){
    
    print(paste0("returning tree with scalar"))
    
    #return tree
    return(trees[[i]])
    
  } else {
    #Calculate final likelihoods for each starting scalar tree
    
    #Vector to store
    final.liks <- c()
    
    for(i in 1:length(trees)){
      
      print(paste0("getting final likelihood for starting scalar ", start.bins[i]))
      
      #Get tree to temporary object
      temp.tree <- trees[[i]]
      
      #Multiply temp.tree edge.lengths by scalar
      temp.tree$edge.length <- temp.tree$edge.length * temp.tree$scalar
      
      #Get logLik
      final.liks[i] <- fitMkNew(temp.tree, 
                             x, 
                             model=model, 
                             fixedQ = Q)$logLik
    }
    
    #Get maximum likelihood
    max.final.lik <- max(final.liks)
    
    #Get index assoc. with maximum likelihood
    max.index <- which(final.liks == max(final.liks))
    
    #Get final start state
    final.start <- bins[max.index]
    
    #Extract final tree
    final.tree <- trees[[max.index]]
    
    print(paste0("the best starting scalar is ", final.start, ", associated logLik is ", max.final.lik))
    
    class(final.tree) <- c("phylo","phylo.scaled")
    
    print(paste0("returning tree with scalar"))
    
    return(final.tree)
  }
}





plotScaledTree <- function (tree, colors = NULL, fsize = 1, ftype = "reg", lwd = 2, 
                            pts = FALSE, node.numbers = FALSE, mar = NULL, add = FALSE, 
                            offset = NULL, direction = "rightwards", type = "phylogram", 
                            setEnv = TRUE, part = 1, xlim = NULL, ylim = NULL, nodes = "intermediate", 
                            tips = NULL, maxY = NULL, hold = TRUE, split.vertical = FALSE, 
                            lend = 2, asp = NA, outline = FALSE, plot = TRUE) 
{
  if (inherits(tree, "multiPhylo")) {
    par(ask = TRUE)
    for (i in 1:length(tree)) plotSimmap(tree[[i]], colors = colors, 
                                         fsize = fsize, ftype = ftype, lwd = lwd, pts = pts, 
                                         node.numbers = node.numbers, mar, add, offset, direction, 
                                         type, setEnv, part, xlim, ylim, nodes, tips, maxY, 
                                         hold, split.vertical, lend, asp, outline, plot)
  }
  else {
    if (!inherits(tree, "phylo")) 
      stop("tree should be object of class \"phylo\"")
    if (is.null(tree$scalar)) 
      stop("tree should contain mapped states on edges.")
    ftype <- which(c("off", "reg", "b", "i", "bi") == ftype) - 
      1
    if (!ftype) 
      fsize = 0
    if (is.null(colors)) {
      st <- sort(unique(unlist(sapply(tree$scalar, names))))
      colors <- palette()[1:length(st)]
      names(colors) <- st
      if (length(st) > 1) {
        cat("no colors provided. using the following legend:\n")
        print(colors)
      }
    }
    tree$tip.label <- gsub("_", " ", tree$tip.label)
    if (is.null(mar)) 
      mar = rep(0.1, 4)
    if (hold) 
      null <- dev.hold()
    if (type == "phylogram") {
      if (direction %in% c("upwards", "downwards")) {
        if (outline) {
          fg <- par()$fg
          par(fg = "transparent")
          black <- colors
          black[] <- fg
          updownPhylogram(tree, colors = black, fsize, 
                          ftype, lwd = lwd + 2, pts, node.numbers, 
                          mar, add, offset, direction, setEnv, xlim, 
                          ylim, nodes, tips, split.vertical, lend, 
                          asp, plot)
          par(fg = fg)
        }
        updownPhylogram(tree, colors, fsize, ftype, 
                        lwd, pts, node.numbers, mar, add = if (outline) 
                          TRUE
                        else add, offset, direction, setEnv, xlim, 
                        ylim, nodes, tips, split.vertical, lend, asp, 
                        plot)
      }
      else {
        if (outline) {
          fg <- par()$fg
          par(fg = "transparent")
          black <- colors
          black[] <- fg
          plotPhylogram.rateMap(tree, colors = black, fsize, 
                        ftype, lwd = lwd + 2, pts, node.numbers, 
                        mar, add, offset, direction, setEnv, xlim, 
                        ylim, nodes, tips, split.vertical, lend, 
                        asp, plot)
          par(fg = fg)
        }
        plotPhylogram.rateMap(tree, colors, fsize, ftype, lwd, 
                      pts, node.numbers, mar, add = if (outline) 
                        TRUE
                      else add, offset, direction, setEnv, xlim, 
                      ylim, nodes, tips, split.vertical, lend, asp, 
                      plot)
      }
    }
    else if (type == "fan") {
      if (outline) {
        fg <- par()$fg
        par(fg = "transparent")
        black <- colors
        black[] <- fg
        plotFan(tree, colors = black, fsize, ftype, 
                lwd = lwd + 2, mar, add, part, setEnv, xlim, 
                ylim, tips, maxY, lend, plot, offset)
        par(fg = fg)
      }
      plotFan(tree, colors, fsize, ftype, lwd, mar, add = if (outline) 
        TRUE
        else add, part, setEnv, xlim, ylim, tips, maxY, 
        lend, plot, offset)
    }
    else if (type == "cladogram") {
      if (outline) {
        fg <- par()$fg
        par(fg = "transparent")
        black <- colors
        black[] <- fg
        plotCladogram(tree, colors = black, fsize, ftype, 
                      lwd = lwd + 2, mar, add, offset, direction, 
                      xlim, ylim, nodes, tips, lend, asp, plot)
        par(fg = fg)
      }
      plotCladogram(tree, colors, fsize, ftype, lwd, mar, 
                    add = if (outline) 
                      TRUE
                    else add, offset, direction, xlim, ylim, nodes, 
                    tips, lend, asp, plot)
    }
    if (hold) 
      null <- dev.flush()
  }
}

     

## functions plot stochastic character mapped trees
## written by Liam Revell 2011-2023

# function to plot simmap tree in type "phylogram"
# written by Liam J. Revell 2011-2023
plotPhylogram.rateMap<-function(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
                        add,offset,direction,setEnv,xlim,ylim,placement,tips,split.vertical,lend,
                        asp,plot,underscore){
  if(split.vertical&&!setEnv){
    cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
    spit.vertical<-FALSE
  }
  # set offset fudge (empirically determined)
  offsetFudge<-1.37
  # reorder
  cw<-reorderRateMap(tree)
  pw<-reorderRateMap(tree,order = "postorder")
  # count nodes and tips
  n<-Ntip(cw)
  m<-cw$Nnode
  # Y coordinates for nodes
  Y<-matrix(NA,m+n,1)
  # first, assign y coordinates to all the tip nodes
  if(is.null(tips)) Y[cw$edge[cw$edge[,2]<=n,2]]<-1:n
  else Y[cw$edge[cw$edge[,2]<=n,2]]<-if(is.null(names(tips))) 
    tips[sapply(1:Ntip(cw),function(x,y) which(y==x),y=cw$edge[cw$edge[,2]<=n,2])]
  else if(!underscore) tips[gsub(" ","_",cw$tip.label)]
  # get Y coordinates of the nodes
  nodes<-unique(pw$edge[,1])
  for(i in 1:m){
    if(placement=="intermediate"){ 
      desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
      Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
    } else if(placement=="centered"){
      desc<-getDescendants(pw,nodes[i])
      desc<-desc[desc<=Ntip(pw)]
      Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
    } else if(placement=="weighted"){
      desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
      n1<-desc[which(Y[desc]==min(Y[desc]))]
      n2<-desc[which(Y[desc]==max(Y[desc]))]
      v1<-pw$edge.length[which(pw$edge[,2]==n1)]
      v2<-pw$edge.length[which(pw$edge[,2]==n2)]
      Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
    } else if(placement=="inner"){
      desc<-getDescendants(pw,nodes[i])
      desc<-desc[desc<=Ntip(pw)]
      mm<-which(abs(Y[desc]-median(Y[1:Ntip(pw)]))==min(abs(Y[desc]-
                                                              median(Y[1:Ntip(pw)]))))
      if(length(mm>1)) mm<-mm[which(Y[desc][mm]==min(Y[desc][mm]))]
      Y[nodes[i]]<-Y[desc][mm]
    }
  }
  # compute node heights
  H<-nodeHeights(cw)
  # open plot
  par(mar=mar)
  if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
  if(!add) plot.new()
  ###
  if(is.null(xlim)){
    pp<-par("pin")[1]
    sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
      offsetFudge*fsize*strwidth("W",units="inches")
    alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    xlim<-if(direction=="leftwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
  }
  if(is.null(ylim)) ylim=range(Y)
  if(direction=="leftwards") H<-max(H)-H
  plot.window(xlim=xlim,ylim=ylim,asp=asp)
  if(plot){
    ####
    if(!split.vertical){
      for(i in 1:m) lines(H[which(cw$edge[,1]==nodes[i]),1],
                          Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],col=colors[which(names(colors) == cw$scalar[match(nodes[i],
                                                                                                     cw$edge[,1])][1])],lwd=lwd,lend=lend)
    }
    for(i in 1:nrow(cw$edge)){
      x<-H[i,1]
        if(direction=="leftwards")
          lines(c(x,x-cw$edge.length[i]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
                col=colors[which(names(colors) == (cw$scalar[i]))],lwd=lwd,lend=lend)
        else lines(c(x,x+cw$edge.length[i]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
                   col=colors[which(names(colors) == (cw$scalar[i]))],lwd=lwd,lend=lend)
        if(pts) points(c(x,x+cw$edge.length[i]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
                       pch=20,lwd=(lwd-1))
        x<-x+if(direction=="leftwards") -cw$scalar[i] else cw$scalar[i]
      }
    if(node.numbers){
      symbols(if(direction=="leftwards") max(H) else 0,
              mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
              rectangles=matrix(c(1.2*fsize*strwidth(as.character(Ntip(cw)+1)),
                                  1.4*fsize*strheight(as.character(Ntip(cw)+1))),1,2),inches=FALSE,
              bg="white",add=TRUE)
      text(if(direction=="leftwards") max(H) else 0,
           mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),Ntip(cw)+1,
           cex=fsize)
      for(i in 1:nrow(cw$edge)){
        x<-H[i,2]
        if(cw$edge[i,2]>Ntip(cw)){
          symbols(x,Y[cw$edge[i,2]],
                  rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),
                                      1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,
                  bg="white",add=TRUE)
          text(x,Y[cw$edge[i,2]],cw$edge[i,2],cex=fsize)
        }
      }
    }
    if(direction=="leftwards") pos<-if(par()$usr[1]>par()$usr[2]) 4 else 2
    if(direction=="rightwards") pos<-if(par()$usr[1]>par()$usr[2]) 2 else 4
    for(i in 1:n) if(ftype) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,
                                 offset=offset,cex=fsize,font=ftype)
  }
  if(setEnv){
    PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
             show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
             font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
             x.lim=xlim,y.lim=ylim,
             direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
             edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
                                    function(x,y,z) y[match(x,z)],y=H,z=cw$edge),yy=Y[,1])
    assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
  }
  if(plot) if(split.vertical) splitEdgeColor(cw,colors,lwd)
}

reorderRateMap <- function (tree, order = "cladewise", index.only = FALSE, ...) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  ii <- reorder.phylo(tree, order, index.only = TRUE, ...)
  if (!index.only) {
    if (inherits(ii, "phylo")) 
      ii <- whichorder(ii$edge[, 2], tree$edge[, 2])
    tree$edge <- tree$edge[ii, ]
    tree$edge.length <- tree$edge.length[ii]
    if (!is.null(tree$scalar)) {
      tree$scalar <- tree$scalar[ii]
    }
    attr(tree, "order") <- order
    return(tree)
  }
  else return(ii)
}

fitMkNew <- function (tree, x, model = "SYM", fixedQ = NULL, ...) 
{
  if (hasArg(output.liks)) 
    output.liks <- list(...)$output.liks
  else output.liks <- FALSE
  if (hasArg(q.init)) 
    q.init <- list(...)$q.init
  else q.init <- length(unique(x))/sum(tree$edge.length)
  if (hasArg(opt.method)) 
    opt.method <- list(...)$opt.method
  else opt.method <- "nlminb"
  if (hasArg(min.q)) 
    min.q <- list(...)$min.q
  else min.q <- 1e-12
  if (hasArg(max.q)) 
    max.q <- list(...)$max.q
  else max.q <- max(nodeHeights(tree)) * 100
  if (hasArg(logscale)) 
    logscale <- list(...)$logscale
  else logscale <- FALSE
  N <- Ntip(tree)
  M <- tree$Nnode
  if (is.matrix(x)) {
    x <- x[tree$tip.label, ]
    m <- ncol(x)
    states <- colnames(x)
  }
  else {
    x <- to.matrix(x, sort(unique(x)))
    m <- ncol(x)
    states <- colnames(x)
  }
  if (hasArg(pi)) 
    pi <- list(...)$pi
  else pi <- "equal"
  if (is.numeric(pi)) 
    root.prior <- "given"
  if (pi[1] == "equal") {
    pi <- setNames(rep(1/m, m), states)
    root.prior <- "flat"
  }
  else if (pi[1] == "estimated") {
    pi <- if (!is.null(fixedQ)) 
      statdist(fixedQ)
    else statdist(summary(fitMk(tree, x, model), quiet = TRUE)$Q)
    cat(paste("Using pi estimated from the stationary", 
              "distribution of Q assuming a flat prior.\npi =\n"))
    print(round(pi, 6))
    cat("\n")
    root.prior <- "stationary"
  }
  else if (pi[1] == "fitzjohn") 
    root.prior <- "nuisance"
  if (is.numeric(pi)) {
    pi <- pi/sum(pi)
    if (is.null(names(pi))) 
      pi <- setNames(pi, states)
    pi <- pi[states]
  }
  if (is.null(fixedQ)) {
    if (is.character(model)) {
      rate <- matrix(NA, m, m)
      if (model == "ER") {
        k <- rate[] <- 1
        diag(rate) <- NA
      }
      else if (model == "ARD") {
        k <- m * (m - 1)
        rate[col(rate) != row(rate)] <- 1:k
      }
      else if (model == "SYM") {
        k <- m * (m - 1)/2
        ii <- col(rate) < row(rate)
        rate[ii] <- 1:k
        rate <- t(rate)
        rate[ii] <- 1:k
      }
    }
    else {
      if (ncol(model) != nrow(model)) 
        stop("model is not a square matrix")
      rate <- model
      m <- ncol(rate)
      states <- as.character(1:ncol(rate))
      k <- max(rate)
      
      if(length(states) != ncol(x)){
        missing <- states[which(!states %in% colnames(x))]
        x <- cbind(x,matrix(data=0,nrow=nrow(x),ncol=length(missing),
                            dimnames = list(rownames(x),
                                            missing)))
        x <- x[,sort(colnames(x))]
        
      }
    }
    Q <- matrix(0, m, m)
  }
  else {
    m <- ncol(fixedQ)
    states <- as.character(1:ncol(fixedQ))
    rate <- matrix(NA, m, m)
    k <- m * (m - 1)
    rate[col(rate) != row(rate)] <- 1:k
    Q <- fixedQ
    if(ncol(fixedQ) != ncol(x)){
      missing <- states[which(!states %in% colnames(x))]
      x <- cbind(x,matrix(data=0,nrow=nrow(x),ncol=length(missing),
                          dimnames = list(rownames(x),
                                          missing)))
    }
  }
  index.matrix <- rate
  tmp <- cbind(1:m, 1:m)
  rate[tmp] <- 0
  rate[rate == 0] <- k + 1
  liks <- rbind(x, matrix(0, M, m, dimnames = list(1:M + N, 
                                                   states)))
  pw <- reorder(tree, "pruningwise")
  lik <- function(Q, output.liks = FALSE, pi, ...) {
    if (hasArg(output.pi)) 
      output.pi <- list(...)$output.pi
    else output.pi <- FALSE
    if (is.Qmatrix(Q)) 
      Q <- unclass(Q)
    if (any(is.nan(Q)) || any(is.infinite(Q))) 
      return(1e+50)
    comp <- vector(length = N + M, mode = "numeric")
    parents <- unique(pw$edge[, 1])
    root <- min(parents)
    for (i in 1:length(parents)) {
      anc <- parents[i]
      ii <- which(pw$edge[, 1] == parents[i])
      desc <- pw$edge[ii, 2]
      el <- pw$edge.length[ii]
      v <- vector(length = length(desc), mode = "list")
      for (j in 1:length(v)) {
        v[[j]] <- EXPM(Q * el[j]) %*% liks[desc[j], 
        ]
      }
      if (anc == root) {
        if (is.numeric(pi)) 
          vv <- Reduce("*", v)[, 1] * pi
        else if (pi[1] == "fitzjohn") {
          D <- Reduce("*", v)[, 1]
          pi <- D/sum(D)
          vv <- D * D/sum(D)
        }
      }
      else vv <- Reduce("*", v)[, 1]
      comp[anc] <- sum(vv)
      liks[anc, ] <- vv/comp[anc]
    }
    if (output.liks) 
      return(liks[1:M + N, , drop = FALSE])
    else if (output.pi) 
      return(pi)
    else {
      logL <- -sum(log(comp[1:M + N]))
      if (is.na(logL)) 
        logL <- Inf
      return(logL)
    }
  }
  if (is.null(fixedQ)) {
    if (length(q.init) != k) 
      q.init <- rep(q.init[1], k)
    q.init <- if (logscale) 
      log(q.init)
    else q.init
    if (opt.method == "optim") {
      fit <- if (logscale) 
        optim(q.init, function(p) lik(makeQ(m, exp(p), 
                                            index.matrix), pi = pi), method = "L-BFGS-B", 
              lower = rep(log(min.q), k), upper = rep(log(max.q), 
                                                      k))
      else optim(q.init, function(p) lik(makeQ(m, p, index.matrix), 
                                         pi = pi), method = "L-BFGS-B", lower = rep(min.q, 
                                                                                    k), upper = rep(max.q, k))
    }
    else if (opt.method == "none") {
      fit <- list(objective = lik(makeQ(m, q.init, index.matrix), 
                                  pi = pi), par = q.init)
    }
    else {
      fit <- if (logscale) 
        nlminb(q.init, function(p) lik(makeQ(m, exp(p), 
                                             index.matrix), pi = pi), lower = rep(log(min.q), 
                                                                                  k), upper = rep(log(max.q), k))
      else nlminb(q.init, function(p) lik(makeQ(m, p, 
                                                index.matrix), pi = pi), lower = rep(0, k), 
                  upper = rep(max.q, k))
    }
    if (logscale) 
      fit$par <- exp(fit$par)
    if (pi[1] == "fitzjohn") 
      pi <- setNames(lik(makeQ(m, fit$par, index.matrix), 
                         FALSE, pi = pi, output.pi = TRUE), states)
    obj <- list(logLik = if (opt.method == "optim") -fit$value else -fit$objective, 
                rates = fit$par, index.matrix = index.matrix, states = states, 
                pi = pi, method = opt.method, root.prior = root.prior)
    if (output.liks) 
      obj$lik.anc <- lik(makeQ(m, obj$rates, index.matrix), 
                         TRUE, pi = pi)
  }
  else {
    fit <- lik(Q, pi = pi)
    if (pi[1] == "fitzjohn") 
      pi <- setNames(lik(Q, FALSE, pi = pi, output.pi = TRUE), 
                     states)
    obj <- list(logLik = -fit, rates = Q[sapply(1:k, function(x, 
                                                              y) which(x == y), index.matrix)], index.matrix = index.matrix, 
                states = states, pi = pi, root.prior = root.prior)
    if (output.liks) 
      obj$lik.anc <- lik(makeQ(m, obj$rates, index.matrix), 
                         TRUE, pi = pi)
  }
  lik.f <- function(q) -lik(q, output.liks = FALSE, pi = if (root.prior == 
                                                             "nuisance") 
    "fitzjohn"
    else pi)
  obj$lik <- lik.f
  class(obj) <- "fitMk"
  return(obj)
}

is.Qmatrix<-function(x) "Qmatrix" %in% class(x)

makeQ<-function(m,q,index.matrix){
  Q<-matrix(0,m,m)
  Q[]<-c(0,q)[index.matrix+1]
  diag(Q)<-0
  diag(Q)<--rowSums(Q)
  Q
}

EXPM<-function(x,...){
  e_x<-if(isSymmetric(x)) matexpo(x) else expm(x,...)
  dimnames(e_x)<-dimnames(x)
  e_x
}


