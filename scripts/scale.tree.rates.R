#### LOAD PACKAGES ####
library(phytools)
library(ape)


#### FUNCTIONs ####
scaleTreeRates <- function(tree,tip.states,max.ratios,nbins,max.transition,
                           type="discrete",model="ARD",use.eigen=F,
                           use.expm=F,mp=1e+50,print.rates=F,
                           ...){
  
  
  #Get max deviation from parent
  set.range <- max.transition
  
  #Set bins
  bins <- c(1/seq(from = max.ratios[1], to = 1,length.out = nbins[1] + 1),
            seq(from = 1, to = max.ratios[2], length.out = nbins[2] + 1))
  bins <- unique(bins)
  
  #Set scalar
  tree$scalar <- rep(1,times = length(tree$edge.length))
  
  #Set tip.states
  
  #get rates
  if(hasArg(rates)){
    
    rates <- list(...)$rates
  }
  
  x <- tip.states
  
  trees <- list()
  
  #Loop through starting scalars and edges
  for(i in 1:length(bins)){
    
    print(paste0("starting scalar: ",bins[i]))
    #Get start state
    start.scalar <- bins[i]
    
    for(j in 1:length(tree$edge.length)){
      
      print(paste0("starting scalar ", bins[i] ,", edge ",j))
      
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
      
      temp.tree$edge.length <- temp.tree$edge.length * temp.tree$scalar
      
      #Fit our unscaled model
      if(hasArg(rates) == F){
      
        fit.unscaled <- ace2(x=x, phy=temp.tree, type=type, model=model, 
                            use.eigen=use.eigen,use.expm=use.expm,mp=mp,
                            print.rates=print.rates)
      } else {
        
        fit.unscaled <- ace2(x=x, phy=temp.tree, type=type, model=model, 
                             use.eigen=use.eigen,use.expm=use.expm,mp=mp,
                             print.rates=print.rates,rates=rates)
      }
      
      if(type == "discrete"){
        #Save likelihood
        lik.best <- fit.unscaled$loglik
      } else {
        lik.best <- fit.unscaled$resloglik
      }
      
      for(k in 1:length(poss.scalars)){
        
        #Set current scalar
        cur.scalar <- poss.scalars[k]
        
        #Reassign temp.tree
        temp.tree <- tree
        
        #Assign current scalar to temporary tree
        temp.tree$scalar[j] <- cur.scalar
        
        temp.tree$edge.length <- temp.tree$edge.length * temp.tree$scalar
        
        if(hasArg(rates) == F){
        
          #Fit our scaled model
          fit.scaled <- ace2(x=x, phy=temp.tree, type=type, model=model, 
                             use.eigen=use.eigen,use.expm=use.expm,mp=mp,
                             print.rates=print.rates)
        } else {
          
          fit.scaled <- ace2(x=x, phy=temp.tree, type=type, model=model, 
                             use.eigen=use.eigen,use.expm=use.expm,mp=mp,
                             print.rates=print.rates,rates = rates)
          
        }
        
        #Extract loglik
        if(type == "discrete"){
          #Save likelihood
          lik.scaled <- fit.scaled$loglik
        } else {
          lik.scaled <- fit.scaled$resloglik
        }
        
        #Check if new likelihood is better than current best
        if(lik.scaled >= lik.best){
          
          lik.best <- lik.scaled
          
          index.best <- k
          
        }
      }
      
      #Extract best scalar
      best.scalar <- poss.scalars[index.best]
      
      #Assign best scalar to scalar object in tree
      tree$scalar[j] <- best.scalar
      
    }
    
    #Assign to list of trees
    trees[[i]] <- tree
  }
  
  #Calculate final likelihoods for each starting scalar tree
  
  #Vector to store
  final.liks <- c()
  
  for(i in 1:length(trees)){
    
    print(paste0("getting final likelihood for starting scalar ", bins[i]))
    
    #Get tree to temporary object
    temp.tree <- trees[[i]]
    
    #Multiply temp.tree edge.lengths by scalar
    temp.tree$edge.length <- temp.tree$edge.length * temp.tree$scalar
    
    if(hasArg(rates)){
    
      #Get logLik
      fit <- ace2(x=x, phy=temp.tree, type=type, model=model, 
                  use.eigen=use.eigen,use.expm=use.expm,mp=mp,
                  print.rates=print.rates)
    } else {
      
      #Get logLik
      fit <- ace2(x=x, phy=temp.tree, type=type, model=model, 
                  use.eigen=use.eigen,use.expm=use.expm,mp=mp,
                  print.rates=print.rates,rates=rates)
    }
    
    if(type == "discrete"){
      final.liks[i] <- fit$loglik
    } else {
      final.liks[i] <- fit$resloglik
      
    }
    
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


## function to plot simmap tree in type "phylogram"
## written by Liam J. Revell 2011-2023
updownPhylogram<-function(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
                          add,offset,direction,setEnv,xlim,ylim,placement,tips,split.vertical,lend,
                          asp,plot,underscore){
  if(split.vertical&&!setEnv){
    cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
    spit.vertical<-FALSE
  }
  # set offset fudge (empirically determined)
  offsetFudge<-1.37
  # reorder
  cw<-reorderSimmap(tree)
  pw<-reorderSimmap(tree,"postorder")
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
  if(is.null(ylim)){
    pp<-par("pin")[2]
    sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
      offsetFudge*fsize*strwidth("W",units="inches")
    alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    ylim<-if(direction=="downwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
  }
  if(is.null(xlim)) xlim=range(Y)
  if(direction=="downwards") H<-max(H)-H
  plot.window(xlim=xlim,ylim=ylim,asp=asp)
  ####
  if(plot){
    if(!split.vertical){
      for(i in 1:m) lines(Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],
                          H[which(cw$edge[,1]==nodes[i]),1],
                          col=colors[names(cw$maps[[match(nodes[i],
                                                          cw$edge[,1])]])[1]],lwd=lwd,lend=lend)
    }
    for(i in 1:nrow(cw$edge)){
      x<-H[i,1]
      for(j in 1:length(cw$maps[[i]])){
        if(direction=="downwards")
          lines(c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),c(x,x-cw$maps[[i]][j]),
                col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
        else lines(c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),c(x,x+cw$maps[[i]][j]),
                   col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
        if(pts) points(c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),c(x,x+cw$maps[[i]][j]),
                       pch=20,lwd=(lwd-1))
        x<-x+if(direction=="downwards") -cw$maps[[i]][j] else cw$maps[[i]][j]
        j<-j+1
      }
    }
    if(node.numbers){
      symbols(mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
              if(direction=="downwards") max(H) else 0,
              rectangles=matrix(c(1.2*fsize*strwidth(as.character(Ntip(cw)+1)),
                                  1.4*fsize*strheight(as.character(Ntip(cw)+1))),1,2),inches=FALSE,
              bg="white",add=TRUE)
      text(mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
           if(direction=="downwards") max(H) else 0,Ntip(cw)+1,
           cex=fsize)
      for(i in 1:nrow(cw$edge)){
        x<-H[i,2]
        if(cw$edge[i,2]>Ntip(cw)){
          symbols(Y[cw$edge[i,2]],x,
                  rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),
                                      1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,
                  bg="white",add=TRUE)
          text(Y[cw$edge[i,2]],x,cw$edge[i,2],cex=fsize)
        }
      }
    }
    if(direction=="downwards") pos<-if(par()$usr[3]>par()$usr[4]) 2 else 4
    if(direction=="upwards") pos<-if(par()$usr[3]>par()$usr[4]) 2 else 4
    for(i in 1:n){
      shift<-offset*fsize*strwidth("W")*(diff(par()$usr[3:4])/diff(par()$usr[1:2]))
      if((direction=="downwards"&&diff(par()$usr[3:4])>0) ||
         (direction=="upwards"&&diff(par()$usr[3:4])<0)) shift<--shift
      if(ftype){
        text(labels=cw$tip.label[i],Y[i],
             H[which(cw$edge[,2]==i),2]+shift,
             pos=pos,offset=0,cex=fsize,font=ftype,
             srt=if(direction=="downwards") 270 else 90)
      }
    }
  }
  if(setEnv){
    PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
             show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
             font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
             x.lim=xlim,y.lim=ylim,
             direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
             edge=cw$edge,xx=Y[,1],yy=sapply(1:(Ntip(cw)+cw$Nnode),
                                             function(x,y,z) y[match(x,z)],y=H,z=cw$edge))
    assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
  }
  if(plot) if(split.vertical) splitEdgeColor(cw,colors,lwd)
}

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

# function to plot simmap tree in type "fan"
# written by Liam J. Revell 2013-2017
plotFan<-function(tree,colors,fsize,ftype,lwd,mar,add,part,setEnv,xlim,ylim,tips,maxY,lend,plot,offset){
  if(!plot) cat("plot=FALSE option is not permitted for type=\"fan\". Tree will be plotted.\n")
  if(is.null(offset)) offset<-1
  # reorder
  cw<-reorder(tree)
  pw<-reorder(tree,"pruningwise")
  # count nodes and tips
  n<-Ntip(cw)
  m<-cw$Nnode 
  # get Y coordinates on uncurved space
  Y<-vector(length=m+n)
  if(is.null(tips)) tips<-1:n
  if(part<1.0) Y[cw$edge[cw$edge[,2]<=n,2]]<-0:(n-1)
  else Y[cw$edge[cw$edge[,2]<=n,2]]<-tips
  nodes<-unique(pw$edge[,1])
  for(i in 1:m){
    desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
    Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
  }
  if(is.null(maxY)) maxY<-max(Y)
  Y<-setNames(Y/maxY*2*pi,1:(n+m))
  Y<-part*cbind(Y[as.character(cw$edge[,2])],Y[as.character(cw$edge[,2])])
  R<-nodeHeights(cw)
  # now put into a circular coordinate system
  x<-R*cos(Y)
  y<-R*sin(Y)
  # optimize x & y limits
  par(mar=mar)
  offsetFudge<-1.37 # empirically determined
  OFFSET<-0
  pp<-par("pin")[1]
  sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
    offsetFudge*OFFSET*fsize*strwidth("W",units="inches") 
  alp<-optimize(function(a,H,sw,pp) (2*a*1.04*max(H)+2*sw-pp)^2,H=R,sw=sw,pp=pp,
                interval=c(0,1e6))$minimum
  if(part<=0.25) x.lim<-y.lim<-c(0,max(R)+sw/alp)
  else if(part>0.25&&part<=0.5){ 
    x.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
    y.lim<-c(0,max(R)+sw/alp)
  } else x.lim<-y.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
  if(is.null(xlim)) xlim<-x.lim
  if(is.null(ylim)) ylim<-y.lim
  # plot tree
  if(!add) plot.new()
  plot.window(xlim=xlim,ylim=ylim,asp=1)
  # plot radial lines (edges)
  ## first, the lines emerging from the root (if there are only two):
  jj<-which(cw$edge[,1]==(Ntip(cw)+1))
  if(length(jj)==2){
    m.left<-cumsum(cw$maps[[jj[1]]])/sum(cw$maps[[jj[1]]])
    xx.left<-c(x[jj[1],1],x[jj[1],1]+(x[jj[1],2]-x[jj[1],1])*m.left)
    yy.left<-c(y[jj[1],1],y[jj[1],1]+(y[jj[1],2]-y[jj[1],1])*m.left)
    m.right<-cumsum(cw$maps[[jj[2]]])/sum(cw$maps[[jj[2]]])
    xx.right<-c(x[jj[2],1],x[jj[2],1]+(x[jj[2],2]-x[jj[2],1])*m.right)
    yy.right<-c(y[jj[2],1],y[jj[2],1]+(y[jj[2],2]-y[jj[2],1])*m.right)
    xx<-c(xx.left[length(xx.left):1],xx.right[2:length(xx.right)])
    yy<-c(yy.left[length(yy.left):1],yy.right[2:length(yy.right)])
    col<-colors[c(names(m.left)[length(m.left):1],names(m.right))]
    segments(xx[2:length(xx)-1],yy[2:length(yy)-1],xx[2:length(xx)],yy[2:length(yy)],
             col=col,lwd=lwd,lend=lend)
  } else jj<-NULL
  for(i in 1:nrow(cw$edge)){
    if(i%in%jj==FALSE){
      maps<-cumsum(cw$maps[[i]])/sum(cw$maps[[i]])
      xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
      yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
      for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colors[names(maps)[i]],
                                       lwd=lwd,lend=lend)
    }
  }
  # plot circular lines
  for(i in 1:m+n){
    r<-R[match(i,cw$edge)]
    a1<-min(Y[which(cw$edge==i)])
    a2<-max(Y[which(cw$edge==i)])
    draw.arc(0,0,r,a1,a2,lwd=lwd,col=colors[names(cw$maps[[match(i,cw$edge[,1])]])[1]])
  }
  # plot labels
  for(i in 1:n){
    ii<-which(cw$edge[,2]==i)
    aa<-Y[ii,2]/(2*pi)*360
    adj<-if(aa>90&&aa<270) c(1,0.25) else c(0,0.25)
    tt<-if(aa>90&&aa<270) paste(cw$tip.label[i],paste(rep(" ",offset),
                                                      collapse=""),sep="") else paste(paste(rep(" ",offset),collapse=""),
                                                                                      cw$tip.label[i],sep="")
    aa<-if(aa>90&&aa<270) 180+aa else aa
    if(ftype) text(x[ii,2],y[ii,2],tt,srt=aa,adj=adj,cex=fsize,font=ftype)
  }
  if(setEnv){
    PP<-list(type="fan",use.edge.length=TRUE,node.pos=1,
             show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
             font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
             x.lim=xlim,y.lim=ylim,direction="rightwards",tip.color="black",
             Ntip=Ntip(cw),Nnode=cw$Nnode,edge=cw$edge,
             xx=c(x[sapply(1:n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2],x[1,1],
                  if(m>1) x[sapply(2:m+n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2] else c()),
             yy=c(y[sapply(1:n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2],y[1,1],
                  if(m>1) y[sapply(2:m+n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2] else c()))
    assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
  }
}

## internal function for slanted cladogram
## written by Liam J. Revell 2017-2023
plotCladogram<-function(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,mar=NULL,
                        add=FALSE,offset=NULL,direction="rightwards",xlim=NULL,ylim=NULL,
                        nodes="intermediate",tips=NULL,lend=2,asp=NA,plot=TRUE,underscore=FALSE){
  placement<-nodes
  # set offset fudge (empirically determined)
  offsetFudge<-1.37
  # reorder
  cw<-reorderSimmap(tree)
  pw<-reorderSimmap(tree,"postorder")
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
    for(i in 1:nrow(cw$edge)){
      x<-H[i,1]
      y<-Y[cw$edge[i,1]]
      m<-(Y[cw$edge[i,2]]-Y[cw$edge[i,1]])/(H[i,2]-H[i,1])
      if(is.finite(m)){
        for(j in 1:length(cw$maps[[i]])){
          if(direction=="leftwards")
            lines(c(x,x-cw$maps[[i]][j]),c(y,y-cw$maps[[i]][j]*m),
                  col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
          else lines(c(x,x+cw$maps[[i]][j]),c(y,y+cw$maps[[i]][j]*m),
                     col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
          x<-x+if(direction=="leftwards") -cw$maps[[i]][j] else cw$maps[[i]][j]
          y<-y+if(direction=="leftwards") -m*cw$maps[[i]][j] else m*cw$maps[[i]][j]
          j<-j+1
        }
      } else {
        lines(rep(x,2),Y[cw$edge[i,]],col=colors[names(cw$maps[[i]])[1]],lwd=lwd,
              lend=lend)
      }
    }
    if(direction=="leftwards") pos<-if(par()$usr[1]>par()$usr[2]) 4 else 2
    if(direction=="rightwards") pos<-if(par()$usr[1]>par()$usr[2]) 2 else 4
    for(i in 1:n) if(ftype) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,
                                 offset=offset,cex=fsize,font=ftype)
  }
  PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
           show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
           font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
           x.lim=xlim,y.lim=ylim,
           direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
           edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
                                  function(x,y,z) y[match(x,z)],y=H,z=cw$edge),yy=Y[,1])
  assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
}


## adds legend to an open stochastic map style plot
## written by Liam J. Revell 2013, 2016, 2017
add.simmap.legend<-function(leg=NULL,colors,prompt=TRUE,vertical=TRUE,...){
  if(hasArg(shape)) shape<-list(...)$shape
  else shape<-"square"
  if(prompt){
    cat("Click where you want to draw the legend\n")
    x<-unlist(locator(1))
    y<-x[2]
    x<-x[1]
  } else {
    if(hasArg(x)) x<-list(...)$x
    else x<-0
    if(hasArg(y)) y<-list(...)$y
    else y<-0
  }
  if(hasArg(fsize)) fsize<-list(...)$fsize
  else fsize<-1.0
  if(is.null(leg)) leg<-names(colors)
  h<-fsize*strheight(LETTERS[1])
  w<-par()$mfcol[2]*h*abs(diff(par()$usr[1:2])/diff(par()$usr[3:4]))
  flipped<-par()$usr[1]>par()$usr[2]
  if(vertical){
    y<-y-0:(length(leg)-1)*1.5*h
    x<-rep(x+w/2,length(y))
    text(x + if(flipped) -w else w,y,leg,pos=4,cex=fsize/par()$cex)
  } else {
    sp<-abs(fsize*max(strwidth(leg)))
    x<-x + if(flipped) w/2-0:(length(leg)-1)*1.5*(sp+w) else -w/2+0:(length(leg)-1)*1.5*(sp+w)
    y<-rep(y+w/2,length(x))
    text(x,y,leg,pos=4,cex=fsize/par()$cex)
  }
  if(shape=="square") symbols(x,y,squares=rep(w,length(x)),bg=colors,add=TRUE,inches=FALSE)
  else if(shape=="circle") nulo<-mapply(draw.circle,x=x,y=y,col=colors,
                                        MoreArgs=list(nv=200,radius=w/2))
  else stop(paste("shape=\"",shape,"\" is not a recognized option.",sep=""))
}

# function plots a tree; in the new version this is just a wrapper for plotSimmap
# written by Liam Revell 2012-2017
plotTree<-function(tree,...){
  if(hasArg(color)) color<-list(...)$color
  else color<-NULL
  if(hasArg(fsize)) fsize<-list(...)$fsize
  else fsize<-1.0
  if(hasArg(ftype)) ftype<-list(...)$ftype
  else ftype<-"reg"
  if(hasArg(lwd)) lwd<-list(...)$lwd
  else lwd<-2
  if(hasArg(pts)) pts<-list(...)$pts
  else pts<-FALSE
  if(hasArg(node.numbers)) node.numbers<-list(...)$node.numbers
  else node.numbers<-FALSE
  if(hasArg(mar)) mar<-list(...)$mar
  else mar<-NULL
  if(hasArg(add)) add<-list(...)$add
  else add<-FALSE
  if(hasArg(offset)) offset<-list(...)$offset
  else offset<-NULL
  if(hasArg(type)) type<-list(...)$type
  else type<-"phylogram"
  if(hasArg(direction)) direction<-list(...)$direction
  else direction<-"rightwards"
  if(hasArg(setEnv)) setEnv<-list(...)$setEnv
  else setEnv<-TRUE
  if(hasArg(part)) part<-list(...)$part
  else part<-1.0
  if(hasArg(xlim)) xlim<-list(...)$xlim
  else xlim<-NULL
  if(hasArg(ylim)) ylim<-list(...)$ylim
  else ylim<-NULL
  if(hasArg(nodes)) nodes<-list(...)$nodes
  else nodes<-"intermediate"
  if(hasArg(tips)) tips<-list(...)$tips
  else tips<-NULL
  if(hasArg(maxY)) maxY<-list(...)$maxY
  else maxY<-NULL
  if(hasArg(hold)) hold<-list(...)$hold
  else hold<-TRUE
  if(hasArg(lend)) lend<-list(...)$lend
  else lend<-2
  if(hasArg(asp)) asp<-list(...)$asp
  else asp<-NA
  if(hasArg(plot)) plot<-list(...)$plot
  else plot<-TRUE
  if(hasArg(underscore)) underscore<-list(...)$underscore
  else underscore=FALSE
  if(inherits(tree,"multiPhylo")){
    par(ask=TRUE)
    if(!is.null(color)) names(color)<-"1"
    for(i in 1:length(tree)) plotTree(tree[[i]],color=color,fsize=fsize,ftype=ftype,
                                      lwd=lwd,pts=pts,node.numbers=node.numbers,mar=mar,add=add,offset=offset,
                                      direction=direction,type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,
                                      nodes=nodes,tips=tips,maxY=maxY,hold=hold,lend=lend,asp=asp,plot=plot,
                                      underscore=underscore)
  } else {
    if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
    tree$maps<-as.list(tree$edge.length)
    for(i in 1:length(tree$maps)) names(tree$maps[[i]])<-c("1")
    if(!is.null(color)) names(color)<-"1"
    plotSimmap(tree,colors=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,
               node.numbers=node.numbers,mar=mar,add=add,offset=offset,direction=direction,
               type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,nodes=nodes,tips=tips,maxY=maxY,
               hold=hold,lend=lend,asp=asp,plot=plot,underscore=underscore)
  }
}

## S3 method for objects of class "simmap" & "multiSimmap"
## added by Liam J. Revell 2015
plot.simmap<-function(x,...) plotSimmap(x,...)
plot.multiSimmap<-function(x,...) plotSimmap(x,...)

## function to split vertical plotted lines by the states of daughter edges
## written by Liam J. Revell 2015
splitEdgeColor<-function(tree,colors,lwd=2){
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  for(i in 1:tree$Nnode+Ntip(tree)){
    daughters<-tree$edge[which(tree$edge[,1]==i),2]
    cols<-vector()
    for(j in 1:length(daughters)){
      jj<-which(tree$edge[,2]==daughters[j])
      cols[j]<-if(tree$maps[[jj]][1]==0&&length(tree$maps[[jj]])>1) colors[names(tree$maps[[jj]])[2]] 
      else colors[names(tree$maps[[jj]])[1]]
    }
    ii<-order(obj$yy[c(i,daughters)])
    jj<-order(obj$yy[daughters])
    x0<-x1<-rep(obj$xx[i],length(daughters))
    y0<-obj$yy[c(i,daughters)][ii][1:length(daughters)]
    y1<-obj$yy[c(i,daughters)][ii][2:(length(daughters)+1)]
    cols<-cols[jj]
    for(j in 1:length(x0)) segments(x0[j],y0[j],x1[j],y1[j],col=cols[j],lwd=lwd,lend=2)
  }
}

ace2 <- function (x, phy, type = "continuous", method = if (type == 
                                                            "continuous") "REML" else "ML", CI = TRUE, model = if (type == 
                                                                                                                   "continuous") "BM" else "ER", scaled = TRUE, kappa = 1, 
                  corStruct = NULL, ip = 0.1, use.expm = FALSE, use.eigen = TRUE, 
                  marginal = FALSE,mp=1e+50,print.rates=F,...) 
{
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  if (is.null(phy$edge.length)) 
    stop("tree has no branch lengths")
  type <- match.arg(type, c("continuous", "discrete"))
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (nb.node != nb.tip - 1) 
    stop("\"phy\" is not rooted AND fully dichotomous.")
  if (length(x) != nb.tip) 
    stop("length of phenotypic and of phylogenetic data do not match.")
  if (!is.null(names(x))) {
    if (all(names(x) %in% phy$tip.label)) 
      x <- x[phy$tip.label]
    else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
  }
  obj <- list()
  if (kappa != 1) 
    phy$edge.length <- phy$edge.length^kappa
  if (type == "continuous") {
    switch(method, REML = {
      minusLogLik <- function(sig2) {
        if (sig2 < 0) return(1e+100)
        V <- sig2 * vcv(phy)
        distval <- mahalanobis(x, center = mu, cov = V)
        logdet <- sum(log(eigen(V, symmetric = TRUE, 
                                only.values = TRUE)$values))
        (nb.tip * log(2 * pi) + logdet + distval)/2
      }
      mu <- rep(ace(x, phy, method = "pic")$ace[1], nb.tip)
      out <- nlm(minusLogLik, 1, hessian = TRUE)
      sigma2 <- out$estimate
      se_sgi2 <- sqrt(1/out$hessian)
      tip <- phy$edge[, 2] <= nb.tip
      minus.REML.BM <- function(p) {
        x1 <- p[phy$edge[, 1] - nb.tip]
        x2 <- numeric(length(x1))
        x2[tip] <- x[phy$edge[tip, 2]]
        x2[!tip] <- p[phy$edge[!tip, 2] - nb.tip]
        -(-sum((x1 - x2)^2/phy$edge.length)/(2 * sigma2) - 
            nb.node * log(sigma2))
      }
      out <- nlm(function(p) minus.REML.BM(p), p = rep(mu[1], 
                                                       nb.node), hessian = TRUE)
      obj$resloglik <- -out$minimum
      obj$ace <- out$estimate
      names(obj$ace) <- nb.tip + 1:nb.node
      obj$sigma2 <- c(sigma2, se_sgi2)
      if (CI) {
        se <- .getSEs(out)
        tmp <- se * qt(0.025, nb.node)
        obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
      }
    }, pic = {
      if (model != "BM") stop("the \"pic\" method can be used only with model = \"BM\".")
      phy <- reorder(phy, "postorder")
      phenotype <- numeric(nb.tip + nb.node)
      phenotype[1:nb.tip] <- if (is.null(names(x))) x else x[phy$tip.label]
      contr <- var.con <- numeric(nb.node)
      ans <- .C(C_pic, as.integer(nb.tip), as.integer(phy$edge[, 
                                                               1]), as.integer(phy$edge[, 2]), as.double(phy$edge.length), 
                as.double(phenotype), as.double(contr), as.double(var.con), 
                as.integer(CI), as.integer(scaled))
      obj$ace <- ans[[5]][-(1:nb.tip)]
      names(obj$ace) <- nb.tip + 1:nb.node
      if (CI) {
        se <- sqrt(ans[[7]])
        tmp <- se * qnorm(0.025)
        obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
      }
    }, ML = {
      if (model == "BM") {
        tip <- phy$edge[, 2] <= nb.tip
        dev.BM <- function(p) {
          if (p[1] < 0) return(1e+100)
          x1 <- p[-1][phy$edge[, 1] - nb.tip]
          x2 <- numeric(length(x1))
          x2[tip] <- x[phy$edge[tip, 2]]
          x2[!tip] <- p[-1][phy$edge[!tip, 2] - nb.tip]
          -2 * (-sum((x1 - x2)^2/phy$edge.length)/(2 * 
                                                     p[1]) - nb.node * log(p[1]))
        }
        out <- nlm(function(p) dev.BM(p), p = c(1, rep(mean(x), 
                                                       nb.node)), hessian = TRUE)
        obj$loglik <- -out$minimum/2
        obj$ace <- out$estimate[-1]
        names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
        se <- .getSEs(out)
        obj$sigma2 <- c(out$estimate[1], se[1])
        if (CI) {
          tmp <- se[-1] * qt(0.025, nb.node)
          obj$CI95 <- cbind(obj$ace + tmp, obj$ace - 
                              tmp)
        }
      }
    }, GLS = {
      if (is.null(corStruct)) stop("you must give a correlation structure if method = \"GLS\".")
      if (class(corStruct)[1] == "corMartins") M <- corStruct[1] * 
          dist.nodes(phy)
      if (class(corStruct)[1] == "corGrafen") phy <- compute.brlen(attr(corStruct, 
                                                                        "tree"), method = "Grafen", power = exp(corStruct[1]))
      if (class(corStruct)[1] %in% c("corBrownian", "corGrafen")) {
        dis <- dist.nodes(attr(corStruct, "tree"))
        MRCA <- mrca(attr(corStruct, "tree"), full = TRUE)
        M <- dis[as.character(nb.tip + 1), MRCA]
        dim(M) <- rep(sqrt(length(M)), 2)
      }
      one2n <- 1:nb.tip
      varAY <- M[-one2n, one2n]
      varA <- M[-one2n, -one2n]
      DF <- data.frame(x)
      V <- corMatrix(Initialize(corStruct, DF), corr = FALSE)
      invV <- solve(V)
      o <- gls(x ~ 1, DF, correlation = corStruct)
      GM <- o$coefficients
      obj$ace <- drop(varAY %*% invV %*% (x - GM) + GM)
      names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
      if (CI) {
        se <- sqrt((varA - varAY %*% invV %*% t(varAY))[cbind(1:nb.node, 
                                                              1:nb.node)])
        tmp <- se * qnorm(0.025)
        obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
      }
    })
  }
  else {
    if (method != "ML") 
      stop("only ML estimation is possible for discrete characters.")
    if (any(phy$edge.length < 0)) 
      stop("some branches have negative length")
    if (!is.factor(x)) 
      x <- factor(x)
    if(is.character(model)){
      nl <- nlevels(x)
      lvls <- levels(x)
      x <- as.integer(x)
    } else {
      nl <- ncol(model)
      lvls <- 1:ncol(model) 
      x <- as.integer(x)
    }
    if (is.character(model)) {
      rate <- matrix(NA, nl, nl)
      switch(model, ER = np <- rate[] <- 1, ARD = {
        np <- nl * (nl - 1)
        rate[col(rate) != row(rate)] <- 1:np
      }, SYM = {
        np <- nl * (nl - 1)/2
        sel <- col(rate) < row(rate)
        rate[sel] <- 1:np
        rate <- t(rate)
        rate[sel] <- 1:np
      })
    }else {
      if (ncol(model) != nrow(model)) 
        stop("the matrix given as 'model' is not square")
      if (ncol(model) != nl) 
        stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
      rate <- model
      np <- max(rate)
    }
    index.matrix <- rate
    tmp <- cbind(1:nl, 1:nl)
    index.matrix[tmp] <- NA
    rate[tmp] <- 0
    rate[rate == 0] <- np + 1
    liks <- matrix(0, nb.tip + nb.node, nl)
    TIPS <- 1:nb.tip
    liks[cbind(TIPS, x)] <- 1
    if (anyNA(x)) 
      liks[which(is.na(x)), ] <- 1
    phy <- reorder(phy, "postorder")
    Q <- matrix(0, nl, nl)
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    EL <- phy$edge.length
    if (use.eigen) {
      dev <- function(p, output.liks = FALSE) {
        if (any(is.nan(p)) || any(is.infinite(p))) 
          return(1e+50)
        comp <- numeric(nb.tip + nb.node)
        Q[] <- c(p, 0)[rate]
        diag(Q) <- -rowSums(Q)
        decompo <- eigen(Q)
        lambda <- decompo$values
        GAMMA <- decompo$vectors
        invGAMMA <- solve(GAMMA)
        for (i in seq(from = 1, by = 2, length.out = nb.node)) {
          j <- i + 1L
          anc <- e1[i]
          des1 <- e2[i]
          des2 <- e2[j]
          v.l <- GAMMA %*% diag(exp(lambda * EL[i])) %*% 
            invGAMMA %*% liks[des1, ]
          v.r <- GAMMA %*% diag(exp(lambda * EL[j])) %*% 
            invGAMMA %*% liks[des2, ]
          v <- v.l * v.r
          comp[anc] <- sum(v)
          liks[anc, ] <- v/comp[anc]
        }
        if (output.liks) 
          return(liks[-TIPS, , drop = FALSE])
        dev <- -2 * sum(log(comp[-TIPS]))
        if (is.na(dev)) 
          Inf
        else dev
      }
    }
    else {
      if (!requireNamespace("expm", quietly = TRUE) && 
          use.expm) {
        warning("package 'expm' not available; using function 'matexpo' from 'ape'")
        use.expm <- FALSE
      }
      E <- if (use.expm) 
        expm::expm
      else matexpo
      dev <- function(p, output.liks = FALSE) {
        if (any(is.nan(p)) || any(is.infinite(p))) 
          return(1e+50)
        if(print.rates == T){
          print(p)
        }
        comp <- numeric(nb.tip + nb.node)
        Q[] <- c(p, 0)[rate]
        diag(Q) <- -rowSums(Q)
        for (i in seq(from = 1, by = 2, length.out = nb.node)) {
          j <- i + 1L
          anc <- e1[i]
          des1 <- e2[i]
          des2 <- e2[j]
          v.l <- E(Q * EL[i]) %*% liks[des1, ]
          v.r <- E(Q * EL[j]) %*% liks[des2, ]
          v <- v.l * v.r
          comp[anc] <- sum(v)
          liks[anc, ] <- v/comp[anc]
        }
        if (output.liks) 
          return(liks[-TIPS, , drop = FALSE])
        dev <- -2 * sum(log(comp[-TIPS]))
        if (is.na(dev)) 
          Inf
        else dev
      }
    }
    if(hasArg(rates) == F){
      out <- nlminb(rep(ip, length.out = np), function(p) dev(p), 
                    lower = rep(0, np), upper = rep(mp, np))
      obj$loglik <- -out$objective/2
      obj$rates <- out$par
    } else {
      p <- list(...)$rates
      if(length(p) != np){
        stop("rates supplied do not match number of rates in matrix")
      }
      out <- dev(p)
      obj$loglik <- -out[1]/2
      obj$rates <- p
    }
    oldwarn <- options("warn")
    options(warn = -1)
    out.nlm <- try(nlm(function(p) dev(p), p = obj$rates, 
                       iterlim = 1, stepmax = 0, hessian = TRUE), silent = TRUE)
    options(oldwarn)
    obj$se <- if (inherits(out.nlm, "try-error")) {
      warning("model fit suspicious: gradients apparently non-finite")
      rep(NaN, np)
    }
    else .getSEs(out.nlm)
    obj$index.matrix <- index.matrix
    if (CI) {
      lik.anc <- dev(obj$rates, TRUE)
      if (!marginal) {
        Q[] <- c(obj$rates, 0)[rate]
        diag(Q) <- -rowSums(Q)
        for (i in seq(to = 1, by = -2, length.out = nb.node)) {
          anc <- e1[i] - nb.tip
          des1 <- e2[i] - nb.tip
          if (des1 > 0) {
            P <- matexpo(Q * EL[i])
            tmp <- lik.anc[anc, ]/(lik.anc[des1, ] %*% 
                                     P)
            lik.anc[des1, ] <- (tmp %*% P) * lik.anc[des1, 
            ]
          }
          j <- i + 1L
          des2 <- e2[j] - nb.tip
          if (des2 > 0) {
            P <- matexpo(Q * EL[j])
            tmp <- lik.anc[anc, ]/(lik.anc[des2, ] %*% 
                                     P)
            lik.anc[des2, ] <- (tmp %*% P) * lik.anc[des2, 
            ]
          }
          lik.anc <- lik.anc/rowSums(lik.anc)
        }
      }
      rownames(lik.anc) <- nb.tip + 1:nb.node
      colnames(lik.anc) <- lvls
      obj$lik.anc <- lik.anc
    }
  }
  obj$call <- match.call()
  class(obj) <- "ace"
  obj
}

.getSEs <- function(out)
{
  h <- out$hessian
  if (any(diag(h) == 0)) {
    warning("The likelihood gradient seems flat in at least one dimension (gradient null):\ncannot compute the standard-errors of the transition rates.\n")
    se <- rep(NaN, nrow(h))
  } else {
    se <- sqrt(diag(solve(h)))
  }
  se
}



