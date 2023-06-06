library(phytools)
install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.0.3.tar.gz")
library(rlang)
install.packages("https://cran.r-project.org/src/contrib/Archive/cli/cli_3.2.0.tar.gz")
library(devtools)
install_github("dpasqualin/sfreemap")
library(sfreemap)
library(parallel)

#### SFREEMAP FUNCTIONS ####
sfreemap <- function(tree, tip_states, Q=NULL, type="standard", model="SYM", method="empirical", ...) {
  
  # Am I running on windows? Windows does not have support for the kind of
  # parallelism we are using
  
  # Should this program run in parallel?
  if (hasArg(parallel)) {
    parallel <- list(...)$parallel
    if (isTRUE(parallel) && !support_parallel_mode()) {
      ## parallel requested but parallel mode not supported
      warning('parallel mode is not available on this machine.', call. = FALSE)
      parallel <- FALSE
    }
  } else {
    parallel <- support_parallel_mode()
  }
  
  # When running in parallel, choose how many cores do use.
  # default to all cores available on the machine
  mc.cores <- detectCores()
  if (hasArg(mc.cores)) {
    tmp <- list(...)$mc.cores
    # ignore things like NULL
    if (is.numeric(tmp)) {
      mc.cores <- tmp
    }
  }
  if (mc.cores == 1) {
    parallel <- FALSE
  }
  
  # how many omp threads should be created?
  omp <- 1
  if (hasArg(omp)) {
    omp <- list(...)$omp
  }
  
  # Defining the prior distribution for the root node of the tree,
  # also known as "pi"
  prior <- "equal"
  if (hasArg(prior)) {
    prior <- list(...)$prior
  }
  
  # available types and models
  dna_models <- names(getModels())
  standard_models <- c('SYM', 'ER', 'ARD')
  valid_models <- list(
    "standard" = standard_models
    , "dna" = dna_models
  )
  
  # check for tree object class
  if (all(!inherits(tree, "phylo"), !inherits(tree, "multiPhylo"))) {
    stop("'tree' should be an object of class 'phylo' or 'multiPhylo'")
  }
  
  # check for type
  if (! type %in% names(valid_models)) {
    stop('Unknown type', type)
  }
  
  # check for model
  if (! model %in% valid_models[[type]]) {
    stop('Unknown model ', model)
  }
  
  if (all(inherits(prior, "list"), !inherits(Q, "list"))) {
    stop("if 'prior' is a list 'Q' should be a list with of same size.")
  }
  
  if (all(inherits(prior, "list"), inherits(Q, "list"), length(prior) != length(Q))) {
    stop("if 'prior' and 'Q' are lists, their number of elements must match.")
  }
  
  if (all(inherits(Q, "list"), inherits(tree, "multiPhylo"), length(Q) != length(tree))) {
    stop("if 'Q' is a list and 'tree' is a 'multiPhylo' object, their number of elements must match.")
  }
  
  # a helper function to call sfreemap multiple times, combining
  # trees with Q rate matrices and priors
  call_multiple <- function(idx, tree, tip_states, Q, prior) {
    if (inherits(tree, "multiPhylo")) {
      tree <- tree[[idx]]
    }
    
    if (all(inherits(tip_states, "matrix"), ncol(tip_states) > 1)) {
      tip_states <- tip_states[,idx]
    }
    
    if (all(!is.null(Q), inherits(Q, "list"))) {
      Q <- Q[[idx]]
    }
    
    if (inherits(prior, "list")) {
      prior <- prior[[idx]]
    }
    
    params <- list(
      "tree" <- tree
      , "tip_states" <- tip_states
      , "Q" = Q
      , "prior" = prior
      , "model" = model
      , "type" = type
      , "method" = method
      , "..." = ...
    )
    
    return (do.call(sfreemap, params))
  }
  
  # helper to decide whether to call 'call_multiple' in serial or parallel
  serial_or_parallel <- function(times, tree, tip_states, Q, prior) {
    if (parallel) {
      mtrees <- mclapply(1:times, call_multiple, tree, tip_states, Q
                         , prior, mc.cores=mc.cores)
    } else {
      mtrees <- lapply(1:times, call_multiple, tree, tip_states, Q, prior)
    }
    return (fix_return(mtrees))
  }
  
  # with some combination of parameters we might have a list of multiPhylo
  # objects (a list of a list), so we need to convert it to a single
  # multiPhylo object
  fix_return <- function(mtrees) {
    tmp <- mtrees[[1]]
    if (inherits(tmp, "multiPhylo") || inherits(tmp, "list")) {
      mtrees <- c(mapply(c, mtrees))
    }
    
    if (length(mtrees) > 1) {
      class(mtrees) <- c("sfreemap", "multiPhylo")
    } else {
      mtrees <- mtrees[[1]]
    }
    
    return (mtrees)
  }
  
  # Everything below these tests assume the program is running on with a
  # single tree, single rate matrix and single tip label. So here we check
  # parameters and call sfreemap multiple times if needed.
  if (inherits(tree, "multiPhylo")) {
    # if 'multiPhylo', call sfreemap for each tree
    return(serial_or_parallel(length(tree), tree, tip_states, Q, prior))
  } else if (inherits(Q, "list")) {
    # if multiple rate matrix, call sfreemap for each one.
    # serial_or_parallell will handle the case when we have an equal number
    # of rate matrices and trees, where sfreemap should match tree 1
    # with rate matrix 1, 2 with 2, and so on..
    return(serial_or_parallel(length(Q), tree, tip_states, Q, prior))
  } else if (inherits(prior, "list")) {
    # we can run multiple priors on trees and Q rate matrices
    return(serial_or_parallel(length(prior), tree, tip_states, Q, prior))
  } else if (all(inherits(tip_states, "matrix"), ncol(tip_states) > 1)) {
    # if single tree but multiple tip_states (dna type), call sfreemap
    # for each set of tip label
    return(serial_or_parallel(ncol(tip_states), tree, tip_states, Q, prior))
  } else if (!is.rooted(tree)) {
    # all trees must be rooted
    stop("'tree' must be rooted")
  }
  
  # tol gives the tolerance for zero elements in Q.
  # (Elements less then tol will be reset to tol)
  tol <- 1e-8
  if (hasArg(tol)) {
    tol <- list(...)$tol
  }
  
  # FIXME: this was not tested yet, not sure it it works
  # prior a list containing alpha and beta parameters for the gamma
  # prior distribution on the transition rates in Q. Note that alpha
  # and beta can be single values or vectors, if different prior
  # are desired for each value in Q
  gamma_prior <- list(alpha=1, beta=1, use.empirical=FALSE)
  if (hasArg(gamma_prior)) {
    pr <- list(...)$gamma_prior
    gamma_prior[names(pr)] <- pr
  }
  
  # burn_in for the MCMC
  burn_in <- 1000
  if (hasArg(burn_in)) {
    burn_in <- list(...)$burn_in
  }
  
  # sample_freq for the MCMC
  sample_freq <- 100
  if (hasArg(sample_freq)) {
    sample_freq <- list(...)$sample_freq
  }
  
  # number os simulations for the MCMC
  n_simulation <- 10
  if (hasArg(n_simulation)) {
    n_simulation <- list(...)$n_simulation
  }
  
  # FIXME: this was not tested yet, not sure it it works
  # A single numeric value or a vector containing the (normal)
  # sampling variances for the MCMC
  vQ <- 0.1
  if (hasArg(vQ)) {
    vQ <- list(...)$vQ
  }
  
  # We set the class here so we can use functions like reorder.sfreemap
  class(tree) <- c("sfreemap", "phylo")
  
  if (all(!is.null(Q), is.matrix(Q))) {
    if(hasArg(poss.states)){
      QP <- Q_matrix(tree, tip_states, Q, model, prior, tol, type, poss.states=list(...)$poss.states)
    } else {
      QP <- Q_matrix(tree, tip_states, Q, model, prior, tol, type)
    }
  } else if (type == 'standard') {
    # standard data type has currently two ways of estimating the rate
    # matrix
    if (method == "empirical") {
      QP <- Q_empirical(trees[[1]], tip_states=tips, prior, model, tol, omp)
    } else if (method == "mcmc") {
      QP <- Q_mcmc(tree, tip_states, prior, model, gamma_prior, tol, burn_in
                   , sample_freq, vQ, n_simulation, omp)
      Q <- lapply(QP, function(x) x$Q)
      prior <- lapply(QP, function(x) x$prior)
      # TODO: how to pass on logL?
      return(serial_or_parallel(length(QP), tree, tip_states, Q, prior))
    }
    # Estimating Q when using nucleotide data
  } else if (all(type == "dna", is.null(Q))) {
    QP <- Q_dna(tip_states, tree, model, tol)
  }
  
  # NOTE: It's important to notice that below this point it is garantee that
  # we are dealing with a single tree (a "phylo" object), a single Q matrix,
  # a single prior and a single character.
  
  if (type == "dna") {
    states <- c("a", "c", "t", "g", "-")
  } else {
    if(hasArg(poss.states)){
      states <- list(...)$poss.states
    } else {
      states <- NULL
    }
  }
  
  
  tip_states <- build_states_matrix(tree$tip.label, tip_states, states)
  
  # Set the final value
  Q <- QP$Q
  prior <- QP$prior
  logL <- QP$logL
  
  #Set rownames and colnames
  if(hasArg(poss.states)){
    rownames(Q) <- colnames(Q)<- states
  }
  
  # Vector with rewards
  rewards <- rep(1,nrow(Q))
  if (hasArg(rewards)) {
    rewards <- list(...)$rewards
    if (length(rewards) != nrow(Q)) {
      stop("The rewards vector should represent the states")
    }
  }
  names(rewards) <- colnames(Q)
  
  # Defining the transitions of interest. By default, all transitions.
  QL <- matrix(1, nrow=nrow(Q), ncol=ncol(Q))
  diag(QL) <- 0
  if (hasArg(QL)) {
    QL <- list(...)$QL
    if (!all(dim(Q) == dim(QL))) {
      stop("QL must have same dimensions as Q")
    }
  }
  rownames(QL) <- colnames(QL) <- rownames(Q) 
  
  # Acquire more info about the tree.
  tree_extra <- list(
    states = tip_states
    , n_states = nrow(Q)
    , n_edges = length(tree$edge.length)
    , n_tips = nrow(tip_states)
    , n_nodes = nrow(tip_states) + tree$Nnode
  )
  
  # Reorder the tree so the root is the first row of the matrix.
  # We save the original order to make sure we have the result
  # in same order of the tree;
  tree <- reorder(tree, order='pruningwise')
  # Let's set the elements back to the original tree
  tree[['Q']] <- Q
  tree[['prior']] <- prior
  tree[['logL']] <- logL
  
  # Step 1
  # Compute Eigen values and Eigen vectors for the transition rate matrix
  # Q, assuming Q is symmetrical
  Q_eigen <- eigen(Q, TRUE, only.values = FALSE, symmetric = TRUE)
  # The inverse of Q_eigen vectors
  Q_eigen[['vectors_inv']] <- solve(Q_eigen$vectors)
  
  MAP <- list()
  # Step 2
  # Compute P(tp), the transistion probability, for each edge length t of T
  MAP[['tp']] <- transition_probabilities(Q_eigen, tree$edge.length, omp)
  
  # Step 4 and 5
  MAP[['fl']] <- fractional_likelihoods(tree, tree_extra, Q, Q_eigen
                                        , prior, MAP$tp, tol)
  
  # FIXME: for some unknown reason if I remove this line I get a 'out of
  # bounds' error in posterior_restricted_moment(). This line doesn't need to
  # exist
  MAP[['h']] <- list()
  
  # Build dwelling times
  tree[['mapped.edge']] <- tree[['mapped.edge.lmt']] <- tname <- NULL
  states <- rownames(QL)
  n_states <- tree_extra$n_states
  multiplier <- matrix(0, nrow=n_states, ncol=n_states)
  
  for (i in 1:n_states) {
    for (j in 1:n_states) {
      
      if (i == j) {
        value <- rewards[i] # dwelling times
      } else {
        value <- QL[i,j] * Q[i,j] # number of transitions
      }
      
      if (value == 0) {
        next
      }
      multiplier[i,j] <- value
      
      # Step 3
      h <- func_H(multiplier, Q_eigen, tree, tree_extra, omp)
      prm <- posterior_restricted_moment(h, tree, tree_extra, MAP, omp)
      ev <- expected_value(tree, Q, MAP, prm)
      
      if (i == j) {
        # dwelling times
        tree[['mapped.edge']] <- cbind(tree[['mapped.edge']], ev)
      } else {
        # number of transitions
        state_from <- states[i]
        state_to <- states[j]
        tname <- c(tname, paste(state_from, state_to, sep=','))
        tree[['mapped.edge.lmt']] <- cbind(tree[['mapped.edge.lmt']], ev)
      }
      
      multiplier[i,j] <- 0
    }
  }
  bname <- paste(tree$edge[,1], ",", tree$edge[,2], sep="")
  colnames(tree[['mapped.edge']]) <- names(rewards)
  colnames(tree[['mapped.edge.lmt']]) <- tname
  rownames(tree[['mapped.edge']]) <- rownames(tree[['mapped.edge.lmt']]) <- bname
  
  # Return the tree in the original order
  return (reorder(tree, order='cladewise'))
}

expected_value <- function(tree, Q, map, prm) {
  
  likelihood <- map[['fl']][['L']]
  
  ret <- prm / likelihood
  
  # the rownames of the mapped objects
  names(ret) <- paste(tree$edge[,1], ",", tree$edge[,2], sep="")
  
  return(ret)
}

build_states_matrix <- function(tip_labels, tip_states, possible_states=NULL) {
  if (is.null(possible_states)) {
    possible_states <- tip_states
  }
  res <- to.matrix(tip_states, sort(unique(possible_states)))
  return (res[tip_labels,])
}

getModels <- function() {
  models <- list(
    JC = list(optQ=FALSE, optBf=FALSE,   subs=c(0, 0, 0, 0, 0, 0)),
    F81 = list(optQ=FALSE, optBf=TRUE,   subs=c(0, 0, 0, 0, 0, 0)),
    K80 = list(optQ=TRUE, optBf=FALSE,   subs=c(0, 1, 0, 0, 1, 0)),
    HKY = list(optQ=TRUE, optBf=TRUE,    subs=c(0, 1, 0, 0, 1, 0)),
    TrNe = list(optQ=TRUE, optBf=FALSE,  subs=c(0, 1, 0, 0, 2, 0)),
    TrN = list(optQ=TRUE, optBf=TRUE,    subs=c(0, 1, 0, 0, 2, 0)),
    TPM1 = list(optQ=TRUE, optBf=FALSE,  subs=c(0, 1, 2, 2, 1, 0)),
    K81 = list(optQ=TRUE, optBf=FALSE,   subs=c(0, 1, 2, 2, 1, 0)),
    TPM1u = list(optQ=TRUE, optBf=TRUE,  subs=c(0, 1, 2, 2, 1, 0)),
    TPM2 = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 1, 0, 2, 0)),
    TPM2u = list(optQ=TRUE, optBf=TRUE,  subs=c(1, 2, 1, 0, 2, 0)),
    TPM3 = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 0, 1, 2, 0)),
    TPM3u = list(optQ=TRUE, optBf=TRUE,  subs=c(1, 2, 0, 1, 2, 0)),
    TIM1e = list(optQ=TRUE, optBf=FALSE, subs=c(0, 1, 2, 2, 3, 0)),
    TIM1 = list(optQ=TRUE, optBf=TRUE,   subs=c(0, 1, 2, 2, 3, 0)),
    TIM2e = list(optQ=TRUE, optBf=FALSE, subs=c(1, 2, 1, 0, 3, 0)),
    TIM2 = list(optQ=TRUE, optBf=TRUE,   subs=c(1, 2, 1, 0, 3, 0)),
    TIM3e = list(optQ=TRUE, optBf=FALSE, subs=c(1, 2, 0, 1, 3, 0)),
    TIM3 = list(optQ=TRUE, optBf=TRUE,   subs=c(1, 2, 0, 1, 3, 0)),
    TVMe = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 3, 4, 2, 0)),
    TVM = list(optQ=TRUE, optBf=TRUE,    subs=c(1, 2, 3, 4, 2, 0)),
    SYM = list(optQ=TRUE, optBf=FALSE,   subs=c(1, 2, 3, 4, 5, 0)),
    GTR = list(optQ=TRUE, optBf=TRUE,    subs=c(1, 2, 3, 4, 5, 0))
  )
  return (models)
}

Q_matrix <- function(tree, tip_states, Q, model, prior, tol, type, ...) {
  
  if (type == "dna") {
    states <- c("a", "c", "t", "g", "-")
  } else if (hasArg(poss.states)){
    states <- list(...)$poss.states
  } else {
    states <- NULL
  }
  
  tip_states <- build_states_matrix(tree$tip.label, tip_states, states)
  
  states <- tip_states/rowSums(tip_states)
  
  n_states <- ncol(states)
  
  # Reorder tree and create a copy called bt. Not sure why
  # phytools need this...
  new_tree <- bt <- reorder.phylo(tree, "cladewise")
  
  if(hasArg(poss.states)){
    XX <- getPars(bt, states, model, Q=Q, new_tree, tol, n_states,poss.states=list(...)$poss.states)
  } else {
    XX <- getPars(bt, states, model, Q=Q, new_tree, tol, n_states)
  }
  # NOTE: do we need this?
  L <- XX$L
  logL <- XX$loglik
  
  # Set the priors
  if (prior[1] == "equal") {
    prior <- setNames(rep(1/n_states,n_states), colnames(L)) # set equal
  } else if (prior[1] == "estimated") {
    prior <- statdist(Q) # set from stationary distribution
  } else {
    prior <- prior/sum(prior) # obtain from input
  }
  
  return (list(Q=Q, prior=prior, logL=logL))
}

getPars<-function(bt,xx,model,Q,tree,tol,m,omp=1,liks=TRUE, ...){
  if(hasArg(poss.states)){
    XX<-apeAce(bt,xx,model,omp,fixedQ=Q,output.liks=liks,poss.states=list(...)$poss.states)
  } else {
    XX<-apeAce(bt,xx,model,omp,fixedQ=Q,output.liks=liks)
  }
  N<-length(bt$tip.label)
  II<-XX$index.matrix
  lvls<-XX$states
  if(liks){
    L<-XX$lik.anc
    rownames(L)<-N+1:nrow(L)
    if(!is.binary.tree(tree)){
      ancNames<-matchNodes(tree,bt)
      L<-L[as.character(ancNames[,2]),]
      rownames(L)<-ancNames[,1]
    }
    L<-rbind(xx,L)
    rownames(L)[1:N]<-1:N
  } else L<-NULL
  if(any(XX$rates<tol)){
    #message(paste("\nWarning: some elements of Q not numerically distinct from 0; setting to",tol,"\n"))
    XX$rates[XX$rates<tol]<-tol
  }
  Q<-matrix(XX$rates[II],m,m,dimnames=list(lvls,lvls))
  Q[Q<0.000000001] <- 0
  diag(Q)<--rowSums(Q,na.rm=TRUE)
  return(list(Q=Q,L=L,loglik=XX$loglik))
}

apeAce <- function(tree,x,model,omp,fixedQ=NULL,...){
  if(hasArg(output.liks)) output.liks<-list(...)$output.liks
  else output.liks<-TRUE
  ip<-0.1
  nb.tip<-length(tree$tip.label)
  nb.node<-tree$Nnode
  if(is.matrix(x)){
    x<-x[tree$tip.label,]
    nl<-ncol(x)
    lvls<-colnames(x)
  } else {
    x<-tip.states[tree$tip.label]
    if(hasArgs(poss.states)){
      if(!is.factor(x)) x<-factor(x,levels = list(...)$poss.states) 
    } else {
      if(!is.factor(x)) x<-factor(x)
    }
    nl<-nlevels(x)
    lvls<-levels(x)
    x<-as.integer(x)
  }
  if(is.null(fixedQ)){
    if(is.character(model)){
      rate<-matrix(NA,nl,nl)
      if(model=="ER") np<-rate[]<-1
      if(model=="ARD"){
        np<-nl*(nl-1)
        rate[col(rate)!=row(rate)]<-1:np
      }
      if (model=="SYM") {
        np<-nl*(nl-1)/2
        sel<-col(rate)<row(rate)
        rate[sel]<-1:np
        rate<-t(rate)
        rate[sel]<-1:np
      }
    } else {
      if(ncol(model)!=nrow(model)) stop("the matrix given as 'model' is not square")
      if(ncol(model)!=nl) stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
      rate<-model
      np<-max(rate)
    }
    Q<-matrix(0,nl,nl)
  } else {
    rate<-matrix(NA,nl,nl)
    np<-nl*(nl-1)
    rate[col(rate)!=row(rate)]<-1:np
    Q<-fixedQ
  }
  index.matrix<-rate
  tmp<-cbind(1:nl,1:nl)
  index.matrix[tmp]<-NA
  rate[tmp]<-0
  rate[rate==0]<-np+1
  liks<-matrix(0,nb.tip+nb.node,nl)
  TIPS<-1:nb.tip
  if(is.matrix(x)) liks[TIPS,]<-x
  else liks[cbind(TIPS,x)]<-1
  phy<-reorder.phylo(tree,"pruningwise")
  
  dev<-function(p,output.liks=FALSE,fixedQ=NULL){
    if(any(is.nan(p))||any(is.infinite(p))) return(1e50)
    comp<-numeric(nb.tip+nb.node)
    if(is.null(fixedQ)){
      Q[]<-c(p,0)[rate]
      diag(Q)<--rowSums(Q)
    } else Q<-fixedQ
    
    Q_eigen <- eigen(Q, TRUE, only.values = FALSE)
    Q_eigen[['vectors_inv']] <- solve(Q_eigen$vectors)
    tb <- transition_probabilities(Q_eigen, phy$edge.length, omp)
    
    for(i in seq(from=1,by=2,length.out=nb.node)){
      j<-i+1L
      anc<-phy$edge[i,1]
      des1<-phy$edge[i,2]
      des2<-phy$edge[j,2]
      v.l<-tb[,,i]%*%liks[des1,]
      v.r<-tb[,,j]%*%liks[des2,]
      v<-v.l*v.r
      comp[anc]<-sum(v)
      liks[anc,]<-v/comp[anc]
    }
    if(output.liks) return(liks[-TIPS,])
    dev<--2*sum(log(comp[-TIPS]))
    if(is.na(dev)) Inf else dev
  }
  if(is.null(fixedQ)){
    out<-nlminb(rep(ip,length.out=np),function(p) dev(p),lower=rep(0,np),upper=rep(1e50,np))
    obj<-list()
    obj$loglik<--out$objective/2
    obj$rates<-out$par
    obj$index.matrix<-index.matrix
    if(output.liks){
      obj$lik.anc<-dev(obj$rates,TRUE)
      colnames(obj$lik.anc)<-lvls
    }
    obj$states<-lvls
  } else {
    out<-dev(rep(ip,length.out=np),fixedQ=Q)
    obj<-list()
    obj$loglik<--out/2
    diag(fixedQ) <- NA
    obj$rates<-fixedQ#[sapply(1:nl,function(x,y) which(x==y),index.matrix)]
    obj$rates <- obj$rates[!is.na(obj$rates)]
    obj$index.matrix<-index.matrix
    if(output.liks){
      obj$lik.anc<-dev(unlist(as.vector(obj$rates)),TRUE,fixedQ=Q)
      colnames(obj$lik.anc)<-lvls
    }
    obj$states<-lvls
  }
  return(obj)
}

transition_probabilities <- function(Q_eigen, edges, omp) {
  .Call('sfreemap_transition_probabilities', PACKAGE = 'sfreemap', Q_eigen, edges, omp)
}

func_H <- function(multiplier, Q_eigen, tree, tree_extra, omp) {
  .Call('sfreemap_func_H', PACKAGE = 'sfreemap', multiplier, Q_eigen, tree, tree_extra, omp)
}

posterior_restricted_moment <- function(m, tree, tree_extra, map, omp) {
  .Call('sfreemap_posterior_restricted_moment', PACKAGE = 'sfreemap', m, tree, tree_extra, map, omp)
}

fractional_likelihoods <- function(tree, tree_extra, q, q_eigen, prior, trans_prob, tol) {
  .Call('sfreemap_fractional_likelihoods', PACKAGE = 'sfreemap', tree, tree_extra, q, q_eigen, prior, trans_prob, tol)
}

#### BUILD QMATRIX AND QLMATRIX####

#Load matrix
matrix <- read.csv("../data/transition_matrix/matrix_mammalian_empty.csv",header=F)

#matrix first element = 0
matrix[1,1] <- 0

#Rownames and colnames
rownames(matrix) <- matrix[,1]

colnames(matrix) <- matrix[1,]

#Create QL matrix 
matrix2 <- matrix

#Loop through rows and columns to fill
for(i in 2:nrow(matrix)){
  for(j in 2:ncol(matrix)){
    
    #AA Fisions
    if(matrix[1,j]==matrix[i,1]+1){
      matrix[i,j] <- 57.2580
      matrix2[i,j] <- 1
      
      #AA Fusions  
    } else if(matrix[1,j]+1==matrix[i,1]){
      matrix[i,j] <- 58.4389
      matrix2[i,j] <- 1
      
      #SA Fusions
    } else if(matrix[1,j]==matrix[i,1]+49){
      matrix[i,j] <- 1.018652
      matrix2[i,j] <- 1
      
      #Neo XY -> XY
    } else if(matrix[1,j]==matrix[i,1]-50){
      matrix[i,j] <- 1.778527
      matrix2[i,j] <- 1
      
    } else {
      matrix[i,j] <- 0
      matrix2[i,j] <- 0
    }
  } 
}

#Remove first column and row
matrix <- matrix[-1,-1]
matrix2 <- matrix2[-1,-1]

#Reset colnames/rownames to 1-98
rownames(matrix) <- NULL
colnames(matrix) <- NULL
rownames(matrix2) <- NULL
colnames(matrix2) <- NULL


#Set rowsums equal to 0
diag(matrix) <- -rowSums(matrix)

#Get into right format and names
Q <- as.matrix(matrix)
QL <- as.matrix(matrix2)

#### LOAD DATA ####

data <- read.csv("../data/chromes/mammal_chroms_sim_state.csv",
                 as.is=T)[,c(5,9,11,12,15)]
tree <- read.tree("../data/trees/4705sp_mean.nwk")

#### DATA ORGANIZATION ####

#Remove taxa with uncertain sim states
data <- subset(data,is.na(sim.state) == F)

#Convert species names to right format 
data$tree.name <- tolower(data$tree.name)

#Eliminate spaces used to designate subspecies names and replace with underscores
data$tree.name <- gsub(" ","_",data$tree.name)

#Trim phylogeny and chromosome data
data <- data[data$tree.name %in% tree$tip.label, ]
rownames(data) <- 1:nrow(data)
tree <- keep.tip(tree, data$tree.name)

#Scale tree to unit length
tree$edge.length <- tree$edge.length/max(branching.times(tree))

#Identify and remove multiple sim states for taxa (for now, troubleshoot in future)
data.dup <- data[duplicated(data$tree.name),]

#Remove duplicate elements
data.unique <- data[!duplicated(data$tree.name),]

#Use temp df to reorder tip data to match tree order
data.alt <- data.unique

for(j in 1:length(tree$tip.label)){
  index <- which(data.unique$tree.name == tree$tip.label[j])
  data.alt[j,] <- data.unique[index,]
}

data <- data.unique <- data.alt

rm(data.alt,data.unique)

#Create tip states vector
tip.states <- data$sim.state

names(tip.states) <- data$tree.name



#### RUN SFREEMAP ####
mammals.map <- sfreemap(tree=tree,
                        tip_states=tip.states,
                        Q=Q,
                        poss.states=1:98,
                        QL=QL,
                        parallel=F)

#### SAVE OUTPUTS ####

#Output Q matrix
out.Q <- as.data.frame(mammals.map$Q)

#Output priors
out.pi <- as.data.frame(mammals.map$prior)
out.pi <- as.data.frame(cbind(rownames(out.pi),out.pi[,1]))

colnames(out.pi) <- c("state","priors")
rownames(out.pi) <- NULL

#Output logL
out.logL <- as.data.frame(mammals.map$logL)
colnames(out.logL) <- "logL"

#Output mapped dwelling times
out.times <- as.data.frame(mammals.map$mapped.edge)
out.times <- as.data.frame(cbind(rownames(out.times),out.times))

colnames(out.times)[1] <- "edge"
rownames(out.times) <- NULL

#Output mapped transitions
out.trans <- as.data.frame(mammals.map$mapped.edge.lmt)
out.trans <- as.data.frame(cbind(rownames(out.trans),out.trans))

colnames(out.trans)[1] <- "edge"
rownames(out.trans) <- NULL



#Save all
write.csv(out.Q,"../data/sfreemap_output/Q.csv",row.names = F,quote=F)
write.csv(out.pi,"../data/sfreemap_output/priors.csv",row.names = F,quote=F)
write.csv(out.logL,"../data/sfreemap_output/logL.csv",row.names = F,quote=F)
write.csv(out.times,"../data/sfreemap_output/times.csv",row.names = F,quote=F)
write.csv(out.trans,"../data/sfreemap_output/trans.csv",row.names = F,quote=F)


