make.simmap2 <- function (tree, x, model, nsim, rejmax = NULL, rejint = 1000000, monitor = FALSE, ...){
  if (inherits(tree, "multiPhylo")) {
    ff <- function(yy, x, model, nsim, ...) {
      zz <- make.simmap(yy, x, model, nsim, ...)
      if (nsim > 1) 
        class(zz) <- NULL
      return(zz)
    }
    if (nsim > 1) 
      mtrees <- unlist(sapply(tree, ff, x, model, nsim, 
                              ..., simplify = FALSE), recursive = FALSE)
    else mtrees <- sapply(tree, ff, x, model, nsim, ..., 
                          simplify = FALSE)
    class(mtrees) <- c("multiSimmap", "multiPhylo")
  } else {
    if (hasArg(pi)) 
      pi <- list(...)$pi
    else pi <- "equal"
    if (hasArg(message)) 
      pm <- list(...)$message
    else pm <- TRUE
    if (hasArg(tol)) 
      tol <- list(...)$tol
    else tol <- 0
    if (hasArg(Q)) 
      Q <- list(...)$Q
    else Q <- "empirical"
    if (hasArg(burnin)) 
      burnin <- list(...)$burnin
    else burnin <- 1000
    if (hasArg(samplefreq)) 
      samplefreq <- list(...)$samplefreq
    else samplefreq <- 100
    if (hasArg(vQ)) 
      vQ <- list(...)$vQ
    else vQ <- 0.1
    prior <- list(alpha = 1, beta = 1, use.empirical = FALSE)
    if (hasArg(prior)) {
      pr <- list(...)$prior
      prior[names(pr)] <- pr
    }
    if (!inherits(tree, "phylo")) 
      stop("tree should be object of class \"phylo\".")
    if (!is.matrix(x)) 
      xx <- to.matrix(x, sort(unique(x)))
    else xx <- x
    xx <- xx[tree$tip.label, ]
    xx <- xx/rowSums(xx)
    tree <- bt <- reorder.phylo(tree, "cladewise")
    if (!is.binary(bt)) 
      bt <- multi2di(bt, random = FALSE)
    N <- Ntip(tree)
    m <- ncol(xx)
    root <- N + 1
    if (is.character(Q) && Q == "empirical") {
      XX <- getPars(bt, xx, model, Q = NULL, tree, tol, 
                    m, pi = pi, args = list(...))
      L <- XX$L
      Q <- XX$Q
      logL <- XX$loglik
      pi <- XX$pi
      if (pi[1] == "equal") 
        pi <- setNames(rep(1/m, m), colnames(L))
      else if (pi[1] == "estimated") 
        pi <- statdist(Q)
      else if (pi[1] == "fitzjohn") 
        pi <- "fitzjohn"
      else pi <- pi/sum(pi)
      if (pm) 
        printmessage(Q, pi, method = "empirical")
      mtrees <- replicate(nsim,
                            smap2(tree, x, N, m, root, L, Q, pi, logL, rejmax, rejint, monitor),
                          simplify = FALSE)
    }
    else if (is.character(Q) && Q == "mcmc") {
      if (prior$use.empirical) {
        qq <- fitMk(bt, xx, model)$rates
        prior$alpha <- qq * prior$beta
      }
      get.stationary <- if (pi[1] == "estimated") 
        TRUE
      else FALSE
      if (pi[1] %in% c("equal", "estimated")) 
        pi <- setNames(rep(1/m, m), colnames(xx))
      else if (pi[1] == "fitzjohn") 
        pi <- "fitzjohn"
      else pi <- pi/sum(pi)
      XX <- mcmcQ(bt, xx, model, tree, tol, m, burnin, 
                  samplefreq, nsim, vQ, prior, pi = pi)
      L <- lapply(XX, function(x) x$L)
      Q <- lapply(XX, function(x) x$Q)
      logL <- lapply(XX, function(x) x$loglik)
      pi <- if (get.stationary) 
        lapply(Q, statdist)
      else if (pi[1] == "fitzjohn") 
        lapply(XX, function(x) x$pi)
      else lapply(1:nsim, function(x, y) y, y = pi)
      if (pm) 
        printmessage(Reduce("+", Q)/length(Q), Reduce("+", 
                                                      pi)/length(pi), method = "mcmc")
      mtrees <- if (nsim > 1) 
        mapply(smap2, L = L, Q = Q, pi = pi, logL = logL, 
               MoreArgs = list(tree = tree, x = x, N = N, 
                               m = m, root = root, rejmax = rejmax,
                               rejint = rejint, monitor = monitor), SIMPLIFY = FALSE)
      else list(smap2(tree = tree, x = x, N = N, m = m, 
                     root = root, rejmax = rejmax, rejint = rejint, monitor = monitor,
                     L = L[[1]], Q = Q[[1]], pi = pi[[1]], logL = logL[[1]]))
    }
    else if (is.matrix(Q)) {
      XX <- getPars(bt, xx, model, Q = Q, tree, tol, m, 
                    pi = pi, args = list(...))
      L <- XX$L
      logL <- XX$loglik
      pi <- XX$pi
      if (pi[1] == "equal") 
        pi <- setNames(rep(1/m, m), colnames(L))
      else if (pi[1] == "estimated") 
        pi <- statdist(Q)
      else if (pi[1] == "fitzjohn") 
        pi <- "fitzjohn"
      else pi <- pi/sum(pi)
      if (pm) 
        printmessage(Q, pi, method = "fixed")
      mtrees <- replicate(nsim,
                          c(if(monitor == TRUE){print(paste("simulation", "X", "of", nsim, sep = " "))},
                           smap2(tree, x, N, m, root, L, Q, pi, logL, rejmax, rejint, monitor)),
                          simplify = FALSE)
    }
    if (length(mtrees) == 1) 
      mtrees <- mtrees[[1]]
    else class(mtrees) <- c("multiSimmap", "multiPhylo")
  }
  (if (hasArg(message)) 
    list(...)$message
    else TRUE)
  if ((if (hasArg(message)) 
    list(...)$message
    else TRUE) && inherits(tree, "phylo")) 
    message("Done.")
  return(mtrees)
}

getPars <- function(bt,xx,model,Q,tree,tol,m,liks=TRUE,pi,args=list()){
  if(!is.null(args$pi)) args$pi<-NULL
  args<-c(list(tree=bt,x=xx,model=model,fixedQ=Q,output.liks=liks,pi=pi),args)
  obj<-do.call(fitMk,args)
  N<-length(bt$tip.label)
  pi<-obj$pi
  II<-obj$index.matrix+1
  lvls<-obj$states
  if(liks){
    L<-obj$lik.anc
    rownames(L)<-N+1:nrow(L)
    if(!is.binary(tree)){
      ancNames<-matchNodes(tree,bt)
      L<-L[as.character(ancNames[,2]),]
      rownames(L)<-ancNames[,1]
    }
    L<-rbind(xx,L)
    rownames(L)[1:N]<-1:N
  } else L<-NULL	
  Q<-matrix(c(0,obj$rates)[II],m,m,dimnames=list(lvls,lvls))
  if(any(rowSums(Q,na.rm=TRUE)<tol)){
    message(paste("\nWarning: some rows of Q not numerically distinct from 0; setting to",tol,"\n"))
    ii<-which(rowSums(Q,na.rm=TRUE)<tol)
    for(i in 1:length(ii)) Q[ii[i],setdiff(1:ncol(Q),ii[i])]<-tol/(ncol(Q)-1)
  }
  diag(Q)<--rowSums(Q,na.rm=TRUE)
  return(list(Q=Q,L=L,loglik=logLik(obj),pi=pi))
}

printmessage <- function(Q,pi,method){
  if(method=="empirical"||method=="fixed")
    cat("make.simmap is sampling character histories conditioned on\nthe transition matrix\n\nQ =\n")
  else if(method=="mcmc"){
    cat("make.simmap is simulating with a sample of Q from\nthe posterior distribution\n")
    cat("\nMean Q from the posterior is\nQ =\n")
  }
  print(Q)
  if(method=="empirical") cat("(estimated using likelihood);\n")
  else if(method=="fixed") cat("(specified by the user);\n")
  cat("and (mean) root node prior probabilities\npi =\n")
  if(is.list(pi)) pi<-Reduce("+",pi)/length(pi)
  print(pi)
  flush.console()
}

smap2 <- function(tree,x,N,m,root,L,Q,pi,logL,rejmax,rejint,monitor){
  # create the map tree object
  mtree<-tree
  mtree$maps<-list()
  mtree$mapped.edge<-matrix(0,nrow(tree$edge),m,dimnames=list(paste(tree$edge[,1],",",
                                                                    tree$edge[,2],sep=""),colnames(L)))
  # now we want to simulate the node states & histories by pre-order traversal
  NN<-matrix(NA,nrow(tree$edge),2) # our node values
  NN[which(tree$edge[,1]==root),1]<-rstate(L[as.character(root),]/
                                             sum(L[as.character(root),])) # assign root
  for(j in 1:nrow(tree$edge)){
    # conditioned on the start value, assign end value of node (if internal)
    p<-EXPM(Q*tree$edge.length[j])[NN[j,1],]*L[as.character(tree$edge[j,2]),]
    NN[j,2]<-rstate(p/sum(p))
    NN[which(tree$edge[,1]==tree$edge[j,2]),1]<-NN[j,2]
    # now simulate on the branches
    accept <- FALSE
    counter <- 0
    while(!accept){
      map<-sch(NN[j,1],tree$edge.length[j],Q)
      if (counter == rejmax){
        if (monitor == TRUE){
          print(paste("branch", j, "of", nrow(tree$edge),
                    "has exceeded the rejection limit of", format(rejmax, scientific = FALSE),
                    "and will be skipped", sep = " "))
        }
        map <- NA
        accept = TRUE
      } else
        if(names(map)[length(map)]==NN[j,2]){
        if (monitor == TRUE){
          print(paste("branch", j, "of", nrow(tree$edge),
                    "ACCEPTED after", format(counter, scientific = FALSE),
                    "total rejections", sep = " "))
        }
        accept = TRUE
      } else {
        counter <- counter + 1
        if ((counter/rejint) %% 1 == 0){
          if (monitor == TRUE){  
            print(paste("branch", j, "of", nrow(tree$edge),
                    "rejected with", format(counter, scientific = FALSE),
                    "total rejections", sep = " "))
          }
        }
      }
    }
    mtree$maps[[j]]<-map
    for(k in 1:length(mtree$maps[[j]])){
      mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-
      mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]+mtree$maps[[j]][k]
    }
  }
  for(L in 1:length(mtree$maps)){
    if(is.na(mtree$maps[[L]][1]) == TRUE){
      named.edge <- mtree$edge.length[L]
      names(named.edge) <- "fail"
      mtree$maps[[L]] <- named.edge
    }
  }
  mtree$Q<-Q
  mtree$logL<-logL
  if(!inherits(mtree,"simmap")) class(mtree)<-c("simmap",setdiff(class(mtree),"simmap"))
  attr(mtree,"map.order")<-"right-to-left"
  return(mtree)
}

EXPM <- function(x,...){
  e_x<-if(isSymmetric(x)) matexpo(x) else expm(x,...)
  dimnames(e_x)<-dimnames(x)
  e_x
}

expm <- function (x, method = c("Higham08.b", "Higham08", "AlMohy-Hi09", 
                  "Ward77", "PadeRBS", "Pade", "Taylor", "PadeO", "TaylorO", 
                  "R_Eigen", "R_Pade", "R_Ward77", "hybrid_Eigen_Ward"), order = 8, 
                  trySym = TRUE, tol = .Machine$double.eps, do.sparseMsg = TRUE, 
                  preconditioning = c("2bal", "1bal", "buggy")){
    stopifnot(is.numeric(x) || (isM <- inherits(x, "dMatrix")) || 
                inherits(x, "mpfrMatrix"))
    if (length(d <- dim(x)) != 2) 
      stop("argument is not a matrix")
    if (d[1] != d[2]) 
      stop("matrix not square")
    method <- match.arg(method)
    checkSparse <- !nzchar(Sys.getenv("R_EXPM_NO_DENSE_COERCION"))
    isM <- !is.numeric(x) && isM
    if (isM && checkSparse) {
      if (!(method %in% expm.methSparse) && is(x, "sparseMatrix")) {
        if (do.sparseMsg) 
          message("coercing to dense matrix, as required by method ", 
                  dQuote(method))
        x <- as(x, "denseMatrix")
      }
    }
    switch(method, `AlMohy-Hi09` = expm.AlMoHi09(x, p = order), 
           Higham08.b = expm.Higham08(x, balancing = TRUE),
           Higham08 = expm.Higham08(x, balancing = FALSE), Ward77 = {
            if (!is.numeric(x)) x <- as(x, "matrix")
              switch(match.arg(preconditioning), `2bal` = .Call(do_expm, 
                x, "Ward77"), `1bal` = .Call(do_expm, x, "Ward77_1"), 
                buggy = .Call(do_expm, x, "buggy_Ward77"), stop("invalid 'preconditioning'"))
              }, R_Eigen = {
                  isSym <- if (trySym) isSymmetric.matrix(x) else FALSE
                  z <- eigen(x, symmetric = isSym)
                  V <- z$vectors
                  Vi <- if (isSym) t(V) else solve(V)
                  Re(V %*% (exp(z$values) * Vi))
                  }, hybrid_Eigen_Ward = {
                    if (!is.numeric(x)) x <- as(x, "matrix")
                      .Call(do_expm_eigen, x, tol)
                      }, R_Pade = {
                        stopifnot(order >= 2)
                        expm.s.Pade.s(x, order, n = d[1])
                          }, R_Ward77 = {
                            stopifnot(order >= 2)
                            n <- d[1]
                            trShift <- sum(d.x <- diag(x))
                            if (trShift) {
                              trShift <- trShift/n
                              diag(x) <- d.x - trShift
                              }
                                baP <- balance(x, "P")
                                baS <- balance(baP$z, "S")
                                x <- expm.s.Pade.s(baS$z, order)
                                d <- baS$scale
                                x <- x * (d * rep(1/d, each = n))
                                pp <- as.integer(baP$scale)
                                if (baP$i1 > 1) {
                                 for (i in (baP$i1 - 1):1) {
                                  tt <- x[, i]
                                  x[, i] <- x[, pp[i]]
                                  x[, pp[i]] <- tt
                                  tt <- x[i, ]
                                  x[i, ] <- x[pp[i], ]
                                  x[pp[i], ] <- tt
                                  }
                                }
                                if (baP$i2 < n) {
                                 for (i in (baP$i2 + 1):n) {
                                  tt <- x[, i]
                                  x[, i] <- x[, pp[i]]
                                  x[, pp[i]] <- tt
                                  tt <- x[i, ]
                                  x[i, ] <- x[pp[i], ]
                                  x[pp[i], ] <- tt
                                  }
                                }
                                if (trShift) {
                                  exp(trShift) * x
                                } else x
                              }, PadeRBS = {
                              if (!is.numeric(x)) x <- as(x, "matrix")
                                stopifnot((order <- as.integer(order)) >= 1)
                                if (!is.double(x)) storage.mode(x) <- "double"
                                Fobj <- .Fortran(matexpRBS, order, as.integer(d[1]), 
                                                  T = 1, H = x, iflag = integer(1))[c("H", "iflag")]
                                if (Fobj[["iflag"]] < 0) stop("Unable to determine matrix exponential")
                                    Fobj[["H"]]
                                }, {
                                  if (!is.numeric(x)) x <- as(x, "matrix")
                                  if (!is.double(x)) storage.mode(x) <- "double"
                                    order <- as.integer(order)
                                    ntaylor <- npade <- 0L
                                    if (substr(method, 1, 4) == "Pade") npade <- order else ntaylor <- order
                                      res <- if (identical(grep("O$", method), 1L)) .Fortran(matrexpO, 
                                                  X = x, size = d[1], ntaylor, npade, accuracy = double(1))[c("X", 
                                                  "accuracy")] else .Fortran(matrexp, X = x, size = d[1], 
                                                  ntaylor, npade, accuracy = double(1))[c("X", 
                                                  "accuracy")]
                                                  structure(res$X, accuracy = res$accuracy)
                                    })
  }

sch <- function(start,t,Q){
  tol<-t*1e-12
  dt<-setNames(0,start)
  while(sum(dt)<(t-tol)){
    s<-names(dt)[length(dt)]
    dt[length(dt)]<-if(-Q[s,s]>0) rexp(n=1,rate=-Q[s,s]) else t-sum(dt)
    if(sum(dt)<(t-tol)){
      dt<-c(dt,0)
      if(sum(Q[s,][-match(s,colnames(Q))])>0)
        names(dt)[length(dt)]<-rstate(Q[s,][-match(s,colnames(Q))]/sum(Q[s,][-match(s,colnames(Q))]))
      else names(dt)[length(dt)]<-s
    } else dt[length(dt)]<-dt[length(dt)]-sum(dt)+t
  }
  return(dt)
}

expm.Higham08 <- function (A, balancing = TRUE){
  d <- dim(A)
  if (length(d) != 2 || d[1] != d[2]) 
    stop("'A' must be a square matrix")
  n <- d[1]
  if (n <= 1) 
    return(exp(A))
  if (balancing) {
    baP <- balance(A, "P")
    baS <- balance(baP$z, "S")
    A <- baS$z
  }
  nA <- Matrix::norm(A, "1")
  I <- if (is(A, "Matrix")) 
    Diagonal(n)
  else diag(n)
  if (nA <= 2.1) {
    t <- c(0.015, 0.25, 0.95, 2.1)
    l <- which.max(nA <= t)
    C <- rbind(c(120, 60, 12, 1, 0, 0, 0, 0, 0, 0),
               c(30240, 15120, 3360, 420, 30, 1, 0, 0, 0, 0),
               c(17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1, 0, 0),
               c(17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1))
    A2 <- A %*% A
    P <- I
    U <- C[l, 2] * I
    V <- C[l, 1] * I
    for (k in 1:l) {
      P <- P %*% A2
      U <- U + C[l, (2 * k) + 2] * P
      V <- V + C[l, (2 * k) + 1] * P
    }
    U <- A %*% U
    X <- solve(V - U, V + U)
  }
  else {
    s <- log2(nA/5.4)
    B <- A
    if (s > 0) {
      s <- ceiling(s)
      B <- B/(2^s)
    }
    c. <- c(64764752532480000, 32382376266240000, 7771770303897600, 
            1187353796428800, 129060195264000, 10559470521600, 
            670442572800, 33522128640, 1323241920, 40840800, 
            960960, 16380, 182, 1)
    B2 <- B %*% B
    B4 <- B2 %*% B2
    B6 <- B2 %*% B4
    U <- B %*% (B6 %*% (c.[14] * B6 + c.[12] * B4 + c.[10] * 
                          B2) + c.[8] * B6 + c.[6] * B4 + c.[4] * B2 + c.[2] * 
                  I)
    V <- B6 %*% (c.[13] * B6 + c.[11] * B4 + c.[9] * B2) + 
      c.[7] * B6 + c.[5] * B4 + c.[3] * B2 + c.[1] * I
    X <- solve(V - U, V + U)
    if (s > 0) 
      for (t in 1:s) X <- X %*% X
  }
  if (balancing) {
    d <- baS$scale
    X <- X * (d * rep(1/d, each = n))
    pp <- as.integer(baP$scale)
    if (baP$i1 > 1) {
      for (i in (baP$i1 - 1):1) {
        tt <- X[, i]
        X[, i] <- X[, pp[i]]
        X[, pp[i]] <- tt
        tt <- X[i, ]
        X[i, ] <- X[pp[i], ]
        X[pp[i], ] <- tt
      }
    }
    if (baP$i2 < n) {
      for (i in (baP$i2 + 1):n) {
        tt <- X[, i]
        X[, i] <- X[, pp[i]]
        X[, pp[i]] <- tt
        tt <- X[i, ]
        X[i, ] <- X[pp[i], ]
        X[pp[i], ] <- tt
      }
    }
  }
  X
}

balance <- function (A, job = c("B", "N", "P", "S")){
  .Call("R_dgebal", if (is.numeric(A)) A else as(A, "matrix"), 
        match.arg(job))
}

