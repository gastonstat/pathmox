#' @title Basic PLS-PM algorithm
#' 
#' @description
#' Internal function. \code{get_pls_basic} is called by \code{pathmox}, \code{techmox},
#' \code{fix.pathmox}, \code{fix.techmox}, \code{treemox.pls}, and \code{treemox.boot}. 
#' This function uses internal functions of the super cool package \code{plspm}
#' 
#' @param DT data table
#' @param IDM inner design matrix
#' @param blocks blocks of manifest variables
#' @param modes vector of measurement modes
#' @param scheme inner weighting scheme
#' @param scaled logical indicating whether to scale the data
#' @param tol convergence threshold 
#' @param iter maximum number of iterations
#' @export
#' @keywords internal
get_pls_basic <- 
  function(DT, IDM, blocks, modes, scheme, scaled, tol, iter)
  {   
    ##### PARAMETERS
    DM <- DT
    lvs <- nrow(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
      blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo = rowSums(IDM)
    endo[endo != 0] = 1
    plsr <- FALSE
    ### Variable Names
    lvs.names = colnames(IDM)
    mvs.names = colnames(DM)
    
    # apply the selected scaling
    if (scaled) {
      sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
      X <- scale(DM, scale=sd.X)
    } else {
      X <- scale(DM, scale=FALSE)
    }
    dimnames(X) <- list(rownames(DM), mvs.names)
    
    # ==================== Stage 1: Iterative procedure ==================
    out.ws <- get_weights(X, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(out.ws)) stop("The pls algorithm is non convergent") 
    out.weights <- round(out.ws[[1]], 4)
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA, lvs)
    for (k in 1:lvs) 
      w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Z.lvs <- X %*% out.ws[[2]] %*% diag(w.sig, lvs, lvs)
    Y.lvs <- Z.lvs
    if (!scaled) 
      Y.lvs <- DM %*% out.ws[[2]] %*% diag(w.sig, lvs, lvs)
    dimnames(Y.lvs) <- list(rownames(X), lvs.names)
    dimnames(Z.lvs) <- list(rownames(X), lvs.names)
    loadcomu <- get_loads(X, Y.lvs, blocks)    
    loads <- loadcomu[[1]]
    
    # ============ Stage 2: Path coefficients and total effects ==========
    pathmod <- get_paths(IDM, Y.lvs, plsr)
    innmod <- pathmod[[1]]
    Path <- pathmod[[2]]
    R2 <- pathmod[[3]]
    residuals <- pathmod[[4]]
    
    # ============================= Results ==============================
    model <- list(IDM=IDM, blocks=blocks, scheme=scheme, modes=modes, 
                  scaled=scaled, tol=tol, iter=iter)
    resul = list(out.weights=out.weights, loadings=loads, scores=Y.lvs,   
                 path.coefs=Path, R2=R2, residuals=residuals, model=model)
    resul
  }