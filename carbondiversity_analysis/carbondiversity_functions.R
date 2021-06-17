get_specZ_draws_expand <- function(niter, draws, z_info, pointid, spatial_effect = T){
  Zdata <- get_Z_probs(draws, 1, z_info, spatial_effect = spatial_effect)
  Zprobs <- cbind(Zdata$Z_prob, do.call(cbind, lapply(2:niter, function(x) get_Z_probs(draws, x, z_info)$Z_prob)))
  
  missing_comb <- cbind(subset(as.data.frame(table(Zdata$point, Zdata$species)), Freq == 0))
  expand_comb <- cbind(missing_comb[,c(1,2)], matrix(0, nrow = nrow(missing_comb), ncol = ncol(Zdata)-2))
  names(expand_comb) <- names(Zdata)
  
  Zdata_expand <- rbind(Zdata, expand_comb)
  Zdata_expand <- dplyr::left_join(Zdata_expand, pointid, by = c("point" = "SiteCode"))
  Zprobs_expand <- rbind(Zprobs, matrix(0, nrow = nrow(expand_comb), ncol = ncol(Zprobs)))
  
  Zdraws <- cbind(Zdata_expand[,c("SiteNum","species")],
                  matrix(rbinom(nrow(Zprobs_expand) * niter, 1, as.vector(Zprobs_expand)), nrow = nrow(Zprobs_expand)))
  
  Zdraws <- Zdraws[order(Zdraws$SiteNum),]
  
  return(Zdraws)
}

get_pasture_forest_psi <- function(niter = 10, draws, z_info){
  psicomp <- lapply(1:niter, function(x) get_psi_components(draws, x, z_info))
  pasture_psi <- data.frame(id_sp = psicomp[[1]][,1], do.call(cbind, 
                                                              lapply(psicomp, function(x) exp(x[,2] + x[,3]) / (exp(x[,2] + x[,3]) + 1))))
  forest_psi <- data.frame(id_sp = psicomp[[1]][,1], do.call(cbind, 
                                                             lapply(psicomp, function(x) exp(x[,2] - x[,3]) / (exp(x[,2] - x[,3]) + 1))))
  return(list(pasture_psi = pasture_psi, forest_psi = forest_psi))
}











test <- function (x, plot.layout = c(2, 2), point.cols = "blue", alpha = 0.5, line.cols = "black", ppch = 20, pcex = 0.25,
                  plot.linewidth = 2, ...) {
  options(warn.FPU = FALSE)
  PSAMPLE <- 200
  preddata <- rep(0, times = PSAMPLE)
  thisplot <- 0
  one_page_per_plot <- FALSE
  if((plot.layout[1] == 1) && (plot.layout[2] == 1)){
    one_page_per_plot <- TRUE
  }else{
    par(mfrow = plot.layout)
  }
  
  # Plot observed compositional dissimilarity vs. predicted ecological distance
  plot(unlist(lapply(x, function(y) y$ecological)), unlist(lapply(x, function(y) y$observed)), xlab = "Predicted Ecological Distance", 
       ylab = "Observed Compositional Dissimilarity", 
       type = "n", ylim = c(0, 1))
  for(i in 1:length(x)){
    points(x[[i]]$ecological, x[[i]]$observed, pch = ppch, cex = pcex, col = alpha(point.cols[i], alpha))
    overlayX <- seq(from = min(x[[i]]$ecological), to = max(x[[i]]$ecological), 
                    length = PSAMPLE)
    overlayY <- 1 - exp(-overlayX)
    lines(overlayX, overlayY, lwd = plot.linewidth, col = line.cols[i])
  }
  
  thisplot <- thisplot + 1
  if (one_page_per_plot) {
    dev.new()
    dev.next()
  }
  
  # Plot observed compositional dissimilarity vs. predicted compositional dissimilarity
  plot(unlist(lapply(x, function(y) y$predicted)), unlist(lapply(x, function(y) y$observed)), xlab = "Predicted Compositional Dissimilarity", 
       ylab = "Observed Compositional Dissimilarity", 
       type = "n", ylim = c(0, 1))
  for(i in 1:length(x)){
    points(x[[i]]$predicted, x[[i]]$observed, pch = ppch, cex = pcex, col = alpha(point.cols[i], alpha))
    overlayX <- overlayY <- seq(from = min(x[[i]]$predicted), to = max(x[[i]]$predicted), 
                                length = PSAMPLE)
    lines(overlayX, overlayY, lwd = plot.linewidth, col = line.cols[[i]]) 
  }
  
  thisplot <- thisplot + 1
  
  # for(j in 1:length(x)){
  #   preds <- length(x[[j]]$predictors)
  #   predmax <- 0
  #   splineindex <- 1
  #   z <- list()
  #   for (i in 1:preds){
  #     numsplines <- x[[j]]$splines[i]
  #     if (sum(x[[j]]$coefficients[splineindex:(splineindex + numsplines - 1)]) > 0) {
  #       z[[i]] <- .C("GetPredictorPlotData", pdata = as.double(preddata), 
  #                                            as.integer(PSAMPLE), 
  #                                            as.double(x[[j]]$coefficients[splineindex:(splineindex +  numsplines - 1)]), 
  #                                            as.double(x[[j]]$knots[splineindex:(splineindex + numsplines - 1)]), 
  #                                            as.integer(numsplines), 
  #               PACKAGE = "gdm")
  #       v <- max(z[[i]]$pdata)
  #       if (v > predmax) {
  #         predmax <- v
  #       }
  #     }
  #     splineindex <- splineindex + numsplines
  #   }
  #   splineindex <- 1
  # }
  
  zall <- list()
  for(j in 1:length(x)){
    preds <- x[[j]]$predictors
    z <- list()
    for (i in 1:length(preds)) {
      numsplines <- x[[j]]$splines[i]
      if (sum(x[[j]]$coefficients[splineindex:(splineindex + numsplines - 1)]) > 0) {
        # if (one_page_per_plot) {
        #   dev.new()
        #   dev.next()
        # }else {
        #   thisplot <- thisplot + 1
        #   if (thisplot > (plot.layout[1] * plot.layout[2])) {
        #     thisplot <- 1
        #     par(mfrow = plot.layout)
        #   }
        # }
        # if (x[[j]]$geo & i == 1) {
        #   varNam <- "Geographic Distance"
        # }else {
        #   varNam <- preds[i]
        # }
        
        z[[i]] <- .C("GetPredictorPlotData", pdata = as.double(preddata), 
                     as.integer(PSAMPLE),
                     as.double(x[[j]]$coefficients[splineindex:(splineindex + numsplines - 1)]), 
                     as.double(x[[j]]$knots[splineindex:(splineindex + numsplines - 1)]), 
                     as.integer(numsplines), 
                     PACKAGE = "gdm")
        z[[i]][[length(z[[i]]) + 1]] <- preds[i]
      }
      v <- max(z[[i]]$pdata)
      if (v > predmax) {
        predmax <- v
      }
      splineindex <- splineindex + numsplines
    }
    splineindex <- 1
    zall[[j]] <- z
  }
  
  preds <- unique(unlist(lapply(zall, function(x) lapply(x, function(y) y[[6]]))))
  npreds <- length(preds)
  predlist <- list()
  for(i in 1:npreds){
    predlist[[i]] <- lapply(zall, function(x) if(length(which(unlist(x) == preds[i])) == 1){x[[which(unlist(lapply(x, function(y) y[[6]])) == preds[i])]]}
                            else{NA})
  }
  
  for(i in 1:npreds){
    maxpred <- 0 ; minpred <- 0 ; maxdat <- 0 ; mindat <- 0
    for(j in 1:length(predlist[[i]])){
      maxdat <- ifelse(!is.na(predlist[[i]][[j]]),ifelse(max(predlist[[i]][[j]][[4]]) > maxdat, max(predlist[[i]][[j]][[4]]), maxdat), maxdat)
      mindat <- ifelse(!is.na(predlist[[i]][[j]]),ifelse(min(predlist[[i]][[j]][[4]]) > mindat, min(predlist[[i]][[j]][[4]]), mindat), mindat)
      maxpred <- ifelse(!is.na(predlist[[i]][[j]]),ifelse(max(predlist[[i]][[j]][[1]]) > maxpred, max(predlist[[i]][[j]][[1]]), maxpred), maxpred)
      minpred <- ifelse(!is.na(predlist[[i]][[j]]),ifelse(min(predlist[[i]][[j]][[1]]) > minpred, min(predlist[[i]][[j]][[1]]), minpred), minpred)
    }
    
    plot(c(min(mindat), max(maxdat)), c(min(minpred), max(maxpred)), 
         xlab = preds[i], ylab = paste("f(", preds[i], ")", sep = ""), ylim = c(0, predmax), type = "n")
    for(j in 1:length(predlist[[i]])){
      if(!is.na(predlist[[i]][[j]])){
        points(seq(from = predlist[[i]][[j]][[4]][1], to = predlist[[i]][[j]][[4]][predlist[[i]][[j]][[5]]], length = PSAMPLE),
               predlist[[i]][[j]][[1]], type = 'l', col = point.cols[j])
      }
    }
  }
}