# function

GetDensity <- function(cg.mtbr, kWinSize,ref.length) {
   colnames(cg.mtbr) <- c("chr", "posi", "rC_n", "rC_p", "rT_n", "rT_p")
   posi <- cg.mtbr$posi
   rt <- logical(ref.length)
   rt[posi] <- TRUE
   win <- list(L = as.integer(kWinSize / 2), R = as.integer(kWinSize / 2))
   return(swsCalc(rt, win))
}

GetScore <- function(cg.mtbr, kWinSize, ref.length) {
   ##mtbr score sliding windows
   colnames(cg.mtbr) <- c("chr", "posi", "rC_n", "rC_p", "rT_n", "rT_p")
   
   cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
   cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
   
   rC <- integer(ref.length)
   rC[cg.mtbr$posi] <- cg.mtbr$rC
   rT <- integer(ref.length)
   rT[cg.mtbr$posi] <- cg.mtbr$rT
   win <- list(L = as.integer(kWinSize / 2), R = as.integer(kWinSize / 2))
   rCs <- swsCalc(rC, win)
   rTs <- swsCalc(rT, win)
   score <- rCs/(rCs + rTs)
   score[is.na(score[])] <- 0
   
   return(score)
}

RescaleData <- function(density, score) {
   density.scale <- (density - min(density)) / (max(density) - min(density))
   score.scale <- (score - min(score)) / (max(score) - min(score))
   
   return(data.frame(density = density.scale, score = score.scale))
}

NLCor <- function(vo, vf) {
   sse <- sum((vo - vf)^2)
   sst <- sum((vo - mean(vo))^2)

   return(sqrt(1 - sse / sst))
}

