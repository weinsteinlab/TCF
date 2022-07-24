require(rhdf5)
require(fields)
require(spatstat)
require(ks)
require(scales)
require(stats)

raw.2.bin <- function(data, nbins) {
  n <- length(data)
  vars <- ncol(data[[1]])
  samples <- max(sapply(data, nrow))
  all_data <- do.call(rbind, data)
  all_bin <- bin.data(all_data, nbins)
  return(all_bin)
}

bin.data <- function(data, nbins, min_x=FALSE, max_x=FALSE, min_y = FALSE, max_y = FALSE) {
  if (!min_x) {
    min_x = floor(min(data[,1]))
  }
  if (!max_x) {
    max_x = ceiling(max(data[,1]))
  }
  if (!min_y) {
    min_y = floor(min(data[,2]))
  }
  if (!max_y) {
    max_y = ceiling(max(data[,2]))
  }
  i_states <- seq(min_x, max_x, (max_x-min_x)/(nbins))[1:nbins]
  j_states <- seq(min_y, max_y, (max_y-min_y)/(nbins))[1:nbins]
  t(apply(data, 1, function(x) {c(max(which(x[1] >= i_states)), max(which(x[2] >= j_states)), x[3:ncol(data)])}))
}

# load tica h5 file

load_tica <- function(directory) {
  full <- vector("list", 588)
  for (i in 0:587){
    h5file <- paste(directory, "/tica_files/traj", i, "_on_tica_l20.h5", sep="")
    assignfile <- paste(directory, "/assign_files/assigns_", i, ".txt", sep="")
    tica <- t(h5read(h5file, "arr_0"))
    micro <- read.csv(assignfile, header=FALSE)
    tmp <- cbind(tica[,1:2], micro)
    full[[i+1]] <- tmp
    print(i)
  }  
  full
}

tica_data <- load_tica("~/Desktop/TCF_FOR_GK")
data <- raw.2.bin(tica_data, 100)
weights <- read.table("~/Desktop/TCF_FOR_GK/msm_populations.dat")

bootstrap_sample <- function(tica_data, weights, n_bins) {
   tica_data<- sample(tica_data, length(tica_data), replace=TRUE)  
   data <- raw.2.bin(tica_data, n_bins)
   out <- calculate.tcf.weighted(data, weights, 100, matrix(ncol = 2, nrow = 2, c(1, 0, 0, 1)))
   return(out)
}

calculate.tcf.weighted <- function(data, weights, nbins, varcov) {
  histogram <- matrix(ncol = nbins, nrow = nbins, 0)
  for (i in 1:nrow(data)) {
    print(i)
    histogram[c(data[i,1]), c(data[i,2])] <-  histogram[c(data[i,1]), c(data[i,2])] + weights[(data[i,3]+1),]/length(which(data[,3] == data[i,3]))
  }
  kT_eq <- 0.593*(310/298)
  prob <- histogram
  i_prob <- apply(prob, 1, sum)
  j_prob <- apply(prob, 2, sum)
  TCF <- matrix(ncol = nbins, nrow = nbins, 0)
  for (i in 1:nrow(TCF)){
    for (j in 1:ncol(TCF)) {
      TCF[i,j] <- -kT_eq*(log(prob[i,j]) - log(i_prob[i]) - log(j_prob[j]))
    }
  }
  TCF[prob == 0] <- NA
  A <- -log(prob)
  A[A == Inf] <- NA
  list("TCF" = TCF, "Smooth_TCF" = as.matrix(blur(as.im(TCF), varcov = varcov, bleed = FALSE, normalise = TRUE)), "Free_Energy" = as.matrix(blur(as.im(A), varcov = varcov, bleed = FALSE, normalise = TRUE)), "Probability" = prob)
}

trans_TCF <- calculate.tcf.weighted(data, weights, 100, matrix(ncol = 2, nrow = 2, c(1, 0, 0, 1)))

plot.TCF <- function(x, ax1, ax2, name1, name2, zlim=FALSE) {
  if (zlim == FALSE) {
      zlim = c(-max(abs(x), na.rm = TRUE), max(abs(x), na.rm = TRUE))
  }
  par(mar=c(5,5,5,7))
  image(x, col = colorRampPalette(c("blue", "white", "red"))(100), zlim = zlim, axes = FALSE, xlab = name1, ylab = name2)
  axis(1, at = seq(0, 1, 1/(length(ax1)-1)), labels = ax1)
  axis(2, at = seq(0, 1, 1/(length(ax2)-1)), labels = ax2)
  image.plot(x, legend.only = TRUE, col = colorRampPalette(c("blue", "white", "red"))(100), zlim = zlim, legend.lab = "kT")
}

boots_100 = vector("list", 100)
for (i in 1:20) {
  boots_100[[i]] = bootstrap_sample(tica_data, weights, 100)
}

add.alpha <- function(COLORS, ALPHA){
  if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
  RGB <- col2rgb(COLORS, alpha=TRUE)
  RGB[4,] <- round(RGB[4,]*ALPHA)
  NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
  return(NEW.COLORS)
}

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

# colorRampPaletteAlpha()
colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
  if (interpolate=='linear') {
    l <- approx(a, n=n)
  } else {
    l <- spline(a, n=n)
  }
  l$y[l$y > 255] <- 255 # Clamp if spline is > 255
  cr <- addalpha(cr, l$y/255.0)
  return(cr)
}

filter_boots <- function(boots){
  M <- vector(mode="list", length=length(boots))
  for (i in 1:length(boots)) {
    boot_tmp <- boots[[i]]$Smooth_TCF
    boot_tmp[boot_tmp = NA] <- 0
    M[[i]] <- boot_tmp
  }
  B <- do.call(cbind, M)
  B <- array(B, dim=c(dim(M[[1]]), length(M)))
  ave_v <- apply(B, c(1,2), mean, na.rm=TRUE)
  max_v <- apply(B, c(1,2), quantile, probs=0.95, na.rm=TRUE)
  min_v <- apply(B, c(1,2), quantile, probs=0.05, na.rm=TRUE)
  mask_v <- apply(B, c(1,2), equant, x = 0)
  output <- mask_v
  output[mask_v > 0.5] <- 1 - mask_v[mask_v > 0.5] 
  return(output)
}

equant <- function(data, x) {
  if (sum(is.na(data)) == length(data)){
    return(0)
  } else {
     E <- ecdf(data)
     1 - return(E(x))
  }   
}

boots= vector("list", 20)
for (i in 1:20) {
  boots[[i]] = bootstrap_sample(tica_data, weights, 100)
}

mask = filter_boots(boots)
mask[mask > 0.05 ] <- 1
mask[mask <= 0.05] <- 0
mask = 1 - mask

plot.TCF(trans_TCF$Smooth_TCF, c(2, 4, 6, 8), c(-11, -9, -7, -5, -3), 'tICA 1', 'tICA 2')
plot.TCF(trans_TCF$Smooth_TCF*mask, c(2, 4, 6, 8), c(-11, -9, -7, -5, -3), 'tICA 1', 'tICA 2')
