



source("r/rfuns.r")
testDat = data.table::setDF(data.table::fread("data/optdigits.tes", header = F))
testY = testDat$V65
testDat$V65 = NULL
testX = testDat
  

tmp = data.table::setDF(data.table::fread("data/optdigits.tra", header = F))
# tmp2 = data.table::setDF(data.table::fread("data/optdigits.tes", header = F))
# dat = rbind(tmp, tmp2)
dat = tmp
dat$y = dat$V65
dat$V65 = NULL
tmp = range(unlist(dat[-ncol(dat)]))
for(i in 1:(length(dat) -1L)) dat[[i]] = (dat[[i]] - tmp[1]) / tmp[2]
for(i in 1:(length(testX) -1L)) testX[[i]] = (testX[[i]] - tmp[1]) / tmp[2]
set.seed(123)
dat = dat[sample(nrow(dat)), ]


X = dat[-ncol(dat)]
Y = dat[[ncol(dat)]]
tmp = lapply(aggregate(list(ind = 1:length(Y)), by = list(Y), function(x) x)[[2]], function(x)
{
  t(matrix(colMeans(X[x, ]), nrow = 8))
})
pdf("figure/digitShowAvg.pdf", width = 8, height = 8 * 0.2)
par(mar = c(0.5, 0.5, 0.5, 0.5), family = "serif", mfrow = c(2, 5))
for(i in 1:length(tmp))
{
  image(t(tmp[[i]])[, 8:1], col = colorRampPalette(c("white", "black"))(64), xaxt = "n", yaxt = "n", xlab ="", ylab = "", bty ="n")
}
dev.off()


nullDims = which(unlist(lapply(X, function(x) diff(range(x)) < 1e-10)))
X = X[-nullDims]
testX = testX[-nullDims]
colnames(X) = paste0("V", 1:ncol(X))
colnames(testX) = paste0("V", 1:ncol(testX))
tX = t(X); dimnames(tX) = NULL
ttestX = t(testX); dimnames(ttestX) = NULL




# Kmeans.
if(T)
{
  
  
  # For every k (number of clusters), perform K-means 30 times. Collect all results.
  if(T)
  {
    
    
    set.seed(123)
    rst = list()
    Ks = 3:30
    for(k in Ks)
    {
      cat(k, "")
      rst[[length(rst) + 1]] = list()
      for(i in 1:30)
      {
        rst[[length(rst)]][[i]] = GMKMcharlie::KM(X = tX, centroid = tX[, sample(1:nrow(X), k)], verbose = F, maxCore = 7)
      }
    }
    
    
  }
  
  
  # For every k, select the clustering having the least sum of squares.
  if(T)
  {
    
    
    rstSelected = lapply(rst, function(x)
    {
      tmp = which.min(unlist(lapply(x, function(u)
      {
        sum(unlist(lapply(u, function(v) sum(v$member2centroidDistance ^ 2))))
      })))
      x[[tmp]]
    })
    
    
    clusterings = lapply(rstSelected, function(x) lapply(x, function(y) X[y$clusterMember, ] ))
    clusteringsSS = unlist(lapply(rstSelected, function(x) sum(unlist(lapply(x, function(y) sum(y$member2centroidDistance ^ 2))))))
    
    
    pdf("figure/kmeansSSdigit.pdf", width = 8, height = 8 * 0.618)
    par(mar = c(4.1, 5, 0, 0), family = "serif")
    plot(x = Ks, y = clusteringsSS, lwd = 2, col = "darkblue", bty = "L", xlab = "K", ylab = "In-cluster sum of squares", cex.lab = 2, cex.axis = 1.5, type = "l")
    legend("top", bty = "n", legend = "Digit data", cex = 2)
    dev.off()
    
    
    pdf("figure/kmeansSSDiffDigit.pdf", width = 8, height = 8 * 0.618)
    par(mar = c(4.1, 5, 0, 0), family = "serif")
    plot(x = Ks[1:length(diff(clusteringsSS))], y = diff(clusteringsSS), lwd = 2, col = "darkblue", bty = "L", xlab = "K", ylab = "Sum of square difference", cex.lab = 2, cex.axis = 1.5, type = "l")
    # tmp = Ks[which.min(clusteringsSS)]
    lines(x = c(13, 13), y = c(-1e100, 1e100), type = "l", lty = 2)
    legend("left", bty = "n", legend = "First dip, optimal K = 12", cex = 2)
    legend("right", bty = "n", legend = "Digit data", cex = 2)
    dev.off()
    
    
    bestClustering = rstSelected[[which(unlist(lapply(rstSelected, function(x) length(x))) == 12L)]]
    bestClusteringYs = lapply(bestClustering, function(x) Y[x$clusterMember])
    bestClusteringYs = bestClusteringYs[order(unlist(lapply(bestClusteringYs, function(tmp)
    {
      sort(unique(tmp))[which.max(table(tmp))]
    })))]
    
    
    pdf("figure/kmeansDigit12clusters.pdf", width = 10, height = 10 * 0.3)
    par(mar = c(3, 4, 0, 0), family = "serif", mfrow = c(3, 4))
    tmp = c(7, 3, 6, 2, 5, 4, 0, 1, 9, 9, 8, 1)
    for(i in 1:length(bestClusteringYs))
    {
      tmp = bestClusteringYs[[i]]
      hist(tmp, breaks = seq(-0.5, by = 1, len = 11), col = scales::alpha("darkblue", 0.5), xlab = "", ylab = "", xaxt = "n", border = NA, main = "")
      axis(side = 1, at = 0:9, labels = 0:9, cex.axis = 1.25)
      legend("center", legend = sort(unique(tmp))[which.max(table(tmp))], cex = 2, bty = "n", text.col = "red")
    }
    dev.off()
    
    
    bestClustering = rstSelected[[which(unlist(lapply(rstSelected, function(x) length(x))) == 10L)]]
    bestClusteringYs = lapply(bestClustering, function(x) Y[x$clusterMember])
    bestClusteringYs = bestClusteringYs[order(unlist(lapply(bestClusteringYs, function(tmp)
    {
      sort(unique(tmp))[which.max(table(tmp))]
    })))]
    pdf("figure/kmeansDigit10clusters.pdf", width = 10, height = 10 * 0.3)
    par(mar = c(3, 4, 0, 0), family = "serif", mfrow = c(3, 4))
    tmp = c(7, 3, 6, 2, 5, 4, 0, 1, 9, 9, 8, 1)
    for(i in 1:length(bestClusteringYs))
    {
      tmp = bestClusteringYs[[i]]
      hist(tmp, breaks = seq(-0.5, by = 1, len = 11), col = scales::alpha("darkblue", 0.5), xlab = "", ylab = "", xaxt = "n", border = NA, main = "")
      axis(side = 1, at = 0:9, labels = 0:9, cex.axis = 1.25)
      legend("center", legend = sort(unique(tmp))[which.max(table(tmp))], cex = 2, bty = "n", text.col = "red")
    }
    dev.off()
    
    
  }
  
}




# GMM
if(T)
{
  
  
  library(mclust)
  dir.create("gauResult", showWarnings = F)
  data.table::fwrite(Xgau, file = "data/digitXgau.csv", col.names = F)
  
  
  Ks = 3:20
  cmdlist = paste0("python python/gmm.py data/digitXgau.csv ", Ks, " gauResult/digitXgau-", Ks, ".csv")
  for(x in cmdlist) system(x, wait = F)
  
  
  rstNames = paste0("gauResult/digitXgau-", Ks, ".csv")
  gmmRst = lapply(rstNames, function(x) 
  {
    x = as.numeric(data.table::fread(x, header = F)[[1]])
    list(bic = x[1], lab = as.integer(x[-1]))
  })
  
  
  bics = unlist(lapply(gmmRst, function(x) x$bic))
  pdf("figure/gmmRstBicByNgau.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(x = Ks, y = bics, lwd = 2, col = "darkblue", bty = "L", xlab = "K", ylab = "BIC", cex.lab = 2, cex.axis = 1.5, type = "l")
  legend("top", bty = "n", legend = "Digit data", cex = 2)
  dev.off()
  
  
  bestClustering = gmmRst[[which.min(bics)]]
  tmp = aggregate(list(ind = 1:nrow(X)), by = list(bestClustering$lab), function(x) x)[[2]]
  bestClusteringYs = lapply(tmp, function(x) Y[x])
  bestClusteringYs = bestClusteringYs[order(unlist(lapply(bestClusteringYs, function(tmp)
  {
    sort(unique(tmp))[which.max(table(tmp))]
  })))]
  
  
  pdf("figure/gmmsDigit9clusters.pdf", width = 10, height = 10 * 0.3)
  par(mar = c(3, 4, 0, 0), family = "serif", mfrow = c(3, 3))
  tmp = c(7, 3, 6, 2, 5, 4, 0, 1, 9, 9, 8, 1)
  for(i in 1:length(bestClusteringYs))
  {
    tmp = bestClusteringYs[[i]]
    hist(tmp, breaks = seq(-0.5, by = 1, len = 11), col = scales::alpha("darkblue", 0.5), xlab = "", ylab = "", xaxt = "n", border = NA, main = "")
    axis(side = 1, at = 0:9, labels = 0:9, cex.axis = 1.25)
    legend("center", legend = sort(unique(tmp))[which.max(table(tmp))], cex = 2, bty = "n", text.col = "red")
  }
  dev.off()
}




# PCA
if(T)
{
  
  
  XcolMeans = colMeans(X)
  Xcentered = apply(X, 2, function(x) x - mean(x))
  Xsvd = svd(Xcentered)
  tmp = t(t(testX) - XcolMeans); dimnames(tmp) = NULL
  # svdRst$v * x = test, what is x?
  testXeigenCoor = t(solve(Xsvd$v) %*% t(testX))
  
  
  pdf("figure/digitPCA.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  d = rev(cumsum(rev(Xsvd$d ^ 2)))
  d = d / d[1]
  plot(d, xlab = "Eigen index", ylab = "Cumulative variance", cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  
  
  pdf("figure/digitPCAplotIn2dEigenSpace.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(Xsvd$u[, 1:2] %*% diag(Xsvd$d[1:2]), xlab = "x", ylab = "y", cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  
  
  pdf("figure/digitPCAplotIn3dEigenSpace.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(0, 0, 0, 1.25), family = "serif")
  tmp2 = Xsvd$u[, 1:3] %*% diag(Xsvd$d[1:3])
  plot3D::points3D(x = tmp2[,1], y = tmp2[,2], z = tmp2[,3], colkey = NULL)
  dev.off()
  
  
  
  
  # Take the top compoents that explain 90% of the variance and run NN
  if(T)
  {
    
    
    threshold = 0.9
    tmp = which(d >= 1 - threshold)
    tmp = c(tmp, tmp[length(tmp)])
    # tmp = 1:10
    
    
    Xnew = Xsvd$u[, tmp] %*% diag(Xsvd$d[tmp])
    testXnew = testXeigenCoor[, tmp]
    
    
    model = keras::keras_model_sequential()
    keras::layer_flatten(model, input_shape = ncol(Xnew))
    
    
    keras::layer_dense(model, units = 128, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, units = 64, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, 10, activation = "softmax")
    
    
    opt = keras::optimizer_adam(lr = 0.0005)
    keras::compile(model, loss = "sparse_categorical_crossentropy", optimizer = opt, metrics = "accuracy")
    system.time({fm = keras::fit(model, x = Xnew, y = Y, batch_size = 32, epochs = 100, validation_split = 0.3, verbose = 2, shuffle = T, validation_data = list(x_val = testXnew, y_val = testY))}) # https://keras.rstudio.com/reference/fit.html
    # fm
    # lapply(fm$metrics$loss)
    # Plot
    if(T)
    {
      
      
      lossRange = range(c(range(fm$metrics$loss), range(fm$metrics$val_loss)))
      # accRange = range(c(range(fm$metrics$accuracy), range(fm$metrics$val_accuracy)))
      accRange = c(0, 1)
      lossYaxis = seq(lossRange[1], lossRange[2], len = 5)
      accYaxis = seq(accRange[1], accRange[2], len = 5)
      tmp = fm$metrics
      tmp$loss = (tmp$loss - lossRange[1]) / diff(lossRange)
      tmp$val_loss = (tmp$val_loss - lossRange[1]) / diff(lossRange)
      tmp$accuracy = (tmp$accuracy - accRange[1]) / diff(accRange)
      tmp$val_accuracy = (tmp$val_accuracy - accRange[1]) / diff(accRange)
      ylim = range(unlist(tmp))
      pdf("figure/digitPca10compNeuralnet.pdf", width = 8, height = 8 * 0.618)
      par(mar = c(4.1, 5, 1, 5), family = "serif")
      plot(tmp$loss, type = "l", lwd = 2, col = "darkblue", yaxt = "n", xlab = "Epoch", ylab = "Loss", cex.lab = 2, cex.axis = 1.5, bty = "n", ylim = ylim)
      axis(side = 2, at = seq(0, 1, len = 5), labels = round(lossYaxis, 2), cex.axis = 1.5)
      axis(side = 4, at = seq(0, 1, len = 5), labels = round(accYaxis, 2), cex.axis = 1.5)
      mtext("Accuracy", side = 4, line = 3, cex = 2)
      lines(tmp$val_loss, type = "l", lwd = 2, col = "red")
      lines(tmp$accuracy, type = "l", lwd = 2, col = "skyblue", lty = 1)
      lines(tmp$val_accuracy, type = "l", lwd = 2, col = "olivedrab3", lty = 1)
      legend("right", legend = c("Training loss", "Validation loss", "Training accuracy", "Validation accuracy"), lwd = c(2, 2, 2, 2), col = c("darkblue", "red", "skyblue", "olivedrab3"), bty = "n", cex = 2)
      dev.off()
      
      
    }
    
    
  }
  
  
  
  
  # Run GMM and use the Gaussian probability matrix as new features.
  if(T)
  {
    
    
    # We tried using log probabilities, but the results are really bad.
    # We tried using high dimensional Gaussians --- use 22 PCs, but in the end, the probability matrix mostly have 1s and 0s.
    XcolMeans = colMeans(X)
    Xcentered = apply(X, 2, function(x) x - mean(x))
    Xsvd = svd(Xcentered)
    testXeigenCoor = t(solve(Xsvd$v) %*% t(testX))
    Ncomp = 22L
    Nclust = 10L
    whichPCs2select = 1:Ncomp
    Xnew = Xsvd$u[, whichPCs2select] %*% diag(Xsvd$d[whichPCs2select])
    testXnew = testXeigenCoor[, whichPCs2select]
    
    
    dir.create("gauResult", showWarnings = F)
    trainDataPath = "data/digitXpcaTrain.csv"
    testDataPath = "data/digitXpcaTest.csv"
    trainSavePath = "data/digitXpcaTrainGauScore.csv"
    testSavePath = "data/digitXpcaTestGauScore.csv"
    data.table::fwrite(Xnew, file = trainDataPath, col.names = F)
    data.table::fwrite(testXnew, file = testDataPath, col.names = F)
    system(paste0("python python/gmm.py ", trainDataPath, " ", testDataPath, " ", Nclust, " ", trainSavePath, " ", testSavePath), wait = T)
    
    
    Xnew = as.matrix(data.table::fread(trainSavePath, header = F))
    Xnew = (Xnew - min(Xnew)) / max(Xnew)
    testXnew = as.matrix(data.table::fread(testSavePath, header = F))
    testXnew = (testXnew - min(Xnew)) / max(Xnew)
    
    
    model = keras::keras_model_sequential()
    keras::layer_flatten(model, input_shape = ncol(Xnew))
    
    
    keras::layer_dense(model, units = 128, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, units = 64, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, 10, activation = "softmax")
    
    
    opt = keras::optimizer_adam(lr = 0.0005)
    keras::compile(model, loss = "sparse_categorical_crossentropy", optimizer = opt, metrics = "accuracy")
    system.time({fm = keras::fit(model, x = Xnew, y = Y, batch_size = 32, epochs = 100, validation_split = 0.3, verbose = 2, shuffle = T, validation_data = list(x_val = testXnew, y_val = testY))}) # https://keras.rstudio.com/reference/fit.html
    # fm
    # lapply(fm$metrics$loss)
    # Plot
    if(T)
    {
      
      
      lossRange = range(c(range(fm$metrics$loss), range(fm$metrics$val_loss)))
      # accRange = range(c(range(fm$metrics$accuracy), range(fm$metrics$val_accuracy)))
      accRange = c(0, 1)
      lossYaxis = seq(lossRange[1], lossRange[2], len = 5)
      accYaxis = seq(accRange[1], accRange[2], len = 5)
      tmp = fm$metrics
      tmp$loss = (tmp$loss - lossRange[1]) / diff(lossRange)
      tmp$val_loss = (tmp$val_loss - lossRange[1]) / diff(lossRange)
      tmp$accuracy = (tmp$accuracy - accRange[1]) / diff(accRange)
      tmp$val_accuracy = (tmp$val_accuracy - accRange[1]) / diff(accRange)
      ylim = range(unlist(tmp))
      pdf("figure/digitPca10compClusterNeuralnet.pdf", width = 8, height = 8 * 0.618)
      par(mar = c(4.1, 5, 1, 5), family = "serif")
      plot(tmp$loss, type = "l", lwd = 2, col = "darkblue", yaxt = "n", xlab = "Epoch", ylab = "Loss", cex.lab = 2, cex.axis = 1.5, bty = "n", ylim = ylim)
      axis(side = 2, at = seq(0, 1, len = 5), labels = round(lossYaxis, 2), cex.axis = 1.5)
      axis(side = 4, at = seq(0, 1, len = 5), labels = round(accYaxis, 2), cex.axis = 1.5)
      mtext("Accuracy", side = 4, line = 3, cex = 2)
      lines(tmp$val_loss, type = "l", lwd = 2, col = "red")
      lines(tmp$accuracy, type = "l", lwd = 2, col = "skyblue", lty = 1)
      lines(tmp$val_accuracy, type = "l", lwd = 2, col = "olivedrab3", lty = 1)
      legend("right", legend = c("Training loss", "Validation loss", "Training accuracy", "Validation accuracy"), lwd = c(2, 2, 2, 2), col = c("darkblue", "red", "skyblue", "olivedrab3"), bty = "n", cex = 2)
      dev.off()
    }
    
    
  }
  
  
}




# ICA
if(T)
{
  
  
  set.seed(42)
  icaRst = list()
  for(i in 2:ncol(X))
  {
    cat(i, "")
    icaRst[[length(icaRst) + 1]] = fastICA::fastICA(X, n.comp = i, alg.typ = c("deflation"), fun = c("logcosh"), alpha = 1.0, method = c("C"), row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL)
  }
  for(i in 2:length(icaRst)) icaRst[[i]]$X = NULL
  # save(icaRst, file = "data/icaRst.Rdata")
  kurtoAll = lapply(icaRst, function(x) kurto(x$S))
  abskurtoAll = as.numeric(lapply(kurtoAll, function(x) mean(abs(x))))
  
  
  # Since we don't know which components are important. 
  pdf("figure/digitDataKurtosis.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(abskurtoAll, type = "l", col = "darkblue", xlab = "Number of components", ylab = "Mean absolute excess kurtosis", bty = "L", cex.axis = 1.5, cex.lab= 2, lwd = 2)
  dev.off()
  
  
  # Try clustering on 10 independent components, in the hope that each image (observation) can be fully characterized by these 10 features.
  if(T)
  {
    
    
    tmp = icaRst[[which(unlist(lapply(icaRst, function(x) ncol(x$K))) == 10L)]]
    Xnew = tmp$S
    tXnew = t(Xnew)
    set.seed(123)
    rst = list()
    k = 10
    for(i in 1:100)
    {
      cat(i, "")
      rst[[length(rst) + 1]] = GMKMcharlie::KM(X = tXnew, centroid = tXnew[, sample(1:nrow(Xnew), k)], verbose = F, maxCore = 7)
    }
    rst = rst[[which.min(unlist(lapply(rst, function(x) sum(unlist(lapply(x, function(y) sum(y$member2centroidDistance ^ 2)))))))]]
    
    
    # bestClustering = rstSelected[[which(unlist(lapply(rstSelected, function(x) length(x))) == 10L)]]
    bestClusteringYs = lapply(rst, function(x) Y[x$clusterMember])
    bestClusteringYs = bestClusteringYs[order(unlist(lapply(bestClusteringYs, function(tmp)
    {
      sort(unique(tmp))[which.max(table(tmp))]
    })))]
    pdf("figure/kmeansDigit10clustersICA.pdf", width = 10, height = 10 * 0.3)
    par(mar = c(3, 4, 0, 0), family = "serif", mfrow = c(3, 4))
    # tmp = c(7, 3, 6, 2, 5, 4, 0, 1, 9, 9, 8, 1)
    for(i in 1:length(bestClusteringYs))
    {
      tmp = bestClusteringYs[[i]]
      hist(tmp, breaks = seq(-0.5, by = 1, len = 11), col = scales::alpha("darkblue", 0.5), xlab = "", ylab = "", xaxt = "n", border = NA, main = "")
      axis(side = 1, at = 0:9, labels = 0:9, cex.axis = 1.25)
      legend("center", legend = sort(unique(tmp))[which.max(table(tmp))], cex = 2, bty = "n", text.col = "red")
    }
    dev.off()
    
  }
  
  
  
  
  # Take 22 or 10 components and run NN. 22 or 10 is to match what we did in PCA.
  if(T)
  {
    
    
    Ncomp = 22L
    tmp = icaRst[[Ncomp - 1]]
    Xnew = tmp$S
    tXnew = t(Xnew)
    testXnew = as.matrix(testX) %*% tmp$K %*% tmp$W
    ttestXnew = t(testXnew)
    
    
    model = keras::keras_model_sequential()
    keras::layer_flatten(model, input_shape = ncol(Xnew))
    
    
    keras::layer_dense(model, units = 128, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, units = 64, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, 10, activation = "softmax")
    
    
    opt = keras::optimizer_adam(lr = 0.0005)
    keras::compile(model, loss = "sparse_categorical_crossentropy", optimizer = opt, metrics = "accuracy")
    system.time({fm = keras::fit(model, x = Xnew, y = Y, batch_size = 32, epochs = 100, validation_split = 0.3, verbose = 2, shuffle = T, validation_data = list(x_val = testXnew, y_val = testY))}) # https://keras.rstudio.com/reference/fit.html
    # fm
    # lapply(fm$metrics$loss)
    # Plot
    if(T)
    {
      
      
      lossRange = range(c(range(fm$metrics$loss), range(fm$metrics$val_loss)))
      # accRange = range(c(range(fm$metrics$accuracy), range(fm$metrics$val_accuracy)))
      accRange = c(0, 1)
      lossYaxis = seq(lossRange[1], lossRange[2], len = 5)
      accYaxis = seq(accRange[1], accRange[2], len = 5)
      tmp = fm$metrics
      tmp$loss = (tmp$loss - lossRange[1]) / diff(lossRange)
      tmp$val_loss = (tmp$val_loss - lossRange[1]) / diff(lossRange)
      tmp$accuracy = (tmp$accuracy - accRange[1]) / diff(accRange)
      tmp$val_accuracy = (tmp$val_accuracy - accRange[1]) / diff(accRange)
      ylim = range(unlist(tmp))
      pdf("figure/digitIca22compNeuralnet.pdf", width = 8, height = 8 * 0.618)
      par(mar = c(4.1, 5, 1, 5), family = "serif")
      plot(tmp$loss, type = "l", lwd = 2, col = "darkblue", yaxt = "n", xlab = "Epoch", ylab = "Loss", cex.lab = 2, cex.axis = 1.5, bty = "n")
      axis(side = 2, at = seq(0, 1, len = 5), labels = round(lossYaxis, 2), cex.axis = 1.5)
      axis(side = 4, at = seq(0, 1, len = 5), labels = round(accYaxis, 2), cex.axis = 1.5)
      mtext("Accuracy", side = 4, line = 3, cex = 2)
      lines(tmp$val_loss, type = "l", lwd = 2, col = "red")
      lines(tmp$accuracy, type = "l", lwd = 2, col = "skyblue", lty = 1)
      lines(tmp$val_accuracy, type = "l", lwd = 2, col = "olivedrab3", lty = 1)
      legend("right", legend = c("Training loss", "Validation loss", "Training accuracy", "Validation accuracy"), lwd = c(2, 2, 2, 2), col = c("darkblue", "red", "skyblue", "olivedrab3"), bty = "n", cex = 2)
      dev.off()
      
      
    }  
    
  }
  
  
  
  
  # Run GMM and use the Gaussian probability matrix as new features.
  if(T)
  {
    
    
    # XcolMeans = colMeans(X)
    # Xcentered = apply(X, 2, function(x) x - mean(x))
    # Xsvd = svd(Xcentered)
    # testXeigenCoor = t(solve(Xsvd$v) %*% t(testX))
    # whichPCs2select = 1:10
    # Xnew = Xsvd$u[, whichPCs2select] %*% diag(Xsvd$d[whichPCs2select])
    # testXnew = testXeigenCoor[, whichPCs2select]
    
    
    Ncomp = 10L
    Nclust = 10L
    tmp = icaRst[[Ncomp - 1]]
    Xnew = tmp$S
    testXnew = as.matrix(testX) %*% tmp$K %*% tmp$W
    
    
    dir.create("gauResult", showWarnings = F)
    trainDataPath = "data/digitXicaTrain.csv"
    testDataPath = "data/digitXicaTest.csv"
    trainSavePath = "data/digitXicaTrainGauScore.csv"
    testSavePath = "data/digitXicaTestGauScore.csv"
    data.table::fwrite(Xnew, file = trainDataPath, col.names = F)
    data.table::fwrite(testXnew, file = testDataPath, col.names = F)
    system(paste0("python python/gmm.py ", trainDataPath, " ", testDataPath, " ", Nclust, " ", trainSavePath, " ", testSavePath), wait = T)
    
    
    Xnew = as.matrix(data.table::fread(trainSavePath, header = F))
    Xnew = (Xnew - min(Xnew)) / max(Xnew)
    testXnew = as.matrix(data.table::fread(testSavePath, header = F))
    testXnew = (testXnew - min(Xnew)) / max(Xnew)
    
    
    model = keras::keras_model_sequential()
    keras::layer_flatten(model, input_shape = ncol(Xnew))
    
    
    keras::layer_dense(model, units = 128, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, units = 64, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, 10, activation = "softmax")
    
    
    opt = keras::optimizer_adam(lr = 0.0005)
    keras::compile(model, loss = "sparse_categorical_crossentropy", optimizer = opt, metrics = "accuracy")
    system.time({fm = keras::fit(model, x = Xnew, y = Y, batch_size = 32, epochs = 100, validation_split = 0.3, verbose = 2, shuffle = T, validation_data = list(x_val = testXnew, y_val = testY))}) # https://keras.rstudio.com/reference/fit.html
    # fm
    # lapply(fm$metrics$loss)
    # Plot
    if(T)
    {
      
      
      lossRange = range(c(range(fm$metrics$loss), range(fm$metrics$val_loss)))
      # accRange = range(c(range(fm$metrics$accuracy), range(fm$metrics$val_accuracy)))
      accRange = c(0, 1)
      lossYaxis = seq(lossRange[1], lossRange[2], len = 5)
      accYaxis = seq(accRange[1], accRange[2], len = 5)
      tmp = fm$metrics
      tmp$loss = (tmp$loss - lossRange[1]) / diff(lossRange)
      tmp$val_loss = (tmp$val_loss - lossRange[1]) / diff(lossRange)
      tmp$accuracy = (tmp$accuracy - accRange[1]) / diff(accRange)
      tmp$val_accuracy = (tmp$val_accuracy - accRange[1]) / diff(accRange)
      ylim = range(unlist(tmp))
      pdf("figure/digitIca10compClusterNeuralnet.pdf", width = 8, height = 8 * 0.618)
      par(mar = c(4.1, 5, 1, 5), family = "serif")
      plot(tmp$loss, type = "l", lwd = 2, col = "darkblue", yaxt = "n", xlab = "Epoch", ylab = "Loss", cex.lab = 2, cex.axis = 1.5, bty = "n", ylim = ylim)
      axis(side = 2, at = seq(0, 1, len = 5), labels = round(lossYaxis, 2), cex.axis = 1.5)
      axis(side = 4, at = seq(0, 1, len = 5), labels = round(accYaxis, 2), cex.axis = 1.5)
      mtext("Accuracy", side = 4, line = 3, cex = 2)
      lines(tmp$val_loss, type = "l", lwd = 2, col = "red")
      lines(tmp$accuracy, type = "l", lwd = 2, col = "skyblue", lty = 1)
      lines(tmp$val_accuracy, type = "l", lwd = 2, col = "olivedrab3", lty = 1)
      legend("right", legend = c("Training loss", "Validation loss", "Training accuracy", "Validation accuracy"), lwd = c(2, 2, 2, 2), col = c("darkblue", "red", "skyblue", "olivedrab3"), bty = "n", cex = 2)
      dev.off()
    }
    
    
  }
 
   
}




# Random projection
# Random projection.
# (1 - eps)||u - v||^2 < ||p(u) - p(v)||^2 < (1 + eps)||u - v||^2
findNcomp = function(eps, Ndata) { 4 * log(Ndata) / (eps ^ 2 / 2 - eps ^ 3 / 3) }
findNcomp(1, nrow(X)); findNcomp(0.5, nrow(X)); findNcomp(0.2, nrow(X)); findNcomp(0.1, nrow(X))
# 198; 396; 1903; 7070;
if(T)
{
  
  
  set.seed(42)
  
  
  XorignalPairwiseD = dist(X)
  Ncomponent = c(2L, 10L, 30L, 62L)
  
  
  # Gaussian random projection
  if(T)
  {
    XlowDs = lapply(Ncomponent, function(x) 
    {
      lapply(1:20, function(i)
      {
        R = matrix(rnorm(x * ncol(X)) / x ^ 0.5, ncol = x)
        as.matrix(X) %*% R
      })
    })
  }
  
  
  reconsErr = lapply(XlowDs, function(x) 
  {
    cat(".")
    tmp = as.data.frame(lapply(x, function(u) as.numeric(dist(u) / XorignalPairwiseD)))
    rowMeans(tmp)
  })
  
  
  pdf("figure/digitRPreconErrHist.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif", mfrow = c(2, 2))
  tmp = lapply(reconsErr, function(x) x ^ 2)
  hist(tmp[[1]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "Frequency")
  legend("topright", bty = "n", cex = 2, legend = "N(component) = 2")
  hist(tmp[[2]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "10")
  hist(tmp[[3]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "Projected / original", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "30")
  hist(tmp[[4]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "62")
  dev.off()
  
  
  
  
  # Try K-means one the best lowD embedding.
  if(T)
  {
    
    
    set.seed(42)
    tmp = as.matrix(X)
    XlowDs = lapply(1:50, function(i)
    {
      R = matrix(rnorm(10L * ncol(X)) / 10 ^ 0.5, ncol = 10)
      tmp %*% R
    })
    reconsErr = unlist(lapply(XlowDs, function(x) 
    {
      cat(".")
      mean(abs(dist(x) / XorignalPairwiseD - 1))
    }))
    XlowDs = XlowDs[[which.min(reconsErr)]]
    
    
    Xnew = XlowDs
    tXnew = t(Xnew)
    set.seed(123)
    rst = list()
    k = 10
    for(i in 1:100)
    {
      cat(i, "")
      rst[[length(rst) + 1]] = GMKMcharlie::KM(X = tXnew, centroid = tXnew[, sample(1:nrow(Xnew), k)], verbose = F, maxCore = 7)
    }
    rst = rst[[which.min(unlist(lapply(rst, function(x) sum(unlist(lapply(x, function(y) sum(y$member2centroidDistance ^ 2)))))))]]
    
    
    bestClusteringYs = lapply(rst, function(x) Y[x$clusterMember])
    bestClusteringYs = bestClusteringYs[order(unlist(lapply(bestClusteringYs, function(tmp)
    {
      sort(unique(tmp))[which.max(table(tmp))]
    })))]
    pdf("figure/kmeansDigit10clustersRP10dim.pdf", width = 10, height = 10 * 0.3)
    par(mar = c(3, 4, 0, 0), family = "serif", mfrow = c(3, 4))
    # tmp = c(7, 3, 6, 2, 5, 4, 0, 1, 9, 9, 8, 1)
    for(i in 1:length(bestClusteringYs))
    {
      tmp = bestClusteringYs[[i]]
      hist(tmp, breaks = seq(-0.5, by = 1, len = 11), col = scales::alpha("darkblue", 0.5), xlab = "", ylab = "", xaxt = "n", border = NA, main = "")
      axis(side = 1, at = 0:9, labels = 0:9, cex.axis = 1.25)
      legend("center", legend = sort(unique(tmp))[which.max(table(tmp))], cex = 2, bty = "n", text.col = "red")
    }
    dev.off()
    
  }
  
  
  
  
  # Take 10 components and run NN
  if(T)
  {
    
    
    set.seed(42)
    XorignalPairwiseD = dist(X)
    Ncomponent = c(2L, 10L, 30L, 62L)
    
    
    # Gaussian random projection
    if(T)
    {
      
      
      XlowDs = lapply(Ncomponent, function(x) 
      {
        lapply(1:20, function(i)
        {
          R = matrix(rnorm(x * ncol(X)) / x ^ 0.5, ncol = x)
          X = as.matrix(X) %*% R
          dimnames(X) = NULL
          list(X = X, R = R)
        })
      })
      
      
    }
    
    
    bestLowEmbed = lapply(XlowDs, function(x) 
    {
      cat(".")
      tmp = which.min(unlist(lapply(x, function(u) sum(abs(dist(u$X) - XorignalPairwiseD)))))
      cat(tmp, "")
      x[[tmp]]
    })
    # for(i in 1:length(bestLowEmbed)) dimnames(bestLowEmbed[[i]]) = NULL
    
    
    Ncomp = 10L
    tmp = bestLowEmbed[[which(unlist(lapply(bestLowEmbed, function(x) ncol(x$R) == Ncomp)))]]
    Xnew = tmp$X
    testXnew = as.matrix(testX) %*% tmp$R
    
    
    model = keras::keras_model_sequential()
    keras::layer_flatten(model, input_shape = ncol(Xnew))
    
    
    keras::layer_dense(model, units = 128, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, units = 64, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, 10, activation = "softmax")
    
    
    opt = keras::optimizer_adam(lr = 0.0005)
    keras::compile(model, loss = "sparse_categorical_crossentropy", optimizer = opt, metrics = "accuracy")
    system.time({fm = keras::fit(model, x = Xnew, y = Y, batch_size = 32, epochs = 100, validation_split = 0.3, verbose = 2, shuffle = T, validation_data = list(x_val = testXnew, y_val = testY))}) # https://keras.rstudio.com/reference/fit.html
    # fm
    # lapply(fm$metrics$loss)
    # Plot
    if(T)
    {
      
      
      lossRange = range(c(range(fm$metrics$loss), range(fm$metrics$val_loss)))
      # accRange = range(c(range(fm$metrics$accuracy), range(fm$metrics$val_accuracy)))
      accRange = c(0, 1)
      lossYaxis = seq(lossRange[1], lossRange[2], len = 5)
      accYaxis = seq(accRange[1], accRange[2], len = 5)
      tmp = fm$metrics
      tmp$loss = (tmp$loss - lossRange[1]) / diff(lossRange)
      tmp$val_loss = (tmp$val_loss - lossRange[1]) / diff(lossRange)
      tmp$accuracy = (tmp$accuracy - accRange[1]) / diff(accRange)
      tmp$val_accuracy = (tmp$val_accuracy - accRange[1]) / diff(accRange)
      ylim = range(unlist(tmp))
      pdf("figure/digitRP10compNeuralnet.pdf", width = 8, height = 8 * 0.618)
      par(mar = c(4.1, 5, 1, 5), family = "serif")
      plot(tmp$loss, type = "l", lwd = 2, col = "darkblue", yaxt = "n", xlab = "Epoch", ylab = "Loss", cex.lab = 2, cex.axis = 1.5, bty = "n", ylim = ylim)
      axis(side = 2, at = seq(0, 1, len = 5), labels = round(lossYaxis, 2), cex.axis = 1.5)
      axis(side = 4, at = seq(0, 1, len = 5), labels = round(accYaxis, 2), cex.axis = 1.5)
      mtext("Accuracy", side = 4, line = 3, cex = 2)
      lines(tmp$val_loss, type = "l", lwd = 2, col = "red")
      lines(tmp$accuracy, type = "l", lwd = 2, col = "skyblue", lty = 1)
      lines(tmp$val_accuracy, type = "l", lwd = 2, col = "olivedrab3", lty = 1)
      legend("right", legend = c("Training loss", "Validation loss", "Training accuracy", "Validation accuracy"), lwd = c(2, 2, 2, 2), col = c("darkblue", "red", "skyblue", "olivedrab3"), bty = "n", cex = 2)
      dev.off()
      
      
    }  
    
  }
  
  
  
  
  # Run GMM and use the Gaussian probability matrix as new features.
  if(T)
  {
    
    set.seed(42)
    XorignalPairwiseD = dist(X)
    Ncomponent = c(2L, 10L, 30L, 62L)
    
    
    # Gaussian random projection
    if(T)
    {
      
      
      XlowDs = lapply(Ncomponent, function(x) 
      {
        lapply(1:20, function(i)
        {
          R = matrix(rnorm(x * ncol(X)) / x ^ 0.5, ncol = x)
          X = as.matrix(X) %*% R
          dimnames(X) = NULL
          list(X = X, R = R)
        })
      })
      
      
    }
    
    
    bestLowEmbed = lapply(XlowDs, function(x) 
    {
      cat(".")
      tmp = which.min(unlist(lapply(x, function(u) sum(abs(dist(u$X) - XorignalPairwiseD)))))
      cat(tmp, "")
      x[[tmp]]
    })
    # for(i in 1:length(bestLowEmbed)) dimnames(bestLowEmbed[[i]]) = NULL
    
    
    Ncomp = 10L
    Nclust = 10L
    tmp = bestLowEmbed[[which(unlist(lapply(bestLowEmbed, function(x) ncol(x$R) == Ncomp)))]]
    Xnew = tmp$X
    testXnew = as.matrix(testX) %*% tmp$R
    
    
    dir.create("gauResult", showWarnings = F)
    trainDataPath = "data/digitXrpTrain.csv"
    testDataPath = "data/digitXrpTest.csv"
    trainSavePath = "data/digitXrpTrainGauScore.csv"
    testSavePath = "data/digitXrpTestGauScore.csv"
    data.table::fwrite(Xnew, file = trainDataPath, col.names = F)
    data.table::fwrite(testXnew, file = testDataPath, col.names = F)
    system(paste0("python python/gmm.py ", trainDataPath, " ", testDataPath, " ", Nclust, " ", trainSavePath, " ", testSavePath), wait = T)
    
    
    Xnew = as.matrix(data.table::fread(trainSavePath, header = F))
    Xnew = (Xnew - min(Xnew)) / max(Xnew)
    testXnew = as.matrix(data.table::fread(testSavePath, header = F))
    testXnew = (testXnew - min(Xnew)) / max(Xnew)
    
    
    model = keras::keras_model_sequential()
    keras::layer_flatten(model, input_shape = ncol(Xnew))
    
    
    keras::layer_dense(model, units = 128, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, units = 64, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, 10, activation = "softmax")
    
    
    opt = keras::optimizer_adam(lr = 0.0005)
    keras::compile(model, loss = "sparse_categorical_crossentropy", optimizer = opt, metrics = "accuracy")
    system.time({fm = keras::fit(model, x = Xnew, y = Y, batch_size = 32, epochs = 100, validation_split = 0.3, verbose = 2, shuffle = T, validation_data = list(x_val = testXnew, y_val = testY))}) # https://keras.rstudio.com/reference/fit.html
    # fm
    # lapply(fm$metrics$loss)
    # Plot
    if(T)
    {
      
      
      lossRange = range(c(range(fm$metrics$loss), range(fm$metrics$val_loss)))
      # accRange = range(c(range(fm$metrics$accuracy), range(fm$metrics$val_accuracy)))
      accRange = c(0, 1)
      lossYaxis = seq(lossRange[1], lossRange[2], len = 5)
      accYaxis = seq(accRange[1], accRange[2], len = 5)
      tmp = fm$metrics
      tmp$loss = (tmp$loss - lossRange[1]) / diff(lossRange)
      tmp$val_loss = (tmp$val_loss - lossRange[1]) / diff(lossRange)
      tmp$accuracy = (tmp$accuracy - accRange[1]) / diff(accRange)
      tmp$val_accuracy = (tmp$val_accuracy - accRange[1]) / diff(accRange)
      ylim = range(unlist(tmp))
      pdf("figure/digitRp10compClusterNeuralnet.pdf", width = 8, height = 8 * 0.618)
      par(mar = c(4.1, 5, 1, 5), family = "serif")
      plot(tmp$loss, type = "l", lwd = 2, col = "darkblue", yaxt = "n", xlab = "Epoch", ylab = "Loss", cex.lab = 2, cex.axis = 1.5, bty = "n", ylim = ylim)
      axis(side = 2, at = seq(0, 1, len = 5), labels = round(lossYaxis, 2), cex.axis = 1.5)
      axis(side = 4, at = seq(0, 1, len = 5), labels = round(accYaxis, 2), cex.axis = 1.5)
      mtext("Accuracy", side = 4, line = 3, cex = 2)
      lines(tmp$val_loss, type = "l", lwd = 2, col = "red")
      lines(tmp$accuracy, type = "l", lwd = 2, col = "skyblue", lty = 1)
      lines(tmp$val_accuracy, type = "l", lwd = 2, col = "olivedrab3", lty = 1)
      legend("right", legend = c("Training loss", "Validation loss", "Training accuracy", "Validation accuracy"), lwd = c(2, 2, 2, 2), col = c("darkblue", "red", "skyblue", "olivedrab3"), bty = "n", cex = 2)
      dev.off()
      
      
    }
    
    
  }
  
}




# Diffusion map
if(T)
{
  
  # Change eps, run, save.
  if(F)
  {
    
    
    dmat = as.matrix(dist(X))
    kernelDmat = exp(-dmat ^ 2 * 2)
    P = kernelDmat / colSums(kernelDmat) # Symmetric
    M = P
    
    
    Rcpp::sourceCpp("src/speedup.cpp", verbose = T)
    rst = list()
    for(i in 1:1000)
    {
      lastM = M
      M = matmul(M, M, maxCore = 15) 
      M = M / rowSums(M)
      meanErr = mean(abs(M - lastM))
      cat(i, ", err = ", meanErr, "    ")
      rst[[length(rst) + 1]] = M
      if(meanErr < 1e-17) break 
    }
    stationaryP = rst[[length(rst)]]
    M1 = P
    M1 = apply(M1, 2, function(x) x - mean(x))
    M1svd = svd(M1 / stationaryP ^ 0.5)
    M2 = rst[[1]]
    M2 = apply(M2, 2, function(x) x - mean(x))
    M2svd = svd(M2 / stationaryP ^ 0.5)
    M4 = rst[[2]]
    M4 = apply(M4, 2, function(x) x - mean(x))
    M4svd = svd(M4 / stationaryP ^ 0.5)
    # save(M1svd, file = "data/digit-eps-2.0-M1svd.Rdata")
    # save(M2svd, file = "data/digit-eps-2.0-M2svd.Rdata")
    # save(M4svd, file = "data/digit-eps-2.0-M4svd.Rdata")
  }
  
  
  
  
  load("data/digit-eps-0.5-M4svd.Rdata")
  tmp = aggregate(list(ind = 1:length(Y)), by = list(Y), function(x) x)
  palette = randomcoloR::distinctColorPalette(length(unique(Y)))
  tmpcols = unlist(mapply(function(x, y) rep(x, length(y)), palette, tmp[[2]]))
  tmpcols = tmpcols[order(unlist(tmp[[2]]))]
  
  
  pdf("figure/digitDiffusionMap2d.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(M4svd$u[,1:2], pch = 16, col = tmpcols, cex = 0.5, bty = "L", xlab = "x", ylab = "y", cex.lab = 2, cex.axis = 1.5)
  dev.off()
  
  

  
  # Take the components that express 70% of the variance.
  tmp = rev(cumsum(rev(M4svd$d ^ 2))); tmp = tmp / tmp[1]
  threshold = 0.99
  tmp = which(tmp >= 1 - threshold)
  Xnew = M4svd$u[, tmp] %*% diag(M4svd$d[tmp])
  tXnew = t(Xnew)
  
  
  # Try K-means one the best lowD embedding.
  if(T)
  {
    
    
    set.seed(123)
    rst = list()
    k = 10
    for(i in 1:100)
    {
      cat(i, "")
      rst[[length(rst) + 1]] = GMKMcharlie::KM(X = tXnew, centroid = tXnew[, sample(1:nrow(Xnew), k)], verbose = F, maxCore = 7)
    }
    rst = rst[[which.min(unlist(lapply(rst, function(x) sum(unlist(lapply(x, function(y) sum(y$member2centroidDistance ^ 2)))))))]]
    
    
    bestClusteringYs = lapply(rst, function(x) Y[x$clusterMember])
    bestClusteringYs = bestClusteringYs[order(unlist(lapply(bestClusteringYs, function(tmp)
    {
      sort(unique(tmp))[which.max(table(tmp))]
    })))]
    pdf("figure/kmeansDigit10clustersDiffusion99percentVariancePreserved.pdf", width = 10, height = 10 * 0.3)
    par(mar = c(3, 4, 0, 0), family = "serif", mfrow = c(3, 4))
    # tmp = c(7, 3, 6, 2, 5, 4, 0, 1, 9, 9, 8, 1)
    for(i in 1:length(bestClusteringYs))
    {
      tmp = bestClusteringYs[[i]]
      hist(tmp, breaks = seq(-0.5, by = 1, len = 11), col = scales::alpha("darkblue", 0.5), xlab = "", ylab = "", xaxt = "n", border = NA, main = "")
      axis(side = 1, at = 0:9, labels = 0:9, cex.axis = 1.25)
      legend("center", legend = sort(unique(tmp))[which.max(table(tmp))], cex = 2, bty = "n", text.col = "red")
    }
    dev.off()
    
    
  }
  
  
}






# T-stne
if(T)
{
  
  
  tsneRst = list()
  for(d in 2:3)
  {
    cat(d, "")
    tsneRst[[length(tsneRst) + 1]] = Rtsne::Rtsne(as.matrix(X), dims = d, initial_dims = 50, perplexity = 30, theta = 0.5, check_duplicates = TRUE, pca = TRUE, partial_pca = FALSE, max_iter = 1000, verbose = getOption("verbose", FALSE), is_distance = FALSE, Y_init = NULL, pca_center = TRUE, pca_scale = FALSE, normalize = F, momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 10)
  }
  save(tsneRst, file = "data/digitTsneRst-2-3.Rdata")
  
 
   
  
  tsne2 = tsneRst[[1]]$Y
  tmp = aggregate(list(ind = 1:length(Y)), by = list(Y), function(x) x)
  palette = randomcoloR::distinctColorPalette(length(unique(Y)))
  tmpcols = unlist(mapply(function(x, y) rep(x, length(y)), palette, tmp[[2]]))
  tmpcols = tmpcols[order(unlist(tmp[[2]]))]
  pdf("figure/digitTsne2d.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  # plot3D::points3D(tsne3[,1], tsne3[,2], tsne3[,3])
  plot(tsne2, pch = 16, col = tmpcols, cex = 0.5, bty = "L", xlab = "x", ylab = "y", cex.lab = 2, cex.axis = 1.5, cex.lab = 2, cex.axis = 1.5)
  legend("bottomleft", bty = "n", col = unique(tmpcols)[order(unique(Y))], legend = sort(unique(Y)), pch = rep(15, 10), cex = 1.5)
  dev.off()
  
  
  
  
  # Try K-means one the best lowD embedding.
  if(T)
  {
    
    
    Xnew = tsneRst[[1]]$Y
    tXnew = t(Xnew)
    set.seed(123)
    rst = list()
    k = 10
    for(k in 3:20)
    {
      iend = 30L
      if(k == 10L) iend = 100L
      rst[[length(rst) + 1]] = list()
      for(i in 1:iend)
      {
        cat(i, "")
        rst[[length(rst)]][[i]] = GMKMcharlie::KM(X = tXnew, centroid = tXnew[, sample(1:nrow(Xnew), k)], verbose = F, maxCore = 7)
      }
    }
    
    
    rst = lapply(rst, function(X) 
    {
      X[[which.min(unlist(lapply(X, function(x) sum(unlist(lapply(x, function(y) sum(y$member2centroidDistance ^ 2)))))))]]
    })
    
    
    inclusterSS = unlist(lapply(rst, function(x) sum(unlist(lapply(x, function(y) sum((y$member2centroidDistance) ^ 2))))))
    
    
    pdf("figure/kmeansDigit10clustersUseTsne2dSS.pdf", width = 8, height = 8 * 0.618)
    par(mar = c(4.1, 5, 0, 0), family = "serif")
    plot(y = inclusterSS, x = 3:20, xlab = "Number of clusters", ylab = "In-cluster sum of square", bty = "L", cex.lab = 2, cex.axis = 1.5, lwd = 2, col = "darkblue", type = "l")
    inclusterSSdiff = diff(inclusterSS)
    inclusterSSdiff = c(inclusterSSdiff[1], inclusterSSdiff)
    inclusterSSdiff = (inclusterSSdiff - min(inclusterSSdiff)) / diff(range(inclusterSSdiff)) * diff(range(inclusterSS)) + min(inclusterSS)
    lines(y = inclusterSSdiff, x = 3:20, lwd = 2, col = "red", type = "l")
    lines(x = 13, y = inclusterSSdiff[11], type = "p", cex = 2)
    legend("right", legend = c("In-cluster sum of square", "Difference"), bty = "n", cex = 2, pch = c(15, 15), col = c("darkblue", "red"))
    dev.off()
    
    
    
    
    rst12cl = rst[[11]]
    bestClusteringYs = lapply(rst12cl, function(x) Y[x$clusterMember])
    bestClusteringYs = bestClusteringYs[order(unlist(lapply(bestClusteringYs, function(tmp)
    {
      sort(unique(tmp))[which.max(table(tmp))]
    })))]
    pdf("figure/kmeansDigit13clustersUseTsne2d.pdf", width = 10, height = 10 * 0.3)
    par(mar = c(3, 4, 0, 0), family = "serif", mfrow = c(3, 5))
    # tmp = c(7, 3, 6, 2, 5, 4, 0, 1, 9, 9, 8, 1)
    for(i in 1:length(bestClusteringYs))
    {
      tmp = bestClusteringYs[[i]]
      hist(tmp, breaks = seq(-0.5, by = 1, len = 11), col = scales::alpha("darkblue", 0.5), xlab = "", ylab = "", xaxt = "n", border = NA, main = "")
      axis(side = 1, at = 0:9, labels = 0:9, cex.axis = 1.25)
      legend("center", legend = sort(unique(tmp))[which.max(table(tmp))], cex = 2, bty = "n", text.col = "red")
    }
    dev.off()
    
    
  }
  
  
  
  
  # Try run NN.
  if(T)
  {
    
    
    tsneRst = list()
    for(d in 2:3)
    {
      cat(d, "")
      tmp = rbind(as.matrix(X), as.matrix(testX))
      tsneRst[[length(tsneRst) + 1]] = Rtsne::Rtsne(tmp, dims = d, initial_dims = 50, perplexity = 30, theta = 0.5, check_duplicates = TRUE, pca = TRUE, partial_pca = FALSE, max_iter = 1000, verbose = getOption("verbose", FALSE), is_distance = FALSE, Y_init = NULL, pca_center = TRUE, pca_scale = FALSE, normalize = F, momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 10)
    }
    
  
    tmp = tsneRst[[1]]
    Xnew = tmp$Y[1:nrow(X), ]
    testXnew = tmp$Y[(nrow(X) + 1L):nrow(tmp$Y), ]
    
    
    model = keras::keras_model_sequential()
    keras::layer_flatten(model, input_shape = ncol(Xnew))
    
    
    keras::layer_dense(model, units = 128, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, units = 64, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, 10, activation = "softmax")
    
    
    opt = keras::optimizer_adam(lr = 0.0005)
    keras::compile(model, loss = "sparse_categorical_crossentropy", optimizer = opt, metrics = "accuracy")
    system.time({fm = keras::fit(model, x = Xnew, y = Y, batch_size = 32, epochs = 100, validation_split = 0.3, verbose = 2, shuffle = T, validation_data = list(x_val = testXnew, y_val = testY))}) # https://keras.rstudio.com/reference/fit.html
    # fm
    # lapply(fm$metrics$loss)
    # Plot
    if(T)
    {
      
      
      lossRange = range(c(range(fm$metrics$loss), range(fm$metrics$val_loss)))
      # accRange = range(c(range(fm$metrics$accuracy), range(fm$metrics$val_accuracy)))
      accRange = c(0, 1)
      lossYaxis = seq(lossRange[1], lossRange[2], len = 5)
      accYaxis = seq(accRange[1], accRange[2], len = 5)
      tmp = fm$metrics
      tmp$loss = (tmp$loss - lossRange[1]) / diff(lossRange)
      tmp$val_loss = (tmp$val_loss - lossRange[1]) / diff(lossRange)
      tmp$accuracy = (tmp$accuracy - accRange[1]) / diff(accRange)
      tmp$val_accuracy = (tmp$val_accuracy - accRange[1]) / diff(accRange)
      ylim = range(unlist(tmp))
      pdf("figure/digitTsne2compNeuralnet.pdf", width = 8, height = 8 * 0.618)
      par(mar = c(4.1, 5, 1, 5), family = "serif")
      plot(tmp$loss, type = "l", lwd = 2, col = "darkblue", yaxt = "n", xlab = "Epoch", ylab = "Loss", cex.lab = 2, cex.axis = 1.5, bty = "n", ylim = ylim)
      axis(side = 2, at = seq(0, 1, len = 5), labels = round(lossYaxis, 2), cex.axis = 1.5)
      axis(side = 4, at = seq(0, 1, len = 5), labels = round(accYaxis, 2), cex.axis = 1.5)
      mtext("Accuracy", side = 4, line = 3, cex = 2)
      lines(tmp$val_loss, type = "l", lwd = 2, col = "red")
      lines(tmp$accuracy, type = "l", lwd = 2, col = "skyblue", lty = 1)
      lines(tmp$val_accuracy, type = "l", lwd = 2, col = "olivedrab3", lty = 1)
      legend("right", legend = c("Training loss", "Validation loss", "Training accuracy", "Validation accuracy"), lwd = c(2, 2, 2, 2), col = c("darkblue", "red", "skyblue", "olivedrab3"), bty = "n", cex = 2)
      dev.off()
    }  
    
  }
  
  
  
  
  # Run NN on clustered lowD
  if(T)
  {
    
    
    tsneRst = list()
    for(d in 2:3)
    {
      cat(d, "")
      tmp = rbind(as.matrix(X), as.matrix(testX))
      tsneRst[[length(tsneRst) + 1]] = Rtsne::Rtsne(tmp, dims = d, initial_dims = 50, perplexity = 30, theta = 0.5, check_duplicates = TRUE, pca = TRUE, partial_pca = FALSE, max_iter = 1000, verbose = getOption("verbose", FALSE), is_distance = FALSE, Y_init = NULL, pca_center = TRUE, pca_scale = FALSE, normalize = F, momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 10)
    }
    
    
    
    
    tmp = tsneRst[[1]]
    Xnew = tmp$Y[1:nrow(X), ]
    testXnew = tmp$Y[(nrow(X) + 1L):nrow(tmp$Y), ]
    
    
    dir.create("gauResult", showWarnings = F)
    trainDataPath = "data/digitXtsneTrain.csv"
    testDataPath = "data/digitXtsneTest.csv"
    trainSavePath = "data/digitXtsneTrainGauScore.csv"
    testSavePath = "data/digitXtsneTestGauScore.csv"
    data.table::fwrite(Xnew, file = trainDataPath, col.names = F)
    data.table::fwrite(testXnew, file = testDataPath, col.names = F)
    system(paste0("python python/gmm.py ", trainDataPath, " ", testDataPath, " ", Nclust, " ", trainSavePath, " ", testSavePath), wait = T)
    
    
    Xnew = as.matrix(data.table::fread(trainSavePath, header = F))
    Xnew = (Xnew - min(Xnew)) / max(Xnew)
    testXnew = as.matrix(data.table::fread(testSavePath, header = F))
    testXnew = (testXnew - min(Xnew)) / max(Xnew)
    
    
    model = keras::keras_model_sequential()
    keras::layer_flatten(model, input_shape = ncol(Xnew))
    
    
    keras::layer_dense(model, units = 128, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, units = 64, activation = NULL) # "relu"
    keras::layer_activation_leaky_relu(model, alpha = 0.01)
    keras::layer_dropout(model, rate = 0.5)
    
    
    keras::layer_dense(model, 10, activation = "softmax")
    
    
    opt = keras::optimizer_adam(lr = 0.0005)
    keras::compile(model, loss = "sparse_categorical_crossentropy", optimizer = opt, metrics = "accuracy")
    system.time({fm = keras::fit(model, x = Xnew, y = Y, batch_size = 32, epochs = 100, validation_split = 0.3, verbose = 2, shuffle = T, validation_data = list(x_val = testXnew, y_val = testY))}) # https://keras.rstudio.com/reference/fit.html
    # fm
    # lapply(fm$metrics$loss)
    # Plot
    if(T)
    {
      
      
      lossRange = range(c(range(fm$metrics$loss), range(fm$metrics$val_loss)))
      # accRange = range(c(range(fm$metrics$accuracy), range(fm$metrics$val_accuracy)))
      accRange = c(0, 1)
      lossYaxis = seq(lossRange[1], lossRange[2], len = 5)
      accYaxis = seq(accRange[1], accRange[2], len = 5)
      tmp = fm$metrics
      tmp$loss = (tmp$loss - lossRange[1]) / diff(lossRange)
      tmp$val_loss = (tmp$val_loss - lossRange[1]) / diff(lossRange)
      tmp$accuracy = (tmp$accuracy - accRange[1]) / diff(accRange)
      tmp$val_accuracy = (tmp$val_accuracy - accRange[1]) / diff(accRange)
      ylim = range(unlist(tmp))
      pdf("figure/digitTsnecompClusterNeuralnet.pdf", width = 8, height = 8 * 0.618)
      par(mar = c(4.1, 5, 1, 5), family = "serif")
      plot(tmp$loss, type = "l", lwd = 2, col = "darkblue", yaxt = "n", xlab = "Epoch", ylab = "Loss", cex.lab = 2, cex.axis = 1.5, bty = "n", ylim = ylim)
      axis(side = 2, at = seq(0, 1, len = 5), labels = round(lossYaxis, 2), cex.axis = 1.5)
      axis(side = 4, at = seq(0, 1, len = 5), labels = round(accYaxis, 2), cex.axis = 1.5)
      mtext("Accuracy", side = 4, line = 3, cex = 2)
      lines(tmp$val_loss, type = "l", lwd = 2, col = "red")
      lines(tmp$accuracy, type = "l", lwd = 2, col = "skyblue", lty = 1)
      lines(tmp$val_accuracy, type = "l", lwd = 2, col = "olivedrab3", lty = 1)
      legend("right", legend = c("Training loss", "Validation loss", "Training accuracy", "Validation accuracy"), lwd = c(2, 2, 2, 2), col = c("darkblue", "red", "skyblue", "olivedrab3"), bty = "n", cex = 2)
      dev.off()
      
      
    }
    
    
    
    
    
  }
  
  
  
}










