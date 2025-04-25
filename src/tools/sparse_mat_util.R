#function for analysis 

sweep_sparse = function(x, margin, stats, fun = "*") {
  #handle error 
  if (!class(x)[1] %in% c("dgTMatrix","dgCMatrix","dgRMatrix","dgeMatrix","matrix")){
    stop("Only several saprse matrix types are supported: dgTMatrix,dgCMatrix,dgRMatrix,dgeMatrix,matrix")
  }
  if((length(stats)!=1) & (length(stats)!=dim(x)[margin])){
    stop("stats doesn't match the dimension of the certain margin")
  }
  #deal with stats
  if (class(stats)[1] == "matrix"){
    stats = c(stats)
  }
  #to deal with a non-sparse matrix
  if (class(x)[1]  == "matrix"){
    return(sweep(x, margin,stats,fun))
  }
  if (class(x)[1] == "dgeMatrix"){
    res <- sweep(as.matrix(x), margin, stats, fun)
    res <- as(res, "dgeMatrix")
    return(res)
  }
  if (!class(x)[1] %in% c("dgTMatrix","dgCMatrix","dgRMatrix")){
    stop("Only several saprse matrix types are supported: dgTMatrix, dgCMatrix, dgRMatrix")
  }
  if (class(x)[1]=="dgRMatrix"){
    i <- rep(1:x@Dim[1], diff(x@p)) -1
  }else{
    i <- x@i
  }
  if (class(x)[1]=="dgCMatrix"){
    j <- rep(1:x@Dim[2], diff(x@p)) -1
  }else{
    j <- x@j
  }
  if (margin == 1) {
    idx <- i + 1
  } else {
    idx <- j + 1
  }
  if (class(stats)[1] == "matrix"){
    stats = c(stats)
  }
  if(length(stats )==1){
    stats <- rep(stats, x@Dim[margin])
  }
  f <- match.fun(fun)
  x@x <- f(x@x, stats[idx])
  return(x)
}

#sparse transformation: such as log+1
transform_sparse = function(x, fun = "log2") {
  if (!class(x)[1] %in% c("dgTMatrix","dgCMatrix","dgRMatrix")){
    stop("Only several saprse matrix types are supported: dgTMatrix, dgCMatrix, dgRMatrix")
  }
  if (class(fun) != "function"){
    f <- match.fun(fun)
  }else{
    f <- fun
  }
  x@x <- f(x@x)
  return(x)
}
