chromSampler <- function(dat = NULL,
                         info.column = NULL,
                         data.column = NULL,
                         is.diploid = T){
  # this function will sample a single chromosome number from a given
  # range of chromosome numbers.
  # chromosome numbers can only be separated by "," or "-"
  # assumes chromosome number increases by 2
  for(i in 1:nrow(dat)){
    if(gregexpr(pattern = "-",dat[i,info.column])[[1]][1] > 0 | gregexpr(pattern = ",",dat[i,info.column])[[1]][1] > 0){
      if(dat[,info.column][i] != ""){
        foo.a <- as.list(unlist(strsplit(dat[,info.column][i], split = ",")))
        foo.x <- vector(mode = "list", length = length(foo.a))
        for(j in 1:length(foo.a)){
          if((0 < gregexpr(pattern = "-", foo.a[[j]])[[1]])){
            x <- unlist(strsplit(x = foo.a[[j]], split = "-"))
            if(is.diploid == T){
              foo.x[[j]] <- seq(from = as.numeric(x[1]),
                                to = as.numeric(x[2]),
                                by = 2)
            }else{
              foo.x[[j]] <- seq(from = as.numeric(x[1]),
                                to = as.numeric(x[2]),
                                by = 1)
            }
          }else{
            foo.x[[j]] <- as.numeric(foo.a[[j]])
          }
          dat[,data.column][i] <-  sample(unlist(foo.x), size = 1)
        }
      }
    }
  }
  return(dat)
}