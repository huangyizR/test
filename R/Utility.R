#' @export

split_by_sequence = function (gene_df) {
  oo = gene_df$V2 %>% as.integer()
  k = diff(oo)
  for (i in 1:length(k) ) {
    if (k[i] == 1){
    } else {
      temp = k[i]
      i = i + 1
      while(  i <= length(k) ){
        if ( k[i] == 1) {
          k[i] = temp
          i = i + 1
        } else {
          break
        }
        if (i > length(k)) {
          break
        }
      }
    }
  }
  k = c(1,k)
  k = split(gene_df,inverse.rle(within.list(rle(k), values <- seq_along(values))))
  return (k)
}

#' @export
get_region = function(x) {
  min = x$V2[1] %>% as.integer()
  max = x$V2[nrow(x)] %>% as.integer()
  result = str_c(x$V1[1],":",min,":",max,":",x$V3[1],":",x$V4[1])
  return (result)
}
#' @export
get_gene_region = function (x){
  x = x[which(!duplicated(x)),]
  if (nrow(x) != 1) {
    x = split_by_sequence(x)
    x = lapply(x,get_region)
  } else {
    y = vector("list",1)
    y[[1]] = get_region(x)
    x = y
  }
  return (x)
}

#' @export
sep_SYMBOL = function (x) {
  x = split(x,x$V3)
  x = lapply(x,get_gene_region)
  return (x)
}

