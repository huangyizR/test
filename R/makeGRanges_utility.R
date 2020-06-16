#' @export

get_gene_info = function(x) {
  x = mcols(x)
  info = str_c(x$LOCATION, ":", x$LOCSTART, ":", x$TXID, ":" , x$GENEID , ":", x$SYMBOL)
  return(info)
}

#' @export

getGRanges <- function (circle_df) {
  seqnames = circle_df$Chromosome
  ranges = IRanges(circle_df$start_coord, end = circle_df$end_coord)
  dis_cord = circle_df$dis_cord
  split_reads = circle_df$split_reads
  circle_score = circle_df$circle_score

  gr <- GRanges(
    seqnames = Rle(seqnames),
    ranges = ranges,
    dis_cord = dis_cord,
    split_reads = split_reads,
    circle_score = circle_score,
    # loc_ID = circle_df$loc_ID
  )

  return (gr)
}


#' @export

temp_f = function (row_of_circle_df) {
  start_coord <- row_of_circle_df[1,"start_coord"]:row_of_circle_df[1,"end_coord"]
  end_coord   <- start_coord
  df <- data.frame ( Chromosome  = row_of_circle_df$Chromosome ,
                     start_coord = start_coord,
                     end_coord   = end_coord,
                     dis_cord    = row_of_circle_df$dis_cord,
                     split_reads = row_of_circle_df$split_reads,
                     circle_score = row_of_circle_df$circle_score,
                     loc_ID   = row_of_circle_df$loc_ID,
                     stringsAsFactors = FALSE)
  return(df)
}


#' @export

getGRanges_single <- function (circle_df) {
  n_rows <- nrow(circle_df)
  loc_ID  = str_c(circle_df$Chromosome,":",circle_df$start_coord,":",circle_df$end_coord)
  circle_df = mutate(circle_df,"loc_ID" = loc_ID)
  circle_df = mutate(circle_df, "number" = 1:n_rows)
  extended_df_list <-  lapply(split(circle_df, circle_df$number), temp_f)
  extended_circle_df <- do.call(rbind, extended_df_list)
  return (extended_circle_df)
}

