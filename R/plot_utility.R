#' @export
get_annotation = function (vector_string) {
  a_s = paste(vector_string , collapse = '') %>% str_detect("BHLHE39")
  return (a_s)
}

#' @export
trim_intron = function(a_gene) {
  a_gene_intron = a_gene[a_gene$region == "intron"]
  if (length(a_gene_intron) == 0) {
    return(a_gene)
  }
  a_gene_non_intron = a_gene[a_gene$region != "intron"]
  match = findOverlaps(a_gene_intron,a_gene_non_intron)
  if (length(match) == 0) {
    return (a_gene)
  }
  match_df <- as.data.frame(match)
  match_df_split <- split(match_df, match_df$queryHits )
  intron_list = vector(mode = "list", length = length(match_df_split))
  trimed_intron = granges(a_gene_intron)[0]
  for (i in 1:length(intron_list)){
    a_intron = match_df_split[[i]][1,1]
    non_intron = match_df_split[[i]]$subjectHits

    a_intron = a_gene_intron[a_intron]
    non_intron = a_gene_non_intron[non_intron]

    trimed_intron = c(trimed_intron,setdiff(a_intron,non_intron))
  }
  if (trimed_intron %>% length() != 0) {
    trimed_intron$region = "intron"
    trimed_intron$SYMBOL = a_gene_intron$SYMBOL[1]
    result = c(trimed_intron, a_gene_non_intron)
  } else {
    result = c(a_gene_non_intron)
  }
  return (result)
}

#' @export
grToCirTrack = function (a_gr_Annotation, trimIntron = TRUE) {
  a_gr_Annotation = str_split(a_gr_Annotation,":", simplify = TRUE)
  a_gr_Annotation = as.data.frame(a_gr_Annotation, stringsAsFactors = FALSE)
  names(a_gr_Annotation) = c("gene","start","end","region","SYMBOL")
  a_gr_Annotation$start = as.integer(a_gr_Annotation$start)
  a_gr_Annotation$end = as.integer(a_gr_Annotation$end)
  a_gr_Annotation_gene_GR = GRanges (seqnames = Rle(a_gr_Annotation$gene),
                                     ranges = IRanges(a_gr_Annotation$start,
                                                      end = a_gr_Annotation$end),
                                     region = a_gr_Annotation$region,
                                     SYMBOL = a_gr_Annotation$SYMBOL)
  a_gr_Annotation_gene_GR_split = split(a_gr_Annotation_gene_GR,a_gr_Annotation_gene_GR$SYMBOL)
  #if ( sum(a_gr_Annotation_gene_GR_split$region == "intron" ) >= 2) {
  #    trimIntron == TRUE;
  #}
  if (trimIntron){
    for (i in 1:length(a_gr_Annotation_gene_GR_split)) {
      a_gr_Annotation_gene_GR_split[[i]] = trim_intron(a_gr_Annotation_gene_GR_split[[i]])
    }
  }
  a_gr_Annotation_gene_GR_split_trim = unlist(a_gr_Annotation_gene_GR_split)
  cirTrack_input = data.frame(gene = seqnames(a_gr_Annotation_gene_GR_split_trim) %>% as.character(),
                              start = start(ranges(a_gr_Annotation_gene_GR_split_trim)),
                              end =  end(ranges(a_gr_Annotation_gene_GR_split_trim)) ,
                              region = a_gr_Annotation_gene_GR_split_trim$region,
                              SYMBOL = a_gr_Annotation_gene_GR_split_trim$SYMBOL,
                              stringsAsFactors = FALSE)
  return (cirTrack_input)
}

#' @export
purning_sort_gr = function(granges) {

  sorted_seqnames = str_c("chr",1:22)
  sorted_seqnames = c(sorted_seqnames, "chrX")
  sorted_seqnames = c(sorted_seqnames, "chrY")
  seqlevels(granges, pruning.mode = "coarse") <- sorted_seqnames
  seqlevels(granges) = sorted_seqnames
  granges = sort(granges)
  return (granges)
}

#' @export
getOverLap_gr = function (template_gr, gr_list) {
  n = length(gr_list)
  temp_m = matrix(0,length(template_gr),n)
  colnames(temp_m) = names(gr_list)
  hit_df = as.data.frame(temp_m, stringsAsFactors = FALSE)
  for (i in 1:length(gr_list)) {
    hit = findOverlaps(gr_list[[i]], template_gr);
    hit_n = subjectHits(hit)[!duplicated(subjectHits(hit))]
    hit_df[,i][hit_n] = 1
  }
  values(template_gr) <- DataFrame(hit_df)
  return(template_gr)
}


#' @export
GRsubsetByRegionLen = function (gr, chr ,start, end, strand = "*" , length = 0) {
  subset_gr = gr[seqnames(gr) == chr & start(gr) > start & end(gr) < end & strand(gr) ==  strand]
  if (length != 0) {
    gr[width(gr) < length]
  }
  return (subset_gr)
}

#' @export
gr2bed = function(gr) {
  df = data.frame(Chromosome  = seqnames(gr) %>% as.character(),
                  start_coord = start(ranges(gr)),
                  end_coord = end(ranges(gr)),
                  dis_cord = values(gr)$dis_cord,
                  split_reads = values(gr)$split_reads,
                  circle_score = values(gr)$circle_score,stringsAsFactors = FALSE)
  return (df)
}

#' @export
parse_mantaID = function (vcf) {
  ID_parse = vcf$MATEID %>% str_split(":",simplify = TRUE)
  ID = str_c(ID_parse[,1],":",ID_parse[,2],":",ID_parse[,3],":", ID_parse[,4], ":",ID_parse[,5])
  return (ID)
}

#' @export
vrtobed = function (vcf_vr) {
  gr = granges(vcf_vr)
  df = data.frame(chr = as.character(seqnames(vcf_vr)),
                  start = ranges(granges(vcf_vr)) %>% start(), end = ranges(granges(vcf_vr)) %>% end(),stringsAsFactors = FALSE)
  return(df)
}

#' @export
gr2bed = function(gr) {
  df = data.frame(Chromosome  = seqnames(gr) %>% as.character(),
                  start_coord = start(ranges(gr)),
                  end_coord = end(ranges(gr)),
                  mathylation = values(gr)$mathylation,stringsAsFactors = FALSE)
  return (df)
}
