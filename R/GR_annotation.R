#' @export

get_GR_annotation <- function(vcf_files) {
  vcf_df = read.table(vcf_files, sep = '\t', stringsAsFactors = FALSE)
  vcf_df_bed = filter_cirBed (vcf_df)
  vcf_df_bed_gr = getGRanges(vcf_df_bed)

  vcf_df_bed_ext = getGRanges_single(vcf_df_bed)
  vcf_df_bed_ext_gr = getGRanges(vcf_df_bed_ext)

  txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
  vcf_ext_gr_annot <- locateVariants( vcf_df_bed_ext_gr, txdb_hg38, AllVariants())
  org_hum <- org.Hs.eg.db
  mcols(vcf_ext_gr_annot)[["SYMBOL"]] <- select(org_hum, keys = mcols(vcf_ext_gr_annot)$GENEID, columns = "SYMBOL", keytype = "ENTREZID" )[,2]
  if ( sum(is.na(vcf_ext_gr_annot$SYMBOL)) != 0 ) {
    vcf_ext_gr_annot[ which( is.na(vcf_ext_gr_annot$SYMBOL) )]$SYMBOL = "UNKNOWN"
  }
  mcols(vcf_ext_gr_annot)[["loc_info"]] = str_c(as.character(seqnames(vcf_ext_gr_annot)),":",
                                                start(ranges(vcf_ext_gr_annot)),":",
                                                str_replace_na(mcols(vcf_ext_gr_annot)$LOCATION),":",
                                                #        str_replace_na(mcols(vcf_ext_gr_annot)$LOCSTART),":",
                                                #        str_replace_na(mcols(vcf_ext_gr_annot)$TXID),":",
                                                #        str_replace_na(mcols(vcf_ext_gr_annot)$GENEID),":",
                                                str_replace_na(mcols(vcf_ext_gr_annot)$SYMBOL) )
  match_df <- as.data.frame(findOverlaps(vcf_df_bed_gr, vcf_ext_gr_annot))

  match_df_split <- split(match_df, match_df$queryHits )

  info_list = vector(mode = "list", length = length(vcf_df_bed_gr))

  for (i in 1:length(match_df_split)) {
    info_list[[as.integer(names(match_df_split)[i])]] = mcols(vcf_ext_gr_annot[match_df_split[[i]]$subjectHits])[["loc_info"]]
  }
  for (i in 1:length(info_list) ) {
    if (is.null(info_list[[i]]) ) {
      info_list[[i]] = "Intergenic:no_txdbCall"
    }
  }
  scaled_region = info_list
  for (i in 1:length(scaled_region) ) {
    if (scaled_region[[i]] != "Intergenic:no_txdbCall") {
      print(i)
      x = scaled_region[[i]]
      x_split = str_split(x, ":", simplify = TRUE) %>% as.data.frame(stringsAsFactors = FALSE)
      x_list = split(x_split, x_split$V4)
      result = lapply(x_list, sep_SYMBOL) %>% unlist() %>% unname()
      scaled_region[[i]] = result
    } else {
      chr = seqnames(vcf_df_bed_gr[i]) %>% as.character()
      start = ranges(vcf_df_bed_gr[i]) %>% start()
      end   = ranges(vcf_df_bed_gr[i]) %>% end()
      tag = str_c(chr,":",start,":",end,":","Intergenic",":","no_txdbCall")
      scaled_region[[i]] = tag
    }
  }
  scaled_region = CharacterList(scaled_region)
  mcols(vcf_df_bed_gr)[["Annotation"]] = scaled_region
  return (vcf_df_bed_gr)
}

