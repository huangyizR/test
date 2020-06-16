#' @export
multi_plot_ecDNA = function (initialize_template, ecDNA_list, save = TRUE) {
  gene_symbol = c()
  for (i in 1:length(ecDNA_list)){
    gene_symbol = c(gene_symbol, ecDNA_list[[i]]$SYMBOL %>% unique())
  }
  gene_symbol = gene_symbol %>% unique()
  gene_text = "KeyGenes";
  for (i in 1:length(gene_symbol) ) {
    gene_text = str_c(gene_text, "_",gene_symbol[i]);
  }
  if (save) {
    title = str_c("POETIC36_cfVSt2_",gene_text,".pdf")
    pdf(title)
  }

  circos.par("track.height" = 0.05)
  circos.genomicInitialize(initialize_template)
  n = 1
  color =  c("#FF000040", "#00FF0040", "#0000FF40")
  for (i in 1:length(ecDNA_list) ) {


    circos.genomicTrack(ecDNA_list[[i]], ylim = c(0, n + 0.5 ),
                        panel.fun = function(region, value, ...) {
                          genes = unique(value$SYMBOL)
                          for (i in seq_along(genes)){
                            a_gene_value = value[value$SYMBOL == genes[i],]
                            print(head(value[value$SYMBOL == genes[i],]))
                            a_gene_region = region[value$SYMBOL == genes[i],]
                            current_gene_start = min(a_gene_region$start)
                            current_gene_end   = max(a_gene_region$end)
                            if (genes[i] == "UNKNOWN"){
                              circos.lines(c(current_gene_start, current_gene_end),
                                           c(n - i/4 + 1/4, n - i/4 + 1/4), col = color[i])
                              circos.genomicRect(a_gene_region[1, , drop = FALSE], ytop = n - i/4 + 1/4 + 0.29,
                                                 ybottom = n - i/4 + 1/4 - 0.29, col = "#A9A9A9", border = NA)
                              next
                            }
                            circos.lines(c(current_gene_start, current_gene_end),
                                         c(n - i/4 + 1/4, n - i/4 + 1/4), col = color[i])
                            test_y = n - i/4 + 1/4 + 0.4
                            test_x = (current_gene_start + current_gene_end) / 2
                            circos.text(test_x, test_y, a_gene_value$SYMBOL[1], cex = 0.6)

                            for (j in 1:nrow(a_gene_value)) {
                              if (a_gene_value$region[j] == "coding") {
                                print("found a coding!")
                                circos.genomicRect(a_gene_region[j, , drop = FALSE], ytop = n - i/4 + 1/4 + 0.60,
                                                   ybottom = n - i/4 + 1/4 - 0.60, col = "#BD381B", border = NA)    #red
                              } else if (a_gene_value$region[j] == "intron") {
                                print("found a intron!")
                                circos.genomicRect(a_gene_region[j, , drop = FALSE], ytop = n - i/4 + 1/4 + 0.39,
                                                   ybottom = n - i/4 + 1/4 - 0.39, col = "#00FF0040", border = NA)  #light green
                              } else if (a_gene_value$region[j] == "promoter") {
                                circos.genomicRect(a_gene_region[j, , drop = FALSE], ytop = n - i/4 + 1/4 + 0.53,
                                                   ybottom = n - i/4 + 1/4 - 0.53, col = "#0000FF40", border = NA)   #dark blue
                              } else if (a_gene_value$region[j] == "spliceSite") {
                                circos.genomicRect(a_gene_region[j, , drop = FALSE], ytop = n - i/4 + 1/4 + 0.65,
                                                   ybottom = n - i/4 + 1/4 - 0.65, col = "#C018C8", border = NA)    #pink
                              } else if (a_gene_value$region[j] == "threeUTR") {
                                circos.genomicRect(a_gene_region[j, , drop = FALSE], ytop = n - i/4 + 1/4 + 0.29,
                                                   ybottom = n - i/4 + 1/4 - 0.29, col = "#FF000040", border = NA)   # yellow
                              }
                            }
                          }
                        }, bg.border = NA)

  }


  text(0,1, gene_text, col = 1)
  legend("topleft",
         c("coding","intron","promoter","spliceSite","threeUTR","intergenic"),
         fill=c("#BD381B","#00FF0040","#0000FF40","#C018C8","#FF000040","#A9A9A9") )
  circos.clear()
  if (save) {
    dev.off()
  }
}
