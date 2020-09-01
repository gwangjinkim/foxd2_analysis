

mkdir -p ~/Downloads/DEanalysis
rsync -rav -e "ssh -p 64168" josephus@132.230.165.168:/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907 ~/Downloads/DEanalysis
rsync -rav -e "ssh -p 64168" josephus@132.230.165.168:/media/josephus/Elements/DEanalysis/foxd2-C3-mk3-201907 ~/Downloads/DEanalysis

rsync -rav -e "ssh -p 64168" josephus@132.230.165.168:/media/josephus/Elements/DEanalysis/foxd2-foxd2-F7-201907 ~/Downloads/DEanalysis






conda create --name xlsxgrep
source activate xlsxgrep
conda install -c conda-forge python
pip install xlsxgrep

cd /media/josephus/Elements/DEanalysis
source activate xlsxgrep
xlsxgrep "kidney morphogenesis" --recursive --with-filename --with-sheetname .


for f in foxd2*;
  do xlsxgrep "kidney morphogenesis" --recursive --with-filename --with-sheetname ${f}
done

# hit!

rsync -rav -e "ssh -p 64168" josephus@132.230.165.168:/media/josephus/Elements/DEanalysis/foxd2-C3-F7-201907 ~/Downloads/DEanalysis










cmputer <- "work"

if (computer == "ss") {
  script_dir <- "/media/daten/arnold/josephus/script"
  out_dir    <- "/media/daten/arnold/josephus/results/DEanalaysis"
  meta_dir   <- "/media/daten/arnold/josephus/metas"
} else {
  script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
  out_dir    <- "/media/josephus/Elements/DEanalysis"
  meta_dir   <- "/media/josephus/Elements/metas"
}
dir.create(out_dir, recursive=TRUE)

source(file.path(script_dir, "190122.rdownstream.central.R"))






# foxd2-C3-F7-201907/c3-vs-f7/GO/DE_sig_foxd2-C3-F7_c3-vs-f7_190801160815_0.05_1_over_GO_BP.xlsx: down.l2FC.srt: 	GO:0001656	metanephros development	23/1882	87/23239	2.73636731704451e-07	4.15678999185964e-05	3.25926598261359e-05	Wnt4/Aph1a/Gdnf/Gata3/Rdh10/Lif/Sox8/Pdgfa/Kif26b/Eya1/Aqp1/Fras1/Pdgfb/Lhx1/Pou3f3/Hnf1b/Pax8/Pax2/Irx2/Wnt7b/Tfap2a/Sall1/Irx1	23.0


# foxd2-C3-F7-201907/c3-vs-f7/GO/DE_sig_foxd2-C3-F7_c3-vs-f7_190801160815_0.05_1_over_GO_BP.xlsx: down.l2FC.srt: 	GO:0016055	Wnt signaling pathway	65/1882	418/23239	2.65747153339829e-07	4.15678999185964e-05	3.25926598261359e-05	Gsk3a/Arntl/Mdfic/Sulf2/Gid8/Porcn/Notch1/Stk4/Aes/Wnt4/Trpm4/Rac1/Nfkb1/Wls/Nphp4/Src/Dkk3/Gata3/Xiap/Kdm6a/Atp6ap2/Ccnd1/Zfp703/Bicc1/Jade1/Wnt5a/Wnt9a/Prkaa2/Dvl1/Zbtb33/Vangl2/Fam53b/Egfr/Ccne1/Nrarp/Ror1/Gprc5b/Fzd4/Pitx2/Fzd6/Met/Daam2/Axin2/Lgr6/Prickle1/Sema5a/Wnt5b/Mgat3/Wnt10a/Scel/Nkx2-5/Folr1/Cdh1/Cav1/Sfrp1/Mir196a-1/Hnf1b/Ptpru/Nog/Celsr1/Wnt6/Ccdc88c/Wnt7b/Celsr2/Sall1	65.0


# read-in the GO table

# for a givn GO term slurp the genes

# read-in the sig-DE table

# select the genes

# build a heatmap - first without scaling - just absolute

# then with scaling and centering









###########################################
# which files I are dealing with?
###########################################

other_go_analysis_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-C3-F7-201907/c3-vs-f7/GO/DE_sig_foxd2-C3-F7_c3-vs-f7_190801160815_0.05_1_over_GO_BP.xlsx'

right_go_analysis_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907/f7-vs-wt/GO/DE_sig_foxd2-mk4_f7-vs-wt_190801155819_0.05_1_over_GO_BP.xlsx'

right_counts_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907/f7-vs-wt/count-table/nrm-counts-foxd2-mk4-f7-vs-wt.xlsx'

right_DE_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907/f7-vs-wt/DE-table/DE_sig_foxd2-mk4_f7-vs-wt_190801155819_0.05_1.xlsx'


###########################################
# which genes I am dealing with?
###########################################

kidney.morphogenesis <- c("Wnt4", "Tgfb1", "Gdnf", "Gata3", "Hs2st1", "Lif", "Vangl2", "Adamts16", "Hoxb7", "Ahr", "Sox8", "Lama5", "Wnk4", "Kif26b", "Fgf1", "Eya1", "Fras1", "Lhx1", "Npnt", "Hnf1b", "Nog", "Pax8", "Pax2", "Wnt6", "Irx2", "Wnt7b", "Sall1", "Irx1")


metanephros.development <- c("Wnt4", "Aph1a", "Gdnf", "Gata3", "Rdh10", "Lif", "Sox8", "Pdgfa", "Kif26b", "Eya1", "Aqp1", "Fras1", "Pdgfb", "Lhx1", "Pou3f3", "Hnf1b", "Pax8", "Pax2", "Irx2", "Wnt7b", "Tfap2a", "Sall1", "Irx1")


wnt.signalling.pathway <- c("Gsk3a", "Arntl", "Mdfic", "Sulf2", "Gid8", "Porcn", "Notch1", "Stk4", "Aes", "Wnt4", "Trpm4", "Rac1", "Nfkb1", "Wls", "Nphp4", "Src", "Dkk3", "Gata3", "Xiap", "Kdm6a", "Atp6ap2", "Ccnd1", "Zfp703", "Bicc1", "Jade1", "Wnt5a", "Wnt9a", "Prkaa2", "Dvl1", "Zbtb33", "Vangl2", "Fam53b", "Egfr", "Ccne1", "Nrarp", "Ror1", "Gprc5b", "Fzd4", "Pitx2", "Fzd6", "Met", "Daam2", "Axin2", "Lgr6", "Prickle1", "Sema5a", "Wnt5b", "Mgat3", "Wnt10a", "Scel", "Nkx2-5", "Folr1", "Cdh1", "Cav1", "Sfrp1", "Mir196a-1", "Hnf1b", "Ptpru", "Nog", "Celsr1", "Wnt6", "Ccdc88c", "Wnt7b", "Celsr2", "Sall1")



###########################################
# read-in each table
###########################################

require(xlsx2dfs)
other_go_dfs <- xlsx2dfs(other_go_analysis_fpath)
right_go_dfs <- xlsx2dfs(right_go_analysis_fpath)
right_counts_dfs <- xlsx2dfs(right_counts_fpath)
right_DE_dfs <- xlsx2dfs(right_DE_fpath)

other_go_df <- other_go_dfs[["all.l2FC.srt"]]
right_go_df <- right_go_dfs[["all.l2FC.srt"]]
right_counts_df <- right_counts_dfs$`nrm-counts`
right_DE_df <- right_DE_dfs$all.l2FC.srt

############################################
# get genes from GO
############################################

get_go_genes <- function(go_df, term) {
  strsplit(x=go_df[grepl(pattern=term, x=go_df$Description), "geneID"], 
           split="/")[[1]]
} # works!

get_go_id <- function(go_df, term) {
  rownames(go_df)[grepl(pattern=term, x=go_df$Description)]
}

show_go_terms <- function(go_df) {
  sort(unique(go_df$Description))
}

show_go_terms(other_go_dfs[[3]])
kidney.morphogenesis <- get_go_genes(other_go_dfs[[3]], "kidney morphogenesis")
metanephros.development <- get_go_genes(other_go_dfs[[3]], "metanephros development")
wnt.sig.pathway <- get_go_genes(other_go_dfs[[3]], "Wnt signaling pathway")

show_go_terms(right_go_dfs[[3]])

kidney.dev <- get_go_genes(right_go_dfs[[3]], "kidney development")


# how about generating for each GO term a heatmap? that is absolutely possible!


###############################################
# generate for the genes a heatmap
###############################################





###############################################
# helper functions
###############################################
time.now <- function() format(Sys.time(), "%y%m%d%H%M%S")

common.string <- function(s1, s2, acc="") {
  # Return common beginning string
  first.s1 = substr(s1, 1, 1)
  first.s2 = substr(s2, 1, 1)
  if (s1 == "" && s2 == "" || first.s1 != first.s2) {
    acc
  } else  { # s1[1] == s2[1]
    common.string(substr(s1, 2, nchar(s1)), substr(s2, 2, nchar(s2)), paste0(acc, first.s1))
  }
} # works!

cdr <- function(vec) if (length(vec) < 2) c() else vec[2:length(vec)]
car <- first <- function(vec) vec[1]
cadr <- second <- function(vec) vec[2]
cons <- function(x, vec) c(x, vec)
nreverse <- rev
null <- function(x) length(x) == 0

grouper <- function(vec, acc=c(), last.common=c(), prev="") {
  if (null(vec)) {
    nreverse(acc)
  } else if (length(vec) == 1) {
    grouper(cdr(vec), cons(last.common, acc))
  } else {
    next.common <- common.string(first(vec), second(vec))
    if (null(last.common)) {
      grouper(cdr(vec), cons(next.common, acc), next.common, car(vec))
    } else {
      prev.common <- common.string(prev, car(vec))
      if (last.common == prev.common) {
        grouper(cdr(vec), cons(prev.common, acc), prev.common, car(vec))
      } else {
        grouper(cdr(vec), cons(next.common, acc), next.common, car(vec))
      }
    }
  }
} # works!

split_by_group <- function(vec, indexes=FALSE) {
  dfs <- split(data.frame(x= if (indexes) 1:length(vec) else vec, 
                          stringsAsFactors=FALSE), 
               grouper(vec))
  lapply(dfs, function(df) df$x)
}

average_df <- function(df, col_names=NULL) {
  # averages df by similarity of columnnames
  if (is.null(col_names)) col_names <- colnames(df)
  column_indexes <- split_by_group(col_names, indexes=TRUE)
  res.df <- Reduce(cbind, lapply(column_indexes, function(idxs) rowMeans(df[, idxs]))) # rowSds
  colnames(res.df) <- names(column_indexes)
  res.df
} # works as expected!


###############################################
# heatmap function
###############################################

require(pheatmap)
generate_heatmap <- function(df, 
                             col_names=NULL,
                             scale=TRUE,
                             non_scale_log2=FALSE,
                             color.1 = "black",
                             color.2 = "red",
                             cluster_rows = F,
                             col_steps = 100,
                             fontsize_row = 5,
                             display_numbers = FALSE,
                             number_color = "white",
                             fontsize_number = 5,
                             outname = NULL) {
  # `outname` should be without ending!
  df.avg <- average_df(df, col_names = col_names)
  if (scale) {
    scaledata <- t(scale(t(df.avg)))
    scaledata <- scaledata[complete.cases(scaledata), ]
  } else {
    scaledata <- if (non_scale_log2) log2(df.avg + .1) else df.avg
  }
  # prepare heatmap
  colfunc <- colorRampPalette(c(color.1, color.2))
  
  if (is.null(outname)) {
    pheatmap(scaledata,
             cluster_rows = cluster_rows,
             cluster_cols = F,
             cellwidth = 40,
             col = colfunc(col_steps),
             fontsize_row = fontsize_row,
             border_color = NA,
             display_numbers = display_numbers,
             number_color = number_color,
             fontsize_number = fontsize_number)
  } else {
  # print heatmap
    {
      setEPS()
      postscript(paste0(outname, ".eps"))
      # svg(paste0(outname, ".svg"))
      pheatmap(scaledata,
               cluster_rows = cluster_rows,
               cluster_cols = F,
               cellwidth = 40,
               col = colfunc(col_steps),
               fontsize_row = fontsize_row,
               border_color = NA,
               display_numbers = display_numbers,
               number_color = number_color,
               fontsize_number = fontsize_number)
      dev.off()
    }

    {
      svg(paste0(outname, ".svg"))
      pheatmap(scaledata,
               cluster_rows = cluster_rows,
               cluster_cols = F,
               cellwidth = 40,
               col = colfunc(col_steps),
               fontsize_row = fontsize_row,
               border_color = NA,
               display_numbers = display_numbers,
               number_color = number_color,
               fontsize_number = fontsize_number)
      dev.off()
    }
  #   # print counts
  #   res.all.gois <- gois[gois %in% rownames(res$all.names.srt)]
  #   res.up.gois <- gois[gois %in% rownames(res$up.names.srt)]
  #   res.down.gois <- gois[gois %in% rownames(res$down.names.srt)]
  #   
  #   res.gois <- list(cnts.nrm[gois.present, ], cnts.gois, res.all.gois, res.up.gois, res.down.gois)
  #   names(res.gois) <- c("goi.counts", "goi.counts.avg", "goi.DE.all", "goi.DE.up", "goi.DE.down")
  #   write.dfs(res.gois, paste0(outname, ".xlsx"))
    dfs2xlsx(withNames("goi_all", scaledata),
            paste0(outname, ".xlsx"))
 }
}

df <- right_counts_df
generate_heatmap(df[kidney.morphogenesis, ], col_names=c(rep("F7", 5), rep("wt", 5)), scale=TRUE)

generate_heatmap(df[kidney.dev, ], col_names=c(rep("F7", 5), rep("wt", 5)), scale=TRUE)







generate_heatmap(df[kidney.morphogenesis, ], 
                 col_names=c(rep("F7", 5), rep("wt", 5)), 
                 scale=FALSE, 
                 fontsize_row = 5,
                 non_scale_log2=TRUE,
                 display_numbers = TRUE)

generate_heatmap(df[kidney.morphogenesis, ], 
                 col_names=NULL, 
                 scale=FALSE, 
                 fontsize_row = 5,
                 non_scale_log2=TRUE,
                 display_numbers = TRUE)





generate_heatmap(df[kidney.dev, ], 
                 col_names=NULL, 
                 scale=FALSE, 
                 fontsize_row = 5,
                 non_scale_log2=TRUE,
                 display_numbers = TRUE) 
                 
generate_heatmap(df[kidney.dev, ], 
                 col_names=NULL, 
                 scale=TRUE, 
                 fontsize_row = 5,
                 non_scale_log2=TRUE,
                 display_numbers = TRUE) # that looks good


generate_heatmap(df[kidney.dev, ], 
                 col_names=c(rep("F7", 5), rep("wt", 5)), 
                 scale=TRUE, 
                 fontsize_row = 5,
                 non_scale_log2=TRUE,
                 display_numbers = TRUE) 

generate_heatmap(df[kidney.dev, ], 
                 col_names=c(rep("F7", 5), rep("wt", 5)), 
                 scale=FALSE, 
                 fontsize_row = 5,
                 non_scale_log2=FALSE,
                 display_numbers = TRUE) # the most realistic picture

generate_heatmap(df[kidney.dev, ], 
                 col_names=c(rep("F7", 5), rep("wt", 5)), 
                 scale=FALSE, 
                 fontsize_row = 5,
                 non_scale_log2=TRUE,
                 display_numbers = TRUE) #

generate_heatmap(df[kidney.dev, ], 
                 col_names=c(rep("F7", 5), rep("wt", 5)), 
                 scale=TRUE, 
                 fontsize_row = 5,
                 non_scale_log2=TRUE,
                 display_numbers = TRUE) # most unrealistic black and white picture





generate_heatmap(df[wnt.sig.pathway, ], 
                 col_names=c(rep("F7", 5), rep("wt", 5)), 
                 scale=FALSE, 
                 fontsize_row = 3,
                 non_scale_log2=FALSE,
                 display_numbers = TRUE,
                 fontsize_number = 3) 
                 
# save all these variants for 
#  each genes/go terms
#  for each combination
#  in different folders
#  in svg and eps and jpg format
# and send them out!
# 
# if gene name length > 65 then go to 3






############################################
# test the function
############################################

out_dir = "/media/josephus/Elements/DEanalysis/tmp"
dir.create(out_dir, recursive = TRUE)

go_terms <- show_go_terms(other_go_dfs[[3]])

go_term <- go_terms[[1]]
go <- gsub(" ", "_", go_term)
genes <- get_go_genes(other_go_dfs[[3]], go_term)
size <- if (length(genes) <= 65) 5 else 3

mode <- "scaled_f7_vs_wt"
# 'DE_sig_foxd2-C3-F7_c3-vs-f7_190801160815_0.05_1_over_GO_BP.xlsx'
# '/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907/f7-vs-wt/count-table/nrm-counts-foxd2-mk4-f7-vs-wt.xlsx'

symbol2string <- function(x) deparse(substitute(x))
symbol2fstring <- function(x) gsub("\\.|\\$", "_", deparse(substitute(x)))

fname <- paste0(go, "_", mode, "_", time.now())
out_fpath <- file.path(out_dir, fname)
dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)


generate_heatmap(df[kidney.dev, ], 
                 col_names=c(rep("F7", 5), rep("wt", 5)), 
                 scale=TRUE, 
                 fontsize_row = size,
                 non_scale_log2=TRUE,
                 display_numbers = FALSE,
                 fontsize_number = size,
                 outname = out_fpath) #



plot_go_term_heatmap_with_numbers <- function(go_term, df, out_dir, go_df) {
  go <- gsub(" ", "_", go_term)
  genes <- get_go_genes(go_df, go_term)
  size <- if (length(genes) <= 65) 5 else 3
  
  mode <- "scaled_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=c(rep("F7", 5), rep("wt", 5)), 
                   scale=TRUE, 
                   fontsize_row = size,
                   non_scale_log2=TRUE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #

  mode <- "unscaled_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=c(rep("F7", 5), rep("wt", 5)), 
                   scale=FALSE, 
                   fontsize_row = size,
                   non_scale_log2=FALSE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #

  mode <- "unscaledlog2_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=c(rep("F7", 5), rep("wt", 5)), 
                   scale=FALSE, 
                   fontsize_row = size,
                   non_scale_log2=TRUE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #
  
}


df <- right_counts_df
go_df <- other_go_dfs[[3]]
out_dir = "/media/josephus/Elements/DEanalysis/heatmap_foxd2-F7-mk4-201907/f7-vs-wt/heatmaps_from_foxd2-C3-F7_c3-vs-f7_GO_BP"
plot_go_term_heatmap_with_numbers("kidney morphogenesis", df, out_dir, go_df)
plot_go_term_heatmap_with_numbers("metanephros development", df, out_dir, go_df)
plot_go_term_heatmap_with_numbers("Wnt signaling pathway", df, out_dir, go_df)

out_dir = "/media/josephus/Elements/DEanalysis/tmp"
unlink(out_dir, recursive=TRUE)



df <- right_counts_df
go_df <- right_go_dfs[[3]] # downregulated
out_dir = "/media/josephus/Elements/DEanalysis/heatmap_foxd2-F7-mk4-201907/f7-vs-wt/heatmaps_from_foxd2-mk4_f7-vs-wt_GO_BP_down"
go_terms <- show_go_terms(go_df)

for (go_term in go_terms) {
  plot_go_term_heatmap_with_numbers(go_term, df, out_dir, go_df)
}


# rsync -rav -e "ssh -p 64168" josephus@132.230.165.168:/media/josephus/Elements/DEanalysis/heatmap_foxd2-F7-mk4-201907 ~/Downloads/DEanalysis/


df <- right_counts_df
go_df <- right_go_dfs[[2]] # upregulated
out_dir = "/media/josephus/Elements/DEanalysis/heatmap_foxd2-F7-mk4-201907/f7-vs-wt/heatmaps_from_foxd2-mk4_f7-vs-wt_GO_BP_up"
go_terms <- show_go_terms(go_df)

for (go_term in go_terms) {
  plot_go_term_heatmap_with_numbers(go_term, df, out_dir, go_df)
}





################################################################################



###########################################
# retrive all genes associated to a GO term
###########################################
# https://www.biostars.org/p/52101/


other_go_analysis_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-C3-F7-201907/c3-vs-f7/GO/DE_sig_foxd2-C3-F7_c3-vs-f7_190801160815_0.05_1_over_GO_BP.xlsx'

require(xlsx2dfs)
other_go_dfs <- xlsx2dfs(other_go_analysis_fpath)

get_go_id <- function(go_df, term, exact = TRUE) {
  rownames(go_df)[grepl(pattern=if (exact) paste0("^", term, "$") else term, x=go_df$Description)]
}

get_go_term <- function(go_df, term) {
  go_df$Description[grepl(pattern=term, x=go_df$Description)]
}

go_terms <- c("kidney morphogenesis", "metanephros development", "Wnt signaling pathway")
go_ids <- sapply(go_terms, function(gt) get_go_id(other_go_dfs[[3]], gt))
go_ids_more <- sapply(go_terms, function(gt) get_go_id(other_go_dfs[[3]], gt, exact=FALSE))
go_terms_found <- sapply(go_terms, function(gt) get_go_term(other_go_dfs[[3]], gt))




require(biormaRt)

# go_id2genes_initialize <- function() {
#    ensembl <- NULL
#    df <- NULL
#    function(go_id) {
#      if (is.null(ensembl)) ensembl <<- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#      if (is.null(df)) df <<- getBM(attributes = c('external_gene_name', 'go_id'),
#                       mart = ensembl)
#      df$external_gene_name[df$go_id == go_id]
#    }
# }
# 
# go_id2genes <- go_id2genes_initialize()
# go_id2genes(go_ids[[1]])
# # Error in listMarts(host = host, path = path, port = port, includeHosts = TRUE,  : 
# #   Unexpected format to the list of available marts.
# # Please check the following URL manually, and try ?listMarts for advice.
# # http://www.ensembl.org:80/biomart/martservice?type=registry&requestid=biomaRt
# # 
# 
# 
# # ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# # listAttributes(ensembl)
# 

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
df <- getBM(attributes = c('external_gene_name', 'go_id'),
                      mart = ensembl)
genes <- list()
for (go_id in go_ids) {
  genes[[go_id]] <- df$external_gene_name[df$go_id == go_id]
}

right_go_analysis_fpath <- '/home/josephus/Downloads/DEanalysis/foxd2-foxd2-F7-20200403/f7-vs-wt_foxd2/GO/DE_sig_foxd2-foxd2-F7-corr_f7-vs-wt_foxd2_200403113317_0.05_1_over_GO_BP.xlsx' 






###############################################
# generate for the genes a heatmap
###############################################





###############################################
# helper functions
###############################################
time.now <- function() format(Sys.time(), "%y%m%d%H%M%S")

common.string <- function(s1, s2, acc="") {
  # Return common beginning string
  first.s1 = substr(s1, 1, 1)
  first.s2 = substr(s2, 1, 1)
  if (s1 == "" && s2 == "" || first.s1 != first.s2) {
    acc
  } else  { # s1[1] == s2[1]
    common.string(substr(s1, 2, nchar(s1)), substr(s2, 2, nchar(s2)), paste0(acc, first.s1))
  }
} # works!

cdr <- function(vec) if (length(vec) < 2) c() else vec[2:length(vec)]
car <- first <- function(vec) vec[1]
cadr <- second <- function(vec) vec[2]
cons <- function(x, vec) c(x, vec)
nreverse <- rev
null <- function(x) length(x) == 0

grouper <- function(vec, acc=c(), last.common=c(), prev="") {
  if (null(vec)) {
    nreverse(acc)
  } else if (length(vec) == 1) {
    grouper(cdr(vec), cons(last.common, acc))
  } else {
    next.common <- common.string(first(vec), second(vec))
    if (null(last.common)) {
      grouper(cdr(vec), cons(next.common, acc), next.common, car(vec))
    } else {
      prev.common <- common.string(prev, car(vec))
      if (last.common == prev.common) {
        grouper(cdr(vec), cons(prev.common, acc), prev.common, car(vec))
      } else {
        grouper(cdr(vec), cons(next.common, acc), next.common, car(vec))
      }
    }
  }
} # works!

split_by_group <- function(vec, indexes=FALSE) {
  dfs <- split(data.frame(x= if (indexes) 1:length(vec) else vec, 
                          stringsAsFactors=FALSE), 
               grouper(vec))
  lapply(dfs, function(df) df$x)
}

average_df <- function(df, col_names=NULL) {
  # averages df by similarity of columnnames
  if (is.null(col_names)) col_names <- colnames(df)
  column_indexes <- split_by_group(col_names, indexes=TRUE)
  res.df <- Reduce(cbind, lapply(column_indexes, function(idxs) rowMeans(df[, idxs]))) # rowSds
  colnames(res.df) <- names(column_indexes)
  res.df
} # works as expected!

get_sample_names <- function(df) {
  sample_names <- grouper(colnames(right_counts_df))
  gsub("(_|-|\ )+$", "", sample_names)
}


###############################################
# heatmap function
###############################################

plot_go_term_heatmap_with_numbers <- function(go_term, df, out_dir, go_df = NULL, genes = NULL) {
  go <- gsub(" ", "_", go_term)
  if (is.null(genes)) genes <- get_go_genes(go_df, go_term)
  size <- if (length(genes) <= 65) 5 else 3
  
  mode <- "scaled_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=get_sample_names(df), 
                   scale=TRUE, 
                   fontsize_row = size,
                   non_scale_log2=TRUE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #

  mode <- "unscaled_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=get_sample_names(df), 
                   scale=FALSE, 
                   fontsize_row = size,
                   non_scale_log2=FALSE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #

  mode <- "unscaledlog2_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=get_sample_names(df), 
                   scale=FALSE, 
                   fontsize_row = size,
                   non_scale_log2=TRUE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #
  
}

###########################################
# which files I are dealing with?
###########################################
right_counts_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907/f7-vs-wt/count-table/nrm-counts-foxd2-mk4-f7-vs-wt.xlsx'

## correction: you have to take new analysis' file!
right_counts_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-foxd2-F7-20200403/f7-vs-wt_foxd2/f7-vs-wt_foxd2/count-table/nrm-counts-foxd2-foxd2-F7-corr-f7-vs-wt_foxd2.xlsx'

###########################################
# read-in each table
###########################################

require(xlsx2dfs)
right_counts_dfs <- xlsx2dfs(right_counts_fpath)
right_counts_df <- right_counts_dfs$`nrm-counts`

###########################################
# which genes I am dealing with?
###########################################

kidney.morphogenesis <- c("Wnt4", "Tgfb1", "Gdnf", "Gata3", "Hs2st1", "Lif", "Vangl2", "Adamts16", "Hoxb7", "Ahr", "Sox8", "Lama5", "Wnk4", "Kif26b", "Fgf1", "Eya1", "Fras1", "Lhx1", "Npnt", "Hnf1b", "Nog", "Pax8", "Pax2", "Wnt6", "Irx2", "Wnt7b", "Sall1", "Irx1")


metanephros.development <- c("Wnt4", "Aph1a", "Gdnf", "Gata3", "Rdh10", "Lif", "Sox8", "Pdgfa", "Kif26b", "Eya1", "Aqp1", "Fras1", "Pdgfb", "Lhx1", "Pou3f3", "Hnf1b", "Pax8", "Pax2", "Irx2", "Wnt7b", "Tfap2a", "Sall1", "Irx1")


wnt.signalling.pathway <- c("Gsk3a", "Arntl", "Mdfic", "Sulf2", "Gid8", "Porcn", "Notch1", "Stk4", "Aes", "Wnt4", "Trpm4", "Rac1", "Nfkb1", "Wls", "Nphp4", "Src", "Dkk3", "Gata3", "Xiap", "Kdm6a", "Atp6ap2", "Ccnd1", "Zfp703", "Bicc1", "Jade1", "Wnt5a", "Wnt9a", "Prkaa2", "Dvl1", "Zbtb33", "Vangl2", "Fam53b", "Egfr", "Ccne1", "Nrarp", "Ror1", "Gprc5b", "Fzd4", "Pitx2", "Fzd6", "Met", "Daam2", "Axin2", "Lgr6", "Prickle1", "Sema5a", "Wnt5b", "Mgat3", "Wnt10a", "Scel", "Nkx2-5", "Folr1", "Cdh1", "Cav1", "Sfrp1", "Mir196a-1", "Hnf1b", "Ptpru", "Nog", "Celsr1", "Wnt6", "Ccdc88c", "Wnt7b", "Celsr2", "Sall1")


out_dir = "/media/josephus/Elements/DEanalysis/heatmap_foxd2-foxd2-F7-20200416/f7-vs-wt_foxd2/heatmaps_GO_genes"
plot_go_term_heatmap_with_numbers("kidney morphogenesis", right_counts_df, out_dir, genes=kidney.morphogenesis)
plot_go_term_heatmap_with_numbers("metanephros development", right_counts_df, out_dir, genes=metanephros.development)
plot_go_term_heatmap_with_numbers("Wnt signaling pathway", right_counts_df, out_dir, genes=wnt.signalling.pathway)

# rsync -rav -e "ssh -p 64168" josephus@132.230.165.168:/media/josephus/Elements/DEanalysis/heatmap_foxd2-foxd2-F7-20200403 ~/Downloads/DEanalysis/
# rsync -rav -e "ssh -p 64168" josephus@132.230.165.168:/media/josephus/Elements/DEanalysis/heatmap_foxd2-foxd2-F7-20200416 ~/Downloads/DEanalysis/

##################################
# with biomaRt
##################################

> genes
$`GO:0060993`
 [1] "Vangl2"   "Wwtr1"    "Nphp3"    "Gcnt4"    "Wnt9b"    "Gcnt3"   
 [7] "Wnt4"     "Agtr2"    "Zmpste24" "Prkx"     "Hnf1b"    "Sox4"    
[13] "Gcnt1"   

$`GO:0001656`
 [1] "Six2"   "Fbn1"   "Bmp7"   "Gdf11"  "Itga8"  "Bcl2"   "Greb1l" "Slit2" 
 [9] "Osr1"   "Nkx3-1" "Ctsh"   "Wnt4"   "Irx1"   "Foxc2"  "Irx2"   "Gli3"  
[17] "Eya1"   "Bmp4"   "Irx3"   "Spry1"  "Pax8"   "Id2"    "Rdh10"  "Aph1a" 
[25] "Fgf10"  "Aph1c"  "Id3"    "Sox17"  "Lhx1"   "Osr2"   "Fgf8"   "Shh"   
[33] "Gdnf"   "Nf1"    "Hoxa11" "Hoxc11"

$`GO:0016055`
  [1] "Tmem88"        "Lgr4"          "Strn"          "Csnk1g1"      
  [5] "Cdc73"         "Rnf43"         "Wif1"          "Gm28635"      
  [9] "Lzts2"         "Csnk1g3"       "Nkd2"          "Daam2"        
 [13] "Axin1"         "Pkd2"          "Csnk1g2"       "Mark1"        
 [17] "Klhl12"        "Apcdd1"        "Wnt10b"        "Ccne1"        
 [21] "Lgr6"          "Axin2"         "Notum"         "Rspo4"        
 [25] "Sost"          "Gsk3b"         "Lmbr1l"        "Slc9a3r1"     
 [29] "Csnk1d"        "Apc"           "Gid8"          "Ldb1"         
 [33] "Fzd10"         "Tspan12"       "Dact1"         "Brd7"         
 [37] "Fzd3"          "Wnt1"          "Vrk3"          "Cul3"         
 [41] "Csnk2a2"       "Nkd1"          "Arl6"          "Rtf1"         
 [45] "Wnt8b"         "Fzd6"          "Nphp3"         "Mesd"         
 [49] "Cpz"           "Grk5"          "Tmem131l"      "Daam1"        
 [53] "Kremen1"       "Csnk1e"        "Csnk2a1"       "Znrf3"        
 [57] "Wnt9b"         "Wnt3"          "Cxxc4"         "Sostdc1"      
 [61] "Csnk1a1"       "Ctnnd1"        "Cyld"          "Ptk7"         
 [65] "Tax1bp3"       "Pias4"         "Ror1"          "Fzd7"         
 [69] "Bcl9"          "Porcn"         "Frat2"         "Celsr2"       
 [73] "Xiap"          "Lrp5"          "Fermt2"        "Tle3"         
 [77] "Ndrg2"         "Gsk3a"         "Amer2"         "Wnt4"         
 [81] "Wnt8a"         "Myc"           "Reck"          "Wnt16"        
 [85] "Zbtb33"        "Zranb1"        "Amer1"         "Fam53b"       
 [89] "Prkaa1"        "Wnt3a"         "Wnt9a"         "Drd2"         
 [93] "Lrp6"          "Ccnd1"         "Vrk2"          "Tnks2"        
 [97] "Pygo1"         "Ubac2"         "Sfrp4"         "Prkaa2"       
[101] "Ctnnbip1"      "Otulin"        "Tnks"          "Fzd5"         
[105] "Rnf138"        "Tcf7l1"        "Chd8"          "Ccar2"        
[109] "Dkk4"          "Sfrp5"         "Fzd1"          "Hhex"         
[113] "Spin1"         "Cdk14"         "Krt6a"         "Dvl3"         
[117] "Wnt7b"         "Nxn"           "Sfrp1"         "Mitf"         
[121] "Tle4"          "Ddx3x"         "Tle1"          "Ryk"          
[125] "Amotl1"        "Usp34"         "Wls"           "Ddb1"         
[129] "Tgfb1i1"       "Tmem198"       "Amer3"         "Senp2"        
[133] "Dixdc1"        "Wnt6"          "Wnt10a"        "Trabd2b"      
[137] "Tcf7l2"        "Amfr"          "Shisa6"        "Kremen2"      
[141] "Dact3"         "Wnt5a"         "Apc2"          "Sox17"        
[145] "Ctnnb1"        "Cpe"           "Zbed3"         "Pitx2"        
[149] "Nlk"           "Adgra2"        "Btrc"          "Ndp"          
[153] "Amotl2"        "Wnt5b"         "Epm2a"         "Ddit3"        
[157] "Etv2"          "Rnf146"        "Lef1"          "Wnt2"         
[161] "Dkk1"          "Fbxw4"         "Frzb"          "Ror2"         
[165] "Vrk1"          "Csnk2b"        "Ctr9"          "Fzd8"         
[169] "Tle5"          "Ccny"          "Tle2"          "Rspo3"        
[173] "Tcf7"          "Fzd4"          "Disc1"         "Dvl1"         
[177] "Dvl2"          "Lrp4"          "Wnt2b"         "Rspo2"        
[181] "Grk6"          "Draxin"        "Ccdc88c"       "Bcl7b"        
[185] "Dkk2"          "Pkd1"          "Fzd2"          "Tnik"         
[189] "Vax2"          "Fzd9"          "Wnt7a"         "Wwox"         
[193] "Cela1"         "Calcoco1"      "Ccn4"          "Wnt11"        
[197] "Paf1"          "Leo1"          "Wdr61"         "Dkk3"         
[201] "Fbxw11"        "2210016L21Rik" "Sfrp2"  

Well, only the Wnt-signalling contains more genes ...
but the others fewer genes ... strange ...
If finding better gene lists



#############################3
# with table
#############################


fpath <- "/media/josephus/archive/go/gene_association.mgi.gz"
go_df <- read.delim(gzfile(fpath), comment.char="!", header=FALSE, sep="\t", quote="")
# V5 is go_id and V3 is gene symbol!

go_df$V3[as.character(go_df$V5) == go_ids["kidney morphogenesis"]]


> go_df$V3[as.character(go_df$V5) == go_ids["kidney morphogenesis"]]
 [1] Agtr2    Ahr      Gcnt1    Gcnt3    Gcnt4    Hnf1b    Nphp3    Nphp3   
 [9] Prkx     Sox4     Tshz3    Vangl2   Wnt4     Wnt9b    Wwtr1    Zmpste24
24836 Levels: 0610005C13Rik 0610006L08Rik 0610009B22Rik ... Zzz3
> go_df$V3[as.character(go_df$V5) == go_ids["metanephros development"]]
 [1] Aph1a  Aph1c  Bcl2   Bcl2   Bmp4   Bmp4   Eya1   Eya1   Fbn1   Foxc2 
[11] Gdf11  Gdnf   Gdnf   Gdnf   Gdnf   Gli3   Greb1l Hoxa11 Hoxc11 Hoxd11
[21] Id2    Itga8  Itga8  Nf1    Osr1   Pax2   Pax2   Pax8   Pds5a  Rdh10 
[31] Robo2  Shh    Six2   Slc5a1 Slit2  Spry1  Wnt4   Wnt4   Wt1    Wt1   
[41] Wt1    Id3    Lhx1   Irx2   Irx3   Irx1   Nkx3-1 Pax8   Tshz3  Osr2  
24836 Levels: 0610005C13Rik 0610006L08Rik 0610009B22Rik ... Zzz3
> go_df$V3[as.character(go_df$V5) == go_ids["Wnt signaling pathway"]]
  [1] 2210016L21Rik Aes           Amer1         Amer2         Amer3        
  [6] Amotl1        Amotl2        Apc           Apc           Apcdd1       
 [11] Arl6          Axin1         Axin1         Axin2         Bcl7b        
 [16] Bcl9          Brd7          Btrc          Calcoco1      Ccar2        
 [21] Ccdc88c       Ccdc88c       Ccnd1         Ccne1         Ccny         
 [26] Cd44          Cdc73         Cdk14         Cela1         Celsr2       
 [31] Chd8          Cpe           Cpz           Csnk1d        Csnk1e       
 [36] Csnk1g1       Csnk1g2       Csnk1g3       Csnk2a2       Csnk2b       
 [41] Ctnnb1        Ctnnbip1      Ctr9          Cul3          Cxxc4        
 [46] Cyld          Daam1         Daam2         Dab2          Dact1        
 [51] Dact3         Ddb1          Ddit3         Ddx3x         Disc1        
 [56] Dixdc1        Dixdc1        Dkk1          Dkk2          Dkk3         
 [61] Dkk4          Draxin        Draxin        Drd2          Dvl1         
 [66] Dvl2          Dvl2          Dvl3          Dvl3          Etv2         
 [71] Fam53b        Fbxw11        Fbxw4         Fermt2        Frat1        
 [76] Frat2         Frat2         Frat2         Frzb          Fzd1         
 [81] Fzd10         Fzd2          Fzd2          Fzd3          Fzd3         
 [86] Fzd4          Fzd4          Fzd4          Fzd5          Fzd6         
 [91] Fzd6          Fzd7          Fzd7          Fzd8          Fzd8         
 [96] Fzd9          Gid8          Grk5          Grk6          Gsk3a        
[101] Gsk3b         Hbp1          Hhex          Hic1          Invs         
[106] Klhl12        Kremen1       Kremen1       Kremen2       Kremen2      
[111] Krt6a         Ldb1          Lef1          Lef1          Leo1         
[116] Lgr4          Lgr6          Lrp4          Lrp5          Lrp6         
[121] Lrrfip2       Lzts2         Macf1         Macf1         Mark1        
[126] Mark2         Mesd          Mir106b       Mir15a        Mir182       
[131] Mir195a       Mir196a-1     Mir196a-2     Mir200a       Mir200b      
[136] Mir20a        Mir221        Mir222        Mir224        Mir27b       
[141] Mir29b-1      Mir32         Mir93         Mitf          Myc          
[146] Ndp           Ndrg2         Nkd1          Nkd2          Nlk          
[151] Notum         Nphp3         Nxn           Otulin        Paf1         
[156] Peg12         Pias4         Pitx2         Porcn         Porcn        
[161] Prkaa1        Prkaa2        Ptk7          Ptk7          Pygo1        
[166] Rnf138        Rnf146        Rnf146        Rnf43         Ror1         
[171] Ror2          Rspo1         Rspo2         Rspo3         Rspo4        
[176] Rtf1          Ryk           Ryk           Senp2         Sfrp1        
[181] Sfrp2         Sfrp4         Sfrp5         Shisa6        Slc9a3r1     
[186] Sost          Sostdc1       Sox17         Spin1         Strn         
[191] Tax1bp3       Tcf7          Tcf7          Tcf7l1        Tcf7l1       
[196] Tcf7l2        Tcf7l2        Tgfb1i1       Tle1          Tle2         
[201] Tle3          Tle4          Tmem131l      Tmem198       Tmem88       
[206] Tnik          Tnks          Tnks          Tnks2         Trabd2b      
[211] Usp34         Vax2          Wdr61         Wif1          Wisp1        
[216] Wls           Wnt1          Wnt10b        Wnt2          Wnt3a        
[221] Wnt5a         Wnt5a         Wnt5a         Wnt6          Wnt6         
[226] Wnt7a         Wnt7b         Wnt9b         Wwox          Xiap         
[231] Zbed3         Zbtb33        Znrf3         Zranb1        Wnt2b        
[236] Wnt5b         Axl           Mst1r         Ryk           Wnt3a        
[241] Wnt3          Csnk1a1       Wnt8a         Mertk         Wnt9a        
[246] Wnt6          Wnt7b         Wnt2          Wnt7a         Wnt8b        
[251] Wnt9b         Tyro3         Wnt4          Wnt16         Met          
[256] Wnt1          Wnt5a         Wnt10a        Wnt10b        Wnt11        
24836 Levels: 0610005C13Rik 0610006L08Rik 0610009B22Rik ... Zzz3



go_tables = list( "kidney morphogenesis" = as.character(go_df$V3[as.character(go_df$V5) == go_ids["kidney morphogenesis"]]),
                  "metanephros development" = as.character(go_df$V3[as.character(go_df$V5) == go_ids["metanephros development"]]),
                  "Wnt signaling pathway" = as.character(go_df$V3[as.character(go_df$V5) == go_ids["Wnt signaling pathway"]]))

genes

go_bps = list( "kidney morphogenesis" = kidney.morphogenesis,
                  "metanephros development" = metanephros.development,
                  "Wnt signaling pathway" = wnt.signalling.pathway)

fused_genes <- Map(function(...) unique(unlist(...)), go_tables, genes, go_bps)

save.image(file="~/foxd2.analysis.20200403.rdata")
####################################
# todo:
####################################

Vielleicht alle diese Gene aufsammeln und dann ein heatmap generieren!



# kill 5181
# kill 31868
# kill 3987
# 
#  8051
#  
# 4273
# 10191
# 4965
# 19202
# kill 4390
# kill 6251








################################################################################

Hi Josephus,


hier nochmal das directory für die RNA seq Daten für die heatmap:


RNAseq_201907_analysis, 

      foxd2-F7-mk4-201907, 

          f7-vs-wt


wäre es möglich 3 verschiedene maps zu machen? alles aus  f7-vs-wt? 


im skype call haben wir besprochen, dass wir gerne eine oder  heatmap auf GO terms basierend ins Papier machen würden und eine die aus Genen besteht die wir aus anderen verschiedenen gründen relevant für den Phänotyp halten. GO term Listen wären 2 gut, eine der terms die am prominentesten runter ist und einer hoch:


    Go term heatmap: aus der Liste der GO terms in f7-vs-wt folder wäre GO:0001822 (kidney development) gut. Das sind diese Gene:Mmp17/Smad9/Fgf1/Col4a4/Tfap2a/Sim1/Wnt2b/Adamts1/Tgfb2/Aqp1/Cys1/Pax2/Npnt/Egr1/Agt/Lrp4/Col4a3/Wnt4/Fgfr2/Id3/Fras1/Gli2/Pygo1/Enpep



2. aus der Liste der GO terms in f7-vs-wt folder GO:0001822 (regulation of MAPK)
Gene sind Sash1/Hgf/Samd5/Flt1/Ephb2/Tgfa/Il1rn/Fgfr1/Dusp9/Dusp4/Sfrp1/Ghr/Ghrl/Cd74/Fgd4/Fzd8/Wnt7b/Pdgfrb/Spred2/Mif/Vegfa/Uchl1/Tiam1/Pdgfb/Ptpn6

3. Genliste nicht GO term basiert sonderen aus Interesse für die Gene:
PAX2, PAX8, wnt6, wnt4, wnt7a, wnt7b, Fras1, Sall1, EYA1 Rab26b, FOXD1, ETV4, ETV5, LAMA5, LRP2, FGFR2, FAT2, RET, MYCN, GFRA1, ITGB3, MAP3K9, MAP3K5, MAPK10,  MAPT, GDNF

als Kontrollgene kann man 
für hetamap 1 schauen ob die gene in f7 niedriger sind als für wt. 
für heatmap 2 sollte alles für f7 höher sein als für wt. 
für heatmap 3 sollten Pax2 und wnt4 in f7 niedriger sein als für wt. 

die Gene aus heatmap 

das wäre total super, danke! 

wenn Du Fragen hast welche Daten ruf mich an, am besten am handy, ich bin meist nicht im Büro. 




################################################################################

# work
# source activate newR

right_go_analysis_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907/f7-vs-wt/GO/DE_sig_foxd2-mk4_f7-vs-wt_190801155819_0.05_1_over_GO_BP.xlsx'
right_counts_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907/f7-vs-wt/count-table/nrm-counts-foxd2-mk4-f7-vs-wt.xlsx'


kidney.development <- c("Mmp17", "Smad9", "Fgf1", "Col4a4", "Tfap2a", "Sim1", "Wnt2b", "Adamts1", "Tgfb2", "Aqp1", "Cys1", "Pax2", "Npnt", "Egr1", "Agt", "Lrp4", "Col4a3", "Wnt4", "Fgfr2", "Id3", "Fras1", "Gli2", "Pygo1", "Enpep")

regulation.of.mapk <- c("Sash1", "Hgf", "Samd5", "Flt1", "Ephb2", "Tgfa", "Il1rn", "Fgfr1", "Dusp9", "Dusp4", "Sfrp1", "Ghr", "Ghrl", "Cd74", "Fgd4", "Fzd8", "Wnt7b", "Pdgfrb", "Spred2", "Mif", "Vegfa", "Uchl1", "Tiam1", "Pdgfb", "Ptpn6")

capitalize <- Vectorize(function(s) paste0(toupper(substring(s, 1, 1)), tolower(substring(s, 2))))
goi <- capitalize(c("PAX2", "PAX8", "wnt6", "wnt4", "wnt7a", "wnt7b", "Fras1", "Sall1", "EYA1", "Rab26", "Rab26os", "FOXD1", "ETV4", "ETV5", "LAMA5", "LRP2", "FGFR2", "FAT2", "RET", "MYCN", "GFRA1", "ITGB3", "MAP3K9", "MAP3K5", "MAPK10", "MAPT", "GDNF"))
# test if goi names are all present in list

# there is no Rab26b only Rab26 or Rab26os!

###############################################
# generate for the genes a heatmap
###############################################

###############################################
# helper functions
###############################################
time.now <- function() format(Sys.time(), "%y%m%d%H%M%S")

common.string <- function(s1, s2, acc="") {
  # Return common beginning string
  first.s1 = substr(s1, 1, 1)
  first.s2 = substr(s2, 1, 1)
  if (s1 == "" && s2 == "" || first.s1 != first.s2) {
    acc
  } else  { # s1[1] == s2[1]
    common.string(substr(s1, 2, nchar(s1)), substr(s2, 2, nchar(s2)), paste0(acc, first.s1))
  }
} # works!

cdr <- function(vec) if (length(vec) < 2) c() else vec[2:length(vec)]
car <- first <- function(vec) vec[1]
cadr <- second <- function(vec) vec[2]
cons <- function(x, vec) c(x, vec)
nreverse <- rev
null <- function(x) length(x) == 0

grouper <- function(vec, acc=c(), last.common=c(), prev="") {
  if (null(vec)) {
    nreverse(acc)
  } else if (length(vec) == 1) {
    grouper(cdr(vec), cons(last.common, acc))
  } else {
    next.common <- common.string(first(vec), second(vec))
    if (null(last.common)) {
      grouper(cdr(vec), cons(next.common, acc), next.common, car(vec))
    } else {
      prev.common <- common.string(prev, car(vec))
      if (last.common == prev.common) {
        grouper(cdr(vec), cons(prev.common, acc), prev.common, car(vec))
      } else {
        grouper(cdr(vec), cons(next.common, acc), next.common, car(vec))
      }
    }
  }
} # works!

split_by_group <- function(vec, indexes=FALSE) {
  dfs <- split(data.frame(x= if (indexes) 1:length(vec) else vec, 
                          stringsAsFactors=FALSE), 
               grouper(vec))
  lapply(dfs, function(df) df$x)
}

average_df <- function(df, col_names=NULL) {
  # averages df by similarity of columnnames
  if (is.null(col_names)) col_names <- colnames(df)
  column_indexes <- split_by_group(col_names, indexes=TRUE)
  res.df <- Reduce(cbind, lapply(column_indexes, function(idxs) rowMeans(df[, idxs]))) # rowSds
  colnames(res.df) <- names(column_indexes)
  res.df
} # works as expected!

get_sample_names <- function(df) {
  sample_names <- grouper(colnames(df))
  gsub("(_|-|\ )+$", "", sample_names)
}

# rtrim <- Vectorize(function(s, pattern) gsub(paste0(pattern, "$"), "", s))
# generate_names <- function(df) rtrim(grouper(colnames(df)), pattern="(_|-|\ )")


###############################################
# heatmap function
###############################################
require(pheatmap)
generate_heatmap <- function(df, 
                             col_names=NULL,
                             scale=TRUE,
                             non_scale_log2=FALSE,
                             color.1 = "black",
                             color.2 = "red",
                             cluster_rows = F,
                             col_steps = 100,
                             fontsize_row = 5,
                             display_numbers = FALSE,
                             number_color = "white",
                             fontsize_number = 5,
                             outname = NULL) {
  # `outname` should be without ending!
  df.avg <- average_df(df, col_names = col_names)
  if (scale) {
    scaledata <- t(scale(t(df.avg)))
    scaledata <- scaledata[complete.cases(scaledata), ]
  } else {
    scaledata <- if (non_scale_log2) log2(df.avg + .1) else df.avg
  }
  # prepare heatmap
  colfunc <- colorRampPalette(c(color.1, color.2))
  
  if (is.null(outname)) {
    pheatmap(scaledata,
             cluster_rows = cluster_rows,
             cluster_cols = F,
             cellwidth = 40,
             col = colfunc(col_steps),
             fontsize_row = fontsize_row,
             border_color = NA,
             display_numbers = display_numbers,
             number_color = number_color,
             fontsize_number = fontsize_number)
  } else {
  # print heatmap
    {
      setEPS()
      postscript(paste0(outname, ".eps"))
      # svg(paste0(outname, ".svg"))
      pheatmap(scaledata,
               cluster_rows = cluster_rows,
               cluster_cols = F,
               cellwidth = 40,
               col = colfunc(col_steps),
               fontsize_row = fontsize_row,
               border_color = NA,
               display_numbers = display_numbers,
               number_color = number_color,
               fontsize_number = fontsize_number)
      dev.off()
    }

    {
      svg(paste0(outname, ".svg"))
      pheatmap(scaledata,
               cluster_rows = cluster_rows,
               cluster_cols = F,
               cellwidth = 40,
               col = colfunc(col_steps),
               fontsize_row = fontsize_row,
               border_color = NA,
               display_numbers = display_numbers,
               number_color = number_color,
               fontsize_number = fontsize_number)
      dev.off()
    }
  #   # print counts
  #   res.all.gois <- gois[gois %in% rownames(res$all.names.srt)]
  #   res.up.gois <- gois[gois %in% rownames(res$up.names.srt)]
  #   res.down.gois <- gois[gois %in% rownames(res$down.names.srt)]
  #   
  #   res.gois <- list(cnts.nrm[gois.present, ], cnts.gois, res.all.gois, res.up.gois, res.down.gois)
  #   names(res.gois) <- c("goi.counts", "goi.counts.avg", "goi.DE.all", "goi.DE.up", "goi.DE.down")
  #   write.dfs(res.gois, paste0(outname, ".xlsx"))
    dfs2xlsx(withNames("goi_all", scaledata),
            paste0(outname, ".xlsx"))
 }
}


plot_go_term_heatmap_with_numbers <- function(go_term, df, out_dir, go_df = NULL, genes = NULL) {
  go <- gsub(" ", "_", go_term)
  if (is.null(genes)) genes <- get_go_genes(go_df, go_term)
  size <- if (length(genes) <= 65) 5 else 3
  
  mode <- "scaled_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=get_sample_names(df), 
                   scale=TRUE, 
                   fontsize_row = size,
                   non_scale_log2=TRUE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #

  mode <- "unscaled_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=get_sample_names(df), 
                   scale=FALSE, 
                   fontsize_row = size,
                   non_scale_log2=FALSE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #

  mode <- "unscaledlog2_f7_vs_wt"
  fname <- paste0(go, "_", mode, "_", time.now())
  out_fpath <- file.path(out_dir, fname)
  dir.create(dirname(out_fpath), recursive=TRUE, showWarnings=FALSE)
  generate_heatmap(df[genes, ],
                   col_names=get_sample_names(df), 
                   scale=FALSE, 
                   fontsize_row = size,
                   non_scale_log2=TRUE,
                   display_numbers = TRUE,
                   fontsize_number = size,
                   outname = out_fpath) #
  
}




###########################################
# read-in each table
###########################################

require(xlsx2dfs)
right_go_dfs <- xlsx2dfs(right_go_analysis_fpath)
right_counts_dfs <- xlsx2dfs(right_counts_fpath)

right_go_df <- right_go_dfs[["all.l2FC.srt"]]
right_counts_df <- right_counts_dfs$`nrm-counts`


all(goi %in% rownames(right_counts_df)) # True

out_dir = "/media/josephus/Elements/DEanalysis/heatmap_foxd2-F7-mk4-201907/f7-vs-wt/heatmaps_foxd2-mk4_f7-vs-wt"
plot_go_term_heatmap_with_numbers("kidney development", right_counts_df, out_dir, genes=kidney.development)
plot_go_term_heatmap_with_numbers("regulation of mapk", right_counts_df, out_dir, genes=regulation.of.mapk)
plot_go_term_heatmap_with_numbers("gois", right_counts_df, out_dir, genes=goi)

rsync -rav -e "ssh -p 64168" josephus@132.230.165.168:/media/josephus/Elements/DEanalysis/heatmap_foxd2-F7-mk4-201907/ ~/Downloads/DEanalysis/heatmap_foxd2-F7-mk4-201907/


## best is because of the error of sample names to run the entire analysis one day again!



Verwendet wurde das File:
right_counts_fpath <- '/media/josephus/Elements/DEanalysis/foxd2-F7-mk4-201907/f7-vs-wt/count-table/nrm-counts-foxd2-mk4-f7-vs-wt.xlsx'
(sample names korrigiert)

das zugehörige Metafile sieht so aus:
('num' bedeutet 'numerator': zaehler im Bruch (A/B gelesen A vs. B -> A)
 'denom' bedeutet 'denominator': nenner im Bruch (A/B gelesen A vs. B -> B)
also wurde hier mk4-F7 vs. mk4-ctrl untersucht)

nano meta-foxd2-201907p-mk4.txt

sampleName	fileName	condition	testing
wt_mk4_1	/media/josephus/archive_big/count/foxd2_201907p/external-MN-11-mK4-ctrl-rep1-16878-sym-fcount.tab	wt	denom
wt_mk4_2	/media/josephus/archive_big/count/foxd2_201907p/external-MN-12-mK4-ctrl-rep2-16877-sym-fcount.tab	wt	denom
wt_mk4_3	/media/josephus/archive_big/count/foxd2_201907p/external-MN-13-mK4-ctrl-rep3-16876-sym-fcount.tab	wt	denom
wt_mk4_4	/media/josephus/archive_big/count/foxd2_201907p/external-MN-14-mK4-ctrl-rep4-16875-sym-fcount.tab	wt	denom
wt_mk4_5	/media/josephus/archive_big/count/foxd2_201907p/external-MN-15-mK4-ctrl-rep5-16874-sym-fcount.tab	wt	denom
F7_mk4_1	/media/josephus/archive_big/count/foxd2_201907p/external-MN-16-mK4-F7-rep1-16873-sym-fcount.tab	f7	num
F7_mk4_2	/media/josephus/archive_big/count/foxd2_201907p/external-MN-17-mK4-F7-rep2-16872-sym-fcount.tab	f7	num
F7_mk4_3	/media/josephus/archive_big/count/foxd2_201907p/external-MN-18-mK4-F7-rep3-16871-sym-fcount.tab	f7	num
F7_mk4_4	/media/josephus/archive_big/count/foxd2_201907p/external-MN-19-mK4-F7-rep4-16870-sym-fcount.tab	f7	num
F7_mk4_5	/media/josephus/archive_big/count/foxd2_201907p/external-MN-20-mK4-F7-rep5-16869-sym-fcount.tab	f7	num

Es wurden die Gene geplottet:

kidney.development <- c("Mmp17", "Smad9", "Fgf1", "Col4a4", "Tfap2a", "Sim1", "Wnt2b", "Adamts1", "Tgfb2", "Aqp1", "Cys1", "Pax2", "Npnt", "Egr1", "Agt", "Lrp4", "Col4a3", "Wnt4", "Fgfr2", "Id3", "Fras1", "Gli2", "Pygo1", "Enpep")

regulation.of.mapk <- c("Sash1", "Hgf", "Samd5", "Flt1", "Ephb2", "Tgfa", "Il1rn", "Fgfr1", "Dusp9", "Dusp4", "Sfrp1", "Ghr", "Ghrl", "Cd74", "Fgd4", "Fzd8", "Wnt7b", "Pdgfrb", "Spred2", "Mif", "Vegfa", "Uchl1", "Tiam1", "Pdgfb", "Ptpn6")

capitalize <- Vectorize(function(s) paste0(toupper(substring(s, 1, 1)), tolower(substring(s, 2))))
goi <- capitalize(c("PAX2", "PAX8", "wnt6", "wnt4", "wnt7a", "wnt7b", "Fras1", "Sall1", "EYA1", "Rab26", "Rab26os", "FOXD1", "ETV4", "ETV5", "LAMA5", "LRP2", "FGFR2", "FAT2", "RET", "MYCN", "GFRA1", "ITGB3", "MAP3K9", "MAP3K5", "MAPK10", "MAPT", "GDNF"))
# goi: gene of interest

Jede Gengruppe wurde 3x geplotet: 
_scaled_        zentriert um den Mittelwert und std dev als 1 gesetzt
_unscaled_      Rohwerte
_unscaledlog2_  log2 der Rohwerte um Unterschiede auch im kleinen Bereich besser zu sehen








