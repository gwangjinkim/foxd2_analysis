
# 190831.rdownstream.fox2.analsyis.R


##############################
# fonts problem
##############################

# library(extrafont)
# extrafont::font_import()


####################################################################
# source 
####################################################################
# this loads the functions necessary

computer <- "ss"

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



data_names  <-  c("foxd2-mk3",
                  "foxd2-mk4",
                  "foxd2-mk4_foxd2",
                  "foxd2-mk3_mk4_wts",
                  "foxd2-C3-F7",
                  "foxd2-empty-mk4",
                  "foxd2-empty-mk3",
                  "foxd2-foxd2-F7",
                  "foxd2-foxd2-C3")
outdir_names <- c("foxd2-C3-mk3-201907",
                  "foxd2-F7-mk4-201907",
                  "foxd2-foxd2-empty-201907",
                  "foxd2-mk3-mk4-wts-201907",
                  "foxd2-C3-F7-201907",
                  "foxd2-empty-wt-mk4-201907",
                  "foxd2-empty-wt-mk3-201907",
                  "foxd2-foxd2-F7-201907",
                  "foxd2-foxd2-C3-201907")
meta_fnames <- c("meta-foxd2-201907p-mk3.txt", 
                "meta-foxd2-201907p-mk4.txt", 
                "meta-foxd2-201907p-mk4_foxd2.txt", 
                "meta-foxd2-201907p-mk3-mk4-wts.txt", 
                "meta-foxd2-201907p-mk3_C3-mk4_F7-mutants.txt", 
                "meta-foxd2-201907p-mk4_empty-mk4.txt", 
                "meta-foxd2-201907p-mk4_empty-mk3.txt", 
                "meta-foxd2-201907p-mk4_foxd2-mk4_F7.txt", 
                "meta-foxd2-201907p-mk4_foxd2-mk4_C3.txt")

####################################################################
# DE settings
####################################################################
# set here tresholds
# important are lFC (log2FoldChange for DE analysis 1 or 2.5 we use)

alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 1       # log2FC limit
GO.alpha <- 0.01
GO.FDR   <- 0.05
GO.N     <- 30
GO.cl.N  <- 10

########################
# gskb install
########################

# wget https://bioconductor.org/packages/release/data/experiment/src/contrib/gskb_1.16.0.tar.gz
# install.packages("~/R_sources/gskb_1.16.0.tar.gz", type="source", dependencies=TRUE)


if (!require(foreach)) {
  install.packages("foreach")
  require(foreach)
}

if (!require(doParallel)) {
  install.packages("doParallel")
  require(doParallel)
}

numCores <- 30
registerDoParallel(numCores)

foreach (i = 1:9) %dopar% {
  
  #########################################
  # set variables here which are common to all runs
  #########################################
  
  indirpath <- ""
  
  #########################################
  # define variables specific for each run
  #########################################
  
  dataname <- data_names[i]
  outdirpath <- file.path(out_dir, outdir_names[i])
  dir.create(outdirpath, recursive=TRUE)
  metapath <-   file.path(meta_dir, meta_fnames[i])
  
  ####################################################################
  # get meta information
  ####################################################################
  # this section reads in the meta file and makes it available
  # no need to change anything here
  
  meta.df <- read.table(metapath, sep = '\t',
                        header = TRUE,
                        stringsAsFactors = FALSE)
  meta.df$condition <- factor(meta.df$condition, # ensures preferred order
                              levels = unique(meta.df$condition))
  
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')
  
  ####################################################################
  # inferred paths
  ####################################################################
  # prepares subfolders - no changes to make
  
  outdirpath <- file.path(outdirpath, core.name(meta.df))
  meta.outdir <- file.path(outdirpath, "meta")
  cnts.outdir <- file.path(outdirpath, "count-table")
  DE.outdir <- file.path(outdirpath, "DE-table")
  plot.outdir <- file.path(outdirpath, "glimma-plots")
  hm.outdir <- file.path(outdirpath, "heatmap")
  hm.final.outdir <- file.path(hm.outdir, "final")
  hm.up.outdir <- file.path(hm.outdir, "up")
  hm.down.outdir <- file.path(hm.outdir, "down")
  hm.all.outdir <- file.path(hm.outdir, "all")
  hm.up.list.outdir <- file.path(hm.outdir, "up-list")
  hm.down.list.outdir <- file.path(hm.outdir, "down-list")
  hm.all.list.outdir <- file.path(hm.outdir, "all-list")
  GO.outdir <- file.path(outdirpath, "GO")
  pathway.outdir <- file.path(outdirpath, "pathway")
  goi.outdir <- file.path(outdirpath, "goi")
  lgenes.outdir <- file.path(outdirpath, "listgenes")
  gskb.outdir <- file.path(outdirpath, "gskb")
  PCA.outdirpath <- file.path(outdirpath, "pca")
  
  #####################################################################
  # ensure existence of output paths
  #####################################################################
  # creates  subfolders - no changes to make
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.all.list.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=lgenes.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=gskb.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=PCA.outdirpath, recursive = TRUE, showWarnings = FALSE)
  
  # https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
  #######################################################################
  # Create meta table
  #######################################################################
  # creates the meta file in subfolder 'meta' - no changes to make
  
  write.table(meta.df,
              file = file.path(meta.outdir, basename(metapath)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  #######################################################################
  # Create DESeq2.obj and DESeq2.obj.disp
  #######################################################################
  # does the DESeq2 analysis - no changes to make
  
  DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
  DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)
  
  #######################################################################
  # Create the count tables
  #######################################################################
  # prints out count tables - no changes to make
  
  cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = FALSE, averaged = FALSE,
                        sheetName = "raw.all")
  cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = TRUE, averaged = FALSE,
                        sheetName = "normalized.all")
  cnts.avg.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                            dataname = dataname,
                            printp = TRUE, normalized = TRUE, averaged = TRUE,
                            sheetName = "avg.normalized.all")

  #######################################################################
  # Create DE table
  #######################################################################
  # prints out DE tables - no changes to make
  
  res <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                    dataname = dataname,
                    printp = TRUE)
  resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                       dataname = dataname,
                       printp = TRUE,
                       filterp = TRUE, alpha = alpha, lFC = lFC)
  print.cnts.DE.sortings(cnts.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create plots
  #######################################################################
  # creates interactive plots in folder "glimma" - no changes to make
  
  meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
                dataname, top=300, launch = FALSE)
  
  meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
                dataname, alpha = alpha, lFC = lFC,
                launch = FALSE)
  
  #######################################################################
  # Create heatmap
  #######################################################################
  # creates heatmap with scaled values in 'heatmap' folder
  # you can change 'k =' to what you want. e.g. 12, 24, etc.
  # set the seed for reproducibility
  
  set.seed(123) # for reproducibility of heatmaps set seeds to a value
  try(scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm$`nrm-counts-avg`, resSig,
                                              outdirpath = hm.all.outdir, selected.genes = NULL, name.add = "all",
                                              dataname = dataname,
                                              k = 8, printp=TRUE,
                                              alpha=alpha, lFC=lFC, filterp=TRUE,
                                              xlsxp=TRUE, csvp=FALSE, tsvp=FALSE))
  
  
  #######################################################################
  # Create GO
  #######################################################################
  # GO analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # inDir is outdirpath
  # outBase is GO.outdir
  last <- function(l) l[length(l)]
  DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))

  
  try(do_all_GO_enrichment(DEfpath, GO.outdir))
  
  #######################################################################
  # KEGG
  #######################################################################
  # KEGG analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # DEfpath is like above
  # outBase <- pathway.outdir
  try(do_KEGG_enrichment(DEfpath, pathway.outdir))

  #######################################################################
  # gskb
  #######################################################################
  # gskb analysis for mouse - actually no changes needed
  # but if you don't want it (quicker analysis) - comment this section out
  
  # DEfpath is like above
  # outBase <- gskb.outdir
  try(do_all_gskb_enrichments(DEfpath,
                          gskb.outdir,
                          all_ez_grouped))
  
  #######################################################################
  # nicer MDS/PCA plots using plotly
  #######################################################################
  # creates nicer non-interactive PCA plots (for publication)
  # - no changes needed
  # - but if PCA plot not important - comment this section out
  
  # request DE folder for significant genes path in previous results
  last <- function(l) l[length(l)]
  cntsDEfpath <- last(dir(DE.outdir, pattern = "DE-cnts-sig", full.names = TRUE))
  
  try(meta2plyMDSplot(meta.df, 
                  DESeq2.obj.disp, 
                  PCA.outdirpath,
                  dataname, 
                  top=500, 
                  launch = FALSE))
  
  try(meta2plyMDSplotDE(cntsDEfpath = cntsDEfpath, 
                    outdirpath = PCA.outdirpath,
                    dataname = dataname, 
                    corename = core.name(meta.df), 
                    top=500, 
                    launch = FALSE))
  
  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)

}



for (i in 1:6) {
  
  #########################################
  # set variables here which are common to all runs
  #########################################
  
  indirpath <- ""
  
  #########################################
  # define variables specific for each run
  #########################################
  
  # here, I have 2 analysis
  if (i == 1) {
    dataname <- "dyn2h1ko-ud"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-ud-201901"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-ud.txt"
  }
  
  if (i == 2) {
    dataname <- "ift-ud"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift-201901"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-ud.txt"
  }
  
  if (i == 3) {
    dataname <- "dyn2h1ko-ift-ud"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift-dft140ko-201901"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dft140ko-ud.txt"
  }
  
  if (i == 4) {
    dataname <- "dyn2h1ko-15d"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-15d-201901"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-15d.txt"
  }

  if (i == 5) {
    dataname <- "ift-15d"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift-15d-201901"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-15d.txt"
  }
  
  
  if (i == 6) {
    dataname <- "dyn2h1ko-ift-15d"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-ift-15d-201901"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dyn2h1ko-15d.txt"
  }


  
  ####################################################################
  # get meta information
  ####################################################################
  # this section reads in the meta file and makes it available
  # no need to change anything here
  
  meta.df <- read.table(metapath, sep = '\t',
                        header = TRUE,
                        stringsAsFactors = FALSE)
  meta.df$condition <- factor(meta.df$condition, # ensures preferred order
                              levels = unique(meta.df$condition))
  
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')
  
  ####################################################################
  # inferred paths
  ####################################################################
  # prepares subfolders - no changes to make
  
  outdirpath <- file.path(outdirpath, core.name(meta.df))
  meta.outdir <- file.path(outdirpath, "meta")
  cnts.outdir <- file.path(outdirpath, "count-table")
  DE.outdir <- file.path(outdirpath, "DE-table")
  plot.outdir <- file.path(outdirpath, "glimma-plots")
  hm.outdir <- file.path(outdirpath, "heatmap")
  hm.final.outdir <- file.path(hm.outdir, "final")
  hm.up.outdir <- file.path(hm.outdir, "up")
  hm.down.outdir <- file.path(hm.outdir, "down")
  hm.all.outdir <- file.path(hm.outdir, "all")
  hm.up.list.outdir <- file.path(hm.outdir, "up-list")
  hm.down.list.outdir <- file.path(hm.outdir, "down-list")
  hm.all.list.outdir <- file.path(hm.outdir, "all-list")
  GO.outdir <- file.path(outdirpath, "GO")
  pathway.outdir <- file.path(outdirpath, "pathway")
  goi.outdir <- file.path(outdirpath, "goi")
  lgenes.outdir <- file.path(outdirpath, "listgenes")
  gskb.outdir <- file.path(outdirpath, "gskb")
  PCA.outdirpath <- file.path(outdirpath, "pca")
  
  #####################################################################
  # ensure existence of output paths
  #####################################################################
  # creates  subfolders - no changes to make
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.all.list.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=lgenes.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=gskb.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=PCA.outdirpath, recursive = TRUE, showWarnings = FALSE)
  
  # https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
  #######################################################################
  # Create meta table
  #######################################################################
  # creates the meta file in subfolder 'meta' - no changes to make
  
  write.table(meta.df,
              file = file.path(meta.outdir, basename(metapath)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  #######################################################################
  # Create DESeq2.obj and DESeq2.obj.disp
  #######################################################################
  # does the DESeq2 analysis - no changes to make
  
  DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
  DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)
  
  #######################################################################
  # Create the count tables
  #######################################################################
  # prints out count tables - no changes to make
  
  cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = FALSE, averaged = FALSE,
                        sheetName = "raw.all")
  cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = TRUE, averaged = FALSE,
                        sheetName = "normalized.all")
  cnts.avg.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                            dataname = dataname,
                            printp = TRUE, normalized = TRUE, averaged = TRUE,
                            sheetName = "avg.normalized.all")
  
  #######################################################################
  # Create DE table
  #######################################################################
  # prints out DE tables - no changes to make
  
  res <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                    dataname = dataname,
                    printp = TRUE)
  resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                       dataname = dataname,
                       printp = TRUE,
                       filterp = TRUE, alpha = alpha, lFC = lFC)
  print.cnts.DE.sortings(cnts.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create plots
  #######################################################################
  # creates interactive plots in folder "glimma" - no changes to make
  
  meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
                dataname, top=300, launch = FALSE)
  
  meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
                dataname, alpha = alpha, lFC = lFC,
                launch = FALSE)
  
  #######################################################################
  # Create heatmap
  #######################################################################
  # creates heatmap with scaled values in 'heatmap' folder
  # you can change 'k =' to what you want. e.g. 12, 24, etc.
  # set the seed for reproducibility
  
  set.seed(123) # for reproducibility of heatmaps set seeds to a value
  scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm$`nrm-counts-avg`, resSig,
                                              outdirpath = hm.all.outdir, selected.genes = NULL, name.add = "all",
                                              dataname = dataname,
                                              k = 8, printp=TRUE,
                                              alpha=alpha, lFC=lFC, filterp=TRUE,
                                              xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)
  
  
  #######################################################################
  # Create GO
  #######################################################################
  # GO analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # inDir is outdirpath
  # outBase is GO.outdir
  last <- function(l) l[length(l)]
  DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))

  
  do_all_GO_enrichment(DEfpath, GO.outdir)
  
  #######################################################################
  # KEGG
  #######################################################################
  # KEGG analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # DEfpath is like above
  # outBase <- pathway.outdir
  do_KEGG_enrichment(DEfpath, pathway.outdir)

  #######################################################################
  # gskb
  #######################################################################
  # gskb analysis for mouse - actually no changes needed
  # but if you don't want it (quicker analysis) - comment this section out
  
  # DEfpath is like above
  # outBase <- gskb.outdir
  do_all_gskb_enrichments(DEfpath,
                          gskb.outdir,
                          all_ez_grouped)
  
  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)
}

























for (i in 1:6) {
  
  #########################################
  # set variables here which are common to all runs
  #########################################
  
  indirpath <- ""
  
  #########################################
  # define variables specific for each run
  #########################################
  
  # here, I have 2 analysis
  if (i == 1) {
    dataname <- "dyn2h1ko-ud-4rm-3rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-ud-201901-4rm-3rm-test"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-ud-4rm-3rm.txt"
  }
  
  if (i == 2) {
    dataname <- "ift-ud-15rm-3rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift-201901-15rm-3rm-test"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-ud-15rm-3rm.txt"
  }
  
  if (i == 3) {
    dataname <- "dyn2h1ko-ift-ud-15rm-4rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift-dft140ko-201901-15rm-4rm-test"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dft140ko-ud-15rm-4rm.txt"
  }
  
  if (i == 4) {
    dataname <- "dyn2h1ko-15d-45rm-12rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-15d-201901-45rm-12rm-test"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-15d-45rm-12rm.txt"
  }

  if (i == 5) {
    dataname <- "ift-15d-5rm-12rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift-15d-201901-5rm-12rm-test"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-15d-5rm-12rm.txt"
  }
  
  
  if (i == 6) {
    dataname <- "dyn2h1ko-ift-15d-5rm-45rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-ift-15d-201901-5rm-45rm-test"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dyn2h1ko-15d-5rm-45rm.txt"
  }


  
  ####################################################################
  # get meta information
  ####################################################################
  # this section reads in the meta file and makes it available
  # no need to change anything here
  
  meta.df <- read.table(metapath, sep = '\t',
                        header = TRUE,
                        stringsAsFactors = FALSE)
  meta.df$condition <- factor(meta.df$condition, # ensures preferred order
                              levels = unique(meta.df$condition))
  
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')
  
  ####################################################################
  # inferred paths
  ####################################################################
  # prepares subfolders - no changes to make
  
  outdirpath <- file.path(outdirpath, core.name(meta.df))
  meta.outdir <- file.path(outdirpath, "meta")
  cnts.outdir <- file.path(outdirpath, "count-table")
  DE.outdir <- file.path(outdirpath, "DE-table")
  plot.outdir <- file.path(outdirpath, "glimma-plots")
  hm.outdir <- file.path(outdirpath, "heatmap")
  hm.final.outdir <- file.path(hm.outdir, "final")
  hm.up.outdir <- file.path(hm.outdir, "up")
  hm.down.outdir <- file.path(hm.outdir, "down")
  hm.all.outdir <- file.path(hm.outdir, "all")
  hm.up.list.outdir <- file.path(hm.outdir, "up-list")
  hm.down.list.outdir <- file.path(hm.outdir, "down-list")
  hm.all.list.outdir <- file.path(hm.outdir, "all-list")
  GO.outdir <- file.path(outdirpath, "GO")
  pathway.outdir <- file.path(outdirpath, "pathway")
  goi.outdir <- file.path(outdirpath, "goi")
  lgenes.outdir <- file.path(outdirpath, "listgenes")
  gskb.outdir <- file.path(outdirpath, "gskb")
  PCA.outdirpath <- file.path(outdirpath, "pca")
  
  #####################################################################
  # ensure existence of output paths
  #####################################################################
  # creates  subfolders - no changes to make
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.all.list.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=lgenes.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=gskb.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=PCA.outdirpath, recursive = TRUE, showWarnings = FALSE)
  
  # https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
  #######################################################################
  # Create meta table
  #######################################################################
  # creates the meta file in subfolder 'meta' - no changes to make
  
  write.table(meta.df,
              file = file.path(meta.outdir, basename(metapath)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  #######################################################################
  # Create DESeq2.obj and DESeq2.obj.disp
  #######################################################################
  # does the DESeq2 analysis - no changes to make
  
  DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
  DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)
  
  #######################################################################
  # Create the count tables
  #######################################################################
  # prints out count tables - no changes to make
  
  cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = FALSE, averaged = FALSE,
                        sheetName = "raw.all")
  cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = TRUE, averaged = FALSE,
                        sheetName = "normalized.all")
  cnts.avg.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                            dataname = dataname,
                            printp = TRUE, normalized = TRUE, averaged = TRUE,
                            sheetName = "avg.normalized.all")
  
  #######################################################################
  # Create DE table
  #######################################################################
  # prints out DE tables - no changes to make
  
  res <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                    dataname = dataname,
                    printp = TRUE)
  resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                       dataname = dataname,
                       printp = TRUE,
                       filterp = TRUE, alpha = alpha, lFC = lFC)
  print.cnts.DE.sortings(cnts.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create plots
  #######################################################################
  # creates interactive plots in folder "glimma" - no changes to make
  
  meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
                dataname, top=300, launch = FALSE)
  
  meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
                dataname, alpha = alpha, lFC = lFC,
                launch = FALSE)
  
  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)
}








for (i in c(1, 2, 4, 5, 6)) {
  
  #########################################
  # set variables here which are common to all runs
  #########################################
  
  indirpath <- ""
  
  #########################################
  # define variables specific for each run
  #########################################
  
  # here, I have 2 analysis
  if (i == 1) {
    dataname <- "dyn2h1ko-ud-4rm-32rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-ud-201901-4rm-32rm-test-1"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-ud-4rm-32rm.txt"
  }
  
  if (i == 2) {
    dataname <- "ift140ko-ud-15rm-32rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift140ko-201901-15rm-32rm-test-1"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-ud-15rm-32rm.txt"
  }
  
#   if (i == 3) {
#     dataname <- "dyn2h1ko-ift140ko-ud-15rm-4rm"
#     outdirpath <- "/media/josephus/Elements/DEanalysis/ift140ko-dft140ko-201901-15rm-4rm-test"
#     metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dft140ko-ud-15rm-4rm.txt"
#   }

  
  if (i == 4) {
    dataname <- "dyn2h1ko-15d-45rm-45rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-15d-201901-45rm-45rm-test-1"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-15d-45rm-45rm.txt"
  }

  if (i == 5) {
    dataname <- "ift140ko-15d-54rm-45rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift140ko-15d-201901-54rm-45rm-test-1"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-15d-54rm-45rm.txt"
  }
  
  
  if (i == 6) {
    dataname <- "dyn2h1ko-ift140ko-15d-54rm-4rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-ift140ko-15d-201901-54rm-4rm-test-1"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dyn2h1ko-15d-54rm-4rm.txt"
  }

  if (i == 7) {
    dataname <- "dyn2h1ko-wt-15d-45rm-345rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn2h1ko-wt-15d-201901-45rm-345rm-test-2"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-15d-45rm-345rm.txt"
  }

  if (i == 8) {
    dataname <- "ift140ko-wt-15d-54rm-345rm"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift140ko-wt-15d-201901-54rm-345rm-test-2"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-15d-54rm-345rm.txt"
  }


  
  ####################################################################
  # get meta information
  ####################################################################
  # this section reads in the meta file and makes it available
  # no need to change anything here
  
  meta.df <- read.table(metapath, sep = '\t',
                        header = TRUE,
                        stringsAsFactors = FALSE)
  meta.df$condition <- factor(meta.df$condition, # ensures preferred order
                              levels = unique(meta.df$condition))
  
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')
  
  ####################################################################
  # inferred paths
  ####################################################################
  # prepares subfolders - no changes to make
  
  outdirpath <- file.path(outdirpath, core.name(meta.df))
  meta.outdir <- file.path(outdirpath, "meta")
  cnts.outdir <- file.path(outdirpath, "count-table")
  DE.outdir <- file.path(outdirpath, "DE-table")
  plot.outdir <- file.path(outdirpath, "glimma-plots")
  hm.outdir <- file.path(outdirpath, "heatmap")
  hm.final.outdir <- file.path(hm.outdir, "final")
  hm.up.outdir <- file.path(hm.outdir, "up")
  hm.down.outdir <- file.path(hm.outdir, "down")
  hm.all.outdir <- file.path(hm.outdir, "all")
  hm.up.list.outdir <- file.path(hm.outdir, "up-list")
  hm.down.list.outdir <- file.path(hm.outdir, "down-list")
  hm.all.list.outdir <- file.path(hm.outdir, "all-list")
  GO.outdir <- file.path(outdirpath, "GO")
  pathway.outdir <- file.path(outdirpath, "pathway")
  goi.outdir <- file.path(outdirpath, "goi")
  lgenes.outdir <- file.path(outdirpath, "listgenes")
  gskb.outdir <- file.path(outdirpath, "gskb")
  PCA.outdirpath <- file.path(outdirpath, "pca")
  
  #####################################################################
  # ensure existence of output paths
  #####################################################################
  # creates  subfolders - no changes to make
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.all.list.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=lgenes.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=gskb.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=PCA.outdirpath, recursive = TRUE, showWarnings = FALSE)
  
  # https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
  #######################################################################
  # Create meta table
  #######################################################################
  # creates the meta file in subfolder 'meta' - no changes to make
  
  write.table(meta.df,
              file = file.path(meta.outdir, basename(metapath)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  #######################################################################
  # Create DESeq2.obj and DESeq2.obj.disp
  #######################################################################
  # does the DESeq2 analysis - no changes to make
  
  DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
  DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)
  
  #######################################################################
  # Create the count tables
  #######################################################################
  # prints out count tables - no changes to make
  
  cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = FALSE, averaged = FALSE,
                        sheetName = "raw.all")
  cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = TRUE, averaged = FALSE,
                        sheetName = "normalized.all")
  cnts.avg.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                            dataname = dataname,
                            printp = TRUE, normalized = TRUE, averaged = TRUE,
                            sheetName = "avg.normalized.all")
  
  #######################################################################
  # Create DE table
  #######################################################################
  # prints out DE tables - no changes to make
  
  res <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                    dataname = dataname,
                    printp = TRUE)
  resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                       dataname = dataname,
                       printp = TRUE,
                       filterp = TRUE, alpha = alpha, lFC = lFC)
  print.cnts.DE.sortings(cnts.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create plots
  #######################################################################
  # creates interactive plots in folder "glimma" - no changes to make
  
  meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
                dataname, top=300, launch = FALSE)
  
  meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
                dataname, alpha = alpha, lFC = lFC,
                launch = FALSE)
  

  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)
}


































##################################
# full runs
##################################

source("/home/josephus/Dropbox/a_myRNAseq/rdownstream/190220.rdownstream.central.R")

####################################################################
# DE settings
####################################################################
# set here tresholds
# important are lFC (log2FoldChange for DE analysis 1 or 2.5 we use)

alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 1       # log2FC limit
GO.alpha <- 0.01
GO.FDR   <- 0.05
GO.N     <- 30
GO.cl.N  <- 10

for (i in 4:19) {
  
  
  
  #########################################
  # set variables here which are common to all runs
  #########################################
  
  indirpath <- ""
  
  #########################################
  # define variables specific for each run
  #########################################
  
  if (i == 1) {
    dataname <- "dyn_ud_1234_wt_ud_12345"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn_ud_1234-vs-wt_ud_12345"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-ud.txt"
  }
  
  if (i == 2) {
    dataname <- "ift_ud_12345_wt_ud_12345"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_ud_12345-vs-wt_ud_12345"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-ud.txt"
  }
  
  if (i == 3) {
    dataname <- "ift_ud_12345_dyn_ud_1234"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_ud_12345-vs-dyn_ud_1234"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dft140ko-ud.txt" # right
  }
  
  if (i == 4) {
    dataname <- "dyn_15d_12345_wt_15d_12345"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn_15d_12345-vs-wt_15d_12345"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-15d.txt" # right
  }

  if (i == 5) {
    dataname <- "ift_15d_12345_wt_15d_12345"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_15d_12345-vs-wt_15d_12345"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-15d.txt" # right
  }
  
  # changed # to be!
  if (i == 6) {
    dataname <- "ift_15d_12345_dyn_15d_12345"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_15d_12345-vs-dyn_15d_12345"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dyn2h1ko-15d.txt" # nano
  }
  

  if (i == 7) {
    dataname <- "dyn_ud_123_wt_ud_1245"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn_ud_123-wt_ud_1245"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-ud-4rm-3rm.txt" # right
  }
  
  if (i == 8) {
    dataname <- "ift_ud_234_wt_ud_1245"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_ud_234-wt_ud_1245"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-ud-15rm-3rm.txt" # right
  }
  
  if (i == 9) {
    dataname <- "ift_ud_234_dyn_ud_123"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_ud_234-vs-dyn_ud_123"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dyn2h1ko-ud-15rm-4rm.txt"
  }
  
  if (i == 10) {
    dataname <- "dyn_15d_123_wt_15d_345"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn_15d_123-vs-wt_15d_345"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-15d-45rm-12rm.txt" # right
  }

  if (i == 11) {
    dataname <- "ift_15d_1234_wt_15d_345"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_15d_1234-vs-wt_15d_345"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-15d-5rm-12rm.txt" # right
  }
  
  if (i == 12) {
    dataname <- "ift_15d_1234_dyn_15d_123"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_15d_1234-vs-dyn_15d_123"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-dyn2h1ko-15d-5rm-45rm.txt" # right

  } # redo

  if (i == 13) {
    dataname <- "dyn_ud_123_wt_ud_145"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn_ud_123-vs-wt_ud_145"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-ud-4rm-32rm.txt" # right
  }
  
  if (i == 14) {
    dataname <- "ift_ud_234_wt_ud_145"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_ud_234-vs-wt_ud_145"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-ud-15rm-32rm.txt" # right
  }
  

  if (i == 15) {
    dataname <- "dyn_15d_123_wt_15d_123"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn_15d_123-vs-wt_15d_123"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-15d-45rm-45rm.txt" # right
  }

  if (i == 16) {
    dataname <- "ift_15d_123_wt_15d_123"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_15d_123-vs-wt_15d_123"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-15d-54rm-45rm.txt" # right
  }
  
  
  if (i == 17) {
    dataname <- "ift_15d_123_dyn_15d_123"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_15d_123-vs-dyn_15d_123"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-ift140ko-15d-54rm-4rm.txt"
  }

  if (i == 18) {
    dataname <- "dyn_15d_123_wt_15d_12"
    outdirpath <- "/media/josephus/Elements/DEanalysis/dyn_15d_123-vs-wt_15d_12"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-dyn2h1ko-wt-15d-45rm-345rm.txt" # right
  }

  if (i == 19) {
    dataname <- "ift_15d_123_wt_15d_12"
    outdirpath <- "/media/josephus/Elements/DEanalysis/ift_15d_123-vs-wt_15d_12"
    metapath   <- "/media/josephus/Elements/metas/meta-foxd2-201901-ift140ko-wt-15d-54rm-345rm.txt" # right
  }


  ####################################################################
  # get meta information
  ####################################################################
  # this section reads in the meta file and makes it available
  # no need to change anything here
  
  meta.df <- read.table(metapath, sep = '\t',
                        header = TRUE,
                        stringsAsFactors = FALSE)
  meta.df$condition <- factor(meta.df$condition, # ensures preferred order
                              levels = unique(meta.df$condition))
  
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')
  
  ####################################################################
  # inferred paths
  ####################################################################
  # prepares subfolders - no changes to make
  
  outdirpath <- file.path(outdirpath, core.name(meta.df))
  meta.outdir <- file.path(outdirpath, "meta")
  cnts.outdir <- file.path(outdirpath, "count-table")
  DE.outdir <- file.path(outdirpath, "DE-table")
  plot.outdir <- file.path(outdirpath, "glimma-plots")
  hm.outdir <- file.path(outdirpath, "heatmap")
  hm.final.outdir <- file.path(hm.outdir, "final")
  hm.up.outdir <- file.path(hm.outdir, "up")
  hm.down.outdir <- file.path(hm.outdir, "down")
  hm.all.outdir <- file.path(hm.outdir, "all")
  hm.up.list.outdir <- file.path(hm.outdir, "up-list")
  hm.down.list.outdir <- file.path(hm.outdir, "down-list")
  hm.all.list.outdir <- file.path(hm.outdir, "all-list")
  GO.outdir <- file.path(outdirpath, "GO")
  pathway.outdir <- file.path(outdirpath, "pathway")
  goi.outdir <- file.path(outdirpath, "goi")
  lgenes.outdir <- file.path(outdirpath, "listgenes")
  gskb.outdir <- file.path(outdirpath, "gskb")
  PCA.outdirpath <- file.path(outdirpath, "pca")
  
  #####################################################################
  # ensure existence of output paths
  #####################################################################
  # creates  subfolders - no changes to make
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.all.list.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=lgenes.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=gskb.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=PCA.outdirpath, recursive = TRUE, showWarnings = FALSE)
  
  # https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
  #######################################################################
  # Create meta table
  #######################################################################
  # creates the meta file in subfolder 'meta' - no changes to make
  
  write.table(meta.df,
              file = file.path(meta.outdir, basename(metapath)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  #######################################################################
  # Create DESeq2.obj and DESeq2.obj.disp
  #######################################################################
  # does the DESeq2 analysis - no changes to make
  
  DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
  DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)
  
  #######################################################################
  # Create the count tables
  #######################################################################
  # prints out count tables - no changes to make
  
  cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = FALSE, averaged = FALSE,
                        sheetName = "raw.all")
  cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = TRUE, averaged = FALSE,
                        sheetName = "normalized.all")
  cnts.avg.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                            dataname = dataname,
                            printp = TRUE, normalized = TRUE, averaged = TRUE,
                            sheetName = "avg.normalized.all")

  #######################################################################
  # Create DE table
  #######################################################################
  # prints out DE tables - no changes to make
  
  res <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                    dataname = dataname,
                    printp = TRUE)
  resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                       dataname = dataname,
                       printp = TRUE,
                       filterp = TRUE, alpha = alpha, lFC = lFC)
  print.cnts.DE.sortings(cnts.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create plots
  #######################################################################
  # creates interactive plots in folder "glimma" - no changes to make
  
  meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
                dataname, top=300, launch = FALSE)
  
  meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
                dataname, alpha = alpha, lFC = lFC,
                launch = FALSE)
  
  #######################################################################
  # Create heatmap
  #######################################################################
  # creates heatmap with scaled values in 'heatmap' folder
  # you can change 'k =' to what you want. e.g. 12, 24, etc.
  # set the seed for reproducibility
  
  set.seed(123) # for reproducibility of heatmaps set seeds to a value
  try(scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm$`nrm-counts-avg`, resSig,
                                              outdirpath = hm.all.outdir, selected.genes = NULL, name.add = "all",
                                              dataname = dataname,
                                              k = 8, printp=TRUE,
                                              alpha=alpha, lFC=lFC, filterp=TRUE,
                                              xlsxp=TRUE, csvp=FALSE, tsvp=FALSE))
  
  
  #######################################################################
  # Create GO
  #######################################################################
  # GO analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # inDir is outdirpath
  # outBase is GO.outdir
  last <- function(l) l[length(l)]
  DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))

  
  try(do_all_GO_enrichment(DEfpath, GO.outdir))
  
  #######################################################################
  # KEGG
  #######################################################################
  # KEGG analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # DEfpath is like above
  # outBase <- pathway.outdir
  try(do_KEGG_enrichment(DEfpath, pathway.outdir))

  #######################################################################
  # gskb
  #######################################################################
  # gskb analysis for mouse - actually no changes needed
  # but if you don't want it (quicker analysis) - comment this section out
  
  # DEfpath is like above
  # outBase <- gskb.outdir
  try(do_all_gskb_enrichments(DEfpath,
                          gskb.outdir,
                          all_ez_grouped))
  
  #######################################################################
  # nicer MDS/PCA plots using plotly
  #######################################################################
  # creates nicer non-interactive PCA plots (for publication)
  # - no changes needed
  # - but if PCA plot not important - comment this section out
  
  # request DE folder for significant genes path in previous results
  last <- function(l) l[length(l)]
  cntsDEfpath <- last(dir(DE.outdir, pattern = "DE-cnts-sig", full.names = TRUE))
  
  try(meta2plyMDSplot(meta.df, 
                  DESeq2.obj.disp, 
                  PCA.outdirpath,
                  dataname, 
                  top=500, 
                  launch = FALSE))
  
  try(meta2plyMDSplotDE(cntsDEfpath = cntsDEfpath, 
                    outdirpath = PCA.outdirpath,
                    dataname = dataname, 
                    corename = core.name(meta.df), 
                    top=500, 
                    launch = FALSE))
  
  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)

}



################################################################################

###########################################
# foxd2 f7-vs-wt correction 20200403
###########################################

# nano /media/josephus/Elements/metas/meta-foxd2-20200403-mk4_wt_foxd2-mk4_f7.txt

sampleName	fileName	condition	testing
foxd2_gw320_1	/media/josephus/archive_big/count/foxd2_201907p/external-MN-26-mK4-Gw320-foxd2-rep1-16863-sym-fcount.tab	wt_foxd2	denom
foxd2_gw320_2	/media/josephus/archive_big/count/foxd2_201907p/external-MN-27-mK4-Gw320-foxd2-rep2-16862-sym-fcount.tab	wt_foxd2	denom
foxd2_gw320_3	/media/josephus/archive_big/count/foxd2_201907p/external-MN-28-mK4-Gw320-foxd2-rep3-16861-sym-fcount.tab	wt_foxd2	denom
foxd2_gw320_4	/media/josephus/archive_big/count/foxd2_201907p/external-MN-29-mK4-Gw320-foxd2-rep4-16860-sym-fcount.tab	wt_foxd2	denom
foxd2_gw320_5	/media/josephus/archive_big/count/foxd2_201907p/external-MN-30-mK4-Gw320-foxd2-rep5-16859-sym-fcount.tab	wt_foxd2	denom
F7_mk4_1	/media/josephus/archive_big/count/foxd2_201907p/external-MN-16-mK4-F7-rep1-16873-sym-fcount.tab	f7	num
F7_mk4_2	/media/josephus/archive_big/count/foxd2_201907p/external-MN-17-mK4-F7-rep2-16872-sym-fcount.tab	f7	num
F7_mk4_3	/media/josephus/archive_big/count/foxd2_201907p/external-MN-18-mK4-F7-rep3-16871-sym-fcount.tab	f7	num
F7_mk4_4	/media/josephus/archive_big/count/foxd2_201907p/external-MN-19-mK4-F7-rep4-16870-sym-fcount.tab	f7	num
F7_mk4_5	/media/josephus/archive_big/count/foxd2_201907p/external-MN-20-mK4-F7-rep5-16869-sym-fcount.tab	f7	num





####################################################################
# source 
####################################################################
# this loads the functions necessary

computer <- "ss"

if (computer == "work") {
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

computer <- "ss"

if (computer == "work") {
  script_dir <- "/media/daten/arnold/josephus/script"
  out_dir    <- "/media/daten/arnold/josephus/results/DEanalaysis"
  meta_dir   <- "/media/daten/arnold/josephus/metas"
} else {
  script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
  out_dir    <- "/media/josephus/Elements/DEanalysis"
  meta_dir   <- "/media/josephus/Elements/metas"
}
dir.create(out_dir, recursive=TRUE)


data_names  <-  c("foxd2-foxd2-F7-corr")
outdir_names <- c("foxd2-foxd2-F7-20200403")
meta_fnames <- c("meta-foxd2-20200403-mk4_wt_foxd2-mk4_f7.txt")

####################################################################
# DE settings
####################################################################
# set here tresholds
# important are lFC (log2FoldChange for DE analysis 1 or 2.5 we use)

alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 1       # log2FC limit
GO.alpha <- 0.01
GO.FDR   <- 0.05
GO.N     <- 30
GO.cl.N  <- 10

# ####################################################################
# # Use functions for analysis
# ####################################################################
# this big for loop handles the runs
# I change 1:2 to the n needed
# like 1:7 etc
# if you want to jump over other - just give the sequence
# for (i in c(1, 4, 7)) { ... }

########################
# gskb install
########################

# wget https://bioconductor.org/packages/release/data/experiment/src/contrib/gskb_1.16.0.tar.gz
# install.packages("~/R_sources/gskb_1.16.0.tar.gz", type="source", dependencies=TRUE)


if (!require(foreach)) {
  install.packages("foreach")
  require(foreach)
}

if (!require(doParallel)) {
  install.packages("doParallel")
  require(doParallel)
}

numCores <- 1
registerDoParallel(numCores)

# foreach (i = 1:1) %dopar% {
i <- 1
  #########################################
  # set variables here which are common to all runs
  #########################################
  
  indirpath <- ""
  
  #########################################
  # define variables specific for each run
  #########################################
  
  dataname <- data_names[i]
  outdirpath <- file.path(out_dir, outdir_names[i])
  dir.create(outdirpath, recursive=TRUE)
  metapath <-   file.path(meta_dir, meta_fnames[i])
  
  ####################################################################
  # get meta information
  ####################################################################
  # this section reads in the meta file and makes it available
  # no need to change anything here
  
  meta.df <- read.table(metapath, sep = '\t',
                        header = TRUE,
                        stringsAsFactors = FALSE)
  meta.df$condition <- factor(meta.df$condition, # ensures preferred order
                              levels = unique(meta.df$condition))
  
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')
  
  ####################################################################
  # inferred paths
  ####################################################################
  # prepares subfolders - no changes to make
  
  outdirpath <- file.path(outdirpath, core.name(meta.df))
  meta.outdir <- file.path(outdirpath, "meta")
  cnts.outdir <- file.path(outdirpath, "count-table")
  DE.outdir <- file.path(outdirpath, "DE-table")
  plot.outdir <- file.path(outdirpath, "glimma-plots")
  hm.outdir <- file.path(outdirpath, "heatmap")
  hm.final.outdir <- file.path(hm.outdir, "final")
  hm.up.outdir <- file.path(hm.outdir, "up")
  hm.down.outdir <- file.path(hm.outdir, "down")
  hm.all.outdir <- file.path(hm.outdir, "all")
  hm.up.list.outdir <- file.path(hm.outdir, "up-list")
  hm.down.list.outdir <- file.path(hm.outdir, "down-list")
  hm.all.list.outdir <- file.path(hm.outdir, "all-list")
  GO.outdir <- file.path(outdirpath, "GO")
  pathway.outdir <- file.path(outdirpath, "pathway")
  goi.outdir <- file.path(outdirpath, "goi")
  lgenes.outdir <- file.path(outdirpath, "listgenes")
  gskb.outdir <- file.path(outdirpath, "gskb")
  PCA.outdirpath <- file.path(outdirpath, "pca")
  
  #####################################################################
  # ensure existence of output paths
  #####################################################################
  # creates  subfolders - no changes to make
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.all.list.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=lgenes.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=gskb.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=PCA.outdirpath, recursive = TRUE, showWarnings = FALSE)
  
  # https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
  #######################################################################
  # Create meta table
  #######################################################################
  # creates the meta file in subfolder 'meta' - no changes to make
  
  write.table(meta.df,
              file = file.path(meta.outdir, basename(metapath)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  #######################################################################
  # Create DESeq2.obj and DESeq2.obj.disp
  #######################################################################
  # does the DESeq2 analysis - no changes to make
  
  DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
  DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)
  
  #######################################################################
  # Create the count tables
  #######################################################################
  # prints out count tables - no changes to make
  
  cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = FALSE, averaged = FALSE,
                        sheetName = "raw.all")
  cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = TRUE, averaged = FALSE,
                        sheetName = "normalized.all")
  cnts.avg.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                            dataname = dataname,
                            printp = TRUE, normalized = TRUE, averaged = TRUE,
                            sheetName = "avg.normalized.all")

  #######################################################################
  # Create DE table
  #######################################################################
  # prints out DE tables - no changes to make
  
  res <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                    dataname = dataname,
                    printp = TRUE)
  resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                       dataname = dataname,
                       printp = TRUE,
                       filterp = TRUE, alpha = alpha, lFC = lFC)
  print.cnts.DE.sortings(cnts.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create plots
  #######################################################################
  # creates interactive plots in folder "glimma" - no changes to make
  
  meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
                dataname, top=300, launch = FALSE)
  
  meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
                dataname, alpha = alpha, lFC = lFC,
                launch = FALSE)
  
  #######################################################################
  # Create heatmap
  #######################################################################
  # creates heatmap with scaled values in 'heatmap' folder
  # you can change 'k =' to what you want. e.g. 12, 24, etc.
  # set the seed for reproducibility
  
  set.seed(123) # for reproducibility of heatmaps set seeds to a value
  try(scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm$`nrm-counts-avg`, resSig,
                                              outdirpath = hm.all.outdir, selected.genes = NULL, name.add = "all",
                                              dataname = dataname,
                                              k = 8, printp=TRUE,
                                              alpha=alpha, lFC=lFC, filterp=TRUE,
                                              xlsxp=TRUE, csvp=FALSE, tsvp=FALSE))
  
  
  #######################################################################
  # Create GO
  #######################################################################
  # GO analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # inDir is outdirpath
  # outBase is GO.outdir
  last <- function(l) l[length(l)]
  DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))

  
  try(do_all_GO_enrichment(DEfpath, GO.outdir))
  
  #######################################################################
  # KEGG
  #######################################################################
  # KEGG analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # DEfpath is like above
  # outBase <- pathway.outdir
  try(do_KEGG_enrichment(DEfpath, pathway.outdir))

  #######################################################################
  # gskb
  #######################################################################
  # gskb analysis for mouse - actually no changes needed
  # but if you don't want it (quicker analysis) - comment this section out
  
  # DEfpath is like above
  # outBase <- gskb.outdir
  try(do_all_gskb_enrichments(DEfpath,
                          gskb.outdir,
                          all_ez_grouped))
  
  #######################################################################
  # nicer MDS/PCA plots using plotly
  #######################################################################
  # creates nicer non-interactive PCA plots (for publication)
  # - no changes needed
  # - but if PCA plot not important - comment this section out
  
  # request DE folder for significant genes path in previous results
  last <- function(l) l[length(l)]
  cntsDEfpath <- last(dir(DE.outdir, pattern = "DE-cnts-sig", full.names = TRUE))
  
  try(meta2plyMDSplot(meta.df, 
                  DESeq2.obj.disp, 
                  PCA.outdirpath,
                  dataname, 
                  top=500, 
                  launch = FALSE))
  
  try(meta2plyMDSplotDE(cntsDEfpath = cntsDEfpath, 
                    outdirpath = PCA.outdirpath,
                    dataname = dataname, 
                    corename = core.name(meta.df), 
                    top=500, 
                    launch = FALSE))
  
  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)

# }


















