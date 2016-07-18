#!/usr/bin/Rscript
#Matti Ruuskanen, Jul 2016
#Lopez R, Cowley A, Li W, McWilliam H.(2013) Using EMBL-EBI Services via Web Interface and Programmatically via Web Services. Curr Protoc Bioinformatics Volume 48 (2014) p.3.12.1-3.12.50. DOI: 10.1002/0471250953.bi0312s48 
#McWilliam H., Li W., Uludag M., Squizzato S., Park Y.M., Buso N., Cowley A.P., Lopez R.(2013) Analysis Tool Web Services from the EMBL-EBI Nucleic Acids Research 41: W597-W600. PubMed Id: 23671338 Abstract DOI: 10.1093/nar/gkt376
#This program eats a list of accession numbers as an input and fetches their sequences from the database using DbFetch
#Output is a tab limited .csv file with metadata and protein sequences, and a .fa file of the CDS (also possible to get up and downstream sequences in a specified range)

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("optparse", "parallel", "doParallel", "foreach")
ipak(packages)

ncores <- detectCores()

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL, 
              help = "dataset file name", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = getwd(),
              help = "output folder [default = %default]", metavar = "character"),
  make_option(c("-p", "--processors"), type = "numeric", default = ncores,
              help = "number of cores to be used [default = all cores]", metavar = "numeric"),
  make_option(c("-u", "--upstream"), type = "numeric", default = 0,
              help = "extract sequence n nucleotides upstream [default = %default]", metavar = "numeric"),
  make_option(c("-d", "--downstream"), type = "numeric", default = 0,
              help = "extract sequence n nucleotides downstream [default = %default]", metavar = "numeric")
); 

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call. = FALSE)
}

if (opt$processors > ncores) {
  opt$processors <- ncores
  }
registerDoParallel(cores <- opt$processors)

setwd(opt$out)

if (opt$upstream < 0) {
  opt$upstream <- 0
  }
if (opt$downstream < 0) {
  opt$downstream <- 0
  }

#read in the file
df <- read.table(opt$file, header = F)

sequences <- list()
list <- list()
failed <- list()
dirs<-c("sequence","list","failed")
for(i in dirs){
   if(!dir.exists(i)){
     dir.create(i)
   }  
}

targetseq <- df[,1]
targetseqfiles <- paste("sequence_", targetseq, ".csv", sep = "")
todoseq <- targetseq[!targetseqfiles %in% dir(dirs[1])]
parresult <- foreach (i = unique(todoseq), .combine = rbind) %dopar% {
  source("extractAcc.R")
  out <- extractAcc(i, opt$upstream, opt$downstream)
  if (out == "obsolete") {
    write.table(unlist(i), paste("failed/failed_", i, ".txt", sep = ""), row.names = F, quote = F, col.names = F)
  } else {
  write.table(unlist(out[[1]]), paste("sequence/sequence_", i, ".txt", sep = ""), row.names = F, quote = F, col.names = F)
  write.csv(out[[2]],paste("list/list_",i,".csv",sep = ""),row.names = F)
  }
}

j <- 0
for( i in dir("sequence")){
  j <- j+1
  sequences[[j]] <- readLines(paste("sequence/",i,sep = ""))
}
sequences <- unlist(sequences)
sequences <- sequences[complete.cases(sequences)]
write.table(sequences,file = paste(opt$out, "/out.fa", sep = ""), row.names = F, quote = F, col.names = F)


j <- 0
for( i in dir("list")){
  j <- j+1
  list[[j]] <- read.csv(paste("list/", i, sep = ""))
}
list <- do.call("rbind", list)
list <- list[complete.cases(list$EMBL.id),]
write.csv(list, file = paste(opt$out, "/out.csv", sep = ""), row.names = F)

if (!is.null(failed)) {
  j <- 0
  for( i in dir("failed")){
    j <- j+1
    failed[[j]] <- readLines(paste("failed/", i, sep = ""))
  }
  failed <- unlist(failed)
  write.table(failed, file = paste(opt$out, "/out_failed.txt", sep = ""), row.names = F, quote = F, col.names = F)
}
unlink(c("sequence", "list", "failed"), recursive = T)
