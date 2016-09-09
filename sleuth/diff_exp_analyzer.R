library("sleuth")
library("optparse")
library("files")

#command line arguments
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-d", "--directory"), type='character', default=character(0),
              help="Please specify the directory"),
  make_option(c("-ge", "--genovar"), type='character', default=character(0),
              help="Please specify the genotype variable name"),
  make_option(c("-ba", "--batch"), action="store_true", default=FALSE,
              help="Batch correction method"),
  make_option(c("-bv", "--batchvar"), type='character', default=character(0),
              help="Please specify the batch variable name"),
  make_option(c("-s", "--shiny"), action='store_true', default=FALSE,
              help="Command to open shiny console")
)

opt = parse_args(OptionParser(option_list=option_list))

try (if(length(opt$d) == 0) stop('Directory cannot be empty'))

try (if(!file.exists(opt$d)) stop('Directory must exist'))
# try (if(file.exists(opt$d)) setwd(opt$d))

if(length(opt$g) == 0){
  genovar = 'genotypezmt'
} else {
  genovar = paste('genotype', opt$ge, sep='')
}

if(length(opt$batchvar) == 0){
  batchvar = 'batchb'
} else {
  batchvar = paste('batch', opt$batchvar, sep='')
}




#gene info for sleuth
print("Fetching bioMart info:")
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
print('#renaming genes:')
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


#point to your directory
base_dir <- opt$d

#get ids
sample_id <- dir(file.path(base_dir, "results"))
print(sample_id)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
print(kal_dirs)
s2c <- read.table(file.path(base_dir, "rna_seq_info.txt"), header = TRUE, stringsAsFactors= FALSE)
print(s2c)

if (opt$batch == TRUE){
  s2c <- dplyr::select(s2c, sample= experiment, genotype, batch)
  s2c <- dplyr::mutate(s2c, path = kal_dirs)
  so <- sleuth_prep(s2c, ~ genotype+batch, target_mapping= t2g)
  so <- sleuth_fit(so,~ genotype+batch, fit_name = 'full')
  so <- sleuth_fit(so, ~batch, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
} else {
  s2c <- dplyr::select(s2c, sample = experiment, genotype)
  s2c <- dplyr::mutate(s2c, path = kal_dirs)
  so <- sleuth_prep(s2c, ~ genotype, target_mapping= t2g)
  so <- sleuth_fit(so,~ genotype, fit_name = 'full')
}

#print(s2c)

#prepend and make object, state maximum model here
#so <- sleuth_prep(s2c, ~ genotype, target_mapping= t2g)

#fit the model
#so <- sleuth_fit(so,~ genotype, fit_name = 'full')

#Wald test implementations
#no interactions
so <- sleuth_wt(so, which_beta = genovar, which_model = 'full')


#write results to tables
results_table <- sleuth_results(so, genovar, 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betas.csv', sep='/'))

results_table <- sleuth_results(so, genovar, 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betas.csv', sep='/'))

if (opt$batch == TRUE) {
  so <- sleuth_wt(so, which_beta = '(Intercept)', which_model = 'full')
  so <- sleuth_wt(so, which_beta = batchvar, which_model = 'full')
  results_table <- sleuth_results(so, '(Intercept)','full', test_type= 'wt')
  write.csv(results_table, paste(base_dir, 'intercept.csv', sep='/'))
  results_table <- sleuth_results(so, batchvar,'full', test_type= 'wt')
  write.csv(results_table, paste(base_dir, 'batch.csv', sep='/'))
  sr <- sleuth_results(so, 'reduced:full', 'lrt')
  write.csv(sr, paste(base_dir, 'batch_lrt.csv', sep='/'))
}


#if you want to look at shiny
if (opt$shiny){
  sleuth_live(so)
}
