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

#gene info for sleuth
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
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
s2c <- dplyr::select(s2c, sample = experiment, genotype)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

#prepend and make object, state maximum model here
so <- sleuth_prep(s2c, ~ genotype, target_mapping= t2g)

#fit the model
so <- sleuth_fit(so,~ genotype, fit_name = 'full')

#Wald test implementations
#no interactions
so <- sleuth_wt(so, which_beta = genovar, which_model = 'full')

#if you want to look at shiny
if (opt$shiny){
  sleuth_live(so)
}

#write results to tables
results_table <- sleuth_results(so, genovar,'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betas.csv', sep='/'))