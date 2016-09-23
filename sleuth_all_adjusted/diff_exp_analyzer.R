library("sleuth")
library("optparse")
library("files")

#gene info for sleuth
print("Fetching bioMart info:")
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
print('#renaming genes:')
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


#point to your directory
base_dir <- 'kallisto'

#get ids
sample_id <- dir(file.path(base_dir, "results"))
print(sample_id)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
print(kal_dirs)
s2c <- read.table(file.path(base_dir, "rna_seq_info.txt"), header = TRUE, stringsAsFactors= FALSE)
print(s2c)

s2c <- dplyr::select(s2c, sample= experiment, genotype, batch)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
so <- sleuth_prep(s2c, ~ genotype+batch, target_mapping= t2g)
so <- sleuth_fit(so,~ genotype+batch, fit_name = 'full')
so <- sleuth_fit(so, ~batch, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

#print(s2c)

#prepend and make object, state maximum model here
#so <- sleuth_prep(s2c, ~ genotype, target_mapping= t2g)

#fit the model
#so <- sleuth_fit(so,~ genotype, fit_name = 'full')

#Wald test implementations
#no interactions
so <- sleuth_wt(so, which_beta = 'genotypebB', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'genotypebA', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'genotypebC', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'genotypebD', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'genotypebE', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'genotypebF', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'genotypebG', which_model = 'full')


#write results to tables
results_table <- sleuth_results(so, 'genotypebB', 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betasB.csv', sep='/'))
results_table <- sleuth_results(so, 'genotypebA', 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betasA.csv', sep='/'))
results_table <- sleuth_results(so, 'genotypebC', 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betasC.csv', sep='/'))
results_table <- sleuth_results(so, 'genotypebD', 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betasD.csv', sep='/'))
results_table <- sleuth_results(so, 'genotypebE', 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betasE.csv', sep='/'))
results_table <- sleuth_results(so, 'genotypebF', 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betasF.csv', sep='/'))
results_table <- sleuth_results(so, 'genotypebG', 'full', test_type= 'wt')
write.csv(results_table, paste(base_dir, 'betasG.csv', sep='/'))

#so <- sleuth_wt(so, which_beta = '(Intercept)', which_model = 'full')
#so <- sleuth_wt(so, which_beta = batchvar, which_model = 'full')

#results_table <- sleuth_results(so, '(Intercept)','full', test_type= 'wt')
#write.csv(results_table, paste(base_dir, 'intercept.csv', sep='/'))

#results_table <- sleuth_results(so, batchvar,'full', test_type= 'wt')
#write.csv(results_table, paste(base_dir, 'batch.csv', sep='/'))

sr <- sleuth_results(so, 'reduced:full', 'lrt')
write.csv(sr, paste(base_dir, 'batch_lrt.csv', sep='/'))



sleuth_live(so)

