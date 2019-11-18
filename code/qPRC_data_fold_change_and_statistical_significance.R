# Rscript qPRC_data_fold_change_and_statistical_significance.R <input_file> > <output_file>
# Rscript qPRC_data_fold_change_and_statistical_significance.R example_data.txt > example_output.txt

# Please note that the following are currently hard-coded and need to be converted into input parameters so as to generalise this script:
# treatment vs control cell-lines
# target genes vs control gene

library(ggplot2)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
library(plyr)
options(width=200)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'example_data.txt' ) # for testing

args_infile = args[1] # example_data.txt

set.seed(1)

qpcr = read.table(infile, sep = "\t", header=TRUE)

# 41 cycles is considered to show that target is not present.
qpcr$threshold_cycle = ifelse( (qpcr$threshold_cycle=="Undetermined"), 41, qpcr$threshold_cycle )
qpcr$Ct_mean = ifelse( (qpcr$Ct_mean=="Undetermined"), 41, qpcr$Ct_mean )
qpcr$Ct_SD = ifelse( is.na(qpcr$Ct_SD), 0, qpcr$Ct_SD )
qpcr$threshold_cycle = as.numeric(qpcr$threshold_cycle)
qpcr$Ct_mean = as.numeric(qpcr$Ct_mean)
qpcr$Ct_SD = as.numeric(qpcr$Ct_SD)

cell_line = c( "H28-BR0.1", "H28-BR0.3", "MeT5A-BR0.1", "MeT5A-BR0.3", "REN-BR0.1", "REN-BR0.3", "VMC23-BR0.3" )
control = c( "H28", "H28", "MeT5A", "MeT5A", "REN", "REN", "VMC23" )
target_gene = c( "ATP5G1V2", "ATP5G2", "ATP5G2V3", "ATP5G3", "ATP5B", "ATP5IF1", "RNA18S1" )
reference_gene = rep.int( c("ATP5G1"), length(target_gene) )

fold_change_target_gene = rep( target_gene, each=length(cell_line) )
fold_change_reference_gene = rep( reference_gene, each=length(cell_line) )
fold_change_cell_line = rep( cell_line, length(target_gene) )
fold_change_control = rep( control, length(target_gene) )

temp = cbind( fold_change_target_gene, fold_change_reference_gene )
fold_change = temp
temp = cbind( temp, fold_change_cell_line )
fold_change = temp
temp = cbind( temp, fold_change_control )
fold_change = temp
colnames(fold_change) = c( "target_gene", "reference_gene", "cell_line", "control" )
fold_change = as.data.frame(fold_change)
rm(temp)
fold_change_template = fold_change

list_ct_values = function( expression_gene, cell_line, all_values )
{
	all_values[ ((all_values$gene==expression_gene) & (all_values$cell_line==cell_line)), c("threshold_cycle") ]
}

calculate_fold_change = function( params, all_values )
{
	target = params[1]
	ref = params[2]
	cell = params[3]
	control = params[4]
	ct_values_for_target_gene_in_sample = list_ct_values( target, cell, all_values=all_values )
	ct_values_for_reference_gene_in_sample = list_ct_values( ref, cell, all_values=all_values )
	ct_values_for_target_gene_in_control = list_ct_values( target, control, all_values=all_values )
	ct_values_for_reference_gene_in_control = list_ct_values( ref, control, all_values=all_values )
	mean_of_ct_values_for_target_gene_in_sample = mean(ct_values_for_target_gene_in_sample)
	mean_of_ct_values_for_reference_gene_in_sample = mean(ct_values_for_reference_gene_in_sample)
	mean_of_ct_values_for_target_gene_in_control = mean(ct_values_for_target_gene_in_control)
	mean_of_ct_values_for_reference_gene_in_control = mean(ct_values_for_reference_gene_in_control)
	mean_of_ct_values_for_target_gene_in_sample = ifelse( (mean_of_ct_values_for_target_gene_in_sample>=36), 41, mean_of_ct_values_for_target_gene_in_sample )
	mean_of_ct_values_for_reference_gene_in_sample = ifelse( (mean_of_ct_values_for_reference_gene_in_sample>=36), 41, mean_of_ct_values_for_reference_gene_in_sample )
	mean_of_ct_values_for_target_gene_in_control = ifelse( (mean_of_ct_values_for_target_gene_in_control>=36), 41, mean_of_ct_values_for_target_gene_in_control )
	mean_of_ct_values_for_reference_gene_in_control = ifelse( (mean_of_ct_values_for_reference_gene_in_control>=36), 41, mean_of_ct_values_for_reference_gene_in_control )
	delta_CQ_for_sample = mean_of_ct_values_for_target_gene_in_sample - mean_of_ct_values_for_reference_gene_in_sample
	delta_CQ_for_control = mean_of_ct_values_for_target_gene_in_control - mean_of_ct_values_for_reference_gene_in_control
	delta_CQ_for_sample = ifelse( (mean_of_ct_values_for_target_gene_in_sample<41), delta_CQ_for_sample, 41 )
	delta_CQ_for_control = ifelse( (mean_of_ct_values_for_target_gene_in_control<41), delta_CQ_for_control, 41 )
	delta_delta_CQ = delta_CQ_for_sample - delta_CQ_for_control
	delta_delta_CQ = ifelse( (delta_CQ_for_sample<41), delta_delta_CQ, 100 )
	fold_change = 2^(-1 * delta_delta_CQ)
	fold_change
}

fold_change$fold_change = apply( fold_change, 1, calculate_fold_change, all_values=qpcr )
fold_change$fold_change = sprintf( "%3.3f", fold_change$fold_change )
fold_change

# 	   target_gene reference_gene   cell_line control fold_change
# 	1     ATP5G1V2         ATP5G1   H28-BR0.1     H28       1.662
# 	2     ATP5G1V2         ATP5G1   H28-BR0.3     H28       0.204
# 	3     ATP5G1V2         ATP5G1 MeT5A-BR0.1   MeT5A       0.097
# 	4     ATP5G1V2         ATP5G1 MeT5A-BR0.3   MeT5A       0.111
# 	5     ATP5G1V2         ATP5G1   REN-BR0.1     REN       1.137
# 	6     ATP5G1V2         ATP5G1   REN-BR0.3     REN       0.094
# 	7     ATP5G1V2         ATP5G1 VMC23-BR0.3   VMC23       0.132
# 	8       ATP5G2         ATP5G1   H28-BR0.1     H28       0.000
# 	9       ATP5G2         ATP5G1   H28-BR0.3     H28       0.000
# 	10      ATP5G2         ATP5G1 MeT5A-BR0.1   MeT5A       0.000
# 	11      ATP5G2         ATP5G1 MeT5A-BR0.3   MeT5A       0.000
# 	12      ATP5G2         ATP5G1   REN-BR0.1     REN       1.600
# 	13      ATP5G2         ATP5G1   REN-BR0.3     REN       0.086
# 	14      ATP5G2         ATP5G1 VMC23-BR0.3   VMC23       0.282
# 	15    ATP5G2V3         ATP5G1   H28-BR0.1     H28       0.668
# 	16    ATP5G2V3         ATP5G1   H28-BR0.3     H28       0.000
# 	17    ATP5G2V3         ATP5G1 MeT5A-BR0.1   MeT5A       0.000
# 	18    ATP5G2V3         ATP5G1 MeT5A-BR0.3   MeT5A       0.000
# 	19    ATP5G2V3         ATP5G1   REN-BR0.1     REN       0.643
# 	20    ATP5G2V3         ATP5G1   REN-BR0.3     REN       0.481
# 	21    ATP5G2V3         ATP5G1 VMC23-BR0.3   VMC23       0.000
# 	22      ATP5G3         ATP5G1   H28-BR0.1     H28       8.605
# 	23      ATP5G3         ATP5G1   H28-BR0.3     H28       0.669
# 	24      ATP5G3         ATP5G1 MeT5A-BR0.1   MeT5A       0.688
# 	25      ATP5G3         ATP5G1 MeT5A-BR0.3   MeT5A       2.426
# 	26      ATP5G3         ATP5G1   REN-BR0.1     REN       0.804
# 	27      ATP5G3         ATP5G1   REN-BR0.3     REN       0.199
# 	28      ATP5G3         ATP5G1 VMC23-BR0.3   VMC23       0.347
# 	29       ATP5B         ATP5G1   H28-BR0.1     H28       0.000
# 	30       ATP5B         ATP5G1   H28-BR0.3     H28       0.000
# 	31       ATP5B         ATP5G1 MeT5A-BR0.1   MeT5A       0.000
# 	32       ATP5B         ATP5G1 MeT5A-BR0.3   MeT5A       0.000
# 	33       ATP5B         ATP5G1   REN-BR0.1     REN       0.869
# 	34       ATP5B         ATP5G1   REN-BR0.3     REN       0.000
# 	35       ATP5B         ATP5G1 VMC23-BR0.3   VMC23       0.321
# 	36     ATP5IF1         ATP5G1   H28-BR0.1     H28       0.895
# 	37     ATP5IF1         ATP5G1   H28-BR0.3     H28       0.294
# 	38     ATP5IF1         ATP5G1 MeT5A-BR0.1   MeT5A       0.007
# 	39     ATP5IF1         ATP5G1 MeT5A-BR0.3   MeT5A       0.001
# 	40     ATP5IF1         ATP5G1   REN-BR0.1     REN       0.821
# 	41     ATP5IF1         ATP5G1   REN-BR0.3     REN       0.152
# 	42     ATP5IF1         ATP5G1 VMC23-BR0.3   VMC23       0.453
# 	43     RNA18S1         ATP5G1   H28-BR0.1     H28       0.721
# 	44     RNA18S1         ATP5G1   H28-BR0.3     H28       0.224
# 	45     RNA18S1         ATP5G1 MeT5A-BR0.1   MeT5A       0.004
# 	46     RNA18S1         ATP5G1 MeT5A-BR0.3   MeT5A       0.000
# 	47     RNA18S1         ATP5G1   REN-BR0.1     REN       0.704
# 	48     RNA18S1         ATP5G1   REN-BR0.3     REN       0.321
# 	49     RNA18S1         ATP5G1 VMC23-BR0.3   VMC23       0.397

calculate_randomised_fold_change = function( params, all_values )
{
	target = params[1]
	ref = params[2]
	cell = params[3]
	control = params[4]
	ct_values_for_target_gene_in_sample = list_ct_values( target, cell, all_values=all_values )
	ct_values_for_reference_gene_in_sample = list_ct_values( ref, cell, all_values=all_values )
	ct_values_for_target_gene_in_control = list_ct_values( target, control, all_values=all_values )
	ct_values_for_reference_gene_in_control = list_ct_values( ref, control, all_values=all_values )
	num_values_for_target_gene_in_sample = length(ct_values_for_target_gene_in_sample)
	num_values_for_reference_gene_in_sample = length(ct_values_for_reference_gene_in_sample)
	num_values_for_target_gene_in_control = length(ct_values_for_target_gene_in_control)
	num_values_for_reference_gene_in_control = length(ct_values_for_reference_gene_in_control)
	randomised_ct_values_for_target_gene_in_sample_and_control = sample( c( ct_values_for_target_gene_in_sample, ct_values_for_target_gene_in_control ) )
	randomised_ct_values_for_reference_gene_in_sample_and_control = sample( c( ct_values_for_reference_gene_in_sample, ct_values_for_reference_gene_in_control ) )
	ct_values_for_target_gene_in_sample = randomised_ct_values_for_target_gene_in_sample_and_control[ 1 : num_values_for_target_gene_in_sample ]
	ct_values_for_target_gene_in_control = randomised_ct_values_for_target_gene_in_sample_and_control[ (num_values_for_target_gene_in_sample+1) : (num_values_for_target_gene_in_sample+num_values_for_target_gene_in_control) ]
	ct_values_for_reference_gene_in_sample = randomised_ct_values_for_reference_gene_in_sample_and_control[ 1 : num_values_for_reference_gene_in_sample ]
	ct_values_for_reference_gene_in_control = randomised_ct_values_for_reference_gene_in_sample_and_control[ (num_values_for_reference_gene_in_sample+1) : (num_values_for_reference_gene_in_sample+num_values_for_reference_gene_in_control) ]
	mean_of_ct_values_for_target_gene_in_sample = mean(ct_values_for_target_gene_in_sample)
	mean_of_ct_values_for_reference_gene_in_sample = mean(ct_values_for_reference_gene_in_sample)
	mean_of_ct_values_for_target_gene_in_control = mean(ct_values_for_target_gene_in_control)
	mean_of_ct_values_for_reference_gene_in_control = mean(ct_values_for_reference_gene_in_control)
	mean_of_ct_values_for_target_gene_in_sample = ifelse( (mean_of_ct_values_for_target_gene_in_sample>=36), 41, mean_of_ct_values_for_target_gene_in_sample )
	mean_of_ct_values_for_reference_gene_in_sample = ifelse( (mean_of_ct_values_for_reference_gene_in_sample>=36), 41, mean_of_ct_values_for_reference_gene_in_sample )
	mean_of_ct_values_for_target_gene_in_control = ifelse( (mean_of_ct_values_for_target_gene_in_control>=36), 41, mean_of_ct_values_for_target_gene_in_control )
	mean_of_ct_values_for_reference_gene_in_control = ifelse( (mean_of_ct_values_for_reference_gene_in_control>=36), 41, mean_of_ct_values_for_reference_gene_in_control )
	delta_CQ_for_sample = mean_of_ct_values_for_target_gene_in_sample - mean_of_ct_values_for_reference_gene_in_sample
	delta_CQ_for_control = mean_of_ct_values_for_target_gene_in_control - mean_of_ct_values_for_reference_gene_in_control
	delta_CQ_for_sample = ifelse( (mean_of_ct_values_for_target_gene_in_sample<41), delta_CQ_for_sample, 41 )
	delta_CQ_for_control = ifelse( (mean_of_ct_values_for_target_gene_in_control<41), delta_CQ_for_control, 41 )
	delta_delta_CQ = delta_CQ_for_sample - delta_CQ_for_control
	delta_delta_CQ = ifelse( (delta_CQ_for_sample<41), delta_delta_CQ, 100 )
	fold_change = 2^(-1 * delta_delta_CQ)
	fold_change
}

number_of_sampling_for_randomisation_test = 10000
fold_change_randomisation_test = fold_change
fold_change_randomisation_test$fold_change_randomisation_test = NA
for (i in seq(1, number_of_sampling_for_randomisation_test)) {
	randomised_fold_change = fold_change_template
	randomised_fold_change$fold_change = apply( randomised_fold_change, 1, calculate_randomised_fold_change, all_values=qpcr )
	randomised_fold_change$fold_change = sprintf( "%3.3f", randomised_fold_change$fold_change )
	old_fold_change_randomisation_test = fold_change_randomisation_test$fold_change_randomisation_test
	if (is.na(old_fold_change_randomisation_test)) {
		fold_change_randomisation_test$fold_change_randomisation_test = randomised_fold_change$fold_change
	} else {
		new_fold_change_randomisation_test = cbind( fold_change_randomisation_test$fold_change_randomisation_test, randomised_fold_change$fold_change )
		fold_change_randomisation_test$fold_change_randomisation_test = new_fold_change_randomisation_test
	}
}

count_number_of_greater_than = function( observed_fold_change, fold_change_randomisation_test )
{
	num_greater = 0
	for (i in seq( 1, length(fold_change_randomisation_test) )) {
		this_randomisation_test = fold_change_randomisation_test[i]
		if (observed_fold_change > 1) {
			if (this_randomisation_test > observed_fold_change) {
				num_greater = num_greater + 1
			}
		} else if (observed_fold_change < 1) {
			if (this_randomisation_test < observed_fold_change) {
				num_greater = num_greater + 1
			}
		}
	}
	total_num = length(fold_change_randomisation_test)
	c( num_greater, total_num )
}

determine_significance = function( params )
{
	observed_fold_change = params[5]
	idx_start = 6
	idx_end = length(params)
	fold_change_randomisation_test = as.numeric(params[idx_start:idx_end])
	how_many_results_greater_than = count_number_of_greater_than( observed_fold_change, fold_change_randomisation_test )
	how_many_greater = how_many_results_greater_than[1]
	how_many_total = how_many_results_greater_than[2]
	pvalue = how_many_greater/how_many_total
	up_or_down = ''
	if (observed_fold_change > 1) {
		up_or_down = 'UP'
	} else if (observed_fold_change < 1) {
		up_or_down = 'DOWN'
	}
	significant_up_or_down = ''
	if (pvalue <= 0.05) {
		significant_up_or_down = up_or_down
	}
	significance_result = c( up_or_down, significant_up_or_down, pvalue )
}

fold_change_significance = fold_change
significance_result = apply( fold_change_randomisation_test, 1, determine_significance )
significance_result_transposed = t(significance_result)
colnames(significance_result_transposed) = c( "up_down", "signif_up_down", "pvalue" )
significance_result_transposed = as.data.frame(as.matrix(significance_result_transposed))
fold_change_significance$up_down = significance_result_transposed$up_down
fold_change_significance$signif_up_down = significance_result_transposed$signif_up_down
fold_change_significance$pvalue = as.numeric(significance_result_transposed$pvalue)
fold_change_significance$pvalue = sprintf( "%3.4f", fold_change_significance$pvalue )

fold_change_significance

# 	   target_gene reference_gene   cell_line control fold_change up_down signif_up_down pvalue
# 	1     ATP5G1V2         ATP5G1   H28-BR0.1     H28       1.662      UP                0.1072
# 	2     ATP5G1V2         ATP5G1   H28-BR0.3     H28       0.204    DOWN           DOWN 0.0095
# 	3     ATP5G1V2         ATP5G1 MeT5A-BR0.1   MeT5A       0.097    DOWN           DOWN 0.0457
# 	4     ATP5G1V2         ATP5G1 MeT5A-BR0.3   MeT5A       0.111    DOWN                0.0562
# 	5     ATP5G1V2         ATP5G1   REN-BR0.1     REN       1.137      UP                0.3482
# 	6     ATP5G1V2         ATP5G1   REN-BR0.3     REN       0.094    DOWN           DOWN 0.0336
# 	7     ATP5G1V2         ATP5G1 VMC23-BR0.3   VMC23       0.132    DOWN           DOWN 0.0446
# 	8       ATP5G2         ATP5G1   H28-BR0.1     H28       0.000    DOWN                0.8458
# 	9       ATP5G2         ATP5G1   H28-BR0.3     H28       0.000    DOWN                0.9452
# 	10      ATP5G2         ATP5G1 MeT5A-BR0.1   MeT5A       0.000    DOWN                0.7939
# 	11      ATP5G2         ATP5G1 MeT5A-BR0.3   MeT5A       0.000    DOWN                0.7971
# 	12      ATP5G2         ATP5G1   REN-BR0.1     REN       1.600      UP                0.1772
# 	13      ATP5G2         ATP5G1   REN-BR0.3     REN       0.086    DOWN           DOWN 0.0393
# 	14      ATP5G2         ATP5G1 VMC23-BR0.3   VMC23       0.282    DOWN                0.0819
# 	15    ATP5G2V3         ATP5G1   H28-BR0.1     H28       0.668    DOWN                0.2004
# 	16    ATP5G2V3         ATP5G1   H28-BR0.3     H28       0.000    DOWN                0.5020
# 	17    ATP5G2V3         ATP5G1 MeT5A-BR0.1   MeT5A       0.000    DOWN                0.9463
# 	18    ATP5G2V3         ATP5G1 MeT5A-BR0.3   MeT5A       0.000    DOWN                0.9489
# 	19    ATP5G2V3         ATP5G1   REN-BR0.1     REN       0.643    DOWN                0.2217
# 	20    ATP5G2V3         ATP5G1   REN-BR0.3     REN       0.481    DOWN                0.1247
# 	21    ATP5G2V3         ATP5G1 VMC23-BR0.3   VMC23       0.000    DOWN                0.5497
# 	22      ATP5G3         ATP5G1   H28-BR0.1     H28       8.605      UP             UP 0.0202
# 	23      ATP5G3         ATP5G1   H28-BR0.3     H28       0.669    DOWN                0.2116
# 	24      ATP5G3         ATP5G1 MeT5A-BR0.1   MeT5A       0.688    DOWN                0.2932
# 	25      ATP5G3         ATP5G1 MeT5A-BR0.3   MeT5A       2.426      UP                0.0571
# 	26      ATP5G3         ATP5G1   REN-BR0.1     REN       0.804    DOWN                0.2840
# 	27      ATP5G3         ATP5G1   REN-BR0.3     REN       0.199    DOWN           DOWN 0.0362
# 	28      ATP5G3         ATP5G1 VMC23-BR0.3   VMC23       0.347    DOWN                0.0803
# 	29       ATP5B         ATP5G1   H28-BR0.1     H28       0.000    DOWN                1.0000
# 	30       ATP5B         ATP5G1   H28-BR0.3     H28       0.000    DOWN                1.0000
# 	31       ATP5B         ATP5G1 MeT5A-BR0.1   MeT5A       0.000    DOWN                1.0000
# 	32       ATP5B         ATP5G1 MeT5A-BR0.3   MeT5A       0.000    DOWN                1.0000
# 	33       ATP5B         ATP5G1   REN-BR0.1     REN       0.869    DOWN                0.3811
# 	34       ATP5B         ATP5G1   REN-BR0.3     REN       0.000    DOWN                0.5092
# 	35       ATP5B         ATP5G1 VMC23-BR0.3   VMC23       0.321    DOWN                0.1615
# 	36     ATP5IF1         ATP5G1   H28-BR0.1     H28       0.895    DOWN                0.3262
# 	37     ATP5IF1         ATP5G1   H28-BR0.3     H28       0.294    DOWN           DOWN 0.0104
# 	38     ATP5IF1         ATP5G1 MeT5A-BR0.1   MeT5A       0.007    DOWN           DOWN 0.0462
# 	39     ATP5IF1         ATP5G1 MeT5A-BR0.3   MeT5A       0.001    DOWN                0.0817
# 	40     ATP5IF1         ATP5G1   REN-BR0.1     REN       0.821    DOWN                0.3347
# 	41     ATP5IF1         ATP5G1   REN-BR0.3     REN       0.152    DOWN           DOWN 0.0369
# 	42     ATP5IF1         ATP5G1 VMC23-BR0.3   VMC23       0.453    DOWN                0.1632
# 	43     RNA18S1         ATP5G1   H28-BR0.1     H28       0.721    DOWN                0.2336
# 	44     RNA18S1         ATP5G1   H28-BR0.3     H28       0.224    DOWN           DOWN 0.0107
# 	45     RNA18S1         ATP5G1 MeT5A-BR0.1   MeT5A       0.004    DOWN           DOWN 0.0414
# 	46     RNA18S1         ATP5G1 MeT5A-BR0.3   MeT5A       0.000    DOWN                0.0504
# 	47     RNA18S1         ATP5G1   REN-BR0.1     REN       0.704    DOWN                0.1624
# 	48     RNA18S1         ATP5G1   REN-BR0.3     REN       0.321    DOWN           DOWN 0.0389
# 	49     RNA18S1         ATP5G1 VMC23-BR0.3   VMC23       0.397    DOWN                0.1197

