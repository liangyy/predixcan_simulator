library(optparse)

option_list <- list(
    make_option(c("-p", "--pred_expr"), type="character", default=NULL,
                help="Predicted expression (output of Predict.py)",
                metavar="character"),
    make_option(c("-e", "--pve"), type="numeric", default=NULL,
                help="Proportion of variation explained by predicted expression",
                metavar="character"),
    make_option(c("-b", "--beta_dist"), type="character", default=NULL,
                help="Effect size distribution ('infinitesimal', 'spike_n_slab_N', with N be the fraction of non-zeros)",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output simulated phenotype",
                metavar="character"),
    make_option(c("-o", "--output_beta"), type="character", default=NULL,
                help="Output corresponding gene effect sizes",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

.is_spike_n_slab = function(str) {
  ! is.na(stringr::str_match(str, 'spike_n_slab_')[1,1])
}
.get_spike_n_slab_frac = function(str) {
  as.numeric(stringr::str_remove(str, 'spike_n_slab_'))
}

simulate_beta = function(distr, n) {
  if(distr == 'infinitesimal') {
    return(rnorm(n))
  } else if(.is_spike_n_slab(distr)){
    beta = rep(0, n)
    non_zero_ind = runif(n) < .get_spike_n_slab_frac(distr)
    beta[non_zero_ind] = rnorm(sum(non_zero_ind))
    return(beta)
  } else {
    message('beta_dist = ', distr, ' not implemented.')
    quit()
  }
}

calc_var_error = function(ygene, pve) {
  var(ygene) / pve * (1 - pve)
}

library(data.table)
pred_expr = fread(opt$pred_expr, header = TRUE, sep = '\t', data.table = FALSE)
x_mat = pred_expr[, c(-1, -2)]  # predicted expression numbers  
indiv_mat = pred_expr[, 1 : 2]  # individual IDs

# simulate effect sizes
beta = simulate_beta(opt$beta_dist, ncol(x_mat))

# calculate poly gene effect
y_gene = as.matrix(x_mat) %*% beta

# calculate error term
var_error = calc_var_error(y_gene, opt$pve)

# simulate phenotype
y = y_gene + rnorm(nrow(x_mat), sd = sqrt(var_error))

# save y
y = cbind(indiv_mat, data.frame(pheno = y))
write.table(y, opt$output, sep = '\t', quote = FALSE, row = FALSE, col = TRUE)

# save gene-level effect size
df_beta = data.frame(gene_id = colnames(x_mat), effect_size = beta)
write.table(df_beta, opt$output_beta, sep = '\t', quote = FALSE, row = FALSE, col = TRUE)
