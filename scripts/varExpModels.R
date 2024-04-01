library(variancePartition)
library(ggplot2)

create_random_effects_formula <- function(vectors) {
  # Combine the vectors into a formula string

  formula_str <- paste("(1 | ", vectors, ")")
  formula_str <- paste(formula_str, collapse = " + ")
  # Convert the string to a formula object
  #formula_obj <- as.formula(formula_str)
  return(formula_str)
}

RowVar <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
pca.df <- read.csv(snakemake@input[[1]][1], row.names=1)
libstats.df <- read.csv(snakemake@input[[2]][1], row.names=1)
#Rescale total bp such that it will allow the model to converge downstream -- should not affect results
libstats.df$Z_total_bp <- (libstats.df$Total_bp - mean(libstats.df$Total_bp) ) / sd(libstats.df$Total_bp)
lme.df <- cbind(pca.df, libstats.df)

#get meta-fields:
meta.field <- c()
pc.field <- c()
for (meta in colnames(pca.df )){
  if (!startsWith(meta, "PC") ){
    meta.field <- c(meta.field, meta)
    lme.df[,meta] <- as.factor(lme.df[,meta])
  } else {pc.field <- c(meta, pc.field)} }

y.df <- t(pca.df[,pc.field]) #transform PC df to right format

#we need to filter out response variables with variance < 1e-05 so the model can converge
variance_results <- RowVar(y.df) > 1e-05
y.df <- y.df[variance_results, ]

#call the mixed-effect model and calculte variance explained by each component:
re.form <- create_random_effects_formula(meta.field)
form <- as.formula(paste("~ Z_total_bp + ", re.form))

#use variance partition to analyze each component
varPart <- fitExtractVarPartModel(y.df, form, lme.df)
write.csv(varPart, snakemake@output[[2]][1])
# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart, decreasing=TRUE)

# Figure 1a
# Bar plot of variance fractions for all PCs
barplot <- plotPercentBars(vp)
ggsave(snakemake@output[[1]][1], dpi=300)