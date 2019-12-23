#' Add A1 Allele and Genetic Model in ordinalGWAS Results
#'
#' By default, the SNP names generated from plink --recode A/AD are formatted as SNP_A1 (additive model) or SNP_HET (dominant model). This function splits this column into SNP and A1 in the results. For AD model, a new Model column is added in the results.
#'
#' @param results The result data frame from runAnalysis().
#'
#' @return A new result data frame with A1 allele (for A and AD models) and genetic model (for AD model only). Columns include:
#' \itemize{
#'   \item Phenotype - Phenotype analyzed.
#'   \item SNP - SNP name.
#'   \item A1 - The effect allele.
#'   \item BETA - Beta value from regression analysis.
#'   \item SE - Standard error of the beta value.
#'   \item Tvalue - Test statistic.
#'   \item P - P value.
#'   \item OR - Odds ratio.
#'   \item L95 - Lower boundary of 0.95 confidence interval for odds ratio.
#'   \item U95 - Upper boundary of 0.95 confidence interval for odds ratio.
#'   \item Model - Genetic model (for input data of AD model only): Additive or Dominant.
#' }
#'
#' @examples
#' File.A <- system.file("extdata", "example.A.raw", package="ordinalGWAS")
#' File.AD <- system.file("extdata", "example.AD.raw", package="ordinalGWAS")
#' File.pheno <- system.file("extdata", "example.pheno.txt", package="ordinalGWAS")
#' File.covar <- system.file("extdata", "example.covar.txt", package="ordinalGWAS")
#' File.pheno.covar <- system.file("extdata", "example.pheno.covar.txt", package="ordinalGWAS")
#'
#' # Run analysis for a genotype file of additive model. Split SNP name into SNP and A1.
#' myObj <- loadData( geno.file=File.A, pheno.file=File.pheno, covar.file=File.covar )
#' checkedObj <- checkData( obj=myObj, pheno.name=c("PHENOTYPE1","PHENOTYPE2"),  covar.name=c("COVARIATE1","COVARIATE2") )
#' results <- runAnalysis( obj=checkedObj )
#' new.results <- splitSNP(results)
#'
#' # Run analysis with a genotype file of AD model. Add A1 and genetic model columns.
#' myObj <- loadData( geno.file=File.AD, pheno.file=File.pheno, covar.file=File.covar )
#' checkedObj <- checkData( obj=myObj, pheno.name=c("PHENOTYPE1","PHENOTYPE2"),  covar.name=c("COVARIATE1","COVARIATE2") )
#' results <- runAnalysis( obj=checkedObj )
#' new.results <- splitSNP(results)
#'
#' @export
splitSNP <- function(results){

x <- results$SNP
x <- gsub( "^X", "", x )
# if SNP format is chr:bp (eg. 15:10001), there is a leading X when reading raw file
temp <- strsplit( x, "_" )
results$SNP <- sapply( temp, "[[", 1 )
a1 <- sapply( temp, "[[", 2 )

if( ! "HET" %in% a1 ){	# additive model, raw file from plink --recode A
	results$A1 <- a1
	results <- results[ , c( "Phenotype", "SNP", "A1", 
		"BETA", "SE", "Tvalue", "P", "OR", "L95", "U95" ) ]
}else{	# additive and dominant model, raw file from plink --recode AD
	isDominant <- a1=="HET"
	isAdditive <- !isDominant
	a1 <- a1[isAdditive]
	results$A1 <- c(rbind(a1,a1))
	results$Model <- "Additive"
	results$Model[isDominant] <- "Dominant"
	results <- results[ , c( "Phenotype", "SNP", "A1", 
		"BETA", "SE", "Tvalue", "P", "OR", "L95", "U95", "Model" ) ]	
}	# Model column is added in AD model results

return(results)

}
