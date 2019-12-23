#' Run Ordered Logicstic Regression Analysis
#'
#' This function performs ordered logistic regresssion analysis with data and analysis options from checkData().
#'
#' @param obj The data object created by checkData() containing input data and analysis options (phenotypes and covariates).
#'
#' @return A data frame of analysis results, including the following columns: Phenotype, SNP, BETA, SE, Tvalue, P, OR, L95, U95.
#' \itemize{
#'   \item Phenotype - Phenotype analyzed.
#'   \item SNP - SNP name. By default, the SNP name is formatted as SNP_A1 (additive model) or SNP_HET (dominant model).
#'   \item BETA - Beta value from regression analysis.
#'   \item SE - Standard error of the beta value.
#'   \item Tvalue - Test statistic.
#'   \item P - P value.
#'   \item OR - Odds ratio.
#'   \item L95 - Lower boundary of 0.95 confidence interval for odds ratio.
#'   \item U95 - Upper boundary of 0.95 confidence interval for odds ratio.
#' }
#'
#' @examples
#' File.A <- system.file("extdata", "example.A.raw", package="ordinalGWAS")
#' File.pheno <- system.file("extdata", "example.pheno.txt", package="ordinalGWAS")
#' File.covar <- system.file("extdata", "example.covar.txt", package="ordinalGWAS")
#' File.pheno.covar <- system.file("extdata", "example.pheno.covar.txt", package="ordinalGWAS")
#'
#' # Run analysis with default phenotype in the genotype file
#' myObj <- loadData( geno.file=File.A )
#' checkedObj <- checkData( obj=myObj )
#' results <- runAnalysis( obj=checkedObj )
#'
#' # Run analysis with a pheno file and include all phenotypes
#' checkedObj <- checkData( obj=myObj, all.pheno=T )
#' results <- runAnalysis( obj=checkedObj )
#'
#' # Run analysis pheno and covar included in the same file
#' myObj <- loadData( geno.file=File.A, pheno.file=File.pheno.covar, same.pheno.covar.file=T )
#' checkedObj <- checkData( obj=myObj, pheno.name="PHENOTYPE1", covar.name="COVARIATE1" )
#' results <- runAnalysis( obj=checkedObj )
#'
#' # Run analysis with a pheno file and a covar file
#' myObj <- loadData( geno.file=File.A, pheno.file=File.pheno, covar.file=File.covar )
#' checkedObj <- checkData( obj=myObj, pheno.name=c("PHENOTYPE1","PHENOTYPE2"),  covar.name=c("COVARIATE1","COVARIATE2") )
#' results <- runAnalysis( obj=checkedObj )
#'
#' @export
runAnalysis <- function(obj){

# unpack obj
SNPs <- obj[["SNPs"]]
pheno.name <- obj[["pheno.name"]]
covar.name <- obj[["covar.name"]]
merged <- obj[["merged"]]

# final.result: a list of df. each df is result of one PHEN
final.result <- list()

# result0: initialize the result matrix for one PHEN
zero <- numeric( length(SNPs) )
result0 <- data.frame( BETA=zero, SE=zero, Tvalue=zero, P=zero, 
    OR=zero, L95=zero, U95=zero )    # 7 columns
rownames(result0) <- SNPs    # add row names: SNPs

for( PHEN in pheno.name ){	# loop starts for one phenotype

	result <- result0

	# loop starts for one snp
	for( snp in SNPs ){

	# build model
	fmla <- paste( PHEN, "~", snp )
	if( !identical(covar.name,NA) ){
	to.add <- paste( covar.name, collapse=" + " )
	fmla <- paste( fmla, "+", to.add )
	}
	myModel <- MASS::polr( fmla, data=merged, Hess=TRUE)

	# coef and P
	ctable <- coef(summary(myModel))    # coef table
	temp <- ctable[1,]  # a named vector. first row is the snp data: Value, Std. Error, t value
	p <- pnorm(abs(temp["t value"]), lower.tail = FALSE) * 2
	names(p) <- "P" # compute P
	temp <- c(temp,p)   # 4 items: Value, Std. Error, t value, P

	# OR
	OR <- exp( temp["Value"] )
	names(OR) <- "OR"

	# 95% CI for OR
	ci <- confint.default(myModel, parm=snp)    # this is 95% CI for BETA
	# use R base version of confint instead of MASS version # note there is a space in parm
	ciOR <- exp( ci )   # this is 95% CI for OR
	# see notes below on confint function

	# combine and write to result
	temp <- c( temp, OR, ciOR ) # Value, Std. Error, t value, P, OR, L95, U95
	temp <- unname(temp)    # 7 items
	result[ snp, ] <- temp  # add to result

	} # loop ends for one snp

result$SNP <- rownames(result)
result <- result[,c(8,1:7)]	# new order: SNP, and result data
final.result[[PHEN]] <- result

} # loop ends for PHEN

to.return <- plyr::ldply(final.result, rbind, .id="Phenotype")
# merge all tables in the list into one df, add the first column of Phenotype (names of the list)
return( to.return )

}
