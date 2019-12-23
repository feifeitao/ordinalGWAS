#' Load Input Data Files for ordinalGWAS
#'
#' This function loads input genotype, phenotype, and covariate files and creates a data object for ordinalGWAS.
#'
#' @param geno.file Path to the input genotype file (.raw). The file is created by PLINK --recode A or AD options.
#' @param pheno.file Path to the input phenotype file. Optional. If pheno.file is not specified, the phenotype from genotype file will be used for analysis.
#' @param covar.file Path to the input covariate file. Optional. If covar.file is not specified & same.pheno.covar.file is False, the analysis will not include covariates.
#' @param same.pheno.covar.file Logical. Default is False. Set to True if phenotypes and covariates are included in one file. In this case, the pheno.file option will be used for both phenotype and covariate files, and the covar.file option will not be used.
#' 
#' @return The function returns a data object of input datasets. The data object is a list of 4 components: geno, pheno, covar, same.pheno.covar.file.
#' \itemize{
#'   \item geno - Data frame of genotypes. First column is sample ID (FID and IID joined by a space). Second column is the default phenotype from genotype file. Other columns are genotypes.
#'   \item pheno - Data frame of ID and phenotypes. Value is NA if there is no input phenotype file.
#'   \item covar - Data frame of ID and covariates. Value is NA if there is no input covariate file, or if same.pheno.covar.file is set to True.
#'   \item same.pheno.covar.file - Logical. True if phenotypes and covariates are included in one file.
#' }
#'
#' @examples
#' File.A <- system.file("extdata", "example.A.raw", package="ordinalGWAS")
#' File.pheno <- system.file("extdata", "example.pheno.txt", package="ordinalGWAS")
#' File.covar <- system.file("extdata", "example.covar.txt", package="ordinalGWAS")
#' File.pheno.covar <- system.file("extdata", "example.pheno.covar.txt", package="ordinalGWAS")
#'
#' # load a genotype file
#' myObj <- loadData( geno.file=File.A )
#'
#' # load a genotype file and a phenotype file
#' myObj <- loadData( geno.file=File.A, pheno.file=File.pheno )
#'
#' # load a genotype file, and a file containing both phenotypes and covariates
#' myObj <- loadData( geno.file=File.A, pheno.file=File.pheno.covar, same.pheno.covar.file=T )
#'
#' # load a genotype file, a phenotype file, and a covariate file
#' myObj <- loadData( geno.file=File.A, pheno.file=File.pheno, covar.file=File.covar )
#'
#' @export
loadData <- function( geno.file, pheno.file=NA, covar.file=NA, same.pheno.covar.file=F ){

geno <- read.table( geno.file, sep=" ", header=T, stringsAsFactors=F )
N <- ncol(geno)	# 1-6: FID, IID, PAT, MAT, SEX, PHENOTYPE. followed by genotypes
geno$ID <- paste( geno$FID, geno$IID, by=":" )
geno <- geno[ , c( N+1, 6:N) ]	# ID, PHENOTYPE, and genotypes
temp <- list( geno=geno )

if( is.na(pheno.file) ){
temp[["pheno"]] <- NA
}else{
pheno <- read.table( pheno.file, header=T, stringsAsFactors=F )
N <- ncol(pheno)
pheno$ID <- paste( pheno$FID, pheno$IID, by=":" )
pheno <- pheno[ , c( N+1, 3:N ) ]	# ID, phenotypes
temp[["pheno"]] <- pheno
}

if( same.pheno.covar.file==T | is.na(covar.file) ){
temp[["covar"]] <- NA	# if same file for pheno and covar, or covar file not specified
}else{	# covar file is specified
covar <- read.table( covar.file, header=T, stringsAsFactors=F )
N <- ncol(covar)
covar$ID <- paste( covar$FID, covar$IID, by=":" )
covar <- covar[ , c( N+1, 3:N ) ]	# ID, covariates
temp[["covar"]] <- covar
}

temp[["same.pheno.covar.file"]] <- same.pheno.covar.file

return( temp )	# geno, pheno, covar, same.pheno.covar.file

}
