#' Load Analysis Options and Check Input Data for ordinalGWAS
#'
#' This function reads analysis options (phenotypes/covariates), and checks the input data object from loadData().
#'
#' @param obj The data object created by loadData().
#' @param pheno.name A character vector of phenotypes for analysis. Optional.
#' @param covar.name A character vector of covariates for analysis. Optional.
#' @param all.pheno Logical. Default is False. If set to True, all phenotypes in the phenotype file will be analyzed.
#' @param all.covar Logical. Default is False. If set to True, all covariates in the covariate file will be included in the analysis.
#' 
#' @return The function returns a data object containing input data and analysis options. The data object is a list of 4 components: SNPs, pheno.name, covar.name, merged.
#' \itemize{
#'   \item SNPs - A character vector of SNPs in the input data.
#'   \item pheno.name - A character vector of phenotypes for analysis.
#'   \item covar.name - A character vector of covariates for analysis.
#'   \item merged - The full data frame of sample IDs, genotypes, phenotypes and covariates.
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
#' # check data
#' checkedObj <- checkData( obj=myObj )
#'
#' # load a genotype file, and a file for both phenotypes and covariates
#' myObj <- loadData( geno.file=File.A, pheno.file=File.pheno.covar, same.pheno.covar.file=T )
#' # select phenotype and covariate for analysis, and check data
#' checkedObj <- checkData( obj=myObj, pheno.name="PHENOTYPE1", covar.name="COVARIATE1" )
#' 
#' # load a genotype file, a phenotype file, and a covariate file
#' myObj <- loadData( geno.file=File.A, pheno.file=File.pheno, covar.file=File.covar )
#' # select phenotypes and covariates for analysis, and check data
#' checkedObj <- checkData( obj=myObj, pheno.name=c("PHENOTYPE1","PHENOTYPE2"), covar.name=c("COVARIATE1","COVARIATE2") )
#'
#' @export
checkData <- function( obj, pheno.name=NA, covar.name=NA, all.pheno=F, all.covar=F ){

# unpack obj
geno <- obj[["geno"]]
pheno <- obj[["pheno"]]
covar <- obj[["covar"]]
same.pheno.covar.file <- obj[["same.pheno.covar.file"]]
SNPs <- colnames(geno)[3:ncol(geno)]

# merge geno and pheno (with selected phenotypes), set pheno.name
if( identical(pheno,NA) ){
	pheno.name <- "PHENOTYPE"
	merged <- geno
} else if ( !(identical(pheno.name,NA)) & all.pheno==T ){
	print( "Options pheno.name and all.pheno cannot be used together." )
	return(NULL)
} else if ( !(identical(pheno.name,NA)) & same.pheno.covar.file==T ){
	all.names <- c(pheno.name,covar.name)
	subpheno <- pheno[ , c( "ID", all.names ) ]
	merged <- merge( geno, subpheno, by="ID", all.x=F, all.y=F )
} else if ( !(identical(pheno.name,NA)) ){
	subpheno <- pheno[ , c( "ID", pheno.name ) ]
	merged <- merge( geno, subpheno, by="ID", all.x=F, all.y=F )
} else if ( all.pheno==T ){
	pheno.name <- colnames(pheno)[2:ncol(pheno)]
	merged <- merge( geno, pheno, by="ID", all.x=F, all.y=F )
}

# add covar into merged
if( same.pheno.covar.file==T | identical(covar,NA) ){
	abc <- 0	# no changes to merged or covar.name
} else if ( !(identical(covar.name,NA)) & all.covar==T ){
	print( "Options covar.name and all.covar cannot be used together." )
	return(NULL)
} else if ( !(identical(covar.name,NA)) ){
	if( length(intersect(covar.name,pheno.name))>0 ){
		print( "Phenotypes and covariates cannot have same variable names." )
		return(NULL)
	} else {
		subcovar <- covar[ , c("ID",covar.name) ]
		merged <- merge( merged, subcovar, by="ID", all.x=F, all.y=F )
	}
} else if ( all.covar==T ){
	covar.name <- colnames(covar)[2:ncol(covar)]
	if( length(intersect(covar.name,pheno.name))>0 ){
		print( "Phenotypes and covariates cannot have same variable names." )
		return(NULL)
	} else { merged <- merge( merged, covar, by="ID", all.x=F, all.y=F ) }
}

# phenotypes: change -9 to NA, and change to ordered factor
for( x in pheno.name ){
check <- merged[,x]==-9
check[ is.na(check) ] <- F
merged[check,x] <- NA
merged[,x] <- as.ordered(merged[,x])
}

# covariates: change -9 to NA
if( !identical(covar.name,NA) & same.pheno.covar.file==F ){
for( x in covar.name ){
	check <- merged[,x]==-9
	check[ is.na(check) ] <- F
	merged[check,x] <- NA
}}

to.return <- list( SNPs=SNPs, pheno.name=pheno.name, covar.name=covar.name, merged=merged )
return(to.return)

}
