#' Wrapper function to add annotated germline calls from HaplotypeCaller (GATK4)
#' @param somatic somatic data object to modify
#' @param fn.HC.germline path annotated germline text file (e.q. SAMPLE_ID_snpeff_dbsnp_clinvar_gnomAD.txt)
#' @export
#' library( VariantAnnotation )
#' library( stringr )
#' library( plyr )
#' library( poppy )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Re-purposed read_germline_pathogenic function
# function name read_germline_HC_pathogenic 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


read_germline_HC_pathogenic = function( somatic, fn_path ){
    gm = read.table( fn.HC.germline, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    names(gm) = c("chrom", "pos", "ref", "alt", "symbol", "ensembl_transcript", "effect", "aa", "clinvar", "ad_ref", "ad_alt")
    gm = gm[gm$clinvar=="Pathogenic"| gm$clinvar=="Pathogenic/Likely_pathogenic"| gm$clinvar=="Pathogenic|_protective",]
    positions = unique(gm$pos)
    chroms = rep("", length(positions))
    symbols = rep("", length(positions))
    effects = rep("", length(positions))
    clinvars = rep("", length(positions))
    somatic = c(somatic, list( "germline_HC_pathogenic" = gm ) )
    somatic
}
