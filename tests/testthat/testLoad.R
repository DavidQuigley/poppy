#library( poppy )
#library( GenomicRanges )
context("loading code")

test_that("somatic data load works", {
    pkg = "poppy"

    sample_id = "CA071-01_WGS"
    sample_base = "CA-071"
    dir_data_root=test_path( sample_id )
    fn_segs = paste0( dir_data_root, "/", sample_id, ".purple.cnv.somatic.tsv")
    fn_purity = paste0( dir_data_root, "/", sample_id, ".purple.purity.tsv" )
    fn_manta = paste0(dir_data_root,'/somaticSV_with_INV.vcf')
    fn_snp_strelka = paste0( dir_data_root, "/somatic.snvs_PASS_snpeff_dbsnp_clinvar.txt")
    fn_indel_strelka = paste0( dir_data_root, "/somatic.indels_PASS_snpeff_dbsnp_clinvar.txt")
    fn_snp_mutect = paste0( dir_data_root, "/",sample_id,"_mutect_filtered_PASS_snpeff_dbsnp_clinvar.txt")
    fn_cna_by_gene = paste0( dir_data_root, "/",sample_id,".purple.cnv.gene.tsv")

    expect_equal( file.exists( fn_segs ), TRUE )
    expect_equal( file.exists( fn_purity ), TRUE )
    expect_equal( file.exists( fn_manta ), TRUE )
    expect_equal( file.exists( fn_snp_strelka ), TRUE )
    expect_equal( file.exists( fn_indel_strelka ), TRUE )
    expect_equal( file.exists( fn_snp_mutect ), TRUE )
    expect_equal( file.exists( fn_cna_by_gene ), TRUE )

    library_type = "WGS"

    expect_error( poppy::Somatic( sample_base=sample_base, sample_id=sample_id, library_type="bogus") )
    somatic = poppy::Somatic( sample_base=sample_base, sample_id=sample_id, library_type=library_type)
    expect_equal( somatic$sample_base, sample_base )
    expect_equal( somatic$sample_id, sample_id )

    min_reads = 5
    somatic = poppy::add_manta(somatic, fn_manta, MIN_READS = min_reads)
    expect_equal( length( somatic$manta ), 662 )
    nt_for_microhomology = 3
    somatic = poppy::calculate_manta_summary( somatic, nt_for_microhomology )
    expect_equal( dim(somatic$list_sv)[1], 662 )
    expect_equal( sum( somatic$list_sv$svtype=="DEL" ), 413 )
    somatic = poppy::add_SNV_strelka( somatic,
                                      fn_strelka_clinvar_snv = fn_snp_strelka,
                                      fn_strelka_clinvar_indel = fn_indel_strelka )
    somatic = poppy::add_SNV_mutect( somatic, fn_mutect_clinvar_snv = fn_snp_mutect )
    somatic = poppy::calculate_strelka_properties( somatic )
    expect_equal( dim(somatic$strelka)[1], 303052 )
    somatic$strelka[somatic$strelka$gene=="TP53",]
    somatic = poppy::somatic_add_CNA_by_gene_purple(somatic, fn_cna_by_gene )
    expect_equal( dim(somatic$CNA_genes)[1], 25417 )
    expect_equal( somatic$CNA_genes$maxCopyNumber[ dimnames( somatic$CNA_genes)[[1]] == "BRCA2" ], 2.0005)
    somatic = poppy::somatic_add_CNA_purple( somatic, fn_segs, fn_purity )
    expect_equal(dim(somatic$segments)[1], 871)
} )
