## set up command-line options
option_list = list(
  optparse::make_option(c("-o", "--outdir"), type="character", default=NULL,
                       help="Output directory"),
  optparse::make_option("--true_ibd", type="character", default=NULL,
                       help="True IBD segments TSV"),
  optparse::make_option("--inf_ibd", type="character", default=NULL,
                       help="Inferred IBD segments TSV"),
  optparse::make_option("--true_sfs", type="character", default=NULL,
                       help="True SFS TSV"),
  optparse::make_option("--est_sfs", type="character", default=NULL,
                       help="Estimated SFS TSV"),
  optparse::make_option("--true_ne", type="character", default=NULL,
                       help="True Ne TSV"),
  optparse::make_option("--est_ne", type="character", default=NULL,
                       help="Estimated Ne TSV"),
  optparse::make_option("--min_overlap_frac", type="character", default="0.5",
                       help="Minimum segment overlap fraction [default=%default]"),
  # Transmission chain analysis options
  optparse::make_option("--true_pairs", type="character", default=NULL,
                       help="True transmission pairs TSV"),
  optparse::make_option("--inf_pairs", type="character", default=NULL,
                       help="Inferred transmission pairs TSV"),
  optparse::make_option("--directed", type="character", default="FALSE",
                       help="Whether transmission pairs are directed (TRUE/FALSE) [default=%default]"),
  # Outbreak cluster analysis options  
  optparse::make_option("--true_clusters", type="character", default=NULL,
                       help="True outbreak clusters TSV"),
  optparse::make_option("--inf_clusters", type="character", default=NULL,
                       help="Inferred outbreak clusters TSV"),
  # Genetic distance calibration options
  optparse::make_option("--inf_links", type="character", default=NULL,
                       help="Inferred links with genetic distances TSV")
);