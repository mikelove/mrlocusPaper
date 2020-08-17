# Simulation using twas_sim and Snakemake

* `Snakemake` - provides all code used to run simulation and the
  methods evaluated
* `sims/simulation_files_240.tgz` - a gzipped tar file containing all
  the 240 simulation files
    - These simulation files are empty `.sim` files which specify
      simulation settings and random seed by their filename
	- E.g. `1_7r_500nq_1pctm_0.1h2_0.01ve` has random seed 7, n=500
      for the eQTL, 1% of SNPs are causal (eQTL) in the locus, the
      gene heritability is 10% and the percent variance explained of
      the trait by the gene is 1%.
	- The simulations are numbered 1-9 and null1, null2, null3.
	  For a mapping of the numbers to the letters used in the paper,
	  see `sim_review.R`
* `final_plot.R` - code used to make Figure 2 and Supplementary
  Figures of simulation performance metrics
* `null_plot.R` - code used to assess the 3 null simulation settings
* `mrlocus.R`, `twmr.R`, `pmr_summary_egger.R`, `lda_mr_egger.R` - for
  R-based methods, the usage of the R functions are provided in these
  files
* `sim_review.R` - code used for making the diagram explaining the simulations
* `timing.R` - code used for making the timing plot
