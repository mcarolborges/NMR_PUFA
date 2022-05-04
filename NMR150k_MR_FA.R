############################################################################
#
#                      MR of fatty acids on CVD
#
############################################################################

############################################################################
#                                Set - UP                                  #
############################################################################


##
## Set working directory 
setwd(Sys.getenv('FA_CVD'))


##
## Clear the work environment
rm(list = ls())


##
## Setting repository (UoB)
options(repos = c(CRAN ="http://www.stats.bris.ac.uk/R/"))


##
## Setting digits
options(digits = 10)


##
## Library directory
.libPaths()


##
## Updating packages 
#update.packages(ask = "FALSE")


##
## Install packages
#install.packages(c("data.table", "purrr", "devtools", "xlsx", "dplyr", "ggplot2"))
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments")
#devtools::install_github("NightingaleHealth/ggforestplot")
#library(remotes)
#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)


##
## Load library
library(data.table)
library(purrr)
library(TwoSampleMR)
library(MRInstruments)
library(xlsx)
library(dplyr)
library(ggplot2)
library(ggforestplot)


##
## Key directories

dat.dir <- paste0(Sys.getenv('MR_FA'), "/data/")


############################################################################
#                  Extract SNPs and SNP-exposure data
############################################################################

##
## Fatty acids traits in UKBB and Kettunen GWAS

ao <- available_outcomes()

fa_ls <- list(
  DHA     = list(UKBB = 'ukb-DHA',     Kettunen = 'met-c-852'),
  Omega_3 = list(UKBB = 'ukb-Omega_3', Kettunen = 'met-c-855'),
  LA      = list(UKBB = 'ukb-LA',      Kettunen = 'met-c-893'),
  Omega_6 = list(UKBB = 'ukb-Omega_6', Kettunen = 'met-c-856')
)


##
## Function to extract SNP-exposure data 

extr_expdat <- function(fa_trait, excl.fads) {
  
  ## Discovery sample
  dis <- "UKBB"
  
  ## Fatty acid trait ID in discovery sample
  id.dis <- fa_ls[[fa_trait]][[dis]]
  
  ## Select SNPs and extract association data in discovery sample
  f <- paste0(fa_trait, "_int.imputed.txt.gz")
  dis.dat <- fread(paste0(dat.dir, "nmr_ukb_150k/", f)) %>%
    mutate(P_BOLT_LMM_INF = as.numeric(P_BOLT_LMM_INF)) %>%
    filter(., INFO > 0.8 & P_BOLT_LMM_INF < 5e-8 & A1FREQ < 0.99 & A1FREQ > 0.01) %>%
    format_data(., 
                type = "exposure",
                snp_col = "SNP",
                beta_col = "BETA",
                se_col = "SE",
                eaf_col = "A1FREQ",
                effect_allele_col = "ALLELE1",
                other_allele_col = "ALLELE0",
                pval_col = "P_BOLT_LMM_INF",
                info_col = "INFO",
                chr_col = "CHR",
                pos_col = "BP"
    ) %>%
    clump_data(., clump_kb = 10000, clump_r2 = 0.001) %>%
    mutate(id.exposure = id.dis, samplesize.exposure = 114999)
  
  ## Exclude FADS SNPs if excl.fads=="yes"
  if(excl.fads=="yes") {
    # Identify FADS SNPs
    fads.snps <- filter(dis.dat, chr.exposure == 11 & pos.exposure > 61567097-500000 & pos.exposure < 61634826+500000) %>% pull(SNP)
    # Remove FADS SNPs from association data in discovery sample
    dis.dat <- filter(dis.dat, ! SNP %in% fads.snps) 
  }								
  
  ## Indicate if FADS SNPs were excluded from the dataset
  dis.dat <- mutate(dis.dat, exposure = fa_trait, sample = dis, exclfads = excl.fads, chr.exposure = as.character(chr.exposure))
}


##
## Extract SNP-exposure data from UKBB/Kett (w/o and w FADS SNPs) 

# List of arguments for extr_expdat
expdat_args <- cross_df(list(fa_trait = names(fa_ls), excl.fads = c("no", "yes")))

# Extract exposure data
expdat <- pmap(expdat_args, extr_expdat) %>% bind_rows


############################################################################
#            Check R2/F of UKBB SNPs on UKBB Kettunen data 
############################################################################

##
## Function to approximate R2 and F statistics in discovery and replication samples

check_ivrep <- function(fa_trait, excl.fads) {
  
  ## Discovery sample
  dis <- "UKBB"
  
  ## Replication sample
  rep <- "Kettunen"
  
  ## Fatty acid trait ID in discovery sample
  id.dis <- fa_ls[[fa_trait]][[dis]]
  
  ## Fatty acid trait ID in replication sample
  id.rep <- fa_ls[[fa_trait]][[rep]]
  
  ## Extract association data in discovery sample for selected SNPs
  dis.dat <- filter(expdat, id.exposure == id.dis & exclfads == excl.fads)
  
  ## Extract association data in replication sample for selected SNPs
  rep.dat <- extract_outcome_data(dis.dat$SNP, outcomes = id.rep)
  
  ## Harmonise discovery and replication data
  dis.rep.dat <- harmonise_data(dis.dat, rep.dat) 
  # P.S.: assuming forward strand as ambiguous SNPs incorrectly harmonised should not affect R2/F estimation
  
  ## Check if there are SNPs with harmonisation problems (e.g. incompatible alleles)
  print(table(dis.rep.dat$mr_keep == F))
  
  ## Calculate R2 and F for each SNP for discovery and replication
  dis.rep.dat <- mutate(dis.rep.dat,
                        r2_dis = 2 * (beta.exposure^2) * eaf.exposure * (1-eaf.exposure),
                        f_dis  = ((beta.exposure/se.exposure)^2),
                        r2_rep = 2 * (beta.outcome^2) * eaf.outcome * (1-eaf.outcome),
                        f_rep  = ((beta.outcome/se.outcome)^2)
  )
  # P.S.: betas and SEs in SD units for both discovery and replication datasets
  
  ## Calculate total R2 and mean F across SNPs for discovery and replication
  iv.sumdat <- summarise(dis.rep.dat,
                         id_discovery   =  unique(id.exposure),
                         id_replication =  unique(id.outcome),
                         Nsnps_available = n(), 
                         R2_Discovery   =  round(sum(r2_dis), 3), 
                         R2_Replication =  round(sum(r2_rep), 3), 
                         F_Discovery    =  round(mean(f_dis), 0), 
                         F_Replication  =  round(mean(f_rep), 0)
  ) %>%
    mutate(., 
           Fatty_acid = fa_trait,
           Discovery = dis,
           Replication = rep,
           Nsnps_selected = nrow(dis.dat),
           FADS_excluded = excl.fads
    ) %>%
    arrange(Replication, FADS_excluded) %>%
    select(Fatty_acid, Discovery, Replication, FADS_excluded, Nsnps_selected, Nsnps_available, R2_Discovery, R2_Replication, F_Discovery, F_Replication, id_discovery, id_replication)														
  
}


##
## Approximate R2 and F statistics in discovery and replication samples

# List of arguments for check_ivrep
ivrep_args <- cross_df(list(fa_trait = names(fa_ls), excl.fads = c("no", "yes")))

# Extract R2 and F statistics
ivrep <- pmap(ivrep_args, check_ivrep) %>% bind_rows 


############################################################################
#               Import and format fatty acids data from published GWAS
############################################################################

##
## GWAS consortia

fa.source <- c("charge", "ket", "shin")


##
## Location of Fatty acids GWAS data from "charge", "ket", "shin"

d <- paste0(Sys.getenv('MR_FA'), "/data/", fa.source, "/hg19/")                

f.path <- list.files(path = d, pattern = ".tab$", full.names = T)


##
## Read and format fatty acids GWAS data 

fname <- stringr::str_split(f.path, "/hg19//", simplify = TRUE)[,2] # File names need to be in the same order as in f.path

fa.outdat <- map(f.path, ~read_outcome_data(., 
                                            snps = unique(expdat$SNP), 
                                            sep = "\t", 
                                            snp_col = "snp", 
                                            eaf_col = "effect_allele_freq",
                                            pval_col = "p",
                                            samplesize_col = "n"
)) %>%
  map2(., fname, ~mutate(.x, outcome = .y)) %>%
  bind_rows 


##
## Fatty acids map file

fa.info <- read.xlsx("data/fatty acid instruments_v2.xlsx", sheetName = "FA fingerprint")  
fa.info <- data.table::melt(fa.info, 
                            id.vars = c("abbreviation", "Label", "order", "Formula", "Group", "Configuration", "HMDB.ID", "PubChem.CID", "IUPAC.name"), 
                            measure =c("filename_charge", "filename_ket", "filename_shin"), 
                            variable.name = "source", 
                            value.name = c("filename")
) %>%
  filter(!is.na(filename))


############################################################################
#               Check IV x FA association                                      
############################################################################


harmo_fadat <- function(id, excl.fads) {
  
  print(paste("Doing:", id, "; FADS excluded:", excl.fads))
  
  ## Harmonise effect alleles
  dat <- filter(expdat, id.exposure == id & exclfads == excl.fads) %>%
    harmonise_data(., fa.outdat, action = 2) %>%
    mutate(exclfads = excl.fads)
}


##
## Harmonise data for each SNP subset (1: UKB; 2: UKB-noFADS)

fadat_arg <- cross_df(list(id = unique(expdat$id.exposure), excl.fads = c("no", "yes")))

fa_dat <- pmap(fadat_arg, harmo_fadat) %>% bind_rows


##
## Function to run IVW for each SNP subset (1: UKB; 2: UKB-noFADS)

run_ivw <- function(excl.fads) {
  
  ## Subset according to exclfads
  dat <- filter(fa_dat, exclfads == excl.fads)
  
  ## Run IVW
  res <- mr(dat, method_list = "mr_ivw") %>%
    mutate(z = b / se) 
  
  ## Generate between-SNP Heterogeneity statistics
  res <- mr_heterogeneity(dat, method_list =  "mr_ivw") %>%
    select(id.exposure, id.outcome, Q, Q_df, Q_pval) %>%
    merge(res,., by = c("id.exposure", "id.outcome")) %>%
    mutate(exclfads = excl.fads)
  
  ## Create indicator variable for strength of SNP-fatty acids association
  res$stars <- cut(res$pval, breaks=c(-Inf, 5e-8, 5e-5, 5e-2, Inf), label=c("***", "**", "*", "")) 
  
  ## Merge data with fatty acids information file			
  res <- res %>%
    base::merge(., fa.info, by.x = "outcome", by.y = "filename") %>%
    arrange(order) %>%
    mutate(
      abbreviation = factor(abbreviation, levels = rev(factor(unique(abbreviation)))),
      xlab = case_when(source=="filename_charge" ~ "CHARGE\nPlasma FA\n(N = 8,631)\n[GC]",
                       source=="filename_ket" ~  "Kettunen2016\nPlasma FA\n(N = 13,516)\n[NMR]",
                       source=="filename_shin" ~  "Shin2014\nPlasma FA\n(N = 7,352)\n[MS]",
      )
    )
}		


##
## Run IVW for each SNP subset (1: UKB; 2: UKB-noFADS)

fa_res <- map(c("no", "yes"), ~run_ivw(.)) %>% bind_rows


##
## Function to create heatmaps IV vs FA

mk_heatmap <- function(excl.fads, id, grp1, grp2) {
  
  ## Exposure name
  exp <- filter(fa_res, id.exposure == id) %>% pull(exposure) %>% as.character %>% unique
  
  ## Sample (UKBB vs Kettunen)
  samp <- "UKBB"
  
  ## FA class
  fa_grp <- "PUFA"
  
  ## Subset data
  df <- filter(fa_res, exclfads == excl.fads & id.exposure == id & (Group == grp1 | Group == grp2))
  
  ## Make plot			
  p <- ggplot(df, aes(x = xlab, y = abbreviation, fill = z)) +
    geom_tile() +
    colorspace::scale_fill_continuous_diverging("Blue-Red 3", limits = c(-10, 25)) +
    theme_minimal() +
    theme(
      legend.title = element_text(size=16),
      legend.text = element_text(size=11),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      strip.text.y = element_text(size=12),
      panel.grid.major = element_line(colour = "white", size=3),
      plot.margin = unit(c(1, 0.5, 1, 0.5), "cm")) + 
    geom_text(aes(label=stars), color="black", size=6, angle = 180) + 
    facet_grid(Group ~ ., scales = "free") 		
  
  ## Save plot
  ggsave(paste0("./output/graphs/heatmap_150kIV_", samp, "_", exp, "_vs_", fa_grp, "_exclFADS-", excl.fads, ".png"), plot = p)
}


##
## Heatmap for each IV vs fatty acids

# List of arguments for heatmap
heatmap_args <- cross_df(list(excl.fads = c("no", "yes"), id = unique(expdat$id.exposure)))

# IVs vs PUFAs
walk2(heatmap_args[[1]], heatmap_args[[2]], ~mk_heatmap(excl.fads = .x, id = .y, grp1 = "n-3 PUFA", grp2 = "n-6 PUFA")) 


##
## Save results
writexl::write_xlsx(list(iv_strength = ivrep, 
                         exp_dat = expdat,
                         iv_fa = fa_dat, 
                         iv_fa_res = fa_res), 
                    path = "./output/IV_nmr150k.xlsx")


############################################################################
#               Import and format CVD GWAS data
############################################################################

##
## Information file for CVD endpoints
cvd.info <- read.xlsx(paste0(Sys.getenv('SCRATCH'), "/gwas_sum_data/cvd/cvd_info.xlsx"), sheetName = "CVD") %>%
  filter(type_outcome == "primary")


##
## GWAS traits
gwas.traits <- cvd.info %>%
  filter(!is.na(file) & (source == "GWAS" | source == "GWAS+UKB")) %>%
  pull(abbrev) %>% 
  as.character 

##
## UKBB traits
ukb.traits <- cvd.info %>%
  filter(!is.na(file) & source == "UKB") %>%
  pull(abbrev) %>% 
  as.character				


##
## FinnGen traits
finngen.traits <- cvd.info %>%
  filter(source == "FinnGen") %>%
  pull(abbrev) %>% 
  as.character	

##
## Metanalysed traits
pool.traits <- cvd.info %>%
  filter(!is.na(file) & source == "Pooled") %>%
  pull(abbrev) %>% 
  as.character	

##
## GWAS files		
gwas.files <- paste0(dat.dir, "cvd/", "gwas/", gwas.traits, "/", gwas.traits, "_gwas.tab.gz")


##
## UKBB files
ukb.files <- paste0(dat.dir, "cvd/", "ukb/", ukb.traits, "/", ukb.traits, "_ukb.tab.gz")


##
## FinnGen files
finngen.files <- paste0(dat.dir, "finngen/", finngen.traits, "/", finngen.traits, "_finngen.tab.gz")


##
## Metanalysed files
pool.files <- paste0(dat.dir, "cvd_pooled/", pool.traits, "/", pool.traits, "_pooled.tab.gz")


##
## Extract SNP-outcome summary data for relevant SNPs

extr.out.dat <- function(f, y, z) {
  
  df <- fread(f) %>%
    format_data(.,
                snps = unique(expdat$SNP),
                type = "outcome", 
                phenotype_col = "trait")
  
  info <- filter(cvd.info, abbrev == y & source %in% z)
  
  df <- merge(df, info, by.x="outcome", by.y="abbrev")
  
}

gwas.dat <- map2(gwas.files, gwas.traits, ~extr.out.dat(.x, .y, z = c("GWAS", "GWAS+UKB"))) %>%
  bind_rows 

ukb.dat <- map2(ukb.files, ukb.traits, ~extr.out.dat(.x, .y, z = "UKB")) %>%
  bind_rows 

finngen.dat <- map2(finngen.files, finngen.traits, ~extr.out.dat(.x, .y, z = "FinnGen")) %>%
  bind_rows 

pool.dat <- map2(pool.files, pool.traits, ~extr.out.dat(.x, .y, z = "Pooled")) %>%
  bind_rows 


############################################################################
#               Combine data and harmonise
############################################################################

##
## Combine outcome data

outdat <- bind_rows(gwas.dat, ukb.dat, finngen.dat, pool.dat)


##
## Harmonise SNP-exposure and SNP-outcome data

uvmr_dat <- filter(expdat, exclfads == "no") %>%
  harmonise_data(., outdat, action = 2) 


############################################################################
#                     Multivariable regression estimates
############################################################################

##
## Import estimates from multivariable regression

mv_res <- read.csv("output/MV_nmr150k.csv") %>%
  mutate(exposure = sub("_int", "", exposure)) 


############################################################################
#                     Univariable MR estimates
############################################################################

##
## Outcome info file (links id to outcome name)

out.info <- outdat %>%
  select(id.outcome, trait, type_outcome, order, source) %>%
  distinct


##
## Run UVMR

## UVMR estimates
eff <- mr(uvmr_dat, method_list = c("mr_ivw", "mr_egger_regression")) %>%
  mutate(id.model = paste0(id.exposure, "_X_", id.outcome)) 

## Between SNP heterogeneity			
het <- mr_heterogeneity(uvmr_dat) %>%
  mutate(id.model = paste0(id.exposure, "_X_", id.outcome)) %>%
  select(id.model, method, Q, Q_df, Q_pval)

## MR-Egger intercept
int <- mr_pleiotropy_test(uvmr_dat) %>%
  mutate(id.model = paste0(id.exposure, "_X_", id.outcome)) %>%
  rename(egger_intercept_SE = se, egger_intercept_pval = pval) %>%
  select(id.model, egger_intercept, egger_intercept_SE, egger_intercept_pval )

## Combine results from UVMR				
uvmr_res <- merge(eff, het, by = c("id.model", "method"), all.x = T) %>%
  merge(., int, by = "id.model") %>%
  merge(., out.info, by = "id.outcome") %>%
  mutate(outcome = trait, type_model = "UVMR") %>%
  select(id.model, id.exposure, id.outcome, exposure, outcome, order, type_outcome, source, type_model, method, 
         nsnp, b, se, pval, Q, Q_df, Q_pval, egger_intercept, egger_intercept_SE, egger_intercept_pval)


##
## Leave-one-out analysis

uvmr_dat2 <- mutate(uvmr_dat, exposure = paste0(exposure, " (UKBB)"), outcome = paste0(trait, " (", source, ")"))
uvmr_res2 <- mutate(uvmr_res, exposure = paste0(exposure, " (UKBB)"), outcome = paste0(outcome, " (", source, ")"))

loo <- mr_leaveoneout(uvmr_dat2) 
pdf(file="./output/graphs/loo_MR_nmr150k.pdf")
mr_leaveoneout_plot(loo)
dev.off()


############################################################################
#                     Multivariable MR estimates
############################################################################

##
## Exposures for MVMR

mv.exp <- unique(expdat$exposure)


##
## Covariables for MVMR

mv.cov <- c("Total_FA", "Total_TG", "Clinical_LDL_C", "ApoB")


##
## Select IVs for covariables

covdat <- map(mv.cov, ~fread(paste0(dat.dir, "nmr_ukb_150k/", ., "_int.imputed.txt.gz"))) %>%
  map2(., mv.cov, ~mutate(.x, exposure = sub("_int[.]imputed[.]txt[.]gz$", "", .y))) %>%
  map(., ~mutate_at(., vars(P_BOLT_LMM_INF), ~as.numeric(.))) %>%
  map(., ~filter(., INFO > 0.8 
                 & 
                   P_BOLT_LMM_INF < 5e-08
                 & 
                   A1FREQ < 0.99
                 & 
                   A1FREQ > 0.01							  
  )) %>%
  map(., ~format_data(., 
                      type = "exposure",
                      phenotype_col = "exposure",
                      snp_col = "SNP",
                      beta_col = "BETA",
                      se_col = "SE",
                      eaf_col = "A1FREQ",
                      effect_allele_col = "ALLELE1",
                      other_allele_col = "ALLELE0",
                      pval_col = "P_BOLT_LMM_INF",
                      info_col = "INFO",
                      chr_col = "CHR",
                      pos_col = "BP"
  )) %>%
  map(., ~clump_data(., clump_kb = 10000, clump_r2 = 0.001)) %>%
  bind_rows 


##
## Select unique IVs for exposures/covariables

exp.snps <- filter(expdat, exposure %in% mv.exp) %>% pull(SNP) %>% unique
cov.snps <- filter(covdat, exposure %in% mv.cov) %>% pull(SNP) %>% unique
mv.snps <- c(exp.snps, cov.snps) %>% unique


##
## Extract SNP-exposure and SNP-covariable data for all MVMR SNPs

mv.expdat <- c(mv.exp, mv.cov) %>%
  map(., ~paste0(., "_int.imputed.txt.gz")) %>%
  map(., ~fread(paste0(dat.dir, "nmr_ukb_150k/", .))) %>%
  map(., ~format_data(., 
                      type = "exposure",
                      snps = mv.snps,
                      phenotype_col = "exposure",
                      snp_col = "SNP",
                      beta_col = "BETA",
                      se_col = "SE",
                      eaf_col = "A1FREQ",
                      effect_allele_col = "ALLELE1",
                      other_allele_col = "ALLELE0",
                      pval_col = "P_BOLT_LMM_INF",
                      info_col = "INFO",
                      chr_col = "CHR",
                      pos_col = "BP"
  )) %>%
  map2(., c(mv.exp, mv.cov), ~mutate(.x, exposure = .y)) %>%
  bind_rows 

##
## Extract SNP-outcome data for all MVMR SNPs

extr.mv.outdat <- function(f, y, z) {
  
  df <- fread(f) %>%
    format_data(.,
                snps = mv.snps,
                type = "outcome", 
                phenotype_col = "trait")
  
  info <- filter(cvd.info, abbrev == y & source %in% z)
  
  df <- merge(df, info, by.x="outcome", by.y="abbrev")
  
}

mv.gwasdat <- map2(gwas.files, gwas.traits, ~extr.mv.outdat(.x, .y, z = c("GWAS", "GWAS+UKB"))) %>%
  bind_rows 

mv.ukbdat <- map2(ukb.files, ukb.traits, ~extr.mv.outdat(.x, .y, z = "UKB")) %>%
  bind_rows 

mv.finngendat <- map2(finngen.files, finngen.traits, ~extr.mv.outdat(.x, .y, z = "FinnGen")) %>%
  bind_rows 

mv.pooldat <- map2(pool.files, pool.traits, ~extr.mv.outdat(.x, .y, z = "Pooled")) %>%
  bind_rows 

mv.outdat <- bind_rows(mv.gwasdat, mv.ukbdat, mv.finngendat, mv.pooldat)


##
## Outcome linking file 

mvout.info <- mv.outdat %>%
  select(id.outcome, trait, type_outcome, order, source) %>%
  distinct


##		
## Harmonise all to be on the same effect allele

mv.dat <- harmonise_data(mv.expdat, mv.outdat, action = 2) 


##
## Import phenotypic covariance matrix for NMR data

pcov <- read.table("./data/pcov_nmr_dat_all.txt", row.names = 1) %>% as.matrix	


##
## Function to format and run MVMR
run_mvmr <- function(exp1, exp2, exp3 = NULL, idout) {
  
  # Print model
  print(paste("Doing MVMR:", idout, "~", exp1, "+", exp2, "+", exp3))
  
  ################ SNP selection ######################
  
  # Select SNPs associated with at least one exposure at GWAS threshold (P < 0.0000005)							
  if(is.null(exp3)) {
    gw.snps <- filter(mv.dat, (exposure == exp1 |  exposure == exp2) & id.outcome == idout) %>% 
      filter(pval.exposure < 5e-08) %>%
      pull(SNP) %>%
      unique %>%
      as.character
  } else {			
    gw.snps <- filter(mv.dat, (exposure == exp1 |  exposure == exp2 |  exposure == exp3) & id.outcome == idout) %>% 
      filter(pval.exposure < 5e-08) %>%
      pull(SNP) %>%
      unique %>%
      as.character
  }
  
  print(paste("N SNPs associated with at least one exposure:", length(gw.snps)))
  
  # Clump these to avoid variants in LD between exposures (use p-value from first exposure)
  nold.snps <- filter(mv.dat, exposure == exp1 & id.outcome == idout) %>%
    filter(SNP %in% gw.snps) %>% 
    clump_data %>% 
    pull(SNP) %>%
    as.character
  
  print(paste("N SNPs not in LD associated with at least one exposure:", length(nold.snps)))
  
  ################ Format data ######################
  
  # Select dataset for index exposure
  exp1.dat <- filter(mv.dat, exposure == exp1 & id.outcome == idout & SNP %in% nold.snps) %>% 
    select(SNP, beta.exposure, se.exposure)
  
  # Select dataset for adjusted exposure
  exp2.dat <- filter(mv.dat, exposure == exp2 & id.outcome == idout & SNP %in% nold.snps) %>% 
    select(SNP, beta.exposure, se.exposure)
  
  # Select dataset for second adjusted exposure
  if(is.null(exp3)==F) {
    exp3.dat <- filter(mv.dat, exposure == exp3 & id.outcome == idout & SNP %in% nold.snps) %>% 
      select(SNP, beta.exposure, se.exposure)
  }
  
  # Select dataset for one outcome (keeping only independent SNPs))
  out1.dat <- filter(mv.dat, exposure == exp1 & id.outcome == idout & SNP %in% nold.snps) %>% 
    select(SNP, beta.outcome, se.outcome)
  
  # Merge datasets
  if(is.null(exp3)) {
    mv.raw <- merge(exp1.dat, exp2.dat, by = "SNP") %>%
      merge(., out1.dat, by = "SNP")
  } else {
    mv.raw <- merge(exp1.dat, exp2.dat, by = "SNP") %>%
      merge(., out1.dat, by = "SNP") %>%
      merge(., exp3.dat, by = "SNP")
  }
  
  # Format data for MVMR
  if(is.null(exp3)) {
    mv.input <- MVMR::format_mvmr(
      BXGs = mv.raw[,c(2,4)],
      BYG = mv.raw[,6],
      seBXGs = mv.raw[,c(3,5)],
      seBYG = mv.raw[,7],
      RSID = mv.raw[,1])
  } else {		
    mv.input <- MVMR::format_mvmr(
      BXGs = mv.raw[,c(2,4,8)],
      BYG = mv.raw[,6],
      seBXGs = mv.raw[,c(3,5,9)],
      seBYG = mv.raw[,7],
      RSID = mv.raw[,1])
  }
  
  # Number of SNPs included
  nsnps <- length(unique(mv.input$SNP))
  
  ################ Approximate covariance matrix for SNP-exp1 and SNP-exp2  ######################
  
  # Select variance-covariance matrix for exp1 and exp2
  if(is.null(exp3)) {
    pcov_sel <- pcov[c(exp1, exp2), c(exp1, exp2)]	
  } else {
    pcov_sel <- pcov[c(exp1, exp2, exp3), c(exp1, exp2, exp3)]	
  }
  print("Pheno covar matrix: ")
  print(pcov_sel)
  
  # Convert to covariance matrices for SNP-exposure effects 
  if(is.null(exp3)) {
    gcov <- MVMR::phenocov_mvmr(pcov_sel, mv.input[,c("sebetaX1", "sebetaX2")])	
  } else {
    gcov <- MVMR::phenocov_mvmr(pcov_sel, mv.input[,c("sebetaX1", "sebetaX2", "sebetaX3")])	
  }
  
  ################ Run MVMR  ######################
  
  # Test for weak instruments 
  if(is.null(exp3)) {
    mv.weakiv <- MVMR::strength_mvmr(r_input = mv.input, gencov = gcov) %>%
      as.data.table(.) %>%
      rename(cF_exp1 = exposure1, cF_exp2 = exposure2) %>%
      mutate(model = paste0(exp1, "_x_", exp2, "_x_", exp3, "_x_", idout))	
  } else {
    mv.weakiv <- MVMR::strength_mvmr(r_input = mv.input, gencov = gcov) %>%
      as.data.table(.) %>%
      rename(cF_exp1 = exposure1, cF_exp2 = exposure2, cF_exp3 = exposure3) %>%
      mutate(model = paste0(exp1, "_x_", exp2, "_x_", exp3, "_x_", idout))
  }
  
  
  # Test for horizontal pleiotropy 
  mv.pleio <- MVMR::pleiotropy_mvmr(r_input = mv.input, gencov = gcov) %>%
    as.data.table(.) %>%
    mutate(model = paste0(exp1, "_x_", exp2, "_x_", exp3, "_x_", idout)) 
  
  # Run MVMR to estimate direct effects
  mv.effect <- MVMR::ivw_mvmr(r_input = mv.input) %>%
    as.data.table(.) %>%
    mutate(type_exp = paste0("exp", row.names(.)),
           model = paste0(exp1, "_x_", exp2, "_x_", exp3, "_x_", idout)
    ) 
  
  if(is.null(exp3)) {
    mv.effect <- mutate(mv.effect, exp = case_when(type_exp == "exp1" ~ exp1, type_exp == "exp2" ~ exp2))
  } else {
    mv.effect <- mutate(mv.effect, exp = case_when(type_exp == "exp1" ~ exp1, type_exp == "exp2" ~ exp2, type_exp == "exp3" ~ exp3))
  }											
  
  mv.effect <- rename(mv.effect, SE = "Std. Error", t = "t value", pval = "Pr(>|t|)") %>%								
    dcast(., model ~ type_exp, value.var = c("exp", "Estimate", "SE", "t", "pval"))											
  
  # Combine all results
  mv.res <- merge(mv.weakiv, mv.pleio, by = "model") %>% 
    merge(., mv.effect, by = "model") %>%
    mutate(nSNPs = nsnps, 
           type_model = "MVMR",
           id.outcome = idout
    ) %>%
    merge(., mvout.info, by = "id.outcome")
  
  return(mv.res)
  
}

# List of arguments for MVMR
mvmr_args_2exp <- cross_df(list(exp1 = mv.exp, exp2 = mv.cov, idout = unique(mvout.info$id.outcome)))
mvmr_args_3exp <- cross_df(list(exp1 = mv.exp, exp2 = "Total_TG", exp3 = "Clinical_LDL_C", idout = unique(mvout.info$id.outcome)))

# Apply MVMR

mvmr.out_2exp <- pmap(mvmr_args_2exp, run_mvmr) %>% 
  bind_rows %>%				
  mutate(method = paste0(type_model, " (+", exp_exp2, ")"), exp_source = "UKB")

mvmr.out_3exp <- pmap(mvmr_args_3exp, run_mvmr) %>% 
  bind_rows %>%
  mutate(method = paste0(type_model, " (+", exp_exp2, "+", exp_exp3, ")"), exp_source = "UKB")

mvmr.out <- bind_rows(mvmr.out_2exp, mvmr.out_3exp)

# Save all results
writexl::write_xlsx(list(
  sum_dat          = uvmr_dat,
  uvmr_res         = uvmr_res,
  loo              = loo, 
  MVMR             = mvmr.out
),
path = "./output/MR_nmr150k.xlsx")

# Combine UVMR + MVMR results
all.res <- mvmr.out %>%
  rename(exposure = exp_exp1,
         b = Estimate_exp1,
         se = SE_exp1,
         pval = pval_exp1,
         id.model = model,
         outcome = trait,
         F = cF_exp1,
         nsnp = nSNPs
  ) %>%
  select(id.model, exposure, id.outcome, outcome, type_outcome, order, source, method, type_model, exp_source, b, se, pval, F, nsnp) %>%
  bind_rows(., uvmr_res) %>%
  bind_rows(., mv_res)


############################################################################
#                     Make plot comparing methods (study specific)
############################################################################

##
## Subsets data for plots - ALL data sources

make_plot1 <- function(expname, typeout) {
  
  # Subset data for plot
  df <- all.res %>%
    filter(method == "Inverse variance weighted" & exposure %in% expname & type_outcome %in% typeout) %>%
    arrange(order) %>%
    mutate(outcome = factor(outcome, levels = rev(factor(unique(outcome)))),
           source = factor(source, levels = rev(c("GWAS", "UKB", "GWAS+UKB", "FinnGen", "Pooled", "MR-Base")))
    )  
  
  # Generate plot (1 FA, 1 outcome type, fads y/n)
  p <- ggforestplot::forestplot(
    df = df,
    name = outcome,
    estimate = b,
    se = se,
    pvalue = pval,
    psignif = 0.05/9,
    xlab = "OR (95% CI)",
    title = expname,
    colour = source,
    shape = source,
    logodds = T	
  ) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))
  
  # Save plot
  ggsave(paste0("./output/graphs/forest_MR_NMR150k_", "UKB", "-SNPs_", expname, "_vs_", typeout, ".png"), plot = p, width = 10, height = 10)									
}						


# List of arguments for MVMR
plot_args <- cross_df(list(expname = mv.exp, typeout = as.character(unique(all.res$type_outcome))))

# Generate and save plots
pwalk(plot_args, make_plot1) 


############################################################################
#                     Make plot comparing methods (pooled data)
############################################################################	

##
## Subsets data for plots - *** POOLED data only ***

make_plot2 <- function(expname, typeout) {
  
  # Subset data for plot
  df <- all.res %>%
    filter(! method == "MV (crude)") %>%
    filter( (source == "Pooled") | (source == "UKB" & outcome == "Aortic valve stenosis") | (source == "UKB" & method == "MV (adjusted)") ) %>%
    arrange(order) %>%
    mutate(
      outcome = factor(outcome, levels = rev(factor(unique(outcome)))),
      method = factor(method, levels = rev(c("MV (adjusted)", "Inverse variance weighted", "MR Egger",
                                             "MVMR (+Total_FA)", "MVMR (+Total_TG)", "MVMR (+Clinical_LDL_C)", 
                                             "MVMR (+ApoB)", "MVMR (+Total_TG+Clinical_LDL_C)")))
    ) %>%
    filter(exposure %in% expname & type_outcome %in% typeout) 	
  
  # Generate plot (1 FA, 1 outcome type)
  p <- ggforestplot::forestplot(
    df = df,
    name = outcome,
    estimate = b,
    se = se,
    pvalue = pval,
    psignif = 0.05/9,
    xlab = "OR (95% CI)",
    title = expname,
    colour = method,
    #shape = source,
    logodds = T	
  ) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18)) +
    scale_x_continuous(limits = c(0.2, 2.5), breaks = c(0.5, 1, 1.5, 2, 2.5))
  
  # Save plot
  ggsave(paste0("./output/graphs/forest_MR_NMR150k_", "UKB", "-SNPs_", expname, "_vs_", typeout, "_POOLEDonly", ".png"), plot = p, width = 10, height = 10)
}					

# Generate and save plots
pwalk(plot_args, make_plot2) 

q("no")