Jayachandran_2021_tenofovir_p24 <- function() {
  description <- paste(
    "Mechanistic ex vivo HIV-1 viral-dynamics PK/PD model linking rectal mucosal",
    "mononuclear cell (MMC) tenofovir-diphosphate (TFVdp) concentration to",
    "suppression of cumulative p24 antigen expression in an ex vivo rectal-tissue",
    "explant challenge, from the RMP-02/MTN-006 pre-exposure-prophylaxis trial.",
    "An exponential HIV-1 growth compartment with a first-order death rate feeds a",
    "delayed p24 antigen expression compartment; the drug inhibits viral growth",
    "through a linear effect E = slope * CP, where CP is the MMC TFVdp",
    "concentration measured 30 minutes postdose and degrading exponentially over",
    "the ex vivo assay time course. The MMC TFVdp concentration is supplied as an",
    "exogenous covariate (CONC_TFVDP_FMOLMC) from the companion multicompartment",
    "PK model Jayachandran_2021_tenofovir, whose parameters were fixed when this",
    "PK/PD layer was fitted. Parameters are for the total-MMC cell type, the",
    "cell type the authors used for their target-concentration simulations.",
    sep = " "
  )
  reference <- paste(
    "Jayachandran P, Garcia-Cremades M, Vucicevic K, Bumpus NN, Anton P,",
    "Hendrix C, Savic R. A Mechanistic In Vivo/Ex Vivo",
    "Pharmacokinetic-Pharmacodynamic Model of Tenofovir for HIV Prevention.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10(3):179-187.",
    "doi:10.1002/psp4.12583. Viral-dynamics equations (Eqs 5-6) from Methods;",
    "final total-MMC parameter estimates from Table 2; drug-effect and",
    "drug-degradation forms from Methods and the Supplementary Model Code 2",
    "NONMEM control stream. MMC TFVdp concentrations come from the companion PK",
    "model; see modellib('Jayachandran_2021_tenofovir').",
    sep = " "
  )
  vignette <- "Jayachandran_2021_tenofovir"
  units <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  covariateData <- list(
    CONC_TFVDP_FMOLMC = list(
      description = "MMC tenofovir-diphosphate concentration measured 30 minutes postdose, driving the ex vivo drug effect",
      units = "fmol/million cells",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Static per-explant exogenous concentration supplied from the companion",
        "multicompartment PK model Jayachandran_2021_tenofovir (MMC TFVdp biophase",
        "state). Enters the linear drug effect E = slope * CONC_TFVDP_FMOLMC *",
        "exp(-kdeg * t), where t is the ex vivo assay time. Set to 0 for a",
        "baseline (no-treatment) explant, which makes E = 0 and reproduces the",
        "untreated viral-growth trajectory; the source Supplementary Model Code 2",
        "gates the effect with IF(VISIT.NE.2), i.e. zero at the baseline visit.",
        "In the source control stream the column is CP. Sweep it across 0 to",
        "11,000 fmol/million cells to reproduce the target-suppression simulation",
        "(Figure 5); a concentration of about 9,000 fmol/million cells suppresses",
        "apparent viral replication to 1%.",
        sep = " "
      ),
      source_name = "CP (MMC TFVdp concentration, fmol/million cells)"
    )
  )

  # virus is the canonical free-virion pool; p24 (cumulative p24 antigen
  # expression, the delayed effect compartment) has no canonical and is declared
  # paper-specific. Both states are latent PD quantities with no biological
  # specimen matrix.
  paper_specific_compartments <- c("p24")

  compartmentData <- list(
    virus = list(analyte = "HIV-1 virions", units = "virions/mL", specimen = "not applicable", verified = TRUE),
    p24 = list(analyte = "HIV-1 p24 capsid antigen", units = "pg/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human (ex vivo rectal-tissue explant challenge)",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "22-66 years",
    sex_female_pct = 22.2,
    disease_state = "HIV-1 seronegative healthy adults; rectal-tissue explants challenged ex vivo with HIV-1",
    dose_range = "Single oral 300 mg TDF; single and multiple (7 daily) rectal 1% TFV gel (see the companion PK model)",
    regions = "USA (Los Angeles, CA and Pittsburgh, PA)",
    notes = paste(
      "RMP-02/MTN-006 (NCT00984971). Rectal-tissue explants were collected at",
      "baseline and postdose, infected ex vivo with HIV-1 within 1-2 hours of",
      "harvest, and supernatant p24 antigen quantified on assay days 0/1, 4, 7,",
      "11 and 14. Measurements across four biopsies were averaged and accumulated",
      "to give a cumulative p24 antigen level at each post-challenge timepoint",
      "(682 cumulative p24 observations, 24% below the 10 pg/mL limit of",
      "quantification, Table S1). The model was fitted for CD4-, CD4+ and total",
      "MMC cell types; only slope and the residual error differed between them",
      "(Table 2), and the total-cell-type parameters shipped here are those used",
      "for the target-concentration simulations.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # HIV-1 viral dynamics (Methods Eqs 5-6, total-MMC cell type Table 2).
    #   d virus / dt = kgrow * virus * (1 - E) - kdeath * virus
    #   d p24   / dt = ke0  * (ppc * virus - p24)
    # The p24 compartment is a delayed effect compartment (Methods); its
    # first-order expression rate maps to the canonical effect-compartment
    # equilibration rate ke0 (source k_p), and the virus-to-p24 conversion
    # ratio to the pseudo-partition coefficient ppc (source R). Cell growth
    # and p24 expression rates were fixed to literature values, the death
    # rate and ratio fixed from the baseline (no-treatment) model.
    # ------------------------------------------------------------------
    lkgrow <- fixed(log(0.0320)); label("HIV-1 viral growth rate (1/h)")                       # Table 2 'k growth' = 0.0320 /h (fixed to literature value)
    lkdeath <- fixed(log(0.0193)); label("HIV-1 viral death rate (1/h)")                       # Table 2 'k death' = 0.0193 /h (fixed from baseline model)
    ke0 <- fixed(log(0.00400)); label("p24 antigen expression rate (1/h)")                     # Table 2 'k p24' = 0.00400 /h (fixed to literature value); log scale
    ppc <- fixed(log(0.0404)); label("Virus-to-p24 antigen conversion ratio (pg/virions)")     # Table 2 'Ratio' = 0.0404 pg/virions (fixed from baseline model); log scale
    lkdeg <- fixed(log(0.0018)); label("Ex vivo drug degradation rate (1/h)")                  # Table 2 'k degradation' = 0.0018 /h (fixed; TFV degradation kinetics, ref 29)
    lslope <- log(0.00011); label("Linear drug-effect slope ((pg/mL)/(fmol/million cells)), total MMC") # Table 2 TOTAL 'Slope' = 0.00011 (RSE 11%)

    # Between-subject variability on the p24 expression rate only; the
    # source $OMEGA gives the variance (0.525) directly and the Table 2
    # %CV of 72.5% is sqrt(0.525). Fixed for modelling the treatment effect.
    etake0 ~ fixed(0.525)  # variance 0.525; sqrt(0.525) = 0.725 = Table 2 'IIV k p24 72.5%'

    # Combined additive + proportional residual error (Table 2, TOTAL).
    propSd <- 0.876; label("Proportional residual error, cumulative p24 (fraction)")  # Table 2 TOTAL 'Proportional error' 87.6% CV (RSE 6%)
    addSd <- 3.61; label("Additive residual error, cumulative p24 (pg/mL)")           # Table 2 TOTAL 'Additive error' 3.61 pg/mL (RSE 11%)
  })

  model({
    kgrow <- exp(lkgrow)
    kdeath <- exp(lkdeath)
    kp24 <- exp(ke0 + etake0)
    ratio <- exp(ppc)
    kdeg <- exp(lkdeg)
    slope <- exp(lslope)

    # Linear drug effect. The MMC TFVdp concentration is measured 30 minutes
    # postdose and degrades exponentially over the ex vivo assay time t; a
    # baseline explant carries CONC_TFVDP_FMOLMC = 0, giving E = 0.
    edrug <- slope * CONC_TFVDP_FMOLMC * exp(-kdeg * t)

    # Initial conditions: 10^4 virions, no p24 (Methods).
    virus(0) <- 10000
    p24(0) <- 0

    d/dt(virus) <- kgrow * virus * (1 - edrug) - kdeath * virus
    d/dt(p24) <- kp24 * (ratio * virus - p24)

    Cp24 <- p24
    Cp24 ~ add(addSd) + prop(propSd)
  })
}
