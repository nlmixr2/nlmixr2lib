Yamazaki_2008_crizotinib_mouse <- function() {
  description <- paste(
    "Preclinical (athymic mouse; GTL16 gastric carcinoma or U87MG glioblastoma",
    "xenograft). Integrated PK + cMet phosphorylation (effect-compartment link",
    "model) + exponential tumor-growth-inhibition (TGI) model for orally",
    "administered crizotinib (PF02341066), an ATP-competitive cMet receptor",
    "tyrosine kinase inhibitor. PK is one-compartment first-order absorption",
    "with a fixed 0.8 h lag, fitted by naive-pooled analysis (dose-group-specific",
    "estimates due to nonlinear kinetics; the encoded set is Study 2 at 50 mg/kg).",
    "The cMet phosphorylation response is the Sheiner 1979 link model with E0,",
    "Emax, and Hill coefficient all fixed at 1 (Imax 1/(1 + Ce/EC50) form). The",
    "tumor-growth model is exponential, with the growth rate inhibited by the",
    "plasma concentration via a sigmoidal Imax 1/(1 + Cc/EC50_tumor) function",
    "(Emax fixed at 1; the saturable tumor-volume capacity term TG50 was",
    "rejected by the authors as TG50 >> Tmax). Default TGI parameters reproduce",
    "the GTL16 fit; the U87MG variant (kin_tumor=0.0134, kout_tumor=0.00236,",
    "EC50_tumor=94.1 ng/mL) is documented in population$notes and demonstrated",
    "in the validation vignette."
  )
  reference <- paste(
    "Yamazaki S, Skaptason J, Romero D, Lee JH, Zou HY, Christensen JG,",
    "Koup JR, Smith BJ, Koudriakova T. Pharmacokinetic-pharmacodynamic",
    "modeling of biomarker response and tumor growth inhibition to an",
    "orally available cMet kinase inhibitor in human tumor xenograft",
    "mouse models. Drug Metab Dispos. 2008;36(7):1267-1274.",
    "doi:10.1124/dmd.107.019711. PMID 18381487."
  )
  vignette <- "Yamazaki_2008_crizotinib_mouse"
  units <- list(time = "hour", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species       = "mouse (athymic nude; GTL16 human gastric carcinoma or U87MG human glioblastoma subcutaneous xenograft)",
    n_subjects    = 200L,
    n_studies     = 3L,
    sex_female_pct = NA_real_,
    disease_state = "GTL16 gastric-carcinoma or U87MG glioblastoma xenograft tumors",
    dose_range    = "PF02341066 (crizotinib) 3.13-50 mg/kg PO once daily for 9-11 days; plasma+cMet phosphorylation sampled at 1, 4, 8, 24 h after last dose (n=3/timepoint)",
    regions       = "Preclinical (Pfizer Global Research and Development, La Jolla, CA)",
    notes         = paste(
      "Three repeated-dose studies: (1) GTL16 xenograft at 8.5, 17, 34 mg/kg;",
      "(2) GTL16 xenograft at 3.13, 6.25, 12.5, 25, 50 mg/kg; (3) U87MG xenograft",
      "at 3.13, 6.25, 12.5, 25, 50 mg/kg. Subjects n=3/timepoint were humanely",
      "euthanized for plasma and tumor PD sampling; tumor volume was measured",
      "longitudinally during the treatment period. PK was estimated by",
      "naive-pooled analysis with dose-group-specific parameters because of",
      "nonlinear kinetics (saturation of hepatic/intestinal clearance at higher",
      "doses; CL/F and V/F tended higher at lower doses, Table 1 footnote and",
      "Discussion). The PK values encoded here are Study 2 (GTL16) at 50 mg/kg",
      "(post-saturation, linear-PK regime, full-efficacy dose): ka=0.331 1/h,",
      "CL/F=1.80 L/h/kg, V/F=5.56 L/kg. The full 11-row Table 1 of per-dose-group",
      "PK estimates is reproduced in the validation vignette. PD-link model",
      "parameters (Table 2 Link Model column) come from a joint fit across the",
      "GTL16 cohort only: ke0=0.135 1/h, EC50_phospho=18.5 ng/mL. TGI exponential",
      "growth parameters (Table 3) were fit separately to each cell line:",
      "GTL16 (encoded default) kin_tumor=0.0130 1/h, kout_tumor=0.00672 1/h,",
      "EC50_tumor=213 ng/mL; U87MG kin_tumor=0.0134 1/h, kout_tumor=0.00236 1/h,",
      "EC50_tumor=94.1 ng/mL. The IDR-only and IDR-with-effect-compartment PD",
      "alternatives (Table 2) were superseded by the simpler link model in the",
      "paper's Discussion (combined model collapses to the link model because",
      "kout >> ke0, indicating near-instantaneous phospho equilibration)."
    )
  )

  ini({
    # ---- PK: 1-compartment, first-order absorption with fixed 0.8 h lag ----
    # Naive-pooled fit; dose-group-specific values (Table 1). Encoded set is
    # Study 2 (GTL16) at 50 mg/kg -- post-saturation, linear-PK regime,
    # full-efficacy dose. Full Table 1 is in the vignette Source-trace section.
    lka     <- log(0.331);  label("Absorption rate constant ka (1/h)")                 # Table 1 PKPD 2 at 50 mg/kg (SE 0.018)
    lcl     <- log(1.80);   label("Apparent oral clearance CL/F (L/h/kg)")             # Table 1 PKPD 2 at 50 mg/kg (SE 0.17)
    lvc     <- log(5.56);   label("Apparent oral volume of distribution V/F (L/kg)")   # Table 1 PKPD 2 at 50 mg/kg (SE 0.68)
    ltlag   <- fixed(log(0.8)); label("Absorption lag time tlag (h) -- fixed")         # Results paragraph 1 ("fixed absorption lag time of 0.8 h")

    # ---- PD link model: cMet phosphorylation inhibition (GTL16 cohort) ----
    # Equations 1 (effect compartment) and 2 (Imax response) of Methods.
    # E0, Emax, and Hill (gamma) are all fixed at 1 (Methods + Table 2).
    lke0          <- log(0.135); label("Effect-compartment equilibration rate ke0 (1/h)")     # Table 2 Link Model (SE 0.020)
    lec50_phospho <- log(18.5);  label("EC50 for cMet phosphorylation inhibition (ng/mL)")    # Table 2 Link Model (SE 2.65)

    # ---- TGI: exponential growth + Imax inhibition by plasma Cc (GTL16 default) ----
    # Simplified Equation 8 of Methods (TG50 term dropped because Discussion:
    # TG50 estimates exceeded the observed maximum tumor volume by >10000 mm^3
    # so the saturable capacity collapses to 1). Emax_tumor fixed at 1.
    lkin_tumor    <- log(0.0130);  label("Tumor first-order growth rate kin (1/h)")           # Table 3 GTL16 (SE 0.00214)
    lkout_tumor   <- log(0.00672); label("Tumor first-order loss rate kout (1/h)")            # Table 3 GTL16 (SE 0.00243)
    lec50_tumor   <- log(213);     label("EC50 for tumor growth inhibition (ng/mL)")          # Table 3 GTL16 (SE 123)

    # ---- Baseline tumor volume (operator-derived) ----
    # Methods: initial condition was the "individual tumor volume (cubic
    # millimeters)" per mouse; a typical-value population baseline is not
    # reported numerically. The 150 mm^3 below was read off Figures 5 (GTL16)
    # and 6 (U87MG) day-0 tumor volumes (~100-200 mm^3 across vehicle and
    # treated groups) for use as a simulation default; documented in the
    # validation vignette under Assumptions and deviations.
    lrbase_tumor  <- log(150); label("Typical baseline tumor volume V_T0 (mm^3) -- figure-derived")  # operator-extracted from Figures 5 and 6 (paper does not report a numerical typical baseline)

    # ---- Inter-individual variability ----
    # Paper Methods: "an interanimal variability for kout was estimated using
    # an exponential variance model" in the TGI exponential model. The
    # magnitude is NOT reported in Table 3 or anywhere else in the paper.
    # Encoded fixed at 0 per missing-value convention; downstream users
    # running stochastic VPCs must supply their own omega^2.
    etalkout_tumor ~ fixed(0)  # Methods: kout IIV reported as estimated, magnitude not given (see vignette Errata)

    # ---- Residual error ----
    # PK: Results section ("Residual variability was estimated to be 28%")
    # implies a proportional residual on Cc.
    propSd           <- 0.28;       label("Proportional residual error on plasma PF02341066 concentration (fraction)")  # Results paragraph 1 ("residual variability was estimated to be 28%")
    # PD (phospho) and TGI (tumor_vol): Methods states "Residual variability
    # was characterized by a proportional error model" but no numerical
    # magnitudes are reported in Tables 2 or 3. Encoded fixed(0) per
    # missing-value convention; documented in the vignette Errata.
    propSd_phospho   <- fixed(0);   label("Proportional residual error on cMet phosphorylation (fraction) -- magnitude not reported")    # Methods (proportional model used; magnitude not in Table 2)
    propSd_tumor_vol <- fixed(0);   label("Proportional residual error on tumor volume (fraction) -- magnitude not reported")            # Methods (proportional model used; magnitude not in Table 3)
  })

  model({
    # 1. Individual parameters
    ka            <- exp(lka)
    cl            <- exp(lcl)
    vc            <- exp(lvc)
    tlag          <- exp(ltlag)
    ke0           <- exp(lke0)
    ec50_phospho  <- exp(lec50_phospho)
    kin_tumor     <- exp(lkin_tumor)
    kout_tumor    <- exp(lkout_tumor + etalkout_tumor)
    ec50_tumor    <- exp(lec50_tumor)
    rbase_tumor   <- exp(lrbase_tumor)

    # 2. PK ODEs: one-compartment oral with first-order absorption
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    alag(depot)   <- tlag
    # central is mg/kg (dose amt enters as mg/kg); vc is L/kg; central/vc is
    # mg/L = ug/mL; the * 1000 converts to the ng/mL scale used by EC50_phospho
    # and EC50_tumor (Tables 2 and 3) and declared in units$concentration.
    Cc <- central / vc * 1000                                # ng/mL

    # 3. Effect compartment for the link model (Sheiner 1979).
    # dCe/dt = ke0 * (Cp - Ce); Ce drives the cMet phosphorylation inhibition
    # (Equation 1 of Methods). effect(0) defaults to 0 (no drug pre-dose).
    d/dt(effect) <- ke0 * (Cc - effect)
    Ce <- effect

    # 4. cMet phosphorylation inhibition (link model, Equation 2 of Methods).
    # E = E0 - Emax * Ce^gamma / (EC50^gamma + Ce^gamma), with E0 = Emax = gamma = 1
    # collapses to E = 1 - Ce / (EC50 + Ce). Baseline (Ce = 0) gives E = 1
    # (ratio to control animals); near-complete inhibition (Ce >> EC50) gives
    # E -> 0.
    phospho <- 1 - Ce / (ec50_phospho + Ce)

    # 5. Tumor-growth ODE (simplified Equation 8 of Methods, with Imax = 1
    # and the TG50 saturable-capacity term dropped per the Results /
    # Discussion simplification). Drug effect is driven by plasma Cc, not
    # the effect-site Ce (Equation 8 references the plasma concentration of
    # PF02341066 directly).
    tgi_inhib <- 1 - Cc / (ec50_tumor + Cc)
    d/dt(tumor_vol) <- kin_tumor * tgi_inhib * tumor_vol - kout_tumor * tumor_vol
    tumor_vol(0)    <- rbase_tumor                            # mm^3 starting tumor volume (typical value; per-subject override via initial-condition cmt event)

    # 6. Observations (multi-output: three independent endpoints).
    # Downstream fitting with a DVID column requires the data to encode the
    # DVID values matching the order of these residual statements (1=Cc,
    # 2=phospho, 3=tumor_vol). Simulation via rxSolve produces all three
    # columns regardless of DVID.
    Cc        ~ prop(propSd)
    phospho   ~ prop(propSd_phospho)
    tumor_vol ~ prop(propSd_tumor_vol)
  })
}
attr(Yamazaki_2008_crizotinib_mouse, "message") <-
  "Default TGI parameters reproduce the GTL16 xenograft fit (Table 3). To switch to the U87MG xenograft variant, override the TGI parameters: kin_tumor=0.0134 1/h, kout_tumor=0.00236 1/h, EC50_tumor=94.1 ng/mL (Table 3 U87MG column). The validation vignette shows how to do this with model$ini updates."
Yamazaki_2008_crizotinib_mouse
