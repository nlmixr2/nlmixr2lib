Jiao_2009_sirolimus <- function() {
  description <- "One-compartment population PK model for oral sirolimus in Chinese adult de novo renal transplant recipients on triple immunosuppression with ciclosporin and corticosteroids (Jiao 2009). First-order absorption with ka fixed at the literature value 0.752 1/h. Covariate effects on apparent clearance: linear-deviation effects of total cholesterol and whole-blood ciclosporin trough concentration centred on the cohort medians, multiplicative power-form effects of concomitant silymarin and glycyrrhizin co-therapy in hepatically impaired patients, and a power-form effect of the current sirolimus daily dose centred at 2 mg. Apparent volume of distribution carries a linear-deviation effect of ciclosporin trough concentration."
  reference <- paste(
    "Jiao Z, Shi XJ, Li ZD, Zhong MK.",
    "Population pharmacokinetics of sirolimus in de novo Chinese adult",
    "renal transplant patients.",
    "Br J Clin Pharmacol. 2009;68(1):47-54.",
    "doi:10.1111/j.1365-2125.2009.03392.x.",
    sep = " "
  )
  vignette <- "Jiao_2009_sirolimus"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    TCHOL = list(
      description        = "Fasting total cholesterol (whole-blood lipid panel).",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on apparent sirolimus CL/F centred at the cohort median 5.66 mmol/L: typical CL/F decreases by 0.662 L/h per 1 mmol/L increase in TCHOL above the median. Cohort range 0.82-10.9 mmol/L per Jiao 2009 Table 1. Per the Discussion: higher cholesterol implies more lipoprotein-binding capacity and a lower unbound fraction available for hepatic metabolism, reducing apparent CL/F.",
      source_name        = "TC"
    ),
    CP_CSA_NGML = list(
      description        = "Whole-blood ciclosporin (CsA) trough concentration (predose, end of dosing interval).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-record whole-blood CsA C0 supplied as a time-varying perpetrator-drug covariate. Linear-deviation effects on both CL/F (slope -0.00417 L/h per ng/mL) and V/F (slope -7.27 L per ng/mL) centred at the cohort median 104 ng/mL. Cohort range 0-508 ng/mL per Jiao 2009 Table 1; reflects therapeutic drug monitoring of CsA across the CsA-dose-reduction (Phase I) and CsA-discontinuation (Phase II) arms. Captures the well-documented in-vivo CYP3A4 / P-gp inhibition of sirolimus by ciclosporin.",
      source_name        = "C0"
    ),
    DOSE = list(
      description        = "Current sirolimus daily dose (DDS) supplied as a per-record covariate.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on CL/F: cl_eff *= (DOSE / 2)^0.479 with reference DDS = 2 mg/day (cohort median per Jiao 2009 Table 1). Cohort range 0.5-6.0 mg/day. The dose-dependence reflects a combination of low oral bioavailability (~17% sirolimus) and possible therapeutic-drug-monitoring-induced bias (subjects with higher CL/F titrate to higher doses), as discussed in the source. Time-varying per subject because the maintenance dose was titrated to target trough levels.",
      source_name        = "DDS"
    ),
    CONMED_SILYMARIN = list(
      description        = "Concomitant silymarin (milk thistle, Silybum marianum) co-therapy at the PK observation.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant silymarin)",
      notes              = "1 = subject is receiving silymarin co-therapy at the PK record, 0 = no concomitant silymarin. In Jiao 2009 silymarin was prescribed in patients with abnormal liver-function indices (all 15 silymarin-treated subjects had elevated ALT), so the modelled effect reflects combined hepatic impairment + silymarin CYP3A4 / P-gp inhibition. Multiplicative power-form effect on CL/F: cl_eff *= 0.65^CONMED_SILYMARIN, i.e. ~35% reduction in apparent sirolimus CL/F.",
      source_name        = "SLM"
    ),
    CONMED_GLYCYRRHIZIN = list(
      description        = "Concomitant glycyrrhizin (licorice, Glycyrrhiza uralensis) co-therapy at the PK observation.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant glycyrrhizin)",
      notes              = "1 = subject is receiving glycyrrhizin co-therapy at the PK record, 0 = no concomitant glycyrrhizin. Like silymarin, glycyrrhizin in this cohort (19 subjects) was prescribed for hepatic protection in patients with abnormal liver-function indices, so the modelled effect reflects combined hepatic impairment + glycyrrhizin CYP3A4 inhibition. Multiplicative power-form effect on CL/F: cl_eff *= 0.661^CONMED_GLYCYRRHIZIN, i.e. ~34% reduction in apparent sirolimus CL/F.",
      source_name        = "GLZ"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 112L,
    n_observations   = 804L,
    n_studies        = 1L,
    n_centres        = 5L,
    age_range        = "19-68 years",
    age_median       = "42 years (mean 42 +/- 9.9)",
    weight_range     = "44.5-91.0 kg",
    weight_median    = "60.0 kg (mean 60.4 +/- 9.43)",
    sex_female_pct   = 30.4,
    race_ethnicity   = "Chinese (single-country multicentre cohort across Guangzhou, Shanghai and Beijing).",
    disease_state    = "Chinese adult de novo primary renal allograft recipients on triple immunosuppression. Patients with active major infection, hepatitis B / C, HIV, decreased platelets, severe dyslipidaemia (fasting triglycerides >= 4.6 mmol/L or fasting total cholesterol >= 7.8 mmol/L) or multiple organ transplants were excluded.",
    dose_range       = "Loading 6 mg PO within 24 h of transplant surgery, followed by 0.5-6.0 mg/day maintenance titrated to whole-blood trough C0 of 6-12 ng/mL (months 1-3 and Phase I), 12-20 ng/mL (months 4-6 of Phase II) and 10-20 ng/mL (months 7-12 of Phase II). Cohort median DDS 2.0 mg/day.",
    regions          = "China (first Affiliated Hospital of Sun Yat-Sen University, Changzheng Hospital of Second Military Medical University, Shanghai No.1 People's Hospital of Jiaotong University, Chinese People's Liberation Army General Hospital, and Friendship Hospital of Capital Medical University).",
    transplant_type  = "Primary renal allograft; cadaveric n = 104 (93%), living donor n = 8 (7%).",
    co_medication    = "Ciclosporin (CsA, twice daily, sirolimus dosed 4 h after morning CsA), corticosteroids (no withdrawal permitted), and a wide range of other medications. Co-medications taken in > 10% of the cohort included antihypertensives (amlodipine, atenolol, captopril, clonidine, diltiazem, felodipine, metoprolol, nifedipine, losartan), diuretics (furosemide), lipid-lowering agents (fluvastatin, fenofibrate), aspirin, erythropoietin, H2-receptor antagonists (famotidine, ranitidine), proton-pump inhibitors (omeprazole), induction-therapy antibodies (antilymphocyte immunoglobulin, basiliximab, daclizumab) and herbal formulations (silymarin, cultured cordyceps mycelia, scutellarin, glycyrrhizin).",
    trial_design     = "Phase I (CsA dose reduction, n = 71) and Phase II (CsA discontinuation, n = 41) sub-cohorts; both arms received induction therapy with CsA + sirolimus + corticosteroids in the first 3 months before diverging. Phase I tapered CsA to a C0 target of 50-100 ng/mL by end of month 4; Phase II discontinued CsA by end of month 4. Trial registered as NCT00257387.",
    notes            = "Eight hundred and four sirolimus trough concentrations from 112 subjects; median 6 (range 1-18) concentrations per subject. All concentrations sampled at the end of the dosing interval with a median time of 24 h (range 22.8-25 h) post-dose. Data collected over a mean of 167 post-transplant days (median 172, range 4-388). Three external laboratories quantified sirolimus by two HPLC-UV methods (linear ranges 2.5-60 and 2.5-70 ng/mL, intra-/inter-day precision < 9.2% / 7.6%) and one LC-MS/MS method (linear range 0.2-50 ng/mL, intra-/inter-day precision < 15.4% / 11.6%). CsA was quantified by AxSYM (Shanghai sites) or TDx (other sites) automated fluorescence polarization with monoclonal antibodies; no assay-normalisation was applied to the original C0 values in the reported final model (sensitivity analysis showed assay-rescaling shifted CsA-related parameters < 20% with < 2 OFV units alternation)."
  )

  ini({
    # Final-model NONMEM estimates from Jiao 2009 Table 2 (NONMEM column).
    # Bootstrap median estimates (969 successful bootstrap fits out of 1000)
    # agreed with NONMEM point estimates within 5% across all parameters.

    # Apparent clearance and volume reference values are typical-individual
    # estimates AT THE COHORT MEDIAN COVARIATE VALUES: TCHOL = 5.66 mmol/L,
    # CP_CSA_NGML = 104 ng/mL, DOSE = 2 mg/day, CONMED_SILYMARIN = 0,
    # CONMED_GLYCYRRHIZIN = 0. The covariate-effect terms in model() shift
    # the typical value away from the reference for off-median covariates.
    lcl <- log(10.1);  label("Apparent clearance, CL/F, at cohort-median covariates (L/h)")  # Table 2 theta1 = 10.1 L/h (3.0% RSE)
    lvc <- log(3670);  label("Apparent central volume of distribution, V/F, at the CsA-trough reference (L)")  # Table 2 theta7 = 3670 L (9.3% RSE)

    # Absorption rate constant is FIXED at the literature value from Zahir
    # et al. (2006), reference [32] in Jiao 2009 [paper Methods: "Because no
    # data from the absorption phase were available, Ka was fixed at 0.752
    # l h-1 according to the literature and its interindividual variability
    # (IIV) was not estimated"]. The Results note that varying ka by 5- and
    # 10-fold changed final estimates by no more than 10% / 4%, confirming
    # the fit is insensitive to the exact ka value.
    lka <- fixed(log(0.752)); label("First-order absorption rate constant, ka (1/h; literature-fixed)")  # Jiao 2009 Methods (paper ref 32 = Zahir 2006)

    # Covariate effects on CL/F. Linear-deviation forms (TCHOL, CP_CSA_NGML)
    # carry the signed coefficient from Eq. 9; the power-form effects
    # (CONMED_SILYMARIN, CONMED_GLYCYRRHIZIN) are stored on the natural-log
    # scale so model() can write `exp(e_conmed_silymarin_cl * CONMED_SILYMARIN)`
    # to reproduce 0.65^SLM. The DOSE exponent enters directly as a power.
    e_tchol_cl              <- -0.662;       label("Linear-deviation effect of TCHOL on CL/F (L/h per mmol/L deviation from 5.66)")              # Table 2 theta2 = 0.662 (25.1% RSE); enters with a minus sign in the bracket (higher TCHOL -> lower CL/F)
    e_cp_csa_ngml_cl        <- -0.00417;     label("Linear-deviation effect of CP_CSA_NGML on CL/F (L/h per ng/mL deviation from 104)")         # Table 2 theta3 = 4.17e-3 (44.4% RSE); enters with a minus sign in the bracket (higher CsA C0 -> lower CL/F)
    e_conmed_silymarin_cl   <- log(0.650);   label("Log-effect of concomitant silymarin on CL/F (unitless)")                                    # Table 2 theta4 = 0.650 (16.3% RSE); 35% reduction in CL/F when CONMED_SILYMARIN = 1
    e_conmed_glycyrrhizin_cl <- log(0.661);  label("Log-effect of concomitant glycyrrhizin on CL/F (unitless)")                                 # Table 2 theta5 = 0.661 (25.4% RSE); 34% reduction in CL/F when CONMED_GLYCYRRHIZIN = 1
    e_dose_cl               <- 0.479;        label("Power-form exponent of (DOSE / 2) on CL/F (unitless)")                                       # Table 2 theta6 = 0.479 (16.6% RSE); (DDS/2)^0.479 multiplier on CL/F

    # Covariate effect on V/F.
    e_cp_csa_ngml_vc        <- -7.27;        label("Linear-deviation effect of CP_CSA_NGML on V/F (L per ng/mL deviation from 104)")            # Table 2 theta8 = 7.27 (14.9% RSE); enters with a minus sign in the V/F equation (higher CsA C0 -> lower V/F)

    # Inter-individual variability. Jiao 2009 Methods Eq. 1: P_i = P_typ * exp(eta_i)
    # (exponential / log-normal IIV). Table 2 reports interindividual variability
    # as a percentage (CV%); convert to the natural-log variance via
    # omega^2 = log(1 + CV^2). ka has no IIV by paper design (Methods, paper
    # ref 43). Off-diagonal correlations are not reported.
    etalcl ~ 0.0551   # 23.8% CV; omega^2 = log(1 + 0.238^2) = 0.0551          (Table 2 omega_CL/F = 23.8%, 26.7% RSE)
    etalvc ~ 0.2789   # 56.7% CV; omega^2 = log(1 + 0.567^2) = 0.2789          (Table 2 omega_V/F = 56.7%, 24.8% RSE)

    # Residual unexplained variability: exponential error model (Methods Eq. 3,
    # Y = IPRED * exp(eps)). Table 2 reports 29.9% (8.0% RSE); under the
    # exponential / log-normal residual the reported value is the SD of the
    # log-residual eps, so it maps directly to expSd. nlmixr2's `~ lnorm(expSd)`
    # implements Y ~ lnorm(IPRED, expSd) on the log-concentration scale.
    expSd <- 0.299;  label("Log-normal residual error SD (unitless)")  # Table 2 sigma = 29.9% reported as exponential-error magnitude
  })

  model({
    # Reference values from Jiao 2009 final-model Eq. 9-10 (cohort medians).
    tchol_ref   <- 5.66     # mmol/L (Table 1 cohort median TC)
    cp_csa_ref  <- 104      # ng/mL (Table 1 cohort median CsA C0)
    dose_ref    <- 2        # mg/day (Table 1 cohort median sirolimus DDS)

    # Apparent CL/F: linear-deviation TCHOL + CP_CSA_NGML bracket multiplied
    # by herb power-form factors and the DOSE power-form factor. IIV applied
    # multiplicatively on the full typical value per Methods Eq. 1. The
    # bracket reproduces Jiao 2009 Eq. 9 verbatim:
    #   typical CL/F = [theta1 + e_tchol_cl*(TCHOL - 5.66)
    #                   + e_cp_csa_ngml_cl*(CP_CSA_NGML - 104)]
    #                  * exp(e_conmed_silymarin_cl * CONMED_SILYMARIN)
    #                  * exp(e_conmed_glycyrrhizin_cl * CONMED_GLYCYRRHIZIN)
    #                  * (DOSE / dose_ref)^e_dose_cl
    cl_bracket <- exp(lcl) +
      e_tchol_cl * (TCHOL - tchol_ref) +
      e_cp_csa_ngml_cl * (CP_CSA_NGML - cp_csa_ref)
    cl_typ <- cl_bracket *
      exp(e_conmed_silymarin_cl * CONMED_SILYMARIN) *
      exp(e_conmed_glycyrrhizin_cl * CONMED_GLYCYRRHIZIN) *
      (DOSE / dose_ref)^e_dose_cl
    cl <- cl_typ * exp(etalcl)

    # Apparent V/F: linear-deviation CP_CSA_NGML factor (Eq. 10), with
    # IIV applied multiplicatively per Methods Eq. 1.
    vc_typ <- exp(lvc) + e_cp_csa_ngml_vc * (CP_CSA_NGML - cp_csa_ref)
    vc <- vc_typ * exp(etalvc)

    # Absorption rate constant: fixed for all subjects (no IIV).
    ka <- exp(lka)

    # Micro-rate constants.
    kel <- cl / vc

    # ODE system: oral one-compartment with first-order absorption (NONMEM
    # ADVAN2 / TRANS2 per Methods). Dose enters depot in mg, central state
    # in mg, vc in L -> central/vc gives mg/L; multiply by 1000 to convert
    # to ng/mL (reporting units of the sirolimus assay).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
