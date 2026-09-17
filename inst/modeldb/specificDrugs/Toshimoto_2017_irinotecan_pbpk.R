Toshimoto_2017_irinotecan_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, five-unit tandem dispersion liver, segregated-flow",
    "intestine, enterohepatic circulation). Irinotecan (CPT-11) and its four",
    "metabolites SN-38, SN-38G, NPC and APC after a 600 mg 90-minute",
    "intravenous infusion in adult cancer patients. Five structurally",
    "identical PBPK modules -- one per chemical species -- are coupled by",
    "hepatic and intestinal metabolic clearances, so the whole system is one",
    "jointly fitted model with 117 ODE states. Each module carries a central",
    "(blood) compartment, perfusion-limited muscle, skin, adipose and gut",
    "serosa, a permeability-limited liver resolved as five hepatic",
    "extracellular (sinusoidal) sub-compartments in series each exchanging",
    "with its own hepatocyte sub-compartment, a three-compartment biliary",
    "transit chain, an intestinal lumen, an enterocyte and a mucosal-blood",
    "compartment, and faecal and urinary sinks. The 46 biochemical parameters",
    "were fitted with the Cluster Newton method (CNM), which returns a family",
    "of 30 equally well-fitting parameter vectors rather than a single",
    "optimum; the values below are set ID 2, the set the authors themselves",
    "used for the paper's headline 1,000,000-virtual-patient simulation and",
    "one of the six sets satisfying more than five of their seven",
    "clinical-reproduction criteria; the vignette shows how to substitute",
    "any of the other 29. Between-subject variability is the lognormal enzyme /",
    "transporter / transit variability of Supplementary Table 3B and the",
    "normal physiological variability of Supplementary Table 3A; six",
    "germline polymorphisms (UGT1A1*28, SLCO1B1 c.521T>C and c.388A>G,",
    "ABCG2 c.421C>A, ABCB1 c.3435C>T, ABCC2 c.-24C>T) scale the enzyme and",
    "transporter activities they control, with the per-process contribution",
    "fractions of the Methods. Two integrator states accumulate the unbound",
    "SN-38 AUC in plasma and in the enterocyte, which are the paper's",
    "surrogates for neutropenia and for delayed diarrhoea. The paper reports",
    "no residual-error model (CNM is a fixed-effects least-squares method),",
    "so the propSd terms are placeholders. See the vignette Errata for the",
    "three mass-balance corrections applied to the published equations."
  )
  reference <- paste(
    "Toshimoto K, Tomaru A, Hosokawa M, Sugiyama Y. Virtual Clinical Studies",
    "to Examine the Probability Distribution of the AUC at Target Tissues",
    "Using Physiologically-Based Pharmacokinetic Modeling: Application to",
    "Analyses of the Effect of Genetic Polymorphism of Enzymes and",
    "Transporters on Irinotecan Induced Side Effects.",
    "Pharm Res. 2017;34(8):1584-1600. doi:10.1007/s11095-017-2153-z.",
    "The ODE system and the hybrid-to-elementary parameter conversions are",
    "transcribed from the Supplementary Text (ESM_1). Fixed physiological and",
    "physicochemical constants are Supplementary Table 1 (ESM_7); the",
    "genotype activity ratios and allele frequencies are Supplementary Table 2",
    "(ESM_8); the between-subject variability is Supplementary Table 3",
    "(ESM_9); the 30 CNM parameter sets are Supplementary Table 4 (ESM_10).",
    "Equations 1-5 are the article's own numbered equations, read from the",
    "publisher-supplied renderings Article_Equ1.gif to Article_Equ5.gif that",
    "ship in the same supplementary archive.",
    sep = " "
  )
  vignette <- "Toshimoto_2017_irinotecan_pbpk"

  # Compartment vocabulary. `central`, `muscle`, `skin`, `adipose`,
  # `gut_lumen`, `a_feces` and `a_urine` are canonical
  # (inst/references/compartment-names.md); `sn38` is the registered
  # metabolite suffix and `sn38g`, `npc` and `apc` are registered alongside
  # it by this extraction. The remaining stems follow the two registered
  # models of this same laboratory and model family,
  # Aoki_2024_bosentan_pbpk.R and Tsuchitani_2024_telmisartan_pbpk.R:
  # `is_liver<n>` is the hepatic extracellular (sinusoidal) space,
  # `int_liver<n>` the hepatocyte space, `serosa` the gut serosal tissue and
  # `ehc<n>` the biliary transit chain. `intestine_ent` (absorptive,
  # metabolising enterocyte) and `intestine_muc` (mucosal blood) are the
  # unsegmented analogues of that model's `<segment>_ent` / `<segment>_muc`;
  # the canonical `enterocyte` is deliberately NOT reused because its
  # register entry defines it as a terminal sink whose contents never reach
  # plasma, which is the opposite of the role here. Each of these stems now
  # has a second independent paper, so all are candidates for promotion to
  # canonical; per the standing ruling that promotion is a separate step,
  # they are declared paper-specific here.
  paper_specific_compartments <- c(
    "serosa",
    "is_liver1",
    "is_liver2",
    "is_liver3",
    "is_liver4",
    "is_liver5",
    "int_liver1",
    "int_liver2",
    "int_liver3",
    "int_liver4",
    "int_liver5",
    "ehc1",
    "ehc2",
    "ehc3",
    "intestine_ent",
    "intestine_muc",
    "serosa_sn38",
    "is_liver1_sn38",
    "is_liver2_sn38",
    "is_liver3_sn38",
    "is_liver4_sn38",
    "is_liver5_sn38",
    "int_liver1_sn38",
    "int_liver2_sn38",
    "int_liver3_sn38",
    "int_liver4_sn38",
    "int_liver5_sn38",
    "ehc1_sn38",
    "ehc2_sn38",
    "ehc3_sn38",
    "intestine_ent_sn38",
    "intestine_muc_sn38",
    "serosa_sn38g",
    "is_liver1_sn38g",
    "is_liver2_sn38g",
    "is_liver3_sn38g",
    "is_liver4_sn38g",
    "is_liver5_sn38g",
    "int_liver1_sn38g",
    "int_liver2_sn38g",
    "int_liver3_sn38g",
    "int_liver4_sn38g",
    "int_liver5_sn38g",
    "ehc1_sn38g",
    "ehc2_sn38g",
    "ehc3_sn38g",
    "intestine_ent_sn38g",
    "intestine_muc_sn38g",
    "serosa_npc",
    "is_liver1_npc",
    "is_liver2_npc",
    "is_liver3_npc",
    "is_liver4_npc",
    "is_liver5_npc",
    "int_liver1_npc",
    "int_liver2_npc",
    "int_liver3_npc",
    "int_liver4_npc",
    "int_liver5_npc",
    "ehc1_npc",
    "ehc2_npc",
    "ehc3_npc",
    "intestine_ent_npc",
    "intestine_muc_npc",
    "serosa_apc",
    "is_liver1_apc",
    "is_liver2_apc",
    "is_liver3_apc",
    "is_liver4_apc",
    "is_liver5_apc",
    "int_liver1_apc",
    "int_liver2_apc",
    "int_liver3_apc",
    "int_liver4_apc",
    "int_liver5_apc",
    "ehc1_apc",
    "ehc2_apc",
    "ehc3_apc",
    "intestine_ent_apc",
    "intestine_muc_apc",
    "auc_u_sn38",
    "auc_u_ent_sn38"
  )

  # Time in hours. Every volume, flow and clearance in Supplementary Tables 1
  # and 4 is expressed per kilogram of body weight and is multiplied by WT in
  # model(), so the tissue states hold concentrations in umol/L and the lumen,
  # bile, faecal and urinary states hold amounts in umol. Irinotecan doses are
  # supplied in umol (irinotecan hydrochloride trihydrate MW 677.19 g/mol, so
  # the 600 mg reference dose is 600000 / 677.19 = 886 umol; Supplementary
  # Table 1B states that dose as 12.66 umol/kg at 70 kg, which reproduces
  # 886 / 70 = 12.66 exactly and fixes both the salt form and the reference
  # weight).
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The tissue states are concentrations (umol/L); the
  # luminal, biliary, faecal and urinary states are amounts (umol).
  compartmentData <- list(
    central = list(analyte = "irinotecan", units = "umol/L", specimen = "whole blood", verified = TRUE),
    muscle = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    serosa = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver1 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver2 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver3 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver4 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver5 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver1 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver2 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver3 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver4 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver5 = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    ehc1 = list(analyte = "irinotecan", units = "umol", specimen = "bile", verified = TRUE),
    ehc2 = list(analyte = "irinotecan", units = "umol", specimen = "bile", verified = TRUE),
    ehc3 = list(analyte = "irinotecan", units = "umol", specimen = "bile", verified = TRUE),
    gut_lumen = list(analyte = "irinotecan", units = "umol", specimen = "administration site", verified = TRUE),
    intestine_ent = list(analyte = "irinotecan", units = "umol/L", specimen = "tissue", verified = TRUE),
    intestine_muc = list(analyte = "irinotecan", units = "umol/L", specimen = "whole blood", verified = TRUE),
    a_feces = list(analyte = "irinotecan", units = "umol", specimen = "faeces", verified = TRUE),
    a_urine = list(analyte = "irinotecan", units = "umol", specimen = "urine", verified = TRUE),
    central_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "whole blood", verified = TRUE),
    muscle_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    skin_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    adipose_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    serosa_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver1_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver2_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver3_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver4_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver5_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver1_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver2_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver3_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver4_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver5_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    ehc1_sn38 = list(analyte = "SN-38", units = "umol", specimen = "bile", verified = TRUE),
    ehc2_sn38 = list(analyte = "SN-38", units = "umol", specimen = "bile", verified = TRUE),
    ehc3_sn38 = list(analyte = "SN-38", units = "umol", specimen = "bile", verified = TRUE),
    gut_lumen_sn38 = list(analyte = "SN-38", units = "umol", specimen = "administration site", verified = TRUE),
    intestine_ent_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "tissue", verified = TRUE),
    intestine_muc_sn38 = list(analyte = "SN-38", units = "umol/L", specimen = "whole blood", verified = TRUE),
    a_feces_sn38 = list(analyte = "SN-38", units = "umol", specimen = "faeces", verified = TRUE),
    a_urine_sn38 = list(analyte = "SN-38", units = "umol", specimen = "urine", verified = TRUE),
    central_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "whole blood", verified = TRUE),
    muscle_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    skin_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    adipose_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    serosa_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver1_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver2_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver3_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver4_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver5_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver1_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver2_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver3_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver4_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver5_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    ehc1_sn38g = list(analyte = "SN-38G", units = "umol", specimen = "bile", verified = TRUE),
    ehc2_sn38g = list(analyte = "SN-38G", units = "umol", specimen = "bile", verified = TRUE),
    ehc3_sn38g = list(analyte = "SN-38G", units = "umol", specimen = "bile", verified = TRUE),
    gut_lumen_sn38g = list(analyte = "SN-38G", units = "umol", specimen = "administration site", verified = TRUE),
    intestine_ent_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "tissue", verified = TRUE),
    intestine_muc_sn38g = list(analyte = "SN-38G", units = "umol/L", specimen = "whole blood", verified = TRUE),
    a_feces_sn38g = list(analyte = "SN-38G", units = "umol", specimen = "faeces", verified = TRUE),
    a_urine_sn38g = list(analyte = "SN-38G", units = "umol", specimen = "urine", verified = TRUE),
    central_npc = list(analyte = "NPC", units = "umol/L", specimen = "whole blood", verified = TRUE),
    muscle_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    skin_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    adipose_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    serosa_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver1_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver2_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver3_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver4_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver5_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver1_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver2_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver3_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver4_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver5_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    ehc1_npc = list(analyte = "NPC", units = "umol", specimen = "bile", verified = TRUE),
    ehc2_npc = list(analyte = "NPC", units = "umol", specimen = "bile", verified = TRUE),
    ehc3_npc = list(analyte = "NPC", units = "umol", specimen = "bile", verified = TRUE),
    gut_lumen_npc = list(analyte = "NPC", units = "umol", specimen = "administration site", verified = TRUE),
    intestine_ent_npc = list(analyte = "NPC", units = "umol/L", specimen = "tissue", verified = TRUE),
    intestine_muc_npc = list(analyte = "NPC", units = "umol/L", specimen = "whole blood", verified = TRUE),
    a_feces_npc = list(analyte = "NPC", units = "umol", specimen = "faeces", verified = TRUE),
    a_urine_npc = list(analyte = "NPC", units = "umol", specimen = "urine", verified = TRUE),
    central_apc = list(analyte = "APC", units = "umol/L", specimen = "whole blood", verified = TRUE),
    muscle_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    skin_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    adipose_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    serosa_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver1_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver2_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver3_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver4_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver5_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver1_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver2_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver3_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver4_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver5_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    ehc1_apc = list(analyte = "APC", units = "umol", specimen = "bile", verified = TRUE),
    ehc2_apc = list(analyte = "APC", units = "umol", specimen = "bile", verified = TRUE),
    ehc3_apc = list(analyte = "APC", units = "umol", specimen = "bile", verified = TRUE),
    gut_lumen_apc = list(analyte = "APC", units = "umol", specimen = "administration site", verified = TRUE),
    intestine_ent_apc = list(analyte = "APC", units = "umol/L", specimen = "tissue", verified = TRUE),
    intestine_muc_apc = list(analyte = "APC", units = "umol/L", specimen = "whole blood", verified = TRUE),
    a_feces_apc = list(analyte = "APC", units = "umol", specimen = "faeces", verified = TRUE),
    a_urine_apc = list(analyte = "APC", units = "umol", specimen = "urine", verified = TRUE),
    auc_u_sn38 = list(analyte = "SN-38", units = "umol/L*h", specimen = "whole blood", verified = TRUE),
    auc_u_ent_sn38 = list(analyte = "SN-38", units = "umol/L*h", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = paste(
        "Total body weight. Every volume, blood flow and clearance in this",
        "model is tabulated per kilogram of body weight and is multiplied",
        "by WT inside model(), so WT sets the absolute size of the whole",
        "physiological system."
      ),
      units = "kg",
      type = "continuous",
      source_name = "Body weight",
      notes = paste(
        "Supplementary Table 3A gives the virtual-population distribution",
        "as normal with mean 74.87 kg and CV 15.2 percent, back-calculated",
        "by the authors from the reported body-surface-area distribution",
        "(Methods, Generation of Virtual Patients, item 1). Supplementary",
        "Table 1 footnote a uses 70 kg as the reference weight for the per-",
        "kilogram flows and volumes and for the 12.66 umol/kg reference",
        "dose."
      )
    ),
    SNP_UGT1A1_RS8175347_HET = list(
      description = paste(
        "Binary germline genotype indicator for UGT1A1 *28 (TA7 promoter",
        "repeat) (rs8175347): 1 if the subject carries the heterozygous",
        "genotype, 0 otherwise. Paired with SNP_UGT1A1_RS8175347_HOM; both",
        "indicators 0 is the wild-type homozygote reference stratum. Time-",
        "fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_UGT1A1_RS8175347_HOM also 0)",
      source_name = "UGT1A1 *28 (TA7 promoter repeat)",
      notes = paste(
        "Supplementary Table 2A gives the heterozygous activity ratio",
        "relative to wild type as 60.2 percent. In model() that activity",
        "scales the UGT1A1-mediated glucuronidation of SN-38 to SN-38G in",
        "the liver and in the enterocyte (contribution 100 percent to",
        "each)."
      )
    ),
    SNP_UGT1A1_RS8175347_HOM = list(
      description = paste(
        "Binary germline genotype indicator for UGT1A1 *28 (TA7 promoter",
        "repeat) (rs8175347): 1 if the subject carries the homozygous-",
        "variant genotype, 0 otherwise. Paired with",
        "SNP_UGT1A1_RS8175347_HET; both indicators 0 is the wild-type",
        "homozygote reference stratum. Time-fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_UGT1A1_RS8175347_HET also 0)",
      source_name = "UGT1A1 *28 (TA7 promoter repeat)",
      notes = paste(
        "Supplementary Table 2A gives the homozygous-variant activity ratio",
        "relative to wild type as 32.2 percent. In model() that activity",
        "scales the UGT1A1-mediated glucuronidation of SN-38 to SN-38G in",
        "the liver and in the enterocyte (contribution 100 percent to",
        "each)."
      )
    ),
    SNP_SLCO1B1_RS4149056_HET = list(
      description = paste(
        "Binary germline genotype indicator for SLCO1B1 c.521T>C (V174A)",
        "(rs4149056): 1 if the subject carries the heterozygous genotype, 0",
        "otherwise. Paired with SNP_SLCO1B1_RS4149056_HOM; both indicators",
        "0 is the wild-type homozygote reference stratum. Time-fixed per",
        "subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_SLCO1B1_RS4149056_HOM also 0)",
      source_name = "SLCO1B1 c.521T>C (V174A)",
      notes = paste(
        "Supplementary Table 2A gives the heterozygous activity ratio",
        "relative to wild type as 63.4 percent. In model() that activity",
        "scales the OATP1B1-mediated hepatic uptake of SN-38 and of SN-38G",
        "(contribution 100 percent to each), multiplicatively with the",
        "rs2306283 activity of the same gene."
      )
    ),
    SNP_SLCO1B1_RS4149056_HOM = list(
      description = paste(
        "Binary germline genotype indicator for SLCO1B1 c.521T>C (V174A)",
        "(rs4149056): 1 if the subject carries the homozygous-variant",
        "genotype, 0 otherwise. Paired with SNP_SLCO1B1_RS4149056_HET; both",
        "indicators 0 is the wild-type homozygote reference stratum. Time-",
        "fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_SLCO1B1_RS4149056_HET also 0)",
      source_name = "SLCO1B1 c.521T>C (V174A)",
      notes = paste(
        "Supplementary Table 2A gives the homozygous-variant activity ratio",
        "relative to wild type as 26.8 percent. In model() that activity",
        "scales the OATP1B1-mediated hepatic uptake of SN-38 and of SN-38G",
        "(contribution 100 percent to each), multiplicatively with the",
        "rs2306283 activity of the same gene."
      )
    ),
    SNP_SLCO1B1_RS2306283_HET = list(
      description = paste(
        "Binary germline genotype indicator for SLCO1B1 c.388A>G (N130D)",
        "(rs2306283): 1 if the subject carries the heterozygous genotype, 0",
        "otherwise. Paired with SNP_SLCO1B1_RS2306283_HOM; both indicators",
        "0 is the wild-type homozygote reference stratum. Time-fixed per",
        "subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_SLCO1B1_RS2306283_HOM also 0)",
      source_name = "SLCO1B1 c.388A>G (N130D)",
      notes = paste(
        "Supplementary Table 2A gives the heterozygous activity ratio",
        "relative to wild type as 161.0 percent. In model() that activity",
        "scales the OATP1B1-mediated hepatic uptake of SN-38 and of SN-38G",
        "(contribution 100 percent to each), multiplicatively with the",
        "rs4149056 activity of the same gene."
      )
    ),
    SNP_SLCO1B1_RS2306283_HOM = list(
      description = paste(
        "Binary germline genotype indicator for SLCO1B1 c.388A>G (N130D)",
        "(rs2306283): 1 if the subject carries the homozygous-variant",
        "genotype, 0 otherwise. Paired with SNP_SLCO1B1_RS2306283_HET; both",
        "indicators 0 is the wild-type homozygote reference stratum. Time-",
        "fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_SLCO1B1_RS2306283_HET also 0)",
      source_name = "SLCO1B1 c.388A>G (N130D)",
      notes = paste(
        "Supplementary Table 2A gives the homozygous-variant activity ratio",
        "relative to wild type as 221.0 percent. In model() that activity",
        "scales the OATP1B1-mediated hepatic uptake of SN-38 and of SN-38G",
        "(contribution 100 percent to each), multiplicatively with the",
        "rs4149056 activity of the same gene."
      )
    ),
    SNP_ABCG2_RS2231142_HET = list(
      description = paste(
        "Binary germline genotype indicator for ABCG2 c.421C>A (Q141K)",
        "(rs2231142): 1 if the subject carries the heterozygous genotype, 0",
        "otherwise. Paired with SNP_ABCG2_RS2231142_HOM; both indicators 0",
        "is the wild-type homozygote reference stratum. Time-fixed per",
        "subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_ABCG2_RS2231142_HOM also 0)",
      source_name = "ABCG2 c.421C>A (Q141K)",
      notes = paste(
        "Supplementary Table 2A gives the heterozygous activity ratio",
        "relative to wild type as 87.5 percent. In model() that activity",
        "scales the BCRP-mediated biliary excretion of SN-38 and the BCRP-",
        "mediated efflux of SN-38 from the enterocyte into the intestinal",
        "lumen (contribution 33 percent to each)."
      )
    ),
    SNP_ABCG2_RS2231142_HOM = list(
      description = paste(
        "Binary germline genotype indicator for ABCG2 c.421C>A (Q141K)",
        "(rs2231142): 1 if the subject carries the homozygous-variant",
        "genotype, 0 otherwise. Paired with SNP_ABCG2_RS2231142_HET; both",
        "indicators 0 is the wild-type homozygote reference stratum. Time-",
        "fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_ABCG2_RS2231142_HET also 0)",
      source_name = "ABCG2 c.421C>A (Q141K)",
      notes = paste(
        "Supplementary Table 2A gives the homozygous-variant activity ratio",
        "relative to wild type as 51.3 percent. In model() that activity",
        "scales the BCRP-mediated biliary excretion of SN-38 and the BCRP-",
        "mediated efflux of SN-38 from the enterocyte into the intestinal",
        "lumen (contribution 33 percent to each)."
      )
    ),
    ABCB1_C3435T_HET = list(
      description = paste(
        "Binary germline genotype indicator for ABCB1 c.3435C>T",
        "(rs1045642): 1 if the subject carries the heterozygous genotype, 0",
        "otherwise. Paired with ABCB1_C3435T_MUT; both indicators 0 is the",
        "wild-type homozygote reference stratum. Time-fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with ABCB1_C3435T_MUT also 0)",
      source_name = "ABCB1 c.3435C>T",
      notes = paste(
        "Supplementary Table 2A gives the heterozygous activity ratio",
        "relative to wild type as 66.6 percent. In model() that activity",
        "scales the MDR1-mediated biliary excretion and enterocyte-to-lumen",
        "efflux of irinotecan (contribution 50 percent to each) and of",
        "SN-38 (33 percent to each)."
      )
    ),
    ABCB1_C3435T_MUT = list(
      description = paste(
        "Binary germline genotype indicator for ABCB1 c.3435C>T",
        "(rs1045642): 1 if the subject carries the homozygous-variant",
        "genotype, 0 otherwise. Paired with ABCB1_C3435T_HET; both",
        "indicators 0 is the wild-type homozygote reference stratum. Time-",
        "fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with ABCB1_C3435T_HET also 0)",
      source_name = "ABCB1 c.3435C>T",
      notes = paste(
        "Supplementary Table 2A gives the homozygous-variant activity ratio",
        "relative to wild type as 33.3 percent. In model() that activity",
        "scales the MDR1-mediated biliary excretion and enterocyte-to-lumen",
        "efflux of irinotecan (contribution 50 percent to each) and of",
        "SN-38 (33 percent to each)."
      )
    ),
    SNP_ABCC2_RS717620_HET = list(
      description = paste(
        "Binary germline genotype indicator for ABCC2 c.-24C>T (rs717620):",
        "1 if the subject carries the heterozygous genotype, 0 otherwise.",
        "Paired with SNP_ABCC2_RS717620_HOM; both indicators 0 is the wild-",
        "type homozygote reference stratum. Time-fixed per subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_ABCC2_RS717620_HOM also 0)",
      source_name = "ABCC2 c.-24C>T",
      notes = paste(
        "Supplementary Table 2A gives the heterozygous activity ratio",
        "relative to wild type as 66.6 percent. In model() that activity",
        "scales the MRP2-mediated biliary excretion and enterocyte-to-lumen",
        "efflux of irinotecan (contribution 50 percent to each), of SN-38",
        "(33 percent to each) and of SN-38G (100 percent to each)."
      )
    ),
    SNP_ABCC2_RS717620_HOM = list(
      description = paste(
        "Binary germline genotype indicator for ABCC2 c.-24C>T (rs717620):",
        "1 if the subject carries the homozygous-variant genotype, 0",
        "otherwise. Paired with SNP_ABCC2_RS717620_HET; both indicators 0",
        "is the wild-type homozygote reference stratum. Time-fixed per",
        "subject."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type homozygote, with SNP_ABCC2_RS717620_HET also 0)",
      source_name = "ABCC2 c.-24C>T",
      notes = paste(
        "Supplementary Table 2A gives the homozygous-variant activity ratio",
        "relative to wild type as 33.3 percent. In model() that activity",
        "scales the MRP2-mediated biliary excretion and enterocyte-to-lumen",
        "efflux of irinotecan (contribution 50 percent to each), of SN-38",
        "(33 percent to each) and of SN-38G (100 percent to each)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = paste(
      "The pharmacokinetic parameters were fitted to the mean blood",
      "concentration-time profiles of irinotecan, SN-38, SN-38G, NPC and",
      "APC reported by van der Bol et al. 2011 (Eur J Cancer 47:831-8), a",
      "cross-over omeprazole interaction study in adult cancer patients.",
      "Each virtual clinical study was sized to the 127 advanced or",
      "metastatic cancer patients of the target study, Teft et al. 2015",
      "(Br J Cancer 112:857-65), and the final probability simulation",
      "used 1,000,000 virtual patients."
    ),
    disease_state = "advanced or metastatic solid tumours treated with irinotecan-based chemotherapy",
    dose_range = paste(
      "600 mg irinotecan hydrochloride trihydrate (886 umol; 12.66",
      "umol/kg at 70 kg) as a single 90-minute intravenous infusion."
    ),
    weight_range = "normal, mean 74.87 kg, CV 15.2 percent (Supplementary Table 3A)",
    regions = paste(
      "Allele frequencies are of mixed provenance: UGT1A1*28 and ABCC2",
      "c.-24C>T from the Teft 2015 target-study cohort, SLCO1B1 c.521T>C",
      "and c.388A>G (with their linkage-disequilibrium joint frequencies)",
      "from the global SLCO1B1 survey of Pasanen 2008, ABCG2 c.421C>A",
      "from Zamber 2003, and ABCB1 c.3435C>T from the 1000 Genomes",
      "European population (Supplementary Table 2A and its footnote c)."
    ),
    notes = paste(
      "Ages, sexes and baseline laboratory values are not reported by",
      "Toshimoto 2017; they belong to the two underlying clinical",
      "studies. Carboxylesterase genotype was deliberately excluded from",
      "the virtual population: the authors' own anti-CES1A",
      "immunodepletion experiment (Results, Figure 3) showed that",
      "removing enough CES1A activity to cut temocapril hydrolysis to",
      "45.1 percent of control left irinotecan hydrolysis at 95.3 percent",
      "of control, establishing that CES2 rather than CES1 forms SN-38 in",
      "human liver."
    )
  )
  ini({
    # ---------------------------------------------------------------
    # 1. Cluster Newton method (CNM) parameters -- Supplementary Table 4,
    #    parameter set ID 2. The CNM returns 30 parameter vectors that fit
    #    the observed profiles equally well (average weighted sum of squares
    #    0.082-0.107 across the 30); ID 2 is the set the authors used for
    #    the 1,000,000-virtual-patient probability simulation reported in
    #    the Discussion, and is one of the six sets satisfying more than
    #    five of their seven clinical-reproduction criteria (Figure 8).
    #    These are estimates, not fixed constants, so they carry no
    #    fixed() wrapper. The vignette tabulates all 30 sets.
    # ---------------------------------------------------------------

    # -- Irinotecan (Supplementary Table 4A, row ID 2). Irinotecan is the
    #    one compound whose hepatic parameters were optimised as
    #    elementary rather than hybrid quantities (Methods, Development of
    #    the PBPK Model).
    lv_central <- log(0.119)
    label("Irinotecan central (blood) volume (L/kg)") # Suppl Table 4A ID 2: Vcentral 0.119
    lka <- log(0.512)
    label("Irinotecan lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table 4A ID 2: ka 0.512
    lps_dif_inf_h <- log(1.045)
    label("Irinotecan sinusoidal passive-diffusion influx clearance (L/h/kg)") # Suppl Table 4A ID 2: PSdif,inf,h 1.045
    lkbile <- log(0.098)
    label("Irinotecan biliary transit rate constant (1/h)") # Suppl Table 4A ID 2: kbile 0.098
    lkfeces <- log(3.047)
    label("Irinotecan faecal transit rate constant (1/h)") # Suppl Table 4A ID 2: kfeces 3.047
    lcl_sn38_h <- log(0.197)
    label("Hepatic intrinsic clearance, irinotecan to SN-38 (L/h/kg)") # Suppl Table 4A ID 2: CL SN-38,h 0.197
    lcl_npc_h <- log(5.726)
    label("Hepatic intrinsic clearance, irinotecan to NPC (L/h/kg)") # Suppl Table 4A ID 2: CL NPC,h 5.726
    lcl_apc_h <- log(8.957)
    label("Hepatic intrinsic clearance, irinotecan to APC (L/h/kg)") # Suppl Table 4A ID 2: CL APC,h 8.957
    lcl_others_h <- log(0.149)
    label("Hepatic intrinsic clearance, irinotecan to other products (L/h/kg)") # Suppl Table 4A ID 2: CL others,h 0.149
    lcl_bile <- log(0.038)
    label("Irinotecan biliary intrinsic clearance (L/h/kg)") # Suppl Table 4A ID 2: CL bile 0.038
    lr_met_ces <- log(1.379)
    label("Hepatic-to-intestinal ratio of CES-mediated SN-38 formation (unitless)") # Suppl Table 4A ID 2: CL SN-38,h / CL SN-38,ent 1.379
    lps_dif_eff_ent <- log(0.085)
    label("Irinotecan enterocyte basolateral passive-diffusion clearance (L/h/kg)") # Suppl Table 4A ID 2: PSdif,eff,ent 0.085
    rdif_ent <- 0.399
    label("Irinotecan enterocyte active-to-passive apical efflux ratio (unitless)") # Suppl Table 4A ID 2: Rdif,ent 0.399
    lr_met_cyp3a <- log(1.106)
    label("Hepatic-to-intestinal ratio of CYP3A-mediated NPC and APC formation (unitless)") # Suppl Table 4A ID 2: CL NPC,h / CL NPC,ent 1.106

    # -- SN-38 (Supplementary Table 4B, row ID 2).
    lv_central_sn38 <- log(0.078)
    label("SN-38 central (blood) volume (L/kg)") # Suppl Table 4B ID 2: Vcentral 0.078
    lka_sn38 <- log(6.454)
    label("SN-38 lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table 4B ID 2: ka 6.454
    rdif_h_sn38 <- 0.188
    label("SN-38 hepatic passive-to-active uptake ratio Rdif,h (unitless)") # Suppl Table 4B ID 2: Rdif,h 0.188
    inv_beta_sn38 <- 1.388
    label("SN-38 reciprocal of the hepatic elimination fraction beta (unitless)") # Suppl Table 4B ID 2: 1/beta 1.388
    lclint_all_sn38 <- log(12.078)
    label("SN-38 hepatic overall intrinsic clearance (L/h/kg)") # Suppl Table 4B ID 2: CLint,all 12.078
    inv_fglu_sn38 <- 1.319
    label("SN-38 reciprocal of the glucuronidated fraction fglu (unitless)") # Suppl Table 4B ID 2: 1/fglu 1.319; Suppl Table 4 footnote a: fglu = 1 - fbile
    lkbile_sn38 <- log(0.053)
    label("SN-38 biliary transit rate constant (1/h)") # Suppl Table 4B ID 2: kbile 0.053
    lkfeces_sn38 <- log(1.053)
    label("SN-38 faecal transit rate constant (1/h)") # Suppl Table 4B ID 2: kfeces 1.053
    rdif_ent_sn38 <- 7.277
    label("SN-38 enterocyte active-to-passive apical efflux ratio (unitless)") # Suppl Table 4B ID 2: Rdif,ent 7.277
    lr_ugt_h_ent_sn38 <- log(1.852)
    label("Hepatic-to-intestinal ratio of UGT1A1-mediated SN-38G formation (unitless)") # Suppl Table 4B ID 2: CL SN-38G,h / CL SN-38G,ent 1.852

    # -- SN-38G (Supplementary Table 4C, row ID 2).
    lv_central_sn38g <- log(0.289)
    label("SN-38G central (blood) volume (L/kg)") # Suppl Table 4C ID 2: Vcentral 0.289
    lka_sn38g <- log(1.527)
    label("SN-38G lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table 4C ID 2: ka 1.527
    rdif_h_sn38g <- 0.924
    label("SN-38G hepatic passive-to-active uptake ratio Rdif,h (unitless)") # Suppl Table 4C ID 2: Rdif,h 0.924
    inv_beta_sn38g <- 9.95
    label("SN-38G reciprocal of the hepatic elimination fraction beta (unitless)") # Suppl Table 4C ID 2: 1/beta 9.950
    lclint_all_sn38g <- log(0.156)
    label("SN-38G hepatic overall intrinsic clearance (L/h/kg)") # Suppl Table 4C ID 2: CLint,all 0.156
    lkbile_sn38g <- log(0.156)
    label("SN-38G biliary transit rate constant (1/h)") # Suppl Table 4C ID 2: kbile 0.156
    lkfeces_sn38g <- log(0.224)
    label("SN-38G faecal transit rate constant (1/h)") # Suppl Table 4C ID 2: kfeces 0.224
    rdif_ent_sn38g <- 3.103
    label("SN-38G enterocyte active-to-passive apical efflux ratio (unitless)") # Suppl Table 4C ID 2: Rdif,ent 3.103
    lkdec_sn38g <- log(0.083)
    label("Luminal deconjugation rate constant, SN-38G to SN-38 (1/h)") # Suppl Table 4C ID 2: kdec 0.083

    # -- NPC (Supplementary Table 4D, row ID 2). NPC is the only
    #    metabolite with BOTH hepatic metabolism (back to SN-38, via CES)
    #    and biliary excretion, which is why it is the only compound for
    #    which Table I lists 1/fbile as an unknown to be optimised.
    lv_central_npc <- log(0.147)
    label("NPC central (blood) volume (L/kg)") # Suppl Table 4D ID 2: Vcentral 0.147
    lka_npc <- log(2.147)
    label("NPC lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table 4D ID 2: ka 2.147
    inv_beta_npc <- 3.631
    label("NPC reciprocal of the hepatic elimination fraction beta (unitless)") # Suppl Table 4D ID 2: 1/beta 3.631
    lclint_all_npc <- log(13.652)
    label("NPC hepatic overall intrinsic clearance (L/h/kg)") # Suppl Table 4D ID 2: CLint,all 13.652
    inv_fbile_npc <- 1.441
    label("NPC reciprocal of the biliary fraction fbile (unitless)") # Suppl Table 4D ID 2: 1/fbile 1.441
    lkbile_npc <- log(0.169)
    label("NPC biliary transit rate constant (1/h)") # Suppl Table 4D ID 2: kbile 0.169
    lkfeces_npc <- log(0.045)
    label("NPC faecal transit rate constant (1/h)") # Suppl Table 4D ID 2: kfeces 0.045

    # -- APC (Supplementary Table 4E, row ID 2).
    lv_central_apc <- log(0.087)
    label("APC central (blood) volume (L/kg)") # Suppl Table 4E ID 2: Vcentral 0.087
    lka_apc <- log(0.253)
    label("APC lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table 4E ID 2: ka 0.253
    inv_beta_apc <- 16.932
    label("APC reciprocal of the hepatic elimination fraction beta (unitless)") # Suppl Table 4E ID 2: 1/beta 16.932
    lclint_all_apc <- log(1.428)
    label("APC hepatic overall intrinsic clearance (L/h/kg)") # Suppl Table 4E ID 2: CLint,all 1.428
    lkbile_apc <- log(4.031)
    label("APC biliary transit rate constant (1/h)") # Suppl Table 4E ID 2: kbile 4.031
    lkfeces_apc <- log(2.236)
    label("APC faecal transit rate constant (1/h)") # Suppl Table 4E ID 2: kfeces 2.236

    # ---------------------------------------------------------------
    # 2. Fixed physiological constants -- Supplementary Table 1A, scaled
    #    per kilogram of body weight. Identical for all five species.
    # ---------------------------------------------------------------
    lq_liver <- fixed(log(1.242))
    label("Hepatic blood flow (L/h/kg)") # Suppl Table 1A: Liver blood flow 1.242
    lq_muscle <- fixed(log(0.642))
    label("Muscle blood flow (L/h/kg)") # Suppl Table 1A: Muscle blood flow 0.642
    lq_skin <- fixed(log(0.257))
    label("Skin blood flow (L/h/kg)") # Suppl Table 1A: Skin blood flow 0.257
    lq_adipose <- fixed(log(0.223))
    label("Adipose blood flow (L/h/kg)") # Suppl Table 1A: Adipose blood flow 0.223
    lq_mucosa <- fixed(log(0.257))
    label("Intestinal mucosal blood flow (L/h/kg)") # Suppl Table 1A: Mucosa (Enterocyte) blood flow 0.257
    lq_serosa <- fixed(log(0.274))
    label("Intestinal serosal blood flow (L/h/kg)") # Suppl Table 1A: Serosa blood flow 0.274
    lv_liver <- fixed(log(0.0241))
    label("Total liver volume (L/kg)") # Suppl Table 1A: Liver volume 0.0241
    lv_muscle <- fixed(log(0.429))
    label("Muscle volume (L/kg)") # Suppl Table 1A: Muscle volume 0.429
    lv_skin <- fixed(log(0.111))
    label("Skin volume (L/kg)") # Suppl Table 1A: Skin volume 0.111
    lv_adipose <- fixed(log(0.143))
    label("Adipose volume (L/kg)") # Suppl Table 1A: Adipose volume 0.143
    lv_mucosa <- fixed(log(0.0075))
    label("Total intestinal mucosa volume (L/kg)") # Suppl Table 1A: Mucosa (Enterocyte) volume 0.0075
    lv_serosa <- fixed(log(0.0089))
    label("Intestinal serosa volume (L/kg)") # Suppl Table 1A: Serosa volume 0.0089
    frac_liver_ex <- fixed(0.278)
    label("Extracellular fraction of the liver volume (unitless)") # Suppl Table 1A: Liver extracellular space 0.278
    frac_mucosa_ex <- fixed(0.118)
    label("Extracellular (mucosal blood) fraction of the mucosa volume (unitless)") # Suppl Table 1A: Mucosa extracellular space 0.118
    ar <- fixed(20)
    label("Enterocyte apical-to-basolateral surface-area ratio (unitless)") # Suppl Table 1A: Apical/Basolateral area ratio (AR) 20
    tinf <- fixed(1.5)
    label("Intravenous infusion duration (h)") # Methods and Results: SN-38 concentration at the end of infusion (90 min)

    # ---------------------------------------------------------------
    # 3. Fixed physicochemical constants -- Supplementary Table 1B.
    #    fb is the unbound fraction in BLOOD (fp/Rb), fh and fgut the
    #    unbound fractions in hepatocyte and enterocyte, and Kp the
    #    tissue-to-blood concentration ratios. SN-38G is unbound in both
    #    intracellular compartments (fh = fgut = 1), and its blood unbound
    #    fraction exceeds 1 because Rb 0.55 is below unity -- fb = fp/Rb =
    #    0.691/0.55 = 1.26 as tabulated.
    # ---------------------------------------------------------------
    fb <- fixed(0.307)
    label("irinotecan unbound fraction in blood (unitless)") # Suppl Table 1B: fb 0.307
    fh <- fixed(0.0356)
    label("irinotecan unbound fraction in the hepatocyte (unitless)") # Suppl Table 1B: fh 0.0356
    fgut <- fixed(0.0563)
    label("irinotecan unbound fraction in the enterocyte (unitless)") # Suppl Table 1B: fgut 0.0563
    lkp_muscle <- fixed(log(3.953))
    label("irinotecan muscle-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,muscle 3.953
    lkp_skin <- fixed(log(4.277))
    label("irinotecan skin-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,skin 4.277
    lkp_adipose <- fixed(log(2.831))
    label("irinotecan adipose-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,adipose 2.831
    lkp_gut <- fixed(log(6.291))
    label("irinotecan gut-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,gut 6.291
    lcl_renal <- fixed(log(0.0977))
    label("irinotecan renal clearance from blood (L/h/kg)") # Suppl Table 1B: CLr 0.0977
    fb_sn38 <- fixed(0.022)
    label("SN-38 unbound fraction in blood (unitless)") # Suppl Table 1B: fb 0.022
    fh_sn38 <- fixed(0.227)
    label("SN-38 unbound fraction in the hepatocyte (unitless)") # Suppl Table 1B: fh 0.227
    fgut_sn38 <- fixed(0.123)
    label("SN-38 unbound fraction in the enterocyte (unitless)") # Suppl Table 1B: fgut 0.123
    lkp_muscle_sn38 <- fixed(log(0.0799))
    label("SN-38 muscle-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,muscle 0.0799
    lkp_skin_sn38 <- fixed(log(0.19))
    label("SN-38 skin-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,skin 0.19
    lkp_adipose_sn38 <- fixed(log(0.16))
    label("SN-38 adipose-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,adipose 0.16
    lkp_gut_sn38 <- fixed(log(0.19))
    label("SN-38 gut-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,gut 0.19
    lcl_renal_sn38 <- fixed(log(0.0475))
    label("SN-38 renal clearance from blood (L/h/kg)") # Suppl Table 1B: CLr 0.0475
    fb_sn38g <- fixed(1.26)
    label("SN-38G unbound fraction in blood (unitless)") # Suppl Table 1B: fb 1.26
    fh_sn38g <- fixed(1)
    label("SN-38G unbound fraction in the hepatocyte (unitless)") # Suppl Table 1B: fh 1
    fgut_sn38g <- fixed(1)
    label("SN-38G unbound fraction in the enterocyte (unitless)") # Suppl Table 1B: fgut 1
    lkp_muscle_sn38g <- fixed(log(0.488))
    label("SN-38G muscle-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,muscle 0.488
    lkp_skin_sn38g <- fixed(log(0.774))
    label("SN-38G skin-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,skin 0.774
    lkp_adipose_sn38g <- fixed(log(0.129))
    label("SN-38G adipose-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,adipose 0.129
    lkp_gut_sn38g <- fixed(log(0.671))
    label("SN-38G gut-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,gut 0.671
    lcl_renal_sn38g <- fixed(log(0.129))
    label("SN-38G renal clearance from blood (L/h/kg)") # Suppl Table 1B: CLr 0.129
    fb_npc <- fixed(0.4)
    label("NPC unbound fraction in blood (unitless)") # Suppl Table 1B: fb 0.4
    fh_npc <- fixed(0.0352)
    label("NPC unbound fraction in the hepatocyte (unitless)") # Suppl Table 1B: fh 0.0352
    fgut_npc <- fixed(0.0614)
    label("NPC unbound fraction in the enterocyte (unitless)") # Suppl Table 1B: fgut 0.0614
    lkp_muscle_npc <- fixed(log(5.08))
    label("NPC muscle-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,muscle 5.08
    lkp_skin_npc <- fixed(log(4.25))
    label("NPC skin-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,skin 4.25
    lkp_adipose_npc <- fixed(log(1.23))
    label("NPC adipose-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,adipose 1.23
    lkp_gut_npc <- fixed(log(7.54))
    label("NPC gut-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,gut 7.54
    cl_renal_npc <- fixed(0)
    label("NPC renal clearance from blood (L/h/kg)") # Suppl Table 1B: CLr 0 (tabulated as 0, footnote 'Assumption (No information)')
    fb_apc <- fixed(0.331)
    label("APC unbound fraction in blood (unitless)") # Suppl Table 1B: fb 0.331
    fh_apc <- fixed(0.042)
    label("APC unbound fraction in the hepatocyte (unitless)") # Suppl Table 1B: fh 0.042
    fgut_apc <- fixed(0.072)
    label("APC unbound fraction in the enterocyte (unitless)") # Suppl Table 1B: fgut 0.072
    lkp_muscle_apc <- fixed(log(1.48))
    label("APC muscle-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,muscle 1.48
    lkp_skin_apc <- fixed(log(1.29))
    label("APC skin-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,skin 1.29
    lkp_adipose_apc <- fixed(log(0.371))
    label("APC adipose-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,adipose 0.371
    lkp_gut_apc <- fixed(log(2.19))
    label("APC gut-to-blood partition coefficient (unitless)") # Suppl Table 1B: Kp,gut 2.19
    lcl_renal_apc <- fixed(log(0.0621))
    label("APC renal clearance from blood (L/h/kg)") # Suppl Table 1B: CLr 0.0621

    # ---------------------------------------------------------------
    # 4. Genotype activity ratios -- Supplementary Table 2A, column
    #    "Activity ratio to wild type (%)", divided by 100. Each is the
    #    multiplicative activity of the affected enzyme or transporter
    #    relative to the wild-type homozygote, so 1.0 would be no effect.
    #    They are experimental / assumed inputs to the virtual-population
    #    generator rather than quantities estimated from the PK data, so
    #    they are fixed.
    # ---------------------------------------------------------------
    e_ugt1a1_het_clglu <- fixed(0.602)
    label("Relative UGT1A1-mediated SN-38 glucuronidation activity, UGT1A1 *28 heterozygote (unitless)") # Suppl Table 2A row 'UGT1A1 *28': heterozygote 60.2 percent of wild type; reference (1)a, activity ratio for UGT1A1 *6 used
    e_ugt1a1_hom_clglu <- fixed(0.322)
    label("Relative UGT1A1-mediated SN-38 glucuronidation activity, UGT1A1 *28 homozygous variant (unitless)") # Suppl Table 2A row 'UGT1A1 *28': homozygote 32.2 percent of wild type; reference (1)a, activity ratio for UGT1A1 *6 used
    e_slco1b1_521_het_psact <- fixed(0.634)
    label("Relative OATP1B1-mediated hepatic uptake activity, SLCO1B1 521T>C heterozygote (unitless)") # Suppl Table 2A row 'SLCO1B1 521T>C': heterozygote 63.4 percent of wild type; reference in-house, unpublished (footnote b)
    e_slco1b1_521_hom_psact <- fixed(0.268)
    label("Relative OATP1B1-mediated hepatic uptake activity, SLCO1B1 521T>C homozygous variant (unitless)") # Suppl Table 2A row 'SLCO1B1 521T>C': homozygote 26.8 percent of wild type; reference in-house, unpublished (footnote b)
    e_slco1b1_388_het_psact <- fixed(1.61)
    label("Relative OATP1B1-mediated hepatic uptake activity, SLCO1B1 388A>G heterozygote (unitless)") # Suppl Table 2A row 'SLCO1B1 388A>G': heterozygote 161 percent of wild type; reference in-house, unpublished (footnote b)
    e_slco1b1_388_hom_psact <- fixed(2.21)
    label("Relative OATP1B1-mediated hepatic uptake activity, SLCO1B1 388A>G homozygous variant (unitless)") # Suppl Table 2A row 'SLCO1B1 388A>G': homozygote 221 percent of wild type; reference in-house, unpublished (footnote b)
    e_abcg2_het_efflux <- fixed(0.875)
    label("Relative BCRP-mediated efflux activity, ABCG2 421C>A heterozygote (unitless)") # Suppl Table 2A row 'ABCG2 421C>A': heterozygote 87.5 percent of wild type; reference in-house, unpublished (footnote b)
    e_abcg2_hom_efflux <- fixed(0.513)
    label("Relative BCRP-mediated efflux activity, ABCG2 421C>A homozygous variant (unitless)") # Suppl Table 2A row 'ABCG2 421C>A': homozygote 51.3 percent of wild type; reference in-house, unpublished (footnote b)
    e_abcb1_het_efflux <- fixed(0.666)
    label("Relative MDR1-mediated efflux activity, ABCB1 3435C>T heterozygote (unitless)") # Suppl Table 2A row 'ABCB1 3435C>T': heterozygote 66.6 percent of wild type; reference Assumption
    e_abcb1_hom_efflux <- fixed(0.333)
    label("Relative MDR1-mediated efflux activity, ABCB1 3435C>T homozygous variant (unitless)") # Suppl Table 2A row 'ABCB1 3435C>T': homozygote 33.3 percent of wild type; reference Assumption
    e_abcc2_het_efflux <- fixed(0.666)
    label("Relative MRP2-mediated efflux activity, ABCC2 -24C>T heterozygote (unitless)") # Suppl Table 2A row 'ABCC2 -24C>T': heterozygote 66.6 percent of wild type; reference Assumption
    e_abcc2_hom_efflux <- fixed(0.333)
    label("Relative MRP2-mediated efflux activity, ABCC2 -24C>T homozygous variant (unitless)") # Suppl Table 2A row 'ABCC2 -24C>T': homozygote 33.3 percent of wild type; reference Assumption

    # ---------------------------------------------------------------
    # 5. Transporter contribution fractions -- Methods, Generation of
    #    Virtual Patients. These say what share of each biliary or
    #    enterocyte-efflux clearance a given transporter carries, and so
    #    how far a genotype or a transporter-specific random effect can
    #    move that clearance. The three 33 percent shares of the SN-38
    #    efflux processes are the paper's own rounding of one third and
    #    sum to 0.99; the residual 1 percent is carried at wild-type
    #    activity by the (1 - sum) term in model().
    # ---------------------------------------------------------------
    frac_abcb1_iri <- fixed(0.5)
    label("MDR1 share of irinotecan biliary excretion and enterocyte efflux (unitless)") # Methods, Generation of Virtual Patients: ABCB1 ... contributions are 50% for irinotecan
    frac_abcc2_iri <- fixed(0.5)
    label("MRP2 share of irinotecan biliary excretion and enterocyte efflux (unitless)") # Methods, Generation of Virtual Patients: ABCC2 ... contributions are 50% for irinotecan
    frac_abcg2_sn38 <- fixed(0.33)
    label("BCRP share of SN-38 biliary excretion and enterocyte efflux (unitless)") # Methods, Generation of Virtual Patients: ABCG2 ... contributions for both are 33%
    frac_abcb1_sn38 <- fixed(0.33)
    label("MDR1 share of SN-38 biliary excretion and enterocyte efflux (unitless)") # Methods, Generation of Virtual Patients: ABCB1 ... 33% for SN-38
    frac_abcc2_sn38 <- fixed(0.33)
    label("MRP2 share of SN-38 biliary excretion and enterocyte efflux (unitless)") # Methods, Generation of Virtual Patients: ABCC2 ... 33% for SN-38
    frac_abcc2_sn38g <- fixed(1.0)
    label("MRP2 share of SN-38G biliary excretion and enterocyte efflux (unitless)") # Methods, Generation of Virtual Patients: ABCC2 ... 100% for SN-38G

    # ---------------------------------------------------------------
    # 6. Intestinal clearance of the irinotecan "other products"
    #    pathway. The pathway exists in the Supplementary Text
    #    enterocyte equation as CL others,ent(irinotecan), but Table I
    #    lists hepatic-to-intestinal ratios for only two of the four
    #    irinotecan metabolic routes (CES to SN-38 and CYP3A to NPC and
    #    APC), so no value for it was ever optimised. It is therefore
    #    fixed at zero, which is the only reading consistent with the
    #    published unknown-parameter list. See the vignette Errata.
    # ---------------------------------------------------------------
    cl_others_ent <- fixed(0)
    label("Intestinal intrinsic clearance, irinotecan to other products (L/h/kg)") # not optimised: Table I lists no hepatic/intestinal ratio for CL others

    # ---------------------------------------------------------------
    # 7. Between-subject variability -- Supplementary Table 3.
    #    Table 3B parameters are LOGNORMAL, so the omega below is
    #    log(CV^2 + 1) and the eta enters model() as exp(eta). Table 3A
    #    physiological parameters are NORMAL, so their omega is CV^2 and
    #    the eta enters as (1 + eta). One random effect is drawn per
    #    PROCESS, not per compound, so a subject with high OATP1B1
    #    activity takes up both SN-38 and SN-38G faster -- that shared
    #    draw is what makes the transporter genotypes interpretable.
    # ---------------------------------------------------------------
    etaloatp1b1 ~ 0.0644423 # Suppl Table 3B row 'PSact,inf,h by OATP1B1': lognormal, CV 25.8 percent; omega = log(1 + 0.258^2). IIV on OATP1B1-mediated hepatic uptake (PSact,inf,h)
    etalps_dif_h ~ 0.0099503 # Suppl Table 3B row 'PSdif,inf,h and PSdif,eff,h': lognormal, CV 10 percent; omega = log(1 + 0.1^2). IIV on hepatic sinusoidal passive diffusion (PSdif,inf,h and PSdif,eff,h)
    etalcyp3a_h ~ 0.1033685 # Suppl Table 3B row 'CLmet,h by CYP3A': lognormal, CV 33 percent; omega = log(1 + 0.33^2). IIV on hepatic CYP3A metabolism
    etalugt1a1_h ~ 0.0861777 # Suppl Table 3B row 'CLmet,h by UGT1A1': lognormal, CV 30 percent; omega = log(1 + 0.3^2). IIV on hepatic UGT1A1 glucuronidation
    etalces_h ~ 0.0861777 # Suppl Table 3B row 'CLmet,h by CES': lognormal, CV 30 percent; omega = log(1 + 0.3^2). IIV on hepatic carboxylesterase metabolism
    etalmdr1_bile ~ 0.0644423 # Suppl Table 3B row 'CLbile by MDR1': lognormal, CV 25.8 percent; omega = log(1 + 0.258^2). IIV on MDR1-mediated biliary excretion
    etalbcrp_bile ~ 0.0644423 # Suppl Table 3B row 'CLbile by BCRP1': lognormal, CV 25.8 percent; omega = log(1 + 0.258^2). IIV on BCRP-mediated biliary excretion
    etalmrp2_bile ~ 0.0644423 # Suppl Table 3B row 'CLbile by MRP2': lognormal, CV 25.8 percent; omega = log(1 + 0.258^2). IIV on MRP2-mediated biliary excretion
    etalcl_bile_other ~ 0.0644423 # Suppl Table 3B row 'CLbile of NPC, APC': lognormal, CV 25.8 percent; omega = log(1 + 0.258^2). IIV on biliary excretion of NPC and APC
    etalkbile ~ 0.0861777 # Suppl Table 3B row 'kbile': lognormal, CV 30 percent; omega = log(1 + 0.3^2). IIV on biliary transit rate constant
    etalka ~ 0.0861777 # Suppl Table 3B row 'ka': lognormal, CV 30 percent; omega = log(1 + 0.3^2). IIV on lumen-to-enterocyte absorption rate constant
    etalkfeces ~ 0.0861777 # Suppl Table 3B row 'kfeces': lognormal, CV 30 percent; omega = log(1 + 0.3^2). IIV on faecal transit rate constant
    etalmdr1_ent ~ 0.0644423 # Suppl Table 3B row 'PSact,eff,ent by MDR1': lognormal, CV 25.8 percent; omega = log(1 + 0.258^2). IIV on MDR1-mediated enterocyte apical efflux
    etalbcrp_ent ~ 0.0644423 # Suppl Table 3B row 'PSact,eff,ent by BCRP': lognormal, CV 25.8 percent; omega = log(1 + 0.258^2). IIV on BCRP-mediated enterocyte apical efflux
    etalmrp2_ent ~ 0.0644423 # Suppl Table 3B row 'PSact,eff,ent by MRP2': lognormal, CV 25.8 percent; omega = log(1 + 0.258^2). IIV on MRP2-mediated enterocyte apical efflux
    etalps_dif_ent ~ 0.0099503 # Suppl Table 3B row 'PSdif,eff,ent and PSdif,inf,ent': lognormal, CV 10 percent; omega = log(1 + 0.1^2). IIV on enterocyte basolateral passive diffusion
    etalcyp3a_ent ~ 0.1033685 # Suppl Table 3B row 'CLmet,ent by CYP3A': lognormal, CV 33 percent; omega = log(1 + 0.33^2). IIV on intestinal CYP3A metabolism
    etalugt1a1_ent ~ 0.0861777 # Suppl Table 3B row 'CLmet,ent by UGT1A1': lognormal, CV 30 percent; omega = log(1 + 0.3^2). IIV on intestinal UGT1A1 glucuronidation
    etalces_ent ~ 0.0861777 # Suppl Table 3B row 'CLmet,ent by CES': lognormal, CV 30 percent; omega = log(1 + 0.3^2). IIV on intestinal carboxylesterase metabolism

    etav_liver ~ 0.012996 # Suppl Table 3A row 'Vh': normal, CV 11.4 percent; omega = 0.114^2. IIV on liver volume
    etav_muscle ~ 0.012996 # Suppl Table 3A row 'Vmuscle': normal, CV 11.4 percent; omega = 0.114^2. IIV on muscle volume
    etav_skin ~ 0.012996 # Suppl Table 3A row 'Vskin': normal, CV 11.4 percent; omega = 0.114^2. IIV on skin volume
    etav_adipose ~ 0.012996 # Suppl Table 3A row 'Vadipose': normal, CV 11.4 percent; omega = 0.114^2. IIV on adipose volume
    etav_mucosa ~ 0.012996 # Suppl Table 3A row 'Vmucosa': normal, CV 11.4 percent; omega = 0.114^2. IIV on intestinal mucosa volume
    etav_serosa ~ 0.012996 # Suppl Table 3A row 'Vserosa': normal, CV 11.4 percent; omega = 0.114^2. IIV on intestinal serosa volume
    etaq_liver ~ 0.014884 # Suppl Table 3A row 'Qh': normal, CV 12.2 percent; omega = 0.122^2. IIV on hepatic blood flow
    etaq_muscle ~ 0.014884 # Suppl Table 3A row 'Qmuscle': normal, CV 12.2 percent; omega = 0.122^2. IIV on muscle blood flow
    etaq_skin ~ 0.014884 # Suppl Table 3A row 'Qskin': normal, CV 12.2 percent; omega = 0.122^2. IIV on skin blood flow
    etaq_adipose ~ 0.014884 # Suppl Table 3A row 'Qadipose': normal, CV 12.2 percent; omega = 0.122^2. IIV on adipose blood flow
    etaq_mucosa ~ 0.014884 # Suppl Table 3A row 'Qmucosa': normal, CV 12.2 percent; omega = 0.122^2. IIV on mucosal blood flow
    etaq_serosa ~ 0.014884 # Suppl Table 3A row 'Qserosa': normal, CV 12.2 percent; omega = 0.122^2. IIV on serosal blood flow

    # ---------------------------------------------------------------
    # 8. Residual error. The CNM minimises the weighted sum of squares of
    #    Equation 6, sum((y - yhat)^2 / y^2), a fixed-effects criterion
    #    that estimates no residual variance, and no residual-error model
    #    is reported anywhere in the paper or its supplement. nlmixr2
    #    requires a residual term per endpoint, so the five propSd values
    #    below are syntactic placeholders and must NOT be read as
    #    estimates. Same convention as Aoki_2024_bosentan_pbpk.R and
    #    Tsuchitani_2024_telmisartan_pbpk.R.
    # ---------------------------------------------------------------
    propSd <- fixed(0.1)
    label("Irinotecan proportional residual error placeholder (fraction)") # not reported in Toshimoto 2017; placeholder only
    propSd_sn38 <- fixed(0.1)
    label("SN-38 proportional residual error placeholder (fraction)") # not reported in Toshimoto 2017; placeholder only
    propSd_sn38g <- fixed(0.1)
    label("SN-38G proportional residual error placeholder (fraction)") # not reported in Toshimoto 2017; placeholder only
    propSd_npc <- fixed(0.1)
    label("NPC proportional residual error placeholder (fraction)") # not reported in Toshimoto 2017; placeholder only
    propSd_apc <- fixed(0.1)
    label("APC proportional residual error placeholder (fraction)") # not reported in Toshimoto 2017; placeholder only
  })

  model({
    # ===============================================================
    # 1. Back-transform the log-scale parameters.
    # ===============================================================
    v_central <- exp(lv_central)
    ka <- exp(lka)
    ps_dif_inf_h <- exp(lps_dif_inf_h)
    kbile <- exp(lkbile)
    kfeces <- exp(lkfeces)
    cl_sn38_h <- exp(lcl_sn38_h)
    cl_npc_h <- exp(lcl_npc_h)
    cl_apc_h <- exp(lcl_apc_h)
    cl_others_h <- exp(lcl_others_h)
    cl_bile <- exp(lcl_bile)
    r_met_ces <- exp(lr_met_ces)
    ps_dif_eff_ent <- exp(lps_dif_eff_ent)
    r_met_cyp3a <- exp(lr_met_cyp3a)
    v_central_sn38 <- exp(lv_central_sn38)
    ka_sn38 <- exp(lka_sn38)
    clint_all_sn38 <- exp(lclint_all_sn38)
    kbile_sn38 <- exp(lkbile_sn38)
    kfeces_sn38 <- exp(lkfeces_sn38)
    r_ugt_h_ent_sn38 <- exp(lr_ugt_h_ent_sn38)
    kp_muscle <- exp(lkp_muscle)
    kp_muscle_sn38 <- exp(lkp_muscle_sn38)
    kp_muscle_sn38g <- exp(lkp_muscle_sn38g)
    kp_muscle_npc <- exp(lkp_muscle_npc)
    kp_muscle_apc <- exp(lkp_muscle_apc)
    kp_skin <- exp(lkp_skin)
    kp_skin_sn38 <- exp(lkp_skin_sn38)
    kp_skin_sn38g <- exp(lkp_skin_sn38g)
    kp_skin_npc <- exp(lkp_skin_npc)
    kp_skin_apc <- exp(lkp_skin_apc)
    kp_adipose <- exp(lkp_adipose)
    kp_adipose_sn38 <- exp(lkp_adipose_sn38)
    kp_adipose_sn38g <- exp(lkp_adipose_sn38g)
    kp_adipose_npc <- exp(lkp_adipose_npc)
    kp_adipose_apc <- exp(lkp_adipose_apc)
    kp_gut <- exp(lkp_gut)
    kp_gut_sn38 <- exp(lkp_gut_sn38)
    kp_gut_sn38g <- exp(lkp_gut_sn38g)
    kp_gut_npc <- exp(lkp_gut_npc)
    kp_gut_apc <- exp(lkp_gut_apc)
    cl_renal <- exp(lcl_renal)
    cl_renal_sn38 <- exp(lcl_renal_sn38)
    cl_renal_sn38g <- exp(lcl_renal_sn38g)
    cl_renal_apc <- exp(lcl_renal_apc)
    v_central_sn38g <- exp(lv_central_sn38g)
    ka_sn38g <- exp(lka_sn38g)
    clint_all_sn38g <- exp(lclint_all_sn38g)
    kbile_sn38g <- exp(lkbile_sn38g)
    kfeces_sn38g <- exp(lkfeces_sn38g)
    kdec_sn38g <- exp(lkdec_sn38g)
    v_central_npc <- exp(lv_central_npc)
    ka_npc <- exp(lka_npc)
    clint_all_npc <- exp(lclint_all_npc)
    kbile_npc <- exp(lkbile_npc)
    kfeces_npc <- exp(lkfeces_npc)
    v_central_apc <- exp(lv_central_apc)
    ka_apc <- exp(lka_apc)
    clint_all_apc <- exp(lclint_all_apc)
    kbile_apc <- exp(lkbile_apc)
    kfeces_apc <- exp(lkfeces_apc)

    # ===============================================================
    # 2. Physiological system, scaled to the individual body weight and
    #    perturbed by the NORMAL between-subject variability of
    #    Supplementary Table 3A. Supplementary Tables 1A and 3A are
    #    per-kilogram, so every volume, flow and clearance below is
    #    multiplied by WT to give absolute L, L/h and L/h.
    # ===============================================================
    v_liver <- exp(lv_liver) * WT * (1 + etav_liver)
    v_muscle <- exp(lv_muscle) * WT * (1 + etav_muscle)
    v_skin <- exp(lv_skin) * WT * (1 + etav_skin)
    v_adipose <- exp(lv_adipose) * WT * (1 + etav_adipose)
    v_mucosa <- exp(lv_mucosa) * WT * (1 + etav_mucosa)
    v_serosa <- exp(lv_serosa) * WT * (1 + etav_serosa)
    q_liver <- exp(lq_liver) * WT * (1 + etaq_liver)
    q_muscle <- exp(lq_muscle) * WT * (1 + etaq_muscle)
    q_skin <- exp(lq_skin) * WT * (1 + etaq_skin)
    q_adipose <- exp(lq_adipose) * WT * (1 + etaq_adipose)
    q_mucosa <- exp(lq_mucosa) * WT * (1 + etaq_mucosa)
    q_serosa <- exp(lq_serosa) * WT * (1 + etaq_serosa)

    #    Liver: extracellular (sinusoidal) space versus hepatocytes, and
    #    mucosa: mucosal blood versus enterocytes (Supplementary Table 1A,
    #    "Fraction of volume"). Each of the five tandem liver units holds
    #    one fifth of the respective total volume, and receives one fifth
    #    of every liver-wide transport and elimination clearance; fdisp is
    #    that one fifth.
    fdisp <- 0.2
    v_liver_ex <- v_liver * frac_liver_ex
    v_liver_cell <- v_liver * (1 - frac_liver_ex)
    v_muc <- v_mucosa * frac_mucosa_ex
    v_ent <- v_mucosa * (1 - frac_mucosa_ex)

    # ===============================================================
    # 3. Hybrid-to-elementary parameter conversion. The CNM optimised
    #    hybrid quantities for every compound except irinotecan
    #    (Methods); the Supplementary Text "Optional equation" block and
    #    Equations 1-5 invert them. Working through them here, rather
    #    than hard-coding the elementary values, keeps every number in
    #    ini() a value that is printed in the paper.
    #
    #      Eq 1  CLint,all = (PSact,inf,h + PSdif,inf,h) *
    #                        (CLmet,h + CLbile) /
    #                        (PSdif,eff,h + CLmet,h + CLbile)
    #      Eq 2  Rdif,h    = PSdif,inf,h / PSact,inf,h
    #      Eq 3  beta      = (CLmet,h + CLbile) /
    #                        (PSdif,eff,h + CLmet,h + CLbile)
    #      Eq 4  fbile     = CLbile / (CLmet,h + CLbile)
    #      Eq 5  Rdif,ent  = PSact,eff,ent / (AR * PSdif,eff,ent)
    #            -- printed INVERTED in the article; see below.
    #
    #    so that PSact,inf,h = CLint,all / beta / (1 + Rdif,h),
    #    PSdif,inf,h = PSdif,eff,h = PSact,inf,h * Rdif,h, and
    #    CLmet,h + CLbile = PSdif,eff,h * beta / (1 - beta).
    #    Passive diffusion is assumed symmetric in both the liver and the
    #    enterocyte (Methods), i.e. PSdif,inf = PSdif,eff in each.
    # ===============================================================

    #    Irinotecan: no active hepatic uptake transporter (Supplementary
    #    Text, "PSact,inf,h = 0 for irinotecan, NPC, and APC"), and every
    #    hepatic clearance is tabulated directly.
    ps_act_inf_h_pop <- 0
    ps_dif_inf_h_pop <- ps_dif_inf_h
    ps_dif_eff_h_pop <- ps_dif_inf_h
    #    Rdif,ent: the article's typeset Equation 5 and the sentence that
    #    introduces it ("the ratio of passive diffusional efflux to the
    #    active efflux in the enterocytes") both give
    #    Rdif,ent = AR * PSdif,eff,ent / PSact,eff,ent, but the
    #    Supplementary Text "Optional equation" block gives the reciprocal,
    #    Rdif,ent = PSact,eff,ent / (AR * PSdif,eff,ent). The two on-disk
    #    sources contradict each other and the choice changes enterocyte
    #    SN-38 exposure 2.4-fold, so it is settled against two quantities
    #    the authors themselves published:
    #      (i) the ID 2 unbound plasma SN-38 AUC threshold of 26.35 nM*h
    #          (Discussion, Perspective of VCS) is reproduced to 2.6 percent
    #          by the supplement form and missed by 17 percent by Equation 5;
    #      (ii) under the supplement form a reduced-function ABCG2, ABCB1 or
    #          ABCC2 genotype RAISES enterocyte SN-38 exposure by 9-14
    #          percent, which is the direction and the magnitude needed to
    #          produce the significant efflux-transporter associations with
    #          diarrhoea in Figure 8; under Equation 5 the same genotypes
    #          LOWER it by 0.6-3.5 percent, which would falsify that figure.
    #    The supplement form is therefore what the authors ran, and is used
    #    here. See the vignette Errata.
    ps_act_eff_ent_pop <- rdif_ent * ar * ps_dif_eff_ent
    #    SAR, the common hepatic-to-intestinal ratio of passive
    #    diffusional permeability (Supplementary Text, "Optional
    #    equation"), is shared by all five compounds. It is identified
    #    from irinotecan, the only compound for which both PSdif,eff,h
    #    and PSdif,eff,ent are tabulated, and then fixes PSdif,eff,ent
    #    for the other four.
    sar <- ps_dif_eff_h_pop / ps_dif_eff_ent
    #    R met,CES and R met,CYP3A (Supplementary Text) carry the hepatic
    #    clearances to their intestinal counterparts.
    cl_sn38_ent_pop <- cl_sn38_h / r_met_ces
    cl_npc_ent_pop <- cl_npc_h / r_met_cyp3a
    cl_apc_ent_pop <- cl_apc_h / r_met_cyp3a

    #    SN-38: active OATP1B1 uptake, hepatic glucuronidation and
    #    biliary excretion. fglu is the glucuronidated share of the
    #    hepatic elimination and fbile the biliary share (Supplementary
    #    Table 4 footnote a: fglu = 1 - fbile).
    beta_sn38 <- 1 / inv_beta_sn38
    fglu_sn38 <- 1 / inv_fglu_sn38
    ps_act_inf_h_sn38_pop <- clint_all_sn38 / beta_sn38 / (1 + rdif_h_sn38)
    ps_dif_eff_h_sn38_pop <- ps_act_inf_h_sn38_pop * rdif_h_sn38
    cl_elim_h_sn38 <- ps_dif_eff_h_sn38_pop * beta_sn38 / (1 - beta_sn38)
    cl_glu_h_pop <- cl_elim_h_sn38 * fglu_sn38
    cl_bile_sn38_pop <- cl_elim_h_sn38 * (1 - fglu_sn38)
    ps_dif_eff_ent_sn38_pop <- ps_dif_eff_h_sn38_pop / sar
    ps_act_eff_ent_sn38_pop <- rdif_ent_sn38 * ar * ps_dif_eff_ent_sn38_pop
    cl_glu_ent_pop <- cl_glu_h_pop / r_ugt_h_ent_sn38

    #    SN-38G: active OATP1B1 uptake, biliary excretion only. It is not
    #    metabolised (Supplementary Text: CLmet,h = 0 for SN-38G), so
    #    fbile = 1 by Equation 4 and is not an optimised unknown -- which
    #    is why Table I leaves the 1/fbile cell for SN-38G empty.
    beta_sn38g <- 1 / inv_beta_sn38g
    ps_act_inf_h_sn38g_pop <- clint_all_sn38g / beta_sn38g / (1 + rdif_h_sn38g)
    ps_dif_eff_h_sn38g_pop <- ps_act_inf_h_sn38g_pop * rdif_h_sn38g
    cl_bile_sn38g_pop <- ps_dif_eff_h_sn38g_pop * beta_sn38g / (1 - beta_sn38g)
    ps_dif_eff_ent_sn38g_pop <- ps_dif_eff_h_sn38g_pop / sar
    ps_act_eff_ent_sn38g_pop <- rdif_ent_sn38g * ar * ps_dif_eff_ent_sn38g_pop

    #    NPC: no active uptake, so Equation 1 collapses to CLint,all =
    #    PSdif,inf,h * beta. NPC is both hydrolysed to SN-38 and excreted
    #    in bile, and 1/fbile splits its hepatic elimination between the
    #    two. See the vignette Errata: the Supplementary Text lists
    #    CLmet,h as 0 for NPC, which cannot be right because SN-38 is
    #    formed from NPC in the same document and because Table I
    #    optimises 1/fbile for NPC and for no other compound.
    beta_npc <- 1 / inv_beta_npc
    fbile_npc <- 1 / inv_fbile_npc
    ps_dif_eff_h_npc_pop <- clint_all_npc / beta_npc
    cl_elim_h_npc <- ps_dif_eff_h_npc_pop * beta_npc / (1 - beta_npc)
    cl_bile_npc_pop <- cl_elim_h_npc * fbile_npc
    cl_sn38_h_npc_pop <- cl_elim_h_npc * (1 - fbile_npc)
    cl_sn38_ent_npc_pop <- cl_sn38_h_npc_pop / r_met_ces
    ps_dif_eff_ent_npc_pop <- ps_dif_eff_h_npc_pop / sar

    #    APC: no active uptake and no metabolism, so fbile = 1.
    beta_apc <- 1 / inv_beta_apc
    ps_dif_eff_h_apc_pop <- clint_all_apc / beta_apc
    cl_bile_apc_pop <- ps_dif_eff_h_apc_pop * beta_apc / (1 - beta_apc)
    ps_dif_eff_ent_apc_pop <- ps_dif_eff_h_apc_pop / sar

    # ===============================================================
    # 4. Genotype and between-subject activity of each enzyme and
    #    transporter. A locus activity is 1 in the wild-type homozygote,
    #    the heterozygote ratio in a heterozygote and the homozygote
    #    ratio in a homozygous-variant carrier. The two SLCO1B1 loci act
    #    on the same protein and are combined multiplicatively, which is
    #    the haplotype reading -- see the vignette Errata.
    # ===============================================================
    a_ugt1a1 <- 1 + (e_ugt1a1_het_clglu - 1) * SNP_UGT1A1_RS8175347_HET +
      (e_ugt1a1_hom_clglu - 1) * SNP_UGT1A1_RS8175347_HOM
    a_oatp521 <- 1 + (e_slco1b1_521_het_psact - 1) * SNP_SLCO1B1_RS4149056_HET +
      (e_slco1b1_521_hom_psact - 1) * SNP_SLCO1B1_RS4149056_HOM
    a_oatp388 <- 1 + (e_slco1b1_388_het_psact - 1) * SNP_SLCO1B1_RS2306283_HET +
      (e_slco1b1_388_hom_psact - 1) * SNP_SLCO1B1_RS2306283_HOM
    a_abcg2 <- 1 + (e_abcg2_het_efflux - 1) * SNP_ABCG2_RS2231142_HET +
      (e_abcg2_hom_efflux - 1) * SNP_ABCG2_RS2231142_HOM
    a_abcb1 <- 1 + (e_abcb1_het_efflux - 1) * ABCB1_C3435T_HET +
      (e_abcb1_hom_efflux - 1) * ABCB1_C3435T_MUT
    a_abcc2 <- 1 + (e_abcc2_het_efflux - 1) * SNP_ABCC2_RS717620_HET +
      (e_abcc2_hom_efflux - 1) * SNP_ABCC2_RS717620_HOM

    #    Protein-level activity = genotype activity times the lognormal
    #    random effect of Supplementary Table 3B.
    act_oatp1b1 <- a_oatp521 * a_oatp388 * exp(etaloatp1b1)
    act_mdr1_bile <- a_abcb1 * exp(etalmdr1_bile)
    act_bcrp_bile <- a_abcg2 * exp(etalbcrp_bile)
    act_mrp2_bile <- a_abcc2 * exp(etalmrp2_bile)
    act_mdr1_ent <- a_abcb1 * exp(etalmdr1_ent)
    act_bcrp_ent <- a_abcg2 * exp(etalbcrp_ent)
    act_mrp2_ent <- a_abcc2 * exp(etalmrp2_ent)

    #    A clearance carried by several transporters scales as the
    #    contribution-weighted mean of their activities, with any share
    #    not attributed to a named transporter held at wild-type activity.
    scl_bile_iri <- frac_abcb1_iri * act_mdr1_bile + frac_abcc2_iri * act_mrp2_bile +
      (1 - frac_abcb1_iri - frac_abcc2_iri)
    scl_ent_iri <- frac_abcb1_iri * act_mdr1_ent + frac_abcc2_iri * act_mrp2_ent +
      (1 - frac_abcb1_iri - frac_abcc2_iri)
    scl_bile_sn38 <- frac_abcg2_sn38 * act_bcrp_bile + frac_abcb1_sn38 * act_mdr1_bile +
      frac_abcc2_sn38 * act_mrp2_bile +
      (1 - frac_abcg2_sn38 - frac_abcb1_sn38 - frac_abcc2_sn38)
    scl_ent_sn38 <- frac_abcg2_sn38 * act_bcrp_ent + frac_abcb1_sn38 * act_mdr1_ent +
      frac_abcc2_sn38 * act_mrp2_ent +
      (1 - frac_abcg2_sn38 - frac_abcb1_sn38 - frac_abcc2_sn38)
    scl_bile_sn38g <- frac_abcc2_sn38g * act_mrp2_bile + (1 - frac_abcc2_sn38g)
    scl_ent_sn38g <- frac_abcc2_sn38g * act_mrp2_ent + (1 - frac_abcc2_sn38g)

    # ===============================================================
    # 5. Individual elementary parameters. The population elementary
    #    values of section 3 are perturbed here, one random effect and
    #    one genotype activity per elementary process, exactly as the
    #    Methods describe ("Hybrid parameters ... used for CNM
    #    optimization were converted to each parameter representing
    #    elementary processes", after which the Monte Carlo draw is made
    #    per parameter). Doing the conversion first and the perturbation
    #    second matters: perturbing a hybrid would propagate one random
    #    effect into several elementary processes at once.
    #
    #    Every permeability-surface-area product and every intrinsic
    #    clearance in Supplementary Tables 1B and 4 is per kilogram
    #    (L/h/kg), so each is multiplied by WT here to match the absolute
    #    L/h blood flows of section 2. The rate constants ka, kbile,
    #    kfeces and kdec are 1/h and are not scaled.
    # ===============================================================

    #    Irinotecan.
    ps_act_inf_h <- ps_act_inf_h_pop * WT
    ps_dif_inf_h_i <- ps_dif_inf_h_pop * WT * exp(etalps_dif_h)
    ps_dif_eff_h_i <- ps_dif_inf_h_i
    cl_sn38_h_i <- cl_sn38_h * WT * exp(etalces_h)
    cl_npc_h_i <- cl_npc_h * WT * exp(etalcyp3a_h)
    cl_apc_h_i <- cl_apc_h * WT * exp(etalcyp3a_h)
    cl_met_h_iri <- cl_sn38_h_i + cl_npc_h_i + cl_apc_h_i + cl_others_h * WT
    cl_bile_i <- cl_bile * WT * scl_bile_iri
    ps_dif_eff_ent_i <- ps_dif_eff_ent * WT * exp(etalps_dif_ent)
    ps_act_eff_ent_i <- ps_act_eff_ent_pop * WT * scl_ent_iri
    cl_sn38_ent_i <- cl_sn38_ent_pop * WT * exp(etalces_ent)
    cl_npc_ent_i <- cl_npc_ent_pop * WT * exp(etalcyp3a_ent)
    cl_apc_ent_i <- cl_apc_ent_pop * WT * exp(etalcyp3a_ent)
    cl_met_ent_iri <- cl_sn38_ent_i + cl_npc_ent_i + cl_apc_ent_i + cl_others_ent * WT
    ka_i <- ka * exp(etalka)
    kbile_i <- kbile * exp(etalkbile)
    kfeces_i <- kfeces * exp(etalkfeces)
    v_central_i <- v_central * WT

    #    SN-38.
    ps_act_inf_h_sn38 <- ps_act_inf_h_sn38_pop * WT * act_oatp1b1
    ps_dif_eff_h_sn38 <- ps_dif_eff_h_sn38_pop * WT * exp(etalps_dif_h)
    cl_glu_h <- cl_glu_h_pop * WT * a_ugt1a1 * exp(etalugt1a1_h)
    cl_bile_sn38_i <- cl_bile_sn38_pop * WT * scl_bile_sn38
    ps_dif_eff_ent_sn38 <- ps_dif_eff_ent_sn38_pop * WT * exp(etalps_dif_ent)
    ps_act_eff_ent_sn38 <- ps_act_eff_ent_sn38_pop * WT * scl_ent_sn38
    cl_glu_ent <- cl_glu_ent_pop * WT * a_ugt1a1 * exp(etalugt1a1_ent)
    ka_sn38_i <- ka_sn38 * exp(etalka)
    kbile_sn38_i <- kbile_sn38 * exp(etalkbile)
    kfeces_sn38_i <- kfeces_sn38 * exp(etalkfeces)
    v_central_sn38_i <- v_central_sn38 * WT

    #    SN-38G.
    ps_act_inf_h_sn38g <- ps_act_inf_h_sn38g_pop * WT * act_oatp1b1
    ps_dif_eff_h_sn38g <- ps_dif_eff_h_sn38g_pop * WT * exp(etalps_dif_h)
    cl_bile_sn38g_i <- cl_bile_sn38g_pop * WT * scl_bile_sn38g
    ps_dif_eff_ent_sn38g <- ps_dif_eff_ent_sn38g_pop * WT * exp(etalps_dif_ent)
    ps_act_eff_ent_sn38g <- ps_act_eff_ent_sn38g_pop * WT * scl_ent_sn38g
    ka_sn38g_i <- ka_sn38g * exp(etalka)
    kbile_sn38g_i <- kbile_sn38g * exp(etalkbile)
    kfeces_sn38g_i <- kfeces_sn38g * exp(etalkfeces)
    kdec <- kdec_sn38g
    v_central_sn38g_i <- v_central_sn38g * WT

    #    NPC.
    ps_dif_eff_h_npc <- ps_dif_eff_h_npc_pop * WT * exp(etalps_dif_h)
    cl_sn38_h_npc <- cl_sn38_h_npc_pop * WT * exp(etalces_h)
    cl_bile_npc_i <- cl_bile_npc_pop * WT * exp(etalcl_bile_other)
    ps_dif_eff_ent_npc <- ps_dif_eff_ent_npc_pop * WT * exp(etalps_dif_ent)
    cl_sn38_ent_npc <- cl_sn38_ent_npc_pop * WT * exp(etalces_ent)
    ka_npc_i <- ka_npc * exp(etalka)
    kbile_npc_i <- kbile_npc * exp(etalkbile)
    kfeces_npc_i <- kfeces_npc * exp(etalkfeces)
    v_central_npc_i <- v_central_npc * WT

    #    APC.
    ps_dif_eff_h_apc <- ps_dif_eff_h_apc_pop * WT * exp(etalps_dif_h)
    cl_bile_apc_i <- cl_bile_apc_pop * WT * exp(etalcl_bile_other)
    ps_dif_eff_ent_apc <- ps_dif_eff_ent_apc_pop * WT * exp(etalps_dif_ent)
    ka_apc_i <- ka_apc * exp(etalka)
    kbile_apc_i <- kbile_apc * exp(etalkbile)
    kfeces_apc_i <- kfeces_apc * exp(etalkfeces)
    v_central_apc_i <- v_central_apc * WT

    #    Renal clearances are per kilogram in Supplementary Table 1B and
    #    carry no random effect: the only renal variability the paper
    #    reports is on CLint,sec (Supplementary Table 3B, CV 34.2 percent),
    #    and mapping that intrinsic quantity onto CLr needs the renal
    #    dispersion model of Equations 12-16, whose Qrtb and ERPF inputs
    #    are not tabulated anywhere on disk. See the vignette Errata.
    cl_renal_i <- cl_renal * WT
    cl_renal_sn38_i <- cl_renal_sn38 * WT
    cl_renal_sn38g_i <- cl_renal_sn38g * WT
    cl_renal_npc_i <- cl_renal_npc * WT
    cl_renal_apc_i <- cl_renal_apc * WT

    # ===============================================================
    # 6. ODE system -- Supplementary Text, "Ordinary differential
    #    equations", written out once per chemical species. The tissue
    #    states hold CONCENTRATIONS (umol/L); gut_lumen, ehc1-3, a_feces
    #    and a_urine hold AMOUNTS (umol), exactly as the published
    #    equations do.
    #
    #    Three departures from the published equations are made, all of
    #    them forced by mass balance and all documented in the vignette
    #    Errata:
    #      (a) the hepatic-extracellular equations print the uptake term
    #          "+ fb * ((PSact,inf,h + PSdif,inf,h)/5) * CHE,i"; the sign
    #          must be negative, because the identical term appears with
    #          a plus in the hepatocyte equation and drug cannot be
    #          created in both compartments at once;
    #      (b) the metabolite formation terms Xi(t) and Z(t) are printed
    #          without the unbound fraction and, in the liver, without
    #          the one-fifth dispersion weight that the matching loss
    #          term in the parent equation carries; both are restored
    #          here so that what leaves the parent is exactly what
    #          arrives in the metabolite;
    #      (c) CLmet,h is printed as 0 for NPC although SN-38 is formed
    #          from NPC elsewhere in the same document; it is taken as
    #          CL SN-38,h(NPC), which is what Table I optimising 1/fbile
    #          for NPC alone implies.
    # ===============================================================

    # ---------------- Irinotecan ----------------
    d/dt(central) <- (q_liver * is_liver5 - q_liver * central +
      q_muscle * (muscle / kp_muscle - central) +
      q_skin * (skin / kp_skin - central) +
      q_adipose * (adipose / kp_adipose - central) -
      cl_renal_i * central) / v_central_i
    d/dt(muscle) <- q_muscle * (central - muscle / kp_muscle) / v_muscle
    d/dt(skin) <- q_skin * (central - skin / kp_skin) / v_skin
    d/dt(adipose) <- q_adipose * (central - adipose / kp_adipose) / v_adipose
    d/dt(serosa) <- q_serosa * (central - serosa / kp_gut) / v_serosa

    #    Hepatic extracellular (sinusoidal) chain. Unit 1 receives the
    #    hepatic-arterial inflow plus the serosal and mucosal venous
    #    return; units 2-5 receive their predecessor's outflow.
    d/dt(is_liver1) <- ((q_liver - q_serosa - q_mucosa) * central +
      q_serosa * serosa / kp_gut + q_mucosa * intestine_muc -
      q_liver * is_liver1 -
      fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver1 +
      fh * fdisp * ps_dif_eff_h_i * int_liver1) / (v_liver_ex * fdisp)
    d/dt(is_liver2) <- (q_liver * (is_liver1 - is_liver2) -
      fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver2 +
      fh * fdisp * ps_dif_eff_h_i * int_liver2) / (v_liver_ex * fdisp)
    d/dt(is_liver3) <- (q_liver * (is_liver2 - is_liver3) -
      fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver3 +
      fh * fdisp * ps_dif_eff_h_i * int_liver3) / (v_liver_ex * fdisp)
    d/dt(is_liver4) <- (q_liver * (is_liver3 - is_liver4) -
      fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver4 +
      fh * fdisp * ps_dif_eff_h_i * int_liver4) / (v_liver_ex * fdisp)
    d/dt(is_liver5) <- (q_liver * (is_liver4 - is_liver5) -
      fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver5 +
      fh * fdisp * ps_dif_eff_h_i * int_liver5) / (v_liver_ex * fdisp)

    #    Hepatocyte chain: uptake in, passive efflux plus metabolism
    #    plus biliary excretion out, formation from the precursor in.
    d/dt(int_liver1) <- (fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver1 -
      fh * fdisp * (ps_dif_eff_h_i + cl_met_h_iri + cl_bile_i) * int_liver1) / (v_liver_cell * fdisp)
    d/dt(int_liver2) <- (fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver2 -
      fh * fdisp * (ps_dif_eff_h_i + cl_met_h_iri + cl_bile_i) * int_liver2) / (v_liver_cell * fdisp)
    d/dt(int_liver3) <- (fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver3 -
      fh * fdisp * (ps_dif_eff_h_i + cl_met_h_iri + cl_bile_i) * int_liver3) / (v_liver_cell * fdisp)
    d/dt(int_liver4) <- (fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver4 -
      fh * fdisp * (ps_dif_eff_h_i + cl_met_h_iri + cl_bile_i) * int_liver4) / (v_liver_cell * fdisp)
    d/dt(int_liver5) <- (fb * fdisp * (ps_act_inf_h + ps_dif_eff_h_i) * is_liver5 -
      fh * fdisp * (ps_dif_eff_h_i + cl_met_h_iri + cl_bile_i) * int_liver5) / (v_liver_cell * fdisp)

    #    Biliary transit chain (amounts). ehc1 collects the biliary
    #    efflux of all five hepatocyte units.
    d/dt(ehc1) <- fh * fdisp * cl_bile_i * (int_liver1 + int_liver2 +
      int_liver3 + int_liver4 + int_liver5) - kbile_i * ehc1
    d/dt(ehc2) <- kbile_i * (ehc1 - ehc2)
    d/dt(ehc3) <- kbile_i * (ehc2 - ehc3)

    #    Intestinal lumen (amount): bile in, apical efflux from the
    #    enterocyte in, absorption and faecal transit out.
    d/dt(gut_lumen) <- kbile_i * ehc3 +
      fgut * (ar * ps_dif_eff_ent_i + ps_act_eff_ent_i) * intestine_ent -
      (ka_i + kfeces_i) * gut_lumen

    #    Enterocyte (concentration) and mucosal blood (concentration).
    d/dt(intestine_ent) <- (ka_i * gut_lumen +
      fb * ps_dif_eff_ent_i * intestine_muc -
      fgut * ((ar + 1) * ps_dif_eff_ent_i + cl_met_ent_iri + ps_act_eff_ent_i) *
      intestine_ent) / v_ent
    d/dt(intestine_muc) <- (q_mucosa * central +
      fgut * ps_dif_eff_ent_i * intestine_ent - q_mucosa * intestine_muc -
      fb * ps_dif_eff_ent_i * intestine_muc) / v_muc

    #    Excretion sinks (amounts).
    d/dt(a_feces) <- kfeces_i * gut_lumen
    d/dt(a_urine) <- cl_renal_i * central

    # ---------------- SN-38 ----------------
    d/dt(central_sn38) <- (q_liver * is_liver5_sn38 - q_liver * central_sn38 +
      q_muscle * (muscle_sn38 / kp_muscle_sn38 - central_sn38) +
      q_skin * (skin_sn38 / kp_skin_sn38 - central_sn38) +
      q_adipose * (adipose_sn38 / kp_adipose_sn38 - central_sn38) -
      cl_renal_sn38_i * central_sn38) / v_central_sn38_i
    d/dt(muscle_sn38) <- q_muscle * (central_sn38 - muscle_sn38 / kp_muscle_sn38) / v_muscle
    d/dt(skin_sn38) <- q_skin * (central_sn38 - skin_sn38 / kp_skin_sn38) / v_skin
    d/dt(adipose_sn38) <- q_adipose * (central_sn38 - adipose_sn38 / kp_adipose_sn38) / v_adipose
    d/dt(serosa_sn38) <- q_serosa * (central_sn38 - serosa_sn38 / kp_gut_sn38) / v_serosa

    #    Hepatic extracellular (sinusoidal) chain. Unit 1 receives the
    #    hepatic-arterial inflow plus the serosal and mucosal venous
    #    return; units 2-5 receive their predecessor's outflow.
    d/dt(is_liver1_sn38) <- ((q_liver - q_serosa - q_mucosa) * central_sn38 +
      q_serosa * serosa_sn38 / kp_gut_sn38 + q_mucosa * intestine_muc_sn38 -
      q_liver * is_liver1_sn38 -
      fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver1_sn38 +
      fh_sn38 * fdisp * ps_dif_eff_h_sn38 * int_liver1_sn38) / (v_liver_ex * fdisp)
    d/dt(is_liver2_sn38) <- (q_liver * (is_liver1_sn38 - is_liver2_sn38) -
      fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver2_sn38 +
      fh_sn38 * fdisp * ps_dif_eff_h_sn38 * int_liver2_sn38) / (v_liver_ex * fdisp)
    d/dt(is_liver3_sn38) <- (q_liver * (is_liver2_sn38 - is_liver3_sn38) -
      fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver3_sn38 +
      fh_sn38 * fdisp * ps_dif_eff_h_sn38 * int_liver3_sn38) / (v_liver_ex * fdisp)
    d/dt(is_liver4_sn38) <- (q_liver * (is_liver3_sn38 - is_liver4_sn38) -
      fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver4_sn38 +
      fh_sn38 * fdisp * ps_dif_eff_h_sn38 * int_liver4_sn38) / (v_liver_ex * fdisp)
    d/dt(is_liver5_sn38) <- (q_liver * (is_liver4_sn38 - is_liver5_sn38) -
      fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver5_sn38 +
      fh_sn38 * fdisp * ps_dif_eff_h_sn38 * int_liver5_sn38) / (v_liver_ex * fdisp)

    #    Hepatocyte chain: uptake in, passive efflux plus metabolism
    #    plus biliary excretion out, formation from the precursor in.
    d/dt(int_liver1_sn38) <- (fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver1_sn38 -
      fh_sn38 * fdisp * (ps_dif_eff_h_sn38 + cl_glu_h + cl_bile_sn38_i) * int_liver1_sn38 +
      fh * fdisp * cl_sn38_h_i * int_liver1 +
      fh_npc * fdisp * cl_sn38_h_npc * int_liver1_npc) / (v_liver_cell * fdisp)
    d/dt(int_liver2_sn38) <- (fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver2_sn38 -
      fh_sn38 * fdisp * (ps_dif_eff_h_sn38 + cl_glu_h + cl_bile_sn38_i) * int_liver2_sn38 +
      fh * fdisp * cl_sn38_h_i * int_liver2 +
      fh_npc * fdisp * cl_sn38_h_npc * int_liver2_npc) / (v_liver_cell * fdisp)
    d/dt(int_liver3_sn38) <- (fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver3_sn38 -
      fh_sn38 * fdisp * (ps_dif_eff_h_sn38 + cl_glu_h + cl_bile_sn38_i) * int_liver3_sn38 +
      fh * fdisp * cl_sn38_h_i * int_liver3 +
      fh_npc * fdisp * cl_sn38_h_npc * int_liver3_npc) / (v_liver_cell * fdisp)
    d/dt(int_liver4_sn38) <- (fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver4_sn38 -
      fh_sn38 * fdisp * (ps_dif_eff_h_sn38 + cl_glu_h + cl_bile_sn38_i) * int_liver4_sn38 +
      fh * fdisp * cl_sn38_h_i * int_liver4 +
      fh_npc * fdisp * cl_sn38_h_npc * int_liver4_npc) / (v_liver_cell * fdisp)
    d/dt(int_liver5_sn38) <- (fb_sn38 * fdisp * (ps_act_inf_h_sn38 + ps_dif_eff_h_sn38) * is_liver5_sn38 -
      fh_sn38 * fdisp * (ps_dif_eff_h_sn38 + cl_glu_h + cl_bile_sn38_i) * int_liver5_sn38 +
      fh * fdisp * cl_sn38_h_i * int_liver5 +
      fh_npc * fdisp * cl_sn38_h_npc * int_liver5_npc) / (v_liver_cell * fdisp)

    #    Biliary transit chain (amounts). ehc1 collects the biliary
    #    efflux of all five hepatocyte units.
    d/dt(ehc1_sn38) <- fh_sn38 * fdisp * cl_bile_sn38_i * (int_liver1_sn38 + int_liver2_sn38 +
      int_liver3_sn38 + int_liver4_sn38 + int_liver5_sn38) - kbile_sn38_i * ehc1_sn38
    d/dt(ehc2_sn38) <- kbile_sn38_i * (ehc1_sn38 - ehc2_sn38)
    d/dt(ehc3_sn38) <- kbile_sn38_i * (ehc2_sn38 - ehc3_sn38)

    #    Intestinal lumen (amount): bile in, apical efflux from the
    #    enterocyte in, absorption and faecal transit out.
    #    SN-38 additionally gains the luminal deconjugation of SN-38G.
    d/dt(gut_lumen_sn38) <- kbile_sn38_i * ehc3_sn38 +
      fgut_sn38 * (ar * ps_dif_eff_ent_sn38 + ps_act_eff_ent_sn38) * intestine_ent_sn38 -
      (ka_sn38_i + kfeces_sn38_i) * gut_lumen_sn38 +
      kdec * gut_lumen_sn38g

    #    Enterocyte (concentration) and mucosal blood (concentration).
    d/dt(intestine_ent_sn38) <- (ka_sn38_i * gut_lumen_sn38 +
      fb_sn38 * ps_dif_eff_ent_sn38 * intestine_muc_sn38 -
      fgut_sn38 * ((ar + 1) * ps_dif_eff_ent_sn38 + cl_glu_ent + ps_act_eff_ent_sn38) *
      intestine_ent_sn38 +
      fgut * cl_sn38_ent_i * intestine_ent +
      fgut_npc * cl_sn38_ent_npc * intestine_ent_npc) / v_ent
    d/dt(intestine_muc_sn38) <- (q_mucosa * central_sn38 +
      fgut_sn38 * ps_dif_eff_ent_sn38 * intestine_ent_sn38 - q_mucosa * intestine_muc_sn38 -
      fb_sn38 * ps_dif_eff_ent_sn38 * intestine_muc_sn38) / v_muc

    #    Excretion sinks (amounts).
    d/dt(a_feces_sn38) <- kfeces_sn38_i * gut_lumen_sn38
    d/dt(a_urine_sn38) <- cl_renal_sn38_i * central_sn38

    # ---------------- SN-38G ----------------
    d/dt(central_sn38g) <- (q_liver * is_liver5_sn38g - q_liver * central_sn38g +
      q_muscle * (muscle_sn38g / kp_muscle_sn38g - central_sn38g) +
      q_skin * (skin_sn38g / kp_skin_sn38g - central_sn38g) +
      q_adipose * (adipose_sn38g / kp_adipose_sn38g - central_sn38g) -
      cl_renal_sn38g_i * central_sn38g) / v_central_sn38g_i
    d/dt(muscle_sn38g) <- q_muscle * (central_sn38g - muscle_sn38g / kp_muscle_sn38g) / v_muscle
    d/dt(skin_sn38g) <- q_skin * (central_sn38g - skin_sn38g / kp_skin_sn38g) / v_skin
    d/dt(adipose_sn38g) <- q_adipose * (central_sn38g - adipose_sn38g / kp_adipose_sn38g) / v_adipose
    d/dt(serosa_sn38g) <- q_serosa * (central_sn38g - serosa_sn38g / kp_gut_sn38g) / v_serosa

    #    Hepatic extracellular (sinusoidal) chain. Unit 1 receives the
    #    hepatic-arterial inflow plus the serosal and mucosal venous
    #    return; units 2-5 receive their predecessor's outflow.
    d/dt(is_liver1_sn38g) <- ((q_liver - q_serosa - q_mucosa) * central_sn38g +
      q_serosa * serosa_sn38g / kp_gut_sn38g + q_mucosa * intestine_muc_sn38g -
      q_liver * is_liver1_sn38g -
      fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver1_sn38g +
      fh_sn38g * fdisp * ps_dif_eff_h_sn38g * int_liver1_sn38g) / (v_liver_ex * fdisp)
    d/dt(is_liver2_sn38g) <- (q_liver * (is_liver1_sn38g - is_liver2_sn38g) -
      fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver2_sn38g +
      fh_sn38g * fdisp * ps_dif_eff_h_sn38g * int_liver2_sn38g) / (v_liver_ex * fdisp)
    d/dt(is_liver3_sn38g) <- (q_liver * (is_liver2_sn38g - is_liver3_sn38g) -
      fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver3_sn38g +
      fh_sn38g * fdisp * ps_dif_eff_h_sn38g * int_liver3_sn38g) / (v_liver_ex * fdisp)
    d/dt(is_liver4_sn38g) <- (q_liver * (is_liver3_sn38g - is_liver4_sn38g) -
      fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver4_sn38g +
      fh_sn38g * fdisp * ps_dif_eff_h_sn38g * int_liver4_sn38g) / (v_liver_ex * fdisp)
    d/dt(is_liver5_sn38g) <- (q_liver * (is_liver4_sn38g - is_liver5_sn38g) -
      fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver5_sn38g +
      fh_sn38g * fdisp * ps_dif_eff_h_sn38g * int_liver5_sn38g) / (v_liver_ex * fdisp)

    #    Hepatocyte chain: uptake in, passive efflux plus metabolism
    #    plus biliary excretion out, formation from the precursor in.
    d/dt(int_liver1_sn38g) <- (fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver1_sn38g -
      fh_sn38g * fdisp * (ps_dif_eff_h_sn38g + 0 + cl_bile_sn38g_i) * int_liver1_sn38g +
      fh_sn38 * fdisp * cl_glu_h * int_liver1_sn38) / (v_liver_cell * fdisp)
    d/dt(int_liver2_sn38g) <- (fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver2_sn38g -
      fh_sn38g * fdisp * (ps_dif_eff_h_sn38g + 0 + cl_bile_sn38g_i) * int_liver2_sn38g +
      fh_sn38 * fdisp * cl_glu_h * int_liver2_sn38) / (v_liver_cell * fdisp)
    d/dt(int_liver3_sn38g) <- (fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver3_sn38g -
      fh_sn38g * fdisp * (ps_dif_eff_h_sn38g + 0 + cl_bile_sn38g_i) * int_liver3_sn38g +
      fh_sn38 * fdisp * cl_glu_h * int_liver3_sn38) / (v_liver_cell * fdisp)
    d/dt(int_liver4_sn38g) <- (fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver4_sn38g -
      fh_sn38g * fdisp * (ps_dif_eff_h_sn38g + 0 + cl_bile_sn38g_i) * int_liver4_sn38g +
      fh_sn38 * fdisp * cl_glu_h * int_liver4_sn38) / (v_liver_cell * fdisp)
    d/dt(int_liver5_sn38g) <- (fb_sn38g * fdisp * (ps_act_inf_h_sn38g + ps_dif_eff_h_sn38g) * is_liver5_sn38g -
      fh_sn38g * fdisp * (ps_dif_eff_h_sn38g + 0 + cl_bile_sn38g_i) * int_liver5_sn38g +
      fh_sn38 * fdisp * cl_glu_h * int_liver5_sn38) / (v_liver_cell * fdisp)

    #    Biliary transit chain (amounts). ehc1 collects the biliary
    #    efflux of all five hepatocyte units.
    d/dt(ehc1_sn38g) <- fh_sn38g * fdisp * cl_bile_sn38g_i * (int_liver1_sn38g + int_liver2_sn38g +
      int_liver3_sn38g + int_liver4_sn38g + int_liver5_sn38g) - kbile_sn38g_i * ehc1_sn38g
    d/dt(ehc2_sn38g) <- kbile_sn38g_i * (ehc1_sn38g - ehc2_sn38g)
    d/dt(ehc3_sn38g) <- kbile_sn38g_i * (ehc2_sn38g - ehc3_sn38g)

    #    Intestinal lumen (amount): bile in, apical efflux from the
    #    enterocyte in, absorption and faecal transit out.
    #    SN-38G loses that same deconjugation flux.
    d/dt(gut_lumen_sn38g) <- kbile_sn38g_i * ehc3_sn38g +
      fgut_sn38g * (ar * ps_dif_eff_ent_sn38g + ps_act_eff_ent_sn38g) * intestine_ent_sn38g -
      (ka_sn38g_i + kfeces_sn38g_i) * gut_lumen_sn38g -
      kdec * gut_lumen_sn38g

    #    Enterocyte (concentration) and mucosal blood (concentration).
    d/dt(intestine_ent_sn38g) <- (ka_sn38g_i * gut_lumen_sn38g +
      fb_sn38g * ps_dif_eff_ent_sn38g * intestine_muc_sn38g -
      fgut_sn38g * ((ar + 1) * ps_dif_eff_ent_sn38g + 0 + ps_act_eff_ent_sn38g) *
      intestine_ent_sn38g +
      fgut_sn38 * cl_glu_ent * intestine_ent_sn38) / v_ent
    d/dt(intestine_muc_sn38g) <- (q_mucosa * central_sn38g +
      fgut_sn38g * ps_dif_eff_ent_sn38g * intestine_ent_sn38g - q_mucosa * intestine_muc_sn38g -
      fb_sn38g * ps_dif_eff_ent_sn38g * intestine_muc_sn38g) / v_muc

    #    Excretion sinks (amounts).
    d/dt(a_feces_sn38g) <- kfeces_sn38g_i * gut_lumen_sn38g
    d/dt(a_urine_sn38g) <- cl_renal_sn38g_i * central_sn38g

    # ---------------- NPC ----------------
    d/dt(central_npc) <- (q_liver * is_liver5_npc - q_liver * central_npc +
      q_muscle * (muscle_npc / kp_muscle_npc - central_npc) +
      q_skin * (skin_npc / kp_skin_npc - central_npc) +
      q_adipose * (adipose_npc / kp_adipose_npc - central_npc) -
      cl_renal_npc_i * central_npc) / v_central_npc_i
    d/dt(muscle_npc) <- q_muscle * (central_npc - muscle_npc / kp_muscle_npc) / v_muscle
    d/dt(skin_npc) <- q_skin * (central_npc - skin_npc / kp_skin_npc) / v_skin
    d/dt(adipose_npc) <- q_adipose * (central_npc - adipose_npc / kp_adipose_npc) / v_adipose
    d/dt(serosa_npc) <- q_serosa * (central_npc - serosa_npc / kp_gut_npc) / v_serosa

    #    Hepatic extracellular (sinusoidal) chain. Unit 1 receives the
    #    hepatic-arterial inflow plus the serosal and mucosal venous
    #    return; units 2-5 receive their predecessor's outflow.
    d/dt(is_liver1_npc) <- ((q_liver - q_serosa - q_mucosa) * central_npc +
      q_serosa * serosa_npc / kp_gut_npc + q_mucosa * intestine_muc_npc -
      q_liver * is_liver1_npc -
      fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver1_npc +
      fh_npc * fdisp * ps_dif_eff_h_npc * int_liver1_npc) / (v_liver_ex * fdisp)
    d/dt(is_liver2_npc) <- (q_liver * (is_liver1_npc - is_liver2_npc) -
      fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver2_npc +
      fh_npc * fdisp * ps_dif_eff_h_npc * int_liver2_npc) / (v_liver_ex * fdisp)
    d/dt(is_liver3_npc) <- (q_liver * (is_liver2_npc - is_liver3_npc) -
      fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver3_npc +
      fh_npc * fdisp * ps_dif_eff_h_npc * int_liver3_npc) / (v_liver_ex * fdisp)
    d/dt(is_liver4_npc) <- (q_liver * (is_liver3_npc - is_liver4_npc) -
      fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver4_npc +
      fh_npc * fdisp * ps_dif_eff_h_npc * int_liver4_npc) / (v_liver_ex * fdisp)
    d/dt(is_liver5_npc) <- (q_liver * (is_liver4_npc - is_liver5_npc) -
      fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver5_npc +
      fh_npc * fdisp * ps_dif_eff_h_npc * int_liver5_npc) / (v_liver_ex * fdisp)

    #    Hepatocyte chain: uptake in, passive efflux plus metabolism
    #    plus biliary excretion out, formation from the precursor in.
    d/dt(int_liver1_npc) <- (fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver1_npc -
      fh_npc * fdisp * (ps_dif_eff_h_npc + cl_sn38_h_npc + cl_bile_npc_i) * int_liver1_npc +
      fh * fdisp * cl_npc_h_i * int_liver1) / (v_liver_cell * fdisp)
    d/dt(int_liver2_npc) <- (fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver2_npc -
      fh_npc * fdisp * (ps_dif_eff_h_npc + cl_sn38_h_npc + cl_bile_npc_i) * int_liver2_npc +
      fh * fdisp * cl_npc_h_i * int_liver2) / (v_liver_cell * fdisp)
    d/dt(int_liver3_npc) <- (fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver3_npc -
      fh_npc * fdisp * (ps_dif_eff_h_npc + cl_sn38_h_npc + cl_bile_npc_i) * int_liver3_npc +
      fh * fdisp * cl_npc_h_i * int_liver3) / (v_liver_cell * fdisp)
    d/dt(int_liver4_npc) <- (fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver4_npc -
      fh_npc * fdisp * (ps_dif_eff_h_npc + cl_sn38_h_npc + cl_bile_npc_i) * int_liver4_npc +
      fh * fdisp * cl_npc_h_i * int_liver4) / (v_liver_cell * fdisp)
    d/dt(int_liver5_npc) <- (fb_npc * fdisp * (0 + ps_dif_eff_h_npc) * is_liver5_npc -
      fh_npc * fdisp * (ps_dif_eff_h_npc + cl_sn38_h_npc + cl_bile_npc_i) * int_liver5_npc +
      fh * fdisp * cl_npc_h_i * int_liver5) / (v_liver_cell * fdisp)

    #    Biliary transit chain (amounts). ehc1 collects the biliary
    #    efflux of all five hepatocyte units.
    d/dt(ehc1_npc) <- fh_npc * fdisp * cl_bile_npc_i * (int_liver1_npc + int_liver2_npc +
      int_liver3_npc + int_liver4_npc + int_liver5_npc) - kbile_npc_i * ehc1_npc
    d/dt(ehc2_npc) <- kbile_npc_i * (ehc1_npc - ehc2_npc)
    d/dt(ehc3_npc) <- kbile_npc_i * (ehc2_npc - ehc3_npc)

    #    Intestinal lumen (amount): bile in, apical efflux from the
    #    enterocyte in, absorption and faecal transit out.
    d/dt(gut_lumen_npc) <- kbile_npc_i * ehc3_npc +
      fgut_npc * (ar * ps_dif_eff_ent_npc + 0) * intestine_ent_npc -
      (ka_npc_i + kfeces_npc_i) * gut_lumen_npc

    #    Enterocyte (concentration) and mucosal blood (concentration).
    d/dt(intestine_ent_npc) <- (ka_npc_i * gut_lumen_npc +
      fb_npc * ps_dif_eff_ent_npc * intestine_muc_npc -
      fgut_npc * ((ar + 1) * ps_dif_eff_ent_npc + cl_sn38_ent_npc + 0) *
      intestine_ent_npc +
      fgut * cl_npc_ent_i * intestine_ent) / v_ent
    d/dt(intestine_muc_npc) <- (q_mucosa * central_npc +
      fgut_npc * ps_dif_eff_ent_npc * intestine_ent_npc - q_mucosa * intestine_muc_npc -
      fb_npc * ps_dif_eff_ent_npc * intestine_muc_npc) / v_muc

    #    Excretion sinks (amounts).
    d/dt(a_feces_npc) <- kfeces_npc_i * gut_lumen_npc
    d/dt(a_urine_npc) <- cl_renal_npc_i * central_npc

    # ---------------- APC ----------------
    d/dt(central_apc) <- (q_liver * is_liver5_apc - q_liver * central_apc +
      q_muscle * (muscle_apc / kp_muscle_apc - central_apc) +
      q_skin * (skin_apc / kp_skin_apc - central_apc) +
      q_adipose * (adipose_apc / kp_adipose_apc - central_apc) -
      cl_renal_apc_i * central_apc) / v_central_apc_i
    d/dt(muscle_apc) <- q_muscle * (central_apc - muscle_apc / kp_muscle_apc) / v_muscle
    d/dt(skin_apc) <- q_skin * (central_apc - skin_apc / kp_skin_apc) / v_skin
    d/dt(adipose_apc) <- q_adipose * (central_apc - adipose_apc / kp_adipose_apc) / v_adipose
    d/dt(serosa_apc) <- q_serosa * (central_apc - serosa_apc / kp_gut_apc) / v_serosa

    #    Hepatic extracellular (sinusoidal) chain. Unit 1 receives the
    #    hepatic-arterial inflow plus the serosal and mucosal venous
    #    return; units 2-5 receive their predecessor's outflow.
    d/dt(is_liver1_apc) <- ((q_liver - q_serosa - q_mucosa) * central_apc +
      q_serosa * serosa_apc / kp_gut_apc + q_mucosa * intestine_muc_apc -
      q_liver * is_liver1_apc -
      fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver1_apc +
      fh_apc * fdisp * ps_dif_eff_h_apc * int_liver1_apc) / (v_liver_ex * fdisp)
    d/dt(is_liver2_apc) <- (q_liver * (is_liver1_apc - is_liver2_apc) -
      fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver2_apc +
      fh_apc * fdisp * ps_dif_eff_h_apc * int_liver2_apc) / (v_liver_ex * fdisp)
    d/dt(is_liver3_apc) <- (q_liver * (is_liver2_apc - is_liver3_apc) -
      fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver3_apc +
      fh_apc * fdisp * ps_dif_eff_h_apc * int_liver3_apc) / (v_liver_ex * fdisp)
    d/dt(is_liver4_apc) <- (q_liver * (is_liver3_apc - is_liver4_apc) -
      fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver4_apc +
      fh_apc * fdisp * ps_dif_eff_h_apc * int_liver4_apc) / (v_liver_ex * fdisp)
    d/dt(is_liver5_apc) <- (q_liver * (is_liver4_apc - is_liver5_apc) -
      fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver5_apc +
      fh_apc * fdisp * ps_dif_eff_h_apc * int_liver5_apc) / (v_liver_ex * fdisp)

    #    Hepatocyte chain: uptake in, passive efflux plus metabolism
    #    plus biliary excretion out, formation from the precursor in.
    d/dt(int_liver1_apc) <- (fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver1_apc -
      fh_apc * fdisp * (ps_dif_eff_h_apc + 0 + cl_bile_apc_i) * int_liver1_apc +
      fh * fdisp * cl_apc_h_i * int_liver1) / (v_liver_cell * fdisp)
    d/dt(int_liver2_apc) <- (fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver2_apc -
      fh_apc * fdisp * (ps_dif_eff_h_apc + 0 + cl_bile_apc_i) * int_liver2_apc +
      fh * fdisp * cl_apc_h_i * int_liver2) / (v_liver_cell * fdisp)
    d/dt(int_liver3_apc) <- (fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver3_apc -
      fh_apc * fdisp * (ps_dif_eff_h_apc + 0 + cl_bile_apc_i) * int_liver3_apc +
      fh * fdisp * cl_apc_h_i * int_liver3) / (v_liver_cell * fdisp)
    d/dt(int_liver4_apc) <- (fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver4_apc -
      fh_apc * fdisp * (ps_dif_eff_h_apc + 0 + cl_bile_apc_i) * int_liver4_apc +
      fh * fdisp * cl_apc_h_i * int_liver4) / (v_liver_cell * fdisp)
    d/dt(int_liver5_apc) <- (fb_apc * fdisp * (0 + ps_dif_eff_h_apc) * is_liver5_apc -
      fh_apc * fdisp * (ps_dif_eff_h_apc + 0 + cl_bile_apc_i) * int_liver5_apc +
      fh * fdisp * cl_apc_h_i * int_liver5) / (v_liver_cell * fdisp)

    #    Biliary transit chain (amounts). ehc1 collects the biliary
    #    efflux of all five hepatocyte units.
    d/dt(ehc1_apc) <- fh_apc * fdisp * cl_bile_apc_i * (int_liver1_apc + int_liver2_apc +
      int_liver3_apc + int_liver4_apc + int_liver5_apc) - kbile_apc_i * ehc1_apc
    d/dt(ehc2_apc) <- kbile_apc_i * (ehc1_apc - ehc2_apc)
    d/dt(ehc3_apc) <- kbile_apc_i * (ehc2_apc - ehc3_apc)

    #    Intestinal lumen (amount): bile in, apical efflux from the
    #    enterocyte in, absorption and faecal transit out.
    d/dt(gut_lumen_apc) <- kbile_apc_i * ehc3_apc +
      fgut_apc * (ar * ps_dif_eff_ent_apc + 0) * intestine_ent_apc -
      (ka_apc_i + kfeces_apc_i) * gut_lumen_apc

    #    Enterocyte (concentration) and mucosal blood (concentration).
    d/dt(intestine_ent_apc) <- (ka_apc_i * gut_lumen_apc +
      fb_apc * ps_dif_eff_ent_apc * intestine_muc_apc -
      fgut_apc * ((ar + 1) * ps_dif_eff_ent_apc + 0 + 0) *
      intestine_ent_apc +
      fgut * cl_apc_ent_i * intestine_ent) / v_ent
    d/dt(intestine_muc_apc) <- (q_mucosa * central_apc +
      fgut_apc * ps_dif_eff_ent_apc * intestine_ent_apc - q_mucosa * intestine_muc_apc -
      fb_apc * ps_dif_eff_ent_apc * intestine_muc_apc) / v_muc

    #    Excretion sinks (amounts).
    d/dt(a_feces_apc) <- kfeces_apc_i * gut_lumen_apc
    d/dt(a_urine_apc) <- cl_renal_apc_i * central_apc

    # ---------------- Exposure integrators ----------------
    #    The paper drives neutropenia off the unbound SN-38 AUC in plasma
    #    and delayed diarrhoea off the unbound SN-38 AUC in the intestinal
    #    epithelium (Methods, Performing VCSs; Discussion, Perspective of
    #    VCS, where the ID 2 thresholds are 26.35 and 53.60 nM*h).
    d/dt(auc_u_sn38) <- fb_sn38 * central_sn38
    d/dt(auc_u_ent_sn38) <- fgut_sn38 * intestine_ent_sn38

    # ===============================================================
    # 7. Dosing. `central` is a CONCENTRATION state, so an amount dose in
    #    umol is divided by the central volume on the way in. The
    #    infusion is given with rate = -2 in the event table, which is
    #    what makes rxode2 honour dur() rather than silently delivering a
    #    bolus.
    # ===============================================================
    f(central) <- 1 / v_central_i
    dur(central) <- tinf

    # ===============================================================
    # 8. Observations. Each state already is a blood concentration in
    #    umol/L.
    # ===============================================================
    Cc <- central
    Cc_sn38 <- central_sn38
    Cc_sn38g <- central_sn38g
    Cc_npc <- central_npc
    Cc_apc <- central_apc
    Cc ~ prop(propSd)
    Cc_sn38 ~ prop(propSd_sn38)
    Cc_sn38g ~ prop(propSd_sn38g)
    Cc_npc ~ prop(propSd_npc)
    Cc_apc ~ prop(propSd_apc)
  })
}
