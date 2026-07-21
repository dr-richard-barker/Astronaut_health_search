#!/usr/bin/env Rscript
# ==============================================================================
# 05_nutraceutical_flags.R
# astronaut-opposite-forcing pipeline
#
# Flags nutraceutical / food-derived compounds among the top LINCS reversal hits
# (T7_top_opposite_forcing.csv) as preparation for a deferred vegan-recipe
# optimization step.
#
# Three independent evidence streams are combined:
#   (1) Curated nutraceutical lexicon (polyphenols, flavonoids, vitamins,
#       phytochemicals, endogenous metabolites with dietary sources).
#   (2) LINCS L1000 single-drug-perturbation GMT names that match the lexicon
#       (independent transcriptomic evidence the compound has a measurable
#        perturbation signature).
#   (3) Broad Drug Repurposing Hub context (clinical_phase == "Launched" +
#       nutraceutical-class MoA strings: "antioxidant", "anti-inflammatory",
#       "flavonoid", "vitamin", "phytoestrogen", etc.).
#
# Outputs:
#   results/tables/T12_nutraceutical_flags.csv   -- per-hit nutraceutical flags
#   results/tables/T12b_nutraceutical_summary.csv -- aggregate counts + classes
#   results/figures/F10_nutraceutical_flags.svg / .png
#
# Richard Barker -- 2026
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

repo <- "/mnt/shared-workspace/astronaut-opposite-forcing"
fig_dir <- file.path(repo, "results", "figures")
tab_dir <- file.path(repo, "results", "tables")
proc_dir <- file.path(repo, "data", "processed")
gmt_path <- "/mnt/datalake/LINCS1000/RNAseq_transcriptomics_genesets/single_drug_perturbations-v1.0.gmt"
broad_path <- file.path(proc_dir, "broad_phase_moa_target.csv")
t7_path <- file.path(tab_dir, "T7_top_opposite_forcing.csv")
t6_path <- file.path(tab_dir, "T6_lincs_reversal_scores.csv")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Curated nutraceutical lexicon (explicit per-row construction)
# ------------------------------------------------------------------------------
# Each entry: canonical name + class + primary dietary source(s).
# Sources: USDA phytochemical databases, Linus Pauling Institute Micronutrient
# Information Center, MSigDB chemical perturbations, published reviews of
# dietary polyphenols and spaceflight-relevant antioxidants.
# Built as an explicit list of (name, class, source) triples to avoid
# vector-length mismatch from rep() alignment.
.lex <- list(
  # Flavonoids
  list("apigenin","flavonoid","parsley, celery, chamomile, citrus, onions, apples, berries, green tea, soy"),
  list("luteolin","flavonoid","parsley, celery, chamomile, citrus, onions, apples, berries, green tea, soy"),
  list("quercetin","flavonoid","parsley, celery, chamomile, citrus, onions, apples, berries, green tea, soy"),
  list("kaempferol","flavonoid","parsley, celery, chamomile, citrus, onions, apples, berries, green tea, soy"),
  list("myricetin","flavonoid","parsley, celery, chamomile, citrus, onions, apples, berries, green tea, soy"),
  list("fisetin","flavonoid","parsley, celery, chamomile, citrus, onions, apples, berries, green tea, soy"),
  list("catechin","flavonoid","green tea, cocoa, apples, berries"),
  list("epicatechin","flavonoid","green tea, cocoa, apples, berries"),
  list("epigallocatechin gallate","flavonoid","green tea, cocoa, apples, berries"),
  list("egcg","flavonoid","green tea, cocoa, apples, berries"),
  list("epigallocatechin","flavonoid","green tea, cocoa, apples, berries"),
  list("theaflavin","flavonoid","black tea, green tea"),
  list("theasinensin a","flavonoid","green tea, black tea"),
  list("naringenin","flavonoid","citrus (grapefruit, oranges), tomatoes"),
  list("hesperetin","flavonoid","citrus (oranges, lemons)"),
  list("hesperidin","flavonoid","citrus (oranges, lemons)"),
  list("daidzein","flavonoid","soy, legumes"),
  list("genistein","flavonoid","soy, legumes"),
  list("biochanin a","flavonoid","chickpeas, red clover, soy"),
  list("formononetin","flavonoid","red clover, chickpeas"),
  list("coumestrol","flavonoid","alfalfa, clover, soy sprouts"),
  list("chrysin","flavonoid","honey, propolis, passionflower"),
  list("galangin","flavonoid","honey, propolis, galangal"),
  list("baicalein","flavonoid","skullcap (Scutellaria baicalensis)"),
  list("wogonin","flavonoid","skullcap (Scutellaria baicalensis)"),
  list("isoginkgetin","flavonoid","ginkgo biloba"),
  list("licochalcone a","flavonoid","licorice root (Glycyrrhiza)"),
  list("licochalcone","flavonoid","licorice root (Glycyrrhiza)"),
  list("isoliquiritigenin","flavonoid","licorice root (Glycyrrhiza)"),
  list("glycyrrhizin","flavonoid","licorice root (Glycyrrhiza)"),
  list("glycyrrhetinic acid","flavonoid","licorice root (Glycyrrhiza)"),
  # Anthocyanins
  list("cyanidin","anthocyanin","berries, red cabbage, grapes, eggplant, plums"),
  list("delphinidin","anthocyanin","berries, red cabbage, grapes, eggplant, plums"),
  list("malvidin","anthocyanin","berries, red cabbage, grapes, eggplant, plums"),
  list("pelargonidin","anthocyanin","berries, red cabbage, grapes, eggplant, plums"),
  list("peonidin","anthocyanin","berries, red cabbage, grapes, eggplant, plums"),
  list("petunidin","anthocyanin","berries, red cabbage, grapes, eggplant, plums"),
  list("malvidin-3-glucoside","anthocyanin","berries, red cabbage, grapes, eggplant, plums"),
  # Stilbenes
  list("resveratrol","stilbene","grapes, blueberries, peanuts, red wine"),
  list("pterostilbene","stilbene","blueberries, grapes"),
  list("piceatannol","stilbene","grapes, passion fruit"),
  list("oxyresveratrol","stilbene","mulberry, artocarpus"),
  # Curcuminoids
  list("curcumin","curcuminoid","turmeric (Curcuma longa)"),
  list("demethoxycurcumin","curcuminoid","turmeric (Curcuma longa)"),
  list("bisdemethoxycurcumin","curcuminoid","turmeric (Curcuma longa)"),
  list("tetrahydrocurcumin","curcuminoid","turmeric (Curcuma longa)"),
  # Phenolic acids / tannins
  list("gallic acid","phenolic_acid_tannin","berries, pomegranate, coffee, rosemary, green tea, whole grains"),
  list("ellagic acid","phenolic_acid_tannin","berries, pomegranate, coffee, rosemary, green tea, whole grains"),
  list("ellagitannin","phenolic_acid_tannin","berries, pomegranate, coffee, rosemary, green tea, whole grains"),
  list("urolithin a","phenolic_acid_tannin","pomegranate, walnuts (gut-microbiota metabolite)"),
  list("urolithin b","phenolic_acid_tannin","pomegranate, walnuts (gut-microbiota metabolite)"),
  list("rosmarinic acid","phenolic_acid_tannin","rosemary, oregano, lemon balm, mint"),
  list("carnosic acid","phenolic_acid_tannin","rosemary, sage"),
  list("carnosol","phenolic_acid_tannin","rosemary, sage"),
  list("chlorogenic acid","phenolic_acid_tannin","coffee, sunflower seeds, artichoke"),
  list("caffeic acid","phenolic_acid_tannin","coffee, fruits, vegetables"),
  list("ferulic acid","phenolic_acid_tannin","whole grains, rice bran, oats"),
  list("coumaric acid","phenolic_acid_tannin","fruits, vegetables, whole grains"),
  list("sinapic acid","phenolic_acid_tannin","mustard seed, rapeseed, bran"),
  list("protocatechuic acid","phenolic_acid_tannin","berries, olives, tea"),
  list("vanillic acid","phenolic_acid_tannin","vanilla, berries, tea"),
  list("syringic acid","phenolic_acid_tannin","berries, wine, spices"),
  # Glucosinolates / isothiocyanates
  list("sulforaphane","glucosinolate_isothiocyanate","cruciferous vegetables (broccoli, Brussels sprouts, watercress)"),
  list("sulforaphane n-acetylcysteine","glucosinolate_isothiocyanate","cruciferous vegetables (broccoli, Brussels sprouts, watercress)"),
  list("indole-3-carbinol","glucosinolate_isothiocyanate","cruciferous vegetables (broccoli, Brussels sprouts, watercress)"),
  list("diindolylmethane","glucosinolate_isothiocyanate","cruciferous vegetables (broccoli, Brussels sprouts, watercress)"),
  list("phenethyl isothiocyanate","glucosinolate_isothiocyanate","watercress, radish"),
  list("allyl isothiocyanate","glucosinolate_isothiocyanate","mustard, horseradish, wasabi"),
  # Carotenoids
  list("lycopene","carotenoid","tomato, watermelon, guava, pink grapefruit"),
  list("beta-carotene","carotenoid","carrot, sweet potato, leafy greens, pumpkin"),
  list("alpha-carotene","carotenoid","carrot, sweet potato, leafy greens"),
  list("lutein","carotenoid","leafy greens (kale, spinach), egg yolk, marigold"),
  list("zeaxanthin","carotenoid","leafy greens, corn, egg yolk, orange pepper"),
  list("astaxanthin","carotenoid","salmon, shrimp, krill, algae"),
  list("canthaxanthin","carotenoid","mushrooms, algae, salmon"),
  list("fucoxanthin","carotenoid","brown seaweed, algae"),
  list("cryptoxanthin","carotenoid","papaya, mango, citrus, pumpkin"),
  # Vitamins fat-soluble
  list("vitamin d","vitamin_fat_soluble","fish, eggs, fortified dairy, leafy greens, vegetable oils, nuts"),
  list("vitamin d3","vitamin_fat_soluble","fish, eggs, fortified dairy, leafy greens, vegetable oils, nuts"),
  list("cholecalciferol","vitamin_fat_soluble","fish, eggs, fortified dairy, UV-exposed mushrooms"),
  list("ergocalciferol","vitamin_fat_soluble","UV-exposed mushrooms, fortified foods"),
  list("calcitriol","vitamin_fat_soluble","endogenous synthesis from vitamin D"),
  list("alfacalcidol","vitamin_fat_soluble","synthetic vitamin D analog"),
  list("vitamin a","vitamin_fat_soluble","liver, dairy, eggs, leafy greens, orange vegetables"),
  list("retinol","vitamin_fat_soluble","liver, dairy, eggs, fortified foods"),
  list("retinoic acid","vitamin_fat_soluble","endogenous metabolite of vitamin A"),
  list("all-trans retinoic acid","vitamin_fat_soluble","endogenous metabolite of vitamin A"),
  list("13-cis retinoic acid","vitamin_fat_soluble","endogenous metabolite of vitamin A"),
  list("vitamin e","vitamin_fat_soluble","vegetable oils, nuts, seeds, leafy greens"),
  list("alpha-tocopherol","vitamin_fat_soluble","vegetable oils, nuts, seeds, leafy greens"),
  list("gamma-tocopherol","vitamin_fat_soluble","vegetable oils, nuts, seeds"),
  list("tocotrienol","vitamin_fat_soluble","palm oil, rice bran, barley"),
  list("vitamin k","vitamin_fat_soluble","leafy greens, broccoli, Brussels sprouts"),
  list("phylloquinone","vitamin_fat_soluble","leafy greens, broccoli, Brussels sprouts"),
  list("menaquinone","vitamin_fat_soluble","fermented foods, animal products, natto"),
  # Vitamins water-soluble
  list("vitamin c","vitamin_water_soluble","citrus, peppers, berries, leafy greens"),
  list("ascorbic acid","vitamin_water_soluble","citrus, peppers, berries, leafy greens"),
  list("ascorbate","vitamin_water_soluble","citrus, peppers, berries, leafy greens"),
  list("thiamine","vitamin_water_soluble","whole grains, legumes, pork, yeast"),
  list("riboflavin","vitamin_water_soluble","dairy, eggs, leafy greens, whole grains"),
  list("niacin","vitamin_water_soluble","meat, fish, poultry, peanuts, whole grains"),
  list("nicotinamide","vitamin_water_soluble","meat, fish, poultry, peanuts, whole grains"),
  list("niacinamide","vitamin_water_soluble","meat, fish, poultry, peanuts, whole grains"),
  list("pantothenic acid","vitamin_water_soluble","meat, eggs, whole grains, legumes, mushrooms"),
  list("pyridoxine","vitamin_water_soluble","poultry, fish, potatoes, bananas, chickpeas"),
  list("biotin","vitamin_water_soluble","eggs, liver, nuts, seeds, sweet potato"),
  list("folate","vitamin_water_soluble","leafy greens, legumes, citrus, fortified grains"),
  list("folic acid","vitamin_water_soluble","fortified grains, leafy greens, legumes"),
  list("cobalamin","vitamin_water_soluble","animal products, fortified foods, nutritional yeast"),
  list("vitamin b12","vitamin_water_soluble","animal products, fortified foods, nutritional yeast"),
  # Minerals
  list("magnesium","mineral","nuts, seeds, legumes, leafy greens, whole grains"),
  list("zinc","mineral","pumpkin seeds, legumes, nuts, whole grains, meat"),
  list("selenium","mineral","brazil nuts, seafood, meat, whole grains"),
  list("iodine","mineral","seaweed, iodized salt, seafood, dairy"),
  list("potassium","mineral","bananas, potatoes, leafy greens, beans, avocado"),
  # Endogenous / diet-derived metabolites
  list("creatine","endogenous_metabolite","red meat, fish, endogenous synthesis"),
  list("l-carnitine","endogenous_metabolite","red meat, fish, dairy, endogenous synthesis"),
  list("carnitine","endogenous_metabolite","red meat, fish, dairy, endogenous synthesis"),
  list("taurine","endogenous_metabolite","seafood, meat, endogenous synthesis"),
  list("betaine","endogenous_metabolite","beets, spinach, whole grains, quinoa"),
  list("coenzyme q10","endogenous_metabolite","organ meats, fatty fish, whole grains, endogenous synthesis"),
  list("coq10","endogenous_metabolite","organ meats, fatty fish, whole grains, endogenous synthesis"),
  list("ubiquinone","endogenous_metabolite","organ meats, fatty fish, whole grains, endogenous synthesis"),
  list("alpha-lipoic acid","endogenous_metabolite","organ meats, spinach, broccoli, endogenous synthesis"),
  list("lipoic acid","endogenous_metabolite","organ meats, spinach, broccoli, endogenous synthesis"),
  list("melatonin","endogenous_metabolite","tart cherries, walnuts, oats, endogenous synthesis"),
  list("n-acetylcysteine","endogenous_metabolite","supplement (cysteine precursor)"),
  list("glutathione","endogenous_metabolite","avocado, asparagus, spinach, endogenous synthesis"),
  list("s-adenosylmethionine","endogenous_metabolite","endogenous synthesis (methionine metabolism)"),
  list("nicotinamide riboside","endogenous_metabolite","dairy milk, yeast, endogenous NAD+ precursor"),
  list("nicotinamide mononucleotide","endogenous_metabolite","endogenous NAD+ precursor, vegetables (trace)"),
  # Organosulfur / terpenoid phytochemicals
  list("allicin","organosulfur_terpenoid","garlic, onion, leek"),
  list("diallyl disulfide","organosulfur_terpenoid","garlic"),
  list("diallyl trisulfide","organosulfur_terpenoid","garlic"),
  list("ajoene","organosulfur_terpenoid","garlic (aged extract)"),
  list("gingerol","organosulfur_terpenoid","ginger (Zingiber officinale)"),
  list("shogaol","organosulfur_terpenoid","ginger (dried/heated)"),
  list("paradol","organosulfur_terpenoid","grains of paradise, ginger"),
  list("piperine","organosulfur_terpenoid","black pepper (Piper nigrum)"),
  list("capsaicin","organosulfur_terpenoid","chili peppers (Capsicum)"),
  list("cannabidiol","organosulfur_terpenoid","hemp, cannabis (non-psychoactive)"),
  list("cbd","organosulfur_terpenoid","hemp, cannabis (non-psychoactive)"),
  list("caryophyllene","organosulfur_terpenoid","black pepper, cloves, cinnamon, hops"),
  list("humulene","organosulfur_terpenoid","hops, sage, ginseng"),
  list("menthol","organosulfur_terpenoid","peppermint, mint"),
  list("carvacrol","organosulfur_terpenoid","oregano, thyme"),
  list("thymol","organosulfur_terpenoid","thyme, oregano"),
  list("eugenol","organosulfur_terpenoid","cloves, cinnamon, basil"),
  list("cinnamaldehyde","organosulfur_terpenoid","cinnamon (Cinnamomum)"),
  # Lignans
  list("secoisolariciresinol","lignan","flaxseed, sesame, whole grains, legumes"),
  list("matairesinol","lignan","flaxseed, sesame, whole grains, sunflower seeds"),
  list("enterolactone","lignan","flaxseed, sesame (gut-microbiota metabolite)"),
  list("enterodiol","lignan","flaxseed, sesame (gut-microbiota metabolite)"),
  # Saponins / withanolides
  list("withaferin a","saponin_withanolide","ashwagandha (Withania somnifera)"),
  list("withanolide","saponin_withanolide","ashwagandha (Withania somnifera)"),
  list("guggulsterone","saponin_withanolide","guggul (Commiphora mukul)"),
  list("ginsenoside","saponin_withanolide","ginseng (Panax)"),
  list("silymarin","saponin_withanolide","milk thistle (Silybum marianum)"),
  list("silibinin","saponin_withanolide","milk thistle (Silybum marianum)"),
  # Omega-3 / lipids
  list("epa","omega3_lipid","fatty fish, fish oil, algae, krill"),
  list("dha","omega3_lipid","fatty fish, fish oil, algae, krill"),
  list("alpha-linolenic acid","omega3_lipid","flaxseed, chia, walnuts, hemp seed"),
  list("linoleic acid","omega3_lipid","safflower, sunflower, corn, soybean oils"),
  # Alkaloids / misc phytochemicals
  list("berberine","alkaloid_misc","barberry, goldenseal, Oregon grape, Coptis"),
  list("columbamine","alkaloid_misc","Coptis, barberry"),
  list("palmatine","alkaloid_misc","Coptis, goldenseal, barberry"),
  list("berbamine","alkaloid_misc","barberry, magnolia"),
  list("honokiol","alkaloid_misc","magnolia (Magnolia officinalis)"),
  list("magnolol","alkaloid_misc","magnolia (Magnolia officinalis)"),
  list("oxyberberine","alkaloid_misc","barberry, Coptis"),
  list("indole-3-propionic acid","alkaloid_misc","gut-microbiota metabolite of tryptophan"),
  list("indole-3-acetic acid","alkaloid_misc","gut-microbiota metabolite, cruciferous vegetables"),
  # Short-chain fatty acids
  list("butyrate","short_chain_fatty_acid","fermented dietary fiber (gut microbiota)"),
  list("sodium butyrate","short_chain_fatty_acid","fermented dietary fiber (gut microbiota)"),
  list("propionate","short_chain_fatty_acid","fermented dietary fiber (gut microbiota)"),
  list("acetate","short_chain_fatty_acid","fermented dietary fiber (gut microbiota)")
)
nutra_lexicon <- rbindlist(lapply(.lex, function(x)
  data.table(canonical = x[[1]], class = x[[2]], source_food = x[[3]])))
rm(.lex)
# Deduplicate (some canonical names appear in multiple class groupings above)
nutra_lexicon <- unique(nutra_lexicon, by = "canonical")

cat(sprintf("[nutra] curated lexicon: %d unique compounds across %d classes\n",
            nrow(nutra_lexicon), length(unique(nutra_lexicon$class))))

# ------------------------------------------------------------------------------
# 2. LINCS GMT drug names (independent transcriptomic evidence)
# ------------------------------------------------------------------------------
lincs_drugs <- data.table()
if (file.exists(gmt_path)) {
  gmt_lines <- readLines(gmt_path, n = 1e5)
  # GMT format: each line is "NAME\tDESCRIPTION\tgene1\tgene2\t..."
  # The drug name is the first tab-separated field. Names ending in -up or -dn
  # are the two directions of a single drug perturbation; strip the suffix.
  first_fields <- sub("\t.*$", "", gmt_lines)
  first_fields <- first_fields[nzchar(first_fields)]
  drug_names <- unique(sub("-(up|dn)$", "", first_fields))
  lincs_drugs <- data.table(lincs_name = drug_names,
                            lincs_name_lower = tolower(drug_names))
  cat(sprintf("[nutra] LINCS GMT: %d unique drug perturbation signatures\n",
              nrow(lincs_drugs)))
} else {
  cat("[nutra] WARNING: LINCS GMT not found at", gmt_path, "\n")
}

# ------------------------------------------------------------------------------
# 3. Broad Drug Repurposing Hub context
# ------------------------------------------------------------------------------
broad <- data.table()
if (file.exists(broad_path)) {
  broad <- fread(broad_path)
  cat(sprintf("[nutra] Broad repurposing hub: %d drugs\n", nrow(broad)))
} else {
  cat("[nutra] WARNING: Broad CSV not found at", broad_path, "\n")
}

# ------------------------------------------------------------------------------
# 4. Load top reversal hits (T7) and full reversal scores (T6)
# ------------------------------------------------------------------------------
t7 <- fread(t7_path)
t6 <- fread(t6_path)
cat(sprintf("[nutra] T7 top reversal hits: %d rows; T6 full scores: %d rows\n",
            nrow(t7), nrow(t6)))

# Normalize the drug-name field for matching. T7's `drug_name` is the
# extract_drug_name() output (often a GSE ID fallback). The full signature
# string carries the real compound name, so we mine both.
# IMPORTANT: lowercase FIRST, then strip non-alphanumeric -- otherwise
# `gsub("[^a-z0-9]", ...)` would delete uppercase letters before tolower().
normalize_name <- function(x) {
  x <- tolower(x)
  gsub("[^a-z0-9]", " ", x)
}
sig_lower <- normalize_name(t7$signature)
drug_lower <- normalize_name(t7$drug_name)

# ------------------------------------------------------------------------------
# 5. Match each T7 hit against the nutraceutical lexicon
# ------------------------------------------------------------------------------
# For every lexicon entry, search for it as a whole-word token in either the
# extracted drug_name or the full signature string. Record which evidence
# streams fired.
flag_rows <- list()
for (i in seq_len(nrow(nutra_lexicon))) {
  nm <- nutra_lexicon$canonical[i]
  cls <- nutra_lexicon$class[i]
  src <- nutra_lexicon$source_food[i]
  # word-boundary regex (handles "apigenin" inside "apigenin and luteolin")
  pat <- sprintf("(^|[^a-z0-9])%s([^a-z0-9]|$)", gsub(" ", "[ _-]?", nm))
  hit_sig <- grepl(pat, sig_lower, perl = TRUE)
  hit_drug <- grepl(pat, drug_lower, perl = TRUE)
  hit_lincs <- nrow(lincs_drugs) > 0 &&
               any(grepl(pat, lincs_drugs$lincs_name_lower, perl = TRUE))
  hit_broad <- nrow(broad) > 0 &&
               any(grepl(pat, normalize_name(broad$pert_iname), perl = TRUE))
  hits <- which(hit_sig | hit_drug)
  if (length(hits) == 0) next
  for (h in hits) {
    flag_rows[[length(flag_rows) + 1]] <- data.table(
      t7_row = h,
      signature = t7$signature[h],
      drug_name_extracted = t7$drug_name[h],
      tau = t7$tau[h],
      wcs = t7$wcs[h],
      nutraceutical = nutra_lexicon$canonical[i],
      class = cls,
      source_food = src,
      evidence_signature = hit_sig[h],
      evidence_drug_name = hit_drug[h],
      evidence_lincs_gmt = hit_lincs,
      evidence_broad = hit_broad,
      n_evidence_streams = hit_sig[h] + hit_drug[h] + hit_lincs + hit_broad
    )
  }
}

flags <- if (length(flag_rows)) rbindlist(flag_rows) else data.table()
cat(sprintf("[nutra] matched %d nutraceutical flag rows across %d T7 hits\n",
            nrow(flags), length(unique(flags$t7_row))))

# ------------------------------------------------------------------------------
# 6. Also scan the FULL T6 (all 420 scored drugs) for nutraceuticals that did
#    not make the top-50 cutoff -- useful for the deferred recipe step.
# ------------------------------------------------------------------------------
# T6 schema is (signature, wcs, tau) -- no extracted drug_name column, so we
# mine the full signature string only.
t6_sig_lower <- normalize_name(t6$signature)
t7_sigs <- t7$signature
t6_flag_rows <- list()
for (i in seq_len(nrow(nutra_lexicon))) {
  nm <- nutra_lexicon$canonical[i]
  cls <- nutra_lexicon$class[i]
  src <- nutra_lexicon$source_food[i]
  pat <- sprintf("(^|[^a-z0-9])%s([^a-z0-9]|$)", gsub(" ", "[ _-]?", nm))
  hit_sig <- grepl(pat, t6_sig_lower, perl = TRUE)
  hits <- which(hit_sig)
  if (length(hits) == 0) next
  for (h in hits) {
    t6_flag_rows[[length(t6_flag_rows) + 1]] <- data.table(
      signature = t6$signature[h],
      tau = t6$tau[h],
      wcs = t6$wcs[h],
      nutraceutical = nm,
      class = cls,
      source_food = src,
      in_top50 = t6$signature[h] %in% t7_sigs
    )
  }
}
t6_flags <- if (length(t6_flag_rows)) rbindlist(t6_flag_rows) else data.table()
# Deduplicate by signature + nutraceutical
t6_flags <- unique(t6_flags, by = c("signature", "nutraceutical"))
cat(sprintf("[nutra] full T6 scan: %d nutraceutical-bearing signatures (%d in top-50)\n",
            nrow(t6_flags), sum(t6_flags$in_top50)))

# ------------------------------------------------------------------------------
# 7. Write T12 (per-hit flags) and T12b (summary)
# ------------------------------------------------------------------------------
# Order T12 by tau (most negative = strongest reversal first), then by evidence
if (nrow(flags) > 0) {
  flags <- flags[order(tau, -n_evidence_streams)]
}
fwrite(flags, file.path(tab_dir, "T12_nutraceutical_flags.csv"))
cat(sprintf("[nutra] wrote %s (%d rows)\n",
            "T12_nutraceutical_flags.csv", nrow(flags)))

# Summary by class
if (nrow(flags) > 0) {
  summary_dt <- flags[, .(
    n_hits = .N,
    n_unique_compounds = uniqueN(nutraceutical),
    mean_tau = mean(tau, na.rm = TRUE),
    min_tau = min(tau, na.rm = TRUE),
    max_evidence_streams = max(n_evidence_streams)
  ), by = class][order(n_hits, decreasing = TRUE)]
} else {
  summary_dt <- data.table(class = character(), n_hits = integer(),
                           n_unique_compounds = integer(),
                           mean_tau = numeric(), min_tau = numeric(),
                           max_evidence_streams = integer())
}
fwrite(summary_dt, file.path(tab_dir, "T12b_nutraceutical_summary.csv"))
cat(sprintf("[nutra] wrote T12b_nutraceutical_summary.csv (%d classes)\n",
            nrow(summary_dt)))

# Also write the full-T6 nutraceutical scan for the deferred recipe step
t12c_path <- file.path(tab_dir, "T12c_nutraceutical_full_scan.csv")
if (file.exists(t12c_path) && nrow(t6_flags) == 0) file.remove(t12c_path)
fwrite(t6_flags, t12c_path)
cat(sprintf("[nutra] wrote T12c_nutraceutical_full_scan.csv (%d rows)\n",
            nrow(t6_flags)))

# ------------------------------------------------------------------------------
# 8. Figure F10 -- nutraceutical flags among top reversal hits
# ------------------------------------------------------------------------------
# Panel A: bar chart of nutraceutical class counts among top-50 hits.
# Panel B: dot plot of tau vs nutraceutical, colored by class, sized by evidence.
svg_path <- file.path(fig_dir, "F10_nutraceutical_flags.svg")
png_path <- file.path(fig_dir, "F10_nutraceutical_flags.png")

# Build a combined plot only if we have flags; otherwise emit an empty-panel
# placeholder so the figure file always exists for the manuscript.
if (nrow(flags) > 0) {
  # Panel A: class counts
  class_counts <- flags[, .N, by = class][order(N, decreasing = TRUE)]
  class_counts[, class := factor(class, levels = class_counts$class)]
  pA <- ggplot(class_counts, aes(x = N, y = class, fill = class)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = N), hjust = -0.2, size = 3.5) +
    scale_fill_manual(values = rep("#0279EE", nrow(class_counts))) +
    labs(title = "Nutraceutical classes among top-50 reversal hits",
         subtitle = sprintf("%d compounds flagged across %d classes",
                            nrow(flags), nrow(class_counts)),
         x = "Number of flagged hits", y = NULL) +
    theme_minimal(base_family = "Liberation Sans") +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    expand_limits(x = max(class_counts$N) * 1.2)

  # Panel B: tau per nutraceutical compound (most negative = strongest reversal)
  comp_dt <- unique(flags[, .(nutraceutical, class, tau, n_evidence_streams)],
                    by = c("nutraceutical", "class"))
  comp_dt <- comp_dt[order(tau)]
  comp_dt[, nutraceutical := factor(nutraceutical,
                                    levels = comp_dt$nutraceutical)]
  pB <- ggplot(comp_dt, aes(x = tau, y = nutraceutical,
                            color = class, size = n_evidence_streams)) +
    geom_point() +
    geom_vline(xintercept = -50, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(
      "flavonoid" = "#75A025", "anthocyanin" = "#FD9BED",
      "stilbene" = "#FF9400", "curcuminoid" = "#E9ED4C",
      "phenolic_acid_tannin" = "#0279EE", "glucosinolate_isothiocyanate" = "#75A025",
      "carotenoid" = "#FF9400", "vitamin_fat_soluble" = "#0279EE",
      "vitamin_water_soluble" = "#0279EE", "mineral" = "#000000",
      "endogenous_metabolite" = "#FD9BED", "organosulfur_terpenoid" = "#75A025",
      "lignan" = "#E9ED4C", "saponin_withanolide" = "#FF9400",
      "omega3_lipid" = "#0279EE", "alkaloid_misc" = "#FD9BED",
      "short_chain_fatty_acid" = "#000000"
    )) +
    scale_size_continuous(range = c(2, 6), name = "Evidence streams") +
    labs(title = "Reversal strength (tau) of flagged nutraceuticals",
         subtitle = "Dashed line = tau = -50 reversal threshold",
         x = "LINCS tau-analog (more negative = stronger reversal)",
         y = NULL, color = "Class") +
    theme_minimal(base_family = "Liberation Sans") +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          legend.position = "right")

  # Combine with cowplot-free patchwork via gridExtra if available, else print
  combined <- tryCatch({
    if (requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      pA + pB + plot_layout(ncol = 2, widths = c(1, 1.4))
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      library(gridExtra)
      gridExtra::grid.arrange(pA, pB, ncol = 2, widths = c(1, 1.4))
    } else {
      pB
    }
  }, error = function(e) pB)

  ggsave(svg_path, combined, width = 12, height = 7, dpi = 150)
  ggsave(png_path, combined, width = 12, height = 7, dpi = 150)
} else {
  # Placeholder
  p <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                           label = "No nutraceuticals flagged in top-50 reversal hits",
                           size = 6) +
    theme_void() + theme(plot.title = element_text(face = "bold"))
  ggsave(svg_path, p, width = 10, height = 6, dpi = 150)
  ggsave(png_path, p, width = 10, height = 6, dpi = 150)
}
cat(sprintf("[nutra] wrote F10_nutraceutical_flags.svg/png\n"))

# ------------------------------------------------------------------------------
# 9. Console summary
# ------------------------------------------------------------------------------
cat("\n========== NUTRACEUTICAL FLAG SUMMARY ==========\n")
cat(sprintf("Top-50 reversal hits scanned: %d\n", nrow(t7)))
cat(sprintf("Curated lexicon size: %d compounds, %d classes\n",
            nrow(nutra_lexicon), length(unique(nutra_lexicon$class))))
cat(sprintf("Flagged nutraceutical hits in top-50: %d rows, %d unique compounds\n",
            nrow(flags), uniqueN(flags$nutraceutical)))
if (nrow(flags) > 0) {
  cat("\nTop flagged nutraceuticals (by reversal strength):\n")
  print(head(unique(flags[, .(nutraceutical, class, tau, n_evidence_streams)],
                   by = c("nutraceutical", "class"))[order(tau)], 15))
}
cat(sprintf("\nFull T6 scan: %d nutraceutical-bearing signatures total\n",
            nrow(t6_flags)))
cat("================================================\n")

cat("\n[05_nutraceutical_flags.R] DONE\n")
