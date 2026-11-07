# Mock GO enrichment table: 10 parent terms, each with 4-9 child terms,
# children grouped into named subcategories of at most 3 terms.
# Columns mimic a clusterProfiler-style result plus parent / subcategory.
#   p.adjust   -> will map to colour (as -log10)
#   Count      -> will map to point size
# GO IDs are random placeholders, not real annotations.

make_mock_go <- function(seed = 42) {
  set.seed(seed)

  # parent -> subcategory -> child terms
  tree <- list(
    "immune system process" = list(
      "Lymphocytes" = c("T cell activation", "B cell proliferation",
                                "antigen presentation"),
      "Innate immunity"     = c("innate immune response",
                                "complement activation"),
      "Inflammation"        = c("cytokine production", "leukocyte migration")),
    "cell cycle" = list(
      "Mitosis"             = c("mitotic nuclear division",
                                "chromosome segregation",
                                "spindle assembly checkpoint"),
      "Transitions"   = c("G1/S transition", "G2/M transition",
                                "DNA replication initiation"),
      "Cell division"       = c("cytokinesis", "centrosome duplication")),
    "metabolic process" = list(
      "Lipids"    = c("fatty acid oxidation",
                                "cholesterol biosynthesis"),
      "Energy"   = c("glycolytic process", "TCA cycle"),
      "Building blocks" = c("amino acid catabolism",
                                      "nucleotide biosynthesis")),
    "signal transduction" = list(
      "Kinases"     = c("MAPK cascade", "JAK-STAT cascade",
                                "Ras protein signaling"),
      "Growth signaling" = c("Wnt signaling pathway",
                                    "Notch signaling pathway",
                                    "Hippo signaling"),
      "Receptors"  = c("GPCR signaling", "TGF-beta signaling",
                                "calcium-mediated signaling")),
    "cell death" = list(
      "Apoptosis"           = c("intrinsic apoptotic signaling",
                                "extrinsic apoptotic signaling"),
      "Non-apoptotic" = c("necroptotic process", "autophagic cell death",
                                "pyroptosis")),
    "development" = list(
      "Cardiovascular"      = c("heart development", "angiogenesis"),
      "Organogenesis"       = c("neuron differentiation", "kidney development",
                                "skeletal system development"),
      "Patterning"          = c("limb morphogenesis", "somite development",
                                "gastrulation")),
    "transport" = list(
      "Membrane"  = c("ion transmembrane transport",
                                "lipid transport"),
      "Trafficking" = c("vesicle-mediated transport",
                                      "protein import into nucleus")),
    "response to stimulus" = list(
      "Physical stress"     = c("response to heat", "response to osmotic stress",
                                "response to UV"),
      "Redox & nutrients" = c("response to oxidative stress",
                                     "response to hypoxia",
                                     "response to starvation"),
      "Damage response" = c("response to DNA damage",
                                  "response to wounding",
                                  "response to cytokine")),
    "gene expression" = list(
      "Transcription"       = c("transcription elongation",
                                "chromatin remodeling", "histone methylation"),
      "RNA & translation" = c("mRNA splicing", "rRNA processing",
                                         "translation initiation")),
    "cytoskeleton organization" = list(
      "Filaments" = c("actin filament polymerization",
                                  "microtubule nucleation",
                                  "intermediate filament organization"),
      "Adhesion & cilia"    = c("focal adhesion assembly", "cilium assembly"))
  )

  parent_ids <- sprintf("GO:%07d", sample(1e5:9e5, length(tree)))

  do.call(rbind, lapply(seq_along(tree), function(i) {
    terms <- unlist(tree[[i]], use.names = FALSE)
    subs  <- rep(names(tree[[i]]), lengths(tree[[i]]))
    n <- length(terms)
    base <- runif(1, 2, 12)  # parent-level strength: children loosely correlated
    data.frame(
      parent_id   = parent_ids[i],
      parent_term = names(tree)[i],
      subcategory = subs,
      ID          = sprintf("GO:%07d", sample(1e5:9e5, n)),
      Description = terms,
      GeneRatio   = NA_character_,
      Count       = sample(5:80, n),
      p.adjust    = 10^-pmax(1.5, rnorm(n, base, 2.5)),
      stringsAsFactors = FALSE
    )
  })) -> df

  df$GeneRatio <- sprintf("%d/%d", df$Count, 600)
  rownames(df) <- NULL
  df
}

if (sys.nframe() == 0) {
  go <- make_mock_go()
  print(table(go$parent_term))
  print(head(go, 12))
  write.csv(go, "mock_go.csv", row.names = FALSE)
}
