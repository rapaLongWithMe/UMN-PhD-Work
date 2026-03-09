# Packages
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

setwd("~/Desktop/Git")

# ---- input: your Araport11 GFF3 (gz or plain) ----
gff_file <- "Araport11_GFF3_genes_transposons.20250813.gff.gz"


# ---- Load GFF ----
gff <- read_tsv(
  gff_file,
  comment = "#",
  col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"),
  col_types = "ccciicccc"
)

# ---- Helpers ----
grab_attr <- function(x, key) {
  # returns the first match of key=value from the attributes field
  m <- str_match(x, paste0("(^|;)", key, "=([^;]+)"))
  ifelse(is.na(m[,3]), NA_character_, m[,3])
}

mode_symbol <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ---- Basic coverage you computed ----
n_gene <- sum(gff$type == "gene")
n_gene_with_symbol <- gff |>
  filter(type == "gene") |>
  summarise(n = sum(str_detect(attributes, "symbol="))) |>
  pull(n)

n_with_symbol_any <- sum(str_detect(gff$attributes, "symbol="))
cat("Total gene features:", n_gene, "\n")
cat("Gene features with symbol (gene rows):", n_gene_with_symbol, "\n")
cat("Any features (gene/mRNA/etc.) with symbol:", n_with_symbol_any, "\n")

# ---- Extract locus-level info from 'gene' rows ----
genes <- gff %>%
  filter(type == "gene") %>%
  transmute(
    GeneID = grab_attr(attributes, "ID"),                   
    Symbol_gene = grab_attr(attributes, "symbol"),
    Desc1 = grab_attr(attributes, "curator_summary"),
    Desc2 = grab_attr(attributes, "computational_description"),
    Desc3 = grab_attr(attributes, "full_name"),
    Desc4 = grab_attr(attributes, "Note"),
    Description = coalesce(Desc1, Desc2, Desc3, Desc4)) %>%
  distinct(GeneID, .keep_all = TRUE)

# ---- Extract symbols from mRNA and map back to locus (Parent = ATxGxxxxx) ----
mrna <- gff %>%
  filter(type == "mRNA") %>%
  transmute(
    Parent    = grab_attr(attributes, "Parent"),            # locus ID (ATxGxxxxx)
    Symbol_m  = grab_attr(attributes, "symbol")) %>%
  filter(!is.na(Parent)) %>%
  # we may have multiple transcripts with differing/missing symbols;
  # pick a single symbol per Parent using majority vote (mode)
  group_by(Parent) %>%
  summarise(Symbol_mrna = mode_symbol(Symbol_m), .groups = "drop")

# ---- Merge: prefer gene-level symbol, else mRNA-level symbol ----
out <- genes %>%
  left_join(mrna, by = c("GeneID" = "Parent")) %>%
  mutate(Symbol = coalesce(Symbol_gene, Symbol_mrna)) %>%
  select(ID = GeneID, Symbol, Description) %>%
  distinct()

# ---- Coverage after rescuing from mRNA ----
n_with_symbol_final <- sum(!is.na(out$Symbol) & out$Symbol != "")
cat("Genes with a symbol after merging mRNA symbols:", n_with_symbol_final, "\n")
cat("Proportion with symbol:", sprintf("%.1f%%", 100 * n_with_symbol_final / nrow(out)), "\n")

# ---- Quick look at remaining missing symbols ----
#missing_syms <- out %>% filter(is.na(Symbol) | Symbol == "")
#cat("Remaining loci without a symbol:", nrow(missing_syms), "\n")
#print(head(missing_syms, 10))

# ---- Write output ----
write_csv(out, "Araport11_geneID_symbol_description.csv")

