library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)
library(readr)


# 1. Data input & phenotype mapping
mutation_raw <- read.table(
  "/data/sdd/Lanyiheng_2023/Cosmic-TCGA-RAWDATA-analysis/Cosmic_GenomeScreensMutant_v100_GRCh38.tsv",
  header = TRUE, sep = "\t", na.strings = TRUE, fill = TRUE
)


phenotype_info <- read_tsv(
  "/data/sdd/Lanyiheng_2023/Cosmic-TCGA-RAWDATA-analysis/gain_of_m6A/m6A_gain_of_m6A_V100_V2/Cosmic_Classification_v100_GRCh38.tsv"
)

# 2. Focus on gain of A mutations
mutation_set <- mutation_raw %>%
  filter(grepl(">A", MUTATION_CDS))

hgvsg_type <- mutation_set %>%
  group_by(HGVSG) %>%
  summarise(annotation_types = paste(unique(MUTATION_DESCRIPTION), collapse = ";"))

mutation_set <- mutation_set %>%
  left_join(hgvsg_type, by = "HGVSG")

mutation_set <- mutation_set %>%
  left_join(
    phenotype_info %>% select(COSMIC_PHENOTYPE_ID, PRIMARY_SITE),
    by = "COSMIC_PHENOTYPE_ID"
  )

# 3. Genome reference
genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens


# 4. Define annotation function for group-specific workflow
annotate_mutations <- function(mut_df, type_pattern, out_prefix, genome) {
  # Filter group
  df <- mut_df %>%
    filter(grepl(type_pattern, MUTATION_DESCRIPTION, ignore.case = TRUE)) %>%
    distinct(COSMIC_SAMPLE_ID, HGVSG, .keep_all = TRUE) %>%       # QC
    add_count(HGVSG, name = "mutation_count") %>%                 # recurrence
    filter(mutation_count > 1)                                    
  
  # Remove NA for position/strand
  df <- df %>% filter(!is.na(GENOME_START), !is.na(CHROMOSOME), !is.na(STRAND))
  
  # Standardize chromosome format
  df <- df %>%
    mutate(
      CHROMOSOME = paste0("chr", gsub("^chr", "", as.character(CHROMOSOME))),
      GENOME_START = as.integer(GENOME_START),
      STRAND = as.character(STRAND)
    ) %>%
    filter(!CHROMOSOME %in% c("chrMT", "chrM", "MT"))
  
  strand_vec <- df$STRAND
  
  # Get flanking sequences
  gr_5bp <- GRanges(seqnames = df$CHROMOSOME,
                    ranges = IRanges(start = df$GENOME_START - 2,
                                     end = df$GENOME_START + 2))
  seqs_5bp <- getSeq(genome, gr_5bp)
  seqs_5bp <- ifelse(strand_vec == "-", as.character(reverseComplement(seqs_5bp)), as.character(seqs_5bp))
  
  gr_up <- GRanges(seqnames = df$CHROMOSOME,
                   ranges = IRanges(start = df$GENOME_START - 2,
                                    end = df$GENOME_START - 1))
  seqs_up <- getSeq(genome, gr_up)
  
  gr_down <- GRanges(seqnames = df$CHROMOSOME,
                     ranges = IRanges(start = df$GENOME_START + 1,
                                      end = df$GENOME_START + 2))
  seqs_down <- getSeq(genome, gr_down)
  
  seqs_up_final <- ifelse(strand_vec == "-", as.character(reverseComplement(seqs_down)), as.character(seqs_up))
  seqs_down_final <- ifelse(strand_vec == "-", as.character(reverseComplement(seqs_up)), as.character(seqs_down))
  
  df$`5bp_flanking` <- seqs_5bp
  df$`2bp-_flanking` <- seqs_up_final
  df$`2bp+_flanking` <- seqs_down_final
  
  df <- df %>% filter(substr(`2bp+_flanking`, 1, 1) == "C")
  
  # DRACH motif annotation
  df <- df %>%
    mutate(
      DRACH_motif =
        substr(`5bp_flanking`, 1, 1) %in% c("A", "G", "T") &
        substr(`5bp_flanking`, 2, 2) %in% c("A", "G") &
        substr(`5bp_flanking`, 4, 4) == "C" &
        substr(`5bp_flanking`, 5, 5) %in% c("A", "C", "T")
    )
  
  # Output
  write.table(df, paste0(out_prefix, ".txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
  
  return(df)
}


# 5. Process SYN/MIS/UTR separately with annotation & QC
syn_df <- annotate_mutations(mutation_set, "synonymous", "GAIN-OF-M6A-syn", genome)
mis_df <- annotate_mutations(mutation_set, "missense", "GAIN-OF-M6A-mis", genome)
utr_df <- annotate_mutations(mutation_set, "UTR", "GAIN-OF-M6A-UTR", genome)
total_df <- bind_rows(syn_df, mis_df, utr_df)
write.table(
  total_df,
  "GAIN-OF-M6A-total.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
# 6. Downstream simple counting and statistical analyses were performed in Excel. 
#When a given HGVSG mutation has multiple functional annotations, missense mutations were prioritized for annotation, followed by UTR mutations, and synonymous mutations were assigned the lowest priority.




























