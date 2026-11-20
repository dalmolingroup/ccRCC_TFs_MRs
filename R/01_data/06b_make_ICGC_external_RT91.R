### Arquivos da validação externa no ICGC-RECA: 91 tumores primários (RT) ###
# Revisão CSBJ-D-26-01039 - Revisor 2, Comentários 1 e 3
#
# Gera as entradas do ml_pipeline_revision.py (external_validation) SÓ com tumores:
#   0_data/expr_MR_ICGC_RT91.csv     amostras (91 RT) x genes, mesmos IDs/colunas do expr_MR_KIRC_full.csv
#   0_data/metadata_ICGC_RT91.csv    amostras x (tissue_type, submitted_donor_id, pathologic_m, ...)
#
# Normalização: voom nas contagens dos 91 tumores (mesma do 21_survival_ICGC.R). O pipeline em Python
# ainda padroniza cada gene (z-score) dentro de cada coorte antes de aplicar o modelo do TCGA.
# Os RA ficam de fora: herdariam o estágio M do doador e distorceriam o z-score.

suppressPackageStartupMessages(library(dplyr))

ICGC_RDATA <- "~/Kidney/0_data/ICGC_DEA_Data.RData"     # count_final_ICGC, donor_data_final
TCGA_EXPR  <- "~/Kidney/0_data/expr_MR_KIRC_full.csv"    # só para copiar os IDs das colunas
OUT_EXPR   <- "~/Kidney/0_data/expr_MR_ICGC_RT91.csv"
OUT_META   <- "~/Kidney/0_data/metadata_ICGC_RT91.csv"
N_EXPECTED <- 91

pe <- path.expand
strip_v <- function(x) sub("\\..*", "", x)

## ---------------------------------------------------------------- tumores RT, 1 por doador
load(pe(ICGC_RDATA))
don <- as.data.frame(donor_data_final)
don$sample <- paste0(don$submitted_donor_id, don$tissue_type)
don <- don[don$tissue_type == "RT" & don$sample %in% colnames(count_final_ICGC), ]
don <- don[!duplicated(don$submitted_donor_id), ]
cat(sprintf("Tumores RT: %d (esperado %d)\n", nrow(don), N_EXPECTED))
stopifnot(nrow(don) == N_EXPECTED)

## ---------------------------------------------------------------- expressão (voom nos 91 RT)
cnt <- count_final_ICGC
rownames(cnt) <- strip_v(cnt$gene_id)
cnt <- as.matrix(cnt[, don$sample])
expr <- GDCRNATools::gdcVoomNormalization(cnt, filter = FALSE)      # genes x amostras

# mesmas colunas (IDs) do arquivo do TCGA usado no treino
tcga_cols <- colnames(read.csv(pe(TCGA_EXPR), row.names = 1, nrows = 1, check.names = FALSE))
keep <- strip_v(tcga_cols) %in% rownames(expr)
cat(sprintf("Genes do TCGA encontrados no ICGC: %d de %d\n", sum(keep), length(tcga_cols)))
if (any(!keep)) cat("Ausentes:", paste(head(tcga_cols[!keep], 20), collapse = ", "), "\n")

X <- t(expr[strip_v(tcga_cols[keep]), , drop = FALSE])               # amostras x genes
colnames(X) <- tcga_cols[keep]                                       # IDs idênticos aos do TCGA
write.csv(X, pe(OUT_EXPR))

## ---------------------------------------------------------------- metadados
meta <- don |> select(sample, submitted_donor_id, tissue_type, pathologic_m, pathologic_t,
                      pathologic_n, any_of(c("age", "gender", "obs.time", "status")))
rownames(meta) <- meta$sample; meta$sample <- NULL
write.csv(meta, pe(OUT_META))

cat("\nEstágio M nos 91 tumores:\n"); print(table(meta$pathologic_m, useNA = "ifany"))
cat("\nSalvos:\n ", OUT_EXPR, "\n ", OUT_META, "\n")
