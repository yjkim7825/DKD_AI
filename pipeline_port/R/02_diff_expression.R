# 02_diff_expression.R ----------------------------------------------------
# 원본: differential expression analysis.R (limma lmFit+eBayes+topTable, 2-group)
# config 기반으로 정리. 입력은 '{샘플}_{그룹}' 라벨 발현행렬.

suppressMessages({ library(limma); library(dplyr); library(ggplot2) })
if (!exists("DATA_DIR")) source("config.R")

# input : 발현행렬(1열=유전자, 헤더='{샘플}_{그룹}')
# ref/alt: 대조/실험군 라벨. tag: 출력 접미사.
run_deg <- function(input, ref = CONTROL_LABEL, alt = CASE_LABEL,
                    tag = paste0(alt, "_vs_", ref),
                    outdir = file.path(RESULT_DIR, "deg"),
                    logfc = LOGFC_FILTER, adjp = ADJP_FILTER, plot = TRUE) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  rt <- read.table(input, header = TRUE, sep = "\t", check.names = FALSE)
  rownames(rt) <- rt[, 1]; exp <- as.matrix(rt[, -1])
  data <- matrix(as.numeric(exp), nrow = nrow(exp),
                 dimnames = list(rownames(exp), colnames(exp)))
  data <- avereps(data)

  Type <- sapply(strsplit(colnames(data), "_"), function(x) tail(x, 1))
  keep <- Type %in% c(ref, alt)
  data <- data[, keep]; Type <- Type[keep]
  data <- data[, order(Type)]; Type <- Type[order(Type)]
  print(table(Type))

  design <- model.matrix(~0 + factor(Type))
  colnames(design) <- levels(factor(Type))
  fit <- lmFit(data, design)
  cont <- makeContrasts(contrasts = paste0(alt, " - ", ref), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))

  allDiff <- topTable(fit2, adjust = "fdr", number = Inf)
  write.table(rbind(id = colnames(allDiff), allDiff),
              file.path(outdir, paste0("all_", tag, ".txt")),
              sep = "\t", quote = FALSE, col.names = FALSE)

  diffSig <- allDiff[abs(allDiff$logFC) > logfc & allDiff$adj.P.Val < adjp, ]
  write.table(rbind(id = colnames(diffSig), diffSig),
              file.path(outdir, paste0("diff_", tag, ".txt")),
              sep = "\t", quote = FALSE, col.names = FALSE)
  up   <- rownames(diffSig[diffSig$logFC > 0, ])
  down <- rownames(diffSig[diffSig$logFC < 0, ])
  writeLines(up,   file.path(outdir, paste0("up_genes_", tag, ".txt")))
  writeLines(down, file.path(outdir, paste0("down_genes_", tag, ".txt")))
  message("[DEG] ", nrow(diffSig), " DEG (up=", length(up), ", down=", length(down), ")")

  if (plot) {
    rt2 <- allDiff
    rt2$Sig <- ifelse(rt2$adj.P.Val < adjp & abs(rt2$logFC) > logfc,
                      ifelse(rt2$logFC > logfc, "Up", "Down"), "Not")
    p <- ggplot(rt2, aes(logFC, -log10(adj.P.Val), col = Sig)) +
      geom_point(size = 0.8) +
      scale_color_manual(values = c(Down = "blue2", Not = "grey", Up = "red2")) +
      labs(title = paste(alt, "vs", ref)) +
      theme(plot.title = element_text(size = 16, hjust = .5, face = "bold"))
    ggsave(file.path(outdir, paste0("vol_", tag, ".pdf")), p, width = 5.5, height = 4.5)
  }
  invisible(diffSig)
}

# ---- 실행 예시 ----
# run_deg(file.path(DATA_DIR, "GSE142025_twoGroups.normalize.txt"),
#         ref = "Control", alt = "DKD", tag = "DKD_vs_Control")
