# ============================================================================
# 05_scrna_stage_compare.R — 단계 의존성 비교 (플랫폼 교락 제거판)
#   ⚠ 두 데이터셋을 절대값으로 섞으면 플랫폼 차이로 교락됨(초기·말기 스케일 다름).
#   → 각 데이터셋 '안에서' PCT vs dPCT 비교 → 논문 단계 서사 직접 검증.
#     초기(GSE131882, Fig6G): dPCT에서 FN1↑ 기대
#     말기(GSE209781, Fig7G): dPCT에서 FN1↓·ALDH2↓ 기대
#   입력 : output/Early_annotated.RDS (01b), Late_annotated.RDS (01)
#   출력 : output/05.stage_PCT_dPCT.csv, 05.stage_FN1_ALDH2.pdf
# ============================================================================
set.seed(123)
suppressMessages({ library(Seurat); library(ggplot2); library(dplyr) })

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
OUT  <- file.path(ROOT, "R-reproduce/06_scrna/output")

# 한 데이터셋 안에서 PCT vs dPCT 의 FN1·ALDH2 평균 + wilcox p
within_pct_dpct <- function(rds, dataset) {
  obj <- readRDS(file.path(OUT, rds)); Idents(obj) <- obj$celltype
  if (!all(c("PCT","dPCT") %in% obj$celltype)) { message("[", dataset, "] PCT/dPCT 없음"); return(NULL) }
  pd <- subset(obj, idents = c("PCT","dPCT"))
  out <- data.frame()
  for (g in c("FN1","ALDH2")) {
    ex <- FetchData(pd, vars = c(g, "celltype"))
    m  <- tapply(ex[[g]], ex$celltype, mean)
    p  <- tryCatch(wilcox.test(ex[[g]] ~ ex$celltype)$p.value, error = function(e) NA)
    out <- rbind(out, data.frame(dataset = dataset, gene = g,
                  PCT = round(m["PCT"],3), dPCT = round(m["dPCT"],3),
                  dir = ifelse(m["dPCT"] < m["PCT"], "↓dPCT", "↑dPCT"),
                  wilcox_p = signif(p,3), row.names = NULL))
  }
  out
}

res <- rbind(
  within_pct_dpct("Early_annotated.RDS", "Early(GSE131882)"),
  within_pct_dpct("Late_annotated.RDS",  "Late(GSE209781)"))
write.csv(res, file.path(OUT, "05.stage_PCT_dPCT.csv"), row.names = FALSE)
cat("\n[단계 의존 — 데이터셋 내부 PCT vs dPCT]\n"); print(res, row.names = FALSE)

# 그림: 데이터셋 × 유전자, PCT vs dPCT 막대
long <- tidyr::pivot_longer(res, c(PCT, dPCT), names_to = "celltype", values_to = "expr")
pdf(file.path(OUT, "05.stage_FN1_ALDH2.pdf"), width = 8, height = 4)
print(ggplot(long, aes(celltype, expr, fill = celltype)) +
        geom_col() + facet_wrap(~ dataset + gene, scales = "free_y", nrow = 2) +
        scale_fill_manual(values = c(PCT = "#00A087", dPCT = "#E69F00")) +
        theme_bw() + labs(title = "PCT vs dPCT within each dataset (FN1 / ALDH2)")); dev.off()

cat("\n★ [05] 완료 (데이터셋 내부 비교 — 플랫폼 교락 제거)\n")
cat("  논문 서사: 초기 dPCT FN1↑ / 말기 dPCT FN1↓·ALDH2↓. dir 열로 실제 방향 확인.\n")
