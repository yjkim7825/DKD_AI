# ============================================================================
# run_ssgsea_heatmap.R — ssGSEA 경로점수 히트맵 (논문 방식)
#   논문: "ssGSEA(GSVA v2.2.0), oxidative stress·inflammatory·EMT·cellular
#          senescence·apoptosis·autophagy. All findings validated in GSE30529."
#     H = GSE142025 (Control/Early/Late) , I = GSE30529 (Control/DKD) ← 논문 검증셋
#     경로 6개: EMT, Apoptosis, Cell aging(=cellular senescence),
#               Inflammatory response, Autophagy, Oxidative stress
#   ※ 저자가 이 코드 공유 안 함 → 논문 본문 그대로 재현
#   방법: GSVA(ssgsea) 로 샘플별 경로점수 → pheatmap (그룹 컬러바)
#   출력: 03_pathway/output/figures/ssGSEA.<데이터>.pdf
# ============================================================================
suppressMessages({ library(GSVA); library(pheatmap); library(clusterProfiler) })

ROOT  <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
PATH  <- file.path(ROOT, "R-reproduce/03_pathway")
NORM  <- file.path(ROOT, "R-reproduce/01_preprocessing/bulk/output/02_normalized")
GMT_H <- file.path(ROOT, "data/4-1. h.all.v2026.1.Hs.symbols.gmt")
OUTF  <- file.path(PATH, "output/figures")
dir.create(OUTF, showWarnings = FALSE, recursive = TRUE)

# ── 경로 6개 유전자셋 만들기 ────────────────────────────────────────────────
# 4개는 Hallmark(로컬 gmt), 2개(Cell aging·Autophagy)는 GO에서 (msigdbr, 없으면 curated)
gmtH <- read.gmt(GMT_H)
getH <- function(name) gmtH$gene[gmtH$term == name]
gs <- list(
  EPITHELIAL_MESENCHYMAL_TRANSITION = getH("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  APOPTOSIS             = getH("HALLMARK_APOPTOSIS"),
  INFLAMMATORY_RESPONSE = getH("HALLMARK_INFLAMMATORY_RESPONSE"),
  OXIDATIVE_STRESS      = getH("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
)
# Cell aging / Autophagy — msigdbr GO:BP 시도, 실패 시 curated 대표 유전자
addGO <- function() {
  ok <- requireNamespace("msigdbr", quietly = TRUE)
  if (ok) {
    m <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)
    if (!is.null(m)) {
      col <- if ("gs_name" %in% names(m)) "gs_name" else "gs_name"
      pick <- function(nm) unique(m$gene_symbol[m[[col]] == nm])
      ca <- pick("GOBP_CELL_AGING"); au <- pick("GOBP_AUTOPHAGY")
      if (length(ca) > 5 && length(au) > 5) return(list(CELL_AGING = ca, AUTOPHAGY = au))
    }
  }
  # fallback: 대표 유전자 (전체 GO셋 아님 — 근사)
  list(
    CELL_AGING = c("CDKN1A","CDKN2A","TP53","LMNB1","GLB1","SERPINE1","IGFBP3","MAP2K3",
                   "MAP2K6","ETS2","TERT","TERF2","SIRT1","FOXO3","RB1"),
    AUTOPHAGY  = c("MAP1LC3B","SQSTM1","BECN1","ATG5","ATG7","ATG12","ATG3","ATG10",
                   "ULK1","ULK2","GABARAP","GABARAPL1","WIPI2","PIK3C3","LAMP2")
  )
}
gs <- c(gs, addGO())
# 논문 행 순서
row_order <- c("EPITHELIAL_MESENCHYMAL_TRANSITION","APOPTOSIS","CELL_AGING",
               "INFLAMMATORY_RESPONSE","AUTOPHAGY","OXIDATIVE_STRESS")

# ── ssGSEA + 히트맵 함수 ────────────────────────────────────────────────────
ssgsea_heatmap <- function(exprFile, tag, group_levels, group_colors) {
  rt  <- read.table(exprFile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  mat <- as.matrix(rt)
  grp <- sapply(strsplit(colnames(mat), "_"), function(x) tail(x, 1))
  ord <- order(factor(grp, levels = group_levels))          # 그룹순 정렬
  mat <- mat[, ord]; grp <- grp[ord]

  # ssGSEA (GSVA 신·구 API 모두 대응)
  sc <- tryCatch(
    GSVA::gsva(GSVA::ssgseaParam(mat, gs)),                  # 신 API (Bioc≥3.18)
    error = function(e) GSVA::gsva(mat, gs, method = "ssgsea", verbose = FALSE)  # 구 API
  )
  sc <- sc[intersect(row_order, rownames(sc)), , drop = FALSE]

  ann <- data.frame(Group = factor(grp, levels = group_levels)); rownames(ann) <- colnames(mat)
  pdf(file.path(OUTF, paste0("ssGSEA.", tag, ".pdf")), width = 11, height = 5)
  pheatmap(sc,
           annotation_col = ann,
           annotation_colors = list(Group = group_colors),
           color = colorRampPalette(c("#4575b4","white","#d73027"))(50),
           scale = "row", cluster_cols = FALSE, cluster_rows = FALSE,
           show_colnames = FALSE,
           gaps_col = cumsum(table(factor(grp, levels = group_levels))),
           main = paste0("ssGSEA pathway score - ", tag),
           fontsize = 10, fontsize_row = 9)
  dev.off()
  cat("저장:", paste0("ssGSEA.", tag, ".pdf"), " (샘플", ncol(mat), " / 경로", nrow(sc), ")\n")
}

# ── 실행 ────────────────────────────────────────────────────────────────────
# H : GSE142025 (Control/Early/Late)
ssgsea_heatmap(file.path(NORM, "GSE142025.normalized.txt"), "GSE142025_H",
               c("Control","Early","Late"),
               c(Control="#4575b4", Early="#f4a582", Late="#d73027"))
# 96804 : GSE96804 (Control/DKD) — 학습셋(비교용)
ssgsea_heatmap(file.path(NORM, "GSE96804.normalized.txt"), "GSE96804",
               c("Control","DKD"),
               c(Control="#4575b4", DKD="#d73027"))
# I : GSE30529 (Control/DKD) — 논문 검증셋
ssgsea_heatmap(file.path(NORM, "GSE30529.normalized.txt"), "GSE30529_I",
               c("Control","DKD"),
               c(Control="#4575b4", DKD="#d73027"))

cat("\n★ ssGSEA 히트맵 완료 — output/figures/ 확인\n")
cat("※ Cell aging·Autophagy 유전자셋: msigdbr 있으면 GO 전체, 없으면 대표 유전자(근사)\n")
