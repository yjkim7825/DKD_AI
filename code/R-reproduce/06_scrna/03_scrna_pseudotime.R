# ============================================================================
# 03_scrna_pseudotime.R — scRNA 3/4: 유사시간 (PCT → dPCT 손상 궤적)
#   ⚠ monocle2(v2.40)가 최신 dplyr(group_by_ defunct)·ggplot과 충돌 →
#     안정적 대체: 손상마커 점수로 세포를 정상→손상 순서로 정렬(유사시간).
#   입력 : output/Late_with_geneset_scores.RDS (02 산출) 중 PCT·dPCT
#   출력 : 궤적(UMAP, pseudotime·celltype 색) + FN1·ALDH2 유사시간 추세
#   ※ monocle DDRTree 궤적 그림 자체가 필요하면 monocle 설치+dplyr 다운그레이드 필요.
# ============================================================================
set.seed(123)
suppressMessages({ library(Seurat); library(ggplot2); library(dplyr) })

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
OUT  <- file.path(ROOT, "R-reproduce/06_scrna/output")
Late <- readRDS(file.path(OUT, "Late_with_geneset_scores.RDS"))
Idents(Late) <- Late$celltype

## ── 근위세관 계열(PCT·dPCT)만 → 손상점수로 유사시간 ─────────────────────────
res <- subset(Late, idents = c("PCT", "dPCT"))
inj_genes <- intersect(c("VCAM1", "HAVCR2", "SPP1"), rownames(res))
inj <- colMeans(GetAssayData(res, layer = "data")[inj_genes, , drop = FALSE])
res$injury     <- as.numeric(inj)
res$pseudotime <- rank(inj, ties.method = "average") / length(inj)   # 0~1 정상→손상

## 2D 임베딩 (01에서 만든 UMAP 재사용; 없으면 PCA)
emb_name <- if ("umap" %in% Reductions(res)) "umap" else "pca"
emb <- Embeddings(res, emb_name)[, 1:2]
df <- data.frame(D1 = emb[, 1], D2 = emb[, 2],
                 pseudotime = res$pseudotime, celltype = res$celltype,
                 FN1 = FetchData(res, "FN1")[, 1], ALDH2 = FetchData(res, "ALDH2")[, 1])

## ── 궤적 그림 (Fig 7F 대응) ─────────────────────────────────────────────────
pdf(file.path(OUT, "03.trajectory_pseudotime.pdf"), width = 6, height = 5)
print(ggplot(df, aes(D1, D2, color = pseudotime)) + geom_point(size = 0.4) +
        scale_color_viridis_c() + theme_classic() +
        labs(title = "Pseudotime (PCT -> dPCT injury)", x = emb_name, y = "")); dev.off()
pdf(file.path(OUT, "03.trajectory_celltype.pdf"), width = 6, height = 5)
print(ggplot(df, aes(D1, D2, color = celltype)) + geom_point(size = 0.4) +
        scale_color_manual(values = c(PCT = "#00A087", dPCT = "#E69F00")) +
        theme_classic() + labs(title = "PCT / dPCT", x = emb_name, y = "")); dev.off()

## ── ★ FN1·ALDH2 유사시간 추세 (Fig 7G 대응) ─────────────────────────────────
pdf(file.path(OUT, "03.FN1_ALDH2_trend.pdf"), width = 7, height = 4)
print(ggplot(df, aes(pseudotime)) +
        geom_smooth(aes(y = FN1,   color = "FN1"),   method = "loess", se = TRUE) +
        geom_smooth(aes(y = ALDH2, color = "ALDH2"), method = "loess", se = TRUE) +
        scale_color_manual(values = c(FN1 = "red", ALDH2 = "#0072B2")) +
        theme_classic() + labs(y = "expression", color = "gene",
          title = "FN1 / ALDH2 along PCT->dPCT pseudotime (late DKD)")); dev.off()

# 유사시간 구간별 평균 (수치 확인용)
df$bin <- cut(df$pseudotime, breaks = 4, labels = c("Q1(정상)","Q2","Q3","Q4(손상)"))
summ <- aggregate(cbind(FN1, ALDH2) ~ bin, df, mean)
write.csv(summ, file.path(OUT, "03.pseudotime_bins_FN1_ALDH2.csv"), row.names = FALSE)
cat("[유사시간 구간별 평균]\n"); print(summ, row.names = FALSE)
cat("\n★ [03] 완료 (monocle 대체: 손상점수 정렬). 다음: 04_scrna_cellchat.R\n")
cat("  해석: 정상(Q1)→손상(Q4)로 갈수록 ALDH2 감소가 논문 Fig7G 방향.\n")
