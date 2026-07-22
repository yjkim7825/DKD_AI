# ============================================================================
# scrna_func.R — 06_scrna 공용 함수
#   · read_dge   : GSE131882 dropEst dgecounts.rds.gz → count matrix
#   · annotate_by_markers : 클러스터 → 세포타입 (마커 z-score 자동 배정)
#   ※ 함수 정의만
# ============================================================================
suppressMessages({ library(Seurat); library(Matrix) })

# Ensembl→심볼 대응표 (GSE209781 features.tsv.gz: col1=ENSG, col2=symbol)
load_ens2sym <- function(features_path) {
  ft <- read.table(features_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  v <- ft[[2]]; names(v) <- ft[[1]]; v
}

# GSE131882 dgecounts.rds.gz → 유전자×세포 count 행렬
#   구조: d$umicount$exon$all (dgCMatrix). ⚠ 이중 gzip → 두 번 해제.
#   ens2sym 주면 Ensembl ID rownames를 심볼로 변환(마커 매칭용).
read_dge <- function(path, ens2sym = NULL) {
  d  <- readRDS(gzcon(gzfile(path, "rb")))           # 이중 gzip 둘 다 해제
  uc <- d$umicount                                    # UMI 카운트 사용
  m  <- if (!is.null(uc$exon$all))   uc$exon$all      # 엑손(표준)
        else if (!is.null(uc$inex$all))   uc$inex$all # 엑손+인트론(대안)
        else if (!is.null(uc$intron$all)) uc$intron$all
        else stop("count matrix 못 찾음 (umicount$exon$all 없음)")
  # 방향 자동정렬: 유전자가 '행'이어야 함. 바코드(ACGT 12+)가 행이면 전치.
  is_bc <- function(x) length(x) && mean(grepl("^[ACGTN]{12,}$", head(x, 200))) > 0.5
  if (is_bc(rownames(m)) && !is_bc(colnames(m))) m <- Matrix::t(m)
  # Ensembl ID → 심볼 변환 (ENSG... 형태면)
  if (!is.null(ens2sym) && mean(grepl("^ENSG", head(rownames(m), 50))) > 0.5) {
    sym  <- ens2sym[sub("\\..*$", "", rownames(m))]           # 버전(.1) 제거 후 매핑
    keep <- !is.na(sym) & sym != "" & !duplicated(sym)        # 매핑됨 & 중복 심볼 첫 것만
    m <- m[keep, , drop = FALSE]; rownames(m) <- sym[keep]
  }
  as(m, "CsparseMatrix")
}

# 클러스터별 마커 z-score로 세포타입 자동 배정 (데이터셋 무관)
#   obj: 클러스터가 Idents에 있는 Seurat / mk: list(celltype = c(마커...))
annotate_by_markers <- function(obj, mk) {
  feats <- intersect(unique(unlist(mk)), rownames(obj))
  avg <- AverageExpression(obj, features = feats, assays = "RNA", layer = "data")$RNA
  colnames(avg) <- sub("^g", "", colnames(avg))    # AverageExpression가 숫자 클러스터에 'g' 접두 → 제거
  z <- t(scale(t(as.matrix(avg))))                 # 유전자별 z (클러스터 간)
  cl2ct <- sapply(colnames(avg), function(cl) {
    s <- sapply(mk, function(gs) {
      g <- intersect(gs, rownames(z)); if (!length(g)) return(-Inf); mean(z[g, cl], na.rm = TRUE)
    })
    if (all(!is.finite(s))) "Unknown" else names(s)[which.max(s)]
  })
  new <- unname(cl2ct[as.character(Idents(obj))]); new[is.na(new)] <- "Unknown"
  names(new) <- colnames(obj)                       # ★ 세포 바코드 이름 부여 (Seurat 매칭용)
  obj$celltype <- factor(new)
  Idents(obj) <- obj$celltype
  obj
}

# 표준 신장 세포 마커 (자동주석용)
KIDNEY_MARKERS <- list(
  PCT = c("LRP2","CUBN","SLC34A1","ALDOB","GATM"),
  dPCT = c("VCAM1","HAVCR2","SPP1"),
  EC = c("PECAM1","FLT1","EMCN"),
  MES = c("PDGFRB","ITGA8","ACTA2","RGS5"),
  `LOH-DCT` = c("UMOD","SLC12A1","SLC12A3","WNK1"),
  CD = c("AQP2","SLC4A1","ATP6V1G3","FOXI1"),
  PODO = c("NPHS1","NPHS2"),
  T = c("CD3D","CD3E","TRAC","NKG7"),
  B = c("CD79A","CD79B","MS4A1"),
  `Mono-Mac` = c("LYZ","CD68","C1QA","CD14","CD163"),
  Mast = c("TPSAB1","CPA3"),
  Neut = c("FCGR3B","S100A8","S100A9"),
  Plasma = c("MZB1","JCHAIN"))
