# ============================================================
# DKD_AI 파이프라인 필수 패키지 일괄 설치 스크립트
# 새 PC에서 이 파일을 통째로 source 하거나 실행하세요.
#   실행: Rscript install_packages.R
# ============================================================

# --- 0. 설치 도구 준비 (BiocManager, remotes) ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes",     quietly = TRUE)) install.packages("remotes")

# --- 1. CRAN 패키지 ---
cran_pkgs <- c(
  "here",       # 경로 자동탐지 (config.R)
  "ggplot2",    # 그림
  "Seurat",     # STEP7 단일세포 (의존성 많음, 설치 오래 걸림)
  "harmony",    # STEP7 배치보정
  "patchwork",  # STEP7 그림 배치
  "ggrepel"     # STEP7 라벨
)
install.packages(cran_pkgs)

# --- 2. Bioconductor 패키지 ---
bioc_pkgs <- c(
  "limma",                                   # STEP1/2 DEG
  "sva",                                      # STEP3 ComBat 배치보정
  "oligo", "affy", "affyio",                  # STEP2 CEL(RMA) 처리
  "GEOquery",                                 # STEP2 GEO 메타데이터
  "org.Hs.eg.db",                             # 유전자 심볼 매핑 (전반)
  "hgu133a2.db",                              # STEP2 GSE30529 칩 어노테이션
  "clusterProfiler", "GSVA", "GSEABase",
  "fgsea", "enrichplot", "DOSE",              # STEP5 GSEA 계열
  "Rsamtools", "VariantAnnotation", "GenomeInfoDb"  # STEP6 MR: VCF 처리
)
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

# --- 3. GitHub 패키지 (TwoSampleMR: STEP6 MR 핵심) ---
# ⚠️ GitHub는 로그인 안 하면 시간당 요청 60회 제한이 있어,
#    여러 개를 연달아 깔면 rate limit(403)에 걸릴 수 있습니다.
#    걸리면 약 1시간 뒤 자동으로 풀리니 그때 이 부분만 다시 실행하세요.
remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never")

message("\n[완료] 패키지 설치 스크립트 끝. 에러가 난 게 있으면 그 패키지만 따로 확인하세요.")
