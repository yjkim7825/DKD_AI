"""
config.py
=========
경로/파라미터 중앙 설정. 원본 R 스크립트의 하드코딩된 setwd("G:\\187geneMR\\...") 를 대체.
로컬 환경에 맞게 이 파일만 수정하면 전체 파이프라인이 동작한다.
"""
from pathlib import Path

# ---- 기본 디렉터리 -----------------------------------------------------
# pipeline_port/ 아래 data/, results/ 를 기본 사용
BASE_DIR   = Path(__file__).resolve().parent.parent      # .../pipeline_port
DATA_DIR   = BASE_DIR / "data"        # 원본/중간 데이터 (GEO 다운로드 위치)
RESULT_DIR = BASE_DIR / "results"     # 산출물
RESULT_DIR.mkdir(parents=True, exist_ok=True)

# ---- 그룹 라벨 규약 ----------------------------------------------------
# 열 이름은 '{샘플ID}_{그룹}' 형식. (data preprocessing 2.R 와 동일)
CONTROL_LABEL = "Control"
CASE_LABEL    = "DKD"          # 2군 비교 시 실험군. Early/Late 통합 시 "DKD" 로 라벨링
EARLY_LABEL   = "Early"
LATE_LABEL    = "Late"

# ---- DEG 필터 (differential expression analysis.R 동일) ----------------
LOGFC_FILTER = 0.585          # |log2FC| 임계 (=1.5배)
ADJP_FILTER  = 0.05           # BH adj.P 임계

# ---- 머신러닝 (machine learning modeling 1.R 동일) ---------------------
RANDOM_SEED  = 123
CV_FOLDS     = 10             # LASSO CV
SVM_RFE_SIZES = [2, 3, 4, 5, 6, 7, 8]   # SVM-RFE 후보 특징 수
ROC_BOOTSTRAP = 2000          # ci.auc(method="bootstrap") 대응 반복수
