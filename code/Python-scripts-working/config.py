"""config.py — 파이프라인 중앙 설정 (경로는 여기 한 곳에서만 정의).

R 쪽 config.R 와 1:1 동일 값. pathlib 로 루트 자동탐지 → 스크립트가 하위 STEP 폴더에
있어도 항상 같은 루트를 참조(../data 깊이 문제 없음).
▶ 사용자는 보통 수정 불필요. 데이터 위치가 다르면 DATA_ROOT 한 줄만 바꾸세요.
"""
from pathlib import Path

# ---- 루트 (자동탐지: 이 파일 위치 기준) ----
CODE_ROOT = Path(__file__).resolve().parent                 # = Python-scripts-working
DATA_ROOT = (CODE_ROOT.parent / "data").resolve()           # = code/data  (R 파이프라인과 공유)

# ---- 산출물/작업 폴더 ----
OUT_DIR = DATA_ROOT / "processed"        # 중간 매트릭스(R RMA 산출물 재사용, 공유)
RES_DIR = CODE_ROOT / "results"          # Python 분석 결과(R results 와 분리)
SCRATCH_DIR = CODE_ROOT / ".scratch"     # 임시 작업(비출력, 하드코딩 없음)
for _d in (OUT_DIR, RES_DIR, SCRATCH_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---- 번호 붙은 원본 RAW 폴더 (읽기 전용) ----
DIR_GSE96804  = DATA_ROOT / "1-1. GSE96804_RAW"
DIR_GSE104948 = DATA_ROOT / "1-2. GSE104948_RAW"
DIR_GSE104954 = DATA_ROOT / "1-3. GSE104954_RAW"
DIR_GSE142025 = DATA_ROOT / "1-4. GSE142025_RAW"
DIR_GSE30529  = DATA_ROOT / "1-5. GSE30529_RAW"
DIR_GSE131882 = DATA_ROOT / "1-6. GSE131882_RAW"
DIR_GSE209781 = DATA_ROOT / "1-7. GSE209781_RAW"
DIR_GSE266146 = DATA_ROOT / "1-8. GSE266146_RAW"
DIR_SUPP      = DATA_ROOT / "6. Article related data"

# ---- 그룹 라벨 규약 (R 과 동일) ----
GROUP_CONTROL = "Control"
GROUP_EARLY   = "Early"
GROUP_LATE    = "Late"
GROUP_DKD     = "DKD"

# ---- DEG 필터 (논문 동일) ----
LOGFC_FILTER = 0.585
ADJP_FILTER  = 0.05

if __name__ == "__main__":
    print(f"[config] CODE_ROOT = {CODE_ROOT}")
    print(f"[config] DATA_ROOT = {DATA_ROOT}")
    print(f"[config] OUT_DIR   = {OUT_DIR}")
    print(f"[config] RES_DIR   = {RES_DIR}")
