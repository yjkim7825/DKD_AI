"""run_all.py — STEP 01~07 순차 실행 마스터 (Python 이식본).

각 STEP .py 를 순서대로 서브프로세스로 실행. 경로는 config.py(pathlib) 자동탐지.
⚠️ 일부 STEP 은 무겁다(STEP6 MR FinnGen 2.1GB, STEP7 scRNA). 또한 일부 라이브러리
   (pycombat/gseapy/scanpy 등) 미설치 시 해당 STEP 은 TODO 안내 후 넘어간다.
   필요한 STEP 만 돌리려면 STEPS 리스트에서 골라 실행하세요.
"""
import subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [
    "01_deg_gse142025/01_to_matrix.py",
    "01_deg_gse142025/02_deg.py",
    "02_rma_microarray/01_gse96804.py",
    "02_rma_microarray/02_gse30529.py",
    "02_rma_microarray/03_gse104948_104954.py",
    "03_merge_combat/01_merge_combat.py",
    "04_ml_roc/01_paper_design_lasso_svm_roc.py",
    "04_ml_roc/02_roc_plots.py",
    "05_gsea/01_gsea.py",
    "06_mr/01_mr.py",            # 무거움 (FinnGen 2.1GB)
    "07_scrna/01_scrna.py",      # 무거움 (scanpy)
]


def main():
    for s in STEPS:
        print(f"\n==================== RUN: {s} ====================")
        subprocess.run([sys.executable, str(ROOT / s)], check=False)
    print("\n[run_all] 전체 STEP 완료")


if __name__ == "__main__":
    main()
