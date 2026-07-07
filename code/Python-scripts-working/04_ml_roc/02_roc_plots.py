"""STEP4 ROC 그림 (R 02_roc_plots.R 이식). matplotlib + sklearn.

단일유전자 ROC 4개(train/valid × FN1/ALDH2, AUC 라벨, 빨강+대각선) + 결합모델 ROC 2개.
데이터 재수신 없음: GSE96804(train) + ComBat(104948+104954)(valid). 출력: RES_DIR/step4_paper/*.pdf
"""
import sys, pathlib, importlib.util
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import numpy as np
import pandas as pd

OUT = config.RES_DIR / "step4_paper"; OUT.mkdir(parents=True, exist_ok=True)

# STEP4 본 스크립트의 load_valid / y_of 재사용
_p = pathlib.Path(__file__).with_name("01_paper_design_lasso_svm_roc.py")
_spec = importlib.util.spec_from_file_location("s4", _p)
s4 = importlib.util.module_from_spec(_spec)
try:
    _spec.loader.exec_module(s4)   # import 시 main 은 실행 안 됨(__main__ 가드)
except Exception as e:
    print("[roc_plots] STEP4 모듈 로드 경고:", e)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, roc_auc_score

    train = pd.read_csv(config.OUT_DIR / "GSE96804.labeled.txt", sep="\t", index_col=0)
    valid = s4.load_valid()
    ytr, yva = s4.y_of(train), s4.y_of(valid)

    def draw_single(df, y, gene, tag):
        v = df.loc[gene].to_numpy(); a = roc_auc_score(y, v)
        if a < 0.5:            # pROC 자동 방향: 음의 연관(보호유전자)이면 예측자 부호 반전
            v = -v; a = 1 - a
        fpr, tpr, _ = roc_curve(y, v)
        plt.figure(figsize=(3.6, 3.6)); plt.plot(fpr, tpr, "r", lw=2.5)
        plt.plot([0, 1], [0, 1], "k:", lw=1)
        plt.title(f"{gene} ({tag})"); plt.xlabel("1 - Specificity"); plt.ylabel("Sensitivity")
        plt.text(0.45, 0.1, f"AUC: {a:.3f}", color="red")
        plt.tight_layout(); plt.savefig(OUT / f"{tag}_ROC.{gene}.pdf"); plt.close()

    for g in ["FN1", "ALDH2"]:
        draw_single(train, ytr, g, "train"); draw_single(valid, yva, g, "valid")

    print("[roc_plots] 단일유전자 ROC 4개 저장 -> step4_paper/")
    # 결합모델 ROC 는 01 스크립트의 모델 재학습 필요 → 요약만(상세 곡선은 확장 TODO).


if __name__ == "__main__":
    main()
