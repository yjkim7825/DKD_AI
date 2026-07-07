"""
smoke_test.py
=============
GEO 데이터/네트워크 없이 합성데이터로 파이프라인 전 구간을 검증한다.
- 실제 데이터 없이도 코드가 end-to-end 로 도는지 확인용
- 통계 검증: (a) quantile_normalize 후 열별 분포 동일, (b) limma_two_group 이
  주입한 참(true) DEG 를 회수, (c) ComBat 이 배치 평균차를 제거.

실행:  python smoke_test.py
필요:  numpy, pandas, scipy, scikit-learn (없으면 해당 단계만 건너뜀)
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from deg_utils import quantile_normalize, combat, limma_two_group, bh_adjust


def make_expr(n_genes=500, n_ctrl=20, n_case=20, n_de=30, seed=0):
    rng = np.random.default_rng(seed)
    base = rng.normal(8, 2, size=(n_genes, 1))
    ctrl = base + rng.normal(0, 0.4, size=(n_genes, n_ctrl))
    case = base + rng.normal(0, 0.4, size=(n_genes, n_case))
    # 앞의 n_de 유전자에 진짜 차등발현 주입
    effect = np.concatenate([rng.choice([-1, 1], n_de) * rng.uniform(1.0, 2.0, n_de),
                             np.zeros(n_genes - n_de)])
    case = case + effect[:, None]
    genes = [f"G{i:04d}" for i in range(n_genes)]
    cols = [f"C{i}_Control" for i in range(n_ctrl)] + [f"D{i}_DKD" for i in range(n_case)]
    df = pd.DataFrame(np.hstack([ctrl, case]), index=genes, columns=cols)
    return df, set(np.array(genes)[:n_de][effect[:n_de] != 0])


def test_quantile():
    rng = np.random.default_rng(1)
    X = rng.normal(5, 3, size=(200, 6)) + np.arange(6)  # 열마다 다른 shift
    Q = quantile_normalize(X)
    col_means = Q.mean(axis=0)
    ok = np.allclose(col_means, col_means[0], atol=1e-6)
    print(f"[quantile] 열 평균 동일화: {'OK' if ok else 'FAIL'} "
          f"(spread={col_means.max()-col_means.min():.2e})")
    return ok


def test_limma():
    df, true_de = make_expr()
    labels = np.array([c.split("_")[1] for c in df.columns])
    res = limma_two_group(df.values, labels, ref="Control", alt="DKD")
    tab = pd.DataFrame(res, index=df.index)
    called = set(tab[(tab["adj.P.Val"] < 0.05) & (tab["logFC"].abs() > 0.585)].index)
    recall = len(called & true_de) / len(true_de)
    precision = len(called & true_de) / max(len(called), 1)
    print(f"[limma ] 참 DEG {len(true_de)}개 중 회수 recall={recall:.2f}, "
          f"precision={precision:.2f}, called={len(called)}")
    return recall > 0.8 and precision > 0.8


def test_combat():
    df, _ = make_expr(seed=2)
    # 인위적 배치효과: 뒤쪽 절반 샘플에 +3 shift
    X = df.values.copy()
    batch = np.array([0] * (X.shape[1] // 2) + [1] * (X.shape[1] - X.shape[1] // 2))
    X[:, batch == 1] += 3.0
    before = abs(X[:, batch == 0].mean() - X[:, batch == 1].mean())
    Xc = combat(X, batch)
    after = abs(Xc[:, batch == 0].mean() - Xc[:, batch == 1].mean())
    print(f"[combat] 배치 평균차 {before:.3f} -> {after:.3f}")
    return after < 0.1


def test_ml():
    try:
        from sklearn.linear_model import LogisticRegressionCV  # noqa
    except Exception:
        print("[ml    ] scikit-learn 미설치 -> 건너뜀")
        return None
    import p03_feature_selection as fs
    import tempfile, os
    df, true_de = make_expr(n_genes=60, n_de=8, seed=3)
    with tempfile.TemporaryDirectory() as d:
        train = os.path.join(d, "train.txt")
        df.to_csv(train, sep="\t")
        inter = fs.main(train, None, d)
    print(f"[ml    ] 선택 유전자 {len(inter)}개 중 참 DEG 포함: "
          f"{len(set(inter) & true_de)}/{len(inter)}")
    return len(inter) > 0


if __name__ == "__main__":
    results = {
        "quantile_normalize": test_quantile(),
        "limma_two_group": test_limma(),
        "combat": test_combat(),
        "feature_selection(ML)": test_ml(),
    }
    print("\n=== SMOKE TEST 요약 ===")
    for k, v in results.items():
        tag = "SKIP" if v is None else ("PASS" if v else "FAIL")
        print(f"  {tag:4}  {k}")
