"""STEP2 공통 : CEL → RMA 는 R 전용(affy/oligo/makecdfenv) — Python 미대응.

이식 원칙에 따라 Python 은 RMA 를 수행하지 않고, R STEP2 산출물
(OUT_DIR/*.labeled.txt)을 '입력'으로 재사용한다. 각 데이터셋 로더는 해당 매트릭스를
읽어 차원·그룹만 확인(스모크). RMA 자체가 필요하면 R 02_rma_microarray/ 를 먼저 실행.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import pandas as pd


def load_labeled(name: str) -> pd.DataFrame:
    """R RMA 산출 매트릭스 로드. (TODO: RMA 는 R 에서만 생성됨)"""
    path = config.OUT_DIR / f"{name}.labeled.txt"
    if not path.exists():
        print(f"[{name}] 매트릭스 없음 -> R STEP2 먼저 실행 필요: {path}")
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    groups = pd.Series([c.split("_")[-1] for c in df.columns]).value_counts().to_dict()
    print(f"[{name}] {df.shape[0]} genes x {df.shape[1]} samples | groups {groups}")
    return df
