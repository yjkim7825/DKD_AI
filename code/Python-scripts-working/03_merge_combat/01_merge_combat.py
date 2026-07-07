"""STEP3 : 병합 + ComBat (R 01_merge_combat.R 이식).

훈련 data.train = ComBat(GSE96804 + GSE104948)  [사구체]
검증 data.test  = ComBat(GSE104954 + GSE30529)  [세뇨관]
공통 유전자 교집합 → 데이터셋 접두사 cbind → ComBat(batch=데이터셋).
라이브러리: ComBat → pycombat (inmoose.pycombat 또는 combat.pycombat).
출력: RES_DIR/data.train.txt, data.test.txt  (R 산출물 미덮어씀 — Python 결과는 RES_DIR)
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import numpy as np
import pandas as pd

# ComBat 구현 선택 (설치된 것 사용). 없으면 TODO.
def _get_combat():
    try:
        from inmoose.pycombat import pycombat_norm            # 최신
        return lambda dat, batch: pycombat_norm(dat, batch)
    except Exception:
        pass
    try:
        from combat.pycombat import pycombat                  # 구 API
        return lambda dat, batch: pycombat(dat, batch)
    except Exception:
        return None

TRAIN_SETS = {"GSE96804": "GSE96804.labeled.txt",  "GSE104948": "GSE104948.labeled.txt"}
TEST_SETS  = {"GSE104954": "GSE104954.labeled.txt", "GSE30529":  "GSE30529.labeled.txt"}


def read_labeled(name, fname):
    df = pd.read_csv(config.OUT_DIR / fname, sep="\t", index_col=0)
    df.columns = [f"{name}_{c}" for c in df.columns]         # 데이터셋 접두사
    return df


def merge_combat(sets, tag):
    mats = [read_labeled(n, f) for n, f in sets.items()]
    common = sorted(set.intersection(*[set(m.index) for m in mats]))
    all_tab = pd.concat([m.loc[common] for m in mats], axis=1)
    batch = np.concatenate([[i] * m.shape[1] for i, m in enumerate(mats)])
    print(f"[{tag}] 공통 유전자 {len(common)} | {all_tab.shape[1]} samples, batch {len(set(batch))}")
    combat = _get_combat()
    if combat is None:
        print(f"[{tag}] TODO: pycombat 미설치 → 'pip install inmoose' 후 재실행. preNorm 만 저장.")
        all_tab.to_csv(config.RES_DIR / f"{tag}.preNorm.txt", sep="\t")
        return
    corrected = combat(all_tab, list(batch))
    corrected = pd.DataFrame(corrected, index=all_tab.index, columns=all_tab.columns)
    corrected.index.name = "geneNames"
    corrected.to_csv(config.RES_DIR / f"{tag}.txt", sep="\t")
    grp = pd.Series([c.split("_")[-1] for c in corrected.columns]).value_counts().to_dict()
    print(f"[{tag}] 저장 -> {tag}.txt | groups {grp}")


def main():
    merge_combat(TRAIN_SETS, "data.train")
    merge_combat(TEST_SETS,  "data.test")


if __name__ == "__main__":
    main()
