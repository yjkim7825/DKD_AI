"""STEP1-1 : GSE142025 per-sample txt → 발현매트릭스 (R 01_to_matrix.R 이식).

36개 샘플별 txt(.gz) (2열: Symbol, value; 이미 log 스케일)를 Symbol 교집합으로 병합.
그룹 = 파일명 접두사 N=Control / B=Early / A=Late. 출력 컬럼 = {sample}_{group}.
출력: OUT_DIR/GSE142025.labeled.txt
"""
import sys, pathlib, gzip, re
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import pandas as pd

PREFIX2GROUP = {"N": config.GROUP_CONTROL, "B": config.GROUP_EARLY, "A": config.GROUP_LATE}


def read_sample(path: pathlib.Path) -> pd.Series:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=0)
    df.columns = ["Symbol", "value"][: df.shape[1]]
    s = df.groupby("Symbol")["value"].mean()          # 중복 심볼 평균
    return s


def main():
    files = sorted(p for p in config.DIR_GSE142025.iterdir()
                   if re.search(r"\.txt(\.gz)?$", p.name))
    print(f"[GSE142025] 파일 수: {len(files)}")
    series, cols = [], []
    for p in files:
        gsm = re.sub(r"\..*$", "", p.name)
        grp = PREFIX2GROUP.get(gsm[0], "NA")
        s = read_sample(p); s.name = f"{gsm}_{grp}"
        series.append(s); cols.append(s.name)
    mat = pd.concat(series, axis=1, join="inner")      # Symbol 교집합
    mat.index.name = "geneNames"
    out = config.OUT_DIR / "GSE142025.labeled.txt"
    mat.to_csv(out, sep="\t")
    print(f"[GSE142025] {mat.shape[0]} genes x {mat.shape[1]} samples -> {out}")


if __name__ == "__main__":
    main()
