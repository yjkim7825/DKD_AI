"""
io_utils.py
============
데이터 입출력 유틸 (원본 R의 read_tsv / read.table / probe 매핑 / avereps 대응).

- load_series_matrix() : GEO series_matrix.txt(.gz) 파서 (하드코딩 skip 대신 자동 탐지)
- load_platform_map()  : GPL 플랫폼 주석 파일에서 probe -> gene symbol 매핑
- avereps()            : 중복 유전자 행 평균 (limma::avereps)
- read_expr_labeled()  : '{샘플}_{그룹}' 라벨 헤더 발현행렬 로드
- write_expr_labeled() : 동일 포맷으로 저장
"""
from __future__ import annotations
import gzip
import io
import re
import numpy as np
import pandas as pd


def _open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def load_series_matrix(path) -> pd.DataFrame:
    """
    GEO series_matrix.txt(.gz) 에서 발현행렬 추출.
    '!series_matrix_table_begin' ~ '!series_matrix_table_end' 블록을 자동 인식.
    반환: index=probe ID(ID_REF), columns=GSM 샘플.
    """
    lines = []
    with _open_maybe_gz(path, "rt") as fh:
        capture = False
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                capture = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if capture:
                lines.append(line)
    if not lines:
        raise ValueError(f"series_matrix 테이블 블록을 찾지 못함: {path}")
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t")
    df = df.rename(columns={df.columns[0]: "ID_REF"})
    df["ID_REF"] = df["ID_REF"].astype(str).str.strip('"')
    df.columns = [c.strip('"') for c in df.columns]
    df = df.set_index("ID_REF")
    return df.apply(pd.to_numeric, errors="coerce")


def load_platform_map(path, id_col="ID", symbol_col="Gene Symbol",
                      skip_prefix=("#", "^", "!"), sep_symbol=" // ",
                      symbol_index=0) -> pd.Series:
    """
    GPL 플랫폼 주석 파일에서 probe ID -> gene symbol 매핑 Series 반환.
    - 주석/헤더 라인(#, ^, !)은 자동 건너뜀
    - 'AAA // BBB // ...' 형태 심볼은 symbol_index 번째 토큰 사용
      (원본 R: str_split(...)[,2] 즉 두번째 -> 필요 시 symbol_index=1 로)
    """
    # 헤더 행 자동 탐지: id_col 이 들어있는 첫 줄
    header_row = None
    with _open_maybe_gz(path, "rt") as fh:
        for i, line in enumerate(fh):
            if line.startswith(skip_prefix):
                continue
            if id_col in line.split("\t"):
                header_row = i
                break
    if header_row is None:
        raise ValueError(f"'{id_col}' 헤더를 플랫폼 파일에서 못 찾음: {path}")

    df = pd.read_csv(path, sep="\t", skiprows=header_row, dtype=str,
                     compression="gzip" if str(path).endswith(".gz") else None)
    df = df[[id_col, symbol_col]].dropna()
    df[symbol_col] = df[symbol_col].apply(
        lambda s: s.split(sep_symbol)[symbol_index].strip()
        if sep_symbol in s else s.strip()
    )
    df = df[df[symbol_col] != ""]
    df[id_col] = df[id_col].astype(str)
    return pd.Series(df[symbol_col].values, index=df[id_col].values)


def avereps(df: pd.DataFrame, gene_ids) -> pd.DataFrame:
    """limma::avereps — 동일 유전자명 행을 평균으로 축약. df: probe x sample."""
    tmp = df.copy()
    tmp.index = pd.Index(gene_ids, name="gene")
    return tmp.groupby(level=0).mean()


def read_expr_labeled(path) -> pd.DataFrame:
    """
    '{샘플}_{그룹}' 라벨 헤더 발현행렬 로드. 1열=유전자명.
    반환: index=gene, columns='{sample}_{Type}'.
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df


def write_expr_labeled(df: pd.DataFrame, path, index_name="geneNames"):
    out = df.copy()
    out.index.name = index_name
    out.to_csv(path, sep="\t")


def labels_from_columns(columns, pattern=r".*_(.+)$") -> np.ndarray:
    """'{sample}_{Type}' 열이름에서 그룹 라벨(Type) 추출."""
    rgx = re.compile(pattern)
    out = []
    for c in columns:
        m = rgx.match(str(c))
        out.append(m.group(1) if m else str(c))
    return np.array(out)
