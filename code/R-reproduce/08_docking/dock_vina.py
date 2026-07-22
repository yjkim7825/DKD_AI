#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dock_vina.py — 분자 도킹 (AutoDock Vina) : 레스베라트롤 ↔ FN1·ALDH2
  논문 재현: 리간드=Resveratrol(PubChem CID 445154),
            표적=FN1(PDB 3m7p) · ALDH2(PDB 8dr9)  (논문 지정 구조)
  CB-Dock2(웹, 논문)와 같은 blind docking을 로컬 Vina로 수행.
  ※ 알고리즘·포켓탐지가 CB-Dock2와 달라 점수 절대값은 다를 수 있음(경향은 동일).
  ※ 네트워크·Vina·OpenBabel 설치가 필요 → 사용자 PC에서 실행 (샌드박스 X).

필요 도구 (설치는 README 참고):
  - AutoDock Vina  (명령어 `vina`)
  - Open Babel     (명령어 `obabel`)
  - Python 3 (표준 라이브러리만 사용)

실행:
  python dock_vina.py
결과:
  output/dock_summary.csv   (표적별 최적 Vina 점수)
  output/<target>_out.pdbqt (도킹 포즈)
"""
import os, sys, subprocess, urllib.request, shutil, re

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "output"); os.makedirs(OUT, exist_ok=True)

LIGAND_CID = 445154                      # Resveratrol
TARGETS = {"FN1": "3m7p", "ALDH2": "8dr9"}   # 논문 지정 PDB
EXHAUST = 16                             # Vina exhaustiveness (높을수록 정밀·느림)

def need(cmd):
    if shutil.which(cmd) is None:
        sys.exit(f"[에러] '{cmd}' 없음 — README의 설치 안내를 따르세요 (vina, obabel).")

def fetch(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        print(f"  (있음) {os.path.basename(path)}"); return
    print(f"  다운로드 {url}")
    urllib.request.urlretrieve(url, path)

def prep_ligand():
    sdf = os.path.join(OUT, f"resveratrol_{LIGAND_CID}.sdf")
    pdbqt = os.path.join(OUT, "ligand.pdbqt")
    fetch(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{LIGAND_CID}/SDF?record_type=3d", sdf)
    # 수소 추가 + 3D + PDBQT 변환 (논문: 수소 추가)
    subprocess.run(["obabel", sdf, "-O", pdbqt, "-h", "--partialcharge", "gasteiger"], check=True)
    return pdbqt

def prep_receptor(pdb_id):
    pdb  = os.path.join(OUT, f"{pdb_id}.pdb")
    clean= os.path.join(OUT, f"{pdb_id}_clean.pdb")
    pdbqt= os.path.join(OUT, f"{pdb_id}.pdbqt")
    fetch(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", pdb)
    # 물(HOH)·헤테로원자 제거, 단백질 ATOM만 (논문: 물 분자 제거)
    with open(pdb) as f, open(clean, "w") as g:
        for line in f:
            if line.startswith("ATOM"):
                g.write(line)
        g.write("END\n")
    # 수소 추가 + 강체 수용체 PDBQT (-xr)
    subprocess.run(["obabel", clean, "-O", pdbqt, "-h", "-xr", "--partialcharge", "gasteiger"], check=True)
    return pdbqt

def box_from_pdbqt(pdbqt, pad=5.0, maxsize=30.0):
    """blind docking: 수용체 전체를 감싸는 박스 (중심=무게중심, 크기=범위+pad, 상한 maxsize)."""
    xs, ys, zs = [], [], []
    with open(pdbqt) as f:
        for l in f:
            if l.startswith(("ATOM", "HETATM")):
                xs.append(float(l[30:38])); ys.append(float(l[38:46])); zs.append(float(l[46:54]))
    cx, cy, cz = (min(xs)+max(xs))/2, (min(ys)+max(ys))/2, (min(zs)+max(zs))/2
    sx = min(max(xs)-min(xs)+pad, maxsize)
    sy = min(max(ys)-min(ys)+pad, maxsize)
    sz = min(max(zs)-min(zs)+pad, maxsize)
    return (cx, cy, cz, sx, sy, sz)

def run_vina(receptor, ligand, box, out_pose, log):
    cx, cy, cz, sx, sy, sz = box
    cmd = ["vina", "--receptor", receptor, "--ligand", ligand,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", f"{sx:.1f}", "--size_y", f"{sy:.1f}", "--size_z", f"{sz:.1f}",
           "--exhaustiveness", str(EXHAUST), "--out", out_pose]
    with open(log, "w") as lg:
        subprocess.run(cmd, check=True, stdout=lg, stderr=subprocess.STDOUT)
    # 최상위 포즈 점수(첫 REMARK VINA RESULT) 파싱
    best = None
    with open(out_pose) as f:
        for l in f:
            m = re.search(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)", l)
            if m: best = float(m.group(1)); break
    return best

def main():
    for c in ("vina", "obabel"): need(c)
    print("[1] 리간드(레스베라트롤) 준비"); lig = prep_ligand()
    rows = [("target", "pdb", "vina_score_kcal_mol", "paper_score")]
    paper = {"FN1": -7.0, "ALDH2": -9.8}
    for gene, pid in TARGETS.items():
        print(f"[2] 수용체 {gene} ({pid}) 준비"); rec = prep_receptor(pid)
        box = box_from_pdbqt(rec)
        print(f"[3] 도킹 {gene} … (blind, exhaustiveness={EXHAUST})")
        pose = os.path.join(OUT, f"{gene}_{pid}_out.pdbqt")
        log  = os.path.join(OUT, f"{gene}_{pid}_vina.log")
        score = run_vina(rec, lig, box, pose, log)
        print(f"    → {gene} Vina 점수: {score} kcal/mol (논문 {paper[gene]})")
        rows.append((gene, pid, score, paper[gene]))
    with open(os.path.join(OUT, "dock_summary.csv"), "w") as f:
        for r in rows: f.write(",".join(map(str, r)) + "\n")
    print("\n★ 도킹 완료 → output/dock_summary.csv, *_out.pdbqt")
    print("  해석: Vina 점수 < -5.0 = 유의미한 결합. 음수 클수록 강함.")
    print("  ※ CB-Dock2(논문)와 알고리즘이 달라 절대값 차이 가능 — '강한 결합' 경향이 재현되면 성공.")

if __name__ == "__main__":
    main()
