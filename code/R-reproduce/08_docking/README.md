# 08_docking — 분자 도킹 (레스베라트롤 ↔ FN1·ALDH2)

논문: 레스베라트롤(PubChem **CID 445154**)을 **FN1(PDB 3m7p)** · **ALDH2(PDB 8dr9)** 에 도킹 → CB-Dock2로 Vina 점수 **FN1 −7.0 / ALDH2 −9.8 kcal/mol** (< −5.0 = 유의미한 결합).

두 가지 재현 방법을 모두 제공한다.

---

## 방법 ① AutoDock Vina (로컬 코드 재현) — `dock_vina.py`

CB-Dock2와 같은 blind docking을 오픈소스 Vina로 로컬 실행. **코드로 기록·재현 가능**.
※ 절대 점수는 CB-Dock2와 다를 수 있음(알고리즘·포켓탐지 차이) — "강한 결합(−5 미만)" 경향이 재현되면 성공.

### 설치 (한 번만)
```bash
# conda 권장 (Windows/Mac/Linux 공통)
conda create -n dock -c conda-forge python=3.10 vina openbabel -y
conda activate dock
```
- Vina 단독: https://github.com/ccsb-scripps/AutoDock-Vina/releases
- Open Babel: https://openbabel.org

### 실행
```bash
cd 08_docking
python dock_vina.py
```
스크립트가 자동으로:
1. 레스베라트롤 SDF (PubChem CID 445154) 다운로드 → 수소 추가 → PDBQT
2. FN1(3m7p)·ALDH2(8dr9) PDB 다운로드 → 물 제거 → 수소 추가 → PDBQT
3. 수용체 전체를 감싸는 blind docking 박스 계산
4. Vina 도킹 → 최적 점수·포즈 저장

### 출력 (`output/`)
| 파일 | 내용 |
| --- | --- |
| `dock_summary.csv` | 표적별 Vina 점수 + 논문값 대조 |
| `<gene>_<pdb>_out.pdbqt` | 도킹 포즈 (PyMOL/Chimera로 시각화) |
| `<gene>_<pdb>_vina.log` | Vina 로그 |

---

## 방법 ② CB-Dock2 (웹툴, 논문과 동일) — 코드 없음

논문이 실제로 쓴 방법. 웹에서 클릭으로 수행.

1. **구조 준비**
   - 리간드: PubChem에서 Resveratrol(CID 445154) → **3D SDF** 다운로드
   - 수용체: RCSB PDB에서 **3m7p**(FN1), **8dr9**(ALDH2) → PDB 다운로드
2. **CB-Dock2 접속**: https://cadd.labshare.cn/cb-dock2/ (또는 검색 "CB-Dock2")
3. 수용체 PDB 업로드 + 리간드 SDF 업로드 → **Dock** 클릭
4. 결과: 자동 포켓 탐지 + **Vina Score** (음수 클수록 강함) + 3D 결합 포즈
5. 논문값과 대조: FN1 ≈ −7.0, ALDH2 ≈ −9.8 kcal/mol

- 장점: 논문과 정확히 같은 방법 / 단점: 웹 클릭이라 코드 기록 불가

---

## 참고
- 논문 처리: 수소 원자 추가 + 물 분자 제거 (두 방법 다 반영)
- 판정 기준: Vina 점수 **< −5.0** = 유의미한 결합 친화도
- 상세·해석: `../docs/08_도킹_논문vs재현.md`
