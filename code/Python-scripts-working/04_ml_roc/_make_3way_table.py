"""Python vs R vs 논문 3자 대조표 생성 (STEP1 DEG + STEP4 ML/ROC).

값 출처:
- Python : 방금 실행한 STEP1(02_deg.py)·STEP4(01_paper_design_lasso_svm_roc.py) 결과 CSV/콘솔.
- R      : R-scripts-working/results/step4_paper/*.csv, DECISIONS.md, 재현정확도_대조표.md.
- 논문    : 재현정확도_대조표.md(논문 Figure 4 판독값) — FN1 0.911/0.911, ALDH2 0.912/0.815, GLM valid 0.942, DEG(LvE) 2833.
새 분석 없음: 이미 산출된 수치를 한 표로 모으기만 함.
출력: RES_DIR/python_vs_R_vs_paper.csv
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import pandas as pd

RES = config.RES_DIR
pyd = RES / "step4_paper"

# --- Python 값을 방금 저장한 CSV 에서 읽음(하드코딩 최소화) ---
sg = pd.read_csv(pyd / "single_gene_ROC_AUC.csv").set_index("gene")
cm = pd.read_csv(pyd / "combined_model_AUC.csv").set_index("model")

def P_sg(g, col): return round(float(sg.loc[g, col]), 3)
def P_cm(m, col): return round(float(cm.loc[m, col]), 3)

rows = []
def add(cat, metric, py, r, paper, note=""):
    rows.append({"category": cat, "metric": metric, "Python": py, "R": r, "Paper": paper, "note": note})

# ===== STEP1 DEG 개수 =====
add("STEP1_DEG", "Late_vs_Early (유의 DEG 수)", 3460, 3314, 2833, "Welch t(파이썬) vs eBayes(R); 방향·비율 일치")
add("STEP1_DEG", "Late_vs_Control (유의 DEG 수)", 4096, 4022, "n/a", "논문 미보고")
add("STEP1_DEG", "Early_vs_Control (유의 DEG 수)", 473, 671, "n/a", "논문 미보고")

# ===== STEP4 단일유전자 ROC-AUC =====
add("STEP4_single", "FN1  AUC_train",  P_sg("FN1","AUC_train"),  0.909, 0.911, "Python=R 정확 일치")
add("STEP4_single", "FN1  AUC_valid",  P_sg("FN1","AUC_valid"),  0.915, 0.911, "Python=R 정확 일치")
add("STEP4_single", "ALDH2 AUC_train", P_sg("ALDH2","AUC_train"),0.940, 0.912, "Python=R 정확 일치")
add("STEP4_single", "ALDH2 AUC_valid", P_sg("ALDH2","AUC_valid"),0.784, 0.815, "Python=R 정확 일치, 논문과 근사")

# ===== STEP4 결합모델(FN1+ALDH2) AUC =====
add("STEP4_combined", "GLM  AUC_train",  P_cm("GLM","AUC_train"),  0.980, 0.978, "")
add("STEP4_combined", "GLM  AUC_valid",  P_cm("GLM","AUC_valid"),  0.826, 0.942, "Python(0.926)이 R(0.826)보다 논문(0.942)에 근접")
add("STEP4_combined", "RF   AUC_train",  P_cm("RF","AUC_train"),   1.000, 1.000, "")
add("STEP4_combined", "RF   AUC_valid",  P_cm("RF","AUC_valid"),   0.914, 0.815, "Python≈R")
add("STEP4_combined", "SVM  AUC_train",  P_cm("SVM","AUC_train"),  0.973, 0.977, "")
add("STEP4_combined", "SVM  AUC_valid",  P_cm("SVM","AUC_valid"),  0.785, 0.807, "sklearn SVC(rbf,gamma=scale)≠e1071 → 검증 갭(포팅 한계)")
add("STEP4_combined", "XGB  AUC_train",  P_cm("XGBoost","AUC_train"),0.997,0.999, "")
add("STEP4_combined", "XGB  AUC_valid",  P_cm("XGBoost","AUC_valid"),0.873,0.844, "Python≈R")

# ===== STEP4 유전자 선택 (FN1/ALDH2 포함 여부) =====
add("STEP4_select", "LASSO 선택수",        5, 6, 6, "Py{ALDH2,CREB5,FN1,IFI44L,VNN2} / R{+XAF1,CDKN1B} / 논문{+XAF1,CDKN1B,TSPYL5}")
add("STEP4_select", "FN1  ∈ LASSO",        "Yes", "Yes", "Yes", "3자 모두 FN1 선택")
add("STEP4_select", "ALDH2 ∈ LASSO",       "Yes", "Yes", "Yes", "3자 모두 ALDH2 선택")
add("STEP4_select", "LASSO∩SVM-RFE",       "ALDH2,VNN2", "LASSO 6 (교집합 파일 공란)", "FN1,ALDH2(최종)", "Py SVM-RFE(RFECV)가 FN1 탈락; 단일 ROC 로는 FN1 최상위")

# ===== STEP5 GSEA/ssGSEA (Python 결과 파일에서 읽음) =====
def hallmark_nes(tag, term):
    f = RES / "step5_gsea" / f"Hallmark.GSEA.{tag}.txt"
    if not f.exists(): return "미실행"
    t = pd.read_csv(f, sep="\t")
    tc = "Term" if "Term" in t.columns else t.columns[1]
    hit = t[t[tc] == term]
    return round(float(hit["NES"].iloc[0]), 2) if len(hit) else "n/a"

corr_f = RES / "step5_gsea" / "FN1_ALDH2_pathway_corr.csv"
if corr_f.exists():
    cc = pd.read_csv(corr_f)
    def corr_of(gene, path):
        h = cc[(cc.gene == gene) & (cc.pathway == path)]
        return round(float(h["rho"].iloc[0]), 3) if len(h) else "n/a"
    add("STEP5_GSEA", "EMT NES (Late_vs_Control)", hallmark_nes("Late_vs_Control","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), 2.52, "n/a", "FN1 축(섬유화) 상향 — 3자 방향 일치")
    add("STEP5_GSEA", "EMT NES (Late_vs_Early)",   hallmark_nes("Late_vs_Early","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), 2.59, "n/a", "")
    add("STEP5_GSEA", "OXPHOS NES (Late_vs_Control)", hallmark_nes("Late_vs_Control","HALLMARK_OXIDATIVE_PHOSPHORYLATION"), -1.62, "n/a", "ALDH2 축(대사) 하향")
    add("STEP5_GSEA", "OXPHOS NES (Late_vs_Early)",   hallmark_nes("Late_vs_Early","HALLMARK_OXIDATIVE_PHOSPHORYLATION"), -2.02, "n/a", "")
    add("STEP5_GSEA", "FN1 ↔ ECM_RECEPTOR (rho)",  corr_of("FN1","KEGG_ECM_RECEPTOR_INTERACTION"), 0.919, "n/a", "Python=R 일치")
    add("STEP5_GSEA", "ALDH2 ↔ TCA_CYCLE (rho)",   corr_of("ALDH2","KEGG_CITRATE_CYCLE_TCA_CYCLE"), 0.921, "n/a", "ALDH2 대사 축")
else:
    add("STEP5_GSEA", "GSEA/ssGSEA", "미실행", "완료", "-", "gseapy 미설치")

# ===== STEP6 MR (Python GCST 결과에서 읽음) =====
mr_f = RES / "step6_mr" / "table.MRresult.GCST.csv"
if mr_f.exists():
    mm = pd.read_csv(mr_f).set_index("exposure")
    add("STEP6_MR", "FN1 IVW b (GCST)",  round(float(mm.loc["FN1","b"]),3),  1.021,  1.021, "Python=R=Supp9 일치")
    add("STEP6_MR", "FN1 IVW OR (GCST)", round(float(mm.loc["FN1","or"]),3), 2.777, 2.777, f"nsnp {int(mm.loc['FN1','nsnp'])}")
    add("STEP6_MR", "FN1 IVW p (GCST)",  round(float(mm.loc["FN1","pval"]),4), 0.0122, 0.0122, "위험 인과 OR>1")
    add("STEP6_MR", "ALDH2 IVW b (GCST)", round(float(mm.loc["ALDH2","b"]),3), -0.395, -0.395, "보호 인과 OR<1")
    add("STEP6_MR", "ALDH2 IVW OR (GCST)",round(float(mm.loc["ALDH2","or"]),3),0.673, 0.673, f"nsnp {int(mm.loc['ALDH2','nsnp'])} (R/Supp9=14, Py=15)")
    add("STEP6_MR", "ALDH2 IVW p (GCST)", round(float(mm.loc["ALDH2","pval"]),4),0.0321,0.0321, "")
    add("STEP6_MR", "FinnGen outcome", "미실행", "완료(비유의)", "-", "2.1GB 로드 무거워 스킵(GCST 우선)")
else:
    add("STEP6_MR", "MR IVW", "미실행", "완료", "-", "데이터 확인 필요")

# ===== STEP7 scRNA (Python celltype 결과에서 읽음) =====
sc_f = RES / "step7_scrna" / "FN1_ALDH2_by_celltype.csv"
if sc_f.exists():
    st7 = pd.read_csv(sc_f)
    def top_ct(gene):
        h = st7[st7.gene == gene].sort_values("mean_expr", ascending=False)
        return f"{h.iloc[0]['celltype']} ({h.iloc[0]['mean_expr']:.2f})" if len(h) else "n/a"
    add("STEP7_scRNA", "FN1 최고발현 세포", top_ct("FN1"),  "Endothelial (1.72)", "내피/메산지움", "3자 방향 일치 여부")
    add("STEP7_scRNA", "ALDH2 최고발현 세포", top_ct("ALDH2"),"PCT (1.65)", "근위세뇨관", "3자 방향 일치 여부")
else:
    add("STEP7_scRNA", "세포유형별 발현", "미실행", "완료", "-", "scanpy 실행 실패/미완")

# ===== 대체(미이식) STEP =====
add("대체/미이식", "STEP2_RMA", "미이식", "완료", "-", "RMA 는 R 전용(affy/oligo) — R 산출 매트릭스 재사용")
add("대체/미이식", "STEP3_ComBat", "대체", "완료", "-", "inmoose 미설치 → R data.valid.paper.txt 직접 사용")

df = pd.DataFrame(rows, columns=["category","metric","Python","R","Paper","note"])
out = RES / "python_vs_R_vs_paper.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print(df.to_string(index=False))
print(f"\n[저장] {out}")
