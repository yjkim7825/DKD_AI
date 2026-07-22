#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_ml — 두 방식 구현·비교 (순수 numpy)
  방식A (논문 Fig4): LASSO(10겹CV, 10→선택) → 학습·검증 AUC>0.8 필터 → GLM/SVM 결합
  방식B (저자 코드): LASSO ∩ SVM-RFE 교집합 → 유전자별 학습/검증 ROC-AUC
입력: data/processed/data.train.txt (GSE96804), data.test.txt (GSE104948/104954)
후보 10개(DEG∩MR): ALDH2 CA2 CDKN1B CREB5 FN1 IFI44L SYTL2 TSPYL5 VNN2 XAF1
"""
import numpy as np, pandas as pd, os, json

np.random.seed(123)
BASE = "/sessions/amazing-sweet-lamport/mnt/DKD_AI/code/data/processed"
# 논문 코호트 정확 일치: 학습=GSE96804, 검증=GSE104948+104954
TRAIN = os.path.join(BASE, "data.train.paper.txt")
TEST  = os.path.join(BASE, "data.valid.paper.txt")
CANDS = ['ALDH2','CA2','CDKN1B','CREB5','FN1','IFI44L','SYTL2','TSPYL5','VNN2','XAF1']

# ---------- 데이터 ----------
def load(path):
    df = pd.read_csv(path, sep='\t', index_col=0)
    genes = [g for g in CANDS if g in df.index]
    sub = df.loc[genes]
    X = sub.T.values.astype(float)                       # 샘플 x 유전자
    y = np.array([0 if c.endswith('Control') else 1 for c in sub.columns])
    return X, y, genes

Xtr_raw, ytr, genes = load(TRAIN)
Xte_raw, yte, genes_te = load(TEST)
assert genes == genes_te, "유전자 순서 불일치"
# 표준화: 학습셋 통계로 학습·검증 둘 다 변환
mu, sd = Xtr_raw.mean(0), Xtr_raw.std(0) + 1e-8
Xtr = (Xtr_raw - mu) / sd
Xte = (Xte_raw - mu) / sd
print(f"[데이터] 학습 {Xtr.shape[0]}명(Control {sum(ytr==0)}/DKD {sum(ytr==1)}) · "
      f"검증 {Xte.shape[0]}명(Control {sum(yte==0)}/DKD {sum(yte==1)}) · 유전자 {len(genes)}개")

# ---------- 공통 함수 ----------
def auc(y, s):
    y = np.asarray(y); s = np.asarray(s, float)
    pos, neg = s[y==1], s[y==0]
    if len(pos)==0 or len(neg)==0: return np.nan
    order = np.argsort(s); ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s)+1)
    # 동점 평균순위
    _, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    csum = np.cumsum(cnt); start = csum - cnt
    avg = (start + csum + 1) / 2.0
    ranks = avg[inv]
    R1 = ranks[y==1].sum()
    return (R1 - len(pos)*(len(pos)+1)/2) / (len(pos)*len(neg))

def sigmoid(z): return 1/(1+np.exp(-np.clip(z,-30,30)))

# LASSO 로지스틱 (proximal gradient / ISTA, L1)
def lasso_logistic(X, y, lam, iters=2000, lr=0.1):
    n,p = X.shape; w = np.zeros(p); b = 0.0
    for _ in range(iters):
        z = X@w + b; pr = sigmoid(z); g = pr - y
        gw = X.T@g/n; gb = g.mean()
        w -= lr*gw; b -= lr*gb
        w = np.sign(w)*np.maximum(np.abs(w)-lr*lam, 0)   # soft-threshold
    return w, b

def cv_lambda(X, y, lambdas, k=10):
    n = len(y); idx = np.arange(n); np.random.shuffle(idx)
    folds = np.array_split(idx, k)
    mean_dev=[]; se_dev=[]
    for lam in lambdas:
        perfold=[]
        for f in folds:
            tr = np.setdiff1d(idx, f)
            w,b = lasso_logistic(X[tr], y[tr], lam)
            pr = np.clip(sigmoid(X[f]@w+b),1e-6,1-1e-6)
            perfold.append(-np.mean(y[f]*np.log(pr)+(1-y[f])*np.log(1-pr)))
        mean_dev.append(np.mean(perfold)); se_dev.append(np.std(perfold)/np.sqrt(k))
    mean_dev=np.array(mean_dev); se_dev=np.array(se_dev)
    i_min=int(np.argmin(mean_dev)); lam_min=lambdas[i_min]
    # 1-SE 규칙: 최소편차+1SE 이내에서 가장 큰 λ(=가장 인색)
    thr=mean_dev[i_min]+se_dev[i_min]
    ok=np.where(mean_dev<=thr)[0]; lam_1se=lambdas[ok.max()]   # lambdas 오름차순
    return lam_min, lam_1se, mean_dev

# 선형 SVM (hinge + L2, 경사하강) → RFE
def linsvm(X, y, C=1.0, iters=1500, lr=0.05):
    n,p = X.shape; w=np.zeros(p); b=0.0; t = 2*y-1     # ±1
    for _ in range(iters):
        m = t*(X@w+b); mask = m<1
        gw = w - C*(X[mask].T@t[mask]); gb = -C*t[mask].sum()
        w -= lr*gw/n; b -= lr*gb/n
    return w,b

def svm_cv_acc(X, y, k=5):
    n=len(y); idx=np.arange(n); np.random.shuffle(idx); folds=np.array_split(idx,k); acc=[]
    for f in folds:
        tr=np.setdiff1d(idx,f); w,b=linsvm(X[tr],y[tr])
        pred=((X[f]@w+b)>0).astype(int); acc.append((pred==y[f]).mean())
    return np.mean(acc)

def svm_rfe(X, y, gene_names, sizes=(2,3,4,5,6,7,8)):
    remaining = list(range(X.shape[1])); best=None
    # 전체에서 시작해 |w| 최소 특징 제거하며 크기별 CV 정확도 기록
    scores = {}
    cur = remaining[:]
    ranking = []
    while len(cur) > 1:
        w,_ = linsvm(X[:,cur], y)
        if len(cur) in sizes:
            scores[len(cur)] = svm_cv_acc(X[:,cur], y)
        worst = int(np.argmin(np.abs(w)))
        ranking.append(cur.pop(worst))
    if len(cur) in sizes: scores[len(cur)] = svm_cv_acc(X[:,cur], y)
    best_size = max(scores, key=scores.get)
    # best_size개 남을 때까지 제거 재현
    cur = list(range(X.shape[1]))
    while len(cur) > best_size:
        w,_ = linsvm(X[:,cur], y); cur.pop(int(np.argmin(np.abs(w))))
    return [gene_names[i] for i in cur], best_size, scores

# GLM (패널티 없는 로지스틱)
def glm(X, y, iters=3000, lr=0.1):
    n,p=X.shape; w=np.zeros(p); b=0.0
    for _ in range(iters):
        g=sigmoid(X@w+b)-y; w-=lr*(X.T@g/n); b-=lr*g.mean()
    return w,b

# ============================================================
# 방식 A — 논문 Figure 4
# ============================================================
print("\n" + "="*60 + "\n방식 A (논문 Fig4): LASSO → ROC필터 → 결합\n" + "="*60)
lambdas = np.round(np.logspace(-3, 0, 25), 5)
lam_min, lam_1se, dev = cv_lambda(Xtr, ytr, lambdas)
wL, bL = lasso_logistic(Xtr, ytr, lam_1se)                    # 인색한 1-SE λ 사용(논문형 소수 선택)
lasso_genes = [genes[i] for i in range(len(genes)) if abs(wL[i])>1e-4]
print(f"[A1 LASSO] λmin={lam_min}, λ1se={lam_1se} → 선택 {len(lasso_genes)}개: {lasso_genes}")

# ROC 필터: 학습·검증 둘 다 AUC>0.8
rows=[]
for i,g in enumerate(genes):
    if g not in lasso_genes: continue
    a_tr=auc(ytr, Xtr[:,i]); a_te=auc(yte, Xte[:,i])
    # 방향: DKD에서 높으면 그대로, 낮으면 부호 반전한 값이 진단력 (AUC는 <0.5면 1-AUC)
    a_tr=max(a_tr,1-a_tr); a_te=max(a_te,1-a_te)
    rows.append((g,a_tr,a_te,a_tr>0.8 and a_te>0.8))
rocA=pd.DataFrame(rows, columns=['gene','AUC_train','AUC_test','pass_0.8'])
finalA=list(rocA[rocA['pass_0.8']]['gene'])
print("[A2 ROC필터]\n", rocA.to_string(index=False))
print(f"→ 최종(AUC>0.8): {finalA}")

# 결합 모델 (GLM, 선형SVM)
idxF=[genes.index(g) for g in finalA]
resA={}
if idxF:
    wg,bg=glm(Xtr[:,idxF],ytr); resA['GLM']=(auc(ytr,Xtr[:,idxF]@wg+bg), auc(yte,Xte[:,idxF]@wg+bg))
    ws,bs=linsvm(Xtr[:,idxF],ytr); resA['SVM(linear)']=(auc(ytr,Xtr[:,idxF]@ws+bs), auc(yte,Xte[:,idxF]@ws+bs))
print("[A3 결합모델] (AUC train / test)")
for k,(a,b) in resA.items(): print(f"   {k:14s} train={a:.3f}  test={b:.3f}")

# ============================================================
# 방식 B — 저자 코드 (LASSO ∩ SVM-RFE)
# ============================================================
print("\n" + "="*60 + "\n방식 B (저자): LASSO ∩ SVM-RFE\n" + "="*60)
svmrfe_genes, best_size, svm_scores = svm_rfe(Xtr, ytr, genes)
print(f"[B1 SVM-RFE] 최적 크기 {best_size} → {svmrfe_genes}")
print(f"           크기별 CV정확도: { {k:round(v,3) for k,v in sorted(svm_scores.items())} }")
inter = [g for g in genes if g in lasso_genes and g in svmrfe_genes]
print(f"[B2 교집합] LASSO({len(lasso_genes)}) ∩ SVM-RFE({len(svmrfe_genes)}) = {inter}")
rowsB = []
for g in inter:
    i = genes.index(g)
    a_tr = max(auc(ytr, Xtr[:, i]), 1 - auc(ytr, Xtr[:, i]))
    a_te = max(auc(yte, Xte[:, i]), 1 - auc(yte, Xte[:, i]))
    rowsB.append((g, a_tr, a_te))
rocB = pd.DataFrame(rowsB, columns=['gene', 'AUC_train', 'AUC_test'])
print("[B3 유전자별 ROC]\n", rocB.to_string(index=False))


# ===== 저장 =====
OUT = "/sessions/amazing-sweet-lamport/mnt/DKD_AI/code/R-reproduce/05_ml/output"
os.makedirs(OUT, exist_ok=True)
rocA.to_csv(OUT + "/A_roc_filter.csv", index=False)
rocB.to_csv(OUT + "/B_intersect_roc.csv", index=False)
A_comb = {}
for k in resA:
    A_comb[k] = [round(resA[k][0], 3), round(resA[k][1], 3)]
summary = {}
summary["candidates"] = genes
summary["lambda_min"] = float(lam_min)
summary["lambda_1se"] = float(lam_1se)
summary["A_lasso"] = lasso_genes
summary["A_final_AUC08"] = finalA
summary["A_combined"] = A_comb
summary["A_roc"] = rocA.round(3).to_dict("records")
summary["B_svmrfe"] = svmrfe_genes
summary["B_best_size"] = int(best_size)
summary["B_intersect"] = inter
summary["B_roc"] = rocB.round(3).to_dict("records")
fh = open(OUT + "/compare_summary.json", "w")
json.dump(summary, fh, ensure_ascii=False, indent=2)
fh.close()
print("[saved]", OUT)
