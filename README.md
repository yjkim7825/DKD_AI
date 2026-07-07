# DKD Multi-omics 재현 (P25: FN1 & ALDH2)

당뇨병성 신장질환(DKD) 다중오믹스 바이오마커 논문
*Lin et al., "Multi-omics and machine learning identify FN1 and ALDH2 as diagnostic
biomarkers and therapeutic targets in early and late DKD", Renal Failure 2025* (P25)
의 파이프라인을 **로컬에서 재현하고 Python으로 이식**한 저장소.

## 재현 결과 (요약)
- STEP4 ROC: FN1 train 0.909 / valid 0.915 (논문 0.911 / 0.911) — 거의 정확 재현
- STEP6 MR: FN1 b=1.021 (위험), ALDH2 b=-0.395 (보호) — 저자 Supp9 완전 일치
- STEP5 GSEA: FN1=EMT/ECM, ALDH2=산화적인산화/대사
- 전체 대조표: `pipeline_port/docs/재현정확도_대조표.md`

## 폴더 구조
- `code/R-scripts-working/`  — 재현 R 파이프라인(STEP1~7, config 중앙화)
- `code/Python-scripts-working/`  — Python 이식본(R과 대칭 구조)
- `code/R-scripts-for-pipeline-reproducibility/`  — 저자 원본 코드(참조)
- `pipeline_port/docs/`  — 데이터 근거표 · 재현정확도 대조표
- `code/data/processed/`  — **축약 데이터**(재현용 매트릭스 3종)
- `code/data/*.gmt`  — Hallmark / KEGG 유전자셋

## 데이터
- **원본 RAW/GWAS(약 12GB)는 저장소에 포함되지 않음.** 출처·다운로드는
  `pipeline_port/docs/데이터_근거표.md` 참조(GEO / OpenGWAS / FinnGen / GWAS Catalog / MSigDB / Nephroseq).
- 저장소에는 **축약 데이터**만 포함: STEP4 ML 입력 매트릭스(train/test/valid) + Hallmark·KEGG gmt + `results/`.
  이 세트만으로 STEP4(ROC)·5(GSEA)·6(MR) 재현 가능.

## 라이선스 / 출처
- 저자 원본 코드: GitHub `ljw71865/R-scripts-for-pipeline-reproducibility`,
  결과·보충자료: figshare DOI 10.6084/m9.figshare.30190087 (CC BY 4.0).
- MSigDB gene set(gmt)은 MSigDB 라이선스를 따름(재배포 시 원 출처 확인 필요).
