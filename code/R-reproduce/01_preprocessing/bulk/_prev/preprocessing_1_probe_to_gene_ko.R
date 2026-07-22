####==================================================================####
####   data preprocessing 1.R  (한글 해석 버전)
####   목적: GEO 마이크로어레이 원본(probe 단위) → 유전자 단위 발현행렬 생성
####   출력: geneMatrix3.txt  (다음 단계 preprocessing 2.R 의 입력)
####==================================================================####


####----기본 설정 (Basic operations)----####
rm(list = ls())          # 메모리에 있는 모든 변수 삭제 → 깨끗한 상태에서 시작
gc()                     # garbage collection: 안 쓰는 메모리 회수
getwd()                  # 현재 작업 폴더 경로 출력 (확인용)
stringsASFactors = FALSE  # ※ 오타. 원래 의도는 options(stringsAsFactors = FALSE)
                          #   문자열을 factor로 자동변환 막으려는 것인데
                          #   이렇게 쓰면 그냥 변수 하나 만들 뿐 효과 없음 (무시해도 됨)

## R 패키지 로드 (필요한 기능들을 불러옴)
library(readr)      # read_tsv 등 빠른 파일 읽기
library(dplyr)      # select, mutate, inner_join 등 데이터 가공
library(tidyverse)  # dplyr/ggplot 등 묶음 (위와 일부 중복)
library(ggplot2)    # 그래프용 (이 파일에선 실제로 안 씀)
library(ggpubr)     # 그래프용 (이 파일에선 실제로 안 씀)
library(limma)      # avereps() 함수 쓰려고 로드 (중복 probe 평균)
library(stringr)    # str_split 등 문자열 처리

setwd("G:\\187geneMR\\03.download\\GSE37263")
# 작업 폴더 지정.
# ★★★ 이 경로는 강의 템플릿 흔적. 네 실제 데이터 폴더 경로로 반드시 수정.


####----발현행렬(expression matrix) 읽기----####
## 표현 행렬 읽기
Sys.setenv("VROOM_CONNECTION_SIZE" = 262144)
# 파일 읽기 버퍼 크기를 2배(약 256KB)로 늘림.
# GEO series_matrix 헤더 줄이 매우 길어서, 기본값이면
# "line too long / connection buffer" 에러가 나는 걸 방지.

exprSet <- read_tsv("GSE37263_series_matrix.txt.gz", skip = 67)
# GEO series_matrix 파일 읽기.
# 앞 67줄은 !Series_... 메타정보라서 skip으로 건너뜀.
# ★ skip 숫자와 파일명은 데이터셋마다 다름 → 실제 표가 시작되는 줄에 맞게 수정.

exprSet <- data.frame(exprSet)   # tibble → 일반 data.frame 형식으로 변환
class(exprSet)                   # 타입이 data.frame 인지 확인 (점검용)

# 행=probe ID, 열=샘플 인 발현행렬 만들기
rownames(exprSet) <- exprSet$ID_REF  # 첫 컬럼(probe ID)을 '행 이름'으로 지정
exprSet <- exprSet[,-1]              # 이제 중복되는 첫 컬럼(ID_REF) 삭제
                                     # → 결과: 행=probe, 열=샘플 인 순수 발현행렬


####----probe ID를 유전자 이름으로 변환할 사전(GPL) 준비----####
library(readr)
library(dplyr)
probe <- read_tsv("GPL5175-3188.txt", skip = 12)
# 플랫폼(GPL) annotation 파일 = "probe ID가 어떤 유전자인지" 알려주는 사전.
# 앞 12줄은 주석이라 skip. ★ 파일명·skip값은 데이터셋(플랫폼)마다 다름.

Probe_ID <- select(probe, c("ID","Gene Symbol"))
# 사전에서 필요한 두 컬럼만 추출: probe ID + 유전자 심볼(Gene Symbol)

ids <- Probe_ID %>%
  mutate(`Gene Symbol` = str_split(`Gene Symbol`, " // ", simplify = TRUE)[, 2])
# 한 probe에 유전자가 여러 개 붙은 경우("AAA // BBB // CCC" 형식) 처리.
# " // " 기준으로 쪼갠 뒤 [ , 2] = 2번째 조각을 대표 유전자로 사용.
# ★ 주의: GPL 형식에 따라 1번째([,1])를 써야 할 수도 있음. 파일 보고 위치 확인.


####----사전을 적용해 probe → 유전자 발현행렬로 변환----####
library(dplyr)
library(tidyverse)

exprSet1 <- exprSet %>%
  rownames_to_column("ID")
# 행 이름(probe ID)을 다시 "ID"라는 실제 컬럼으로 꺼냄 (join 하기 위해)

ids$ID <- as.character(ids$ID)          # join 위해 ID 타입을 문자(character)로 통일
exprset1 <- inner_join(ids, exprSet1, by="ID")
# probe ID를 기준으로 [유전자 이름] + [발현값]을 결합.
# inner_join → 사전과 발현행렬 양쪽에 다 존재하는 probe만 남김.

length(exprset1$ID)              # 남은 probe 개수 확인
length(exprset1$'Gene Symbol')   # 유전자 심볼 개수 확인 (점검용)

exprset1 = avereps(exprset1[,-c(1,2)],       # 1,2열(ID·유전자심볼) 제외한 발현값만 사용
                   ID = exprset1$'Gene Symbol')
# ★핵심 단계: 같은 유전자를 가리키는 여러 probe의 값을 '평균'내서
#   "유전자 1개 = 행 1개"로 정리 (probe 중복 제거). avereps = average replicates.

exprset1 <- data.frame(exprset1)
# 결과를 data.frame으로: 행=유전자 이름, 열=샘플 인 최종 발현행렬


####----결과 저장 (파일로 내보내기)----####
# exprSet 이 데이터프레임이고 행 이름이 유전자명(geneNames)인지 확인
geneMatrix <- exprset1
geneMatrix <- cbind(geneNames = rownames(geneMatrix), geneMatrix)
# 행 이름(유전자명)을 "geneNames"라는 첫 컬럼으로 추가.
# → 텍스트로 저장할 때 유전자 이름이 파일 안에 컬럼으로 들어가게 하려는 것.

# TXT 파일로 내보내기 (탭 구분)
write.table(
  geneMatrix,
  file = "geneMatrix3.txt",  # 저장 파일명
  sep = "\t",                # 탭(TAB)으로 컬럼 구분
  quote = FALSE,             # 값에 따옴표 붙이지 않음
  row.names = FALSE,         # 행 이름은 따로 안 씀 (이미 geneNames 컬럼에 있음)
  col.names = TRUE           # 헤더(샘플명) 기록
)
# 최종 결과물: geneMatrix3.txt
# → 다음 스크립트 preprocessing 2.R 의 입력 파일이 됨.


####==================================================================####
####  ★ 데이터셋마다 반드시 바꿔야 하는 3곳:
####    1) setwd(...)              → 실제 데이터 폴더 경로
####    2) skip = 67 / skip = 12   → 실제 표/주석 시작 줄에 맞게
####    3) str_split(...)[, 2]      → 대표 유전자 위치([,1] 또는 [,2])
####==================================================================####
