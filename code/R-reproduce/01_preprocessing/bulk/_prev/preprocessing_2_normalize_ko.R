####==================================================================####
####   data preprocessing 2.R  (한글 해석 버전)
####   목적: 발현행렬을 log2 자동판정 후 정규화(normalizeBetweenArrays)하고,
####         Control/Early/Late 세 그룹 라벨을 붙여 파일로 저장한다.
####   입력: GSE142025_matrix.txt, s1.txt(Control), s2.txt(Early), s3.txt(Late)
####   출력: GSE142025_threeGroups.normalize.txt
####==================================================================####

#if (!requireNamespace("BiocManager", quietly = TRUE))          # BiocManager(생물정보 패키지 설치 도구)가 없으면
#    install.packages("BiocManager")                            # 설치한다 (필요 시 주석 해제)
#BiocManager::install("limma")                                  # limma 패키지 설치 (필요 시 주석 해제)


# 패키지 불러오기
library(limma)                                                  # limma: 발현행렬 정규화/분석용 패키지 (avereps, normalizeBetweenArrays 사용)

# 파일 경로 설정
inputFile <- "GSE142025_matrix.txt"   # 발현 데이터 파일              # ★ 데이터셋마다 바꿔야 하는 발현행렬 파일명
controlFile <- "s1.txt"               # 대조군(Control) 샘플 정보     # ★ 대조군 샘플 목록 파일명
earlyFile <- "s2.txt"                 # 초기(Early) 그룹 샘플 정보    # ★ 초기 그룹 샘플 목록 파일명
lateFile <- "s3.txt"                  # 후기(Late) 그룹 샘플 정보     # ★ 후기 그룹 샘플 목록 파일명
geoID <- "GSE142025"                                            # ★ GEO 데이터셋 ID (출력 파일명 앞에 붙음)
setwd("G:\\187geneMR\\05.normalize")                            # ★ 작업 디렉터리(입출력 파일이 있는 폴더) 지정

# 발현행렬 읽기 및 전처리
rt <- read.table(inputFile, header=TRUE, sep="\t", check.names=FALSE)  # 탭 구분 텍스트를 데이터프레임으로 읽기(열이름 원본 유지)
rt <- as.matrix(rt)                                             # 행렬 형태로 변환
rownames(rt) <- rt[,1]                                          # 첫 번째 열(유전자명)을 행 이름으로 지정
exp <- rt[,2:ncol(rt)]                                          # 유전자명 열을 제외한 나머지(발현값)만 추출
dimnames <- list(rownames(exp), colnames(exp))                 # 행 이름(유전자)과 열 이름(샘플)을 보관
data <- matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)  # 문자를 숫자형 행렬로 변환(이름 유지)
rt <- avereps(data)                                             # 같은 유전자명이 중복되면 평균으로 합쳐 하나로 만듦

# log2 변환 여부 자동 판단
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  # 전체 값의 분위수(0,25,50,75,99,100%)를 계산
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )         # 값의 크기/범위로 아직 log2 안 됐는지 판정(참이면 미변환)
if(LogC){                                                       # log2 변환이 필요하다고 판단되면
  rt[rt<0]=0                                                    # 음수 값은 0으로 처리
  rt=log2(rt+1)}                                                # log2(x+1) 변환 수행 (0 방지 위해 +1)
data=normalizeBetweenArrays(rt)                                 # 샘플(어레이) 간 분포를 맞춰 정규화

                 # 세 그룹 샘플 정보 읽기
                 sample_control <- read.table(controlFile, header=FALSE, sep="\t", check.names=FALSE)  # 대조군 샘플 목록 읽기(헤더 없음)
                 sample_early <- read.table(earlyFile, header=FALSE, sep="\t", check.names=FALSE)      # 초기 그룹 샘플 목록 읽기
                 sample_late <- read.table(lateFile, header=FALSE, sep="\t", check.names=FALSE)        # 후기 그룹 샘플 목록 읽기

                 # 샘플 이름 정리
                 sampleName_control <- gsub("^ | $", "", as.vector(sample_control[,1]))  # 이름 앞뒤 공백 제거(대조군)
                 sampleName_early <- gsub("^ | $", "", as.vector(sample_early[,1]))      # 이름 앞뒤 공백 제거(초기)
                 sampleName_late <- gsub("^ | $", "", as.vector(sample_late[,1]))        # 이름 앞뒤 공백 제거(후기)

                 # 세 그룹 데이터 추출
                 controlData <- data[, sampleName_control]      # 정규화 행렬에서 대조군 샘플 열만 뽑기
                 earlyData <- data[, sampleName_early]          # 초기 그룹 샘플 열만 뽑기
                 lateData <- data[, sampleName_late]            # 후기 그룹 샘플 열만 뽑기

                 # 세 그룹 데이터 합치기
                 combinedData <- cbind(controlData, earlyData, lateData)  # 세 그룹을 열 방향으로 이어붙이기
                 controlNum <- ncol(controlData)                # 대조군 샘플 개수
                 earlyNum <- ncol(earlyData)                    # 초기 그룹 샘플 개수
                 lateNum <- ncol(lateData)                      # 후기 그룹 샘플 개수

 # 세 그룹 라벨을 붙여 정규화된 전체 유전자 발현량 출력
Type <- c(rep("Control", controlNum), rep("Early", earlyNum), rep("Late", lateNum))  # 각 샘플의 그룹 라벨 벡터 생성
outData <- rbind(id=paste0(colnames(combinedData), "_", Type), combinedData)  # 첫 행에 "샘플명_그룹" 라벨 행 추가
write.table(outData,file=paste0(geoID, "_threeGroups.normalize.txt"), sep="\t",quote=F,col.names=F)  # 탭 구분 파일로 저장(따옴표X, 열이름X)

