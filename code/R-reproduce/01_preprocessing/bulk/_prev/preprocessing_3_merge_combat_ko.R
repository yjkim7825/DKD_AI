####==================================================================####
####   data preprocessing 3.R  (한글 해석 버전)
####   목적: 여러 GEO 데이터셋을 공통 유전자(intersect)로 병합한 뒤,
####         ComBat으로 배치효과(데이터셋/실험실 차이)를 제거한다.
####   입력: 작업 폴더 내 모든 *.txt 발현 파일 (merge.* 두 개는 제외)
####   출력: merge.preNorm.txt(병합), merge.normalize.txt(배치보정 후)
####==================================================================####

######영상 출처: https://ke.biowolf.cn
######생물정보 자습 사이트: https://www.biowolf.cn/
######위챗 공식계정: biowolf_cn
######협력 이메일: biowolf@foxmail.com
######문의 위챗: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))          # BiocManager가 없으면
#    install.packages("BiocManager")                            # 설치한다 (필요 시 주석 해제)
#BiocManager::install("limma")                                  # limma 패키지 설치 (필요 시 주석 해제)

#if (!requireNamespace("BiocManager", quietly = TRUE))          # BiocManager가 없으면
#    install.packages("BiocManager")                            # 설치한다 (필요 시 주석 해제)
#BiocManager::install("sva")                                    # sva 패키지 설치 (ComBat 함수 포함, 필요 시 주석 해제)


#패키지 불러오기
library(limma)                                                  # limma: avereps 등 발현행렬 처리용
library(sva)                                                    # sva: ComBat(배치효과 제거) 함수 제공
setwd("G:\\187geneMR\\06.combat")      #작업 디렉터리 설정      # ★ 병합할 txt 파일들이 모여 있는 폴더 경로

#작업 폴더 내 ".txt"로 끝나는 모든 파일 목록 얻기
files=dir()                                                     # 현재 폴더의 모든 파일/폴더 이름 나열
files=grep("txt$", files, value=T)                             # 그중 이름이 "txt"로 끝나는 것만 선택
geneList=list()                                                # 각 데이터셋의 유전자 목록을 담을 리스트 초기화

#모든 txt 파일에서 유전자 정보를 읽어 geneList에 저장
for(file in files){                                            # 각 txt 파일에 대해 반복
	if(file=="merge.preNorm.txt"){next}                        # 병합 결과 파일이면 건너뛰기
	if(file=="merge.normalize.txt"){next}                      # 배치보정 결과 파일이면 건너뛰기
    rt=read.table(file, header=T, sep="\t", check.names=F)      #입력 파일 읽기(탭 구분, 열이름 원본 유지)
    geneNames=as.vector(rt[,1])      #유전자 이름 추출          # 첫 열(유전자명)을 벡터로 추출
    uniqGene=unique(geneNames)       #유전자 중복 제거          # 중복 유전자명을 하나로 정리
    header=unlist(strsplit(file, "\\.|\\-"))                    # 파일명을 "."과 "-" 기준으로 쪼갬(데이터셋 이름 추출용)
    geneList[[header[1]]]=uniqGene                              # 데이터셋 이름을 key로 유전자 목록 저장
}

#교집합 유전자 구하기
interGenes=Reduce(intersect, geneList)                         # 모든 데이터셋에 공통으로 존재하는 유전자만 남김

#데이터 병합
allTab=data.frame()                                            # 병합 결과를 담을 빈 데이터프레임
batchType=c()                                                  # 각 샘플이 어느 데이터셋(배치)인지 표시할 벡터
for(i in 1:length(files)){                                     # 각 파일에 대해 순번 i로 반복
    inputFile=files[i]                                         # i번째 파일명
	if(file=="merge.preNorm.txt"){next}                        # 병합 결과 파일이면 건너뛰기(주: 원본 코드는 file 변수 사용)
	if(file=="merge.normalize.txt"){next}                      # 배치보정 결과 파일이면 건너뛰기
    header=unlist(strsplit(inputFile, "\\.|\\-"))              # 파일명에서 데이터셋 이름 추출
    #입력 파일을 읽고 정리
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)  # 발현 파일 읽기
    rt=as.matrix(rt)                                          # 행렬로 변환
    rownames(rt)=rt[,1]                                       # 첫 열(유전자명)을 행 이름으로 지정
    exp=rt[,2:ncol(rt)]                                       # 유전자명 열 제외한 발현값만 추출
    dimnames=list(rownames(exp),colnames(exp))               # 행/열 이름 보관
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)  # 숫자형 행렬로 변환
    rt=avereps(data)                                         # 중복 유전자명을 평균으로 합침
    colnames(rt)=paste0(header[1], "_", colnames(rt))        # 샘플명 앞에 데이터셋 이름을 붙여 구분(중복 방지)

    #데이터 병합
    if(i==1){                                                 # 첫 번째 파일이면
    	allTab=rt[interGenes,]                                # 공통 유전자만 뽑아 초기 테이블로 설정
    }else{                                                    # 두 번째부터는
    	allTab=cbind(allTab, rt[interGenes,])                # 공통 유전자 행을 열 방향으로 이어붙임
    }
    batchType=c(batchType, rep(i,ncol(rt)))                  # 이 데이터셋 샘플들에 배치 번호 i 부여
}

#병합된 발현 데이터 출력
outTab=rbind(geneNames=colnames(allTab), allTab)             # 첫 행에 샘플 이름 행을 추가
write.table(outTab, file="merge.preNorm.txt", sep="\t", quote=F, col.names=F)  # 배치보정 전 병합 결과 저장

#병합 데이터에 배치 보정 수행 후, 배치 보정된 발현 데이터 출력
outTab=ComBat(allTab, batchType, par.prior=TRUE)            # ComBat으로 데이터셋 간 배치효과 제거(파라메트릭 사전분포 사용)
outTab=rbind(geneNames=colnames(outTab), outTab)           # 첫 행에 샘플 이름 행을 추가
write.table(outTab, file="merge.normalize.txt", sep="\t", quote=F, col.names=F)  # 배치보정 후 결과 저장


