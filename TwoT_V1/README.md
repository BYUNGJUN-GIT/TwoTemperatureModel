# Two-temperature model (TTM) for Al/SiO2/Si

이 폴더는 기존 `ThreeT_*` 코드 흐름과 유사하게 사용할 수 있는 2온도 모델 세트입니다.

## 파일 구성
- `TwoT_Para_v1_AlSiO2Si.m`: 샘플 geometry/material/pump 조건 정의 + heat absorption + pulse profile + raw data import
- `TwoT_CeLeTdep_Pulse_v1.m`: 1D TTM solver (electron-lattice coupled diffusion)
- `TwoT_MAIN_v1_AlSiO2Si.m`: 실행 스크립트 (계산, raw data 비교 플롯, CSV 저장, sensitivity 옵션)
- `TwoT_CeLeTdep_Pulse_SENS_v1.m`: TE(전자 온도)의 thermal parameter sensitivity 계산

## 사용 방법
1. MATLAB에서 현재 작업 경로를 repository root로 설정
2. `new_folder/TwoT_MAIN_v1_AlSiO2Si.m` 실행
3. 필요 시 `TwoT_MAIN_v1_AlSiO2Si.m`에서 `sens_th = 1`로 설정해 sensitivity 계산
4. 결과 파일 확인
   - `new_folder/Para_AlSiO2Si_default.mat`
   - `new_folder/T_AlSiO2Si_default.csv`

## 옵션 파일 (선택)
- `AbsCal_AlSiO2Si.mat` (변수명: `AbsCal`) : 흡수 프로파일 직접 로드
- `PulseWidth_AlSiO2Si.mat` (변수명: `PulseWidth`) : 펌프-프로브 cross-correlation 로드
- `new_folder/AlSiO2Si_default.mat` (변수명: `RawData`) : 측정 데이터 로드 및 모델 비교

## 출력
- `T_E`: 전자 온도
- `T_L`: 격자(포논) 온도
- (raw data 존재 시) 정규화 데이터와 모델 응답을 함께 저장

> 참고: TTM에서는 보통 격자와 포논을 동일한 thermal bath로 취급합니다.


## 계면 컨덕턴스 입력 형식
- `TwoT_Para_v1_AlSiO2Si.m`의 `G`는 interface당 한 줄(`;`)로 정의합니다.
- 각 줄은 2x2 계면 컨덕턴스 행렬을 펼친 `[G11, G12, G21, G22]` 순서입니다.
- 현재 solver는 `G11`(격자/포논), `G22`(전자) 대각 성분을 사용합니다.
