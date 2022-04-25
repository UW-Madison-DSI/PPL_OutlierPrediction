DetrendedCov.csv- reduced column version of the WasteData_02-09-22_DeTrended.csv with columns: Date, N1, LoessN1, DetrendedN1

WasteData_02-09-22_DeTrended.csv:
Date:	date of when the sample was collected
FlowRate:	amount of flow in the waste water system when it was measured 
BCoV:	percentage of spiked bovine coronavirus recovered
N1:	the measurement of Covid N1 gene concentration
N2:	another measurement of Covid N2 gene concentration
N1Error: Sample error of the N1 measurement when it was observed.
N2Error: Sample error of the N2 measurement when it was observed.
PMMoV: Measurement of Pepper mild mottle virus gene concentration 
LoessN1: smooth version of N1. meant to capture the trend of the data
LoessN2: smooth version of N2. meant to capture the trend of the data
DeTrendedN1: N1-LoessN1 to capture the variance and the noise in the data	
DeTrendedN2: N2-LoessN2 to capture the variance and the noise in the data
