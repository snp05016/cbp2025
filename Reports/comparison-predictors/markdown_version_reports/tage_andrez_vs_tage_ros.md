# tage_andrez vs tage_ros  
Data Quality:
  Valid traces:   100
  Skipped traces: 5 (zero MPKI)

Delta Statistics (Δ = MPKI_TAGE-SC-L - MPKI_TAGE-SCL):
  Average:  +0.0617 MPKI
  Std Dev:  0.1539 MPKI
  Median:   +0.0233 MPKI

Improvement Analysis:
  TAGE-SCL beats TAGE-SC-L: 77/100 traces (77.0%)

Worst-Case MPKI:
  TAGE-SC-L: 23.7948
  TAGE-SCL: 23.3263

i can conclude that there isnt much of a difference of mpki prolly like 0.06 which is very less to even conclude somrthing. they are literally based on the baseline and both provide optimizations and even though they are diffferent optimizations it looks likes they arent doing much overall. if we are worried about the smallest of the smallest improvements then definitely tage andrez wins and it wins on 77% of the benchmarks, so ideally tage andrez very weakly subsumes tage alberto ros. 