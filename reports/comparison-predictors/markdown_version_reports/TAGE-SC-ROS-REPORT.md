# outline of the predictor

- okay i feel like this predictor is not that important to analyse since i already anaylsed a the one by andrez seznec. although, it does have other optimizations that are added on top of tage just like the andrez seznec one. but the optimizations are different.
- first and foremost, it the sequence starts as a quadratic sequence and then goes into a generalized geometric sequence with increasing multipliers. they do since quadratic sequences increase pretty much faster than the the geo sequence in the start and then slow down. so they though why not use qudratic first and then geometric. this also simplifies hardware apparently, since it allows them to use a directly mapped table than a set associative one. this reduces the mpki.
- the improved the "confidecnce" rates so as to use, the statistical correlator over tage. the loop predictor is also treated as the most confident compenent and if the decision is final and cannot be overriden. i guess it makes sense since the 2016 predictor allowed for a "fight" between the sc and loop predictor. moreover the loop predictor is only chosen if the confidence level is >= 2.

## predictor vs the baseline predictor

-![alt text](tage_sc_l_alberto_vs_baseline_plots.png)

- observing the graphs i noticed that there actually isnt much difference between this predictor and andrez seznec.
- nevertheless what i found so interesting is that it essentially has the same improvements as the andrez seznec predictor lowkey nothing changes in terms of the improvements. lowkey if you look athte table below ull find that its almost the same thing if not the worst.

Better performances by framework type (Side-by-Side Comparison):
--------------------------------------------------------------------------------
Framework |       Improve/Total     | Success % | Improve/Total (Table 2) | Success %
                (tage - alberto ros)               (tage - andrez seznec) 
--------------------------------------------------------------------------------
Web       |  25/26                  |   96.2 %  |  25/26                  |   96.2 %
Fp        |  8/14                   |   57.1 %  |  9/14                   |   64.3 %
Int       |  20/37                  |   54.1 %  |  25/37                  |   67.6 %
Infra     |  4/16                   |   25.0 %  |  7/16                   |   43.8 %
Compress  |  1/8                    |   12.5 %  |  1/8                    |   12.5 %
Media     |  3/4                    |   75.0 %  |  3/4                    |   75.0 %
--------------------------------------------------------------------------------
Average   |  61/109                 |   55.96 % |  70/105                 |   66.7 %



- again, even if you compare the, overall mpki differences ull literally find little to no difference here look at the table below,

================================================================================
COMPARISON: TAGE-SC-L-Alberto-Ros vs TAGE-SCL-Andrez-Seznec
================================================================================

Summary Statistics:
--------------------------------------------------------------------------------
Metric                         TAGE-SC-L-Alberto-Ros TAGE-SCL-Andrez-Seznec Difference     
--------------------------------------------------------------------------------
MPKI                           3.3898               3.3918               ↑ 0.06% worse
MR                             2.6697               2.6687               ↓ 0.04% better
IPC                            3.3537               3.5081               ↑ 4.60% better
Cycles                         14106643.7238        14447201.9524        ↑ 2.41% worse
BrPerCyc                       0.4012               0.4257               ↑ 6.10% better
MispBrPerCyc                   0.0089               0.0092               ↑ 2.42% better
CycWPPKI                       151.4364             151.7579             ↑ 0.21% worse


Trace-by-Trace Comparison (MPKI):
--------------------------------------------------------------------------------
TAGE-SC-L-Alberto-Ros wins: 27 traces
TAGE-SCL-Andrez-Seznec wins: 78 traces
Ties: 0 traces


Biggest Improvements (TAGE-SCL-Andrez-Seznec better):
--------------------------------------------------------------------------------
  web_19_trace                               5.0555 →   0.0000 (↓ 5.0555)
  int_1_trace                               11.4406 →  10.5130 (↓ 0.9276)
  int_2_trace                               10.7227 →   9.9876 (↓ 0.7351)
  int_22_trace                               4.8095 →   4.2827 (↓ 0.5268)
  int_21_trace                              23.7948 →  23.3263 (↓ 0.4685)


Biggest Regressions (TAGE-SCL-Andrez-Seznec worse):
--------------------------------------------------------------------------------
  int_9_trace                                0.0000 →   3.6628 (↑ 3.6628)
  web_3_trace                                0.0000 →   3.5523 (↑ 3.5523)
  web_14_trace                               0.0000 →   3.3625 (↑ 3.3625)
  fp_6_trace                                 0.0000 →   0.8550 (↑ 0.8550)
  int_16_trace                               5.5268 →   5.8279 (↑ 0.3011)
