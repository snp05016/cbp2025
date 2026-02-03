# tage-sc 2025 andrez seznec report

- by digging around, and looking at the tage implementation of the same tage-sc predictor of 2016 i found that it had some differences in the MPKI as mentioned in the abstract of the paper. below are the differences,
  - MPKI Difference is that the 2025 TAGE-SC shows a roughly 15.6% lower MPKI (3.363 vs 3.986), a drop of 0.623 MPKI, though direct comparison is approximate due to different trace sets (CBP-5 vs  CBP2025).
  - 2016 Optimizations: Bank-interleaving for TAGE tables, partial associativity on medium histories (2-3% MPKI gain), enhanced neural SC with IMLI counters, global backward     history, and multiple local histories (total SC benefit ~8%).
  - Both use TAGE core (geometric histories, tagged tables), SC for bias correction, but the 2025 predictor focused on more SC tables and new history forms like region/target IMLI for hard-to-predict branches.
​
- the tage 2016 predictor was not designed to be affective in hardware, but only to win a completition, due to the unreasonable number of sc tables in the 2016 implementation. later a realistic tage predictor was presented and as citied 
`"Realistic" meaning that the author estimated that it could be implemented
for an aggressive instruction front-end predicting an instruction
block with up to 4 branches (at most one taken) per cycle."`

```
The CBP2025 TAGE-SC is derived from CBP2016 TAGE-SC-L, and replicates most of the features that would prevent any reasonable direct hardware implementation: huge number of distinct tables, complete table interleaving in TAGE, use of local histories, unrealistic prediction latency, .. It features the new optimizations on allocation/replacement policy on TAGE-SC proposed in [11] as well as the optimizations on the IMLI components in SC
```

- the optimization features of the TAGE 2025 that aren't there in the 2016
  - on a misprediction, a lot of the entries from different tables are allocated at the first time, and by setting the U counter of the first entry, it is protected against replacement. moreover the ucounter is also set directly to 2 which supports faster eviction.  
  - it uses probablistic counters that are determined by the confidence in the prediction ( provided by the longest matching counter ), to filter the allocation of entries. 
  - uses 2 way skewed associativity.
- the structural correlator has also been "improved" or so they say. most of these optimization are done from the article/book or whatever about tage in 2024  [[André Seznec. 2024. TAGE: an engineering cookbook. Technical Report 9561. Inria. 1–73 pages. https://hal.science/hal-04804900]]. something that i found pretty interesting was that they use a tagged IMLI, ie a tagged inner most loop iterator to solve the issue of the fall through mis prediction in the last loop iteration of a loop. its tagged because considering two loops, you would notice that they counters would overwrite each other if it was just a single IMLI, so to solve that the have a tagged IMLI. **something that was interesting to me that they use two IMLI tables, and both of them have differnt purposes** the other that was is the branch context IMLI, it is useful for branches that dont have fixed loop iterations, so they use histories based on the previous anchor branch. anchor branches do cause biases, but due to the corelator in an BrIMLI, it was make a safe biased prediction, its still counting based ratherr than pattern based but its kinda more biased now. 
- the branch predictor still uses 2 global history based components from the GEHL ( geo history length predictor ), 
  - xor pc and ghr
  - xor pc with ( longestmathcprediction and globla history )
- funny thing, they got rid of the loop predictor it did not seem to help a lot apparently. marignal gains and it took space also so they nuked it lmao.
## tage sc 2025 vs tage sc-l 2016
- 