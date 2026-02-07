# performance delta analysis

- Delta Statistics (Δ = MPKI_LVCP - MPKI_RVA-Toru):
  Average:  +0.1664 MPKI
  Std Dev:  0.4111 MPKI
  Median:   +0.0722 MPKI

Improvement Analysis:
  RVA-Toru beats LVCP: 92/105 traces (87.6%)

Worst-Case MPKI:
  LVCP: 23.5688
  RVA-Toru: 21.6842
![alt text](delta_scatter_lvcp_vs_rvatoru.png)
- i can conclude that RVA Toru mostly subsumes LVCP, and rva toru provides superior prediction on 87.6% of workloads. this is a good indictor. and it makes sense why its better. rva-toru uses registers to compute what branches can be predicted, rva hashes multples values and also deals with register renaming and reuse, so it kinda makes sense why it does better. 

- lvcp might perform better in cases, where there are sufficient loads from memory before branches, especially conditional ones where values are compared. and lvcp only tracks memory driven values where as rva-toru takes all register values and their idea of using digests which are 12 bit encodings is pretty neat as well. 
