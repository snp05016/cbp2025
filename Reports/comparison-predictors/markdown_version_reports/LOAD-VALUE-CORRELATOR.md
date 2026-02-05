# LVCP

## details

- the leading philosophy is that TAGE struggles to perdict h2p branches because they are "data dependent". tage relies on branch history and using a sequence of taken or not taken branches. they struggle with data dependent branches because the processor might take the same path based on its historyt but the decision depends on the data value and not the path.
- lvcp solves this by using 3 main ds
  - load marking table -> which is a small cache to mark loads so the predictor knows to track them.
  - load tracking queue -> it stored the prediction inddex by the branch pc, the load pc and the load value.
  - when a branch is encountered it grabs the most recent load value from the tracking queue and combines that vlaue with the branc pc to look up in the correlation table.
  - it finds the highest confidence match, and uses that prediction instead of the tage prediction.
- i mean the essential question is how do you choose the canditate, you choose using the recent loads and their values or in other words their relevance. if a load happened a along time ago, it is moved into the distant load buffer. i guess another important question, is how does it even tie a load and a branch together. the answer to that is that, loads happen before branches and when the branch arrives. it looks back into the queue to find the msot recent load, there is a hash also as i mentioned above. so it has the branch pc ( which branch is asking to be predicted???), load pc ( where did the data come from? ) and load value ( what is the data?? ), if the same combination happened before then it would essentially mean that you would get if the branch was take or not.

## load value correlator vs baseline

![alt text](register_value_aware_predictor_vs_baseline_plots.png)

Better performances in how many frameworks of each type :
--------------------------------------------------------------------------------
Framework | Improve/Total | Success Ratio | Total MPKI Improvement               
Web       |  23/26        |   88.5        |   -10.4414                           
Int       |  20/37        |   54.1        |   -14.4377                           
Fp        |  5/14         |   35.7        |   -2.9306                            
Infra     |  9/16         |   56.3        |   -2.1943                            
Compress  |  1/8          |   12.5        |   -0.1299                            
Media     |  0/4          |   0.0         |   0.0000                             

- so as usual it performs better than the baseline tage predictor. suprisingly the craziest increase is in the web benchmarks. what i think happens is that the way web benchmarks are structured they often take and have a lot of input data, and that would require a lot of access to memory since they require checking a lot of variables for example

```js
if (request.type == GET) { ... }
```

- this would mean that usually how these javascript intensive workloads are is that they would need memory and this would mean that there would be a lot of loads and stores, and the load value correlator is exceptional with loads before branches. here in the exampel there would be a load before the brtanch and this mean the load value correlator would shine.

- the biggest overall improvement would be for the integer frameworks and that also make sense since real time compile benchmark traces would load from memory and then peform the branches.
