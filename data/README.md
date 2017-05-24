Data
====


1. data for figure 3.


2. Tree for RVARI TPS
```
~/bin/standard-RAxML-master/raxmlHPC-PTHREADS  -T 30 -p 12345 -m PROTGAMMAGTR -n tsp-small-best -s tsp.blast.picked.final.mafft
~/bin/standard-RAxML-master/raxmlHPC-PTHREADS -f b -b 12345 -# 100 -T 30 -p 12345 -m PROTGAMMAGTR -n tsp-small -s tsp.blast.picked.final.mafft
~/bin/standard-RAxML-master/raxmlHPC-PTHREADS -m PROTGAMMAGTR -p 12345 -f b -t RAxML_bestTree.tsp-small-best -z RAxML_bootstrap.tsp-small -n tsp-small-merged

```
