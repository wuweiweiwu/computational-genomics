# HW3 Hung-Wei Wu, wuxx1045

I am using python 2.7

```
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
```

### How to run

```
./hw3 <path to sequence file>

ex:

./hw3 hw3.fna
```

### Visualize

```
Rscript hw3-plot-newick.r tree.txt hw3-tip-labels.txt
or
Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt
```


### Bonus: bootstrapping

uncomment code in main() because it takes forever to run

line 338 - 339

```
percentages = bootstrap(root, ids, id_sequences)
write_bootstrap(percentages)
```

Visualize

```
Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt bootstrap.txt
```


### Files


```
edges.txt
tree.txt
bootstrap.txt

newick-tree.pdf
edges-tree.pdf
bootstrap-tree.pdf
```
