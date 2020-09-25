# CIDER: meta-Clustering based on Inter-group Differential Expression

Introduction xxxx

## Installation


* A [potential issue](https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/) about gfortran: `ld: warning: directory not found for option '-L/usr/local/lib/gcc/i686-apple-darwin8/4.2.3/x86_64'`

* gfortran can be installed from [here](https://cran.r-project.org/bin/macosx/tools/)

## Tutorial - clustering



### Step 1: Initial clustering (optional)

```
dist_list <- getDistMat(seu_list)
plotDistMat(dist_list)
```

```
seu_list <- mergeInitialClusters(seu_list, dist_list)
```

### Step 2: Compute inter-group similarity matrix


```
seu$initial_cluster <- paste(seu$Group, seu$Batch, sep = "_")
dist <- getIDEr(seu)
```


Visualisation:

Network Graph

```
net <- plotNetwork(seu, dist)
```




### Step 3: Final clustering



```
seu <- finalClustering(seu, dist)
```


## Compatability





