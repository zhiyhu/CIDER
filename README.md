# CIDEr: meta-Clustering based on Inter-group Differential Expression

Introduction xxxx

## Installation



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





