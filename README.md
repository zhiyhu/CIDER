# CIDER


## Quick Start - asCIDER

**Step 1**: Prepare batch-specific group (initial cluster) annotations

```
seu$initial_cluster <- paste(seu$Group, seu$Batch, sep = "_")
```

See here for more information to de novo create your initial clusters.

```
dist_list <- getDistMat(seu_list)
plotDistMat(dist_list)
seu_list <- mergeInitialClusters(seu_list, dist_list)
```

**Step 2**: Compute inter-group similarity matrix

```
dist <- getIDEr(seu)
```

**Step 3**: Final clustering

```
seu <- finalClustering(seu, dist)
```

## Visualisation options

Network Graph
```
net <- plotNetwork(seu, dist)
```

## Quick Start - evaluation




## Compatability


## Technical notes

This package was developed using R version 3.6.3 (2020-02-29), platform: x86_64-apple-darwin15.6.0 (64-bit).



