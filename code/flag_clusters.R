
flag.clusters <- function(seurat.obj, genes, resolution, quantile) {

  clusters <- as.character(unique(seurat.obj@meta.data[, resolution]))

  batches <- as.character(unique(seurat.obj$run))

  flagged.clusters <- c()

  for (gene in genes) {

    suspicious.clusters <- c()

    for (batch in batches) {

      if (gene == "percent.mt") {
        cluster.mean <- sapply(clusters, FUN = function(x) {
          mean(seurat.obj$percent.mt[which(seurat.obj@meta.data[, resolution] == x &
                                             seurat.obj$run == batch)])
          })
        }

      else {
        cluster.mean <- sapply(clusters, FUN = function(x) {
          mean(seurat.obj@assays$SCT@data[gene, which(seurat.obj@meta.data[, resolution] == x &
                                                        seurat.obj$run == batch)])
          })
        }
      threshold <- quantile(na.omit(cluster.mean), probs = quantile)

      suspicious.clusters <- c(suspicious.clusters, clusters[which(cluster.mean >= threshold)])

      }

    suspicious.clusters.table <- as.data.frame(table(suspicious.clusters))
    colnames(suspicious.clusters.table) <- c("cluster", "freq")

    flag <- as.character(suspicious.clusters.table$cluster[which(suspicious.clusters.table$freq == length(batches))])

    print(paste0("Cluster ", flag, " scored above threshold for ", gene))

    flagged.clusters <- c(flagged.clusters, flag)

  }

  return(unique(flagged.clusters))

}



flag.clusters2 <- function(seurat.obj, resolution, fraction) {

  clusters <- as.character(unique(seurat.obj@meta.data[, resolution]))

  flagged.clusters <- c()

  cluster.table <- as.data.frame(table(seurat.obj@meta.data[, resolution], seurat.obj$run))
  colnames(cluster.table) <- c("cluster", "run", "frequency")

  for (cluster in clusters) {


    idx <- which(cluster.table$cluster == cluster)

    cluster.fraction <- cluster.table$frequency[idx] / sum(cluster.table$frequency[idx])

    if (max(cluster.fraction) > fraction) {

      flag <- cluster

      print(paste0("Cluster ", flag, " consisted of ", round(max(cluster.fraction), digits = 2) * 100,
                   "% cells from one run"))

      flagged.clusters <- c(flagged.clusters, flag)


    }

  }

  return(unique(flagged.clusters))

}

