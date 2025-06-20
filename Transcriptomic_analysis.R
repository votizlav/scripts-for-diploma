#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(tximport)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pathview)
  library(ReactomePA)
  library(EnhancedVolcano)
  library(ggplot2)
  library(pheatmap)
  library(reshape2)
  library(enrichplot)
  library(msigdbr)
  library(apeglm)
})

# 1. Чтение путей от пользователя
gff3_path <- readline(prompt = "Введите полный путь к GFF3: ")
if (!file.exists(gff3_path)) stop("GFF3-файл не найден: ", gff3_path)

meta_path <- readline(prompt = "Введите путь к CSV с метаданными (sample,file,condition): ")
if (!file.exists(meta_path)) stop("CSV с метаданными не найден: ", meta_path)

meta <- tryCatch({
  read.csv(meta_path, stringsAsFactors = FALSE)
}, error = function(e) stop("Невозможно прочитать CSV: ", e$message))

if (!all(c("sample","file","condition") %in% colnames(meta))) {
  stop("CSV должен содержать колонки: sample, file, condition")
}

missing_files <- meta$file[!file.exists(meta$file)]
if (length(missing_files) > 0) {
  warning("Отсутствуют файлы:\n", paste(missing_files, collapse = "\n"))
  ans <- readline(prompt = "Продолжить без них? (y/n): ")
  if (tolower(ans) != "y") stop("Прервано из-за отсутствующих файлов.")
}

conds <- unique(meta$condition)
meta$condition <- factor(meta$condition, levels = conds)

output_base <- readline(prompt = "Введите базовую директорию для вывода: ")
if (!dir.exists(output_base)) {
  dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
}

# 2. Создание tx2gene
message("Создаём tx2gene из GFF3...")
txdb <- tryCatch({
  makeTxDbFromGFF(gff3_path)
}, error = function(e) stop("Ошибка makeTxDbFromGFF: ", e$message))
transcripts_df <- as.data.frame(mcols(transcripts(txdb, columns = c("TXNAME","GENEID"))))
transcripts_df$TXNAME <- as.character(transcripts_df$TXNAME)
transcripts_df$GENEID <- as.character(transcripts_df$GENEID)
tx2gene <- transcripts_df[, c("TXNAME","GENEID")]

# 3. Tximport
message("Импортируем quantification...")
files <- setNames(meta$file, meta$sample)
txi <- tryCatch({
  tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
}, error = function(e) stop("tximport failed: ", e$message))
txi$counts <- as.matrix(txi$counts)

# 4. DESeq2-объект
coldata <- data.frame(row.names = colnames(txi$counts),
                      condition = meta$condition)
dds <- tryCatch({
  DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
}, error = function(e) stop("Ошибка создания DESeqDataSet: ", e$message))

# 5. Запуск DESeq
message("Запускаем DESeq...")
dds <- tryCatch({
  DESeq(dds)
}, error = function(e) stop("DESeq failed: ", e$message))

vsd <- tryCatch({
  vst(dds, blind = FALSE)
}, error = function(e) stop("VST failed: ", e$message))

# 6. Подготовка аннотаций
ensembl_ids <- rownames(dds)
gene_symbols_all <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                           column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_desc_all <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                        column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")

# 7. Общая PCA и корреляция
plotPCA_vsd <- function(vsd_obj, coldata, outfile) {
  p <- tryCatch({
    plotPCA(vsd_obj, intgroup = "condition") + theme(aspect.ratio = 1)
  }, error = function(e) {
    message("PCA plotting failed: ", e$message); return(NULL)
  })
  if (!is.null(p)) {
    png(outfile, width = 800, height = 800, res = 200)
    print(p)
    dev.off()
  }
}

plotSampleCorr <- function(vsd_obj, coldata, outfile) {
  mat <- tryCatch({
    cor(assay(vsd_obj))
  }, error = function(e) {
    message("Correlation calc failed: ", e$message); return(NULL)
  })
  if (!is.null(mat)) {
    annotation_col_df <- data.frame(condition = coldata$condition)
    rownames(annotation_col_df) <- rownames(coldata)
    png(outfile, width = 800, height = 700, res = 200)
    tryCatch({
      pheatmap(mat, main="Sample correlation", annotation_col=annotation_col_df,
               color = colorRampPalette(c("blue","white","red"))(100))
    }, error = function(e) message("Heatmap failed: ", e$message))
    dev.off()
  }
}

plotPCA_vsd(vsd, coldata, file.path(output_base, "PCA_all.png"))
plotSampleCorr(vsd, coldata, file.path(output_base, "SampleCorr_all.png"))

# 8. Ввод контрастов
n_contrasts <- as.integer(readline(prompt = "Сколько контрастов анализировать? "))
if (is.na(n_contrasts) || n_contrasts <= 0) stop("Некорректное число контрастов.")

contrast_list <- vector("list", n_contrasts)
for (i in seq_len(n_contrasts)) {
  cat("Контраст ", i, "/", n_contrasts, "\n", sep = "")
  repeat {
    g1 <- readline(prompt = paste0("  Базовая группа (из: ", paste(levels(meta$condition), collapse=","), "): "))
    g2 <- readline(prompt = paste0("  Сравниваемая группа (из: ", paste(levels(meta$condition), collapse=","), "): "))
    if (g1 %in% levels(meta$condition) && g2 %in% levels(meta$condition) && g1 != g2) break
    cat("  Ошибочный ввод, повторите.\n")
  }
  cname <- paste0(g2, "_vs_", g1)
  contrast_list[[i]] <- list(g1=g1, g2=g2, name=cname)
}

# 9. Функция анализа контраста
run_contrast_analysis <- function(dds_obj, vsd_obj, coldata, contrast, gene_symbols_all, gene_desc_all, outdir) {
  g1 <- contrast$g1; g2 <- contrast$g2; cname <- contrast$name
  message("Анализ контраста: ", cname)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Всегда используем apeglm для shrink
  res <- tryCatch({
    coef_name <- paste0("condition_", g2, "_vs_", g1)
    if (coef_name %in% resultsNames(dds_obj)) {
      lfcShrink(dds_obj, coef = coef_name, type = "apeglm")
    } else {
      message("Коэффициент ", coef_name, " не найден, используем results() без shrink.")
      results(dds_obj, contrast = c("condition", g2, g1))
    }
  }, error = function(e) {
    message("lfcShrink/results failed: ", e$message)
    tryCatch(results(dds_obj, contrast = c("condition", g2, g1)),
             error = function(e2) stop("results also failed: ", e2$message))
  })
  
  res_df <- as.data.frame(res)
  res_df$gene_symbol <- gene_symbols_all[rownames(res_df)]
  res_df$gene_description <- gene_desc_all[rownames(res_df)]
  write.csv(res_df, file = file.path(outdir, paste0("DE_full_", cname, ".csv")), row.names = TRUE)
  
  try({
    png(file.path(outdir, paste0("MA_", cname, ".png")), width=1200, height=1000, res=120)
    plotMA(res, ylim=c(-10,10), main=paste0("MA: ", cname))
    dev.off()
  }, silent = TRUE)
  
  try({
    png(file.path(outdir, paste0("Volcano_", cname, ".png")), width=1200, height=1500, res=120)
    EnhancedVolcano(res_df, lab=res_df$gene_symbol,
                    x="log2FoldChange", y="padj",
                    pCutoff=0.05, FCcutoff=1,
                    title=paste0("Volcano: ", cname),
                    pointSize=2, labSize=4)
    dev.off()
  }, silent = TRUE)
  
  sig <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange) &
                  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
  if (nrow(sig) > 0) {
    sig <- sig[!is.na(sig$gene_symbol) & sig$gene_symbol != "", ]
    write.csv(sig, file = file.path(outdir, paste0("DE_sig_", cname, ".csv")), row.names=TRUE)
  } else {
    message("  Нет значимых генов для контраста ", cname)
  }
  
  if (nrow(sig) > 0) {
    up <- sig[sig$log2FoldChange > 1, , drop=FALSE]
    down <- sig[sig$log2FoldChange < -1, , drop=FALSE]
    run_kegg <- function(genes_df, prefix) {
      ids <- na.omit(mapIds(org.Hs.eg.db, keys=rownames(genes_df),
                            column="ENTREZID", keytype="ENSEMBL", multiVals="first"))
      if (length(ids) == 0) {
        message("    Нет ENTREZID для ", prefix, " в ", cname)
        return(NULL)
      }
      kk <- tryCatch({
        enrichKEGG(gene=ids, organism="hsa", pvalueCutoff=0.05, qvalueCutoff=0.1)
      }, error = function(e) {
        message("    enrichKEGG failed for ", prefix, ": ", e$message); return(NULL)
      })
      if (!is.null(kk) && nrow(as.data.frame(kk))>0) {
        write.csv(as.data.frame(kk), file=file.path(outdir, paste0("KEGG_", prefix, "_", cname, ".csv")), row.names=FALSE)
        try({
          png(file.path(outdir, paste0("KEGG_bar_", prefix, "_", cname, ".png")), width=1000, height=1500, res=120)
          barplot(kk, showCategory=min(50, nrow(as.data.frame(kk))),
                  title=paste0("KEGG ", prefix, ": ", cname))
          dev.off()
        }, silent = TRUE)
        try({
          png(file.path(outdir, paste0("KEGG_dot_", prefix, "_", cname, ".png")), width=1000, height=1500, res=120)
          dotplot(kk, showCategory=min(50, nrow(as.data.frame(kk))),
                  title=paste0("KEGG dot ", prefix, ": ", cname))
          dev.off()
        }, silent = TRUE)
      } else {
        message("    Нет обогащённых путей KEGG для ", prefix, " в ", cname)
      }
      return(kk)
    }
    message("  KEGG up...")
    run_kegg(up, "up")
    message("  KEGG down...")
    run_kegg(down, "down")
  }
  
  message("  Запускаем GSEA Hallmark...")
  res_rn <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$gene_symbol), ]
  if (nrow(res_rn)>0) {
    res_rn <- res_rn[!duplicated(res_rn$gene_symbol), ]
    gene_list <- sort(setNames(res_rn$log2FoldChange, res_rn$gene_symbol), decreasing=TRUE)
    hallmark <- msigdbr(species="Homo sapiens", category="H")[, c("gs_name","gene_symbol")]
    colnames(hallmark) <- c("ID","gene")
    gsea_res <- tryCatch({
      GSEA(geneList=gene_list, TERM2GENE=hallmark, pvalueCutoff=1,
           pAdjustMethod="BH", minGSSize=15, maxGSSize=500, verbose=FALSE)
    }, error = function(e) {
      message("    GSEA failed: ", e$message); return(NULL)
    })
    if (!is.null(gsea_res)) {
      res_all <- gsea_res@result
      res_clean <- subset(res_all, is.finite(NES) & !is.na(p.adjust))
      gsea_res@result <- res_clean
      sig_gsea <- subset(res_clean, p.adjust<=0.25 & abs(NES)>=1)
      write.csv(res_all, file=file.path(outdir, paste0("GSEA_full_", cname, ".csv")), row.names=FALSE)
      if (nrow(sig_gsea)>0) {
        write.csv(sig_gsea, file=file.path(outdir, paste0("GSEA_sig_", cname, ".csv")), row.names=FALSE)
        try({
          png(file.path(outdir, paste0("GSEA_dot_sig_", cname, ".png")), width=2500, height=2500, res=300)
          dotplot(gsea_res, showCategory=min(20, nrow(sig_gsea)),
                  title=paste0("GSEA sig: ", cname))
          dev.off()
        }, silent = TRUE)
        try({
          png(file.path(outdir, paste0("GSEA_ridge_sig_", cname, ".png")), width=2500, height=3000, res=300)
          ridgeplot(gsea_res, showCategory=min(20, nrow(sig_gsea))) +
            labs(title=paste0("GSEA ridge: ", cname))
          dev.off()
        }, silent = TRUE)
        emt_row <- sig_gsea[sig_gsea$Description=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
        if (nrow(emt_row)==1) {
          emt_id <- emt_row$ID
          try({
            png(file.path(outdir, paste0("GSEA_EMT_", cname, ".png")), width=2000, height=1500, res=250)
            gseaplot2(gsea_res, geneSetID=emt_id,
                      title=paste0("EMT GSEA (NES=", round(emt_row$NES,2),
                                   ", FDR=", round(emt_row$p.adjust,3), ") ", cname))
            dev.off()
          }, silent = TRUE)
          le_str <- gsea_res@result$core_enrichment[gsea_res@result$ID==emt_id]
          if (length(le_str)==1 && nzchar(le_str)) {
            le_genes <- strsplit(le_str, "/")[[1]]
            symbol2ens <- mapIds(org.Hs.eg.db, keys=le_genes,
                                 column="ENSEMBL", keytype="SYMBOL", multiVals="first")
            ens_keep <- symbol2ens[!is.na(symbol2ens) & symbol2ens %in% rownames(vsd_obj)]
            if (length(ens_keep)>0) {
              mat <- assay(vsd_obj)[ens_keep, , drop=FALSE]
              rownames(mat) <- names(ens_keep)
              keep_samps <- colnames(vsd_obj)[coldata$condition %in% c(g1,g2)]
              mat2 <- mat[, keep_samps, drop=FALSE]
              if (ncol(mat2)>0) {
                annotation_col2 <- data.frame(condition=coldata$condition[keep_samps])
                rownames(annotation_col2) <- keep_samps
                row_sds <- apply(mat2, 1, sd)
                mat_filt <- mat2[row_sds>0, , drop=FALSE]
                if (nrow(mat_filt)>0) {
                  try({
                    png(file.path(outdir, paste0("GSEA_leadingEdge_EMT_", cname, ".png")),
                        width=2000, height=3200, res=300)
                    pheatmap(mat_filt, scale="row",
                             annotation_col=annotation_col2,
                             main=paste0("Leading Edge EMT ", cname))
                    dev.off()
                  }, silent = TRUE)
                }
              }
            }
          }
        }
      } else {
        message("    Нет значимых путей GSEA для ", cname)
      }
    }
  } else {
    message("    Нечего запускать GSEA: отсутствуют данные")
  }
  message("Контраст ", cname, " обработан.\n")
}

# 10. Запуск по всем контрастам
for (ctr in contrast_list) {
  outdir <- file.path(output_base, ctr$name)
  run_contrast_analysis(dds, vsd, coldata, ctr, gene_symbols_all, gene_desc_all, outdir)
}

message("Все контрасты завершены.")
