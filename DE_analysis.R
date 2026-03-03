library(stringr)
library(tidyverse)
library(NotationConverter)
library(FindMarkersLoupe)
# devtools::install_local("/projectnb/wax-es/00_shinyapp/DEG/SOFT/FindMarkersLoupe", force =TRUE)
## TMP: original location FindMarkersLoupe, mm10 converter and Utils.R: /projectnb/wax-dk/max/RSRC/G190

DEBUG <- TRUE

if(DEBUG){
  seurat_obj <- readRDS("/projectnb/wax-dk/max/ATAC_TCPO/G183_G193_SC/output/module_3_outputs/extract_combine/rds/output.rds")
  s1 = "G183M2"
  s2 = "G183M1"
  c1 = 2
  c2 = 2
}

get_condition_by_sample_id <- function(sid, sobj) {
  sobj@meta.data %>% 
    select(sample_id, condition) %>% 
    filter(sample_id == sid) %>% 
    unique() %>% 
    pull(condition)
}

#' Calculate pairwise clusters comparisons
#'
#' @param seurat_obj 
#' @param s1 sample id 1 
#' @param s2 sample id 2
#' @param c1 cluster id 1
#' @param c2 cluster id 2
#'
#' @return list of the following objects:
#' - segex_output table with differential expression of clusters
#' - segex_filename output file for segex filename
#' - pdf with dotplot
#' 
#' @export
#'
#' @examples
compute_de <- function(seurat_obj, s1,c1,s2,c2, ix) {
  ## create new Idents and use them 
    ## processing of sample/sample
    ## using "sampleid_clusterid"
    new_idents_df <- seurat_obj@meta.data %>% select(CB, sample_id, seurat_clusters) %>% 
      mutate(new_idents = case_when(
        sample_id == s1 & seurat_clusters == c1 ~ "GROUP1",
        sample_id == s2 & seurat_clusters == c2 ~ "GROUP2",
        TRUE ~ "Other"
      )) %>% 
      select(CB,new_idents) 
  
    id.1 <- "GROUP1"
    id.2 <- "GROUP2"
  
  print("Number of cellls in each group")
  print(table(new_idents_df$new_idents))
  
  ## add new Idents to meta.data
  seurat_obj <- AddMetaData(seurat_obj, new_idents_df)
  
  ## activate new Idents
  Idents(seurat_obj) <- "new_idents"
  
  ## PROCESSING
  
  ## findMarkersLoupe (we cannot use the "short" Seurat object, because the entire 
  ## matrix must be used to calculate the intensities)
  markers_short <- FindMarkersLoupe(seurat_obj, id.1 = id.1, id.2 = id.2, formatted = "short")

  ## prepare Segex output data.frame
  segex_output <- exportToSegex(input_df = markers_short, from = "mm10")
  
  s1_label <- get_condition_by_sample_id(s1, seurat_obj)
  s2_label <- get_condition_by_sample_id(s2, seurat_obj)
  str_glue("{ix}_scLoupe_{s1}_{s1_label}_{c1}_vs_{s2}_{s2_label}_{c2}_DiffExp_IntronicMonoExonic.tsv")

  
  ## DOTPLOT
  ## subset of only required clusters (need to create 2 lines DotPlot)
  only_clusters_CB <- WhichCells(input_seurat, idents = c(id.1, id.2))
  ## short 'Seurat' object only CB related to clusters
  only_clusters_seurat <- subset(input_seurat, cells = only_clusters_CB)
  
  ## get top 30 genes from both sides
  ## TODO: double plot top +/-30 genes and +/-30 lncRNA
  ## TODO: need to create histograms which show shift of pvalue if we compare two
  ## different is size clusters (small clusters will not have any significant 
  ## genes by pvalue, but some genes will have good Log2FC)
  
  top60genes <- markers_short %>% 
    #filter( 0.5*((id.2.intensity+1)+ (id.2.intensity+1)) > 1, !grepl("lnc", gname))
    filter(!grepl("lnc", gname)) %>% 
    arrange(desc(log2_fold_change)) %>% 
    mp_extract_n_top_bottom(., n = 30) %>% 
    pull(gname)
  
  top60lncrna <- markers_short %>% 
    filter(grepl("lnc",gname)) %>% 
    arrange(desc(log2_fold_change)) %>% 
    mp_extract_n_top_bottom(.,  n = 30) %>% 
    pull(gname)
  
  ncells <- table(new_idents_df$new_idents)
  
  id.1.ncells <- ncells[[id.1n]]
  id.2.ncells <- ncells[[id.2n]]
  
  cols <- rev(c("#225ea8","#6baed6","#eff3ff","#fbb4b9","#fbb4b9"))
  top60genes_dotplot <- wrap_elements(mp_dotplot(only_clusters_seurat, 
                                                 top60genes, 
                                                 title = str_glue("{id.1}({id.1.ncells} cells) vs {id.2}({id.2.ncells} cells). Top 30 genes with average intensity > 1")) + 
                                        scale_colour_gradientn(colors = cols))
  
  top60lncrna_dotplot <- wrap_elements(mp_dotplot(only_clusters_seurat, 
                                                  top60lncrna, 
                                                  title = str_glue("{id.1}({id.1.ncells} cells) vs {id.2}({id.2.ncells} cells). Top 30 lncRNA")) + 
                                         scale_colour_gradientn(colors = cols))
  
  top60_dotplot <- top60genes_dotplot/top60lncrna_dotplot+plot_layout(heights = c(5,4))
  
  ## back to usual Idents and clean meta.data
  Idents(object = input_seurat) <- input_seurat@meta.data$seurat_clusters
  input_seurat@meta.data$new_idents <- NULL
  
  if (DEBUG) {
    list(segex_output = segex_output$segex,
         lp = markers_short,
         sr = markers_short_seurat,
         segex_filename = segex_fn,
         pdf = top60_dotplot)
  } else {
    list(segex_output = segex_output$segex,
         segex_filename = segex_fn,
         pdf = top60_dotplot)
  }
}

# TODO: average intensity does not work need replacing for something else
# 

# z.t3 <- compute_de(input_seurat, "G190M2", 1, "G190M2", 2, 1)
# res <- left_join(z.t3$lp,z.t3$sr)
# z.t3$pdf

# if (DEBUG){
#   z.t1 <- compute_de(input_seurat, "AGGR", 1, "AGGR", 2, 1)
#   z.t2 <- compute_de(input_seurat, "AGGR", 1, "G190M2", 0, 1)
#   z.t3 <- compute_de(input_seurat, "G190M2", 1, "G190M2", 2, 1)
#   z.t4 <- compute_de(input_seurat, "G190M2", 1, "G183M1", 1, 1)
# }
# 
# View(z.t1$tmp)
# summary(sd(z.t2$tmp$adjusted_p_value))
# #zscore <- (z.t1$tmp$adjusted_p_value-median(z.t1$tmp$adjusted_p_value))/sd(z.t1$tmp$adjusted_p_value)
# zscore <- z.t2$tmp$adjusted_p_value
# zscore <- zscore[zscore != 1.0]
# hist(zscore, breaks = 100)

# z.t2 <- compute_de(input_seurat, "AGGR", 1, "G190M2", 0, 1)
# z.t1$pdf
# z.tmp <- z.t2$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# View(z.tmp)  
#
#
# z.t4 <- compute_de(input_seurat, "G190M2", 1, "G183M1", 1, 1)
# z.tmp <- z.t4$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# 
# View(z.t4$tmp)
#
# 
# View(z.t2$segex_output)
# z.t1 <- compute_de(input_seurat, "AGGR", 0, "AGGR", 2, 1)
# z.tmp <- z.t1$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# View(z.t1$tmp)
# View(z.tmp)


print("Saving Segex TSVs and prepare list of pdfs to save")
pdf_list <- pmap(de_config_validated, function(SAMPLE_ID_1, SAMPLE_ID_2, CLUSTER_ID_1, CLUSTER_ID_2,    ix) {
  
  res_list <- compute_de(input_seurat, SAMPLE_ID_1, CLUSTER_ID_1, SAMPLE_ID_2, CLUSTER_ID_2, ix)
  if (!DEBUG)  write_tsv(res_list$segex_output, res_list$segex_filename, col_names = T)
  
  res_list$pdf
})


print("Saving PDFs")
pdf_list %>% 
  map(function(p){plot_grid(wrap_plots(p))}) %>% 
  marrangeGrob(nrow = 3, ncol = 1) %>%
  ggsave(filename = "dotplots_by_comparision.pdf", width = 15.50, height = 20)








