library(shiny)
library(shiny.router)
library(shiny.tailwind)
library(shinyjs)
library(shinycssloaders)
## Specify your h5 input file and filter in only cell containing droplets
library(Matrix)
library(Seurat)
library(cluster)
library(fitdistrplus)
library(data.table)
library(dplyr)
library(ggrepel)
library(reticulate)
use_python("/home/poqopo/anaconda3/bin/python")

source("/data/project/scitoseq/script/Final/Function/Make_hashtag_obj.R")
source("/data/project/scitoseq/script/Final/Function/Make_resolved_obj.R")
# source("/data/project/scitoseq/script/Final/Function/Unimodal_Integration.R")
# source("/data/project/scitoseq/script/Final/Function/Multimodal_CLR_Integration.R")
# source("/data/project/scitoseq/script/Final/Function/Multimodal_DSB_Integration.R")



source("/data/project/scitoseq/script/Final/Myapp/TutorialPage.R")
source("/data/project/scitoseq/script/Final/Myapp/ContactPage.R")
source("/data/project/scitoseq/script/Final/Myapp/HelpPage.R")
source("/data/project/scitoseq/script/Final/Myapp/HomePage.R")
source("/data/project/scitoseq/script/Final/Myapp/RunPage.R")


ui = fluidPage(
  use_tailwind(),
  useShinyjs(),
  title = "WebSCITO",
  router_ui(
    route("/", root_page),
    route("run", run_page),
    route("tutorial", tutorial_page),
    route("contact", contact_page),
    route("help", help_page)
  )
)

server = function(input, output, session) {
  router_server()
  output$helpmeta <- renderDataTable({
      read.table("/data/project/scitoseq/script/Final/Sample/sample_metadata.csv")
    }, options = list(searching = FALSE, info = FALSE))
  hideElement("shashtagloading")
  hideElement("smultimodalloading")
  hideElement("sresolveloading")
  hideElement("sintegrationloading")
  hideElement("hashtagloading")
  hideElement("multimodalloading")
  hideElement("resolveloading")
  hideElement("integrationloading")

  observeEvent(input$tutorial_start, {
    dir <- input$Samples
    if (dir == "/data/project/scitoseq/script/Final/Sample/Multi_4AB_4Pool/") {
      showElement("smultimodalloading")
      py_run_string("
import muon as mu
mdata = mu.read_h5mu('/data/project/scitoseq/script/Final/Sample/Multi_4AB_4Pool/umap.h5mu')
rna = mdata['rna'].copy()
prot = mdata['prot'].copy()
")
      multi <<- TRUE
    } else {
      showElement("shashtagloading")
      hashtag_list <<- readRDS(paste(dir,"hashtag_list.rds", sep=''))
      batchNum <<- 10 
      abNum <<- 28
      multi <<- FALSE
    }
    if (!multi) {
      output$sdoublet <- renderPlot({
        doublet_df <- as.data.frame(table(hashtag_list$ab_hashtag_obj$HTO_classification.global))
        
        plot <- doublet_df %>% 
          mutate(csum = rev(cumsum(rev(Freq))), 
                 pos = Freq/2 + lead(csum, 1),
                 pos = if_else(is.na(pos), Freq/2, pos),
                 percentage = Freq/sum(Freq)) %>% 
          ggplot(aes(x = "", y = Freq, fill = Var1)) + 
          geom_col(width = 1, color = 1) +
          geom_label_repel(aes(y = pos,
                               label = glue::glue("{Freq} ({scales::percent(percentage)})"),
                               fill = Var1),
                           size = 4,
                           nudge_x = 1,
                           show.legend = FALSE) +
          labs( fill = "Doublet Type" ) +
          coord_polar(theta = "y") +
          ggtitle("Doublet Information")
        theme_void()
        
        plot
      })
      output$sbatch <- renderPlot({
        max_ID_df <- as.data.frame(table(hashtag_list$ab_hashtag_obj$HTO_maxID))
        
        plot <- max_ID_df %>% 
          mutate(csum = rev(cumsum(rev(Freq))), 
                 pos = Freq/2 + lead(csum, 1),
                 pos = if_else(is.na(pos), Freq/2, pos),
                 percentage = Freq/sum(Freq)) %>% 
          ggplot(aes(x = "", y = Freq, fill = Var1)) + 
          geom_col(width = 1, color = 1) +
          geom_label_repel(aes(y = pos,
                               label = glue::glue("{Freq} ({scales::percent(percentage)})"),
                               fill = Var1),
                           size = 4,
                           nudge_x = 1,
                           show.legend = FALSE) +
          labs( fill = "Batch Type" ) +
          coord_polar(theta = "y") +
          ggtitle("Batch Percentage")
        theme_void()
        plot
      })
      hideElement("shashtagloading")
      output$sbasicInfo <- renderUI({
        div(
          div(class = "flex gap-x-10",
              h1(class= "py-10 font-extrabold text-[30px] text-start", "Basic Information of your Input"),
              div(class = "h-fit my-auto",downloadButton('downloadsHashtag', 'Download Hashtag Object'))
          ),
          
          div(class = "w-4/5 mx-auto flex place-content-evenly py-5 gap-x-10",
              plotOutput("sdoublet", height = 400, width = 600) %>% withSpinner(),
              plotOutput("sbatch", height = 400, width = 600) %>% withSpinner()
          )
        )
      })
      output$sridgeplot <- renderPlot({
        RidgePlot(hashtag_list$ab_hashtag_obj, assay = "HTO", features = rownames(hashtag_list$ab_hashtag_obj[["HTO"]])[1:batchNum], ncol = 4)
      })
      output$sridgeplotUI <- renderUI({
        div(
          h1(class= "py-10 font-extrabold text-[30px] text-start", "Ridgeplot by batch"),
          div(class = "mx-auto flex py-5",
              plotOutput("sridgeplot", width = 1200, height = 400 * (batchNum %/% 4)) %>% withSpinner()
          ),
          div(class ="w-full mx-auto my-10 flex place-content-center",
              actionButton(class="mx-auto w-[300px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-xl",
                           inputId = "spca_umap", label = "Draw PCA and UMAP")
          ),
        )
      })
      output$downloadsHashtag <- downloadHandler(
        filename = function() { 
          paste("ab_hashtag_obj", Sys.Date(), ".rds", sep="")
        },
        content = function(file) {
          saveRDS(hashtag_list$ab_hashtag_obj, file)
        })
      
    } else {
      hideElement("smultimodalloading")
      output$sridgeplotUI <- renderUI({
        div(
          div(class ="w-full mx-auto my-10 flex place-content-center",
              actionButton(class="mx-auto w-[300px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-xl",
                           inputId = "spca_umap", label = "Draw PCA and UMAP")
          ),
        )
      })
    }
})
  observeEvent(input$spca_umap, {
    showElement("sresolveloading")
    if(!multi){
      result <- make_resolved_obj(ab = hashtag_list$ab, tx = hashtag_list$tx, discrete = hashtag_list$discrete,
                                  ab_hashtag_obj = hashtag_list$ab_hashtag_obj,
                                  abNum = abNum, batchNum = batchNum)
      selected_ab_cell_matrix <<- after_processing(result$selected_ab_cell_matrix)
      meta_ab_cell_matrix <<- result$meta_ab_cell_matrix
      ptn_names <<- read.table("/data/project/scitoseq/script/Final/Sample/200K_28AB_10Pool/ptn_names.csv")
      hideElement("sresolveloading")
      showElement("sintegrationloading")
      py_run_string("
import re
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import seaborn as sns
import leidenalg
import matplotlib.pyplot as plt

sc.settings.verbosity = 4

adt_2_df = r.selected_ab_cell_matrix
meta_2_df = r.meta_ab_cell_matrix
adt_2_df = adt_2_df.astype('float')
adt_adata = sc.AnnData(adt_2_df.T)
adt_adata.obs = meta_2_df

sc.pp.normalize_per_cell(adt_adata, counts_per_cell_after=1e4)

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split('(\\d+)', text)]


def normalize_byBatch(adt, batch_var='Batch', max_val=10):
    all_batches = list(set(adt.obs[batch_var]))
    all_batches.sort(key=natural_keys)
    ann_dats = []
    for b in all_batches:
        batch_adat = adt[adt.obs[batch_var] == b]
        sc.pp.log1p(batch_adat)
        sc.pp.scale(batch_adat, max_value=max_val)
        ann_dats.append(batch_adat)
    norm_concat = ann_dats[0].concatenate(ann_dats[1:len(ann_dats)+1])
    return norm_concat

adt_adata = normalize_byBatch(adt_adata)

adt_adata.var_names_make_unique()


ptn_names =r.ptn_names

adt_adata.var = ptn_names


sc.tl.pca(adt_adata, svd_solver='arpack', use_highly_variable=None)
sc.pp.neighbors(adt_adata, n_neighbors=10, n_pcs=15)
sc.tl.umap(adt_adata)
sc.tl.leiden(adt_adata)
")
      py_run_string("umapOptions = list(adt_adata.obs.columns)")
      umapOptions <<- py$umapOptions
      umapOptions <<- umapOptions[! umapOptions %in% c("nCount_RNA", "nCount_ADT", "nFeature_ADT", "nFeature_RNA")]
      py_run_string("featureOptions = list(adt_adata.var_names)")
      featureOptions <<- py$featureOptions
    }
    else{
      py_run_string("umapOptions = list(rna.obs.columns)")
      umapOptions <<- py$umapOptions
      umapOptions <<- umapOptions[! umapOptions %in% c("n_genes_by_counts", "total_counts", "total_counts_mt", "pct_counts_mt",
                                                       "total_counts_rb", "pct_counts_rb", "total_counts_hb", 'pct_counts_hb', 
                                                       'complexity', 'log_total_counts', 'log_n_genes_by_counts', 'log_pct_counts_mt', 
                                                       'log_complexity', 'log_pct_counts_hb')]
      py_run_string("featureOptions = list(prot.var_names)")
      featureOptions <<- py$featureOptions
      hideElement("sresolveloading")
    }

    selected_umap_option <- selectInput(inputId = "selected_umap_option",
                                        label = "Choose UMAP View",
                                        choices = umapOptions,
                                        selected = "leiden")
    output$spca <- renderPlot({
      if(multi) {
        py_run_string("sc.pl.pca_variance_ratio(rna, log=True)")
      }
      else {
        py_run_string("sc.pl.pca_variance_ratio(adt_adata, log=True)")
      }
    })
    output$sumap <- renderPlot ({
      if(multi) {
        py_run_string(paste0("sc.pl.umap(rna, color=['",input$selected_umap_option,"'], show = False)"))
        py_run_string("plt.subplots_adjust(right=0.7)")
        py_run_string("plt.show()")
      }
      else {
        py_run_string(paste0("sc.pl.umap(adt_adata, color=['",input$selected_umap_option,"'], show = False)"))
        py_run_string("plt.subplots_adjust(right=0.7)")
        py_run_string("plt.show()")
      }
    })
    hideElement("sintegrationloading")
    output$spcaumap <- renderUI({
      div(
        div(class = 'flex gap-x-10',
            h1(class= "py-10 font-extrabold text-[30px] text-start", "PCA and UMAP"),
            selected_umap_option,
            div(class = "h-fit my-auto",downloadButton('downloadsnResolved', 'Download Resolved Object'))
        ),
        div(class = "w-4/5 mx-auto flex place-content-evenly py-5 gap-x-10",
            plotOutput("spca", height = 400, width = 600) %>% withSpinner(),
            plotOutput("sumap", height = 400, width = 600) %>% withSpinner()
        ),
        div(class = "w-full mx-auto flex place-content-evenly",
            actionButton(class="mx-10 w-full max-w-[200px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-xl",
                         inputId = "sDotPlotBtn", label = "DotPlot show"),
            actionButton(class="mx-10 w-full max-w-[200px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-xl",
                         inputId = "sFeaturePlotBtn", label = "FeaturePlot show")
        )
      )
    })
    output$downloadsnResolved <- downloadHandler(
      filename = function() { 
        paste("dataset-", Sys.Date(), ".h5ad", sep="")
      },
      content = function(file) {
        if(multi) {
          py_run_string("mdata.write(r.file)")
          
        }
        else {
          py_run_string("adt_adata.write_h5ad(r.file)")
        }
      })
    
  })
  observeEvent(input$sDotPlotBtn, {
    selected_dot_option <- selectInput(inputId = "selected_dot_option",
                                       label = "Choose Group View",
                                       choices = umapOptions,
                                       selected = "leiden")
    selected_feature_option <- selectInput(inputId = "selected_feature_option",
                                       label = "Choose Feature View",
                                       choices =  featureOptions,
                                       multiple = TRUE,
                                       selected = featureOptions[1:10])
    if (multi) {
      py_run_string("geneOptions = list(rna.var_names)")
      geneOptions <<- py$geneOptions
      selected_rna_feature_option <- selectInput(inputId = "selected_rna_feature_option",
                                             label = "Choose Feature View",
                                             choices = geneOptions,
                                             multiple = TRUE,
                                             selected = geneOptions[1:10])
    } else {selected_rna_feature_option <- div()}
    
    output$sdotplot <- renderPlot({
      if (multi) {
        py_run_string("
prot.obs['leiden'] = rna.obs['leiden']
ncols=2
nrows=1
figsize=4
wspace=0.5
fig,axs = plt.subplots(nrows=nrows, ncols=ncols,
                       figsize=(ncols*figsize+figsize*wspace*(ncols-1),nrows*figsize))
plt.subplots_adjust(wspace=wspace)
")
        py_run_string(paste0("sc.pl.dotplot(rna, groupby = '",input$selected_dot_option,"', var_names = [",
                             paste(shQuote(input$selected_rna_feature_option), collapse = ', '),"],show=False,ax=axs[0])"))
        py_run_string(paste0("sc.pl.dotplot(prot, groupby = '",input$selected_dot_option,"', var_names = [",
                             paste(shQuote(input$selected_feature_option), collapse = ', '),"],ax=axs[1])"))
        
      }
      else {
        py_run_string(paste0("sc.pl.dotplot(adt_adata, groupby = '",input$selected_dot_option,"', var_names = [",
                             paste(shQuote(input$selected_feature_option), collapse = ', '),"])"))
      }
    })
    
    
    output$sdotplotUI <- renderUI({
      div(
        div(class = 'flex py-10',
            h1(class= " font-extrabold text-[30px] text-start", "DotPlot"),
            div(class = "mx-10", selected_dot_option),
            div(class = "mx-10", selected_feature_option),
            if(multi == TRUE) {
              selected_rna_feature_option
            }
        ),
        plotOutput("sdotplot", width = "100%", height = "500px") %>% withSpinner()
      )
    })
    
  })
  observeEvent(input$sFeaturePlotBtn, {
      selected_adt_feature_option <- selectInput(inputId = "selected_adt_feature_option",
                                                 label = "Choose adt Feature View",
                                                 choices = featureOptions,
                                                 selected =featureOptions[1])
      if (multi == TRUE) {
        py_run_string("geneOptions = list(rna.var_names)")
        geneOptions <<- py$geneOptions
        selected_rna_featureplot_option <- selectInput(inputId = "selected_rna_featureplot_option",
                                                   label = "Choose rna Feature View",
                                                   choices = geneOptions,
                                                   selected = geneOptions[1])
      }
      output$sfeaturePlot <- renderPlot({
        if (multi == TRUE) {
          py_run_string(paste0("
import colorcet as cc
mu.pl.embedding(
    mdata,
    basis='rna:umap',
    color=[",shQuote(input$selected_adt_feature_option),',',shQuote(input$selected_rna_featureplot_option),"],
    cmap=cc.cm.gouldian,
    palette=cc.glasbey_hv)
"))
        } else { 
          py_run_string(paste0("sc.pl.umap(adt_adata, color=['",input$selected_adt_feature_option,"'],vmin=-1,vmax=1 )"))
          
        }
      })
      output$sfeatureplotUI <- renderUI({
        div(
          div(class = 'flex py-10',
              h1(class= " font-extrabold text-[30px] text-start", "FeaturePlot"),
              div(class = "mx-20", selected_adt_feature_option),
              if (multi == TRUE) {
                div(class = "mx-20", selected_rna_featureplot_option)
              }
          ),
          plotOutput("sfeaturePlot", width = 600 * (multi + 1), height = 400) %>% withSpinner()
        )
      })
    })
  
  observeEvent(input$start, {
    ext <- tools::file_ext(input$raw$datapath)
    if (ext == "h5") {
      raw_input <<- Read10X_h5(input$raw$datapath)
    } else if(ext == "rds")  {
      raw_input <<- readRDS(input$raw$datapath)
    } else {
      validate("Please Provide Valid file.")
    }
    abNum <<- input$abNum
    batchNum <<- as.integer(input$batchNum)
    if (!input$multimodal) {
      showElement("hashtagloading")
      hashtag_list <<- make_hashtag_obj(raw_input, abNum, batchNum)
      output$doublet <- renderPlot({
        doublet_df <- as.data.frame(table(hashtag_list$ab_hashtag_obj$HTO_classification.global))
        
        plot <- doublet_df %>% 
          mutate(csum = rev(cumsum(rev(Freq))), 
                 pos = Freq/2 + lead(csum, 1),
                 pos = if_else(is.na(pos), Freq/2, pos),
                 percentage = Freq/sum(Freq)) %>% 
          ggplot(aes(x = "", y = Freq, fill = Var1)) + 
          geom_col(width = 1, color = 1) +
          geom_label_repel(aes(y = pos,
                               label = glue::glue("{Freq} ({scales::percent(percentage)})"),
                               fill = Var1),
                           size = 4,
                           nudge_x = 1,
                           show.legend = FALSE) +
          labs( fill = "Doublet Type" ) +
          coord_polar(theta = "y") +
          ggtitle("Doublet Information")
        theme_void()
        
        plot
      })
      output$batch <- renderPlot({
        max_ID_df <- as.data.frame(table(hashtag_list$ab_hashtag_obj$HTO_maxID))
        
        plot <- max_ID_df %>% 
          mutate(csum = rev(cumsum(rev(Freq))), 
                 pos = Freq/2 + lead(csum, 1),
                 pos = if_else(is.na(pos), Freq/2, pos),
                 percentage = Freq/sum(Freq)) %>% 
          ggplot(aes(x = "", y = Freq, fill = Var1)) + 
          geom_col(width = 1, color = 1) +
          geom_label_repel(aes(y = pos,
                               label = glue::glue("{Freq} ({scales::percent(percentage)})"),
                               fill = Var1),
                           size = 4,
                           nudge_x = 1,
                           show.legend = FALSE) +
          labs( fill = "Batch Type" ) +
          coord_polar(theta = "y") +
          ggtitle("Batch Percentage")
        theme_void()
        plot
      })
      hideElement("hashtagloading")
      output$basicInfo <- renderUI({
        div(
          h1(class= "py-10 font-extrabold text-[30px] text-start", "Basic Information of your Input"),
          div(class = "w-4/5 mx-auto flex place-content-evenly py-5 gap-x-10",
              plotOutput("doublet", height = 400, width = 600) %>% withSpinner(),
              plotOutput("batch", height = 400, width = 600) %>% withSpinner()
          )
        )
      })
      output$ridgeplot <- renderPlot({
        RidgePlot(hashtag_list$ab_hashtag_obj, assay = "HTO", features = rownames(hashtag_list$ab_hashtag_obj[["HTO"]])[1:batchNum], ncol = 4)
      })
      output$ridgeplotUI <- renderUI({
        div(
          h1(class= "py-10 font-extrabold text-[30px] text-start", "Ridgeplot by batch"),
          div(class = "mx-auto flex py-5",
              plotOutput("ridgeplot", width = 1200, height = 400 * (batchNum %/% 4)) %>% withSpinner()
          ),
          div(class ='w-full flex place-content-between',
              div(class ='m-auto w-4/5 flex place-content-evenly gap-x-2',
                  fileInput(inputId="Metadata", label='Meta data file(optional)', width = "30%"), #Need to make accept rule 
                  fileInput(inputId="antibody_names", label='antibody Name file(optional)', width = "30%"), #Need to make accept rule
              ),
              actionButton(class="mx-auto w-[300px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-xl",
                           inputId = "pca_umap", label = "Draw PCA and UMAP")
          )
        )
      })
    } else {
      output$basicInfo <- renderUI({
        div(
          p(class = "font-bold text-[20px]", "For multiome data, we need more information"),
          p(class = "font-bold text-[20px]", "It takes 15~30 minutes to be done."),
          div(class ='w-full flex place-content-between',
              div(class ='mx-auto my-10 w-4/5 flex place-content-evenly gap-x-2',
                  fileInput(inputId="Metadata", label='Batch Meta data', width = "30%"), #Need to make accept rule 
                  fileInput(inputId="antibody_names", label='Antibody Name file(optional)', width = "30%"), #Need to make accept rule
              ),
              actionButton(class="mx-auto w-[300px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-xl",
                           inputId = "pca_umap", label = "start multiome analysis")
          )
        )
      })
    }
    
  }) 
  observeEvent(input$pca_umap, {
    if (is.null(input$Metadata)) {aux <- NULL} else { aux <- read.table(input$Metadata$datapath)}
    if (is.null(input$antibody_names)) {ptn_names <- NULL} else { ptn_names <- read.table(input$antibody_names$datapath)}
    if (input$multimodal) {
      showElement("multimodalloading")
      tx <- raw_input$`Gene Expression`[, which(Matrix::colSums(raw_input$`Gene Expression`) > 0)]
      ab <- raw_input$`Antibody Capture`[, which(Matrix::colSums(raw_input$`Antibody Capture`) > 0)]
      ccd_indices <- intersect(which(Matrix::colSums(tx) > 0), which(Matrix::colSums(ab) > 0))
      
      ab <- ab[, ccd_indices]
      tx <- tx[, ccd_indices]
      
      ab_resolved_list <- list()
      tx_resolved_list <- list()
      
      for (i in 1:as.integer(batchNum)) {
        start <- (1 + abNum*(i-1))
        end <- i*abNum
        
        ab_subset <- ab[start:end, which(meta$Batch == paste("Batch",i,sep=""))]
        ab_resolved_list[[i]] <- ab_subset
        
        tx_subset <- tx[, which(meta$Batch == paste("Batch",i,sep=""))]
        tx_resolved_list[[i]] <- tx_subset
      }
      
      ab_resolved <- do.call(cbind, ab_resolved_list)
      tx_resolved <- do.call(cbind, tx_resolved_list)
      
      cell_names <-make.names(colnames(ab_resolved), unique = TRUE)
      gene_names <- make.names(rownames(tx_resolved), unique = TRUE)
      ptn_names <- ptn_names$Antibody
      py_run_string("
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy import stats
import scanpy.external as sce
import muon as mu
from muon import MuData
from muon import prot as pt
import helper as _hp
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import anndata as ad
import colorcet as cc

tx_df = r.tx_resolved
tx_df = tx_df.astype('float')
rna = sc.AnnData(tx_df.T)
ab_df = r.ab_resolved
ab_df = ab_df.astype('float')
prot = sc.AnnData(ab_df.T)
rna.obs_names= r.cell_names
rna.var_names = r.gene_names
prot.obs_names = r.cell_names
prot.var_names = r.ptn_names
rna.obs['Pool'] = r.meta

og_shape = rna.shape

sc.pp.filter_cells(rna, min_counts=1, inplace=True)
rem_cell_shape = rna.shape
print(f'Removing obsolete cells: {og_shape} > {rem_cell_shape}')

sc.pp.filter_genes(rna, min_counts=1, inplace=True)
rem_cell_gene_shape = rna.shape
print(f'Removing obsolete genes: {rem_cell_shape} > {rem_cell_gene_shape}')
rna.var['mt'] = rna.var_names.str.startswith('MT-')
rna.var['rb'] = rna.var_names.str.startswith(('RPS', 'RPL'))
rna.var['hb'] = rna.var_names.str.contains(('^HB[^(P)]'))
sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=['mt', 'rb', 'hb'],
    percent_top=None,
    log1p=False,
    inplace=True
)
rna.obs['complexity'] = np.log10(rna.obs['n_genes_by_counts']) / np.log10(rna.obs['total_counts'])
rna.obs['complexity'] = rna.obs['complexity'].astype(np.float32)

fil_dict = dict(
    total_counts=(500, 30000),
    n_genes_by_counts=(500, 6000),
    complexity=(0.88,),
    pct_counts_mt=(7.5,),
    pct_counts_hb=(0.2,),
    gene_thres=20
)
print('Performing cell-level filtering...')

mu.pp.filter_obs(rna, 'total_counts', lambda x: x >= fil_dict['total_counts'][0])
print(f'Min UMI counts: {rna.X.shape}')

mu.pp.filter_obs(rna, 'total_counts', lambda x: x < fil_dict['total_counts'][1])
print(f'Max UMI counts: {rna.X.shape}')

mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: x >= fil_dict['n_genes_by_counts'][0])
print(f'Min gene counts: {rna.X.shape}')

mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: x < fil_dict['n_genes_by_counts'][1])
print(f'Max gene counts: {rna.X.shape}')

mu.pp.filter_obs(rna, 'pct_counts_mt', lambda x: x < fil_dict['pct_counts_mt'])
print(f'Max MT gene percentage: {rna.X.shape}')

mu.pp.filter_obs(rna, 'pct_counts_hb', lambda x: x < fil_dict['pct_counts_hb'])
print(f'Max HB gene percentage: {rna.X.shape}')

print('Performing gene-level filtering...')
gene_cell_thres = int(fil_dict['gene_thres'])
sc.pp.filter_genes(rna, min_cells=gene_cell_thres)
print(f'Final filtered data dims: {rna.X.shape}')

print('Recalculating QC metrics...')
sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=['mt', 'rb', 'hb'],
    percent_top=None,
    log1p=False,
    inplace=True
)
rna.obs['complexity'] = np.log10(rna.obs['n_genes_by_counts']) / np.log10(rna.obs['total_counts'])
cols = ['total_counts', 'total_counts_mt', 'total_counts_rb', 'total_counts_hb']
rna.obs[cols] = rna.obs[cols].astype(np.int32)
rna.obs['complexity'] = rna.obs['complexity'].astype(np.float32)
cols = ['total_counts', 'n_cells_by_counts']
rna.var[cols] = rna.var[cols].astype(np.int32)
cols = ['pct_dropout_by_counts']
rna.var[cols] = rna.var[cols].astype(np.float32)
for x in ('n_cells', 'n_counts'):
    if x in rna.var_keys():
        del(rna.var[x])
    if x in rna.obs_keys():
        del(rna.obs[x])
mdata = MuData({'rna': rna, 'prot': prot})
mdata.mod['rna'] = rna
mu.pp.intersect_obs(mdata)
  ")
py_run_string("
sc.pp.filter_genes(prot, min_counts=1, inplace=True)
rem_cell_gene_shape = prot.shape
print(f'Removing obsolete genes: {rem_cell_shape} > {rem_cell_gene_shape}')
pt.pp.clr(prot)
prot.layers['clr_counts'] = prot.X.copy()
prot.obs['total_clr_counts'] = prot.X.sum(axis=1)
mdata.mod['prot'] = prot
rna.layers['counts'] = rna.X.copy()
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
rna.layers['logcounts'] = rna.X.copy()

print('Computing PCA...')
sc.pp.pca(rna)

rna.layers['scaled'] = rna.X.copy()
rna.X = rna.layers['logcounts'].copy()
sc.pp.neighbors(rna, n_neighbors=10, n_pcs=15)
sc.tl.umap(rna)
sc.tl.leiden(rna)
mdata.mod['rna'] = rna
rna = mdata['rna'].copy()
prot = mdata['prot'].copy()
")

py_run_string("
mu.pl.embedding(
    mdata,
    basis='rna:umap',
    color=['CD4','CD4_protein'],
    cmap=cc.cm.gouldian,
    palette=cc.glasbey_hv)
")
    }
    else {
      showElement("resolveloading")
      resolved_obj <- make_resolved_obj(ab = hashtag_list$ab, tx = hashtag_list$tx,aux = aux, discrete = hashtag_list$discrete,
                                         ab_hashtag_obj = hashtag_list$ab_hashtag_obj, ptn_names = ptn_names, abNum = input$abNum, batchNum = input$batchNum)
      selected_ab_cell_matrix <<- result$selected_ab_cell_matrix
      meta_ab_cell_matrix <<- result$meta_ab_cell_matrix
      hideElement("resolveloading")
      showElement("integrationloading")
      py_run_string("
import re
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import seaborn as sns
import leidenalg
import matplotlib.pyplot as plt

sc.settings.verbosity = 4

adt_2_df = r.selected_ab_cell_matrix
meta_2_df = r.meta_ab_cell_matrix
adt_2_df = adt_2_df.astype('float')
adt_adata = sc.AnnData(adt_2_df.T)
adt_adata.obs = meta_2_df

sc.pp.normalize_per_cell(adt_adata, counts_per_cell_after=1e4)

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split('(\\d+)', text)]


def normalize_byBatch(adt, batch_var='Batch', max_val=10):
    all_batches = list(set(adt.obs[batch_var]))
    all_batches.sort(key=natural_keys)
    ann_dats = []
    for b in all_batches:
        batch_adat = adt[adt.obs[batch_var] == b]
        sc.pp.log1p(batch_adat)
        sc.pp.scale(batch_adat, max_value=max_val)
        ann_dats.append(batch_adat)
    norm_concat = ann_dats[0].concatenate(ann_dats[1:len(ann_dats)+1])
    return norm_concat

adt_adata = normalize_byBatch(adt_adata)

adt_adata.var_names_make_unique()


ptn_names =r.ptn_names

adt_adata.var = ptn_names


sc.tl.pca(adt_adata, svd_solver='arpack', use_highly_variable=None)
sc.pp.neighbors(adt_adata, n_neighbors=10, n_pcs=15)
sc.tl.umap(adt_adata)
sc.tl.leiden(adt_adata)
")
      py_run_string("umapOptions = list(adt_adata.obs.columns)")
      umapOptions <<- py$umapOptions
      umapOptions <<- umapOptions[! umapOptions %in% c("nCount_RNA", "nCount_ADT", "nFeature_ADT", "nFeature_RNA")]
      py_run_string("featureOptions = list(adt_adata.var_names)")
      featureOptions <<- py$featureOptions
    }
    selected_umap_option <- selectInput(inputId = "selected_umap_option",
                                        label = "Choose UMAP View",
                                        choices = umapOptions,
                                        selected = "seurat_clusters")
    output$pca <- renderPlot({
      if(input$multimodal) {
        py_run_string("sc.pl.pca_variance_ratio(adt_adata, log=True)")
      }
      else {
        py_run_string("sc.pl.pca_variance_ratio(rna, log=True)")
      }
    })
    output$umap <- renderPlot ({
      if(input$multimodal) {
        py_run_string(paste0("sc.pl.umap(rna, color=['",input$selected_umap_option,"'], show = False)"))
        py_run_string("plt.subplots_adjust(right=0.7)")
        py_run_string("plt.show()")
      }
      else {
        py_run_string(paste0("sc.pl.umap(adt_adata, color=['",input$selected_umap_option,"'], show = False)"))
        py_run_string("plt.subplots_adjust(right=0.7)")
        py_run_string("plt.show()")
      }
    })
    hideElement("multimodalloading")
    hideElement("integrationloading")
    output$pcaumap <- renderUI({
      div(
        div(class = 'flex gap-x-10',
            h1(class= "py-10 font-extrabold text-[30px] text-start", "PCA and UMAP"),
            selected_umap_option,
            div(class = "h-fit my-auto",downloadButton('downloadnResolved', 'Download Resolved Object'))
        ),
        div(class = "w-4/5 mx-auto flex place-content-evenly py-5 gap-x-10",
            plotOutput("pca", height = 400, width = 600) %>% withSpinner(),
            plotOutput("umap", height = 400, width = 600) %>% withSpinner()
        ),
        div(class = "w-full mx-auto flex place-content-evenly",
            actionButton(class="mx-10 w-full max-w-[200px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-xl",
                         inputId = "DotPlotBtn", label = "DotPlot show"),
            actionButton(class="mx-10 w-full max-w-[200px] h-fit p-5 text-white text-[20px] font-bold bg-black hover:bg-slate-50 rounded-xl",
                         inputId = "FeaturePlotBtn", label = "FeaturePlot show")
        )
      )
      
    })
    
    output$downloadnResolved <- downloadHandler(
      filename = function() { 
        paste("n_resolved_obj-", Sys.Date(), ".rds", sep="")
      },
      content = function(file) {
        saveRDS(n_resolved_obj)
      })
    
  })
  observeEvent(input$DotPlotBtn, {
    selected_dot_option <- selectInput(inputId = "selected_dot_option",
                                       label = "Choose Group View",
                                       choices = umapOptions,
                                       selected = "leiden")
    selected_feature_option <- selectInput(inputId = "selected_feature_option",
                                           label = "Choose Feature View",
                                           choices =  featureOptions,
                                           multiple = TRUE,
                                           selected = featureOptions[1:10])
    if (input$multimodal == TRUE) {
      py_run_string("geneOptions = list(rna.var_names)")
      geneOptions <<- py$geneOptions
      selected_rna_feature_option <- selectInput(inputId = "selected_rna_feature_option",
                                                 label = "Choose Feature View",
                                                 choices = geneOptions,
                                                 multiple = TRUE,
                                                 selected = geneOptions[1:10])
    } else {selected_rna_feature_option <- div()}
    
    output$dotplot <- renderPlot({
      if (input$multimodal) {
        py_run_string("
prot.obs['leiden'] = rna.obs['leiden']
ncols=2
nrows=1
figsize=4
wspace=0.5
fig,axs = plt.subplots(nrows=nrows, ncols=ncols,
                       figsize=(ncols*figsize+figsize*wspace*(ncols-1),nrows*figsize))
plt.subplots_adjust(wspace=wspace)
")
        py_run_string(paste0("sc.pl.dotplot(rna, groupby = '",input$selected_dot_option,"', var_names = [",
                             paste(shQuote(input$selected_rna_feature_option), collapse = ', '),"],show=False,ax=axs[0])"))
        py_run_string(paste0("sc.pl.dotplot(prot, groupby = '",input$selected_dot_option,"', var_names = [",
                             paste(shQuote(input$selected_feature_option), collapse = ', '),"],ax=axs[1])"))
        
      }
      else {
        py_run_string(paste0("sc.pl.dotplot(adt_adata, groupby = '",input$selected_dot_option,"', var_names = [",
                             paste(shQuote(input$selected_feature_option), collapse = ', '),"])"))
      }
    })
    output$dotplotUI <- renderUI({
      div(
        div(class = 'flex py-10',
            h1(class= " font-extrabold text-[30px] text-start", "DotPlot"),
            div(class = "mx-10", selected_dot_option),
            div(class = "mx-10", selected_feature_option),
            if(input$multimodal == TRUE) {
              selected_rna_feature_option
            }
        ),
        plotOutput("dotplot", width = "100%", height = "700px") %>% withSpinner()
      )
    })
    
  })
  observeEvent(input$FeaturePlotBtn, {
    selected_adt_feature_option <- selectInput(inputId = "selected_adt_feature_option",
                                               label = "Choose ADT",
                                               choices = featureOptions,
                                               selected = featureOptions[1])
    if (input$multimodal == TRUE) {
      selected_rna_featureplot_option <- selectInput(inputId = "selected_rna_featureplot_option",
                                                     label = "Choose Feature View",
                                                     choices = geneOptions,
                                                     selected = geneOptions[1])
      
    }
    output$featurePlot <- renderPlot({
      if (input$multimodal == TRUE) {
        py_run_string(paste0("
mu.pl.embedding(
    mdata,
    basis='rna:umap',
    color=[",shQuote(input$selected_adt_feature_option),',',shQuote(input$selected_rna_featureplot_option),"],
    cmap=cc.cm.gouldian,
    palette=cc.glasbey_hv)
"))
      } else { 
        py_run_string(paste0("sc.pl.umap(adt_adata, color=['",input$selected_adt_feature_option,"'],vmin=-1,vmax=1 )"))
        
      }
    })
    output$featureplotUI <- renderUI({
      div(
        div(class = 'flex py-10',
            h1(class= " font-extrabold text-[30px] text-start", "FeaturePlot"),
            div(class = "mx-20", selected_adt_feature_option),
            if (input$multimodal == TRUE) {
              div(class = "mx-20", selected_rna_featureplot_option)
            }
        ),
        plotOutput("featurePlot", width = 600 * (input$multimodal + 1), height = 400) %>% withSpinner()
      )
    })
  })
  
}

options(shiny.autoreload = TRUE)
options(shiny.maxRequestSize=2048*(1024*1024))
options(shiny.port = 8989)
shinyApp(ui, server)