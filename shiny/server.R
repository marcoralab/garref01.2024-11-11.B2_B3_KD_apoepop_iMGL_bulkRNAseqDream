library(shiny)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(purrr)
library(variancePartition)
library(forcats)
library(plotly)

load("results.rda", envir = results <- new.env())
res <- results$res
genenames <- results$genenames
mds <- results$mds

mds_vars <- mds$mds |>
  select(-matches("^MD[0-9]{1,2}$"), -sample) |>
  names()

models_contrasts <- res |>
  imap_dfr(\(l, nm) tibble(model = nm, contrast = names(l$inputs$contrasts)))

contrasts_getmod <- function(cont, mod = "foobar", contrast_overlap = TRUE) {
  mc_overlap <- filter(models_contrasts, contrast == cont & model != mod)
  if (contrast_overlap) {
    mc <- mc_overlap
  } else {
    mc_notoverlap <- filter(models_contrasts, model != mod)
    mc <- bind_rows(mc_overlap, mc_notoverlap)
  }
  c(pull(mc, model), "same model")
}

get_contrast2 <- function(cont, mod1, mod2) {
  if (mod2 == "same model") {
    ctrsts <- names(res[[mod1]][["inputs"]][["contrasts"]])
    ctrsts[ctrsts != cont]
  } else {
    ctrsts <- names(res[[mod2]][["inputs"]][["contrasts"]])
    c(cont, ctrsts[ctrsts != cont])
  }
}

plot_volcano <- function(tab, t) {
  if (nrow(tab) > 0) {
    ggplot(data=tab, aes(x=logFC, y=-log10(adj.P.Val), color = diffexpressed, label=label)) +
      geom_point() + 
      scale_color_manual(values=c("#0006b1", "grey", "#bb0c00")) +
      geom_vline(xintercept=c(-1, 1), col="grey", linetype = 'dashed') +
      geom_hline(yintercept=-log10(0.05), col="grey", linetype = 'dashed') +
      theme_minimal() +
      geom_label_repel( max.overlaps = 30) +
      ylab(expression("-log"[10]*"Adjusted p-value")) +
      xlab(expression("log"[2]*"FC")) +
      theme_pubr() +
      theme(legend.position ="none") +
      labs_pubr()+
      labs(title = t) +
      theme(plot.title = element_text(hjust = 0.4))
  } else {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, size = 8,
               label = sprintf("No significant genes for %s", t)) +
      theme_void()
  }
}

plot_gsea <- function(tab, t, st, sigthresh = 0.05) {
  if (is.null(tab)) {
    nogsea <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, size = 8,
               label = sprintf("No GSEA performed for %s", t)) +
      theme_void()
    return(nogsea)
  }
  tab |>
    filter(padj < sigthresh) |>
    arrange(desc(-NES)) |>
    mutate(pathway = fct_reorder(pathway, NES)) |>
    ggplot(aes(x=pathway, y=NES, size=size, color=padj)) + 
      geom_point(alpha = 0.8) + 
      coord_flip() + 
      xlab("") + 
      labs(title = t, 
           subtitle = st) +
      scale_color_gradient(low = "red",  high = "blue") +
      theme_classic()
}

plot_cor <- function(tab1, tab2, nm1, nm2, by, pfilt = 1, valcol, pcol) {
  if (valcol == "NES" && (is.null(tab1) || is.null(tab2))) {
    nogseatab = ifelse(is.null(tab1), nm1, nm2)
    nogsea <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, size = 8,
               label = sprintf("No GSEA performed for %s", nogseatab)) +
      theme_void()
    return(nogsea)
  }
  tab1 <- tab1 |>
    rename(!!nm1 := !!sym(valcol)) |>
    filter(!!sym(pcol) <= pfilt)
  tab2 <- tab2 |>
    rename(!!nm2 := !!sym(valcol)) |>
    filter(!!sym(pcol) <= pfilt)
  
  inner_join(tab1, tab2, by = by) |>
    ggplot(aes(x = !!sym(nm1), y = !!sym(nm2))) +
    geom_point() +
    stat_smooth(method = "lm", col = "blue", linetype = "dashed", se = FALSE,
                formula = y ~ x) +
    stat_cor(method = "pearson", size = 5) +
    labs(
      title = sprintf("%s in %s vs. %s", valcol, nm1, nm2),
      x = sprintf("%s (%s)", valcol, nm1),
      y = sprintf("%s (%s)", valcol, nm2)
    ) +
    theme_minimal()
}

make_dt_dgea <- function(fit, contrast) {
  fit[[contrast]] |>
    select("ENSEMBL ID" = gene_id, "Symbol" = gene_name, logFC,
           "x̄ Exp" = AveExpr, t, "P Value" = P.Value,
           "Q Value" = adj.P.Val, "β" = B) |>
    DT::datatable(
      extensions = 'Buttons',
      options = list(
        dom = 'Bfrtip',  # Include Buttons
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))) |>
    DT::formatSignif(c("logFC", "x̄ Exp", "t", "P Value", "Q Value", "β"), 3)
}

make_dt_gsea <- function(gsea) {
  if (is.null(gsea)) {return(NULL)}
  gsea |>
    select("Geneset" = geneset, "Pathway" = pathway, ES, NES,
           "P Value" = pval, "Q Value" = padj, size) |>
    DT::datatable(
      extensions = 'Buttons',
      options = list(
        dom = 'Bfrtip',  # Include Buttons
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))) |>
    DT::formatSignif(c("ES", "NES", "P Value", "Q Value"), 3)
}

filter_genesets <- function(dgea, ctrst, gset) {
  if (is.null(dgea)) {
    return(NULL)
  } else if (gset == "All genesets") {
    filter(dgea, contrast == ctrst)
  } else {
    filter(dgea, contrast == ctrst, geneset == gset)
  }
}

plot_mds <- function(mds, grp) {
  pct <- round(mds$var_explained * 100)[1:2]
  p <- mds$mds |>
    ggplot(aes(x = MD1, y = MD2, text = sample, color = !!sym(grp))) +
    geom_point() +
    labs(title = grp,
         x = sprintf("Leading logFC dim 1 (%i%%)", pct[1]),
         y = sprintf("Leading logFC dim 2 (%i%%)", pct[2]))
  
  ggplotly(p, tooltip = "text")
}

function(input, output, session) {
  output$notebook <- renderUI({
    tags$iframe(
      seamless="seamless",
      src="index.html",
      style="width:80vw;height:100vh;border: none;margin: 0;padding: 0;")
  })
  output$qc <- renderUI({
    tags$iframe(
      seamless="seamless",
      src="multiqc_report.html",
      style="width:80vw;height:100vh;border: none;margin: 0;padding: 0;")
  })
  updateSelectInput(session, "model", choices = names(res))
  updateSelectInput(session, "mds_var", choices = mds_vars)
  res_use <- reactive({
    req(input$model)
    res[[input$model]]
  })
  observeEvent(input$model, {
    ctrst <- names(res_use()[["inputs"]][["contrasts"]])
    updateSelectInput(session, "contrast", choices = ctrst)
  })
  observeEvent(list(input$contrast, input$contrast_overlap), {
    model2_chs <- contrasts_getmod(input$contrast, input$model,
                                   input$contrast_overlap)
    updateSelectInput(session, "model2", choices = model2_chs)
  })
  observeEvent(list(input$model2, input$contrast), {
    ctrst2_chs <- get_contrast2(input$contrast, input$model, input$model2)
    updateSelectInput(session, "contrast2", choices = ctrst2_chs)
  })
  observeEvent(input$model, {
    genesets_start <- res_use()[["out"]][["gsea"]] 
    if (is.null(genesets_start)) {
      genesets <- c("No GSEA performed")
    } else {
      genesets <- genesets_start |>
        distinct(geneset) |>
        pull(geneset) |>
        c("All genesets")
    }
    updateSelectInput(session, "geneset", choices = genesets)
  })
  res_tab <- reactive({
    req(input$model, input$contrast)
    res_use()[["out"]][["fit"]][[input$contrast]]
  })
  res_tab_mod2 <- reactive({
    req(input$model2, input$contrast2, input$comp_type)
    if (input$model2 == "same model") {
      mdl <- input$model
    } else {
      mdl <- input$model2
    }
    if (input$comp_type == "Differential Expression") {
      res[[mdl]][["out"]][["fit"]][[input$contrast2]]
    } else if (input$comp_type == "Geneset Enrichment") {
      gsea_unfilt <- res[[mdl]][["out"]][["gsea"]]
      if (is.null(gsea_unfilt)) {
        NULL
      } else {
        filter(gsea_unfilt, contrast == input$contrast2)
      }
    }
  })
  res_tab_gsea <- reactive({
    req(input$model, input$contrast)
    gsea_unfilt <- res[[input$model]][["out"]][["gsea"]]
    if (is.null(gsea_unfilt)) {
      NULL
    } else {
      filter_genesets(gsea_unfilt,input$contrast, input$geneset)
    }
  })
  res_tab_mod1 <- reactive({
    req(input$model2, input$contrast, input$comp_type)
    if (input$comp_type == "Differential Expression") {
      res_tab()
    } else if (input$comp_type == "Geneset Enrichment") {
      res_tab_gsea()
    }
  })
  cor_params <- reactive({
    req(input$model2, input$contrast, input$comp_type)
    if (input$comp_type == "Differential Expression") {
      list(by = "gene_name", valcol = "logFC", pcol = "P.Value")
    } else if (input$comp_type == "Geneset Enrichment") {
      list(by = "pathway", valcol = "NES", pcol = "pval")
    }
  })
  res_dt <- reactive(make_dt_dgea(res_use()[["out"]][["fit"]], input$contrast))
  res_dt_gsea <- reactive(make_dt_gsea(res_tab_gsea()))
  plot_contrasts <- reactive(plotContrasts(res_use()[["out"]][["cd"]]))
  plot_title <- reactive(sprintf("%s:\n%s", input$model, input$contrast))
  plot_mds_ <- reactive({
    req(input$mds_var)
    plot_mds(mds, input$mds_var)
  })
  output$model_spec <- renderText({
    res_use()[["inputs"]][["f"]] |>
      as.character() |>
      paste(collapse = " ")
  })
  output$volcano_plot = renderPlot({plot_volcano(res_tab(), plot_title())})
  output$contrast_plot <- renderPlot({plot_contrasts()})
  output$gsea_plot <- renderPlot({plot_gsea(res_tab_gsea(), plot_title(),
                                            input$geneset, input$gsea_qfilt)})
  output$dg_tab <- DT::renderDataTable({ res_dt() }, server = FALSE)
  output$gsea_tab <- DT::renderDataTable({ res_dt_gsea() }, server = FALSE)
  output$cor_plot <- renderPlot({
    plot_cor(res_tab_mod1(), res_tab_mod2(),
             nm1 = input$model, nm2 = input$model2,
             by = cor_params()$by, pfilt = input$comp_pfilt,
             valcol = cor_params()$valcol, pcol = cor_params()$pcol)
  })
  output$mds <- renderPlotly({plot_mds_()})
}
