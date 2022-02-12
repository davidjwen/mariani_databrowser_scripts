#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(tidyverse)
  library(cowplot)
  library(markdown)
  library(shinycssloaders)
  library(data.table)
  library(Matrix)
  library(patchwork)
  library(DT)
})

theme_set(theme_cowplot())
options(datatable.fread.datatable = FALSE)

exprs = readRDS("appdata/monocle_rna_counts.RDS")
umap_data = as_tibble(fread("appdata/monocle_clusters.csv"))
markers = as_tibble(fread("appdata/monocle_markers.csv"))
genes = sort(rownames(exprs))


# sort the factor levels
factor_levels = umap_data %>%
  select(cluster, cluster_num)
ordered_factor_levels = unique(factor_levels[order(factor_levels$cluster_num), ]$cluster)
umap_data$cluster = factor(umap_data$cluster, levels = ordered_factor_levels)

# clean up factor levels
rm(factor_levels)
rm(ordered_factor_levels)


# Define UI for data browser
ui <- fluidPage(

    # Application title
    titlePanel("Single Cell Data Browser"),

    fluidRow(
      column(4, 
        wellPanel(
             selectizeInput("gene_select",
                           "Select Gene",
                           choices = "Col1a1",
                           selected = "Col1a1"),
             radioButtons("split_on_treatment_bool",
                          "Violin Plot - Split on Treatment",
                          choices = c("Combine", "Split"),
                          selected = "Combine")
        ),
        DTOutput("data", width = "75%")
    ),
        
        # Show a plot of the cells
      mainPanel(
        fluidRow(withSpinner(plotOutput(outputId = "gene_plot"))),
        tabsetPanel(
          tabPanel("Clusters",
            fluidRow(withSpinner(plotOutput(outputId = "cluster_plot", click = "cluster_click")))
          ),
          tabPanel("Violin Plots",
                   fluidRow(withSpinner(plotOutput(outputId = "violin_plot"))))
        )
      )
    )
)

# Define server logic required to draw the data browser
server <- function(input, output, session) {
  updateSelectizeInput(session, 
                       "gene_select", 
                       choices = genes,
                       selected = "Col1a1",
                       server = T)
  
  # plots of clusters
  # in theory should never change; need to align clusters
  positions = umap_data %>%
    group_by(cluster) %>%
    summarise(x = median(x = UMAP_1),
              y = median(x = UMAP_2))
  
  plot_text = ggrepel::geom_text_repel(data = positions,
                                       mapping = aes_string(x = "x",
                                                            y = "y",
                                                            label = "cluster"),
                                       size = 3,
                                       inherit.aes = F)
  
  # base plot object
  cell_plot = ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
    labs(x = "UMAP 1", y = "UMAP 2") + 
    facet_grid(cols = vars(sample.class))

  cluster_plot = cell_plot + 
    geom_point(alpha = 0.75, size = 0.1, aes(color = cluster)) +
    labs(color = "cluster") + 
    plot_text + 
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.title = element_blank())


  # plots of genes
  exprs_levels = reactive({
    req(input$gene_select)
    return(exprs[rownames(exprs) == input$gene_select, ])
  })
  
  gene_plot = reactive({
    req(input$gene_select)
    
    if(sum(exprs_levels()) == 0){
      high_color = "#EEEEEE"
    } else{
      high_color = "#5e1ff0"
    }
    gene_plot = cell_plot + 
      geom_point(alpha = 0.5, size = 0.1, na.rm = T, aes(color = log10(exprs_levels() + 0.1))) +
      scale_color_gradient(low = "#EEEEEE", high = high_color) +
      plot_text + 
      labs(color = "log10(expression)")
    
    # need to manually align plots to make the clicky work
    aligned_plots = align_patches(cluster_plot, gene_plot)
    
    return(aligned_plots[[2]])
  })
  
  violin_plot = reactive({
    req(input$gene_select)
    req(input$split_on_treatment_bool)
    
    plot_data = tibble(
      cluster = umap_data$cluster,
      cluster_num = factor(umap_data$cluster_num),
      counts = exprs_levels(),
      exp = umap_data$sample.class
    )
    
    p = ggplot(plot_data, aes(x = cluster_num, y = log10(counts+.1), fill = cluster)) +
      geom_violin(scale = "width", color = "white", alpha=0.75) +
      labs(title = input$gene_select, x = "Cluster", y = "log10(expression + 0.1)") +
      scale_x_discrete() +
      theme(legend.title = element_blank())
    
    if(input$split_on_treatment_bool == "Split"){
      p = p + facet_grid(rows = vars(exp))
    }
    
    return(p)
  })
  

  output$gene_plot = renderPlot({
    gene_plot()
  })

  output$cluster_plot = renderPlot({
    cluster_plot
  })
  
  output$all_plot = renderPlot({
    all_plot()
  })
  
  output$violin_plot = renderPlot(violin_plot())


  ## Set up buffer, to keep the click.  
  click_saved <- reactiveValues(singleclick = NULL)
  
  ## Save the click, once it occurs.
  observeEvent(eventExpr = input$cluster_click, 
                handlerExpr = {
                 click_saved$singleclick = input$cluster_click 
              })
  
  
  # print the table of marker genes
  output_table = reactive({
    selected_clusters = nearPoints(umap_data, click_saved$singleclick) %>%
      select(cluster) %>%
      unique()
    
    cluster_markers = markers %>%
      select(cell_group, gene_short_name) %>%
      filter(cell_group %in% selected_clusters)
    
    # make the name a bit prettier
    names(cluster_markers) = c("Cluster", "Gene")

    return(cluster_markers)
  })
  
  output$data <- renderDT(
    output_table(), server = F, rownames = F, width = 400,
    options = list(autoWidth = T,
                   columnDefs = list(list(width = '200px', targets = "_all")),
                   scrollY = T)
  )

}


# Run the application 
shinyApp(ui = ui, server = server)
