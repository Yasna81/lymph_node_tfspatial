library(shiny)
library(shinythemes)
library(bslib)
library(corrplot)
cor_mat <- readRDS("data/cor_sp.rds")
tf_list <- rownames(cor_mat)
tf_files <- c(
    FOSB = "FOSB.png",
    IRF4 = "IRF4.png",
    ATF3 = "ATF3.png",
    MAFF = "MAFF.png",
    NFE2L2 = "NFE2L2.png"
)

# ===== THEME =====
my_theme <- bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2C3E50",
    secondary = "#18BC9C",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter")
)

ui <- navbarPage(
    theme = my_theme,
    "TF Spatial Explorer",
    
    
    # ---------------------- 1. ABOUT  ----------------------
    tabPanel("About",
             fluidPage(
                 
                 tags$style(HTML("
            .qr-github-container {
                position: fixed;
                top: 15px;
                right: 20px;
                z-index: 1000;
                background: rgba(255, 255, 255, 0.92);
                padding: 8px 12px;
                border-radius: 8px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.15);
                font-size: 13px;
                display: flex;
                align-items: center;
                gap: 10px;
                backdrop-filter: blur(4px);
            }
            .qr-github-container img {
                width: 48px;
                height: 48px;
                border-radius: 4px;
            }
            .qr-github-container a {
                color: #2c3e50;
                text-decoration: none;
                font-weight: 500;
            }
            .qr-github-container a:hover {
                text-decoration: underline;
            }
        ")),
                 div(class = "qr-github-container",
                     tags$a(href = "https://github.com/Yasna81/lymph_node_tfspatial.git", target = "_blank",
                            "github.com/Yasna81/lymph_node_tfspatial"
                     )
                 ),
                 # Header
                 h2("Project: Spatial Regulatory Mapping of the Human Lymph Node",
                    style = "color: #2c3e50; font-weight: bold;"),
                 hr(),
                 
                 fluidRow(
                     # LEFT: Text content
                     column(8,
                            h4("Overview"),
                            p("This interactive dashboard demonstrates regulatory logic of the human lymph node tissue. 
           By integrating single-cell RNA-seq of human PBMCs with spatial transcriptomics (10x Visium) data of lymph node, 
           this tool visualizes how regulatory networks identified in the blood spatially organize to drive 
           distinct immune functions within the tissue architecture."),
                            br(),
                            
                            # === METHODOLOGY - now with custom dark-blue header ===
                            div(style = "border: 1px solid #ddd; border-radius: 8px; overflow: hidden; margin-bottom: 20px; box-shadow: 0 2px 6px rgba(0,0,0,0.1);",
                                div(style = "background-color: #2C3E50; color: white; padding: 12px 20px; font-size: 18px; font-weight: bold;",
                                    "Computational Method"),
                                div(style = "padding: 20px;",
                                    tags$ul(
                                        tags$li(strong("Cell clustering and annotation:"), "  Single-cell RNA-seq data from peripheral blood mononuclear cells (PBMCs) were processed using Seurat (v5). Cells were clustered using graph-based clustering and annotated into major immune cell types based on canonical markers."),
                                        tags$li(strong("Regulon activity inference:"), " Transcription factor (TF) regulon activity was inferred using pySCENIC on the PBMC scRNA-seq dataset. The top 5 most differentially active TFs across immune subsets were selected: ATF3, FOSB, MAFF, NFE2L2, and IRF4."),
                                        tags$li(strong("Spatial projection.Spatial projection of cell types and regulon activities:"), " Using anchor-based integration (Seurat), scRNA-seq-derived cell-type labels and pySCENIC regulon activity scores (AUC values) for the five selected TFs were transferred onto human lymph node spatial transcriptomics data (Visium). This yielded, for each spatial spot, probabilistic cell-type composition and TF regulon activity scores."),
                                        tags$li(strong("Niche definition:"), "Spatial spots were subjected to unsupervised clustering (graph-based) using the projected regulon activity scores of the five TFs (ATF3, IRF4, FOSB, NFE2L2, MAFF) as the primary feature space."),
                                        tags$li(strong("Cluster characterization and validation"), " Each TF-defined spatial niche was characterized by :
                                                (i) dominant transferred cell type(s) , 
                (ii) pathway enrichment using FGSEA on ranked gene lists derived from transferred single-cell signatures, and 
                (iii) identification of niche-specific marker genes through differential expression analysis between clusters.")
                                    )
                                )
                            ),
                            
                            # === KEY FINDINGS  ===
                            div(style = "border: 1px solid #ddd; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 6px rgba(0,0,0,0.1);",
                                div(style = "background-color: #2C3E50; color: white; padding: 12px 20px; font-size: 18px; font-weight: bold;",
                                    "Key Biological Discoveries"),
                                div(style = "padding: 20px;",
                                    tags$ul(
                                        tags$li(strong("Cell type mapping :"), "PBMC label transfer recoverd B cell follicles, T cell zones and DC niches in visium data"),
                                        tags$li(strong("TF Functions :"), " The analysis resolved three distinct functional niches:"),
                                        tags$ul(
                                            tags$li("Differentiation Zone: High IRF4, B-cell dominant (cluster 4)."),
                                            tags$li("Metabolic/stress Zone: High NFE2L2, Memory T-cell dominant (cluster 6)"),
                                            tags$li("Activation Zone: High ATF3/FOSB, Dendritic Cell dominant (clusters 1,2)"),
                                        ),
                                        tags$li(strong("Distinct stress and differentiation niches"), ": Two B-cell clusters showed opposite TF and pathway signatures."),
                                        tags$ul(
                                            tags$li(" Oxidative stress zone(cluster 5) showed elvated NFE2L2 with strong enrichment of G2M checkpoint, glycolysis, IL2-STAT5, interferon-a pathways representing a highly stressed B cell niche"),
                                            tags$li("Differentiation zone(cluster 4) showed higest IRF4 regulation activity and up regulation of CDC20, HLA-DRA, CD37 and CXCL13, markers of antigen representing. ")
                                        )
                                    )
                                )
                            )
                     ),
                     
                     # 
                     column(4,
                            br(), br(),
                            imageOutput("tf_edit_img", height = "auto"),
                            br(),
                            tags$em("TF-defined spatial niches in the human lymph node",
                                    style = "display: block; text-align: center; color: #2C3E50; font-size: 15px;")
                     )
                 ),
                 
                 hr(),
                 
                 # Footer
                 fluidRow(
                     column(12,
                            h5("Developed by Yasna Reihani"),
                            p("Built with R, Shiny."),
                            tags$a(href="https://github.com/Yasna81/lymph_node_tfspatial.git", "View Source Code on GitHub", target="_blank")
                     )
                 )
             )
    ),
    
    # ---------------------- 2. TF MAP ----------------------
    tabPanel("TF Map",
             sidebarLayout(
                 sidebarPanel(
                     selectInput("tf1", "Select TF (left):", choices = c("FOSB", "IRF4", "ATF3", "MAFF", "NFE2L2"), 
                                 selected = "FOSB"),
                     downloadButton("dl_tf1", "Download TF1 Plot"),
                     br(), br(),
                     selectInput("tf2", "Select TF (right):", choices = c("FOSB", "IRF4", "ATF3", "MAFF", "NFE2L2"), 
                                 selected = "IRF4"),
                     downloadButton("dl_tf2", "Download TF2 Plot")
                 ),
                 mainPanel(
                     fluidRow(
                         column(6, imageOutput("tf1_plot")),
                         column(6, imageOutput("tf2_plot"))
                     ),
    
                 br() , br() , br() , br() ,
                 hr(style = "border-top: 2px solid #18BC9C; margin: 40px 0;"),
                 br(), br(),
                 
                 fluidRow(
                     column(10, offset = 1,
                            div(style = "
                                background: linear-gradient(135deg, #f8f9fa 0%, #e9f7f5 100%);
                                border-left: 6px solid #18BC9C;
                                border-radius: 12px;
                                padding: 25px;
                                box-shadow: 0 6px 16px rgba(0,0,0,0.1);
                                font-size: 16px;
                                line-height: 1.8;
                                ",
                                h4("Interpretation of Transcription Factor Activity Maps", 
                                   style = "color: #2C3E50; margin-top: 0; font-weight: bold;"),
                                
                                tags$ul(style = "margin-bottom: 20px;",
                                        tags$li(tags$strong("IRF4"), " – Master regulator of B-cell differentiation and plasma cell fate. High in germinal centers and dark zones."),
                                        tags$li(tags$strong("FOSB & ATF3"), " – Immediate-early stress/response genes (AP-1 family). Mark sites of acute activation, especially dendritic cells and activated T cells."),
                                        tags$li(tags$strong("NFE2L2 (Nrf2)"), " – Central regulator of oxidative stress and metabolism. Enriched in long-lived memory cells and stressed B-cell zones."),
                                        tags$li(tags$strong("MAFF"), " – Often co-expressed with NFE2L2; involved in stress response and redox balance.")
                                ),
                                tags$p(
                                    tags$strong("Spatial Patterns: "),
                                    "These TFs show distinct and often opposing spatial distributions: ",
                                    tags$em("ATF3 and FOSB"), " are co-localized (acute activation zones), while ",
                                    tags$em("IRF4"), " and ", tags$em("NFE2L2"), " display anticorrelated patterns (differentiation vs. stress niches). ",
                                    tags$em("MAFF"), " partially overlaps with NFE2L2 and complements the ATF3/FOSB zones."
                                ),
                                
                                
                                
                                tags$p(
                                    tags$strong("Key Insight: "), 
                                    "These five TFs capture three major functional states in the lymph node: ",
                                    tags$em("differentiation (IRF4)"), ", ",
                                    tags$em("acute activation (ATF3/FOSB)"), ", and ",
                                    tags$em("chronic stress/metabolic adaptation (NFE2L2/MAFF)"),
                                    ". Their spatial patterns reveal how immune programs are organized in tissue."
                                )
                            )
                     )
                 ),
                 
                
             )
             )
    ),
    
    # ---------------------- 3. CLUSTERS ----------------------
    tabPanel("Clusters",
             fluidPage(
                 h3("TF Clusters and Cell Type Maps"),
                 fluidRow(
                     column(
                         6,
                         h4("TF Cluster Map"),
                         downloadButton("dl_tf_cluster", "Download TF Cluster Plot"),
                         br(), br(),
                         imageOutput("tf_cluster_img")
                     ),
                     column(
                         6,
                         h4("Cell Type Map"),
                         downloadButton("dl_celltype_cluster", "Download Cell Type Plot"),
                         br(), br(),
                         imageOutput("celltype_cluster_img")
                     )
                 ),
                 
                 br() ,br(),br(),br(),br(), 
                 
                 # 
                 div(style = "background-color: #f8fff8;
                border-left: 6px solid #2ecc71;
                padding: 18px 22px;
                margin: 200px 0 40px 0;
                border-radius: 10px;
                box-shadow: 0 3px 10px rgba(0,0,0,0.08);",
                     h4(style = "margin-top: 0; color: #1a1a1a;",
                        strong("Key Insight: Each TF Cluster is Dominated by a Specific Immune Cell Type")),
                     p("7 major TF clusters were found. ",
                       "transfering cell types from single cell data reveals that ", strong("every cluster is predominantly composed of one cell type"), ":"),
                     tags$ul(style = "font-size: 15px; line-height: 1.8;",
                             tags$li(strong("Cluster 0 & 1"), " → Dendritic Cells (DC)"),
                             tags$li(strong("Cluster 2"), " → Dendritic Cells (DC)"),
                             tags$li(strong("Cluster 3"), " → Memory CD4 T Cells"),
                             tags$li(strong("Cluster 4"), " → B Cells"),
                             tags$li(strong("Cluster 5"), " → B Cells"),
                             tags$li(strong("Cluster 6"), " → Memory CD4 T Cells"),
                             tags$li(strong("Cluster 7"), " → Dendritic Cells (DC)")
                     ),
                 ),
                 # ===  TF Score per Cluster Heatmap ===
                 br(), br(), hr(), br(),
                 
                 fluidRow(
                     column(
                         width = 6,
                         h4("TF Activity Score per Cluster Heatmap"),
                         downloadButton("dl_tf_score_heatmap", "Download TF Score Heatmap"),
                         br(), br(),
                         imageOutput("tf_score_per_cluster_heatmap", height = "auto")
                     ),
                     column(
                         width = 6,
                         h4("enriched pathways for each cluster"),
                         downloadButton("dl_dotplot", "Download dotplot"),
                         br(), br(),
                         imageOutput("dotplot", height = "auto")
                     )
                 ),
                   br(),br(), 
                 
                 
                 div(style = "background-color: #f8fff8;
                border-left: 6px solid #2ecc71;
                padding: 18px 22px;
                margin: 200px 0 40px 0;
                border-radius: 10px;
                box-shadow: 0 3px 10px rgba(0,0,0,0.08);",
                     h4(style = "margin-top: 0; color: #1a1a1a;",
                        strong("Key Insight: Each TF Cluster is Dominated by a Specific transcription factor")),
                     p("7 major TF clusters were found. ",
                       "transfering TF activity scores from single cell data reveals that ", strong("every cluster is predominantly enriched in one transcription factor"), ":"),
                     tags$ul(style = "font-size: 15px; line-height: 1.8;",
                             tags$li(strong("Cluster 0"), " → IRF4"),
                             tags$li(strong("Cluster 1 & 2"), " → ATF3/FOSB"),
                             tags$li(strong("Cluster 3"), " → MAFF"),
                             tags$li(strong("Cluster 4"), " → IRF4"),
                             tags$li(strong("Cluster 5"), " → NFE2L2"),
                             tags$li(strong("Cluster 6"), " → NFE2L2"),
                             tags$li(strong("Cluster 7"), " → MAFF")
                     ),
                 ),
                 
             )
    ),
    
    # ---------------------- 4. CORRELATION PLOT ----------------------
    tabPanel("Correlation Plot",
             sidebarLayout(
                 sidebarPanel(
                     selectizeInput("corr_tfs", "Select TFs:", choices = tf_list,
                                    multiple = TRUE,
                                    selected = c("IRF4","FOSB","IRF4","ATF3","NFE2L2"),
                                    options = list(placeholder = 'Choose TFs…')),
                     downloadButton("dl_corr", "Download Correlation Plot")
                 ),
                 mainPanel(
                         #plotOutput("cor_plot"),
                         # This shows either the plot OR a helpful message
                         uiOutput("correlation_plot_or_message")
                     
                 )
             )
    )
)

server <- function(input, output, session) {
    
    # -------- TF MAP PLOTS --------
    output$tf1_plot <- renderImage({
        req(input$tf1)
        list(
            src = tf_files[[input$tf1]],
            contentType = 'image/png',
            width = "100%"
        )
    }, deleteFile = FALSE)
    
    output$tf2_plot <- renderImage({
        req(input$tf2)
        list(
            src = tf_files[[input$tf2]],
            contentType = 'image/png',
            width = "100%"
        )
    }, deleteFile = FALSE)
    
    # -------- DOWNLOAD: TF1 --------
    output$dl_tf1 <- downloadHandler(
        filename = function() paste0(input$tf1, "_spatial.png"),
        content = function(file) {
            file.copy(tf_files[[input$tf1]], file)
        }
    )
    
    # -------- DOWNLOAD: TF2 --------
    output$dl_tf2 <- downloadHandler(
        filename = function() paste0(input$tf2, "_spatial.png"),
        content = function(file) {
            file.copy(tf_files[[input$tf2]], file)
        }
    )
    
    # -------- STATIC CLUSTER PLOTS --------
    output$tf_cluster_img <- renderImage({
        list(src = "tf_cluster_img.png",
             contentType = "image/png",
             width = "85%",
             height = "auto",
             style = " display: block; margin: 20px auto; border-radius: 10px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);")
    }, deleteFile = FALSE)
    
    output$celltype_cluster_img <- renderImage({
        list(src = "celltype_cluster_img.png",
             contentType = "image/png",
             width = "85%",
             height = "auto",
             style = "display: block; margin: 20px auto; border-radius: 10px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);")
    }, deleteFile = FALSE)
    
    # -------- DOWNLOAD: TF CLUSTER --------
    output$dl_tf_cluster <- downloadHandler(
        filename = "tf_cluster_img.png",
        content = function(file) {
            file.copy("tf_cluster_img.png", file)
        }
    )
    
    # -------- DOWNLOAD: CELL TYPE CLUSTER --------
    output$dl_celltype_cluster <- downloadHandler(
        filename = "celltype_cluster_img.png",
        content = function(file) {
            file.copy("celltype_cluster_img.png", file)
        }
    )
    
    # -------- CORRELATION PLOT --------
    # Main output
    output$correlation_plot_or_message <- renderUI({
        selected <- input$corr_tfs
        
        if (is.null(selected) || length(selected) < 2) {
            # 
            tags$div(
                style = "text-align: center; padding: 130px 20px; color: #555; 
               padding: 100px 20px; 
               color: #666; 
               font-size: 20px; 
               background-color: #f8f9fa; 
               border: 2px dashed #ccc; 
               border-radius: 12px; 
               margin: 40px;",
                tags$p(tags$strong("No TFs selected yet")),
                tags$p("Please choose at least 2 transcription factors from the left panel",
                       br(),
                       "to display the correlation heatmap."),
                tags$br(),
                tags$i(class = "fa fa-arrow-left", style = "font-size: 28px;")
            )
        } else {
            tagList(
            # Show the actual plot when enough TFs are selected
            plotOutput("cor_plot", height = "750px"),
            br(), br(),
            
            # Beautiful explanation box
            div(style = "background-color: #f8fff8;
                  border-left: 6px solid #18BC9C;
                  padding: 22px 28px;
                  margin: 40px 20px;
                  border-radius: 12px;
                  box-shadow: 0 4px 15px rgba(0,0,0,0.08);
                  font-size: 16px; line-height: 1.8;",
                h4(style = "margin-top: 0; color: #2C3E50;", strong("What does this correlation heatmap show?")),
                tags$ul(
                    tags$li(tags$strong("Blue (+1):"), " Two TFs are highly co-active in the same spatial spots and they likely belong to the same regulatory program."),
                    tags$li(tags$strong("Red (-1):"), " TFs have opposing spatial patterns and they mark mutually exclusive niches (e.g., differentiation vs. stress)."),
                    tags$li(tags$strong("White (~0):"), " No consistent spatial relationship.")
                ),
                p(tags$strong("Biological insight:"), 
                  "Strong positive clusters reveal coordinated TF modules (ATF3 + FOSB = acute activation), while negative correlations highlight functional antagonism in the tissue (IRF4 and NFE2L2), confirming spatial patterns.")
            )
            )
        }
    })
    
    output$cor_plot <- renderPlot({
        req(input$corr_tfs, length(input$corr_tfs) >= 2)   # this line is crucial
        
        tfs <- input$corr_tfs
        mat <- cor_mat[, tfs]   # 
        
        corrplot::corrplot(
            cor(mat, use = "pairwise.complete.obs"),
            method = "color",
            type = "upper",
            order = "hclust",
            tl.cex = 1.1,
            tl.col = "black",
            title = paste(length(tfs), "selected TFs – correlation")
        )
    })
    # -------- DOWNLOAD: CORRELATION PLOT --------
    output$dl_corr <- downloadHandler(
        filename = "correlation_plot.png",
        content = function(file) {
            png(file, width = 1400, height = 1200, res = 150)
            sel <- input$corr_tfs
            submat <- cor_mat[sel, sel, drop = FALSE]
            corrplot(submat, method = "color", type = "upper", tl.cex = 1.1)
            dev.off()
        })
    # Render the tf_clus.png image
    output$tf_edit_img <- renderImage({
        list(src = "tf_edit_img.png",
             contentType = "image/png",
             width = "100%")
    }, deleteFile = FALSE)
    
    output$tf_score_per_cluster_heatmap <- renderImage({
        list(src = "tf_score_per_cluster_heatmap.png",   
             contentType = "image/png",
             width = "60%",
             height = "auto",
             style = "display: block; margin: 20px auto; border-radius: 10px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);",
             alt = "TF score per cluster heatmap")
    }, deleteFile = FALSE)
    
    # Download handler
    output$dl_tf_score_heatmap <- downloadHandler(
        filename = function() { "tf_score_per_cluster_heatmap.png" },
        content = function(file) {
            file.copy("tf_score_per_cluster_heatmap.png", file)
        }
    )
    output$dotplot <- renderImage({
        list(src = "dotplot.png",   
             contentType = "image/png",
             width = "60%",
             height = "auto",
             style = "display: block; margin: 20px auto; border-radius: 10px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);",
             alt = "GSEA in each TF cluster")
    }, deleteFile = FALSE)
    
    # Download handler
    output$dl_dotplot <- downloadHandler(
        filename = function() { "dotplot.png" },
        content = function(file) {
            file.copy("dotplot.png", file)
        }
    )
    
}

shinyApp(ui, server)