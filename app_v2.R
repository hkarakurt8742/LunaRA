library(shiny)
library(shinybusy)
library(DT)
library(ggnet)
library(plotly)
library(clusterProfiler)
library(VariantAnnotation)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(org.Hs.eg.db)
library(myvariant)
library(stringr)
library(shinycssloaders)
library(readr)
library(ggraph)
library(igraph)

# Function to flatten list to character
flatten_list_to_character <- function(input_list) {
  sapply(input_list, function(x) {
    if (is.null(x)) {
      return(NA_character_)
    } else {
      return(paste(x, collapse = ", "))
    }
  })
}

# Define UI for the app
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .bottom-left-logos {
        position: fixed;
        bottom: 25px;
        left: 25px;
        z-index: 1000;
      }
    "))
  ),
  titlePanel(
    div(
      img(src = "idea.png", height = "100px", width = "100px"), 
      "LunaRa"
    )
  ),
  add_busy_gif(src = "https://media2.giphy.com/media/7aqTC87afdWGXF3o8T/giphy.gif?cid=ecf05e475vgr4jui6uha1neqjhfbpsqwsa6fwklob3kzcns1&ep=v1_gifs_related&rid=giphy.gif&ct=s",
               height = 150, width = 150, position = "full-page"),
  sidebarLayout(
    sidebarPanel(
      selectInput("genome", "Select Genome Assembly:", 
                  choices = c("hg19", "hg38"), selected = "hg19"),
      fileInput("vcfFile", "Upload VCF File", accept = c(".vcf")),
      actionButton("process", "Process VCF File"),
      downloadButton("downloadData", "Download Variant Data"),
      selectInput("filter_gene", "Filter by Gene:", 
                  choices = NULL, selected = NULL, multiple = TRUE),
      br(),
      h4("Information"),
      p("LunaRa is a light-weight web application to annotate genomic variants, find pharmacogenomic interactions, and construct a protein-protein interaction network based on genes with variants and their first neighbours."),
      p("Developed by Idea Technology Solutions.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Variant Annotation", DTOutput("variant_table") %>% withSpinner()),
        tabPanel("Protein-Protein Interaction", 
                 plotlyOutput("ppin_network_plot"),
                 plotlyOutput("ppin_go_bp_plot"),
                 plotlyOutput("ppin_kegg_plot")),
        tabPanel("Pharmacogenomic Interactions", uiOutput("pharma_network_plot")),
        tabPanel("Disease Associations", 
                 DTOutput("disease_table") %>% withSpinner(), 
                 plotlyOutput("disease_network_plot")),
        tabPanel("Help", 
                 h4("Help Section"),
                 p("Use this app to analyze genomic variants. Upload a VCF file, select the genome assembly, and click 'Process VCF File'."), 
                 p("You can filter the results by gene and download the data for further analysis."))
      )
    )
  ),
  div(
    class = "bottom-left-logos",
    img(src = "HumanNet.png", height = "50px", width = "250px", style = "margin-right: 10px;"),
    img(src = "myvariant.png", height = "100px", width = "250px")
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$process, {
    req(input$vcfFile)
    
    vcf_path <- input$vcfFile$datapath
    
    # VCF reading with error handling
    vcf <- tryCatch({
      readVcf(vcf_path, genome = input$genome)
    }, error = function(e) {
      showNotification("Error reading VCF file. Please check the file format and contents.", type = "error")
      return(NULL)
    })
    
    if (is.null(vcf)) return()  
    
    txdb <- if (input$genome == "hg19") {
      TxDb.Hsapiens.UCSC.hg19.knownGene
    } else {
      TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    
    hgvs <- formatHgvs(vcf)
    all_variants <- getVariants(hgvs, fields = "all")
    all <- as.data.frame(all_variants)
    
    chromosomes <- seqnames(vcf)
    positions <- start(vcf)
    
    variant_ranges <- GRanges(seqnames = chromosomes, ranges = IRanges(start = positions, end = positions))
    gene_info <- genes(txdb)
    overlaps <- findOverlaps(variant_ranges, gene_info)
    gene_names <- gene_info$gene_id[subjectHits(overlaps)]
    
    result <- data.frame(
      Chromosome = as.character(chromosomes[queryHits(overlaps)]),
      Position = positions[queryHits(overlaps)],
      Gene = gene_names
    )
    
    gene_symbols <- mapIds(org.Hs.eg.db, keys = result$Gene, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    result$GeneSymbol <- gene_symbols
    
    result_combined <- result %>%
      group_by(Chromosome, Position) %>%
      summarise(Gene = first(Gene), 
                GeneSymbols = paste(unique(GeneSymbol), collapse = "; "),
                .groups = 'drop')
    
    query_data <- data.frame(query = all$query, stringsAsFactors = FALSE) %>%
      mutate(
        Chromosome_Number = as.character(str_extract(query, "(?<=chr)\\w+")),
        Position = as.numeric(str_extract(query, "(?<=g\\.)[0-9]+")),
        Change = case_when(
          str_detect(query, "ins") ~ "ins",
          str_detect(query, "del") ~ "del",
          TRUE ~ str_extract(query, "[A-Z]>[A-Z]")
        ),
        Change_Type = case_when(
          str_detect(query, "ins") ~ "ins",
          str_detect(query, "del") ~ "del",
          str_detect(query, "[A-Z]>[A-Z]") ~ "SNV",
          TRUE ~ NA_character_
        )
      )
    
    all$Chromosome_Number <- query_data$Chromosome_Number
    all$Position <- query_data$Position
    all$Change <- query_data$Change
    all$Change_Type <- query_data$Change_Type
    
    result_combined <- result_combined %>%
      mutate(Chromosome = gsub("chr", "", Chromosome))
    
    df1 <- as.data.frame(result_combined[, 1:4])
    df2 <- as.data.frame(all)
    
    all <- df2 %>%
      left_join(df1 %>% dplyr::select(Chromosome, Position, Gene, GeneSymbols), 
                by = c("Chromosome_Number" = "Chromosome", "Position" = "Position"))
    
    df <- all[is.na(all$notfound), ]
    
    new_df <- data.frame(
      variation = df$query,
      Gene_Symbol = df$GeneSymbols,
      Change_Type = df$Change_Type,
      gnomad_exome_rsid = df$gnomad_exome.rsid,
      gnomad_genome_rsid = df$gnomad_genome.rsid,
      dpsnp_resid = df$dbnsfp.eqtlgen.snp_id,
      clinvar_rsid = df$clinvar.rsid,
      dbsnp_rsid = df$dbsnp.rsid,
      stringsAsFactors = FALSE
    )
    new_df <- new_df %>%
      mutate(
        cadd_annotype = flatten_list_to_character(df$cadd.annotype),  
        cadd_consequence = flatten_list_to_character(df$cadd.consequence),
        clinvar_hgvs_coding = flatten_list_to_character(df$clinvar.hgvs.coding),
        dpsnp_transcript_id = flatten_list_to_character(df$dbnsfp.ensembl.transcriptid),
        exac_af = flatten_list_to_character(df$exac.af),
        gnomad_exome_af_af = flatten_list_to_character(df$gnomad_exome.af.af),
        exac_polyphen = flatten_list_to_character(df$cadd.polyphen.cat),
        clinvar_classification = flatten_list_to_character(df$dbnsfp.clinvar.clnsig)
      ) %>%
      mutate(rsid = coalesce(gnomad_exome_rsid, 
                             gnomad_genome_rsid, 
                             clinvar_rsid, 
                             dbsnp_rsid)) %>%
      dplyr::select(-gnomad_exome_rsid, -gnomad_genome_rsid, -clinvar_rsid, -dbsnp_rsid, -dpsnp_resid)
    
    
    updateSelectInput(session, "filter_gene", choices = unique(new_df$Gene_Symbol))
    
    output$variant_table <- renderDT({
      filtered_df <- new_df
      if (!is.null(input$filter_gene) && length(input$filter_gene) > 0) {
        filtered_df <- filtered_df %>% filter(Gene_Symbol %in% input$filter_gene)
      }
      datatable(filtered_df, options = list(pageLength = 10))
    })
    
    # Protein-Protein Interaction
    humannet <- read.table("HumanNet-FN-Symbols.tsv", sep = "\t", header = F)
    subnetwork_gene <- new_df$Gene_Symbol
    
    sample_specific_network <- humannet %>%
      filter(V1 %in% unique(subnetwork_gene) | V2 %in% unique(subnetwork_gene))
    
    sample_specific_network <- na.omit(sample_specific_network)
    
    network <- network::network(sample_specific_network, matrix.type = "edgelist")
    
    ppin_network <- ggnet2(network, node.size = 5, label = TRUE, label.size = 2,
                           label.color = "black", node.color = "chocolate1")
    
    # Protein-Protein Plot
    output$ppin_network_plot <- renderPlotly({
      ggplotly(ppin_network)
    })
    
    # Enrichment Analysis
    network_enrichment <- enrichGO(unique(c(sample_specific_network$V1, sample_specific_network$V2)),
                                   keyType = "SYMBOL", ont = "BP", OrgDb = "org.Hs.eg.db")
    
    network_enrichment_kegg <- enrichKEGG(gene = unique(result_combined$Gene), organism = "hsa")
    
    output$ppin_go_bp_plot <- renderPlotly({
      go_plot <- dotplot(network_enrichment) +
        theme_minimal(base_size = 10) +
        theme(
          aspect.ratio = 1,
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          plot.title = element_text(size = 10, face = "bold"),
        )
      
      ggplotly(go_plot) %>%
        layout(width = 800, height = 400)  
    })
    
    output$ppin_kegg_plot <- renderPlotly({
      kegg_plot <- dotplot(network_enrichment_kegg) +
        theme_minimal(base_size = 10) +
        theme(
          aspect.ratio = 1,
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          plot.title = element_text(size = 10, face = "bold"),
          panel.grid.major = element_line(color = "lightgrey"),
          panel.grid.minor = element_blank()
        )
      
      ggplotly(kegg_plot) %>%
        layout(width = 800, height = 300)  
    })
    
    # Pharmacogenomic network
    pharmgkb <- read.table("PharmGKB_final_network.tsv", sep = "\t", header = TRUE)
    
    sample_pharma_network <- pharmgkb %>%
      filter(Node1 %in% new_df$rsid | Node2 %in% new_df$rsid)
    
    sample_pharma_network <- na.omit(sample_pharma_network)
    
    output$pharma_network_plot <- renderUI({
      if (nrow(sample_pharma_network) == 0) {
        tagList(
          h3("No genomic variants associated with any pharmacogenomic interactions found in the PharmGKB database.")
        )
      } else {
        pharma_network <- network::network(sample_pharma_network, matrix.type = "edgelist")
        pharma_network_plot <- ggnet2(pharma_network, node.size = 5, label = TRUE, label.size = 2,
                                      label.color = "black", node.color = "chocolate1")
        
        ggplotly(pharma_network_plot)
      }
    })
    
    # Disease Associations
    disease_data <- read.table("human_disease_knowledge_filtered.tsv", header=FALSE , sep = "\t")
    disease_data <- disease_data[, c(2, 4)]  
    colnames(disease_data) <- c("Gene_Symbol", "Disease")
    
    # Prepare filtered diseases
    observeEvent(input$process, {
      unique_genes <- unique(unlist(strsplit(na.omit(all$GeneSymbols), "; ")))
      filtered_diseases <- disease_data[disease_data$Gene_Symbol %in% unique_genes, ]
      filtered_diseases <- filtered_diseases %>% distinct(Gene_Symbol, Disease)
      
      # Render the disease associations table
      output$disease_table <- renderDT({
        datatable(filtered_diseases, options = list(pageLength = 10))
      })
      
      # Render the disease network plot
      filtered_diseases <- na.omit(filtered_diseases)
      filtered_diseases <- unique(filtered_diseases)
      
      
      network <- network::network(filtered_diseases, matrix.type = "edgelist")
      
      disease_network <- ggnet2(network, node.size = 5, label = TRUE, label.size = 2,
                             label.color = "black", node.color = "chocolate1")
      
      # Protein-Protein Plot
      output$disease_network_plot <- renderPlotly({
        ggplotly(disease_network)
      })
      
})

    
    # Download
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("variant_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(new_df, file, row.names = FALSE)
      }
    )
  })
}

shinyApp(ui = ui, server = server)
