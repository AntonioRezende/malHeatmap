library(shiny)
library(shinyFiles)
library(bslib)
library(vcfR)
library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

server <- function(input, output, session) {
###############################DIR VCFs#########################  
  volumes<-c(getVolumes()())
  
  shinyDirChoose(
    input,
    'inputdata',
    roots = volumes, session = session, allowDirCreate = FALSE
  )
  
  global <- reactiveValues(datapath = getwd())
  out <- reactive(input$inputdata)
  output$out <- renderText({
    global$datapath
  })
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$inputdata
               },
               handlerExpr = {
                 if (!"path" %in% names(out())) return()
                 home2 <- normalizePath("/")
                 global$datapath <-
                   file.path(home2, paste(unlist(out()$path[-1]), collapse = .Platform$file.sep))
               })
####################################################################
  ##REading VCFs

  vcf_files<-reactive({
    mypaht<-global$datapath
    vcf_file <- list.files(path = mypaht, pattern = "*.vcf", full.names = TRUE)
  })

  
  extract_af_from_gt <- function(vcf_in,mol) {
    
    bed<-read.table("Pfalciparum.genes.bed", header = FALSE, sep = "\t")
    colnames(bed) <- c("chromosome","start","end","gene")
    
    vcf <- read.vcfR(vcf_in)
    
    # Get variant positions (loci) from CHROM, POS fields
    chrom <- vcf@fix[, "CHROM"]
    pos <- vcf@fix[, "POS"]
    loci <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_")
    
    # Extract AF values from GT field
    af <- extract.gt(vcf, element = "AF", as.numeric = TRUE)
    
    #Extract ANN (annotation) from INFO field
    ann<-extract.info(vcf, element = "ANN")
    
    # Create data frame with loci and AF
    #af_data <- data.frame(locus = loci, AF = as.numeric(af), stringsAsFactors = FALSE)
    #return(af_data)
    #temp <- data.frame(chromosome = chrom, locus = loci, af = as.numeric(af), stringsAsFactors = FALSE)
    
    if(mol == "AA"){
      temp <- data.frame(chromosome = chrom, locus = loci, af = as.numeric(af), ann = ann, stringsAsFactors = FALSE)
      
      temp <- temp %>%
        filter(grepl("missense_variant",ann))
      
      temp$ann<-sub(".*\\|p\\.(\\w+\\d+\\w+)\\|.*", "\\1", temp$ann)
      
    
    
    
      final_data<-data.frame()
    
      for(i in seq_along(temp$chromosome)){
        #teste[i,]
        crome <- temp$chromosome[i]
        loc = as.numeric(gsub("Pf.+_","",temp$locus[i],perl = TRUE))
        specficbed<-bed[bed$chromosome==crome,]
        for (j in seq_along(specficbed$chromosome)){
          if (loc>=specficbed$start[j] & loc<=specficbed$end[j]) {
            gene<-specficbed$gene[j]
            #print(gene)
          }
        }
        #final_data<-rbind(final_data,data.frame(chromosome = crome, locus = loc, gene = gene , af = temp$af[i] ,stringsAsFactors = FALSE))
        final_data<-rbind(final_data,data.frame(chromosome = crome, locus = loc, gene = gene , ann = temp$ann[i], af = temp$af[i], stringsAsFactors = FALSE))
      }
      return(final_data)
      
    }else if(mol=="NT"){
      temp <- data.frame(chromosome = chrom, locus = loci, af = as.numeric(af), stringsAsFactors = FALSE)
      
      
      final_data<-data.frame()
      
      for(i in seq_along(temp$chromosome)){
        #teste[i,]
        crome <- temp$chromosome[i]
        loc = as.numeric(gsub("Pf.+_","",temp$locus[i],perl = TRUE))
        specficbed<-bed[bed$chromosome==crome,]
        for (j in seq_along(specficbed$chromosome)){
          if (loc>=specficbed$start[j] & loc<=specficbed$end[j]) {
            gene<-specficbed$gene[j]
            #print(gene)
          }
        }
        final_data<-rbind(final_data,data.frame(chromosome = crome, locus = loc, gene = gene , af = temp$af[i] ,stringsAsFactors = FALSE))
        #final_data<-rbind(final_data,data.frame(chromosome = crome, locus = loc, gene = gene , ann = temp$ann[i], af = temp$af[i], stringsAsFactors = FALSE))
      }
      return(final_data)
    
    }
  }
  
  vcf_processed<-reactive({
    infiles<-vcf_files()
    typemole <-as.character(input$typemol)
    af_list <- lapply(infiles,extract_af_from_gt,typemole)
    
    if(input$typemol=="NT"){
      af_merged <- Reduce(function(x, y, z) {
        full_join(x, y, z, by = c("chromosome", "locus", "gene"))
      }, af_list)
    
      colnames(af_merged)[-c(1, 2, 3)] <- gsub(".SNP.annot.vcf","",basename(infiles))
    
    }else if(input$typemol=="AA"){
      af_merged <- Reduce(function(x, y, z ) {
        full_join(x, y, z,  by = c("chromosome", "locus", "gene", "ann"))
      }, af_list)
      
      colnames(af_merged)[-c(1, 2, 3, 4)] <- gsub(".SNP.annot.vcf","",basename(infiles))
      
      
    }
    
    
    af_merged[is.na(af_merged)] = 0 
    
    
    unique_chromosomes <- unique(af_merged$gene)
    
    return(list(af_merged,unique_chromosomes))
    
  })
  
  
  output$file <-
    renderUI(expr = selectInput(inputId = 'geneid',
                                label = 'Gene ID',
                                choices = vcf_processed()[[2]]))
  
  
  output$plot <- renderPlot({
    chrom<-input$geneid
    af_merged <- vcf_processed()[[1]]
    chrom_data <- af_merged %>%
      filter(gene == chrom)
    
    if(input$typemol=="NT"){
    
    chrom_data<-chrom_data[order(chrom_data$locus),]
    chrom_data <- chrom_data[(rowSums(chrom_data[,-c(1:3) ]))/length(chrom_data[1,-c(1:3)]) >= input$af, ]
    pheatmap(t(as.matrix(chrom_data[,-c(1,2,3)])),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             show_colnames = TRUE,
             labels_col = chrom_data$locus,
             angle_col = 45,
             #color = colorRampPalette(c("blue", "white", "red"))(50),
             color = colorRampPalette((brewer.pal(9,"Spectral")))(50),
             main = paste("Allele Frequency Heatmap for Gene", chrom))
    
    }else if(input$typemol=="AA"){
      chrom_data<-chrom_data[order(chrom_data$locus),]
      chrom_data <- chrom_data[(rowSums(chrom_data[,-c(1:4) ]))/length(chrom_data[1,-c(1:4)]) >= input$af, ]
      pheatmap(t(as.matrix(chrom_data[,-c(1,2,3,4)])),
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               show_colnames = TRUE,
               labels_col = chrom_data$ann,
               angle_col = 45,
               #color = colorRampPalette(c("blue", "white", "red"))(50),
               color = colorRampPalette((brewer.pal(9,"Spectral")))(50),
               main = paste("Allele Frequency Heatmap for Gene", chrom))
      
    }
    
  })
  
  output$comando <- renderText({
    input$typemol
    #global$datapath
  })
  
  
  
  
}
