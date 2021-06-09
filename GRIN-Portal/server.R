################################################
################# GRIN PORTAL ##################
################################################
#
# This is a Shiny web application #
#
#setwd("~/00_GRIN_Portal/App/GRIN_Portal")

################## LIBRARIES ################

library(shiny)
library(plotly)
library(readxl)
library(readr)
library(r3dmol)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(tidyverse)
library(seqinr)
library(bio3d)

############## FUNCTIONS, STYLE ############

numextract <- function(string) {
    as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

log_fun <- function(x) { 
    x = case_when(
        as.numeric(x) > 0 ~ log10(x), #log10 transformation to make the visuals clean 
        as.numeric(x) == 0 ~ 0,
        as.numeric(x) < 0 ~ log10(x * -1) * -1, 
        TRUE ~ -9999)
    return(x)
    
}

plotly_font <- list(
    family = "sans-serif",
    size = 11)

goodbye <- c("zoomIn2d", "zoomOut2d", "hoverClosestGeo", "hoverClosestGl2d",
             "hoverClosestCartesian","lasso2d","select2d","resetScale2d",
             "hoverCompareCartesian", "hoverClosestPie","toggleSpikelines")

# coloring for colorblindness 
lolliplot_fill_scheme <-
    c(
        "signal" = "#CC79A7",
        "terminal" = "#009E73",
        "S" = "#56B4E9",
        "M" = "#E69F00",
        "linker" = "lightgray",
        "GRIN"="lightgray",
        "no" = "#FFFFFF",
        "missense" = "#D55E00",
        "null" = "#0072B2"
    )



GRIN_colors <- c("GRIN1"="#6fbbf7", 
                 "GRIN2A"="#ee6c71", 
                 "GRIN2B"="#ffbc5a", 
                 "GRIN2D"="#bf73cc")

############## DATA ############


# GRIN Registry
GRIN_clin <- read_csv("data/GRIN_clin.csv") %>%  
  mutate(delay_degree = case_when(is.na(delay_degree) | delay_degree == "DD/ID of unknown degree" ~ "Unknown",
                                  TRUE ~ delay_degree)) %>% mutate(delay_degree = factor(delay_degree, 
  levels = c("no DD/ID", "mild DD/ID", "moderate DD/ID", "severe/profound DD/ID","Unknown"))) %>% 
  mutate(seizures = case_when(is.na(seizures) ~ "Unknown",  TRUE ~ seizures)) %>% 
  mutate(seizures = factor(seizures, levels = c("Yes","No","Unknown"))) %>%
  mutate(movement_disorder = case_when(is.na(movement_disorder) ~ "Unknown",
                                       TRUE ~ movement_disorder)) %>%
  mutate(movement_disorder = factor(movement_disorder, levels = c("Yes","No","Unknown")))%>%
    mutate(PMID = ifelse(PMID >0, paste0("<a target=_blank href=https://pubmed.ncbi.nlm.nih.gov/",
                                                      paste0(PMID),'>',"PubMed",'</a>'))) %>% filter(!is.na(ref_c)) %>% 
  filter(origin != "biparental") %>% mutate(ref_p = substr(variant_p, 4,6))

GRIN_disease <- GRIN_clin %>% select(gene, delay, seizures) %>% 
    mutate(disease_type = case_when(is.na(delay) | is.na(seizures) | seizures == "Unknown" ~ "Unknown",
                                    delay=="No" & seizures=="Yes" ~ "Seizures only",
                                    delay=="Yes" & seizures=="Yes" ~ "Seizures & DD/ID",
                                    delay=="Yes" & seizures=="No" ~ "DD/ID only",
                                    TRUE ~ "Unknown")) %>% add_count(gene, disease_type) %>% unique() 

# GRIN ACMG
GRIN_ACMG <- read_excel("data/GRIN_ACMG_variants_new.xlsx")%>% mutate(transcript = case_when(
    gene == "GRIN1" ~ "NM_007327.4",
    gene == "GRIN2A" ~ "NM_001134407.3",
    gene == "GRIN2B" ~ "NM_000834.4",
    gene == "GRIN2D" ~ "NM_000836.2"
))

# Functional
GRIN_functional <- read_csv("data/GRIN_functional_20210602.csv") %>% mutate(CALL = case_when(
  CALL == "Likely LoF (Supporting Pathogenic)" ~ "Likely Loss-of-Function",
  CALL == "Potential LoF (Supporting Pathogenic)" ~ "Potential Loss-of-Function",
  CALL == "Likely GoF (Supporting Pathogenic)" ~ "Likely Gain-of-Function",
  CALL == "Potential GoF (Supporting Pathogenic)" ~ "Potential Gain-of-Function",
  CALL == "No support for Pathogenicity" ~ "No Detectable Effect",
  TRUE ~ CALL
))

GRIN_functional_select <- GRIN_functional %>% select(gene, aa_pos, ref_p, alt_p, Fold_glu, Fold_gly, Fold_mg, Fold_ph, Fold_zn, Fold_znmin, `Fold Popen`,
                                         `Fold_tauW (ms)`, `Fold_Surface/Total`, `Fold_HEKPeakAmplitude pA/pF`, `Fold Synaptic Charge Transfer`,
                                         `Fold Non-Synaptic Charge Transfer`, CALL)

GRIN_functional_select_s <- GRIN_functional %>% select(gene, aa_pos, ref_p, alt_p, CALL)

GRIN_research <- merge(GRIN_clin, GRIN_functional_select, by=c("gene", "aa_pos", "ref_p", "alt_p"), all.x=TRUE, all.y=FALSE) 

GRIN_clin <- merge(GRIN_clin, GRIN_functional_select_s, by=c("gene", "aa_pos", "ref_p", "alt_p"), all.x=TRUE, all.y=FALSE) 

# animal studies
GRIN_animal <- read_excel("data/GRIN_animal.xlsx") %>% mutate(PMID=paste0("<a target=_blank href=https://pubmed.ncbi.nlm.nih.gov/",
                                                                               paste0(PMID),'>',"PubMed",'</a>'))

# GRIN domains
GRIN_domains <- read_excel("data/grin-domains-20210107.xlsx")

# GRIN gnomad
GRIN_gnomad <- read_csv("data/GRIN_gnomad.csv") %>% mutate(c_pos = numextract(c_pos)) %>% mutate(var_type = case_when(is.na(alt_p) ~ "null",
                                                                                                                  TRUE ~ "missense")) %>% filter(gene != "GRIN2C")
# For Variant Analysis
GRIN_var_analysis <- read_csv("data/GRIN_var_analysis3.csv") %>% 
    mutate(c_pos = numextract(variant_c)) %>% mutate(transcript = case_when(
        gene == "GRIN1" ~ "NM_007327.4",
        gene == "GRIN2A" ~ "NM_001134407.3",
        gene == "GRIN2B" ~ "NM_000834.4",
        gene == "GRIN2D" ~ "NM_000836.2"
    )) %>% mutate(CALL = case_when(
      CALL == "Likely LoF (Supporting Pathogenic)" ~ "Likely Loss-of-Function",
      CALL == "Potential LoF (Supporting Pathogenic)" ~ "Potential Loss-of-Function",
      CALL == "Likely GoF (Supporting Pathogenic)" ~ "Likely Gain-of-Function",
      CALL == "Potential GoF (Supporting Pathogenic)" ~ "Potential Gain-of-Function",
      CALL == "No support for Pathogenicity" ~ "No Detectable Effect",
      TRUE ~ CALL))


GRIN_pie_fun <- function(select_gene, select_colors){
    plot_ly(GRIN_disease %>% 
                select(gene, disease_type, n) %>% 
                unique() %>% 
                filter(gene==select_gene), 
            values=~n, labels=~disease_type, type="pie",
            textposition="inside", 
            textinfo="label+percent",
            insidetextfont = list(color = '#333333'),
            hoverinfo = "text",
            text= ~ paste(n, "individuals"),
            marker=list(colors=select_colors,
                        line = list(color = '#FFFFFF', width = 1)), showlegend = FALSE) %>% 
        layout(title=paste(""),font=plotly_font) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE) 
}

GRIN_bar_fun <- function(select_gene, select_pt_new, select_title){
    plot_ly(GRIN_clin %>% rename(select_pt = select_pt_new) %>% filter(gene==select_gene) %>% 
                select(gene, select_pt)  %>%
                add_count(gene, name = "n_gene") %>% add_count(gene, select_pt) %>% distinct(),
            x = ~ select_pt, 
            y = ~ round(n/n_gene, digits = 2), 
            colors = "Set2", 
            color = ~ select_pt, 
            type = "bar", 
            hoverinfo = "text", showlegend = FALSE,
            text= ~ paste0(round(n/n_gene, digits = 2), " (" ,n," individuals)")) %>% 
        layout(title="", 
               font=plotly_font,
               xaxis = list(title=""),
               yaxis = list(title="Proportion of Cohort")) %>%
        config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE) 
}

GRIN_functional_fun <- function(thing, select_function_new) {
    
        thing %>% filter(!is.na(CALL)) %>%
        rename(select_function = select_function_new) %>% group_by(gene) %>% 
        distinct(variant_c, .keep_all = TRUE) %>% ungroup() %>%
        # mutate(select_function = case_when(select_function > 5 ~ 6, select_function < -5 ~ -6, TRUE ~ select_function)) %>%
        # mutate(select_function = case_when(select_function < 0 ~ 1, TRUE ~ select_function)) %>%
    
        mutate(select_function = ifelse(select_function == -1000, min(as.numeric(select_function)) - 5, as.numeric(select_function)),
        select_function = ifelse(select_function > 0, log10(select_function),
                            # log10 transformation to make the visuals clean
                            ifelse(select_function == 0, 0, ifelse(select_function < 0, (log10(select_function * -1) * -1), -1000)))) %>%
    
        filter(!is.na(select_function)) %>% filter(!is.na(delay)) %>%
        mutate(select_function = round(select_function, digits = 2)) %>%
        
        plot_ly(y =  ~ select_function,
                x =  ~ gene,
                type = "violin",
                split =  ~ gene,
                points = "all",
                pointpos = 0,
                hoverinfo = "text",
                text =  ~ paste0(
                    "Variant: ",
                    variant_p,
                    " ,",
                    delay_degree, 
                    ", Potency = ",
                    select_function
                ),
                color =  ~ gene,
                colors = GRIN_colors
        ) %>%
        layout(title = div(select_function_new," fold potency"),
               font = plotly_font,
               xaxis = list(title = ""),
               yaxis = list(title = "Fold Potency Change", range =c(-6, 6)
                            )) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
 
}

GRIN_onset_fun <- function(select_gene) {
    plot_ly(GRIN_clin %>% filter(gene == select_gene), y = ~onset_months/12, x = ~gene, type = "box",
                    boxpoints = "all", jitter = 0.3,
                    pointpos = 0, hoverinfo = "text", showlegend = FALSE,
                    text= ~paste0(round(onset_months/12, digits = 2), " years, ", variant_p)) %>% layout(font=plotly_font,  title="",xaxis = list(title=""),
                                                                                    yaxis = list(title="Seizure onset (years)")) %>%
        config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
}
  
  r3dmol_patho_fun <- function(select_gene, select_data){
    select_data %>% 
      filter(gene == select_gene) %>%  filter(var_type == "missense") %>% 
      mutate(aa_ref = str_sub(variant_p,4,6), label = "pathogenic") %>% 
      group_by(aa_pos,aa_ref,gene) %>% summarise(n_occ = n()) %>% 
      select(aa_pos,aa_ref,n_occ,gene) 
  }
  
  r3dmol_gnomad_fun <- function(select_gene){
      GRIN_gnomad %>% 
      filter(gene == select_gene) %>%
      mutate(aa_ref = ref_p) %>% 
      group_by(aa_pos,aa_ref,gene) %>% 
      summarise(n_occ = n()) %>% 
      select(aa_pos,aa_ref,n_occ,gene) 
  }

  




############### SERVER ###############
shinyServer(function(input, output, session) {

    ############### PANEL DECISION ###############
    observeEvent(input$infoBtn, {
        updateTabsetPanel(session, "TabDisplay", selected = "infoTab")
    })
    
    observeEvent(input$familyBtn, {
        updateTabsetPanel(session, "TabDisplay", selected = "familyTab")
    })
    
    observeEvent(input$variantBtn, {
        updateTabsetPanel(session, "TabDisplay", selected = "variantTab")
    })
    
    observeEvent(input$researchBtn, {
        updateTabsetPanel(session, "TabDisplay", selected = "researchTab")
    })
    
    observeEvent(input$registryBtn, {
        updateTabsetPanel(session, "TabDisplay", selected = "registryTab")
    })
    
    ##### BASIC INFORMATION #####
    
    ### GRIN1
    
    output$GRIN1_pie <- renderPlotly({
        GRIN_pie_fun("GRIN1", c('#6fbbf7', '#8bcb8e', '#e5e5e5')) 
    })
    
    output$GRIN1_delay <- renderPlotly({
        GRIN_bar_fun("GRIN1", "delay_degree", "Developmental Delay/Intellectual Disability")
    })
    
    output$GRIN1_seizures <- renderPlotly({
        GRIN_bar_fun("GRIN1", "seizures", "Seizures")
    })
    
    output$GRIN1_seizure_onset <- renderPlotly({
        GRIN_onset_fun("GRIN1")
    })
    
    output$GRIN1_movement <- renderPlotly({
        GRIN_bar_fun("GRIN1", "movement_disorder", "Movement Disorder")
    })
    
    output$GRIN1_CVI <- renderPlotly({
        GRIN_bar_fun("GRIN1", "CVI", "Cerebral Visual Impairment")
    })
    
    output$GRIN1_MCD <- renderPlotly({
        GRIN_bar_fun("GRIN1", "MCD", "Malformation of Cortical Development")
    })
    
    ### GRIN2A
    output$GRIN2A_pie <- renderPlotly({
        GRIN_pie_fun("GRIN2A", c( '#8bcb8e', '#E69F00','#e5e5e5','#6fbbf7')) 
    })
    
    output$GRIN2A_delay <- renderPlotly({
        GRIN_bar_fun("GRIN2A", "delay_degree", "Developmental Delay/Intellectual Disability")
    })
    
    output$GRIN2A_seizures <- renderPlotly({
        GRIN_bar_fun("GRIN2A", "seizures", "Seizures")
    })
    
    output$GRIN2A_seizure_onset <- renderPlotly({
      
     plot<- plot_ly(GRIN_clin %>% filter(gene == "GRIN2A"), y = ~onset_months/12, x = ~gene, type = "box",
              boxpoints = "all", jitter = 0.3,
              pointpos = 0, hoverinfo = "text", showlegend = FALSE,
              text= ~paste0(round(onset_months/12, digits=2), " years, ", variant_p))
     
       
    plot <- plot %>% layout(font=plotly_font, title="",  xaxis = list(title=""), yaxis = list(title="Seizure onset (years)",
                                                                                   range = c(0, 25))) %>%
        config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    })
    
    output$GRIN2A_seizure_offset <- renderPlotly({
      
      plot<- plot_ly(x="GRIN2A", y=c(10.8, 12, 10, 20, 9, 8, 12, 16, 10), type = "box",
                     hoverinfo = "text",
                     text = c("p.(Arg847*) 10,8 years",	"p.(Trp198*) 12 years",	"p.(Pro415Hisfs*8) 10 years", 
                              "p.(Pro415Hisfs*8) 20 years", "c.1007+1G>A p.? 9 years",
                             "p.(Ala61Glyfs*78) 8 years","del exon 12-14 12 years",	"p.(Ala27Glyfs*112) 16 years", "del 10 years"),
                     boxpoints = "all", jitter = 0.3, 
                     color = "red",
                     pointpos = 0, showlegend = FALSE,
                    )
      
      
      
      plot <- plot %>% layout(font=plotly_font, title="", xaxis = list(title=""), yaxis = list(title="Seizure offset (years)",
              range = c(0, 25))) %>%
        config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
      
    })
    
    output$GRIN2A_movement <- renderPlotly({
        GRIN_bar_fun("GRIN2A", "movement_disorder", "Movement Disorder")
    })
    
    ### GRIN2B
    output$GRIN2B_pie <- renderPlotly({
        GRIN_pie_fun("GRIN2B", c('#6fbbf7', '#e5e5e5', '#8bcb8e')) 
    })
    
    output$GRIN2B_delay <- renderPlotly({
        GRIN_bar_fun("GRIN2B", "delay_degree", "Developmental Delay/Intellectual Disability")
    })
    
    output$GRIN2B_seizures <- renderPlotly({
        GRIN_bar_fun("GRIN2B", "seizures", "Seizures")
    })
    
    output$GRIN2B_seizure_onset <- renderPlotly({
        GRIN_onset_fun("GRIN2B")
    })
    
    output$GRIN2B_movement <- renderPlotly({
        GRIN_bar_fun("GRIN2B", "movement_disorder", "Movement Disorder")
    })
    
    ##### FOR FAMILIES #####
    # nothing to be calculated here # 
 
    ##### VARIANT ANALYSIS ##### 
    
    # choose either GRIN1, GRIN2A, GRIN2B or GRIN2D   
    varGeneInput <- reactive({
        input$var_gene
    })
    
    observeEvent(input$search_c_pos, {
      
      filtered_c_pos <- GRIN_var_analysis %>% filter(c_pos == input$search_c_pos, gene == input$var_gene)
      
      updatePickerInput(session = session, inputId = "search_ref_c",
                        choices = unique(filtered_c_pos$ref_c))
      
      updatePickerInput(session = session, inputId = "search_alt_c",
                        choices = unique(filtered_c_pos$alt_c))
      
    })
    
    observeEvent(input$search_aa_pos, {
      
      filtered_aa_pos <- GRIN_var_analysis %>% filter(aa_pos == input$search_aa_pos, gene == input$var_gene)
      
      updatePickerInput(session = session, inputId = "search_ref_p",
                        choices = sort(unique(filtered_aa_pos$ref_p)))
      
      updatePickerInput(session = session, inputId = "search_alt_p",
                        choices = unique(filtered_aa_pos$alt_p),
                        selected = "Ile")
      
    })
    

    
    varFilterInput <- reactiveValues(data=NULL)
    
    observeEvent(input$search_var_c, {
        varFilterInput$data <- GRIN_var_analysis %>% filter(gene==input$var_gene) %>% filter(c_pos==input$search_c_pos) %>%
            filter(ref_c==input$search_ref_c) %>% filter(alt_c==input$search_alt_c)
        
        print(varFilterInput$data)
    })
    
    observeEvent(input$search_var_p, {
        varFilterInput$data <- GRIN_var_analysis %>% filter(gene==input$var_gene) %>% filter(aa_pos==input$search_aa_pos) %>%
            filter(ref_p==input$search_ref_p) %>% filter(alt_p==input$search_alt_p)
        
        print(varFilterInput$data)
    })
    
    ACMGFilterInput <- reactiveValues(data2=NULL)
    
    observeEvent(input$search_var_c, {
        varFilterInput$data2 <- GRIN_ACMG %>% filter(gene==input$var_gene) %>% filter(c_pos==input$search_c_pos) %>%
            filter(ref_c==input$search_ref_c) %>% filter(alt_c==input$search_alt_c) %>% distinct(alt_c, .keep_all=TRUE)
        
        print(ACMGFilterInput$data2)
    })
    
    observeEvent(input$search_var_p, {
        varFilterInput$data2 <- GRIN_ACMG %>% filter(gene==input$var_gene) %>% filter(aa_pos==input$search_aa_pos) %>%
            filter(ref_p==input$search_ref_p) %>% filter(alt_p==input$search_alt_p) %>% distinct(alt_p, .keep_all=TRUE)
        
        print(ACMGFilterInput$data2)
    })
    

    #### Variant Information ####  
    
    output$grinBox1 <- renderValueBox({
        valueBox(
            paste0("Gene: ",unique(varFilterInput$data$gene)), 
            div("Transcript: ",unique(varFilterInput$data$transcript),br(), "Exon: ",unique(varFilterInput$data$EXON)), icon = icon("dna"),
            color = "yellow"
        )
    })
    
    output$grinBox2 <- renderValueBox({
        
        p_old <- unique(varFilterInput$data$ref_p)
        p_new <- unique(varFilterInput$data$alt_p)
        
        valueBox(
            paste0("Domain: ", unique(varFilterInput$data$domain)), 
            div("Amino Acid Position: ",unique(varFilterInput$data$aa_pos), br(),
                paste0("Amino Acid Change: ", p_old, " (",seqinr::a(p_old), ") "), "-" ,
                paste0(p_new, " (",seqinr::a(p_new),")")), icon = icon("dna"),
            color = "light-blue"
        )
    })
    
    output$grinBox3 <- renderValueBox({
        
       x<- unique(varFilterInput$data$gnomad_allele_count)
       y<- unique(varFilterInput$data$gnomad_allele_frequency)
        
        valueBox(
            paste0("Control variants: ", ifelse(is.na(x), 0, x)), div(paste0("gnomAD Allele Count for this Position"),br(),
                paste0("gnomAD Allele Frequency: ", ifelse(is.na(y), 0, y))), icon = icon("dna"),
            color = "green"
        )
    })
    
    output$grinBox4 <- renderValueBox({
      
      x <- unique(varFilterInput$data$CALL)
      
      valueBox(width = 12,
        paste0(ifelse(is.na(x), 
          "Not available", 
          x)
        ), 
      div(paste0("Functional Consequence")), icon = icon("microscope"),
        color = "maroon"
      )
    })
    
    # ACMG Table
    output$ACMGTable <- DT::renderDataTable({
        validate(
            need(!plyr::empty(varFilterInput$data),
                 "The variant has not been evaluated yet."))
        
        datatable(GRIN_ACMG %>% filter(gene==varGeneInput()) %>% 
                    filter(aa_pos==varFilterInput$data$aa_pos[1]) %>%
                    filter(alt_p == varFilterInput$data$alt_p[1]) %>%
                    filter(alt_c == varFilterInput$data$alt_c[1]) %>%
                    select(gene, variant_c, variant_p, x, ACMG),
                  colnames = c("Gene", "cDNA Level", "Protein Level", "Criteria", "Clinical Significance"),
                  options = list(dom = 't', scrollY = TRUE), escape=FALSE)   
    })
    
    output$ACMGUI_crit <- renderUI({
        
        validate(
            need(!plyr::empty(varFilterInput$data),
                 "The variant has not been evaluated yet."))
        
        sub <- GRIN_ACMG %>% filter(gene==varGeneInput()) %>% 
          filter(aa_pos==varFilterInput$data$aa_pos[1]) %>% 
          filter(alt_c==varFilterInput$data$alt_c[1]) %>%
          filter(alt_p==varFilterInput$data$alt_p[1]) 
        
        print(sub)

       paragraph <- p(align="justify", div(
         h1(strong(sub$ACMG, style="color:#e04c3c;"
           ))), 
         div(strong(paste0(sub$transcript,": ", sub$variant_c," ", sub$variant_p))), br(), fluidRow())

    })
    
    output$ACMGUI_exp <- renderDataTable({
      
      validate(
        need(!plyr::empty(varFilterInput$data),
             "The variant has not been evaluated yet."))
      
      sub <- DT::datatable(GRIN_ACMG %>% filter(gene==varGeneInput()) %>% filter(aa_pos==varFilterInput$data$aa_pos) %>% 
                             filter(alt_c==varFilterInput$data$alt_c)%>% dplyr::distinct("gene", "variant_c", .keep_all=TRUE) %>%
                             select(PVS1, PS1, PS2, PS3, PS4, PM1, PM2, PM5, PM6, PP2, PP3), options = list(dom = 't') 
                           ) 
      return(sub)
      
    })
    
    
    # Table with patient information 
    output$patientTable <- DT::renderDataTable({
        
      validate(
        need(!plyr::empty(varFilterInput$data),
             "There is no data that matches your filters.")) 
      
      datatable(GRIN_clin %>% filter(gene==varGeneInput()) %>% 
                  filter(aa_pos==varFilterInput$data$aa_pos[1]) %>%
                  filter(alt_p == varFilterInput$data$alt_p[1]) %>%
                  select(gene, domain, variant_c, variant_p, origin, seizures, delay_degree, movement_disorder, ASD, ataxia, CVI, MCD, PMID),
                colnames = c("Gene", "Domain", "cDNA level", "Protein level", "Origin", 
                             "Seizures", "DD/ID","MD","ASD","Ataxia","CVI","MCD",   "Link"),
                options = list(dom = 't', scrollY = TRUE), escape=FALSE)
    })
    
    output$grinspider <- renderPlotly({
        
        validate(need(
            !plyr::empty(varFilterInput$data),
            "There is no data that matches your filters."
        ))
        
        
        y <- varFilterInput$data
        
        y <- y %>% mutate(across(19:23, ~case_when(. > 5 ~ 6,
                                                   . < -5 ~ -6, 
                                                   TRUE ~ .)))
        
        fig <- plot_ly(type = 'scatterpolar', fill = 'toself',
                       textfont = list(family = "Arial"))
        
        fig <- fig %>%
            add_trace(r = c(1, 1, 1, 1, 1),
                      theta = colnames(y[, 19:23]),
                      name = "wildtype",
                      fillcolor = "#a0d8e4",
                      opacity = 0.8,
                      marker = list(color = '#6cc3d5',
                                    opacity = 0.6))
        
        fig <- fig %>%
            add_trace(r = as.numeric(y[1, 19:23]),
                      theta = colnames(y[, 19:23]),
                      name = "variant",
                      fillcolor = "#ffa88e",
                      opacity = 0.5,
                      marker = list(color = '#ff7851',
                                    opacity = 0.6))
        
        fig <- fig %>%
            layout(title = paste0(
                    "Functional Consequence: p.", y$ref_p[1], y$aa_pos[1], y$alt_p[1],
                    "<br><i>Surface Expression: ", ifelse(is.na(y$`Fold_Surface/Total`[1]),
                                                          "unknown", paste0(round(y$`Fold_Surface/Total`[1], digits=2), " fold")), 
                    " compared to Wildtype</i>"
                    ),
                   
                font=plotly_font,
                autosize = T, #width = "100%", height ="100%",
                polar = list(radialaxis = list(visible = T, range = c(-6, 6))),
                showlegend = T) %>% 
            config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE) 
        
        fig
        
    })
    
    output$spiderTable <- DT::renderDataTable({
      
      y <- varFilterInput$data
      
      validate(
        need(!plyr::empty(varFilterInput$data),
             "There is no data that matches your filters.")) 
      
      datatable(y %>% select(Fold_glu, Fold_gly, Fold_mg, Fold_zn, Fold_ph) %>% distinct() %>% mutate(across(1:5, ~round(., digits=1))), 
                colnames = c("Fold Glu", "Fold Gly", "Fold Mg", "Fold Zn", "Fold pH"),
                options = list(dom = 't', scrollY = TRUE), escape=FALSE) 
    })
    
    output$spiderSources <- renderText({ 
      
      validate(need(
        !plyr::empty(varFilterInput$data),
        "No Summary available."
      ))
      
      
      y <- varFilterInput$data
      z <- ifelse(is.na(y$CALL[1]), "Not available", paste(y$CALL[1]))
      
    return(z)
      
    })
    
    output$paraTable <- DT::renderDataTable({
        
        validate(
            need(!plyr::empty(varFilterInput$data),
                 "There is no data that matches your filters.")) 
        
        datatable(GRIN_clin %>% filter(gene!=varGeneInput()) %>% 
                      filter(relative_aa_pos==varFilterInput$data$relative_aa_pos) %>%
                    select(gene, domain, variant_c, variant_p, origin, seizures, delay_degree, movement_disorder, ASD, ataxia, CVI, MCD, PMID),
                  colnames = c("Gene", "Domain", "cDNA level", "Protein level", "Origin", 
                               "Seizures", "DD/ID","MD","ASD","Ataxia","CVI","MCD",   "Link"),
                  options = list(dom = 't', scrollY = TRUE), escape=FALSE)
        
    })
    
    # table with similar variants according to user input 
    
    
    output$compareTable <- DT::renderDataTable({
        
        validate(
            need(!plyr::empty(varFilterInput$data),
                 "There is no data that matches your filters.")) 
        
        z <- GRIN_clin %>% 
            filter(case_when(
                input$compareButtons =="Variant Type" ~ var_type=="missense",
                input$compareButtons =="Domain" ~ domain==varFilterInput$data$domain,
                input$compareButtons == "Functional Consequence" ~ CALL == varFilterInput$data$CALL,
                input$compareButtons =="Amino Acid Position" ~ relative_aa_pos==varFilterInput$data$relative_aa_pos)) %>%
          select(gene, domain, variant_c, variant_p, CALL, 
                 origin, seizures, delay_degree, movement_disorder,  ASD, ataxia, CVI, MCD, PMID)
        
        if (input$hide_unknown == TRUE) {
          z <-
            z %>% filter(
              delay_degree != "Unknown",
              seizures != "Unknown",
             CALL != "Unknown",
             CALL != "Insufficient Data"
            )
        }
        
        if (input$same_GRIN == TRUE) {
            z <-
                z %>% filter(gene == input$var_gene)
        }
        
        datatable(z, extensions = "Buttons", 
                  colnames = c("Gene", "Domain", "cDNA level", "Protein level", "Functional Consequence", 
                               "Origin", 
                               "Seizures", "DD/ID","MD", "ASD","Ataxia","CVI","MCD", "Link"),
                  options = list(dom = 'Brtip',
                                 buttons = c('csv', 'excel'), pageLength=100, scrollY = "350px"), escape = FALSE)
    })
    
    output$comparePlot <- renderPlotly({
        
        validate(need(
            !plyr::empty(varFilterInput$data),
            "There is no data that matches your filters."
        ))
        
        z <- GRIN_clin %>%
            filter( case_when(
                    input$compareButtons == "Variant Type" ~ var_type == "missense",
                    input$compareButtons == "Domain" ~ domain == varFilterInput$data$domain,
                    input$compareButtons == "Functional Consequence" ~ CALL == varFilterInput$data$CALL,
                    input$compareButtons == "Amino Acid Position" ~ relative_aa_pos == varFilterInput$data$relative_aa_pos)
            ) %>% select(gene, domain, dataset, variant_c, variant_p, 
                         CALL, 
                         origin, seizures, delay_degree, PMID) %>% 
          mutate(
              CALL = case_when(is.na(CALL) ~ "Unknown",
                               TRUE ~ CALL)) %>%
     mutate(seizures = case_when(seizures == "Yes" ~ "Seizures Present", 
                                              seizures == "No" ~ "Seizures Not Present",
                                              TRUE ~ as.character(seizures)))
        
        if (input$hide_unknown == TRUE) {
            z <-
                z %>% filter(
                    delay_degree != "Unknown",
                    seizures != "Unknown" ,
                    CALL != "Unknown",
                    CALL != "Not enough Data to determine"
                )
        }
        
        if (input$same_GRIN == TRUE) {
            z <-  z %>% filter(gene == input$var_gene) 
        }
        
        z <- z  %>% add_count(seizures, name = "n_seizures") %>% add_count(delay_degree, name = "n_delay") %>% add_count(CALL, name="n_CALL")
        
        
        plotty <- plot_ly()
        
        plotty <- plotty %>%
            add_pie(
                data = z %>% distinct(delay_degree, n_delay),
                labels =  ~ delay_degree,
                values =  ~ n_delay,
                name = "ID",
                textinfo = "label+percent",
                marker = list(
                    colors = c(
                        "#a0d8e4",
                        "#92dec0",
                        "#ffdf9d",
                        "#ffa88e",
                        "#f7bbbe",
                        "#a8d8ca",
                        "#eff2ee",
                        "#9c9c7c"
                    ) ,
                    line = list(color = '#FFFFFF', width = 1)
                ),
                insidetextfont = list(color = '#333333'),
                domain = list(row = 0, column = 0)
            )
        
        
        
        plotty <- plotty %>%
          add_pie(
            data = z %>%
              distinct(seizures, n_seizures),
            labels =  ~ seizures,
            values =  ~ n_seizures,
            name = "Seizures",
            textinfo = "label+percent",
            insidetextfont = list(color = '#333333'),
            marker = list(
              colors = c("#a0d8e4", "#ffdf9d", "#f7bbbe", "#eff2ee"),
              line = list(color = '#FFFFFF', width = 1)
            ),
            domain = list(row = 0, column = 1)
          )
        
       plotty <- plotty %>%
         add_pie(
            data = z %>%
              distinct(CALL, n_CALL),
            labels =  ~ CALL,
            values =  ~ n_CALL,
            name = "Functional Consequence",
            textinfo = "label+percent",
            insidetextfont = list(color = '#333333'),
            marker = list(
              colors = c("#fbb4ae",
                "#b3cde3",
                "#ccebc5",
                "#decbe4",
                "#fed9a6",
                "#ffffcc",
                "#e5d8bd",
               "#fddaec"),
             line = list(color = '#FFFFFF', width = 1)
           ),
          domain = list(row = 0, column = 2)
         )
        
        plotty <-
            plotty %>% layout(
                title = "Information on Intellectual Disability and Seizures and Functional Consequence",
                showlegend = F, font = plotly_font, autosize = T, grid = list(rows = 1, columns = 3),
                xaxis = list(showgrid = FALSE, zeroline = FALSE, howticklabels = FALSE),
                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
            )
        
        return(plotty)
        
        
    })
    
    # update cDNA and protein input according to user 
    
    observeEvent(input$search_var_c, {
      updateNumericInputIcon(
        session = session,
        inputId = c("search_aa_pos"),
        value = c(unique(varFilterInput$data$aa_pos))
      )
      
      updatePickerInput(
        session = session,
        inputId = c("search_ref_p"),
        choices = c(unique(varFilterInput$data$ref_p))
      )
      
      updatePickerInput(
        session = session,
        inputId = c("search_alt_p"),
        choices = c(unique(varFilterInput$data$alt_p))
      )
      
    })
    
    observeEvent(input$search_var_p, {
      updateNumericInputIcon(
        session = session,
        inputId = c("search_c_pos"),
        value = c(unique(varFilterInput$data$c_pos))
      )
      
      updatePickerInput(
        session = session,
        inputId = c("search_ref_c"),
        choices = c(unique(varFilterInput$data$ref_c))
      )
      
      updatePickerInput(
        session = session,
        inputId = c("search_alt_c"),
        choices = c(unique(varFilterInput$data$alt_c))
      )
      
    })
    
    ##### RESEARCH ##### 
    
    # Filter for subset of variants
    res_mod <- callModule(
        module = selectizeGroupServer,
        id = "grin-filters",
        data = GRIN_research,
        vars = c("gene","var_type",  "alt_p", "domain", "delay_degree", "seizures", "CALL")
    )
    
    output$filtered_n <- renderText({
      x <- nrow(res_mod())
      x <- paste("Number of individuals:", x)
      return(x)
    })
    
    # Table with displayed variants 
    output$subsetTable <- DT::renderDataTable({
        req(res_mod())
      
        z <- res_mod() %>% filter(!is.na(variant_p)) %>% filter(origin != "biparental") %>%
          select(gene, domain, variant_c, variant_p, origin, seizures, delay_degree, movement_disorder, ASD, ataxia, CVI, MCD, PMID)
        
        patient_table <- datatable(z, extensions = "Buttons",
                  colnames = c("Gene", "Domain", "cDNA level", "Protein level", "Origin", 
                               "Seizures", "DD/ID","MD", "ASD","Ataxia","CVI","MCD", "Link"), escape=FALSE, 
                  options = list(dom = 'Brtip',buttons = c('csv', 'excel'), pageLength=300, scrollY = "350px"))
        
        x <- res_mod() %>% filter(!is.na(variant_p)) %>% 
          filter(origin != "biparental") %>% select(gene, variant_c, variant_p,
                                                    Fold_glu, Fold_gly, Fold_mg, Fold_ph, Fold_zn,
                                                    `Fold Popen`,	`Fold_tauW (ms)`,	`Fold_Surface/Total`,
                                                    `Fold Synaptic Charge Transfer`, `Fold Non-Synaptic Charge Transfer`, CALL
          ) %>% 
          mutate(across(4:13, ~round(., digits=1))) %>%  
          filter(!is.na(CALL)) %>% unique()
        
        functional_table <- datatable(x, extensions = "Buttons",
                  escape=FALSE, 
                  colnames = c("Gene", "cDNA level", "Protein level", "Fold Glu", "Fold Gly", "Fold Mg",
                               "Fold pH", "Fold Zn", "Fold p(Open)", "Fold tauW (ms)", "Fold Surface/ Total", "Fold Synaptic Charge Transfer",
                               "Fold Non-Synaptic Charge Transfer", "Summary"),
                  options = list(dom = 'Brtip',buttons = c('csv', 'excel'), pageLength=300, scrollX=TRUE, scrollY = "350px"))
        
        if (input$patientFunSwitch == FALSE) {
          return(patient_table)
        } else {
          return(functional_table)
        }
    })
    
    ### Genotype Interface
    
    output$GRIN_overview_plot <- renderPlotly({
        
        # 2D lolliplot with GRIN variants 
        g <- ggplot(data=res_mod()  %>% filter(origin != "biparental") %>% add_count(gene, variant_c, name = "var_count")
        )+
            geom_segment(aes(x=aa_pos, xend=aa_pos, y=4, yend=ifelse(var_type=="missense", 7,8)), colour="black")+
            geom_point(aes(x=aa_pos, y=ifelse(var_type=="missense", 7,8), color=var_type, text=paste(variant_c, variant_p, var_count, "Variant Count")))+
            geom_rect(data=GRIN_domains, aes(xmin=start, xmax=end, ymin=3, ymax=4, fill=domain_gen, text=domain_spec))+
            theme_classic()+
            ylim(c(1,10))+
            scale_color_manual(values = lolliplot_fill_scheme)+
            scale_fill_manual(values = lolliplot_fill_scheme)+
            facet_grid(gene ~ .)+ 
            theme(
                text = element_text(size = 10),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank(),
                legend.position = "none"
            ) 
        
        if (input$gnomad_m == TRUE) {
            g <- g + geom_point(data=GRIN_gnomad %>% filter(var_type=="missense"), 
                                size=1, aes(x=aa_pos, y=2, alpha=0.1*gnomad_allele_count, text=paste0("Position: ",aa_pos,", Allele count: ", gnomad_allele_count)))
        }
        
        #if (input$gnomad_plof == TRUE) {
        #    g <- g + geom_point(data=GRIN_gnomad %>% filter(var_type=="null"), 
        #                        size=1, shape=17, aes(x=aa_pos, y=1, alpha=0.5*gnomad_allele_count,  text=paste("allele count:", gnomad_allele_count)))
        #}
        
        
        g <- ggplotly(g, tooltip = "text") %>% 
            config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  %>% 
            layout(title="",font=plotly_font
            ) 
        
    })
    
    # 3D plot 
    output$threeDmolGRIN1 <- renderR3dmol({
        
      if (input$map_3d_1 == TRUE) {
        
        variant.df <- res_mod() %>% 
          filter(gene == "GRIN2A" | gene == "GRIN1") %>%
          filter(!str_detect(variant_p,"\\?"),
                 !str_detect(variant_p,"\\*"),
                 !str_detect(variant_p,"dup"),
                 !str_detect(variant_p,"ins"),
                 !str_detect(variant_p,"del")) %>% 
          mutate(aa_ref = str_sub(variant_p,4,6),
                 label = "pathogenic") %>% 
          group_by(aa_pos,aa_ref,gene) %>% 
          summarise(n_occ = n()) %>% 
          select(aa_pos,aa_ref,n_occ,gene) 
        
        gnomad.df <- GRIN_gnomad %>% 
          filter(gene == "GRIN2A" | gene == "GRIN1") %>%
          mutate(aa_ref = ref_p) %>% 
          group_by(aa_pos,aa_ref,gene) %>% 
          summarise(n_occ = n()) %>% 
          select(aa_pos,aa_ref,n_occ,gene) 
        
        structure.df <- read_delim("data/pdb/6ira_structure_coordinates.txt",delim = "\t") %>% 
          mutate(Aminoacid = aaa(Aminoacid)) %>% 
          select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain) 
        
      } else {
        
        variant.df <- r3dmol_patho_fun("GRIN1", res_mod())
        
        gnomad.df <- r3dmol_gnomad_fun("GRIN1")  

        structure.df <- read_delim("data/pdb/GRIN1_subunit_coordinates.txt",delim = "\t") %>% 
          mutate(Aminoacid = aaa(Aminoacid)) %>% 
          select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
      }
      
      variant.df <- variant.df %>% 
        left_join(structure.df,by = c("aa_pos" = "Uniprot_position","aa_ref" = "Aminoacid","gene" = "gene")) %>% 
        mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not 
        filter(struc_cov == "yes") 
      
      gnomad.df <- gnomad.df %>% 
        left_join(structure.df,by = c("aa_pos" = "Uniprot_position","aa_ref" = "Aminoacid","gene" = "gene")) %>% 
        mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not 
        filter(struc_cov == "yes") 
      
      sub_color <- c("#ff1919")
      sub_scale <- c(1.2,0.8)
      struc_color <- "wheat"
      
      rot = 270 
      rot_axis = "x"
      spin_axis = "vy"
      
      #Specify yourself- color of the cartoon per subunit 
      subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
      
      #Model for the protein complex 
      
      modelo <- r3dmol(
        viewer_spec = m_viewer_spec(
          cartoonQuality = 10,
          lowerZoomLimit = 50,
          upperZoomLimit = 1000
        )
      ) 
      if (input$map_3d_1 == TRUE) {
        modelo <- modelo %>% m_add_model(data = "data/pdb/GRIN2A_GRIN1.pdb1", format = "pdb") 
      } else {
        modelo <- modelo %>% m_add_model(data = "data/pdb/GRIN1_subunit.pdb", format = "pdb")
      }
      
      if (input$map_3d_1 == TRUE) {
        # Zoom to encompass the whole scene
        modelo <- modelo %>% m_zoom_to() %>% 
          # Set color o cartoon representation
          m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure 
          # Set subunit colors 
          m_set_style(
            sel = m_sel(chain = c("A","C")), ##grin1
            style = m_style_cartoon(color = subunit_color[1])
          ) %>% 
          m_set_style(
            sel = m_sel(chain = c("B","D")), ##grin2a
            style = m_style_cartoon(color = subunit_color[2])
          ) %>% 
          # visualize variants grin1
          m_set_style(
            sel = m_sel(resi = variant.df$Position_in_structure,
                        atom = "CA",
                        chain = c("A","C")),
            style = m_style_sphere(colorScheme = NULL,
                                   color = sub_color[1],
                                   scale = sub_scale[1])
          ) %>% 
          
          # visualize variants grin2a
          m_set_style(
            sel = m_sel(resi = variant.df$Position_in_structure,
                        atom = "CA",
                        chain = c("B","D")),
            style = m_style_sphere(colorScheme = NULL,
                                   color = "#ff1919",
                                   scale = sub_scale[1])
          ) %>% 
          m_set_style(
            sel = m_sel(resi = variant.df$Position_in_structure,
                        atom = "CA",
                        chain = c("B","D")),
            style = m_style_sphere(colorScheme = NULL,
                                   color = "#ff1919",
                                   scale = sub_scale[1])
          )} else {
            modelo <- modelo %>% m_zoom_to() %>% m_set_style(style = m_style_cartoon(color = "wheat")) %>% # select a color of the structure 
              m_set_style(sel = m_sel(chain = "A"),
                          style = m_style_cartoon(color = "wheat")) %>% 
              m_set_style(
                sel = m_sel(resi = variant.df$Position_in_structure,
                            atom = "CA"),
                style = m_style_sphere(colorScheme = NULL,
                                       color = sub_color[1],
                                       scale = sub_scale[1]))
          }
      
      if  (input$gnomad_m == TRUE) {
        
        modelo <- modelo %>% m_set_style(
          sel = m_sel(resi = gnomad.df$Position_in_structure,
                      atom = "CA"),
          style = m_style_sphere(colorScheme = NULL,
                                 color = "#333333",
                                 scale = sub_scale[2]))
        
      }    
      
      
      return(modelo)
        
    })
    
    output$threeDmolGRIN2A <- renderR3dmol({
      if (input$map_3d_2A == TRUE){
        
        variant.df <- res_mod() %>%
          filter(gene == "GRIN2A" | gene == "GRIN1") %>%
          filter(!str_detect(variant_p,"\\?"),
                       !str_detect(variant_p,"\\*"),
                       !str_detect(variant_p,"dup"),
                       !str_detect(variant_p,"ins"),
                       !str_detect(variant_p,"del")) %>% 
                mutate(aa_ref = str_sub(variant_p,4,6),
                       label = "pathogenic") %>% 
                group_by(aa_pos,aa_ref,gene) %>% 
                summarise(n_occ = n()) %>% 
                select(aa_pos,aa_ref,n_occ,gene) 
            
            gnomad.df <- GRIN_gnomad %>% 
              filter(gene == "GRIN2A" | gene == "GRIN1") %>%
              mutate(aa_ref = ref_p) %>% 
              group_by(aa_pos,aa_ref,gene) %>% 
              summarise(n_occ = n()) %>% 
              select(aa_pos,aa_ref,n_occ,gene) 

            structure.df <- read_delim("data/pdb/6ira_structure_coordinates.txt",delim = "\t") %>% 
                mutate(Aminoacid = aaa(Aminoacid)) %>% 
                select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
            
        } else {
          
          variant.df <- r3dmol_patho_fun("GRIN2A", res_mod())
          
          gnomad.df <- r3dmol_gnomad_fun("GRIN2A")
          
          structure.df <- read_delim("data/pdb/GRIN2A_subunit_coordinates.txt",delim = "\t") %>% 
            mutate(Aminoacid = aaa(Aminoacid)) %>% 
            select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
          
        }
            variant.df <- variant.df %>% 
              left_join(structure.df,by = c("aa_pos" = "Uniprot_position","aa_ref" = "Aminoacid","gene" = "gene")) %>%
              mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not 
                filter(struc_cov == "yes") 
            
            gnomad.df <- gnomad.df %>% 
              left_join(structure.df,by = c("aa_pos" = "Uniprot_position","aa_ref" = "Aminoacid","gene" = "gene")) %>% 
              mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not 
              filter(struc_cov == "yes") 
            
            sub_color <- c("#ff1919")
            sub_scale <- c(1.2,0.8)
            struc_color <- "wheat"
            
            rot = 270 
            rot_axis = "x"
            spin_axis = "vy"
            
            #Specify yourself- color of the cartoon per subunit 
            subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
            
            #Model for the protein complex 
            
            modelo <- r3dmol(viewer_spec = m_viewer_spec(
                    cartoonQuality = 10,
                    lowerZoomLimit = 50,
                    upperZoomLimit = 1000
                    )) 
            
            if (input$map_3d_2A == TRUE) {
              modelo <- modelo %>% m_add_model(data = "data/pdb/GRIN2A_GRIN1.pdb1", format = "pdb") 
            } else {
              modelo <- modelo %>% m_add_model(data = "data/pdb/GRIN1_subunit.pdb", format = "pdb")
            }
            
            if (input$map_3d_2A == TRUE) {
              
            modelo <- modelo %>% 
                m_zoom_to() %>% 
                # Set color of cartoon representation
                m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure 
                # Set subunit colors 
                m_set_style(
                    sel = m_sel(chain = c("A","C")), ##grin1
                    style = m_style_cartoon(color = subunit_color[1])
                ) %>% 
                m_set_style(
                    sel = m_sel(chain = c("B","D")), ##grin2a
                    style = m_style_cartoon(color = subunit_color[2])
                ) %>% 
                # visualize variants grin1
                m_set_style(
                    sel = m_sel(resi = variant.df$Position_in_structure,
                                atom = "CA",
                                chain = c("A","C")),
                    style = m_style_sphere(colorScheme = NULL,
                                           color = sub_color[1],
                                           scale = sub_scale[1])
                ) %>% 
                # visualize variants grin2a
                m_set_style(
                    sel = m_sel(resi = variant.df$Position_in_structure,
                                atom = "CA",
                                chain = c("B","D")),
                    style = m_style_sphere(colorScheme = NULL,
                                           color = "#ff1919",
                                           scale = sub_scale[1]))
        } else {
            modelo <- modelo %>% m_zoom_to() %>% 
                # Set color of cartoon representation 
                m_set_style(style = m_style_cartoon(color = "wheat")) %>% # select a color of the structure 
                m_set_style(sel = m_sel(chain = "A"),
                            style = m_style_cartoon(color = "wheat")) %>% 
                m_set_style(
                    sel = m_sel(resi = variant.df$Position_in_structure,
                                atom = "CA"),
                    style = m_style_sphere(colorScheme = NULL,
                                           color = sub_color[1],
                                           scale = sub_scale[1])) 
        }
            
            if  (input$gnomad_m == TRUE) {
              
              modelo <- modelo %>% m_set_style(
                sel = m_sel(resi = gnomad.df$Position_in_structure,
                            atom = "CA"),
                style = m_style_sphere(colorScheme = NULL,
                                       color = "#333333",
                                       scale = sub_scale[2]))
              
            }    


        
        
        return(modelo)
        
    })
    
    output$threeDmolGRIN2B <- renderR3dmol({
        
      variant.df <- r3dmol_patho_fun("GRIN2B", res_mod())
      gnomad.df <- r3dmol_gnomad_fun("GRIN2B")
        
        
        structure.df <- read.pdb("data/pdb/GRIN2B_model.pdb") %>% 
            .$atom %>% 
            distinct(resno,chain,resid) %>% 
            mutate(Aminoacid = paste0(str_sub(resid,1,1),str_sub(resid,2,3) %>% tolower()),
                   struc_cov = "yes",
                   Uniprot_position = resno,
                   Position_in_structure = resno,
                   gene = "GRIN2B",
                   chain = "Dummy") %>% 
            select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
        
        variant.df <- variant.df %>% 
            left_join(structure.df,by = c("aa_pos" = "Uniprot_position","aa_ref" = "Aminoacid","gene" = "gene")) %>% 
            mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not 
            filter(struc_cov == "yes") 
        
        gnomad.df <- gnomad.df %>% 
          left_join(structure.df,by = c("aa_pos" = "Uniprot_position","aa_ref" = "Aminoacid","gene" = "gene")) %>% 
          mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not 
          filter(struc_cov == "yes") 
        
        
        sub_color <- c("#ff1919")
        sub_scale <- c(1.2,0.8)
        struc_color <- "wheat"
        
        rot = 90
        rot_axis = "x"
        spin_axis = "vy"
        subunit_color <- "white"
        
       modelo <- r3dmol(
            viewer_spec = m_viewer_spec(
                cartoonQuality = 10,
                lowerZoomLimit = 50,
                upperZoomLimit = 1000
                )) %>% 
            # Add model to scene
            m_add_model(data = "data/pdb/GRIN2B_model.pdb", format = "pdb") %>% 
            # Zoom to encompass the whole scene
            m_zoom_to() %>% 
            # Set color of cartoon representation 
            m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure 
            m_set_style(sel = m_sel(chain = "A"),
                        style = m_style_cartoon(color = subunit_color)) %>% 
            m_set_style(
                sel = m_sel(resi = variant.df$Position_in_structure,
                            atom = "CA"),
                style = m_style_sphere(colorScheme = NULL,
                                       color = sub_color[1],
                                       scale = sub_scale[1])
            ) %>% 
            m_set_style(
                sel = m_sel(resi =variant.df$Position_in_structure,
                            atom = "CA"),
                style = m_style_sphere(colorScheme = NULL,
                                       color = sub_color[2],
                                       scale = sub_scale[2])
            ) 
        
        if  (input$gnomad_m == TRUE) {
          
          modelo <- modelo %>% m_set_style(
            sel = m_sel(resi = gnomad.df$Position_in_structure,
                        atom = "CA"),
            style = m_style_sphere(colorScheme = NULL,
                                   color = "black",
                                   scale = sub_scale[2])
          )
          
        } 
        
        return(modelo)
        
    })
    

    
    output$threeDmolGRIN2D <- renderR3dmol({
        
        variant.df <- r3dmol_patho_fun("GRIN2D", res_mod())
        gnomad.df <- r3dmol_gnomad_fun("GRIN2D")
        
        structure.df <- read.pdb("data/pdb/GRIN2D_model.pdb") %>% .$atom %>% 
            distinct(resno,chain,resid) %>% 
            mutate(Aminoacid = paste0(str_sub(resid,1,1),str_sub(resid,2,3) %>% tolower()),
                   struc_cov = "yes",
                   Uniprot_position = resno,
                   Position_in_structure = resno,
                   gene = "GRIN2D",
                   chain = "Dummy") %>% 
            select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
        
        variant.df <- variant.df %>% 
          left_join(structure.df,by = c("aa_pos" = "Uniprot_position","aa_ref" = "Aminoacid","gene" = "gene")) %>% 
          mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not 
          filter(struc_cov == "yes") 
        
        gnomad.df <- gnomad.df %>% 
          left_join(structure.df,by = c("aa_pos" = "Uniprot_position","aa_ref" = "Aminoacid","gene" = "gene")) %>% 
          mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not 
          filter(struc_cov == "yes") 
        
        sub_scale <- c(1.2,0.8)
        struc_color <- "wheat"
        subunit_color <- "white"
        sub_color <- "red"
        
       modelo <- r3dmol(
            viewer_spec = m_viewer_spec(
                cartoonQuality = 10,
                lowerZoomLimit = 50,
                upperZoomLimit = 1000
            )
        ) %>% 
            m_add_model(data = "data/pdb/GRIN2D_model.pdb", format = "pdb") %>% 
            m_zoom_to() %>% 
            m_set_style(style = m_style_cartoon(color = struc_color)) %>%
            m_set_style(sel = m_sel(chain = "A"),
                        style = m_style_cartoon(color = subunit_color)) %>% 
            m_set_style(
                sel = m_sel(resi = variant.df$Position_in_structure,
                            atom = "CA"),
                style = m_style_sphere(colorScheme = NULL,
                                       color = sub_color[1],
                                       scale = sub_scale[1])) 
        
        if  (input$gnomad_m == TRUE) {
          
          modelo <- modelo %>% m_set_style(
            sel = m_sel(resi = gnomad.df$Position_in_structure,
                        atom = "CA"),
            style = m_style_sphere(colorScheme = NULL,
                                   color = "#333333",
                                   scale = sub_scale[2]))
        } 
        
        return(modelo)
    })
    
    output$legend_plot <- renderPlot({ 
      

      
      legend <- data.frame(x=c(1,11,21), y=c(1, 1, 1), text=c("Missense", "Truncating", "Control"))
      plot <- ggplot(legend, aes(x=x, y=y, color=text))+
        geom_point(size = 6)+
        scale_color_manual(values = c("Missense"="#D55E00","Truncating"="#0072B2","Control" ="#000000"))+
        ylim(c(0,2))+
        xlim(c(0,40))+
        theme_void()+
        geom_text(aes(label=text), hjust=-0.4, color="black")+
        theme(legend.position = "none")

      
      return(plot)
      
    })
    
    ### Phenotype Interface
    
    output$phenotype_number_plot <- renderPlotly({
        
        plot <- res_mod() %>% group_by(gene) %>% ungroup() %>% select(gene, domain) %>%
          mutate(domain = str_replace(domain,"linker","L-")) %>%  
          mutate(domain = factor(domain, levels = c("SP", "ATD", "L-ATD-S1",
                                                    "S1", "L-S1-M1", "M1",
                                                    "L-M1-M2", "M2", "L-M2-M3",
                                                    "M3", "L-M3-S2", "S2", "L-S2-M4", "M4", "CTD"))) %>%
            group_by(gene, domain) %>% add_count(gene)  %>% ungroup() %>% distinct() %>%
          plot_ly(
            x=~domain, y=~n, split=~gene, type="bar", color=~gene, alpha = 0.8, colors= GRIN_colors) %>% 
            layout(title="Number of Variants with Phenotypic Data per Domain",font=plotly_font, yaxis = list(title = "Number"),
                   xaxis = list(title = "",  tickangle = 45)
            ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE) 
        
        return(plot)
        
    })
    
    output$violin_grin <- renderPlotly({
      plot <- res_mod() %>% filter(!is.na(delay_degree), delay_degree != "Unknown") %>% select(gene, variant_c, variant_p, delay_degree) %>% 
        mutate(delay_number = case_when(delay_degree == "severe/profound DD/ID" ~ 3,
                                        delay_degree == "moderate DD/ID" ~ 2,
                                        delay_degree == "mild DD/ID" ~ 1,
                                        delay_degree == "no DD/ID" ~ 0
                                        )) %>% 
        plot_ly(
          y =  ~ delay_number,
          x =  ~ gene,
          split =  ~ gene,
          type = 'violin',
          points = "all",
          pointpos = 0,
          opacity = 0.8,
          color = ~ gene,
          colors = GRIN_colors,
          jitter = 0.5, 
          opacity = 0.8,
          box = list(visible = T),
          meanline = list(visible = T),
          hoverinfo = "text",
          text =  ~ paste0("Variant: ", variant_c, " ", variant_p)
        ) %>%
        
        layout(
          title = "Degree of Intellectual Disability (ID)",
          font = plotly_font,
          yaxis = list(title = "Degree of ID", range =
                         c(-0.1, 3.1)),
          xaxis = list(title = "")
        ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
      
      return(plot)

        
    })
    
    output$seizure_gene <- renderPlotly({
        
        plot <- res_mod() %>% filter(!is.na(seizures), seizures != "Unknown") %>% add_count(gene, name = "gene_n") %>% 
          add_count(gene, seizures) %>% select(gene, seizures, gene_n, n) %>% unique() %>% mutate(n_new = round(n/gene_n, digits=2)) %>% 
          plot_ly(
                y =  ~ n_new,
                x =  ~ gene,
                split =  ~ seizures,
                type = 'bar',
                opacity = 0.8,
                hoverinfo = "text",
                color = ~ seizures,
                colors = c("Yes" = "#6fbbf7", "No" = "#ee6c71"),
                text =  ~ paste(n, "Individuals")
            ) %>%
            
            layout(
                title = "Proportion of Individuals with Seizures",
                font = plotly_font,
                xaxis = list(title = ""),
                yaxis = list(title = "Proportion of Cohort")
            ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
        
        return(plot)
        
    })
    
    ### Functional Interface
    
    output$functional_number_plot <- renderPlotly({
        
        z <-  res_mod() %>% filter(!is.na(CALL)) %>%
            group_by(gene) %>% distinct(variant_c, .keep_all = TRUE) %>% ungroup() %>% select(gene, domain) %>%
            mutate(domain = str_replace(domain, "linker", "L-")) %>% 
          mutate(domain = factor(domain,
                                 levels = c( "SP", "ATD", "L-ATD-S1",
                                             "S1", "L-S1-M1", "M1",
                                             "L-M1-M2", "M2", "L-M2-M3",
                                             "M3", "L-M3-S2", "S2", "L-S2-M4", "M4", "CTD"))) %>%
            group_by(gene, domain) %>% add_count(gene) %>%  distinct() %>% ungroup() 
        
        plot <-  plot_ly(z,
                         x =  ~ domain,
                         y =  ~ n,
                         split =  ~ gene,
                         type = "bar",
                         color =  ~ gene,
                         alpha = 0.8,
                         colors = GRIN_colors,
                         hoverinfo = "text",
                         text =  ~ paste0(domain, ": ", n, " Variant(s)")
                         ) %>%
          layout(title = "Number of Variants Tested per Domain",
                 font = plotly_font,
                 xaxis = list(title = "",  tickangle = 45),
                 yaxis = list(title = "Number")) %>% 
          config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
        
        return(plot)
    })
    
     output$functional_cons_plot <- renderPlotly({
      
      z <-  res_mod() %>% filter(!is.na(CALL)) %>%
        group_by(gene) %>% distinct(variant_c, .keep_all = TRUE) %>% ungroup() %>% 
       mutate(CALL = factor(CALL, levels = c("Likely Gain-of-Function", "Potential Gain-of-Function",
                                                 "Likely Loss-of-Function", "Potential Loss-of-Function",
                                                 "Complex", "No Detectable Effect", "Insufficient Data")))%>%
        group_by(gene, CALL) %>% add_count(gene) %>%  distinct() %>% ungroup() %>% select(gene, CALL, n) %>% distinct() 
      
      plot <-  plot_ly(z,
                       x =  ~ CALL,
                       y =  ~ n,
                       split =  ~ gene,
                       type = "bar",
                       color =  ~ gene,
                       alpha = 0.8,
                      colors = GRIN_colors,
                      hoverinfo = "text",
                      text =  ~ paste0(CALL, ": ", n, " Variant(s)")
      ) %>%
        layout(title = "Number of Variants per Functional Consequence",
               font = plotly_font,
               xaxis = list(title = "",  tickangle = 45),
               yaxis = list(title = "Number")) %>% 
        config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
      
      return(plot)
      })
    
    
    output$functional_Glu <- renderPlotly({
        x <- GRIN_functional_fun(res_mod(), "Fold_glu") %>% layout(title = "Glutamate")
        return(x)

    })
    
    output$functional_Gly <- renderPlotly({
        x <- GRIN_functional_fun(res_mod(), "Fold_gly") %>% layout(title = "Glycine")
        return(x)
    })
    
    output$functional_Mg <- renderPlotly({
        x <- GRIN_functional_fun(res_mod(), "Fold_mg") %>% layout(title = "Magnesium")
        return(x)
    })
    

    output$functional_Zn <- renderPlotly({
        x <- GRIN_functional_fun(res_mod(), "Fold_zn") %>% layout(title = "Zinc")
        return(x)
    })
    
    output$functional_tauW <- renderPlotly({
        x <- GRIN_functional_fun(res_mod(), "Fold_tauW (ms)") %>% layout(title = "tau W")
        return(x)
    })
    
    output$functional_surface_mean <- renderPlotly({
        x <- GRIN_functional_fun(res_mod(), "Fold_Surface/Total") %>% layout(title = "Surface Expression")
        return(x)
    })
    
    observeEvent(input$info, {
      show_alert(html = TRUE, width="80%",
        title = "Wild Type Values and Methods", showCloseButton = TRUE,
        text = div(style="font-size:12px", align="left", width="100%", div(align="center", img(src="table_oocyte_data.PNG", width="50%")), 
        br(),
        p(strong("L-Glutamate dose response studies"),":  
        The Barth’s solution for these studies contained (in mM):  
        NaCl (90), KCl (1.0), BaCl2 (0.5), HEPES (10), and EDTA (0.01), adjusted to pH 7.4 (at 23 degrees C) with NaOH.  
        The oocyte membrane potential was clamped at -40 mV.  After a steady baseline was obtained, 
        oocytes were perfused with increasing concentrations of L-glutamate for 0.75 min duration each in the continuous presence of 100 uM glycine. 
        Results at each L-glutamate concentration were normalized to the maximum receptor activation levels 
        (defined as 100%) and the EC50 values obtained by fitting concentration-response data with Equation 1 (below)."),
                   
        p(strong("Glycine dose response studies"),":  The Barth’s solution for these studies was the same as for the L-glutamate studies.  
        Oocyte membrane potential was clamped at -40 mV.  After a steady baseline was obtained, oocytes were perfused with 
        increasing concentrations of glycine for 0.75 min duration each in the continuous presence of 100 uM L-glutamate.  
        Results at each glycine concentration were normalized to the maximum receptor activation levels (defined as 100%) 
        and the EC50 values obtained by fitting concentration-response data with Equation 1 (below)."),
        
        p(strong("Mg2+ dose inhibition studies"), ":  The Barth’s solution for Mg2+ studies contained (in mM):  
        NaCl (90), KCl (1.0), BaCl2 (0.5), and HEPES (10) adjusted to pH 7.4 (at 23 degrees C) with NaOH.  
        The oocyte membrane potential was held at -60 mV.  After a steady baseline was obtained, oocytes were 
        maximally activated with 100 uM L-glutamate and 100 uM glycine and, in the continuous presence of maximal 
        L-glutamate and glycine, were perfused with increasing concentrations of Mg2+ (3, 10, 30, 100, 300, and 1000 uM). 
        Results at each Mg2+ concentration were normalized to the maximum receptor activation levels (defined as 100%) 
        and the IC50 values obtained by fitting concentration-inhibition data with Equation 2 (below)."),
        
        p(strong("pH studies"),":  The Barth’s solution for these studies contained (in mM):  NaCl (90), KCl (1.0), BaCl2 (0.5), 
        HEPES (10), and EDTA (0.01), adjusted to either pH 6.8 or pH 7.6 (at 23 degrees C) with HCl or NaOH. 
        The oocyte membrane potential was held at -40 mV.  After a steady baseline was obtained in pH 7.6 recording 
        buffer the oocytes were maximally activated with 100 uM L-glutamate and 100 uM glycine in pH 7.6 buffer.  
        Following a washout period to reestablish baseline, the oocytes were then maximally activated with 100 uM 
        L-glutamate and 100 uM glycine in pH 6.8 buffer.   The % current at pH 6.8 was then determined compared to 
        the current at H 7.6 (defined as 100%)."),
        
        p(strong("Zn2+ dose inhibition studies"),":  Oocytes expressing recombinant human glutamate receptors are perfused with Barth’s 
        solution containing (in mM):  NaCl (90), KCl (1.0), BaCl2 (0.5), Tricine (10), and HEPES (10) adjusted to pH 7.3 
        (at 23 degrees C) with NaOH.  The oocyte membrane potential is clamped at -20 mV.  After a steady baseline is obtained, 
        oocytes are maximally activated with 50 uM L-glutamate and 50 uM glycine, and then in the continuous presence of maximal 
        glutamate and glycine are perfused with increasing concentrations of Zn2+ (concentration varies depending on the specific 
        receptor and variant tested). Results at each Zn2+ concentration are normalized to the maximum receptor activation 
        levels without Zn2+ (defined as 100%) and IC50 values obtained by fitting concentration-inhibition data with Equation 2 (below)."),
        
        p(strong("beta-Lactamase (beta-lac) Assay"),":   HEK cells were plated in 96-well plates and transiently transfected with cDNA encoding 
        GluN-WT or beta-lac- GluN-Variant with a WT partner to form a mature NMDA receptor.  Eight wells were transfected and surface and 
        total expression activities were measured in 4 wells each 24 hr post transfection.  The ratio of surface (unlysed cells) to 
        total (lysed cells) beta-lactamase expression was measured for the  beta-lac-GluN-Variant and compared to b-lac-GluN-WT.  
        See PMID", shiny::a(href="https://pubmed.ncbi.nlm.nih.gov/27839871/", "27839871", target="_blank") ,"for further information."),
        
        p(strong("Determination of tau(weighted)"), ":  Average NMDA receptor-mediated currents were recorded from HEK cells transiently 
        expressing mutant and wild type NMDA receptors in response to rapid application and removal of 1 mM glutamate in the 
        presence of 0.03 mM glycine. Time course following removal of glutamate was fitted by equation 3 and when the weighted 
        tau calculated according to equation 4."),
        
        p(strong("Equation 1"), ":  Response = 100  / ( ( 1 + EC50 / [agonist] ) nH );  where EC50 is the agonist 
        concentration that elicited the half maximal response, and nH is the Hill slope."),
        
        p(strong("Equation 2"),":   Response  =  (100 - minimum) / (1 + ([concentration] / IC50)nH ) + minimum;   
        where minimum is the residual percent response in saturating concentration (constrained to > 0) 
        of the experimental compounds, IC50 is the concentration of antagonist that causes half maximal inhibition, and nH is the Hill slope."),
        
        p(strong("Equation 3"),": Response(time) =  Amplitude1 exp(-time/tau1) + Amplitude2 exp(-time/tau2)"),
        
        p(strong("Equation 4"), ": Tau(weighted) = Ampltidue1 tau1 / (Amplitude1 + Amplitude2) + Ampltidue2 tau2 / (Amplitude1 + Amplitude2)
        Statistical comparisons:  WT vs Variant receptor results were compared using a two-tailed unpaired t-Test (GraphPadPrism 5.0).  
        The log values of the IC50's or EC50's were used for comparison. When the 95% confidence intervals of the results do not overlap a 
        fold effect was then calculated.")),
        
        type = "info",
        position = "center",
        timer = FALSE
      )
    })
    
    
    
    ### Animal Models Tab
    
    output$animalTable <- DT::renderDataTable({
        
        datatable(GRIN_animal %>% select(-DOI, -aa_pos, -ref_p, -alt_p, -reference),
                  colnames = c("Gene", "Variant", "Model", "Reference", "Phenotype", 
                               "Functional Consequence", "Comment"),
                  options = list(dom = 't', scrollY = TRUE), escape=FALSE)
    })
        
    
})
