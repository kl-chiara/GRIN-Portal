################################################
################# GRIN PORTAL ##################
################################################
#
# This is a Shiny web application #
#

################## LIBRARIES ################

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(tippy)
library(r3dmol)
library(DT)
library(readr)
library(tidyverse)
library(vembedr)

################# CSS/STYLE/INFO #################
landing_panel <- "color: #333333;
      height: 200px;
      width: 260px"

spinner_color <- "#2196f3"

sub_style <- "color:gray; font-style: italic; font-size: 14px;"

GRIN1_text <- p(HTML("<em>GRIN1</em>-related disorder"),"is characterized by mild-to-profound developmental delay / 
              intellectual disability in all affected individuals. Other common manifestations are epilepsy, 
              muscular hypotonia, movement disorders, spasticity, feeding difficulties, and behavior problems. 
              A subset of individuals shows a malformation of cortical development consisting of extensive 
              and diffuse bilateral polymicrogyria. More details can be found on", 
                shiny::a(href="https://www.ncbi.nlm.nih.gov/books/NBK542807/", "GeneReviews."))

GRIN2A_text <- p(HTML("<em>GRIN2A</em>-related disorder"), "is characterized by either mild-to-profound developmental 
                delay / intellectual disability in two thirds of affected individuals. Other common 
                manifestations are speech disorders, epilepsy, muscular hypotonia, movement disorders, 
                spasticity, feeding difficulties, and behavior problems. More details can be found on", 
                 shiny::a(href="https://www.ncbi.nlm.nih.gov/books/NBK385627/", "GeneReviews."))

GRIN2B_text <- p(HTML("<em>GRIN2B</em>-related disorder"),"is characterized by mild-to-profound developmental delay / 
                intellectual disability in all affected individuals. Other common manifestations are epilepsy, 
                muscular hypotonia, movement disorders, spasticity, feeding difficulties, and behavior problems. 
                A subset of individuals shows a malformation of cortical development consisting of extensive 
                and diffuse bilateral polymicrogyria. More details can be found on", 
                 shiny::a(href="https://www.ncbi.nlm.nih.gov/books/NBK501979/", "GeneReviews."))

GRIN2D_text <- p(HTML("<em>GRIN2D</em>-related disorder"), "is characterized by moderate-to-profound developmental delay / 
                intellectual disability in all affected individuals. Other common manifestations are epilepsy, muscular hypotonia, 
                movement disorders, spasticity, feeding difficulties, and behavior problems.")

GRIN_clin <- read_csv("data/GRIN_clin.csv") %>%  
  mutate(delay_degree = factor(delay_degree, levels = c("no DD/ID", "mild DD/ID", "moderate DD/ID", "severe/profound DD/ID"))) %>%
  mutate(ref_p = case_when(ref_p == "linker S1-M1" ~ substr(variant_p, 4,6),
                                                               TRUE ~ ref_p))

################# UI #################
shinyUI(
    ##### LANDING PAGE #####
    ui <- 
    navbarPage(#fluid = FALSE, 
        windowTitle = "GRIN Portal",
        id = "TabDisplay",
       theme = "mytheme2.css",
        title = p(icon("dna"), "GRIN Portal"),

        header = tagList(useShinydashboard()),
       
       
        
        tabPanel(
            title = "Welcome",
            value = "welcomeTab",
            tabName = "welcomeTab",
            
            div(style = "font-size:100%", fluidRow(
                column(12, 
                    style = "background-color: white; color: #676767",
                    align = "center",
                    br(), br(),
                    img(src = "banner7.png", width = "100%")
                )
            ),
            
            div(style = 'background-color: #F8FCFE',
                  fluidRow(
                      
                      # Basic Information Tab
                      div(width = "100%", column(2, offset = 1, align = "center",
                          br(),
                         panel(width = 12, 
                              status = "success",
                              heading = "",
                              h2(tags$i(class = "fa fa-dna", style = "color:#676767")),
                              br(),
                              div("GRIN genes, function and associated disorders", 
                             br(), br(),
                              actionBttn(
                                  inputId = "infoBtn",
                                  label = "Basic Information",
                                  color = "success",
                                  block = TRUE,
                                  size = "md",
                                  style = "stretch")),
                              style =  "background-color: #f3faf4;"
                              ), 
                          ),
                      

                      column(2,align = "center",
                          br(),
                          div(panel(
                              status = "warning",
                              heading = "",
                              h2(tags$i(class = "fa fa-child", style = "color:#676767")),
                              br(),
                              div("Foundations, family groups, links to resources and more",br(), 
                              br(),
                              actionBttn(
                                  inputId = "familyBtn",
                                  label = "Families",
                                  color = "warning",
                                  block = TRUE,
                                  style = "stretch"
                              )),
                              style = "background-color: #fff8ef;"
                          ))
                      ),
                      
                      # Variant Analysis Tab
                      column(2, align = "center",
                          br(),
                          div(panel(
                              status = "info",
                              heading = "",
                              h2(tags$i(class = "fa fa-code-branch", style = "color:#676767")),
                              br(),
                              div("Comprehensive information on variant interpretation",br(), 
                              br(),
                              actionBttn(
                                  inputId = "variantBtn",
                                  label = "Variant Analysis",
                                  color = "royal",
                                  block = TRUE,
                                  style = "stretch"
                              )),
                              style = "background-color: #f9f1fa;"
                          ))
                      ),
                      
                      # Research Tab
                      column(2, align = "center",
                          br(),
                          div(panel(
                              status = "danger",
                              heading = "",
                              h2(tags$i(class = "fa fa-microscope", style = "color:#676767")),
                              br(),
                              div("Filter and select a subset of variants for research",br(), 
                              br(),
                              actionBttn(
                                  inputId = "researchBtn",
                                  label = "Research",
                                  color = "danger",
                                  block = TRUE,
                                  style = "stretch"
                              )),
                              style = "background-color: #fdf0f1;"),
                             
                          )
                      ),
                      
                      # Registry Tab
                      column(2, align = "center", br(),
                          div(panel(
                              status = "primary",
                              heading = "",
                              h2(tags$i(class = "fa fa-laptop-code", style = "color:#676767")),
                              br(),
                              div("Information on the GRIN Registry and how to register",br(),
                              br(),
                              actionBttn(
                                  #class = "btn-primary",
                                  inputId = "registryBtn",
                                  label = "GRIN Registry",
                                  color = "primary",
                                  block = TRUE,
                                  style = "stretch"
                              )), style = "background-color: #f1f8fe;")
                            
                          ), br(), br(), br(), br(), br(), br(), br(), br(), br()
                      ) 
                      
                      ))),
            )), # end tab Panel 
        
        ##### BASIC INFORMATION #####
        
        tabPanel(title = "Basic Information", value = "infoTab",
                 
                 fluidRow(
                   column(10, offset = 1,
                     panel(
                         heading = "GRIN Genes and Neurodevelopmental Disorders",
                         status = "primary",
                 
                 # GRIN1 
                 tabsetPanel(
                     type  = "tabs",
                     id = "displaygrin",
                     tabPanel(
                         title = "GRIN1",
                         value = "basic_GRIN1",
                         panel(heading = "Information on GRIN1",
                               fluidRow(column(4, box(width = 12,
                                                 
                                                 title = "History of GRIN1 Research",
                                                 
                                                 timelineBlock(
                                                     reversed = FALSE,
                                                     width = 12,
                                                     
                                                     timelineLabel(2011, color = "teal"),
                                                     timelineItem(
                                                         title = div("First Description of",  HTML("<em>de novo</em>"), "Variants"),
                                                         "in individuals with intellectual disability",
                                                         icon = icon("dna"),
                                                         color = "olive",
                                                         time = shiny::a("Hamdan et al.", href="https://pubmed.ncbi.nlm.nih.gov/21376300", target = '_blank'),
                                                         border = FALSE,
                                                     ),
                                                     
                                                     timelineLabel(2015, color = "aqua"),
                                                     timelineItem(title = "Delineation of the Phenotype",
                                                                  border = FALSE,
                                                                  icon = icon("user"),
                                                                  color = "teal",
                                                                  time = shiny::a("Ohba et al.", href="https://pubmed.ncbi.nlm.nih.gov/25864721", target = '_blank'),
                                                     ),
                                                     
                                                     timelineLabel(2016, color = "yellow"),
                                                     timelineItem(title = "Delineation of the Phenotype",
                                                                  border = FALSE,
                                                                  icon = icon("user"),
                                                                  color = "teal",
                                                                  time = shiny::a("Lemke et al.", href="https://pubmed.ncbi.nlm.nih.gov/27164704", target = '_blank'),
                                                     ),
                                                     
                                                     timelineLabel(2018, color = "red"),
                                                     timelineItem(title = "First Report of Cases with Brain Malformations",
                                                                  border = FALSE,
                                                                  icon = icon("brain"),
                                                                  color = "yellow",
                                                                  time = shiny::a("Fry et al.", href="https://pubmed.ncbi.nlm.nih.gov/29365063", target = '_blank')),
                                                     
                                                 )
                                   )),
                                   column(8, box(title="GRIN1-related Disorders", width = 12,
                                          GRIN1_text)),
                                   column(8,
                                          box(title = "Clinical Information from the GRIN Registry (96 Individuals)", width = 12,
                                              
                                              tabsetPanel(
                                                      
                                              tabPanel(title = "Developmental Delay",
                                                       br(),
                                                       fluidRow(align = "center",
                                                           column(6, box(title=div("GRIN1 Phenotypes", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN1_pie"))),
                                                           column(6, box(title=div("Developmental Delay / Intellectual Disability", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN1_delay"))),
                                                           column(12, align="center", p("Abbreviations: DD/ID, Developmental Delay / Intellectual Disability", style=sub_style)))),
                                              
                                              tabPanel(title = "Seizures",
                                                       br(),
                                                       fluidRow(align = "center",
                                                           column(6, box(title = div("Seizures", style = "font-size: 15px"), width = 12, plotlyOutput(outputId = "GRIN1_seizures"))),
                                                           column(6, box(title = div("Age at Seizure Onset", style = "font-size: 15px"), width = 12, plotlyOutput(outputId = "GRIN1_seizure_onset")))
                                                                  )),
                                              
                                              tabPanel(title="Further Phenotypes",
                                                       br(),
                                                       fluidRow(align = "center", 
                                                           column(6, box(title = div("Movement Disorder", style = "font-size: 15px"), width = 12, plotlyOutput(outputId = "GRIN1_movement"))),
                                                           column(6, box(title = div("Further Phenotypes", style = "font-size: 15px"), width = 12, 
                                                                  img(src = "GRIN1word.png", width = "100%"),
                                                                  footer = div("Frequently described phenotypes in the GRIN Registry. 
                                                                      Further detailed phenotypes will be added soon.", style=sub_style))
                                                                 
                                                                  )
                                                           ))
                                              )
                                              ))
                               ))
                     ),
                     tabPanel(
                         title = "GRIN2A",
                         value = "GRIN2AInfoTab",
                         panel(heading = "Information on GRIN2A",
                               fluidRow(
                                   column(4, box(width = 12,
                                                 
                                                 title = "History of GRIN2A Research",
                                                 
                                                 timelineBlock(
                                                     reversed = FALSE,
                                                     width = 12,
                                                     
                                                     timelineLabel(2010, color = "teal"),
                                                     timelineItem(
                                                         title = div("First Description of",  HTML("<em>de novo</em>"), "Variants"),
                                                         "in individuals with intellectual disability",
                                                         icon = icon("dna"),
                                                         color = "olive",
                                                         time = shiny::a("Endele et al.", href="https://pubmed.ncbi.nlm.nih.gov/20890276", target = '_blank'),
                                                         border = FALSE,
                                                     ),
                                                     
                                                     timelineLabel(2013, color = "aqua"),
                                                     timelineItem(title = "Delineation of the Phenotype",
                                                                  border = FALSE,
                                                                  icon = icon("user"),
                                                                  color = "teal",
                                                                  time = "multiple",
                                                                  div(shiny::a("Carvill et al.", href="https://pubmed.ncbi.nlm.nih.gov/23933818", target = '_blank'),
                                                                      shiny::a("Lemke et al.", href="https://pubmed.ncbi.nlm.nih.gov/23933819", target = '_blank'),
                                                                      shiny::a("Lesca et al.", href="https://pubmed.ncbi.nlm.nih.gov/23933820", target = '_blank')),
                                                     ),
                                                     
                                                     timelineLabel(2014, color = "yellow"),
                                                     timelineItem(title = "First Report on Treatment with Memantine",
                                                                  "in a Gain-of-Function case",
                                                                  border = FALSE,
                                                                  icon = icon("medkit"),
                                                                  color = "blue",
                                                                  time = shiny::a("Pierson et al.", href="https://pubmed.ncbi.nlm.nih.gov/24839611", target = '_blank'),
                                                     ),
                                                     
                                                     timelineLabel(2015, color = "red"),
                                                     timelineItem(title = "Delineation of the Speech Phenotype",
                                                                  border = FALSE,
                                                                  icon = icon("user"),
                                                                  color = "teal",
                                                                  time = shiny::a("Turner et al.", href="https://pubmed.ncbi.nlm.nih.gov/25596506", target = '_blank')),
                                                     
                                                     timelineLabel(2019, color = "green"),
                                                     timelineItem(title = "Comprehensive Genotype-Phenotype Correlation",
                                                                  border = FALSE,
                                                                  icon = icon("user"),
                                                                  color = "teal",
                                                                  time = shiny::a("Strehlow et al.", href="https://pubmed.ncbi.nlm.nih.gov/30544257", target = '_blank')),
                                                     
                                                 )
                                   )),
                                   column(8, box(width=12, title="GRIN2A-related Disorder",
                                          GRIN2A_text)),
                                  
                                   column(8,
                                          box(title = "Clinical Information from the GRIN Registry (158 Individuals)", width = 12,
                                          tabsetPanel(
                                              
                                              tabPanel(title = "Developmental Delay",
                                                       br(),
                                                       fluidRow(align = "center",
                                                           column(6, box(title=div("GRIN2A Phenotypes", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2A_pie"))),
                                                           column(6, box(title=div("Developmental Delay / Intellectual Disability", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2A_delay"))),
                                                           column(12, align="center", p("Abbreviations: DD/ID, Developmental Delay / Intellectual Disability", style=sub_style)))),
                                              
                                              tabPanel(title = "Seizures",
                                                       br(),
                                                       fluidRow(align = "center", 
                                                           column(4, box(title=div("Seizures", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2A_seizures"))),
                                                           column(4, box(title=div("Age at Seizure Onset", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2A_seizure_onset"))),
                                                           column(4, box(title=div("Age at Seizure Offset", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2A_seizure_offset"))))),
                                              
                                              tabPanel(title="Further Phenotypes",
                                                       br(),
                                                       fluidRow(align = "center", 
                                                           column(6, box(title=div("Movement Disorder", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2A_movement"))),
                                                           column(6, box(title=div("Further Phenotypes", style = "font-size: 15px"), width=12, align="center",
                                                                  img(src = "GRIN2Aword.png", width = "100%"),
                                                                  footer = div("Frequently described phenotypes in the GRIN Registry. 
                                                                      Further detailed phenotypes will be added soon.", style=sub_style))
                                                                  
                                                           )
                                                           ))
                                              
                                          ))
                               ))
                     )),
                     
                     tabPanel(
                         title = "GRIN2B",
                         value = "GRIN2BInfoTab",
                         panel(heading = "Information on GRIN2B",
                               fluidRow(column(4, box(width = 12,
                                                      
                                                      title ="History of GRIN2B Research",
                                                      
                                                      timelineBlock(
                                                          reversed = FALSE,
                                                          width = 12,
                                                          
                                                          timelineLabel(2010, color = "teal"),
                                                          timelineItem(
                                                              title = div("First Description of",  HTML("<em>de novo</em>"), "Variants"),
                                                              "in individuals with intellectual disability",
                                                              icon = icon("dna"),
                                                              color = "olive",
                                                              time = shiny::a("Endele et al.", href="https://pubmed.ncbi.nlm.nih.gov/20890276", target = '_blank'),
                                                              border = FALSE,
                                                              
                                                          ),
                                                          
                                                          timelineLabel(2014, color = "aqua"),
                                                          timelineItem(title = "Delineation of the Phenotype",
                                                                       border = FALSE,
                                                                       icon = icon("user"),
                                                                       color = "teal",
                                                                       time = shiny::a("Lemke et al.", href="https://pubmed.ncbi.nlm.nih.gov/23933819", target = '_blank')
                                                          ),
                                                          
                                                          timelineLabel(2017, color = "yellow"),
                                                          timelineItem(title = "Delineation of the Phenotype",
                                                                       "first report of cases with brain malformations, 
                                                                 first report on treatment with memantine in GoF cases",
                                                                       border = FALSE,
                                                                       icon = icon("medkit"),
                                                                       color = "blue",
                                                                       time = shiny::a("Platzer et al.", href="https://pubmed.ncbi.nlm.nih.gov/28377535", target = '_blank'),
                                                          ),
                                                          
                                                          timelineLabel(2019, color = "red"),
                                                          timelineItem(title = "First Report on L-Serine Treatment",
                                                                       "in a Loss-of-Function case",
                                                                       border = FALSE,
                                                                       icon = icon("medkit"),
                                                                       color = "blue",
                                                                       time = shiny::a("Soto et al.", href="https://pubmed.ncbi.nlm.nih.gov/31213567", target = '_blank')),
                                                          
                                                      )
                               )),
                               column(8, 
                                      box(title = "GRIN2B-related Disorders", GRIN2B_text, width = 12)),
                               column(8,
                                      box(title = "Clinical Information from the GRIN Registry (162 Individuals)", width = 12,
                                      tabsetPanel(
                                          
                                          tabPanel(title = "Developmental Delay",
                                                   br(),
                                                   fluidRow(align = "center",
                                                       column(6, box(title=div("GRIN2B Phenotypes", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2B_pie"))),
                                                       column(6, box(title=div("Developmental Delay / Intellectual Disability", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2B_delay"))),
                                                       column(12, align="center", p("Abbreviations: DD/ID, Developmental Delay / Intellectual Disability", style=sub_style)))),
                                          
                                          tabPanel(title = "Seizures",
                                                   br(),
                                                   fluidRow(align = "center",
                                                       column(6, box(title=div("Seizures", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2B_seizures"))),
                                                       column(6, box(title=div("Age at Seizure Onset", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2B_seizure_onset"))))),
                                          
                                          tabPanel(title="Further Phenotypes",
                                                   br(),
                                                   fluidRow(align = "center", 
                                                       column(6, box(title=div("Movement Disorder", style = "font-size: 15px"), width=12, plotlyOutput(outputId = "GRIN2B_movement"))),
                                                       column(6, box(title=div("Further Phenotypes", style = "font-size: 15px"), width=12, align="center",
                                                              div(width = "100%", img(src = "GRIN2Bword.png", width = "100%")),
                                                              footer = div("Frequently described phenotypes in the GRIN Registry. 
                                                                      Further detailed phenotypes will be added soon.", style=sub_style))
                                                              
                                                       )
                                                       ))
                                          
                                      ))
                               ))
                     )),
                     tabPanel(
                         title = "GRIN2D",
                         value = "GRIN2DInfoTab",
                         panel(heading = "Information on GRIN2D",
                               fluidRow(column(5, box(width = 12,
                                                      
                                                      title = "History of GRIN2D Research",
                                                      
                                                      timelineBlock(
                                                          reversed = FALSE,
                                                          width = 12,
                                                          
                                                          timelineLabel(2016, color = "teal"),
                                                          timelineItem(
                                                              title = div("First Description of",  HTML("<em>de novo</em>"), "Variants"),
                                                              "with epileptic encephalopythy, first report on treatment with memantine, 
                                                      ketamine and magnesium in GoF cases",
                                                              icon = icon("dna"),
                                                              color = "olive",
                                                              time = shiny::a("Li et al.", href="https://pubmed.ncbi.nlm.nih.gov/27616483", target = '_blank'),
                                                              border = FALSE,
                                                              #footer = "Here is the footer",
                                                              
                                                          ),
                                                          
                                                          timelineLabel(2019, color = "orange"),
                                                          timelineItem(title = "Delineation of the Phenotype",
                                                                       border = FALSE,
                                                                       icon = icon("user"),
                                                                       color = "teal",
                                                                       time = shiny::a("XiangWei et al.", href="https://pubmed.ncbi.nlm.nih.gov/31504254", target = '_blank')
                                                          ),
                                                      )
                               )),
                               column(7, 
                                      box(title ="GRIN2D-related Disorders", GRIN2D_text, div("Only few individuals have been reported in the
                                                                                              GRIN Registry yet."), width = 12)),
                               column(7)
                               ))
                     )
                     
                 ))))

                     
                 ), # end basic info tab panel
        
        ##### FAMILIES #####
        
        tabPanel(title = "Families", value = "familyTab", 
                 fluidRow(column(10, offset = 1,
                                 panel(heading = "For Families", status = "primary",
                                       tabsetPanel(
                                         tabPanel(title = "What are GRIN-related Disorders?",
                                                  fluidRow(column(5, offset = 1,
                                                  br(), br(),
                                                  box(width = 12, embed_url("https://www.youtube.com/watch?v=45rdftWLdlQ") %>%
                                                    use_rounded() %>% use_align("center") %>% use_bs_responsive(),
                                                  footer = p(div(align="center", "Expert curated subtitles available in English, German, 
                                                  Russian, Spanish, Italian and Romanian.", style=sub_style)))
                                                  ), 
                                                  
                                                  column(5, offset=0, br(), align = "justify",  br(), br(), br(),
                                                         box(title= div("GRIN Genes and Related Disorders", style = "color:white;"), 
                                                                                               width = 12,
                                                                                            status = "warning", solidHeader = TRUE,
                                                  p("If your family has received a diagnosis of GRIN-related disorders, 
                                                  you’ve probably never heard of it before. That’s okay. 
                                                  Your doctors probably haven’t either! But the GRIN Portal community is here to help."), 
                                                  p("This 6-minute", strong("“What are GRIN-related disorders?”") ,"video helps families learn about these 
                                                  disorders in plain language. You’ll learn:"),
                                                  div("• Scientific terms (such as genes, proteins, and variants/mutations) that will help you understand GRIN-related disorders"),
                                                  div("• What causes GRIN-related disorders"),
                                                  div("• How GRIN-related disorders are diagnosed"),
                                                  div("• The symptoms of GRIN-related disorders"), br(),
                                                  div("Thank you to all contributors!"))))),
                                         
                                         tabPanel(title = "Family Organizations",
                                       panel(heading = "International", width = 12,  
                                       fluidRow(column(2, align = "center", offset = 2,
                                                       div(width = 12, shiny::a(
                                               img(src = "https://curegrinfoundation.salsalabs.org/66d65f3a-a897-4ada-b2ce-80fdcbf9a832/21804aa6-eb96-41c9-9213-35a867f086f3.png", width = "100%"),
                                               href = "https://curegrin.org/",
                                               target = '_blank'))
                                           ),
                                           column(2, align = "center", br(), br(), 
                                                  div(width = 12, shiny::a(
                                                   img(src = "https://i0.wp.com/grin2b.com/wp-content/uploads/2017/07/Grin2BFoundation_logo_web-1.png", width =
                                                           "100%"),
                                                   href = "http://grin2b.com/",
                                                   target = '_blank'))
                                           ),
                                           column(2, align = "center", br(), br(), br(), 
                                                  div(width = 12, shiny::a(
                                                      img(src = "https://i0.wp.com/slc6a1connect.org/wp-content/uploads/2019/05/Searchlight_Logo_wSIMONS_RGB.png?resize=750%2C161&ssl=1", width = "100%"),
                                                      href = "https://www.simonssearchlight.org/?s=GRIN",
                                                      target = '_blank'))
                                           ),
                                           column(2, align = "center", br(), 
                                                  div(width = 12, shiny::a(
                                                    img(src = "https://www.gringn.com/wp-content/uploads/2019/03/cropped-logos-5-139x139.png", width = "50%"),
                                                    href = "https://www.gringn.com/",
                                                    target = '_blank')
                                           ))
                                           )
                                       ),
                                       panel(heading = "National", width = 12, 
                                             fluidRow(align = "center",
                                               column(3,
                                                      div(width="100%", panel(heading = "Italy", style = "height: 280px", status = "primary",br(),
                                                     shiny::a(
                                                       img(src = "https://insiemexmaurizio.files.wordpress.com/2020/06/cropped-definitivo-1.png?w=200", width = "50%"),
                                                       href = "https://assmgrin2aitalia.it/",
                                                       target = '_blank'),
                                                       div(width="100%", shiny::a(
                                                 img(src = "grin2b_italia.PNG", width = "50%"),
                                                 href = "http://www.grin2bitalia.it/grin2bitalia-presentazione.html",
                                                 target = '_blank'))
                                               ))),
                                           
                                              column(3, 
                                                     div(width=12, panel(heading= "Germany", style = "height: 280px", width = 3,  status = "primary",
                                                 shiny::a(
                                                 img(src = "gemeinsam_grin.jpg", width = "100%"),
                                                 href = "https://www.facebook.com/groups/208739597319166/",
                                                 target = '_blank')))),
                                           
                                           column(3, 
                                                  div(width = 12,  panel(heading = "Spain", style = "height: 280px", width = 3, br(), status = "primary",
                                               shiny::a(
                                                 img(src = "https://www.loto32malaga.com/Web_2_0/Cabeceras/grinpatias-logo_15_27_7_2020_9_38.png.png", width = "100%"),
                                                 href = "https://www.grinpatias.org/",
                                                 target = '_blank')))),
                                           
                                           column(3, 
                                                  div(width = 12,panel(heading = "Brasil", style = "height: 280px", width = 3, status = "primary",
                                                 br(), br(), br(),
                                                 shiny::a(
                                                 div(img(src = "grin_brasil.PNG", width = "100%")),
                                                 href = "https://grinbr.com/",
                                                 target = '_blank',
                                                 )))),
                                           column(12, div("...and many more! Please", shiny::a("contact us", href="mailto:GRIN@medizin.uni-leipzig.de") ,"if you would like your organization to be listed here.", align = "center"))
                                       ))))
                                       ))),
                 
                 ), # end families tab
        
        ##### VARIANT ANALYSIS #####
        
        tabPanel(title = "Variant Analysis", value = "variantTab",
                 fluidRow(column(10, offset = 1,
                                 fluidRow(panel(status="primary", heading = "Enter variant",
                                                column(2,
                                                       pickerInput(
                                                           inputId = 'var_gene',
                                                           label =  h5(strong('GRIN gene')),
                                                           width = "100%",
                                                           choices = c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2D"),
                                                           options = list(`style` = "default")
                                                       )
                                                ),
                                                column(2,
                                                       radioButtons(
                                                           inputId = "var_origin",
                                                           label = h5(strong("Origin")), 
                                                           choices = c("de novo (confirmed)", "de novo (assumed)", "Inherited", "Unknown/Untested"),
                                                           selected = "Unknown/Untested"
                                                       )
                                                ),
                                                column(2,
                                                       radioButtons(
                                                           inputId = "var_phenotype",
                                                           label = h5(strong("Phenotype")),
                                                           choices = c("Neurodevelopmental Disorder (NDD) + Cortical Malformation", "NDD + Epilepsy or Cerebral Visual Impairment", 
                                                                       "Unspecified NDD", "Other"),
                                                           selected = "Other"
                                                       )),
                                                column(3,
                                                       numericInputIcon(
                                                           inputId = "search_c_pos",
                                                           label = h5(strong("cDNA Position")),
                                                           min = 1,
                                                           max = 10000,
                                                           value = 1949,
                                                           width = "100%"
                                                       ),
                                                       pickerInput(
                                                           inputId = "search_ref_c",
                                                           label = "Ref",
                                                           choices = c("G", "A", "C", "T"), 
                                                           selected = "A"
                                                       ),
                                                       pickerInput(
                                                           inputId = "search_alt_c",
                                                           label = "Alt",
                                                           choices = c("G", "A", "C", "T", "null"), 
                                                           selected = "T"
                                                       ),
                                                       actionButton(inputId = "search_var_c", label = "Search")),
                                                column(3, 
                                                       numericInputIcon(
                                                           inputId = "search_aa_pos",
                                                           label = h5(strong("Amino Acid Position")),
                                                           min = 1,
                                                           max = 1484,
                                                           value = 650,
                                                           width = "100%"
                                                       ),
                                                       pickerInput(
                                                           inputId = "search_ref_p",
                                                           label = "Ref",
                                                           choices = sort(unique(GRIN_clin$ref_p)), 
                                                           selected = "Asn"
                                                       ),
                                                       pickerInput(
                                                           inputId = "search_alt_p",
                                                           label = "Alt",
                                                           choices = sort(unique(GRIN_clin$alt_p)), 
                                                           selected = "Ile"
                                                       ),
                                                       actionButton(inputId = "search_var_p", label = "Search"))
                                 )))),
                 
                 fluidRow(
                   column(10, offset = 1, 
                          fluidRow(
                            panel(
                              status="default", heading = "Variant Information",
                              tabsetPanel(
                                tabPanel("Variant Information",br(),
                                         fluidRow(
                                           div(width = "100%",valueBoxOutput("grinBox1")),
                                           div(width = "100%",valueBoxOutput("grinBox2")),
                                           div(width = "100%",valueBoxOutput("grinBox3")))),
                                tabPanel("ACMG Criteria (Preliminary Use)",br(),
                                         tags$style(
                                           HTML(".tabbable > .nav > li > a[data-value='ACMG Criteria (Preliminary Use)']{background-color: #fcedec;  color:black}")),
                                         fluidRow(
                                           column(10, offset =1, 
                                           box(title = div("Clinical Significance", style ="color:white;"), icon = icon("medkit"), 
                                               solidHeader=TRUE, status = "danger", width = 12, footer = div(alert(status = "danger", strong("Preliminary/incomplete use."), "Criteria for GRIN Registry variants 
                                               were automatically 
                                               annotated according to standard ", shiny::a(href="https://www.acmg.net/docs/Standards_Guidelines_for_the_Interpretation_of_Sequence_Variants.pdf", 
                                                                                 "ACMG Variant Interpretation Guidelines") , "(GRIN1, GRIN2A, GRIN2D) or the ClinGen GRIN Expert Panel 
                                                                                 Specifications to the ACMG/AMP Variant Interpretation Guidelines (GRIN2B, unpublished). They are yet incomplete and not to be 
                                                                                 used in a medical context.", style = sub_style)),
                                                 uiOutput(outputId = "ACMGUI_crit"),
                                               box(title="Explanation", collapsible = TRUE, collapsed = TRUE, width = 12,
                                                   p("The following criteria were automatically annotated:"),
                                                   DTOutput(outputId = "ACMGUI_exp"), br(),
                                                   p("Evidence specifications: VS, very strong; 
                                                     S, strong; M, moderate; P, supporting", style = sub_style), br(),
                                                   p(strong("PVS1"), "- Null variant in a gene where loss of function is a known mechanism of disease. Very Strong evidence."),
                                                   p(strong("PS1"), "- Same amino acid change as a previously established pathogenic variant regardless of nucleotide change. Strong - supporting evidence."),
                                                   p(strong("PS2"), "- De novo (both maternity and paternity confirmed) in a patient with the disease and no family history. Very strong - supporting evidence.
                                                     For specifications see", shiny::a("here.", href="https://www.medizinische-genetik.de/fileadmin/bilder_grafiken_diagramme/allgemeine_informationen/variantenklassifikation/PS2_PM6.jpg")),
                                                   p(strong("PS3"), "- Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product. Strong - moderate evidence."),
                                                   p(strong("PS4"), "- The same variant has been previously identified in two or more/one unrelated affected individuals and has not been reported in gnomAD. Moderate - supporting evidence."),
                                                   p(strong("PM1"), "- Located in an “intolerant 3D microdomain” without  benign variation. Right now only available for GRIN2B. Moderate - supporting evidence."),
                                                   p(strong("PM2"), "-	Completely absent from all population databases (at least from gnomAD). Supporting evidence."),
                                                   p(strong("PM5"), "- Missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before. Moderate - supporting evidence."),
                                                   p(strong("PM6"), "- Assumed de novo, but without confirmation of paternity and maternity. Very strong - supporting evidence.
                                                     For specifications see", shiny::a("here.", href="https://www.medizinische-genetik.de/fileadmin/bilder_grafiken_diagramme/allgemeine_informationen/variantenklassifikation/PS2_PM6.jpg")),
                                                   p(strong("PP2"), "- Missense variant in a gene that has a low rate of benign missense variation and where missense variants are a common mechanism of disease. Supporting evidence."),
                                                   p(strong("PP3"), "- Multiple lines of computational evidence support a deleterious effect on the gene or gene product. Supporting evidence. ")
                                                   )))))))))), 
                 
                 fluidRow(
                   column(10, offset = 1, 
                          fluidRow(
                            panel(
                              status="default", heading = "Registry Information",
                              tabsetPanel(
                                tabPanel("Patient Information",br(), value = "PatientRegistryTab",
                                         box(title = "Individuals with the same variant in the GRIN Registry", width = 12,
                                         DT::dataTableOutput(outputId = "patientTable"),
                                         br(), p("Abbreviations: DD/ID, Developmental Delay/Intellectual Disability;
                                         ASD, Autism Spectrum Disorder; MD, Movement Disorder; MCD, Malformation of Cortical Development; 
                                         CVI, Cerebral Visual Impairment", style=sub_style, align = "center"))),
                                
                                tabPanel("Functional Information", value = "CFERVfunctionalTab",
                                         fluidRow(column(1),column(4, br(), align="justify",
                                                                   box(title = "", width = 12,plotlyOutput(outputId = "grinspider", height = "100%"),
                                                                    br(),
                                                                   
                                                                    footer = p("Functional data represent fold changes in receptor
                                                                      activity by gene variant compared to wild-type. 1 means no
                                                                      change from wild-type. 
                                                                      Fold_glu, Glutamate; Fold_gly, Glycine; Fold_mg, Magnesium; 
                                                                      Fold_zn, Zinc; Fold_ph, pH"
                                                                      ,style=sub_style, align = "justify"))), 
                                                  
                                                column(7, br(), br(), br(), panel(heading = "Functional Summary", align = "justify",
                                                                      
                                                                      h1(strong(textOutput(outputId = "spiderSources"), style="color:#D55E00;")),
                                                                      
                                                                      DT::dataTableOutput(outputId = "spiderTable"), br(),
                                                                      
                                                p(strong("Overview"), "The potential net effect of changes in NMDA receptor function assessed in vitro were categorized 
                                                  by the same numerical threshold used as supporting evidence for pathogenicity in ACMG criteria. These thresholds used 
                                                  the variability of functional measurements in variants observed in the general population to determine confidence intervals 
                                                  with which a patient-derived variant could be judged unlikely to occur by chance in the general population."), 
                                                
                                                p(strong("Likely Gain-of-Function (Supporting Pathogenic): "), "Changes that increase receptor activity for one of the following 
                                                  parameters (without offsetting changes) of >3-fold decrease in glutamate or glycine EC50 value, >2-fold increase in Mg2+ IC50 value, 
                                                  >2-fold increase in the weighted tau deactivation (tauW), >2-fold increase in open probability, >3-fold increase in surface expression, 
                                                  all assays with appropriate biological replicates compared to values determined with the same assay on the same day for wild-type controls."),
                                                
                                                p(strong("Likely Loss-of-Function (Supporting Pathogenic):"),"Changes that decrease receptor activity for one of the following 
                                                  parameters (without offsetting changes) of >3-fold increase in glutamate or glycine EC50 value, >2-fold decrease in Mg2+ IC50 
                                                  value, >2-fold decrease in deactivation tauW, >2-fold decrease in open probability, >3-fold decrease in surface expression, 
                                                  all assays with appropriate biological replicates compared to values determined with the same assay on the same day for wild-type controls."),
                                                
                                                p(strong("Potential Gain-of-Function (Supporting Pathogenic):"),"Changes in two or more parameters that increase receptor 
                                                  activity, including 1.5-3-fold decrease in glutamate or glycine EC50 value, 1.5-2-fold increase in Mg2+ IC50 value, 
                                                  1.5-2-fold increase in deactivation tauW, 1.25-2-fold increase in open probability, or 2-3-fold increase in surface expression, 
                                                  all assays with appropriate biological replicates compared to values determined with the same assay on the same day for wild-type controls"),
                                                
                                                p(strong("Potential Loss-of-Function (Supporting Pathogenic):"), "Changes in two or more parameters that decrease receptor 
                                                  activity, including 1.5-3-fold increase in glutamate or glycine EC50 value, 1.5-2-fold decrease in Mg2+ IC50 value, 1.5-2-fold 
                                                  decrease in deactivation tauW, 1.25-2-fold decrease in open probability, or 2-3-fold decrease in surface expression, all assays 
                                                  with appropriate biological replicates compared to values determined with the same assay on the same day for wild-type controls"),
                                                
                                                p(strong("Complex:"),"Change to more than one parameter in excess of any 
                                                  threshold described above with opposing effects on receptor activity."),
                                                
                                                p(strong("Insufficient Data:"), "Currently there is insufficient data to support a pathogenic 
                                                  classification of this variant, although some parameters appear changed."),
                                                
                                                p(strong("No Support for Pathogenicity:"), "Currently the data do not support a pathogenic classification of this variant.")
                                                       ))
                                                  )),
                                
                                tabPanel("Paralogous Pathogenic Variants", value = "ParalogousPathoTab",
                                         DT::dataTableOutput(outputId = "paraTable"),
                                         br(), p("Abbreviations: DD/ID, Developmental Delay/Intellectual Disability;
                                         ASD, Autism Spectrum Disorder; MD, Movement Disorder; MCD, Malformation of Cortical Development; 
                                         CVI, Cerebral Visual Impairment", style=sub_style, align = "center"))))))),
                 
                 fluidRow(
                   column(10, offset = 1, 
                          fluidRow(
                            panel(
                              status="default", heading = "Comparative Information",
                              strong("Compare the selected variant with other similar variants."),
                              br(),
                              radioGroupButtons(
                                inputId = "compareButtons",
                                label = "Variants with the same:",
                                choices = c("Amino Acid Position", "Domain", #"Functional Consequence", 
                                            "Variant Type"),
                                justified = TRUE,
                                status = "default",
                                checkIcon = list(yes = icon("ok",lib = "glyphicon"))
                              ), 
                              materialSwitch(inputId = "hide_unknown", label = "Hide missing information", status = "default", right=T, inline=T),
                              materialSwitch(inputId = "same_GRIN", label = "Hide other GRIN genes", status = "default", right=T, inline=T),
                              br(),
                              div(width = "100%", plotlyOutput("comparePlot")),
                              fluidRow(
                              column(12, align="center", p("Abbreviations: DD/ID, Developmental Delay / Intellectual Disability", style=sub_style))),
                              DT::dataTableOutput(outputId = "compareTable")
                            ),
                           )))
                 ), # end variant analysis tab
        
        ##### RESEARCH #####
        
        tabPanel(title = "Research", value = "researchTab",
                 fluidRow(
                     column(10, offset=1,
                            panel(heading="Filter GRIN Registry", status="primary",
                                  fluidRow(
                                      column(12,
                                             selectizeGroupUI(
                                                 id = "grin-filters",
                                                 btn_label = "Reset filters",
                                                 params = list(
                                                     gene = list(
                                                         inputId = "gene",
                                                         placeholder="all",
                                                         title = strong("Gene"),
                                                         choices = c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2D"),
                                                         multiple = TRUE
                                                     ),
                                                     varcons = list(
                                                         inputId = "var_type",
                                                         title = strong("Variant Type"),
                                                         placeholder="all",
                                                         choices = c("missense", "null")
                                                     ),
                                                     aachange = list(
                                                         inputId = "alt_p",
                                                         title = strong("Amino Acid Change"),
                                                         placeholder="all",
                                                         choices = unique(GRIN_clin$alt_p)
                                                     ),
                                                     domain = list(
                                                         inputId = "domain",
                                                         title = strong("Domain"),
                                                         placeholder="all",
                                                         choices = unique(GRIN_clin$domain)
                                                     ),
                                                     functional_cons = list(
                                                       inputId = "CALL",
                                                       title = strong("Functional Consequence"),
                                                       placeholder="all",
                                                       choices = c("A", "B"),
                                                       multiple = TRUE
                                                     ),
                                                     phenotype = list(
                                                         inputId = "delay_degree",
                                                         title = strong("DD/ID"),
                                                         placeholder="all",
                                                         choices = c("A", "B"),
                                                         options = list(`selected-text-format` = "count > 3"),
                                                         multiple = TRUE
                                                     ),
                                                     seizures = list(
                                                         inputId = "seizures",
                                                         title = strong("Seizures"),
                                                         placeholder="all",
                                                         choices = c("A", "B"),
                                                         multiple = TRUE
                                                     )
                                                     ))
                                             ))))),
                 fluidRow(

                   column(10, offset = 1,
                          panel(heading = div(strong("Filtered Subset"), textOutput(outputId = "filtered_n")),
                    
                 tabsetPanel(
                   tabPanel(vvalue = "GenotypeInfoLocation",
                     "Genotype Interface",
                     fluidRow(
                       column(
                         8,
                         br(),
                         p(
                           strong("Selected variants are displayed in 2D (lolliplot) and 3D (protein structure).")
                         ),
                         
                         fluidRow(column(
                           7,
                           materialSwitch(
                             inputId = "gnomad_m",
                             label = "Reference population missense variants (gnomAD)",
                             status = "primary",
                             right = T,
                             inline = T
                           )
                           
                         ),
                         column(5,plotOutput(outputId = "legend_plot", height = 30, width = 450) ),
                         #column(
                         #  4,
                         #  materialSwitch(
                         #    inputId = "gnomad_plof",
                         #    label = "gnomAD pLoF variants",
                         #    status = "primary",
                         #    right = T,
                         #    inline = T
                         #  )
                         #)
                         ),
                         addSpinner(plotlyOutput(outputId = "GRIN_overview_plot"), color =
                                      spinner_color),
                         p("The following transcripts were used: GRIN1: NM_007327.4, 
                           GRIN2A: NM_001134407.3, GRIN2B: NM_000834.4, 
                           GRIN2D: NM_000836.2", align="center", style=sub_style)
                       ), 
                       
                       column(
                         4,
                         
                         tags$style(
                           HTML(".tabbable > .nav > li > a[data-value='GRIN1']{background-color: #edf4ed;  color:black}"),
                           HTML(".tabbable > .nav > li > a[data-value='GRIN2A']{background-color: #edf4ed;  color:black}"),
                           HTML(".tabbable > .nav > li > a[data-value='GRIN2B']{background-color: #edf4ed;  color:black}"),
                           HTML(".tabbable > .nav > li > a[data-value='GRIN2D']{background-color: #edf4ed;  color:black}")
                         ),
                         tabsetPanel(
                           type = "pills",
                           tabPanel("GRIN1",
                                    p(fluidRow(column(10, 
                                                      materialSwitch("map_3d_1", label="GRIN1 and GRIN2A Complex Structure", value=TRUE, status="primary", right=TRUE)),
                                               column(2, p(tippy(icon("question-circle"), 
                                                                 tooltip = "PDB Structure: 6IRA (GRIN1 and GRIN2A)",
                                                                 animation = "scale", theme = "light")))
                                               ),
                                      
                                      addSpinner(
                                        r3dmolOutput(
                                          outputId = "threeDmolGRIN1",
                                          width = "100%",
                                          height = "400px"
                                        ),
                                        color = spinner_color, 
                                      ),
                                      div("UniProt GRIN1:", align="center", style=sub_style, shiny::a("Q05586", href="https://www.uniprot.org/uniprot/Q05586", target="_blank")
                                      #  shiny::a("PDB", href="https://www.rcsb.org/structure/6IRA", target="_blank"))
                                      ))),
                           tabPanel("GRIN2A",
                                    p(fluidRow(column(10, 
                                                      materialSwitch("map_3d_2A", label="GRIN1 and GRIN2A Complex Structure", value=TRUE, status="primary", right=TRUE)),
                                               column(2, p(tippy(icon("question-circle"), 
                                                                 tooltip = "PDB Structure: 6IRA (GRIN1 and GRIN2A)",
                                                                 animation = "scale", theme = "light")))),
                                      
                                      addSpinner(
                                        r3dmolOutput(
                                          outputId = "threeDmolGRIN2A",
                                          width = "100%",
                                          height = "400px"
                                        ),
                                        color = spinner_color,
                                      ),
                                      div("UniProt GRIN2A:", align="center", style=sub_style, 
                                          shiny::a("Q12879", href="https://www.uniprot.org/uniprot/Q12879", target="_blank")))),
                                    
                           tabPanel("GRIN2B",
                                    addSpinner(color = spinner_color,
                                               r3dmolOutput(
                                                 outputId = "threeDmolGRIN2B",
                                                 width = "100%",
                                                 height = "400px"
                                               )),
                                               div(br(), br(), "UniProt GRIN2B:", align="center", style=sub_style, 
                                                   shiny::a("Q13224", href="https://www.uniprot.org/uniprot/Q13224", target="_blank"))),

                           tabPanel("GRIN2D",
                                    addSpinner(color = spinner_color,
                                               r3dmolOutput(
                                                 outputId = "threeDmolGRIN2D",
                                                 width = "100%",
                                                 height = "400px"
                                               )),
                                    div(br(), br(), "UniProt GRIN2D:", align="center", style=sub_style, 
                                        shiny::a("O15399", href="https://www.uniprot.org/uniprot/O15399", target="_blank")))
                           
                      
                     
                   )))),
                   
                   tabPanel("Phenotype Interface",br(),
                            fluidRow(
                              column(4,align="justify", plotlyOutput("phenotype_number_plot")
                              ),
                              column(4, plotlyOutput("violin_grin"), align="center",
                                     p("0 - no ID; 1 - mild ID; 2 - moderate ID; 3 - severe/profound ID", style=sub_style)),
                              column(4, plotlyOutput("seizure_gene"))
                            )
                   ),
                   
                   tabPanel("Functional Interface", br(),
                            fluidRow(
                              column(4,align="justify",plotlyOutput("functional_number_plot")
                              ),column(4,
                              tabsetPanel(
                                tabPanel(title = "Glu", br(),
                              plotlyOutput("functional_Glu")),

                              tabPanel(title = "Gly", br(),
                              plotlyOutput("functional_Gly")),
                              tabPanel(title = "Mg", br(),
                              plotlyOutput("functional_Mg")),
                              tabPanel(title = "Zn", br(),
                              plotlyOutput("functional_Zn")),
                              tabPanel(title = "tauW", br(),
                                       plotlyOutput("functional_tauW")),
                              tabPanel(title = "Surface Expr.",  br(),
                                       plotlyOutput("functional_surface_mean"))
                            ),
                            
                            p("Functional data represent fold changes in receptor activity by gene variant compared to wild-type. Zero means no
                              change from wild-type. Logarithmic scaling was used for the plot.",style=sub_style)),
                            
                            column(4, plotlyOutput("functional_cons_plot"), 
                                   br(), br()),
                            column(4, 
                                   
                                   fluidRow(column(12, align="center", "All data from the literature and: ", 
                                   div(width="100%", p(shiny::a(img(src="CFERV-banner.png", width ="60%"), 
                                              href="http://functionalvariants.emory.edu/index.html",
                                              target = '_blank'
                                              )))),
                                   column(4),
                            column(4,br(), align = "center",
                                   actionButton(
                                     inputId = "info",
                                     label = "Methods",
                                     icon = icon("info"),
                                     
                                     style="color: #fff; background-color: #2494f3; border-color: #1e89de"
                                   ))))
                   )),
                   tabPanel("Animal Models", br(),
                            DT::dataTableOutput(outputId = "animalTable"), br(),
                            p("Filtering is not yet available. A more comprehensive list will be updated soon.", style = sub_style)
                   )
                 ), 
                 br(),
                 panel(heading="Displayed Variants",
                       materialSwitch(
                         inputId = "patientFunSwitch",
                         label = "Show Functional Data",
                         status = "primary",
                         right = T,
                         inline = T
                       ),
                       DT::dataTableOutput(outputId = "subsetTable")
                 ))
        ))), # end research tab
        
        ##### GRIN REGISTRY #####
        
        tabPanel(title = "GRIN Registry", value = "registryTab",
                 fluidRow(column(10, offset=1,
                                 panel(heading="How to enter the GRIN Registry",status="primary",
                                       fluidRow(
                                           panel(column(8, panel(heading=
                                                                     "If you are from Europe, Asia, Africa or other",
                                                                 "This chapter of the study is led by geneticist Dr. Johannes Lemke (University of Leipzig). Please fill out the online consent
                                              form and continue with the online clinical questionnaire; you can ask your doctor to help fill it out.
                                              Please keep your return code and return link to re-access your entry at a later time.
                                              The Leipzig team will contact CFERV to conduct functional studies for novel variants entered in the registry.
                                              In case we have questions on your entry, you may be contacted by an administrator.
                                              In case you have questions, please contact GRIN@medizin.uni-leipzig.de."
                                           )),
                                           column(4,align="center", img(src="grin-world-2.PNG", width="95%"), 
                                                  actionBttn(
                                                      inputId="link2",
                                                      label = shiny::a("Participate", href="https://redcap.medizin.uni-leipzig.de/redcap/surveys/?s=PRAEF9N7J7", style="color: white;", target =
                                                                           "_blank"),
                                                      style = "simple",
                                                      color = "warning",
                                                      icon = icon("arrow-right"),
                                                      size = "md",
                                                      block = FALSE,
                                                      no_outline = TRUE
                                                  )))), 
                                       
                                       fluidRow(
                                           panel(
                                               column(8,
                                                      panel(
                                                          heading = "If you are from North or South America or Australia",
                                                          "please contact Jenifer Sargent (study coordinator) at:
                                              Jenifer.Sargent@cuanschutz.edu. This chapter of the study is led by neurologists Drs. Tim Benke and Kristen Park
                                              (University of Colorado) and Dr. Jennifer Bain (", shiny::a("Simons Foundation/Columbia University", href = "https://www.simonssearchlight.org/"
                                              ),")",
                                              
                                              br(),
                                              br(),
                                              "1. Email Jenifer Sargent your de-identified genetic testing results to make sure you qualify.",
                                              br(),
                                              "2. After we review your genetic testing results, Jenifer will send you a consent form; sign and return to her.",
                                              br(),
                                              "3. Jenifer will send you the clinical questionnaire. We are updating the questionnaire and want you to fill out the latest version, even if you've already done it before.",
                                              br(),
                                              "4. Fill out the questionnaire; you can ask your doctor to help fill it out.",
                                              br(),
                                              "5. Return the filled-out questionnaire to Jenifer; Jenifer will contact CFERV to conduct functional studies for novel variants.",
                                              br(),
                                              "6. In the near future, Jenifer will email you for yearly updates. If you haven't heard from us, then please contact us again."
                                                      ),
                                              panel(
                                                p("Both chapters of the registry will be merged, and patients only need to be enrolled at Colorado or Leipzig.
                                              If you have already enrolled but have additional novel information, please get in contact and update your entry."
                                                )
                                              )
                                               ),
                                              column(4,
                                                     align="center", img(src="grin-world-1.PNG", width="95%"),
                                                     actionBttn(
                                                         inputId="link1",
                                                         label = shiny::a("Participate", 
                                                                          href="mailto:Jenifer.Sargent@cuanschutz.edu?subject=GRIN%20Registry&body=Dear%20Jenifer,%0D%0A%0D%0AWe%20would%20like%20to%20participate%20in%20the%20GRIN%20Registry.%20Please%20let%20us%20know%20how%20we%20can%20contribute.", 
                                                                          style="color:white;"
                                                         ),
                                                         icon = icon("arrow-right"),
                                                         style = "simple",
                                                         color = "warning",
                                                         size = "md",
                                                         block = F,
                                                         no_outline = TRUE
                                                     )))),
                                       
                                 ))),
                 
        ),
        
        ##### ABOUT #####
        
        tabPanel(title = "About", value = "aboutTab", 
                 fluidRow(
                     column(10, offset = 1,
                            br(),
                            tabsetPanel(
                              tabPanel("General Information",
                            panel(heading ="About", status = "primary",
                                  p(strong("This is the alpha version of the  GRIN Portal.")),
                                  p("The GRIN Portal is a coalition of investigators seeking to aggregate and harmonize data generated to study 
                                    GRIN-related disorders, and to make summary data interactively accessible for the wider scientific community, 
                                    while providing educational resources for everyone. The goals of this project are: "),br(),
                                  fluidRow(column(6, 
                                  p("- Providing information on GRIN-related disorders"),
                                  p("- Supporting research on GRIN-related disorders"),
                                  p("- Facilitating recruitment of individuals to the global GRIN registry")),
                                  
                                  column(6,
                                         p("- Providing support in variant interpretation and classification"),
                                         p("- Visualizing data from the global GRIN registry"),
                                         p("- Linking researchers, clinicians and families"))),br(),
                                  footer = div("The project is overseen by Johannes Lemke (contact PI), Steve Traynelis (contact PI), Tim Benke and Dennis Lal. 
                                    This GRIN Portal is an ongoing project and interested collaborators are invited to reach out to join the project.")),
                panel(heading = "Teams and People", status = "primary",
                      p("The current version of the GRIN Portal has been developed by an international team of researchers and clinicians:"),
                      fluidRow(#align = "center",
                        
                  column(12, panel(heading = "Team Leaders", 
                                  p(strong("Johannes Lemke"), "(Leipzig, Germany): General concept, clinical & genetic data"), 
                                  p(strong("Tim Benke"), "(Denver, US): Clinical & genetic data"), 
                                  p(strong("Steve Traynelis"), "(Atlanta, US): Molecular data"),
                                  # p(strong("Amy Ramsey"), "(Toronto, Canada): Animal data"),
                                  p(strong("Dennis Lal"), "(Cleveland, US): General concept, web development, bioinformatics, video production"),
                                  )),
                  
                  column(2, panel(style="height: 230px;",heading = "Clinical & Genetic Data",
                                  div(style="height: 100%;",
                                      p(strong("Ilona Krey")),
                                  p(strong("Jenifer Sargent")),
                                  p("Chiara Klöckner"),
                                  p("Vincent Strehlow"),
                                  p("Konrad Platzer"))
                                  )),
                  
                  column(2, panel(style="height: 230px;",heading = "Molecular Data",
                                  div(style="height: 100%",p(strong("Scott Myers")),
                                  p(strong("Hongjie Yuan")))
                  )),
                  
                  column(2, panel(style="height: 230px;",heading = "Animal Data",
                                  div(style="height: 100%",p(strong("Amy Ramsey")))
                  )),
                  
                  column(2, panel(style="height: 230px;",heading = "Web Development",
                                  div(style="height: 100%",p(strong("Chiara Klöckner")),
                                  p("Marie Mcnee"),
                                  p("Tobias Brünger"),
                                  p("Eduardo Perez-Palma"))
                  )),
                  
                  
                  column(2, panel(style="height: 230px;",heading = "Bioinformatics",
                                  div(style="height: 100%",p(strong("Chiara Klöckner")),
                                  p(strong("Tobias Brünger")),
                                  p(strong("Eduardo Perez-Palma")),
                                  p("Marie Mcnee"),
                                  p("Patrick May"))
                  )),
                  
                  column(2, panel(style="height: 230px;", heading = "Video",
                                  div(style="height: 100%",p(strong("Arthur Stefanski")),
                                  p("Chiara Klöckner"),
                                  p("Tobias Brünger"),
                                  p("Marie Mcnee"))
                  )))),
                  
                
                panel(heading = "Impressum", status = "primary", p("We object to any commercial use and disclosure of data."),
                      p(strong("Copyright and use:"), "The authors grants you the right of use to make a private copy for personal purposes. 
                        However, you are not entitled to change and/or pass on the materials or publish them yourself.
                        Upon request, you will receive free information about the personal data stored about you. 
                        To do so, please contact the administrator."),
                      p(strong("No liability:"), "The contents of this web project have been carefully checked and created to 
                        the best of our knowledge. But for the information presented here is no claim to completeness, 
                        timeliness, quality and accuracy. No responsibility can be accepted for any damage caused by reliance 
                        on or use of the contents of this website."))
                     ),
                tabPanel("Terms and Data Information",
                         panel(heading = "Terms of Use", status = "primary",
                               p("All data here are publicly for the benefit of the wider biomedical community. 
                               You can freely explore the data, and we encourage the use and publication of results generated 
                               from these data. However, we encourage you to contact us before embarking on analyses to 
                               check if your proposed analysis overlaps with work currently underway by our team. Further, 
                               we request that you actively acknowledge and give attribution to the GRIN Portal project, and 
                               link back to the relevant page, wherever possible. All users of our data agree to not attempt 
                               to reidentify participants. Our data set has been subjected to extensive quality control, 
                               but may be imperfect so errors may remain."),
                               p("If you spot any results that seem impossible, 
                               or have suggestions for GRIN Portal improvements: ", 
                               shiny::a(href="mailto:GRIN@medizin.uni-leipzig.de; chiara.kloeckner@medizin.uni-leipzig.de", "Contact us"), 
                               "that we can improve.")),
                         panel(heading = "Data Generation", status = "primary",
                               "A full description of the methods used to aggregate and generated in this project will be provided shortly. ")))
                 )
        ))
        

)) # end ui
