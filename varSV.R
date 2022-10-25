library(shiny)
library(dplyr)
library(data.table)


############################################
################ FLUID PAGE ################
############################################

ui <- fluidPage(
      tags$head(
        tags$style(type="text/css", "select { max-width: 270px; }"),
        tags$style(type="text/css", ".span4 { max-width: 270px; }"),
        tags$style(type="text/css", ".well { max-width: 270px; }")
      ),
      
  div(
    style = "display:flex; align-items:flex-start",
    wellPanel( # sidebar
      
  #### INPUT() FUNCTIONS ####
  fileInput(inputId="file1", label="var hit file", multiple = FALSE, width = NULL,
            buttonLabel = "Browse...", placeholder = "No file selected",
            accept = c('text/csv','text/comma-separated-values',
                       'text/tab-separated-values', 'text/plain',
                       '.csv', '.tsv')
            ),

  # region input 
  selectInput(inputId='region',label='Select var region',
              choices=c('chr1_telo_1','chr1_telo_2','chr2_telo_1','chr2_telo_2','chr3_telo_1','chr3_internal_0','chr3_telo_2',
                        'chr4_telo_1','chr4_internal_0', 'chr4_internal_1','chr4_internal_2','chr4_telo_2',
                        'chr5_telo_1','chr5_telo_2','chr6_telo_1','chr6_internal_1','chr6_telo_2',
                        'chr7_telo_1','chr7_internal_1','chr7_telo_2',
                        'chr8_telo_1','chr8_internal_1',  'chr8_telo_2','chr9_telo_1','chr9_telo_2',
                        'chr10_telo_1','chr10_telo_2','chr11_telo_1','chr11_telo_2',
                        'chr12_telo_1','chr12_internal_1','chr12_internal_2','chr12_telo_2',
                        'chr13_telo_1','chr13_telo_2'
                        )
              ),

  # slider input: % ID
  sliderInput(inputId = "pid",
              label = "% Identity for THIS READ",
              value = 90, min = 60, max = 100
  ),
  # buttons to advance reads 
  actionButton("button1", "Previous read"),
  actionButton("button2", "Next read"),
  numericInput("N_read", "Select a read", 1),
  actionButton("go", "Go"),

  # text input: x-axis values 
  textInput("xmin", "X-axis min"),
  textInput("xmax", "X-axis max")

  ),
  
  downloadButton("save", "save"),
  
  #### OUTPUT FUNCTIONS ####
  # * Output() adds a space in the ui for an R object; you must build the object in the server function below
   div( # main panel
     style = "flex-grow:1; resize:horizontal; overflow: hidden",
  # sample plot (experimenting)
  plotOutput("plot", width = "100%", height = "650px"),
  textOutput("N_reads"),
  textOutput("N_hits"),
  textOutput("counter"),
  textOutput('text1'),
  # show the uploaded data in a table 
  tableOutput("table1"),
  tableOutput("table2")
)
)
)


############################################
############## SERVER FUNCTION ##############
############################################

server <- function(input, output) {
  options(shiny.maxRequestSize=4*800*1024^2) # max size in MB of file upload

  # save the .hit data as a reactive to manipulate 
  data <- eventReactive(input$file1, {read.csv(input$file1$datapath, sep='\t',stringsAsFactors = FALSE)
  }) 

  
  ###########################
  #### OTHER INPUT FILES ###
  ###########################
  
  setwd("/Users/emily/Desktop/forgithub/example-chrom12-ANC")
  
  ### Colors and levels for each gene and intergenic sequence 
  segcolors <- read.table("segcol_all_070722.csv", header = TRUE, sep = ",",comment.char = "@",stringsAsFactors = FALSE) # the '#' symbol is interpreted as a comment!
  
  ### Assignments of reads/contigs to regions (specific to sample; be sure to update)
  read.regions <-   read.table("ANC-REF.Chrom-RegionAssignments.csv", header = T, sep = ",", stringsAsFactors = F)
  
  ### Read lengths (specific to sample; be sure to update)
  readlens <- read.table("ANC.30kb.readlengths.txt", header = FALSE, sep = "\t",comment.char = "@",stringsAsFactors = FALSE) 
  
  
  ###########################
   #### HIT FILTERING ####
  ###########################
  
  
  ### Remove hits for segments that aren't in this region & Get colors/levels for segments in this region as a reactive
  region_colors <- reactive({
      a <- subset(segcolors, segcolors$home_region == input$region)
      return(a)
    })
    
  region_segments <- reactive({   
    a <- unique(region_colors()$segments) 
    return(a)
  })

  hits.filt <- reactive({
    input.region.reads <- read.regions %>% filter (assigned_region == input$region)
    input.region.hits <- data()  %>% filter (subject.id %in% input.region.reads$read)
    return(input.region.hits)
  })
  
  # make a reactive with the names of all the reads in this region
  region_read_names <- reactive({
    a<-unique(hits.filt()$subject.id)
    return(a)
  })
  
  ## For pass reads, Insert a hard filter requiring 5 hits with length >200 and identity > 70
  read_names <- reactive({
    
    passreadnames = c()
    for(i in region_read_names()) {
        temp <- hits.filt()  %>% filter (subject.id == i)
        temp2 <- subset(temp, temp$identity > 75) #75 default
        temp3 <- subset(temp2, temp2$len.alignment > 200) # 200 defaul
        if (dim(temp3)[[1]]>=5){ passreadnames=c(passreadnames,i)}  } # 5 passing segments 
    return(passreadnames)
  })
  

  # print how many reads there are for the whole region/without filtering
  output$N_reads <- renderText({ 
    a <- unique(region_read_names())
    paste('Total reads: ',as.character(length(a)))
  })
  
  # print the number of  PASSING reads for the region/identity setting
  output$N_hits <- renderText({
    paste('Pass reads: ',as.character(length(read_names())))
  })

  ################################
  #### CLICKING THROUGH READS ####
  ################################
 
  # Advance to the next read with the press of a button:
  # the value in vars$counter starts at 1 
  vars<-reactiveValues()
  vars = reactiveValues(counter = 1)
  # add one by clicking "next read"
  observe({
    if(!is.null(input$button2)){
      input$button2
      isolate({
        vars$counter <- vars$counter + 1
      })    }  })
  # Subtract one by clicking "Previous read"
  observe({
    if(!is.null(input$button1)){
      input$button1
      isolate({
        vars$counter <- vars$counter - 1
      })    }  })
  # Jump to first read
  observe({
    if(!is.null(input$button3)){
      input$button1
      isolate({
        vars$counter <- 1
      })    }  })
  # Jump to designated read
  observe({
    if(!is.null(input$go)){
      input$go
      isolate({
        vars$counter <- input$N_read
      })    }  })
  output$counter <- renderText({
    paste("Current pass read: ",as.character(vars$counter))
  })

  
  ###################################
  #### GET DATA FOR CURRENT READ ####
  ##################################
  
  ### make plot_data reactive. CONTAINS HITS WHERE subject.id IS THE CURRENT PASS READ ## AND ## query.id is in input$region 
  plot_data <- reactive({
    temp = subset(hits.filt(),hits.filt()$subject.id==read_names()[vars$counter])
    temp2 <- subset(temp,temp$identity >= input$pid) # filter by identity 
    return(temp2)
  })

  ## print read length
  output$text1 <- renderText({
    a <- subset(readlens, readlens$V1 == plot_data()$subject.id)
    paste('Read length: ',as.character(a$V2))
  })
  

  
  ########################
  #### PLOT THIS READ ####
  #########################

  output$plot <-
    renderPlot({
      # get the levels/colors as a nonreactive 
      col_data <- region_colors()
      
      # Set xmin and xmax based on read data, unless there is a text input
      if(!isTruthy(input$xmin)){
        xmin=min(range(plot_data()$s.start),range(plot_data()$s.end))
      } else {
        xmin=as.numeric(input$xmin)
      }
      if(!isTruthy(input$xmax)){
        xmax=max(range(plot_data()$s.start),range(plot_data()$s.end))
      } else {
        xmax=as.numeric(input$xmax)
      }      
      if(!isTruthy(input$ymin)){
        ymin=0
      } else {
        ymin=as.numeric(input$ymin)
      }
      if(!isTruthy(input$ymax)){
        ymax=as.numeric(max(col_data$level))
      } else {
        ymax=as.numeric(input$ymax)
      }      
      
      #####   set up blank plot   
      par(mar = c(5,10,4,2) + 0.1)
      plot(0,0,col='white',xlim=c(xmin-1000,xmax+1000),ylim=c(ymin,ymax),xlab='position',yaxt='n',ylab='',main=read_names()[vars$counter])
      
      ####  determine size of y-axis.
      # need a list of only main (gene) segments, not intergenic.
      # go through all the segments in col_data, removing 'post-', and make the y axis equal to number of genes 
      segpost = c()
      allsegments = col_data$segments
      for (s in 1:length(allsegments)){
        seg =  allsegments[[s]]
        z = strsplit(seg, '-',fixed=TRUE)
        if (length(z[[1]]) > 1){ segpost=c(segpost,z[[1]][-1])
            } else { segpost=c(segpost,z[[1]]) }
      }
      segpost=unique(segpost) 
      axis(2, at=1:length(segpost), labels=segpost, las=2)
      
      #### add one segment to the plot for each row of plot_data
      for (hit in 1:dim(plot_data())[[1]]){
        row = plot_data()[hit,]
        if ( grepl('of', row[1], fixed = TRUE) == TRUE) { # split it
          seg = strsplit(row[[1]], '.',fixed=TRUE)[[1]][[1]] } else{seg=row[[1]]}
        segcolrow = subset(col_data,col_data$segments==seg)
        
        # avoid error for segments w/o colors b/c they aren't in this region 
        if(dim(segcolrow)[[1]]>0) {
                plotcol = as.character(segcolrow[[2]])
                yval = as.numeric(segcolrow[[3]])
                # find start and end points
                xpos1 = min(as.numeric(row[9]),as.numeric(row[10]))
                xpos2 = max(as.numeric(row[9]),as.numeric(row[10]))
                # actually plot
                rect(xleft=xpos1, ybottom=yval-.49, xright=xpos2, ytop=yval+.49, density = -1, angle = 45, col = plotcol,border='black')
                
                ### add a text label for each segment (e.g. 1of15) -- optional 
              #   if ( grepl('.',row$query.id, fixed = TRUE) == TRUE) { # split it
              #    labeltext = strsplit(row$query.id,'.',fixed=TRUE)[[1]][[2]] } else{labeltext='1of1'}
              #   text(x=mean(xpos1,xpos2),y=yval-0.7,labels=labeltext,cex=1,srt=90)
       }} # end hit loop 
      
      
      
      ## recreate the plot as an object to download
      p2 <- reactive({
        par(mar = c(5,10,4,2) + 0.1)
        plot(0,0,col='white',xlim=c(xmin-1000,xmax+1000),ylim=c(ymin,ymax+1),xlab='position',yaxt='n',ylab='',main=read_names()[vars$counter])
        
        for (hit in 1:dim(plot_data())[[1]]){
            row = plot_data()[hit,]
            if ( grepl('of', row[1], fixed = TRUE) == TRUE) { # split it
              seg = strsplit(row[[1]], '.',fixed=TRUE)[[1]][[1]] } else{seg=row[[1]]}
            segcolrow = subset(col_data,col_data$segments==seg)
          
          # avoid error for segments w/o colors b/c they aren't in this region 
          if(dim(segcolrow)[[1]]>0) {
            plotcol = as.character(segcolrow[[2]])
            yval = as.numeric(segcolrow[[3]])
            # find start and end points
            xpos1 = min(as.numeric(row[9]),as.numeric(row[10]))
            xpos2 = max(as.numeric(row[9]),as.numeric(row[10]))
            # actually plot
            rect(xleft=xpos1, ybottom=yval-.49, xright=xpos2, ytop=yval+.49, density = -1, angle = 45, col = plotcol,border=plotcol) 
          }} # end hit loop 
        })
  
      output$save <- downloadHandler(
        file = "save.pdf" , # variable with filename
        content = function(file) {
          pdf(file = file)
          p2()
          dev.off()
        })
    })
  
  
  # display the plot data in case it needs to be examined 
  output$table2 <- renderTable({
    temp <- plot_data()
    if(!isTruthy(input$xmin)){
      xmin=min(range(temp$s.start),range(temp$s.end))
    } else {
      xmin=as.numeric(input$xmin)
    }
    if(!isTruthy(input$xmax)){
      xmax=max(range(temp$s.start),range(temp$s.end))
    } else {
      xmax=as.numeric(input$xmax)
    }
    temp2 <- subset(temp,temp$s.start >= xmin)
    temp3 <- subset(temp2,temp2$s.end <= xmax)
    df <-temp3[order(temp3$s.start),]
    return(df)
  })
} 

shinyApp(ui = ui, server = server)
