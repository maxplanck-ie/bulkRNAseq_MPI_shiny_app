## app.R ##
Rlib="/data/manke/sikora/shiny_apps/Rlibs3.5.0_bioc3.7"
.libPaths(Rlib)
library(shiny,lib.loc=Rlib)
library(shinydashboard,lib.loc=Rlib)
library(rhandsontable,lib.loc=Rlib)

options(shiny.maxRequestSize = 20*1024^2)

ui <- function(request) {dashboardPage(
    dashboardHeader(title = "Dataset selection"),
    ## Sidebar content
    dashboardSidebar(
      selectInput(inputId="genome", label="Select organism", choices=c("PLEASE SELECT AN ORGANISM","Zebrafish [zv9]","Zebrafish [zv10]","Fruitfly [dm6]","Fruitfly [dm3]","Human [hg37]","Human [hg38]","Mouse [mm9]","Mouse [mm10]"), selected = NULL),
        fileInput(inputId="countfile",label="Upload feature counts table.",multiple=FALSE,accept=NULL,buttonLabel = "Browse...", placeholder = "No file selected"),
        textOutput("fileDescription")
        ),
        
    dashboardBody(
        h2("RNAseq counts visualization"),
        uiOutput("resultPanels")
               
    )

 )}

server <- function(input, output, session) {
    
    output$walkThrough<-renderUI(HTML("<ul><li>1.Provide the organism information and upload feature counts count table. Your data will appear in the Input Annotation tab.</li><li>2.If necessary, relabel sample IDs to desired plotting IDs in the corresponding interactive sampleInfo table column, in the Input Annotation tab. Provide group information in the group column. Click on Run analysis. </li><li>3.Provide a text with ensembl gene IDs to analyze in the Data Visualization tab.</li></ul>"))
    output$FAQ<-renderText("Merging data from multiple datasets or batch effect removal are currenlty not supported.\n For questions, bug reports or feature requests, contact sikora@ie-freiburg.mpg.de.\n For reporting issues or pull requests on GitHub, go to https://github.com/maxplanck-ie/bulkRNAseq_MPI_shiny_app .")
    
    output$fileDescription<-renderText("Please provide an integer count table with ensembl gene IDs as rownames and sample IDs as colnames.")
    

##########################read/load processed data from a database##################
    library("data.table",lib.loc=Rlib)
    library("limma",lib.loc=Rlib)
    library("edgeR",lib.loc=Rlib)
    library("car",lib.loc=Rlib)
    library(gplots,lib.loc=Rlib)
    library(RColorBrewer,lib.loc=Rlib)
    library(ggplot2,lib.loc=Rlib)
    library(reshape2,lib.loc=Rlib)
    library(dplyr,lib.loc=Rlib)

    output$sessionInfo <- renderPrint({capture.output(sessionInfo())})
    ######################################################################################
    is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    check_counts<-function(counttable){
      if (all(is.character(counttable$GeneID),sum(grepl("[A-Z].+",counttable$GeneID,ignore.case=TRUE))==nrow(counttable),sum(apply(counttable[,2:ncol(counttable)],2,function(X)is.wholenumber(X)&X>=0))==(ncol(counttable)-1)*nrow(counttable))){
        res<-TRUE
      }else{
        res<-FALSE
      }
      return(res)
    }
    
    check_genes<-function(genefile,counttable){
      if (all(ncol(genefile)==1,nrow(genefile)>2,sum(genefile$V1 %in% counttable$GeneID)==length(genefile$V1))){
        res<-TRUE
      }else{
        res<-FALSE
      }
      return(res)
    }
    
    check_sampleInfo<-function(sampleSheet){
      if (all(sum(sampleSheet$Group %in% "NA")<nrow(sampleSheet),min(as.numeric(table(sampleSheet$Group[!sampleSheet$Group %in% "NA"])))>=2,length(table(sampleSheet$Group[!sampleSheet$Group %in% "NA"]))>=2)){
        res<-TRUE
      }else{
        res<-FALSE
      }
      return(res)
    }
    
   #######################################################################################  
    values <- reactiveValues()
    values$gtab<-""
   #######################################################################################  
    observe({
      if(!input$genome %in% "PLEASE SELECT AN ORGANISM"){
        genome<-isolate(input$genome)
        gv<-c("Zebrafish [zv9]"="Zv9","Zebrafish [zv10]"="GRCz10","Fruitfly [dm6]"="dm6","Fruitfly [dm3]"="dm3","Human [hg37]"="GRCh37","Human [hg38]"="GRCh38","Mouse [mm9]"="GRCm37","Mouse [mm10]"="GRCm38")
        gfile<-system(sprintf("find /data/manke/sikora/shiny_apps/bulkRNAseq_MPI_shiny_app/tables -name \"%s*\"",gv[genome]),intern=TRUE)
        gtab<-fread(gfile,header=FALSE,sep="\t")
        colnames(gtab)<-c("GeneID","GeneSymbol")
        if(grepl("\\.[0-9]{1,2}",gtab$GeneID[1])){gtab$GeneID<-gsub("\\.[0-9]+","",gtab$GeneID)}
        values$gtab<-gtab
        
      }
    })
    
    observeEvent(input$countfile, {
        
        
       if (!is.null(input$countfile)){
             countFile<-isolate(input$countfile)
             datPath<-countFile$datapath  
                          }
        inDr1<-read.table(datPath,header=TRUE,sep="\t",as.is=TRUE,quote="",nrows=1)
        inDr2<-suppressWarnings(fread(datPath,header=TRUE,sep="\t",nrows=1))
        if(ncol(inDr1)==ncol(inDr2)){cnv<-colnames(inDr1)
                                     cnv[1]<-"GeneID"}else{cnv<-c("GeneID",colnames(inDr1))}
        inDat<-fread(input=datPath,header=FALSE,sep="\t",skip=1)
        colnames(inDat)<-cnv
        ###################
        if(!ncol(inDat)>2|!nrow(inDat)>2){showModal(modalDialog(title = "COUNTFILE FORMAT INCORRECT",
                                                 "Please provide an integer counts table with samples as columns and genes as rows!",
                                                 easyClose = TRUE))}
        req(ncol(inDat)>2&nrow(inDat)>2)
        counts_ok<-check_counts(inDat)
        if(!isTruthy(counts_ok)){showModal(modalDialog(title = "COUNTFILE FORMAT INCORRECT",
                                                      "Please provide an integer counts table with samples as columns and genes as rows!",
                                                      easyClose = TRUE))}
        req(isTruthy(counts_ok))
        ######################
        ####handle gencode####
        if(grepl("\\.[0-9]{1,2}",inDat$GeneID[1])){inDat$GeneID<-gsub("\\.[0-9]+","",inDat$GeneID)}
        

###############initiate reactive table to collect sample information ###############

        
        DF<-data.frame(SampleID=colnames(inDat)[2:ncol(inDat)],PlottingID=colnames(inDat)[2:ncol(inDat)],Batch=(rep("NA",(ncol(inDat)-1))),Group=(rep("NA",(ncol(inDat)-1))),stringsAsFactors = F)
        observe({
          if (!is.null(input$hot)) {
            DF = hot_to_r(input$hot)
          } else {
            if (is.null(values[["DF"]]))
              DF <- DF
            else
              DF <- values[["DF"]]
          }
          values[["DF"]] <- DF
        })

        output$hot <- renderRHandsontable({
          DF <- values[["DF"]]
          if (!is.null(DF))
            rhandsontable(DF, stretchH = "all")
       })
        
        
        ###other reactive values
        values$normCounts<-""
        values$genes<-sample(inDat$GeneID,10)

     
        ######################
            observeEvent(input$runanalysis, {
                sampleInfo<-isolate(values[["DF"]])
                sinfo_ok<-check_sampleInfo(sampleInfo)
                if(!isTruthy(sinfo_ok)){showModal(modalDialog(title = "SAMPLE INFORMATION INCOMPLETE",
                                                              "Please provide group information for at least 2 replicates in at least 2 groups.",
                                                              easyClose = TRUE))}
                req(sinfo_ok)
                sampleInfo<-sampleInfo[!sampleInfo$Group %in% "NA",]
                rownames(sampleInfo)<-sampleInfo$PlottingID
                values$sampleInfo<-sampleInfo

                countdata<-as.data.frame(inDat[,match(sampleInfo$SampleID,colnames(inDat)),with=FALSE],stringsAsFactors=FALSE)
                colnames(countdata)<-sampleInfo$PlottingID
                rownames(countdata)<-inDat$GeneID
                output$countDatHead<-renderTable(head(countdata),include.rownames=TRUE,caption="Relabeled unnormalized input data",caption.placement = getOption("xtable.caption.placement", "top"),width=500)

                if(sum(sampleInfo$Batch=="NA")>0){
                design<-model.matrix(~1+Group,data=sampleInfo)} else { 
                design<-model.matrix(~1+Batch+Group,data=sampleInfo)}
                rownames(design)<-sampleInfo$Plotting_ID
                dge <- DGEList(counts=countdata)
                dge <- calcNormFactors(dge)
                v <- voom(dge,design,plot=FALSE)
                expdat<-v$E
                if(sum(sampleInfo$Batch=="NA")==0){
                  expdat<-removeBatchEffect(v$E,batch=sampleInfo$Batch)
                }
                values$normCounts<-expdat
                req(values$normCounts)
                showModal(modalDialog(title = "DATA NORMALISATION COMPLETE",
                                      "Size factors were used to normalize the raw data, which was subsequently transformed to log2CPM.",
                                      easyClose = TRUE))
                
            })#end of observe input$runanalysis
        
        
        ##########################read user input file with gene ids#######################################
        observeEvent(input$file1,{
          inFile<-isolate(input$file1)
          inTab<-read.table(inFile$datapath, header = FALSE,sep="\t",quote="",as.is=TRUE)
          if(grepl("\\.[0-9]{1,2}",inTab[1,1])){inTab[,1]<-gsub("\\.[0-9]+","",inTab[,1])}
          ###################
          genes_ok<-check_genes(inTab,inDat)
          if(!isTruthy(genes_ok)){showModal(modalDialog(title = "GENE ID FILE FORMAT INCORRECT",
                                                        "Please provide a single column headerless file of at least 3 gene IDs matching your countdata!",
                                                        easyClose = TRUE))}
          req(isTruthy(genes_ok))
          values$genes<-inTab$V1
          
        })#end of observe input$file1
        
        
        observe({
        sampleInfo<-values$sampleInfo  
        normCounts<-values$normCounts
        genes<-values$genes
        req(normCounts,genes)
        selNormCounts<-normCounts[match(genes,rownames(normCounts)),]
        selNormCounts<-selNormCounts[complete.cases(selNormCounts),]
        output$downloadNormCounts <- downloadHandler(filename = "NormalizedLog2CPM.xls",content=function(file) {write.table(selNormCounts,file,row.names = TRUE,sep="\t",quote=FALSE,dec = ",")})
        ##annotate with external gene symbol  ----->>> only if organism was selected, otherwise, copy over gene IDs!
        plotdata<-as.data.frame(selNormCounts,stringsAsFactors=FALSE)
        plotdata$GeneID<-rownames(selNormCounts)
        
        if(isTruthy(values$gtab)){
          gtab<-isolate(values$gtab)
          plotdata$GeneSymbol<-gtab$GeneSymbol[match(plotdata$GeneID,gtab$GeneID,nomatch=NA)]
          if(sum(is.na(plotdata$GeneSymbol))==nrow(plotdata)){showModal(modalDialog(title = "ORGANISM INFORMATION NOT MATCHING COUNTFILE",
                                                                                    "None of the Gene IDs in your countfile matched the selected organism annotation, GeneIDs will be used instead of Gene Symbols.",easyClose = TRUE))}
          plotdata$GeneSymbol[is.na(plotdata$GeneSymbol)|plotdata$GeneSymbol==""]<-plotdata$GeneID[is.na(plotdata$GeneSymbol)|plotdata$GeneSymbol==""]
        }else{plotdata$GeneSymbol<-plotdata$GeneID}
        
        if(input$XlabelChoice=="GeneSymbol"&!isTruthy(values$gtab)){showModal(modalDialog(title = "ORGANISM INFORMATION NOT AVAILABLE",
                                                                                         "In order to use Gene Symbol annotation, please select an organism from the list on the sidebar!",easyClose = TRUE))}
        
        ifelse(input$CBfriendly,colpalette<-brewer.pal(4,"Paired")[c(2,4,1,3)],colpalette<-brewer.pal(8,"Dark2"))
        
        output$heatmap<-renderPlot({
          heatmap.2(selNormCounts, scale="row", trace="none", Rowv=FALSE,dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),main=input$projectid, keysize=1.5,margins = c(10,18),labRow=plotdata[,colnames(plotdata) %in% input$XlabelChoice],density.info="none",ColSideColors=colpalette[factor(sampleInfo$Group[match(colnames(selNormCounts),sampleInfo$SampleID)])]) })
        
        x_choice <- reactive({switch(input$XlabelChoice,"GeneID"="GeneID","GeneSymbol"="GeneSymbol")})
        #x_choice<-"GeneID"
        plotdataL<-melt(plotdata,id.vars=c("GeneID","GeneSymbol"),value.name="Log2CPM",variable.name="SampleID")
        plotdataL$Log2CPM<-as.numeric(plotdataL$Log2CPM)
        plotdataL$Group<-sampleInfo$Group[match(plotdataL$SampleID,sampleInfo$PlottingID)]
        plotdataL$GeneID<-factor(plotdataL$GeneID,levels=plotdata$GeneID)
        plotdataL$GeneSymbol<-factor(plotdataL$GeneSymbol,levels=plotdata$GeneSymbol)
        ifelse(input$revlev,plotdataL$Group<-factor(plotdataL$Group,levels=rev(unique(plotdataL$Group))),plotdataL$Group<-factor(plotdataL$Group,levels=unique(plotdataL$Group)))
        output$jplot<-renderPlot({
        ggplot(data=plotdataL,aes(x=eval(as.name(x_choice())),y=Log2CPM,group=Group,colour=Group))+geom_point(size=2,alpha=0.6,position=position_jitterdodge(jitter.width=0.2,jitter.height=0.000001,seed=314))+scale_colour_manual(values=colpalette)+theme(text = element_text(size=16),axis.text = element_text(size=14),axis.text.x=element_text(angle=90,vjust=0),axis.title = element_text(size=14)) +xlab(x_choice())
                               })# end of render jplot
        
        g1<-unique(sampleInfo$Group)[1]
        g2<-unique(sampleInfo$Group)[2]
        
        #res<-{plotdataL %>% group_by(Group,GeneID) %>% summarise(value = list(Log2CPM)) %>% spread(GeneID, value) %>% group_by(GeneID) %>% 
         # mutate(p_value = t.test(unlist(g1), unlist(g2))$p.value, paste0("mean_",g1) = mean(unlist(g1)), paste0("mean_",g2) = mean(unlist(g2)))}
        
        restemp<-summarize(group_by(plotdataL,Group,GeneID),GroupMean=mean(Log2CPM))
        restempW<-dcast(restemp,GeneID~Group)
        restempW$GeneSymbol<-plotdata$GeneSymbol[match(restempW$GeneID,plotdata$GeneID)]
        restempW$pvalue<-unlist(lapply(split(plotdataL,factor(plotdataL$GeneID,levels=unique(plotdataL$GeneID))),FUN=function(X)t.test(Log2CPM~Group,data=X)$p.value))
        restempW$padj<-p.adjust(restempW$pvalue)
        restempW<-as.data.frame(restempW,stringsAsFactors=FALSE)
        
        output$ttest<-renderTable({restempW})
        output$downloadTtest <- downloadHandler(filename = "t.test_results.xls",content=function(file) {write.table(restempW,file,row.names = FALSE,sep="\t",quote=FALSE,dec = ",")})

        })  ##end of observe          
        

        },ignoreInit=TRUE)#end of observe input$submitinput
    
    output$get_vignette <- downloadHandler(
      filename = "bulkRNAseq_app_vignette.html",
      content = function(con) {
        file.copy(from="/data/manke/sikora/shiny_apps/bulkRNAseq_docs/bulkRNAseq_docs.html", to=con, overwrite =TRUE)
      }
    )
    
    output$downloadSessionInfo <- downloadHandler(
      filename = "sessionInfo.txt",
      content = function(con) {
        sink(con)
        print(sessionInfo())
        sink()
      }
    )
   
    
###################################################################################   
    

       
    output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="Walkthrough",
                                                      fluidPage(
                                                          box(title="Walkthrough",uiOutput("walkThrough")),
                                                          box(title="Miscellaneous information",textOutput("FAQ")),
                                                          downloadButton("get_vignette", label = "Download vignette html")
                                                               )
                                                          ),
                                                  tabPanel(title="Input Annotation",
                                                      fluidPage(
                                                          #fluidRow(
                                                              #tableOutput("datHead"),
                                                              #tableOutput("countDatHead")
                                                           #        ),
                                                          fluidRow(
                                                            box(title="Sample description table",rHandsontableOutput("hot",width=800),width=8)
                                                              
                                                                  ),
                                                          fluidRow(
                                                            actionButton(inputId="runanalysis",label="Run analysis",width="200px",style = "color: black;background-color:  	 	#6495ED")
                                                          ),
                                                          fluidRow(
                                                              box(renderText("Please complete missing sample information. Plotting ID will be used to replace Sample ID on plots, Group will be used to construct a design table. Any samples with missing annotation will be removed from analysis."))
                                                                  ) 
                                                              )
                                                          ),
                                                  tabPanel(title="Data Visualization",
                                                      fluidPage(
                                                        fluidRow(
                                                          box(title="Heatmap",plotOutput("heatmap"),height=600),
                                                          box(title="Jitter plot",plotOutput("jplot"),height=600)
                                                        ),
                                                        fluidRow(
                                                          box(title="Gene list",fileInput(inputId="file1", label="Upload gene list.", multiple = FALSE, accept = NULL, width = NULL,buttonLabel = "Browse...", placeholder = "No file selected")),
                                                          box(title = "Plot controls",
                                                              selectInput("XlabelChoice", "Gene label",choices=c("GeneID","GeneSymbol"),selected="GeneID"),
                                                              checkboxInput("CBfriendly", "Colour-blind friendly (up to 4 colours)"),
                                                              checkboxInput("revlev", "Reverse group order"))),
                                                        fluidRow(
                                                          box(title="Method Description",renderText("Genewise Log2 counts per million were mean-centered and scaled across the samples. Average linkage clustering on euclidean distances was performed.")),
                                                          box(title="Method Description",renderText("Genewise Log2 counts per million are plotted using small amounts of jitter as offset.")))
                                                               ),
                                                      downloadButton("downloadNormCounts", label="Download Normalized Counts",style = "color: black;background-color: #6495ED")
                                                          ),
                                                tabPanel(title="Statistical test",
                                                         fluidPage(
                                                           tableOutput("ttest"),
                                                           downloadButton(outputId="downloadTtest", label="Download test results")
                                                         )
                                                         ),
                                                  tabPanel(title="Session Info",
                                                      fluidPage(
                                                          verbatimTextOutput("sessionInfo"),
                                                          downloadButton(outputId="downloadSessionInfo", label="Download session info")
                                                               )
                                                          )
                                                
                                                )


            do.call(tabsetPanel, myTabs)})


}

shinyApp(ui, server)
