library(Rcpp)
library(mixOmics)
library(plotly)
library(clusterGeneration)
library(probFDA) # Probabilistic Fisher Discriminant Analysis
library(kernlab)
library(penalizedLDA)
library(dimRed)
#library(Rdimtools)
library(doParallel)

source("RPluMA.R")
source("RIO.R")


input <- function(inputfile) {
        parameters <<- readParameters(inputfile)


cl<<-makeCluster(parameters["nclust", 2])
registerDoParallel(cl)

ncomp <<- as.integer(parameters["ncomp", 2]) #Number of components for plsda
separationSeq <<- seq(parameters["sepFrom", 2], parameters["sepTo", 2], by = parameters["sepInc", 2]) #Never more than the samples
#separationSeq = seq(-0.9,0.9 , by = 0.1) #Never more than the samples
#samplesSeq = seq(10, 100010, by = 5000)
samplesSeq <<- seq(parameters["sampFrom", 2], parameters["sampTo", 2], by = parameters["sampInc", 2])
#nRepetitions = 100
nRepetitions <<- parameters["nRepetitions", 2]
repetitionsSeq <<- seq(1, nRepetitions, by = 1)
#nSignal = 10
#nNoise = 200
nSignal <<- parameters["nSignals", 2]
nNoise <<- parameters["nNoise", 2]
}

run <- function() {

outputMatrixPCA <<- array(0, dim=c(length(separationSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
outputMatrixPLSDA <<- array(0, dim=c(length(separationSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
outputMatrixpFDA <<- array(0, dim=c(length(separationSeq),length(samplesSeq)))
outputMatrixpLDA <<- array(0, dim=c(length(separationSeq),length(samplesSeq)))
outputMatrixICA <<- array(0, dim=c(length(separationSeq),length(samplesSeq)))
#outputMatrixLDA=array(0, dim=c(length(separationSeq),length(samplesSeq)))
#outputMatrixODP=array(0, dim=c(length(separationSeq),length(samplesSeq)))
#outputMatrixRLDA=array(0, dim=c(length(separationSeq),length(samplesSeq)))
#outputMatrixSPCA=array(0, dim=c(length(separationSeq),length(samplesSeq)))


#We execute with different matrixes to minimize random effects
for (repetition in repetitionsSeq) {
  print(repetition)
  print(repetitionsSeq)
  print(as.numeric(repetition/nRepetitions)*100)
  
  outputMatrixPCAAux = array(0, dim=c(length(separationSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  outputMatrixPLSDAAux = array(0, dim=c(length(separationSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  outputMatrixpFDAAux=array(0,dim=c(length(separationSeq),length(samplesSeq)))
  outputMatrixpLDAAux=array(0,dim=c(length(separationSeq),length(samplesSeq)))
  outputMatrixICAAux=array(0,dim=c(length(separationSeq),length(samplesSeq)))
 # outputMatrixLDAAux=array(0,dim=c(length(separationSeq),length(samplesSeq)))
#outputMatrixODPAux=array(0,dim=c(length(separationSeq),length(samplesSeq)))
#  outputMatrixRLDAAux=array(0,dim=c(length(separationSeq),length(samplesSeq)))
 # outputMatrixSPCAAux=array(0,dim=c(length(separationSeq),length(samplesSeq)))
  
  
  
  start.time <- Sys.time()

  
  #We iterate over different number of features
  separatioCount = 0
  for (nSeparation in separationSeq) {
    separatioCount = separatioCount +1
    # dir.create(paste(nSamples,"x",nSeparation))
    # setwd(paste(nSamples,"x",nSeparation))
    print (nSeparation)
    
    samplesCount = 0
    #We iterate over different number of signals
    for (nSamples in samplesSeq) {
      samplesCount = samplesCount +1

      #We create the data
      c <- genRandomClust(numClust=2, sepVal=nSeparation, numNonNoisy=nSignal, numNoisy=nNoise, 
                          numOutlier=0, numReplicate=1,clustszind= 1,clustSizeEq=(nSamples+250))

      X <- as.data.frame(c$datList$test_1[1:(nSamples*2),1:(nSignal+nNoise)])
      Y <- c$memList$test_1[1:(nSamples*2)]
      
      
      
      plsda <- splsda(X, Y, ncomp = ncomp, scale = TRUE)

      pca <- pca(X, ncomp = ncomp)
      write.csv(X, "x.csv")
      write.csv(Y, "y.csv")
      pfda<-pfda(X,Y,model="all",crit="cv",cv.fold=10,kernel="rbf")
   
      plda<-PenalizedLDA(as.matrix(X),Y,type="standard",K=1,lambda=0.15)
      
    #  ica<-do.ica(X,ndim=1,type='logcosh',tpar=1,maxiter=100)
      
    #  odp<-do.odp(X,Y,ndim=1)
       
    #  rlda<-do.rlda(X,Y,ndim=1)
      
     # spca<-do.spca(X,ndim=1)
      
      
      #We obain the variables that carry the signal
      signalVariables = setdiff(1:ncol(X),c$noisyList$test_1)
      
      #We compute the performance of the models as the percentage of signal features in the important features
      importantVariablesPCA = selectVar(pca, comp = 1)$name[1:nSignal]
      importantVariablesPLSDA = selectVar(plsda, comp = 1)$name[1:nSignal]
      importantVariablespFDA= head(order(pfda$V,decreasing = FALSE),n=nSignal)
      importantVariablespLDA=head(order(plda$discrim,decreasing=FALSE),n=nSignal)
    #  importantVariablesICA=head(order(ica$trfinfo$mean,decreasing=FALSE),n=nSignal)
     # importantVariablesODP=head(order(odp$trfinfo$mean,decreasing=FALSE),n=nSignal)
   #   importantVariablesRLDA=head(order(rlda$trfinfo$mean,decreasing=FALSE),n=nSignal)
     # importantVariablesSPCA=head(order(spca$trfinfo$mean,decreasing=FALSE),n=nSignal) 
      
      for (n in signalVariables) {
        if (paste('x',n,sep = "") %in% importantVariablesPCA){
          outputMatrixPCAAux[separatioCount,samplesCount] = outputMatrixPCAAux[separatioCount,samplesCount] + 1/nSignal
        }
        if (paste('x',n,sep = "") %in% importantVariablesPLSDA){
          outputMatrixPLSDAAux[separatioCount,samplesCount] = outputMatrixPLSDAAux[separatioCount,samplesCount] + 1/nSignal
        }
        if (paste(n,sep = "") %in% importantVariablespFDA){
          outputMatrixpFDAAux[separatioCount,samplesCount] = outputMatrixpFDAAux[separatioCount,samplesCount] + 1/nSignal
        }
        if (paste(n,sep = "") %in% importantVariablespLDA){
          outputMatrixpLDAAux[separatioCount,samplesCount] = outputMatrixpLDAAux[separatioCount,samplesCount] + 1/nSignal
        }
      #  if (paste(n,sep = "") %in% importantVariablesICA){
       #   outputMatrixICAAux[separatioCount,samplesCount] = outputMatrixICAAux[separatioCount,samplesCount] + 1/nSignal
       # }
      #  if (paste(n,sep = "") %in% importantVariablesODP){
      #    outputMatrixODPAux[separatioCount,samplesCount] = outputMatrixODPAux[separatioCount,samplesCount] + 1/nSignal
       # }
       # if (paste(n,sep = "") %in% importantVariablesRLDA){
       #   outputMatrixRLDAAux[separatioCount,samplesCount] = outputMatrixRLDAAux[separatioCount,samplesCount] + 1/nSignal
       # }
       # if (paste(n,sep = "") %in% importantVariablesSPCA){
       #   outputMatrixSPCAAux[separatioCount,samplesCount] = outputMatrixSPCAAux[separatioCount,samplesCount] + 1/nSignal
       # }
    }
  }
  
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  outputMatrixPLSDA <<- outputMatrixPLSDA + outputMatrixPLSDAAux
  outputMatrixPCA <<- outputMatrixPCA + outputMatrixPCAAux
  outputMatrixpFDA <<- outputMatrixpFDA + outputMatrixpFDAAux
  outputMatrixpLDA <<- outputMatrixpLDA + outputMatrixpLDAAux
 # outputMatrixICA = outputMatrixICA + outputMatrixICAAux
 # outputMatrixODP = outputMatrixODP + outputMatrixODPAux
  #outputMatrixRLDA = outputMatrixRLDA + outputMatrixRLDAAux
 # outputMatrixSPCA = outputMatrixSPCA + outputMatrixSPCAAux
}

outputMatrixPLSDA <<- outputMatrixPLSDA/nRepetitions
outputMatrixPCA <<- outputMatrixPCA/nRepetitions
outputMatrixpFDA <<- outputMatrixpFDA/nRepetitions
outputMatrixpLDA <<- outputMatrixpLDA/nRepetitions

}

output <- function(outputfile) {
#outputMatrixICA = outputMatrixICA/nRepetitions
#outputMatrixODP = outputMatrixODP/nRepetitions
#outputMatrixRLDA = outputMatrixRLDA/nRepetitions
#outputMatrixSPCA = outputMatrixSPCA/nRepetitions
write.csv(outputMatrixPLSDA, paste(outputfile, "plsda", "csv", sep="."))
write.csv(outputMatrixPCA, paste(outputfile, "pca", "csv", sep="."))
write.csv(outputMatrixpFDA, paste(outputfile, "pfda", "csv", sep="."))
write.csv(outputMatrixpLDA, paste(outputfile, "plda", "csv", sep="."))


	pcaP <- plot_ly(colors=c("dark blue","orange") ,showlegend=TRUE, name = "PCA", y = separationSeq,x = samplesSeq, z = outputMatrixPCA)%>%  
  layout(title = "Performance of PCA",
         annotations = list(x = 2.17,y = 0.63,text = 'PCA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "Samples"), 
           yaxis = list(title = "# Separation"), 
           zaxis = list(title = "Performance"))) %>% add_surface()
#pcaP

plsdaP <- plot_ly(colors=c("light blue","yellow"),y = separationSeq, x = samplesSeq, z = outputMatrixPLSDA,showlegend=TRUE, name="PLSDA") %>% 
  layout(title = "Performance of PLSDA vs PCA ",
         annotations = list(x = 1.1,y = 0.40,text = 'PLSDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "Samples"), 
           yaxis = list(title = "# Separation"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#plsdaP

pfdaP <- plot_ly(colors=c("pink","brown"),y = separationSeq, x = samplesSeq, z = outputMatrixpFDA,showlegend=TRUE, name="pFDA") %>% 
  layout(title = "Performance of pFDA",
         annotations = list(x = 0.5,y = 0.15,text = 'pFDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "Samples"), 
           yaxis = list(title = "# Separation"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#pFDA

pldaP <- plot_ly(colors=c("light green","dark green"),y = separationSeq, x = samplesSeq, z = outputMatrixpLDA,showlegend=TRUE, name="pFDA") %>% 
  layout(title = "Performance of pLDA vs pFDA",
         annotations = list(x = 0.5,y = 0.05,text = 'pLDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "Samples"), 
           yaxis = list(title = "# Separation"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#pLDA


#icaP <- plot_ly(colors=c("light pink","dark green"),y = separationSeq, x = samplesSeq, z = outputMatrixICA,showlegend=TRUE, name="ICA") %>% 
 # layout(title = "Performance of ICA",
 #        annotations = list(x = 0.5,y = 0.55,text = 'ICA',showarrow = FALSE),
  #       scene = list(
  #         xaxis = list(title = "Samples"), 
  #         yaxis = list(title = "# Separation"), 
   #        zaxis = list(title = "Performance")
 #        )) %>% add_surface()
#ICA


#odpP <- plot_ly(colors=c("grey","dark blue"),y = separationSeq, x = samplesSeq, z = outputMatrixODP,showlegend=TRUE, name="ODP") %>% 
#  layout(title = "Performance of ODP vs ICA",
 #        annotations = list(x = 0.5,y = 0.25,text = 'ODP',showarrow = FALSE),
 #        scene = list(
 #          xaxis = list(title = "Samples"), 
 #          yaxis = list(title = "# Separation"), 
  #         zaxis = list(title = "Performance")
 #        )) %>% add_surface()
#ODP


#rldaP <- plot_ly(colors=c("light yellow","brown"),y = separationSeq, x = samplesSeq, z = outputMatrixRLDA,showlegend=TRUE, name="RLDA") %>% 
 # layout(title = "Performance of RLDA",
   #      annotations = list(x = 0.5,y = 0.05,text = 'RLDA',showarrow = FALSE),
   #      scene = list(
   #        xaxis = list(title = "Samples"), 
   #        yaxis = list(title = "# Separation"), 
  #         zaxis = list(title = "Performance")
  #       )) %>% add_surface()
#RLDA

#spcaP <- plot_ly(colors=c("pink","blue"),y = separationSeq, x = samplesSeq, z = outputMatrixSPCA,showlegend=TRUE, name="SPCA") %>% 
 # layout(title = "Performance of SPCA vs RLDA",
  #       annotations = list(x = 0.5,y = 0.45,text = 'SPCA',showarrow = FALSE),
 #        scene = list(
#           xaxis = list(title = "Samples"), 
#           yaxis = list(title = "# Separation"), 
#           zaxis = list(title = "Performance")
#         )) %>% add_surface()
#SPCA


p <- subplot(pcaP,plsdaP)
p
#htmlwidgets::saveWidget(as_widget(p), outputfile)

p <- subplot(pfdaP,pldaP)
p

#htmlwidgets::saveWidget(as_widget(p), "graph2.html")

#p <- subplot(odpP)
#p

#htmlwidgets::saveWidget(as_widget(p), "graph3.html")

#p <- subplot(rldaP,spcaP)
#p

#htmlwidgets::saveWidget(as_widget(p), "graph4.html")


#title = paste("ClusterLowSamplesSeparation", "nSignal",nSignal,"nNoise",nNoise, "nRep",nRepetitions,".html",sep="")
#Sys.setenv("plotly_username" = "druiz072")
#Sys.setenv("plotly_api_key" = "Gnv8DXsx4cKYFt8uVD90")
#htmlwidgets::saveWidget(p,title)

stopCluster(cl)


}
