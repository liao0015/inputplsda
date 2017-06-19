

maxquant_plsda <- function(proteinInput, experimentInput){


#source("Modulized_functions.r")
suppressMessages(install.packages.auto(mixOmics))

# the first step is to readin the files and organize to generate a tidy formate

# since there is output of the tidied data, this line is only needed to run for the first time, start the analysis from the next step
my_tidy_data<-tidy_protein_groups(file_progeinGroups=proteinInput, 
                                file_experimentDesgin= experimentInput,
                                rows_marked_with = c("Only.identified.by.site", "Reverse","Potential.contaminant"),
                                cols_starts_with = "LFQ.intensity." # the columns to be selected starting with, use the maximum lengh
                                )


# next time you can start from here to save time
my_tidy_data<- readin_tidied_proteingroups(file_progeinGroups="Out_ProteinGroups_filtered.txt",  file_experimentDesgin= "Out_ProteinGroups_grouping.txt")


X <- as.matrix(my_tidy_data$data_matrix)
YL <- my_tidy_data$groups$Group


# # filter and imputate
X<-tidy_IntensityMarix_process(X, Q = 0.75, Imputation = TRUE)


# this part could be done by any method
my.pca <- mixOmics::pca(X, ncomp = 10, center = TRUE, scale = TRUE)

PCA_analysis<-PCA_post_analysis(my.pca, Y = YL)

#this is the plsda part, start a model, then test, how many component is good

my.plsda <- mixOmics::plsda(X, YL, ncomp = 3)

PLSDA_analysis <- PLSDA_post_analysis(my.plsda, X)



# figure test
pdf("plot.pdf")
  # scree plot to eveluate 
  PCA_analysis$plot_Scree
  # PCA plot 
  PCA_analysis$plot_component
  
  # overal component from PLSDA
  PLSDA_analysis$plot_PLSDA_component
  
  # Overall heatmap for all the whole matrix
  PLSDA_analysis$plot_PLSDA_heatmap
  
  # auc curve from the PLSDA 
  PLSDA_analysis$plot_auc
  
  # VIP distribution 
  PLSDA_analysis$plot_vip_distribution
dev.off()
  # heatmap of the selected features
png("plot.png")
  PLSDA_analysis$plot_vip_filtered_variables
print(PLSDA_analysis$plot_vip_filtered_variables)
dev.off()

# finanl task is to return the data.matrix selected
# this table is also output as "Out_VIP_table.txt" in background
taget_matrix <- PLSDA_analysis$vip.filtered.variables

}






#________________________________________________________________________________________

#     install.packages.auto
#________________________________________________________________________________________

# ___Description___: 
# 1: will check packages installed or not, if not will install from CRAN or bioconductor
# 2: the default setting is not update all dependent packages, which can be set ask = TRUE, or just coment this line
# 3: after install, load the package

# ___Arguments___:
# pacakge names, without quote
# 

#____Usage____;
# example: 
# install.packages.auto(qvalue) # from bioconductor
# install.packages.auto(rNMF) # from CRAN


install.packages.auto <- function(x) { 
  
  local({
    r <- getOption("repos")
    r["CRAN"] <- "http://cran.stat.sfu.ca/"
    options(repos = r)
  })
  
  
  x <- as.character(substitute(x)) 
  
  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    update.packages(ask= FALSE) #update installed packages.
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
  }
  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}


#________________________________________________________________________________________
#     color brewer section
#________________________________________________________________________________________


#brewer.pal
#mycol<- colorRampPalette(c("blue", "white", "red"))(100)


#________________________________________________________________________________________
#     PG_filterout_rows
#________________________________________________________________________________________

# ___Application___: 
# 1: parse the dataframe of proteingroups, filtering out the contaminant, reversed, and ided only by site


# ___input___:
# proteingroups is the dataframe read-in directly from proteingroups.txt (as default from Maxquant)
# do not open and edit it in excel
# rows_marked_with: rows marked in these columns as "+" will be deleted


# ___output___:
# a filtered data matrix


PG_filterout_rows<-function(proteinGroups, rows_marked_with = c("Only.identified.by.site", "Reverse","Potential.contaminant")){
  suppressMessages(install.packages.auto(dplyr))
  suppressMessages(install.packages.auto(lazyeval)) # this is very important for passing column names as parameters to the function
  
  proteinGroups_filtered<-proteinGroups
  
  for(filter_name in rows_marked_with)  {
    print(paste("filtering by",filter_name))
    filter_criteria <- lazyeval::interp(quote(x != "+"), x = as.name(filter_name))
    proteinGroups_filtered<-filter(proteinGroups_filtered, filter_criteria) 
    
  }
 
  # notice that the rownames are silently dropped druing filter, even there are row names, 
  # set the column names here, before the 
  rownames(proteinGroups_filtered)<-proteinGroups_filtered$id
  return(proteinGroups_filtered)  
}


#________________________________________________________________________________________
#     PG_filtering
#________________________________________________________________________________________

# ___Application___: 
# 1: parse the dataframe of proteingroups, filtering out the contaminant, reversed, and ided only by site
# 2: only keep the user defined expression columns, and return data matrix

#____input____
# proteingroups is the dataframe read-in directly from proteingroups.txt (as default from Maxquant)
# do not open and edit it in excel
# rows_marked_with: rows marked in these columns as "+" will be deleted
# cols_starts_with: column names start with, used the maxiumn length of the common characters 

# ___output___:
# a data matrix 

PG_filtering<-function(proteinGroups,rows_marked_with = c("Only.identified.by.site", "Reverse","Potential.contaminant"), cols_starts_with = "LFQ.intensity."){
  
  # removing rows marked with 
  proteinGroups_filtered<-PG_filterout_rows(proteinGroups, rows_marked_with)
  # keeping columns starts with
  proteinGroups_filtered <- select(proteinGroups_filtered, starts_with(cols_starts_with))
  #write.table(proteinGroups_filtered,"Out_ProteinGroups_filtered.txt",sep="\t",row.names = TRUE,col.names = NA)  
  return(proteinGroups_filtered)
}




#________________________________________________________________________________________
#     tidy_protein_groups
#________________________________________________________________________________________

# ___Description___: 
# 1: read in proteingroups
# 2: read in experimental desgin
# 3: filtering proteingroups by rows 
# 4: filtering proteingroups by columns
# 5: removing rows with all 0 
# 5: filtering experimental desgin, only keep samples with grouping information
# 6: filtering proteingroups, only keeping samples with grouping information
# 7: tranpose the filtered proteingroups: rows as samples/observations, columns as variables
# 8: reordering the grouping information, to math the ordering of the samples in proteingroups
# 9: double check the of the sample names are matching
# 10: return both the data matrix and experimental design/grouping 

# ___Arguments___:
# file_progeinGroups: full path for progeinGroups 
# row.names: The column name to be the rownames when readin the table, this column cannot be redundant. 
# proteingroups is the dataframe read-in directly from proteingroups.txt (as default from Maxquant)
# do not open and edit it in excel
# rows_marked_with: rows marked in these columns as "+" will be deleted
# cols_starts_with: column names start with, used the maxiumn length of the common characters 


#____Usage____;

# ___Values___:
# a list of three values
# dat_matrix: a tidy data matrix
# groups: a tidy and matched sample grouping file
# a quaility checking data frame to see if the sample order in the table is the same between the datamatrix and the grouping file


tidy_protein_groups <- function(file_progeinGroups="proteinGroups.txt", 
                                file_experimentDesgin= "experimentalDesign.txt",
                                rows_marked_with = c("Only.identified.by.site", "Reverse","Potential.contaminant"),
                                cols_starts_with = "LFQ.intensity." # the columns to be selected starting with, use the maximum lengh
                                )
  {
    # read in experimental desgin 
    #The place I made changes######
    grouping<-read.delim(file_experimentDesgin,header=TRUE, na.strings = "NA")
    #grouping<-ydata
    # clean the experiemntal desgin, in case there are some rawfiles not used in the grouping
    grouping <- replace(grouping, grouping == "", NA)
    temp<-apply(grouping,1,function(x) length(which(is.na(x))))
    grouping<-grouping[which(temp==0),]
    print("Experimental design read in")
    
    # read in protgroups
    print("Reading in ProteinGroups ....")
    #The place I made changes######
    proteinGroups<-read.delim(file_progeinGroups,header=TRUE, blank.lines.skip = TRUE)
    #proteinGroups<-xdata
    print("ProteinGroups read in")
    # filtering proteinGroups, only keeping expression data
    
    proteinGroups_filtered<-PG_filtering(proteinGroups,
                                         rows_marked_with = rows_marked_with, 
                                         cols_starts_with = cols_starts_with)

    #remove rows with all 0 
    
    proteinGroups_filtered<-proteinGroups_filtered[-which(apply(proteinGroups_filtered,1,function(x)all(x == 0))),]
    ##t<-proteinGroups_filtered[apply(proteinGroups_filtered, 1, function(x) !all(x==0)),]
    
    print("ProteinGroups filtered")
    
    #rm(proteinGroups) # release memory  
    proteinGroups_filtered<-t(proteinGroups_filtered)
    rownames(proteinGroups_filtered)<-gsub(cols_starts_with, "", rownames(proteinGroups_filtered))
    
    # only select columns in Experiment desgin:
    proteinGroups_filtered<- proteinGroups_filtered[which(rownames(proteinGroups_filtered) %in% grouping[,1]),]
    grouping<-grouping[match(rownames(proteinGroups_filtered),as.character(grouping[,1])),] # just in case the order is not the same
    alignment_check<-cbind(as.character(grouping[,1]),rownames(proteinGroups_filtered))
    
    write.table(proteinGroups_filtered,"Out_ProteinGroups_filtered.txt",sep="\t",row.names = TRUE,col.names = NA)  
    write.table(grouping,"Out_ProteinGroups_grouping.txt",sep="\t",row.names = TRUE,col.names = NA) 
    
    return(list(data_matrix = proteinGroups_filtered,
                groups = grouping,
                alignment_check = alignment_check
    ))
}


#________________________________________________________________________________________
#     readin_tidied_proteingroups
#________________________________________________________________________________________

# ___Description___: 
# 1: read in tidied proteingroups
# 2: read in tidied experimental desgin


# ___Arguments___:
# file_progeinGroups: full path for tidied progeinGroups ouput 
# file_experimentDesgin: full path for tidied experiment ouput 



#____Usage____;

# ___Values___:
# a list of three values
# dat_matrix: a tidy data matrix
# groups: a tidy and matched sample grouping file
# a quaility checking data frame to see if the sample order in the table is the same between the datamatrix and the grouping file



readin_tidied_proteingroups <- function(file_progeinGroups="Out_ProteinGroups_filtered.txt", 
                                file_experimentDesgin= "Out_ProteinGroups_grouping.txt")
{
  # read in experimental desgin 
  grouping<-read.delim(file_experimentDesgin,header=TRUE, row.names = 1)
  
  # read in protgroups
  proteinGroups<-read.delim(file_progeinGroups,header=TRUE, row.names = 1,  check.names = FALSE) 
  # make sure that check.names = FALSE, otherwise the there will be an X added to the colnames

  # filtering proteinGroups, only keeping expression data
  alignment_check<-cbind(as.character(grouping[,1]),rownames(proteinGroups))
  
  return(list(data_matrix = proteinGroups,
              groups = grouping,
              alignment_check = alignment_check
  ))
}


#________________________________________________________________________________________
#     IntensityMarix_process
#________________________________________________________________________________________

# ___Description___: 
# 1: very basic data process: log10 transformation, normalization on column or row, missing value imputation, etc
# 2: This function uses homemade subfunction "missingvalue_filtering", be sure to load it first

# ___Argument___:
# IntensityMatrix: is the data input, which is intensity matrix selected columns from the proteinGroup file, it could be the LFQ inetnsity, intensity, or intensity for specific labeling state
# Normalize_columns: if set TRUE, normlize column (experiment) means to the same level (to zero) first. Some method, like PCA, have scale function built in, therefore no normalization need before
# Normalize_rows: if set TRUE, normalize row (proteingroyps) means to the same level (tosezro) then. Some method, like PCA, have scale function built in, therefore no normalization need before
# 
# threshold: refer to the function of  "missingvalue_filtering"
#             briefly, 1 as Q100, 0.5 as Q50, 
# Imputation: if set TRUE, do missing value imputation. Only works when there are missing values todo:choose differnt method,

#___Usage___:
#IntensityMarix_process(IntensityMarix, threshold = 1, Imputation = TRUE, Normalize_columns = TRUE, Normalize_rows = TRUE)

# ___Values___:
# 
# 




IntensityMarix_process<-function(IntensityMarix, threshold = 1, Imputation = TRUE, Normalize_columns = TRUE, Normalize_rows = TRUE){
  
  IntensityMarix[IntensityMarix==0]<-NaN # replace the 0 with NaN
  IntensityMarix_log10<-log10(IntensityMarix) # take log10

  # missing value filtering
  tempt_filter_result <- missingvalue_filtering(IntensityMarix_log10, threshold)
  IntensityMarix_log10_NAfiltered <- tempt_filter_result$data_qualified # home made function [missingvalue_filtering]
  IntensityMarix_log10_NAfiltered_filtering_summary <-tempt_filter_result$filtering_summary

  p1<-matrix_ggboxplot(IntensityMarix_log10_NAfiltered, maintitle = "Distribution of NA-filtered")$violinplot
  p2<-matrix_ggboxplot(IntensityMarix_log10_NAfiltered, maintitle = "Distribution of NA-filtered")$boxplot
    
  if(Imputation == "TRUE"){  # missing value imputation
    IntensityMarix<-rrcovNA::impSeqRob(IntensityMarix_log10_NAfiltered)$x
  }
  
  if (Normalize_columns == "TRUE"){ # do column scaling, keep in mind that the scale function in R is scaling by column
    IntensityMarix<-scale(IntensityMarix)
  }
  
  if (Normalize_rows == "TRUE"){ # scaling of each protein
    IntensityMarix<-t(scale(t(IntensityMarix)))
  }
  
  p3<-matrix_ggboxplot(IntensityMarix, maintitle = "Distribution of NA-filtered & Processed")$violinplot
  p4<-matrix_ggboxplot(IntensityMarix, maintitle = "Distribution of NA-filtered & Processed")$boxplot
  
    
  write.table(IntensityMarix,paste("Out_ProteinGroups_NAfiltered_Scaled",".txt",sep=""),sep="\t",row.names = TRUE,col.names = NA)  
  return(list(IntensityMarix_processed = IntensityMarix,
              IntensityMarix_filtering_summary = IntensityMarix_log10_NAfiltered_filtering_summary,
              violinplot.before = p1,
              boxplot.before = p2,
              violinplot.after = p3,
              boxplot.after = p4
              
  ))
}

tidy_IntensityMarix_process<-function(IntensityMarix, Q = 0.75, Imputation = TRUE){
  suppressMessages(install.packages.auto(rrcovNA))
  
  IntensityMarix<-IntensityMarix[,colSums(IntensityMarix == 0) < (nrow(IntensityMarix)*(1-Q))]
  IntensityMarix<-log10(IntensityMarix)
  IntensityMarix[IntensityMarix==-Inf]<-NA
  if(Imputation=="TRUE"){
    IntensityMarix<-t(rrcovNA::impSeqRob(t(IntensityMarix))$x) # now is a inputed(log10(LFQintensity))
  }
  IntensityMarix
}


#________________________________________________________________________________________

#     missingvalue_filtering
#________________________________________________________________________________________

# ___Description___: 
# filter a matrix out rows/clumn with more than NA/infinte preset(inf values could be generated from log transformation or dividing conversion), 
# threshold is the number of the valid values, rows with more valid values than threshold will be kept as qualified
# only do row-wise filtering, transpose first if do column-wise filtering

#__Usage__:
# missingvalue_filtering(data, threshold=3)

# ___Arguments___:
# data: data matrix with missing values, NA/inf
# threshold: how many (percentage) non-missingvalues are required to be in the matrix
#             can be two types, one is the nnumber of the missing values,the other one is the so called Q value, which is the percentage(0 <= Q <= 1) to the number of columns
#             Q value == 1, require no missing values, 
#             Q value == 0, no filtering,
#             Q value will be converted to the number of missing value (1 < number < ncol(data))
#             will report an error if not setup in this range
#             in this function, celing is used to convert the percentage, if not expected, try floor/round etc

# ___Values___:
# qualified data matrix, not.qulified data.matrix, number of rows of quailfied/ not qulified. 
# 



missingvalue_filtering<-function(data.matrix, threshold = 1){ 
  data.matrix[is.infinite(as.matrix(data.matrix))]<-NA # the is.infinite function does not work on data.frame, 
  # in case there are infinte values there
  
  if(threshold < 0|threshold > ncol(data.matrix) ){
    print ("Threshold cannont be smaller than 0 or bigger than the number of columns, please check you threshold setting and rerun!!!")
  }else{

    if(threshold<=1){ # concet the q value to the real missing value number
      threshold <- ceiling(ncol(data.matrix)*threshold)
    }
    
    data_qualified<-data.matrix[which(apply(data.matrix,1,function(x)(sum(is.na(x)))<=(ncol(data.matrix)-threshold))),]
    data_not.qualified<-data.matrix[which(apply(data.matrix,1,function(x)(sum(is.na(x)))>(ncol(data.matrix)-threshold))),]
    return(list(data_qualified=data_qualified, 
                filtering_summary = list(data_not.qualified=data_not.qualified,
                number.qualified=nrow(data_qualified), 
                number.not.qualified=nrow(data_not.qualified))
    ))
    
  }
}





#________________________________________________________________________________________
#     PCA_plot
#________________________________________________________________________________________

# ___Description___: 
# 1: classical 2d plot of


# ___Arguments___:
# data matrix of tidy format
# grouping information is used just for plotting, if not given, there will be no grouping plot on the figure
# grouping information is required to be in a data.frame, with one column named as sample.name, one colum named as Groups

#____Usage____;

# PCA_plot(data_matrix, grouping)$pca.plot

# ___Values___:
# pca.plot as a ggplot2 object, to plot out the PCA plot
# 


PCA_plot<-function(data_matrix, grouping){
  suppressMessages(install.packages.auto(ggplot2))
  suppressMessages(install.packages.auto(ggfortify))# for autoplot
  data_matrix_t<-t(data_matrix)
  data_matrix_t_merge <- merge(grouping, data_matrix_t, by.y=0, by.x = "Sample.Name")  
  row.names(data_matrix_t_merge)<-data_matrix_t_merge$Sample.Name
  p1<-autoplot(prcomp(data_matrix_t), data = data_matrix_t_merge, colour = 'Groups',label = TRUE )
  p1<-p1+labs(title = "PCA Clustering")
  p1<-p1+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
  return(list(pca.plot = p1))
}






#________________________________________________________________________________________
#     PCA_plot_with_ellipse_kmeans
#________________________________________________________________________________________

# ___Description___: 
# 1: do pca on columns, 
# 2: plot a pca plot, with eliipse on the groups
# 3: note that is a indirect strategy, using result of prcomp to plot points, using kmeans grouping result to plot groups
# # level of 0.95 as confidence to draw the ellipse


# ___Argument___:
# a data matrix, with  samples as columns, and features as rows, 

#___Usage___:
# p <- PCA_plot_with_ellipse_kmeans(t(iris[,1:4]),iris$Species)
# p


# ___Values___:
#  a ggplot2 plot, with points and ellipse
# 


PCA_plot_with_ellipse_kmeans <- function(data_matrix, grouping){
  suppressMessages(install.packages.auto(ggplot2))
  pca    <- prcomp(t(data_matrix), retx=T, scale.=T) # do pca
  scores <- pca$x                       # scores for first three PC's
  #scores <- pca$x[,1:3]   
  # k-means clustering [assume 3 clusters]
  number_of_groups <- length(levels(as.factor(grouping))) 
  km  <- kmeans(scores, centers=number_of_groups, nstart=5)
  
  #ggdata <- data.frame(scores, Cluster=km$cluster, Groups=grouping)
  ggdata <- data.frame(scores, Cluster=km$cluster)
  
  
  p<- ggplot(ggdata) +
    geom_point(aes(x=PC1, y=PC2, colour=factor(Cluster)), size=3, shape=20) + 
    geom_text(aes(x=PC1, y=PC2, color=factor(Cluster),label=colnames(data_matrix)))+
    stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)), geom="polygon", level=0.95, alpha=0.2) +
    guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster"))+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5))
  
  
  return(p) 
  
}




#________________________________________________________________________________________

#     PCA_plot_3d_interactive
#________________________________________________________________________________________

# ___Description___: 
# interactive 3d plot of the pca result
# using 3d function of plotly 


# ___Arguments___:
# 
# grouping information is used just for plotting, if not given, there will be no grouping plot on the figure
# grouping information is required to be in a data.frame, with one column named as sample.name, one colum named as Groups

#____Usage____;

# PCA_plot_3d_interactive(data_matrix, grouping)$pca.plot

# ___Values___:
# a general objective of plotly
# 

PCA_plot_3d_interactive<-function(data_matrix, grouping){
  
  data_matrix_t<-t(data_matrix)
  
  
  loading <- prcomp(data_matrix_t)$x
  
  
  loading_merge <- as.data.frame(merge(grouping, loading, by.y=0, by.x = "Sample.Name"))  
  row.names(loading_merge)<-loading_merge$Sample.Name
  
  
  p1 <- plot_ly(loading_merge, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Groups, colors = c('#BF382A', '#0C4B8E', "#1ABC9C")) %>%
    add_markers() %>%
    #add_text(loading_merge, x = ~PC1, y = ~PC2, z = ~PC3,text = ~Sample.Name) %>% 
    add_text(text = ~Sample.Name) %>% 
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))
  
  return(list(pca.plot.3d = p1))
}





#________________________________________________________________________________________

#     PCA_Screeplot
#     PCA_Screeplot_2
#
#________________________________________________________________________________________

# ___Description___: 
# plot PCA screeplot, which shows how much of variance of each Principal Component explains 
# a good pca needs the firt 2/3 components explains most of the varaiance, otherwise, the separtaion is not good
# the pca analysis was done by prcomp
# an alternative is to to use a function of recordPlot to record the last plot as an object, which can be replotted after, where it was invoked. A list of plots can be returned by this way
# a better alternative is to use ggplot2, where all plots are objects

# ___Arguments___:
# PCA_Screeplot: data matrix as input
# PCA_Screeplot_@: the output of prcomp as input. in case you have already finished the pca analysis

#____Usage____;
# PCA_Screeplot(data_matrix)
# PCA_Screeplot_2(prcomp.out)

# ___Values___:
# plot a scree plot
# 

PCA_Screeplot<-function(data_matrix){
  suppressMessages(install.packages.auto(ggplot2))
  pca.output <- prcomp(t(data_matrix), scale.=TRUE, center = TRUE) 
  sd <- pca.output$sdev
  scores <- pca.output$x
  var <- sd^2
  var.percent <- var/sum(var) * 100
  
  #barplot(var.percent, xlab="Principal Component", ylab="Percent of Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="gray", main="Percent of Variance")
  #abline(h=1/nrow(pca.output$rotation)*100, col="red")
  #p1 <- recordPlot()
  
  p1<-ggplot()+geom_bar(aes(x=c(1:length(var.percent)),y=var.percent), stat="identity")
  p1<-p1+geom_hline(yintercept = 1/nrow(pca.output$rotation)*100, colour = "red")
  p1<-p1+labs(x = "Princaple Component Number",y="Percent of Variance",title = "Screeplot of Variance")
  p1<-p1+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
  
  return(list(Scree.plot = p1))
}



PCA_Screeplot_2<-function(prcomp.out){
  suppressMessages(install.packages.auto(ggplot2))
  sd <- pca.output$sdev
  scores <- pca.output$x
  var <- sd^2
  var.percent <- var/sum(var) * 100
  
  #barplot(var.percent, xlab="Principal Component", ylab="Percent of Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="gray", main="Percent of Variance")
  #abline(h=1/nrow(pca.output$rotation)*100, col="red")
  #p1 <- recordPlot()
  
  p1<-ggplot()+geom_bar(aes(x=c(1:length(var.percent)),y=var.percent), stat="identity")
  p1<-p1+geom_hline(yintercept = 1/nrow(pca.output$rotation)*100, colour = "red")
  p1<-p1+labs(x = "Princaple Component Number",y="Percent of Variance",title = "Screeplot of Variance")
  p1<-p1+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
  
  
  return(list(Scree.plot = p1))
}


PCA_post_analysis <- function(result.pca, Y){

  if(missing (result.pca)){
      print("No Input!")
  }else{
    
    install.packages.auto(mixOmics)
    my.pca <-result.pca
    
    
    postscript("temp")
    dev.control('enable')  
    plot(my.pca)
    p1 <- recordPlot()
    dev.off()
    
    postscript("temp")
    dev.control('enable')    
    plotIndiv(my.pca, comp = c(1, 2), ind.names = TRUE,  group = Y, ellipse = TRUE,legend = TRUE, title = 'PCA plot, PCA comp 1 - 2')
    p2 <- recordPlot()
    dev.off()
    
    return(list(plot_Scree = p1,
                plot_component = p2
    ))    
  }
  
  

  
}

#________________________________________________________________________________________
#     
#________________________________________________________________________________________

# ___Description___: 
# 1: 
# 2: 

# ___Arguments___:
# 
# 

#____Usage____;

# ___Values___:
# 
# 


plsda_CrossValication <-function(result_plsda, validation = "Mfold", folds =5, nrepeat =5 ){
  suppressMessages(install.packages.auto(mixOmics))
  print ("Performing crossing validation, will be very slow depending on the settings")
  perf.my.plsda.start <- mixOmics::perf(my.plsda.start, validation = 'Mfold', folds = fold, progressBar = TRUE, auc = TRUE, nrepeat = 2)
  # here add as many plot and return in the list
  #plotIndiv(my.plsda.start, comp = c(1, 2), ind.names = TRUE,  group = Y, ellipse = TRUE,legend = TRUE, title = 'PLSDA plot, Comp 1 - 2')
  plot.distance<-plot(perf.my.plsda.start, overlay = 'measure', sd = TRUE)
  return(list(plot.distance = plot.distance
              
  ))
}

#________________________________________________________________________________________
#  PLSDA_post_analysis   
#________________________________________________________________________________________

# ___Description___: 
# 1: The main target is to 
# 2: after setting up a PLSDA model, 

# ___Arguments___:
# 
# 

#____Usage____;

# ___Values___:
# 
# 


# Postprocess of PSLDA
PLSDA_post_analysis<-function(result.plsda, X){
  # ploting and output
  my.plsda<-result.plsda
  Y<-my.plsda$Y

  postscript("temp") 
  dev.control('enable')  
  plotIndiv(my.plsda, comp = c(1, 2), ind.names = TRUE,  group = Y, ellipse = TRUE,legend = TRUE, title = 'PLSDA plot, Comp 1 - 2')
  p1 <- recordPlot()
  dev.off()

  my.side.color <- color.mixo(as.numeric(Y))
  
  postscript("temp") 
  dev.control('enable') 
  cim(my.plsda, row.sideColors = my.side.color, row.names = Y)
  p2 <- recordPlot()
  dev.off()

   
  #plot an auc
  postscript("temp") 
  dev.control('enable') 
  my.plsda.auroc = mixOmics::auroc(my.plsda, roc.comp = 1)
  p3 <- recordPlot()
  dev.off() 
  
  # extract all the VIPs
  my.vip<-vip(my.plsda)
  write.table(my.vip,"Out_VIP_table.txt",sep="\t",row.names = TRUE,col.names = NA)  
  
  # oveall distribution of all the vips
  my.vip.plot <- matrix_ggboxplot(my.vip, xlabel="Component", ylabel = "VIP Score", maintitle = "VIP Score Across Component")
  p4 <- my.vip.plot$violinplot

  
  # filter vips, keeping >1, and ouput all the corresponding features/variables
  my.vip.filtered <- my.vip[my.vip[,1]>1,]
  my.vip.filtered <- my.vip.filtered[order(my.vip.filtered[,1]),]
  my.vip.filtered.variables <- X[,row.names(my.vip.filtered)]
  write.table(my.vip.filtered.variables,"Out_ProteinGroups_filtered_VIP.txt",sep="\t",row.names = TRUE,col.names = NA)  
  
  # plot the heatmap of the orginal features
  my.vip.filtered.variables.scaled<-scale(my.vip.filtered.variables)
  postscript("temp") 
  dev.control('enable') 
  cim(my.vip.filtered.variables.scaled, row.sideColors = my.side.color, row.names = Y, col.names = FALSE,row.cex = 0.5, scale = TRUE, center =TRUE)
  p5 <- recordPlot()
  dev.off()  
  
  
  return(list( plot_PLSDA_component = p1,
               plot_PLSDA_heatmap = p2,
               plot_auc = p3,
               plot_vip_distribution =p4,
               plot_vip_filtered_variables = p5, 
               vip.filtered.variables = my.vip.filtered.variables
              
  ))
}


#________________________________________________________________________________________
#     matrix_ggboxplot
#________________________________________________________________________________________

# ___Description___: 
# ggplot2 is powerful at boxplot and violinplot, but needs some pre-process of the datamatrix, a bit tricky sometimes
# here, the steps are wrapped up to give out an easy way

# ___Arguments___:
# data_matrix: data matrix
#  xlabel, ylabel, maintitle 

#____Usage____;
# boxplot_ressult <- matrix_ggboxplot(data_matrix, xlabel="Sample", ylabel = "Value", maintitle = "Distribution")
# plot by: boxplot_ressult$boxplot, boxplot_ressult$violinplot

# ___Values___:
# a list of plot, the first object boxplot, second one is violinplot



matrix_ggboxplot<-function(data_matrix, xlabel="Samples", ylabel = "Value", maintitle = "Distribution"){
  suppressMessages(install.packages.auto(ggplot2))
  data_matrix_melt<-reshape2::melt(as.matrix(data_matrix))
  # in data_matrix_melt, Var1 is the orignal row.names, Var2 is the orignial column names, value is the orignial values

  p1<-ggplot(data_matrix_melt, aes(x = Var2, y = value, fill=Var2))+geom_boxplot() 
  p1<-p1+labs(x = xlabel,y=ylabel,title = maintitle, fill = "Samples")
  p1<-p1+theme_bw()+theme(plot.title = element_text(hjust = 0.5)) 
  
  p2<-ggplot(data_matrix_melt, aes(x = Var2, y = value, fill=Var2)) +geom_jitter(shape=21,alpha=0.3) +geom_violin()
  p2<-p2+labs(x = xlabel,y=ylabel,title = maintitle)
  p2<-p2+theme_bw()+theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Samples"))
  
  return(list(boxplot = p1, violinplot = p2))
}


#________________________________________________________________________________________

#     flattenCorrMatrix
#________________________________________________________________________________________

# ___Description___: 
# 1: flatten a matrix from a correlation analysis, into a data.frame with 4 columns
# 2: 

# ___Arguments___:
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

#____Usage____;
# flattenCorrMatrix(res2$r, res2$P)

# Column 1 : row names (variable 1 for the correlation test)
# Column 2 : column names (variable 2 for the correlation test)
# Column 3 : the correlation coefficients
# Column 4 : the p-values of the correlations

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}







#________________________________________________________________________________________

#     correlation_matrix_plot
#________________________________________________________________________________________

# ___Description___: 
# 1: do correlation analysis of columns,
# 2: do visualization of result correlation matrix, with p values
# 3: a simple wrap of two funcitons: Hmisc::rcorr and corrplot::corrplot

# ___Arguments___:
# data matrix:
# a data matrix to do the column correlatin,

#  method_type:  Character, the visualization method of correlation matrix to be used. Currently, 
# it supports seven methods, named "circle" (default), "square", "ellipse", "number", "pie", "shade" and "color". 
# See examples for details.
# The areas of circles or squares show the absolute value of corresponding correlation coefficients. Method "pie" and "shade" came from Michael Friendly's job (with some adjustment about the shade added on), and "ellipse" came from D.J. Murdoch and E.D. Chow's job, see in section References.

# order_type: 
# corresponds the order or corrplot::corrplot

# "original" for original order (default).
# "AOE" for the angular order of the eigenvectors.
# "FPC" for the first principal component order.
# "hclust" for the hierarchical clustering order.
# "alphabet" for alphabetical order.


#____Usage____;
# t <- correlation_matrix_plot(mtcars, order_type = "hclust", plot_type ="color")
# t$corrplot
# t$corrmatrix

# ___Values___:
#  a corrplot and a very detailed matrix plot, see PerformanceAnalytics::chart.Correlation for more details
# 


correlation_matrix_plot <- function(data_matrix, order_type = "hclust", plot_type ="circle" ){
  
  correlation_matrix <- rcorr(as.matrix(data_matrix))
  win.metafile()
  dev.control('enable')
  corrplot(correlation_matrix$r, type="lower", order=order_type,  p.mat = correlation_matrix$P, sig.level = 0.01, insig = "pch", method =plot_type)
  p1 <- recordPlot()
  dev.off()
  
  win.metafile()
  dev.control('enable')
  suppressWarnings(chart.Correlation(data_matrix, histogram=TRUE, pch= "+"))
  p2 <- recordPlot()
  dev.off()
  
  return(list(corrplot = p1, corrmatrix = p2))
}



#________________________________________________________________________________________

#     filter_PSM
#________________________________________________________________________________________

# ___Description___: 
# 1: this function is for MSfragger result parsing and filter, not working for others
# 2: 

# ___Arguments___:
# file input could be peptXML, or tsv
# others could be just as default

#____Usage____;
#filter_PSM(file = "MS_QC_60min.pepXML")
#t<- filter_PSM(file = "MS_QC_60min.tsv")


# ___Values___:
# will write a table out, 
# return a data.frame of filtered table



filter_PSM <- function(file = "file", pepFDR=0.01, score_for_filtering = "hyperscore", decoyprefix="REVERSED_" ){

  suppressMessages(install.packages.auto(pepXMLTab))
  suppressMessages(install.packages.auto(tools))
  
  filetype = file_ext(file) #library(tools)
  
  if (filetype == "pepXML"){ # if this is a a peppxml file, parse it first
    tsv<-pepXML2tab(file)
  } else if(filetype == "tsv"){# if tsv, read in
    tsv <- read.delim(file, header =FALSE) 
    header= c(
      "scanID", "precursor_neutral_mass","retention_time_sec","assumed_charge", "hit_rank","peptide",
      "peptide_prev_aa","peptide_next_aa","protein","num_matched_ions","tot_num_ions of matched theoretical fragment ions",
      "calc_neutral_pep_mass","massdiff","num_tol_term","num_missed_cleavages","modifications",
      "hyperscore","next_score","intercept_of_ep_model","slope_of_pe_mode")
    colnames(tsv) <- header  
  }else{
    print ("Wrong input file type!")
  }
  
  passed <- PSMfilter(tsv, pepFDR=0.01, scorecolumn='hyperscore', hitrank=1, minpeplen=6, decoyprefix='REVERSED_')
  
  filename_base <-file_path_sans_ext(file) #library(tools)
  write.table(passed,paste(filename_base,"_FDR_filtered.txt", sep=""),sep="\t",row.names = TRUE,col.names = NA )
  return(list(filtered = passed))
}



#
function_template<-function(){}
#________________________________________________________________________________________
#     
#________________________________________________________________________________________

# ___Description___: 
# 1: 
# 2: 

# ___Arguments___:
# 
# 

#____Usage____;

# ___Values___:
# 
# 


# todo, matrix filter, filtering by groups
#
#mysubX <- split(X, Y)
#CD <- as.matrix(mysubX$CD)
#Control <- as.matrix(mysubX$Control)
#UC <- as.matrix(mysubX$UC)

#B <- X[,colSums(is.na(X)) < 3]
#B <- X[,colSums(X ==0) < (nrow(X)-3)]
#as.character(levels(Y))[2]
#C<- A[rowSums(is.na(A)) < ncol(A)/2, ]

