#delete memory
rm(list = ls(all = TRUE))
# library(MultiAssayExperiment)
library(SingleCellExperiment)
library(scater)
# library(scran)
library(plyr)
library(dplyr)
library(ggplot2)

all_time_gene_select<-0
special_function=NA
start_time=0
print(paste("start_time:",start_time,sep = ""))

#if you use raw function,this raw_gene_select_function=T;
# raw_gene_select_function='raw_gene_select_function'

raw_gene_select_function='new_gene_select_function'

if(raw_gene_select_function=='raw_gene_select_function'){
  raw_gene_select_function=='raw_gene_select_function'
}else{
  raw_gene_select_function=='new_gene_select_function'
}
# #设置以下值，预处理采用新的固定方法
# preprocessing_function='new'

flag<-"各基因选择方法研究"
acc_top1 <- function(dataname,datatype,pctkeep){                       # dataname的类型为Biase.Rdata
  # flag0 <- paste("Seurat  hvg ) ",sizePer,sep = "")
  
  datatype=datatype
  pctkeep=pctkeep
  cat(dataname)                                       # 在屏幕上输出dataname
  cat("\n")                                      
  dir <- "gene_select_data_confirm/"
  dataname1 <- substr(dataname, 1, nchar(dataname)-6)
  load(paste(dir,"/",dataname,sep = ""))

  # if(any(dataname=="Patel.Rdata",dataname=="Kolod.RData"))
  # {
  #   cat("to log data.....\n")
  #   data <- 2^data
  # 
  # }
  
  #没有基因名，细胞名，设置名称
  if(any(dataname=="Haber.RData",dataname=="Macosko.RData",dataname=="Test_38_correct.RData")) 
  {
    cat("to t data.....\n")
    data <- t(data)
    rownames(data) = paste("g", 1: dim(data)[1], sep = "")
    colnames(data)=paste("c", 1: dim(data)[2], sep = "")
  }
  if(any(dataname=="Cao.Rdata",dataname=="Han.Rdata",dataname=="Klein.Rdata",dataname=="Pbmc68k.Rdata",dataname=="Shekhar.Rdata")) 
  {
    # cat("to t data.....\n")
    # data <- t(data)
    rownames(data) = paste("g", 1: dim(data)[1], sep = "")
    colnames(data)=paste("c", 1: dim(data)[2], sep = "")
  }
  
  if(F){
    raw_gene_select=raw_gene_select_function
    file_name=paste(raw_gene_select,"_cluster.csv",sep="")
    file_name_time=paste(raw_gene_select,"_run_time.csv",sep="")
    file_name_gene_select_counts=paste(raw_gene_select,"_gene_select_counts.csv",sep="")
    # write.table(dataname1,"spearman_4data_log_restore_ne_SIMLR_SIMLR_no_normal_hc_label_cluster.csv",append=T,col.names = FALSE,sep=",",row.names = FALSE)
    write.table(dataname1,file_name,append=T,col.names = FALSE,sep=",",row.names = FALSE)
    write.table(dataname1,file_name_time,append=T,col.names = FALSE,sep=",",row.names = FALSE)
    write.table(dataname1,file_name_gene_select_counts,append=T,col.names = FALSE,sep=",",row.names = FALSE)
  }
  
  #创建对象
  counts<-as.matrix(data)
  sce <- SingleCellExperiment(
    assays = list(counts = counts), 
    rowData=data.frame(gene=rownames(counts)),
    colData=data.frame(cell=colnames(counts),label_true=meta$label)
  )
  
  write.table (flag, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  write.table (dataname1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  gene_length<-paste("原始基因长度:",dim(assay(sce,"counts"))[1],sep="")
  cell_length<-paste("原始细胞长度:",dim(assay(sce,"counts"))[2],sep="")
  write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
  write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  
  
  data_type<-paste("数据类型：",datatype,sep="")
  write.table (data_type, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  
  #原始值最大值最小值
  max_value<-paste("最大值:",max(assay(sce,"counts")),sep="")
  min_value<-paste("最小值:",min(assay(sce,"counts")),sep="")
  write.table (max_value, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
  write.table (min_value, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  
  
  if(raw_gene_select_function=='raw_gene_select_function'){
    if(T)
    {
      gene_select_dir<-"lry_gene_select/"
      gene_select_names<-list.files(gene_select_dir)
      for(i in 1:length(gene_select_names)){
        
        start_time <- Sys.time()
        print(paste("raw_start_time:",start_time,sep = ""))
        select_gene_function <- substr(gene_select_names[i], 1, nchar(gene_select_names[i])-2)
        cat (select_gene_function)
        cat("\n")
        # flag0 <- "a variety of gene select function"
        flag1<-paste("基因选择方法：",select_gene_function,sep="")
        write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        # write.table (flag1, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        
        #加入保存Excel
        raw_gene_select=raw_gene_select_function
        file_name_0=paste(raw_gene_select,"_",sep="")
        file_name_0_0=paste(file_name_0,select_gene_function,sep="")
        file_name=paste(file_name_0_0,"_cluster.csv",sep="")
        file_name_time=paste(file_name_0_0,"_run_time.csv",sep="")
        file_name_gene_select_counts=paste(file_name_0_0,"_gene_select_counts.csv",sep="")
        
        
        # write.table(dataname1,"spearman_4data_log_restore_ne_SIMLR_SIMLR_no_normal_hc_label_cluster.csv",append=T,col.names = FALSE,sep=",",row.names = FALSE)
        write.table(dataname1,file_name,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        write.table(dataname1,file_name_time,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        write.table(dataname1,file_name_gene_select_counts,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        
        
        gene_select_dir_name<-paste(gene_select_dir,gene_select_names[i],sep="")
        source(gene_select_dir_name)
        # browser()
        # start_time_gene_select <- Sys.time() 
        #run gene select
        
        run_select_gene_function<-gene_select(sce,pctkeep)
        
        special_function=select_gene_function
        
        if(special_function=='filterPCA')
        {
          logcounts(sce)<-log2(counts(sce)+1)
          # browser()
          dat <- logcounts(sce)
          TrueClu <- length(table(colData(sce)[['label_true']]))
          pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
          # browser()
          
          #notices:this pca dimension set the same as f1000 nPC30
          # pca <- pca$x[, seq_len(params$nPC), drop = FALSE]
          pca<-pca$x[,seq_len(30),drop=FALSE]
          hcl <- hclust(dist(pca), method = "ward.D2")
          cluster <- cutree(hcl, k = TrueClu)
          label.true <- sce$label_true
          label.cluster<-cluster
          
        }else if(special_function=='filterRAFSIL')
        {
          
          source("gene_select_other_package/RAFSIL/R/RAFSIL.R")
          source("gene_select_other_package/RAFSIL/R/rafsil_FE.R")
          source("gene_select_other_package/RAFSIL/R/rafsil_RF.R")
          source("gene_select_other_package/RAFSIL/R/RAFSIL_utils.R")
          # library(grid)
          # library(gridExtra)
          #install.packages("ClusterR",repos = "https://cran.rstudio.com")
          library(ClusterR)
          #install.packages("randomForest",repos = "https://cran.rstudio.com")
          library(randomForest)
          #pamk()方法
          library(fpc)
          #make.symmetric(X) : could not find function "make.symmetric"
          #but neet package rjags
          # install.packages("rjags",repos = "https://cran.rstudio.com")
          library(rjags)
          #need download:https://sourceforge.net/projects/mcmc-jags/,install,download d/soft/JAGS-4_3_0.exe
          # install.packages("jags",repos = "https://cran.rstudio.com")
          # install.packages("dclone",repos = "https://cran.rstudio.com")
          library(dclone)
          #install.packages("e1071",repos = "https://cran.rstudio.com")
          library(e1071)
          # dat=readRDS('gene_select_other_package/RAFSIL/inst/extdata/example_data.rds')
          # dim(dat$expr)
          # table(dat$labels)
          
          dat=assay(sce,"counts")
          TrueClu <- length(table(colData(sce)[['label_true']]))
          #Next, derive similarities between the cells using RAFSIL1, RAFSIL2, and Spearman correlation.
          
          #- run RAFSIL1 with 50 forests
          res.r1 = RAFSIL(t(dat),nrep = 50, method="RAFSIL1",NumC=TrueClu)
          
          # res.r2 = RAFSIL(t(dat),           method="RAFSIL2")
          # 
          # #- retriev the dissimilarities
          # dis.r1  = res.r1$D
          # dis.r2  = res.r2$D
          # dis.cor = sqrt((1 - cor(dat,method="spearman"))/2)
          # #next plot visual:see :http://www.kostkalab.net/pub_software/RAFSIL/RAFSIL.html
          # hcl <- hclust(dist(dis.cor), method = "ward.D2")
          # cluster <- cutree(hcl, k = TrueClu)
          # label.true <- sce$label_true
          # label.cluster<-cluster
          #李荣远修改只用随机森林聚类1
          label.true <- sce$label_true
          label.cluster<-res.r1$lab
          #this change hc
        }else if(special_function=='filterSC3')
        {
          source("gene_select_other_package/SC3/myMethods.R")
          source("gene_select_other_package/SC3/SC3_OK.R")
          cluster.result <- SC3_OK(inputTags = assay(sce,"counts"), datatype = datatype, SEED=123)
          # cluster.result <- SC3_OK_lrychange(inputTags = sce, datatype = datatype, SEED=123)
          label.true <- sce$label_true
          label.cluster<-cluster.result
          #this change hc
        }else
        {
          flag1<-paste("基因选择%：",pctkeep,sep="")
          write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
          
          flag_not_expressed<-"最终选择的基因和细胞数"
          write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
          gene_length<-paste("基因长度:",dim(assay(run_select_gene_function,"counts"))[1],sep="")
          cell_length<-paste("细胞长度:",dim(assay(run_select_gene_function,"counts"))[2],sep="")
          write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
          write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
          
          write.table(gene_length,file_name_gene_select_counts,append=T,col.names = FALSE,sep=",",row.names = FALSE)
          
          # end_time_gene_select <- Sys.time()
          # all_time_gene_select <-difftime(end_time_gene_select,start_time_gene_select,units="mins")
          
          # TrueClu<-length(levels(factor(phn$phenoid)))
          
          TrueClu <- length(table(colData(sce)[['label_true']]))
          hcl <- hclust(dist(t(assay(run_select_gene_function,"logcounts"))), method = "ward.D2")
          label.cluster <- cutree(hcl, k = TrueClu)
          
          
          # 根据聚类结果画聚类树及类别框
          # plot(result_hc)
          # label.true <- (as.list(meta["label"]))$label
          
          # label.true <- run_select_gene_function$phenoid
          label.true <- run_select_gene_function$label_true
        }
        
        
        # 聚类评估
        # 导入聚类评估函数
        source("evalcluster.R")
        NMI_ARI=evalcluster(label.true,label.cluster)
        
        end_time <- Sys.time()
        print(paste("raw_end_time:",end_time,sep = ""))
        all_time <-difftime(end_time,start_time,units="mins")
        
        # write.table (dataname1, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        # write.table (run_select_gene_function, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        
        write.table(NMI_ARI,file_name,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        write.table(all_time,file_name_time,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        # flag1<-paste("基因选择运行时间(mins)：",all_time_gene_select,sep="")
        # write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        
        write.table (NMI_ARI, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =T, col.names =F, quote =F,append = T,)
        flag1<-paste("总运行时间(mins)：",all_time,sep="")
        write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        
      }
      
    }
    
  }else{
    #delete cell
    
    start_time <- Sys.time()
    print(paste("new_start_time:",start_time,sep = ""))
    
    keep_cell <- colSums(counts(sce) > 0) > 10
    table(keep_cell)
    sce <- sce[,keep_cell]
    
    flag_not_expressed<-"保留列细胞表达值大于0的统计量大于10的细胞"
    write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    gene_length<-paste("基因长度:",dim(assay(sce,"counts"))[1],sep="")
    cell_length<-paste("细胞长度:",dim(assay(sce,"counts"))[2],sep="")
    write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    #delete gene
    # keep_features <- rowSums(counts(sce) > 0) > 0.06*dim(assay(sce,"counts"))[2]
    keep_features <- rowSums(counts(sce) > 0) > 0
    table(keep_features)
    sce <- sce[keep_features, ]
    dim(sce)
    
    # flag_not_expressed<-"统计基因值大于0的个数。该值大于0.06倍细胞，才保留该基因"
    flag_not_expressed<-"删除表达值全0的基因"
    write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    gene_length<-paste("基因长度:",dim(assay(sce,"counts"))[1],sep="")
    cell_length<-paste("细胞长度:",dim(assay(sce,"counts"))[2],sep="")
    write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
    #QC main delete cell
    per.cell <- perCellQCMetrics(sce)
    #new scater,see example
    # per.cell <- perCellQCMetrics(example_sce, 
    #                              subsets=list(Mito=grep("mt-", rownames(example_sce))))
    # colData(example_sce) <- cbind(colData(example_sce), per.cell)
    colData(sce) <- cbind(colData(sce), per.cell)
    # colData(sce) <- cbind(colData(sce), sce)
    
    #去除低于3倍绝对中位差的值。 mad(sce$total_counts)*3。认为
    # 为何设定3倍的MAD？ ：通常把偏离中位数三倍以上的数据作为异常值，和均值标准差方法比，中位数和MAD的计算不受极端异常值的影响，结果更加稳健。
    
    colData(sce)$libsize.drop <- isOutlier(sce$total, nmads = 3, type = "lower", log = TRUE)
    
    # sce$libsize.drop <- isOutlier(sce$total, nmads = 3, type = "lower", log = TRUE)
    
    table(colData(sce)$libsize.drop)
    
    colData(sce)$feature.drop <- isOutlier(sce$detected, nmads = 3, type = "lower", log = TRUE)
    table(colData(sce)$feature.drop)
    
    table(libsize = sce$libsize.drop, feature = sce$feature.drop)
    
    #最后只筛选去掉以上两者离群值的数据
    sce <- sce[, !(sce$libsize.drop | sce$feature.drop)]
    dim(sce)
    
    flag_not_expressed<-"最后只筛选去掉以上两者离群值的数据"
    write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    gene_length<-paste("基因长度:",dim(assay(sce,"counts"))[1],sep="")
    cell_length<-paste("细胞长度:",dim(assay(sce,"counts"))[2],sep="")
    write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
    #normalization
    counts <- assay(sce, "counts")
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)
    assayNames(sce)
    
    end_time_pre_data <- Sys.time()
    all_time_pre_data <-difftime(end_time_pre_data,start_time,units="mins")
    print(paste("pre_data_use_time:",all_time_pre_data,sep = ""))
    
    
    flag_not_expressed<-"固定规范化后值大小，来源scater的规范化"
    write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    max_value<-paste("最大值:",max(assay(sce,"counts")),sep="")
    min_value<-paste("最小值:",min(assay(sce,"counts")),sep="")
    write.table (max_value, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    write.table (min_value, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
    #gene select and cluster
    if(T)
    {
      gene_select_dir<-"lry_gene_select/"
      gene_select_names<-list.files(gene_select_dir)
      for(i in 1:length(gene_select_names)){
        
        start_time <- Sys.time()
        print(paste("new_skip_pre_data_start_time:",start_time,sep = ""))
        select_gene_function <- substr(gene_select_names[i], 1, nchar(gene_select_names[i])-2)
        cat (select_gene_function)
        cat("\n")
        # flag0 <- "a variety of gene select function"
        flag1<-paste("基因选择方法：",select_gene_function,sep="")
        write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        # write.table (flag1, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        
        #加入保存Excel
        
        
        #加入保存Excel
        raw_gene_select=raw_gene_select_function
        file_name_0=paste(raw_gene_select,"_",sep="")
        file_name_0_0=paste(file_name_0,select_gene_function,sep="")
        file_name=paste(file_name_0_0,"_cluster.csv",sep="")
        file_name_time=paste(file_name_0_0,"_run_time.csv",sep="")
        file_name_gene_select_counts=paste(file_name_0_0,"_gene_select_counts.csv",sep="")
        # write.table(dataname1,"spearman_4data_log_restore_ne_SIMLR_SIMLR_no_normal_hc_label_cluster.csv",append=T,col.names = FALSE,sep=",",row.names = FALSE)
        write.table(dataname1,file_name,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        write.table(dataname1,file_name_time,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        write.table(dataname1,file_name_gene_select_counts,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        
        gene_select_dir_name<-paste(gene_select_dir,gene_select_names[i],sep="")
        source(gene_select_dir_name)
        # browser()
        # start_time_gene_select <- Sys.time()
        #run gene select
        
        run_select_gene_function<-gene_select(sce,pctkeep)
        
        special_function=select_gene_function
        
        #if function is PCA,break
        if(special_function=='filterPCA')
        {
          dat <- logcounts(sce)
          TrueClu <- length(table(colData(sce)[['label_true']]))
          pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
          # browser()
          
          #notices:this pca dimension set the same as f1000 nPC30
          # pca <- pca$x[, seq_len(params$nPC), drop = FALSE]
          pca<-pca$x[,seq_len(30),drop=FALSE]
          hcl <- hclust(dist(pca), method = "ward.D2")
          cluster <- cutree(hcl, k = TrueClu)
          label.true <- sce$label_true
          label.cluster<-cluster
          
        }else if(special_function=='filterRAFSIL')
        {
          source("gene_select_other_package/RAFSIL/R/RAFSIL.R")
          source("gene_select_other_package/RAFSIL/R/rafsil_FE.R")
          source("gene_select_other_package/RAFSIL/R/rafsil_RF.R")
          source("gene_select_other_package/RAFSIL/R/RAFSIL_utils.R")
          # library(grid)
          # library(gridExtra)
          #install.packages("ClusterR",repos = "https://cran.rstudio.com")
          library(ClusterR)
          #install.packages("randomForest",repos = "https://cran.rstudio.com")
          library(randomForest)
          #pamk()方法
          library(fpc)
          #make.symmetric(X) : could not find function "make.symmetric"
          #but neet package rjags
          # install.packages("rjags",repos = "https://cran.rstudio.com")
          library(rjags)
          #need download:https://sourceforge.net/projects/mcmc-jags/,install,download d/soft/JAGS-4_3_0.exe
          # install.packages("jags",repos = "https://cran.rstudio.com")
          # install.packages("dclone",repos = "https://cran.rstudio.com")
          library(dclone)
          #install.packages("e1071",repos = "https://cran.rstudio.com")
          library(e1071)
          # dat=readRDS('gene_select_other_package/RAFSIL/inst/extdata/example_data.rds')
          # dim(dat$expr)
          # table(dat$labels)
          
          # dat=logcounts(sce)
          # TrueClu <- length(table(colData(sce)[['label_true']]))
          # #Next, derive similarities between the cells using RAFSIL1, RAFSIL2, and Spearman correlation.
          # 
          # #- run RAFSIL1 with 50 forests
          # res.r1 = RAFSIL(t(dat),nrep = 50, method="RAFSIL1",frq=0)
          # res.r2 = RAFSIL(t(dat),           method="RAFSIL2",frq=0)
          # 
          # #- retriev the dissimilarities
          # dis.r1  = res.r1$D
          # dis.r2  = res.r2$D
          # dis.cor = sqrt((1 - cor(dat,method="spearman"))/2)
          # #next plot visual:see :http://www.kostkalab.net/pub_software/RAFSIL/RAFSIL.html
          # hcl <- hclust(dist(dis.cor), method = "ward.D2")
          # cluster <- cutree(hcl, k = TrueClu)
          # label.true <- sce$label_true
          # label.cluster<-cluster
         
           #this change hc
          
          dat=logcounts(sce)
          TrueClu <- length(table(colData(sce)[['label_true']]))
          #Next, derive similarities between the cells using RAFSIL1, RAFSIL2, and Spearman correlation.
          
          #- run RAFSIL1 with 50 forests
          # res.r1 = RAFSIL(t(dat),nrep = 50, method="RAFSIL1",frq=0)
          res.r1 = RAFSIL(t(dat),nrep = 50, method="RAFSIL1",NumC=TrueClu,frq=0)
          label.true <- sce$label_true
          label.cluster<-res.r1$lab
          
        }else if(special_function=='filterSC3')
        {
          
          dat <- sce
          source("gene_select_other_package/SC3/myMethods.R")
          source("gene_select_other_package/SC3/SC3_OK.R")
          # cluster.result <- SC3_OK(inputTags = data, datatype = datatype, SEED=123)
          cluster.result <- SC3_OK_lrychange(inputTags = dat, datatype = datatype, SEED=123)
          label.true <- sce$label_true
          label.cluster<-cluster.result
          #this change hc
        }else 
        {
          flag1<-paste("基因选择%：",pctkeep,sep="")
          write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
          
          flag_not_expressed<-"最终选择的基因和细胞数"
          write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
          gene_length<-paste("处理最终基因长度:",dim(assay(run_select_gene_function,"counts"))[1],sep="")
          cell_length<-paste("处理最终细胞长度:",dim(assay(run_select_gene_function,"counts"))[2],sep="")
          write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
          write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
          write.table(dim(assay(run_select_gene_function,"counts"))[1],file_name_gene_select_counts,append=T,col.names = FALSE,sep=",",row.names = FALSE)
          
          TrueClu <- length(table(colData(sce)[['label_true']]))
          
          hcl <- hclust(dist(t(assay(run_select_gene_function,"logcounts"))), method = "ward.D2")
          label.cluster <- cutree(hcl, k = TrueClu)
          # 根据聚类结果画聚类树及类别框
          # plot(result_hc)
          # label.true <- (as.list(meta["label"]))$label
          
          # label.true <- run_select_gene_function$phenoid
          label.true <- run_select_gene_function$label_true
        }
        
        # 聚类评估
        # 导入聚类评估函数
        source("evalcluster.R")
        NMI_ARI=evalcluster(label.true,label.cluster)
        
        end_time <- Sys.time()
        all_time_temp <-difftime(end_time,start_time,units="mins")
        print(paste("end_time_gene_select_hc:",all_time_temp,sep = ''))
        all_time <-all_time_temp+all_time_pre_data
        print(paste("all_time:",all_time,sep = ''))
        # write.table (dataname1, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        # write.table (run_select_gene_function, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        
        write.table(NMI_ARI,file_name,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        write.table(all_time,file_name_time,append=T,col.names = FALSE,sep=",",row.names = FALSE)
        
        # 
        # flag1<-paste("基因选择运行时间(mins)：",all_time_gene_select,sep="")
        # write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        
        write.table (NMI_ARI, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =T, col.names =F, quote =F,append = T,)
        flag1<-paste("总运行时间(mins)：",all_time,sep="")
        write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
        
      }
      
    }
  }
  
  
  
  
  
  # #只有运行了下面的函数后才有各式各样的过滤指标。
  # genes<-rownames(rowData(sce))
  # #查看是否有线粒体基因，没有，有92个ERCC序列
  # genes[grepl('^MT-',genes)]
  # genes[grepl('^ERCC-',genes)]
  # Identify the remaining ERCC spike-ins.
  
  #考虑线粒体MT，ERCC基因的质量控制
  if(F){
    is.spike <- grepl("^ERCC", rownames(sce))
    sum(is.spike)
    is.mito<-grepl("^MT",rownames(sce))
    sum(is.mito)
    
    table(is.spike)
    summary(colSums(counts(sce[is.spike, ])))
    isSpike(sce, "ERCC") <- is.spike
    
    flag_not_expressed<-"查看基因中线粒体MT，和ERCC序列个数"
    write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
    ERCC_length<-paste("ERCC个数:",sum(is.spike),sep="")
    MT_length<-paste("MT个数:",sum(is.mito),sep="")
    write.table (ERCC_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    write.table (MT_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
    # Calculate QC metrics
    
    sce <- calculateQCMetrics(sce, feature_controls = list(ERCC=is.spike,MT=is.mito))
    
    #去除低于3倍绝对中位差的值。 mad(sce$total_counts)*3。认为
    # 为何设定3倍的MAD？ ：通常把偏离中位数三倍以上的数据作为异常值，和均值标准差方法比，中位数和MAD的计算不受极端异常值的影响，结果更加稳健。
    colData(sce)$libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "lower", log = TRUE)
    
    table(colData(sce)$libsize.drop)
    ggplot(as.data.frame(colData(sce)), aes(x = total_counts)) + 
      geom_histogram(bins = 20, fill = "grey80") + xlab("Total count") + 
      ylab("Number of cells") + 
      geom_vline(xintercept = min(sce$total_counts[!sce$libsize.drop]), 
                 color = "red", linetype = "dashed") + 
      theme_bw()
    
    colData(sce)$feature.drop <- isOutlier(sce$total_features_by_counts, nmads = 3, type = "lower", log = TRUE)
    table(colData(sce)$feature.drop)
    ggplot(as.data.frame(colData(sce)), aes(x = total_features_by_counts)) + 
      geom_histogram(bins = 20, fill = "grey80") + xlab("Number of detected features") + 
      ylab("Number of cells") + 
      geom_vline(xintercept = min(sce$total_features_by_counts[!sce$feature.drop]), 
                 color = "red", linetype = "dashed") + 
      theme_bw()
    
    table(libsize = sce$libsize.drop, feature = sce$feature.drop)
    
    
    
    #We also filter out cells with a large fraction of ERCC reads.
    
    #设置偏离中位数三倍以上的数据为离群值
    colData(sce)$spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads = 3, type = "higher")
    ggplot(as.data.frame(colData(sce)), aes(x = pct_counts_ERCC)) + 
      geom_histogram(bins = 20, fill = "grey80") + xlab("ERCC proportion (%)") + 
      ylab("Number of cells") + 
      geom_vline(xintercept = max(sce$pct_counts_ERCC[!sce$spike.drop]), 
                 color = "red", linetype = "dashed") + 
      theme_bw()
    
    table(sce$spike.drop)
    
    colData(sce)$mito.drop <- isOutlier(sce$pct_counts_MT, nmads = 3, type = "higher")
    ggplot(as.data.frame(colData(sce)), aes(x = pct_counts_MT)) + 
      geom_histogram(bins = 20, fill = "grey80") + xlab("ERCC proportion (%)") + 
      ylab("Number of cells") + 
      geom_vline(xintercept = max(sce$pct_counts_MT[!sce$mito.drop]), 
                 color = "red", linetype = "dashed") + 
      theme_bw()
    
    table(sce$mito.drop)
    #最后只筛选去掉以上四种离群值的数据
    sce <- sce[, !(sce$libsize.drop | sce$feature.drop | sce$spike.drop | sce$mito.drop)]
    dim(sce)
    
  }
  #不考虑线粒体基因的质量控制，主要是删除细胞
  if(F){
    
    
    # Calculate QC metrics
    
    sce <- calculateQCMetrics(sce)
    
    browser()
    #去除低于3倍绝对中位差的值。 mad(sce$total_counts)*3。认为
    # 为何设定3倍的MAD？ ：通常把偏离中位数三倍以上的数据作为异常值，和均值标准差方法比，中位数和MAD的计算不受极端异常值的影响，结果更加稳健。
    colData(sce)$libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "lower", log = TRUE)
    
    table(colData(sce)$libsize.drop)
    
    colData(sce)$feature.drop <- isOutlier(sce$total_features_by_counts, nmads = 3, type = "lower", log = TRUE)
    table(colData(sce)$feature.drop)
    
    table(libsize = sce$libsize.drop, feature = sce$feature.drop)
    
    #最后只筛选去掉以上四种离群值的数据
    sce <- sce[, !(sce$libsize.drop | sce$feature.drop)]
    dim(sce)
    
    flag_not_expressed<-"最后只筛选去掉以上两者离群值的数据"
    write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    gene_length<-paste("基因长度:",dim(assay(sce,"counts"))[1],sep="")
    cell_length<-paste("细胞长度:",dim(assay(sce,"counts"))[2],sep="")
    write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
  }
  #换成最新版本
  if(F){
    
    
    # Calculate QC metrics
    
    sce <- perCellQCMetrics(sce)
    
    browser()
    #去除低于3倍绝对中位差的值。 mad(sce$total_counts)*3。认为
    # 为何设定3倍的MAD？ ：通常把偏离中位数三倍以上的数据作为异常值，和均值标准差方法比，中位数和MAD的计算不受极端异常值的影响，结果更加稳健。
    colData(sce)$libsize.drop <- isOutlier(sce$total, nmads = 3, type = "lower", log = TRUE)
    
    table(colData(sce)$libsize.drop)
    
    colData(sce)$feature.drop <- isOutlier(sce$detected, nmads = 3, type = "lower", log = TRUE)
    table(colData(sce)$feature.drop)
    
    table(libsize = sce$libsize.drop, feature = sce$feature.drop)
    
    #最后只筛选去掉以上四种离群值的数据
    sce <- sce[, !(sce$libsize.drop | sce$feature.drop)]
    dim(sce)
    
    flag_not_expressed<-"最后只筛选去掉以上两者离群值的数据"
    write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    gene_length<-paste("基因长度:",dim(assay(sce,"counts"))[1],sep="")
    cell_length<-paste("细胞长度:",dim(assay(sce,"counts"))[2],sep="")
    write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
  }
  
  if(F){
    flag_not_expressed<-"几种归一化方法不同大小"
    write.table (flag_not_expressed, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
    # gene_length<-paste("预处理之后直接log2(data+1)，max:",max(log2(counts(sce)+1)),sep="")
    # cell_length<-paste("预处理之后直接log2(data+1)，min:",min(log2(counts(sce)+1)),sep="")
    # write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    # write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    # 
    # gene_length<-paste("f1000里面的规范化方法，max:",max(data2),sep="")
    # cell_length<-paste("f1000里面的规范化方法，min:",min(data2),sep="")
    # write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    # write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
    # gene_length<-paste("seurat里面的规范化方法，max:",max(data2),sep="")
    # cell_length<-paste("seurat里面的规范化方法，min:",min(data2),sep="")
    # write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    # write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    # 
    # gene_length<-paste("M3drop里面的规范化方法，max:",max(data2),sep="")
    # cell_length<-paste("M3drop里面的规范化方法，min:",min(data2),sep="")
    # write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    # write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
    gene_length<-paste("M3drop里面的规范化方法,在log，max:",max(data2),sep="")
    cell_length<-paste("M3drop里面的规范化方法，在log，min:",min(data2),sep="")
    write.table (gene_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\t", row.names =F, col.names =F, quote =F,append = T,)
    write.table (cell_length, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
    
  }
  
  # 
  # HVG <-M3Drop::BrenneckeGetVariableGenes(expr_mat=d0, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.8)
  # M3Drop_hvg <- rownames(HVG)
  # d1= d0[M3Drop_hvg,]
  # write.table (dataname1, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # gene_length<-yy[1]
  # gene_length_info=paste("gene_length:",gene_length,sep="")
  # cell_length<-yy[2]
  # cell_length_info=paste("cell_length:",cell_length,sep="")
  # gene_cell_length=paste(gene_length_info,cell_length_info,sep = "  ")
  # write.table (gene_cell_length, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # M3Drop_hvg_length=paste("M3Drop_hvg_length:",length(M3Drop_hvg),sep = "  ")
  # write.table (M3Drop_hvg_length, file =paste("output/","other_infor.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # 
  # 
  # # cluster.result <- SC3_OK(inputTags = d1, datatype = "count", SEED=123)
  # cluster.result <- SC3_OK(inputTags = d1, datatype = datatype, SEED=123)
  # 
  # end_time <- Sys.time()
  # all_time <-difftime(end_time,start_time,units="mins")
  # 
  # 
  # cat("============\n")
  # cat(dataname)
  # cat("\n")
  # source("evalcluster.R")
  # a <- evalcluster(meta$label,cluster.result)
  
  # write.table(a,file_name,append=T,col.names = FALSE,sep=",",row.names = FALSE)
  # write.table(all_time,file_name_time,append=T,col.names = FALSE,sep=",",row.names = FALSE)
  # 
  # #  flag2 <- "G-SC3 运行结果："
  # flag  <- "===================="
  # flag1 <- "--------------------"
  # #  write.table (flag2, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # write.table (flag, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # write.table (dataname, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # write.table (flag1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # write.table (a, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =T, col.names =F, quote =F,append = T,)
  # write.table ("The time is (mins):", file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # write.table (all_time, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  # 
  
  #  write.table (all_time1, file =paste("output/","ARI.txt",sep=""), sep =" ", eol="\n", row.names =F, col.names =F, quote =F,append = T,)
  
  #  write.table (cluster.result, file =paste("output/","cluster_result.txt",sep=""), sep =" ", eol="\n", row.names =T, col.names =F, quote =F,append = F,)
  #  write.table (meta$label, file =paste("output/","meta_label.txt",sep=""), sep =" ", eol="\n", row.names =T, col.names =F, quote =F,append = F,)   
}

#way =1:均值和方差；way =0:方差；way =2:均值；way =3:随机；
# acc_top1(dataname="Biase.Rdata",datatype="FPKM",sizePer=0.3,way=2)



# #big data
# 
# # acc_top1(dataname="Biase.Rdata",datatype="FPKM",pctkeep=10)
# # browser()
acc_top1(dataname="Goolam.Rdata",datatype="CPM",pctkeep=10)#124cell
# browser()
acc_top1(dataname="Deng.Rdata",datatype="CPM",pctkeep=10)#259cell
acc_top1(dataname="Kumar.Rdata",datatype="TPM",pctkeep=10)#268
acc_top1(dataname="KumarTCC.Rdata",datatype="TPM",pctkeep=10)#268
acc_top1(dataname="Trapnell.RData",datatype="TPM",pctkeep=10)#288
acc_top1(dataname="TrapnellTCC.Rdata",datatype="TPM",pctkeep=10)#288
acc_top1(dataname="Wang.Rdata",datatype="TPM",pctkeep=10)#479
acc_top1(dataname="Wallrapp.Rdata",datatype="count",pctkeep=10)#752
acc_top1(dataname="Patel.Rdata",datatype="TPM",pctkeep=10)#864

#
acc_top1(dataname="Haber.RData",datatype="TPM",pctkeep=10)#1522
acc_top1(dataname="Petropoulos.RData",datatype="TPM",pctkeep=10)#1529
acc_top1(dataname="Tasic.Rdata",datatype="TPM",pctkeep=10)#1784
acc_top1(dataname="Klein.Rdata",datatype="count",pctkeep=10)#2717

acc_top1(dataname="Han.Rdata",datatype="count",pctkeep=10)#2746
acc_top1(dataname="Grun.RData",datatype="count",pctkeep=10)#3083
acc_top1(dataname="Cao.Rdata",datatype="count",pctkeep=10)#4186
acc_top1(dataname="Macosko.RData",datatype="TPM",pctkeep=10)#6418


acc_top1(dataname="Zemmour.Rdata",datatype="count",pctkeep=10)#6106
acc_top1(dataname="Spallanzani.Rdata",datatype="count",pctkeep=10)#5287
acc_top1(dataname="Sala.Rdata",datatype="count",pctkeep=10)#10801
# 
# # 
# # 
acc_top1(dataname="Shekhar.Rdata",datatype="count",pctkeep=10)#27499
# acc_top1(dataname="Pbmc68k.Rdata",datatype="count",pctkeep=10)#68579




# acc_top1(dataname="Biase.Rdata",datatype="FPKM",pctkeep=10)
# browser()
# acc_top1(dataname="Deng.Rdata",datatype="CPM",pctkeep=10)
# browser()
# acc_top1(dataname="Goolam.Rdata",datatype="CPM",pctkeep=10)
# acc_top1(dataname="Haber.RData",datatype="TPM",pctkeep=10)
# acc_top1(dataname="Macosko.RData",datatype="TPM",pctkeep=10)
# acc_top1(dataname="Wang.Rdata",datatype="TPM",pctkeep=10)
# acc_top1(dataname="Kumar.Rdata",datatype="TPM",pctkeep=10)
# acc_top1(dataname="KumarTCC.Rdata",datatype="TPM",pctkeep=10)
# acc_top1(dataname="Trapnell.RData",datatype="TPM",pctkeep=10)
# acc_top1(dataname="TrapnellTCC.Rdata",datatype="TPM",pctkeep=10)
# # acc_top1(dataname="Chung.Rdata",datatype="TPM",pctkeep=10)
# # acc_top1(dataname="Liu.Rdata",datatype="RPKM",pctkeep=10)
# # acc_top1(dataname="Trapnell.Rdata",datatype="FPKM",pctkeep=10)
# # acc_top1(dataname="Usoskin.RData",datatype="TPM",pctkeep=10)
# #注意，以下的datatype可能不对
# acc_top1(dataname="Cao.Rdata",datatype="count",pctkeep=10)
# acc_top1(dataname="Han.Rdata",datatype="count",pctkeep=10)
# acc_top1(dataname="Klein.Rdata",datatype="count",pctkeep=10)
# acc_top1(dataname="Pbmc68k.Rdata",datatype="count",pctkeep=10)
# acc_top1(dataname="Shekhar.Rdata",datatype="count",pctkeep=10)
# acc_top1(dataname="Wallrapp.Rdata",datatype="count",pctkeep=10)
# acc_top1(dataname="Patel.Rdata",datatype="TPM",pctkeep=10)
# acc_top1(dataname="Tasic.Rdata",datatype="TPM",pctkeep=10)
# acc_top1(dataname="Petropoulos.RData",datatype="TPM",pctkeep=10)
# acc_top1(dataname="Grun.RData",datatype="count",pctkeep=10)









