# set working directory
setwd('~/desktop/2/')

library(survival)

# this function can divided sample into 'Low' or 'High' group
# data of RPKM is needed
group_divide <- function(val, treshold) {
    group <- ifelse (val <= treshold, 'Low', 'High')
  return (group)
}

# read data from file
raw_db <- read.table('./os/combined.tumor-RPKM.txt',
                      sep = '\t',
                      header = F,
                      as.is = T)

surv_db <- read.table('./os/survival-time.txt',
                      sep = '\t',
                      header = T,
                      as.is = T)

# reshape the data set to fufiil futher analysis requirement
RPKM_db <- data.frame(raw_db, row.names = raw_db[,1])
RPKM_db <- data.frame(t(RPKM_db[, -1]), stringsAsFactors = F)

# define a set of function to grab infomation from the other file(survival-time.xls)
search_time <- function(ID){
  surv_time <- surv_db[surv_db$ID == ID, 'survial_time']
  return(surv_time)
}
search_cut <- function(ID){
  surv_cut <- surv_db[surv_db$ID == ID, 'cut']
  return(surv_cut)
}
# add survial infomation to RPKM data frame
RPKM_db$surv_time <- sapply(RPKM_db$ID, search_time)
RPKM_db$surv_cut <- sapply(RPKM_db$ID, search_cut)

rm(raw_db)
# initialize a data frame to store final result
result_tab <- data.frame(sample_id = character(),
                        Low_count = numeric(),
                        High_count = numeric(),
                        p_value = numeric(),
                        significance = character())
# used to mapping column
col_count <- length(names(RPKM_db))

# use for circle to analysis each lncRNA's data
for (i in names(RPKM_db)[c(-1, -(col_count), -(col_count - 1))]){
  # extract subset of every lncRNA's RPKM(as data frame form)
  tmp_db <- data.frame(RPKM_db[, c('ID', i, 'surv_time', 'surv_cut')])
  # modify data type to numeric for further analysis
  tmp_db[, i] <- as.numeric(tmp_db[, i])
  # creat 'mid' to store median
  mid <- median(tmp_db[, i])
  # use sapply to handle a colum of data(divide group)
  tmp_db$Group <- sapply(tmp_db[, i], group_divide, treshold = mid)
  # creat a survial object 
  surv_data <- Surv(time = as.numeric(tmp_db$surv_time),
                    event = as.numeric(tmp_db$surv_cut),
                    type = "right"
                    ) 
  # do survial analysis
  surv_test <- survdiff(surv_data~tmp_db$Group, data = tmp_db)
  # prepare output result
  sample_id = i
  p_value <- 1 - pchisq(surv_test$chisq, 1) # df was set to be constant
  p_significance <- ifelse ( p_value <= 0.05, '+', '-')
  L_count <- length(tmp_db[tmp_db$Group == 'Low', i])
  H_count <- length(tmp_db[tmp_db$Group == 'High', i])
  # add result for this lncRNA to result data frame
  result_tab <- rbind(result_tab,data.frame(sample_id, L_count,
                                            H_count, p_value,
                                            p_significance))
  # when p < 0.05 draw a k-m plot for this subset of data
  if (p_value < 0.05) {
    # file was named after the gene ID
    png(filename = paste(i, "png", sep="."), 
        width = 800,
        height = 800)
        plot(survfit(surv_data~tmp_db$Group, data = tmp_db),col=c("red","green"))  
        title(i) 
        legend(1800,0.995,legend=paste('p.value = ',p_value,sep=''),bty='n',cex=1.4)
        legend(max(as.numeric(as.character(tmp_db$surv_time)),na.rm = T)*0.7,0.94,legend=c(paste('High=',H_count),paste('Low=',L_count)),bty='n',cex=1.3,lwd=3,col=c('red','green'))
        
    dev.off()
  }
}
write.table(result_tab,file="result.txt",sep="\t",row.names = FALSE,quote = FALSE)