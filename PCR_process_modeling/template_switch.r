###############################################################################################
# Step1 preparate matrix
###############################################################################################
rm(list=ls())
gc()

library(optparse)
option_list <- list(
  make_option(c("--error_rate"), type = "double", default=FALSE,help="Provide the probability that a template switch may occur "),
  make_option(c("--cycle"), type="integer", default=FALSE,help="Number of PCR cycles"),
  make_option(c("--output"), type="character", default=FALSE,help="Output file name"),
  make_option(c("--len"),type="integer", default=FALSE,help="Length of reads"),
  make_option(c("--reads_nums"),type="integer", default=FALSE,help="Initial number of PCR molecules"),
  make_option(c("--throughput"),type="integer", default=FALSE,help="0:Not Simulation throughput ; 1:Simulation throughpu"),
  make_option(c("--throughput_reads"),type="integer", default=FALSE,help="If simulation throughput,Provide throughput reads nums"),
  make_option(c("--wildtype_ratio"),type="double", default=FALSE,help="Initial ratio of wildtype_ratio"),
  #    make_option(c("--identity"),type="double", default=FALSE,help="Genotype identity"),
  make_option(c("--mode"),type="integer", default=FALSE,help="0:Simulation will calculate each cycle;1:Simulation will calculate only 20/30/40 cycle"))
opt_parser=OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
Cycle=as.integer(opt$cycle)
Output=opt$output
Len=as.integer(opt$len)
Error2=opt$error_rate
Error=Error2*Len
Reads_nums=as.integer(opt$reads_nums)
Wildtype_ratio=opt$wildtype_ratio
Othertype_ratio=1-Wildtype_ratio
Mode=opt$mode
Throughput=opt$throughput
Throughput_reads=opt$throughput_reads
#Identity=opt$identity

kk <- 1 
table <- matrix(NA, nrow = 6, ncol = 40)
mode_table <- matrix(NA, nrow = 6, ncol = 3)
Wildtype_reads_sum_num=Reads_nums*Wildtype_ratio #初始wildtype 分配到的reads总数
Othertype_reads_sum_num=Reads_nums*Othertype_ratio#初始othertype 分配到的reads总数
Othertype_row_num=Len*3 #othertype基因型总数，一种基因型占一行
Per_Othertype_num= round(Othertype_reads_sum_num/Othertype_row_num)#othertype每个基因型初始有多少个reads
Wildtype_row_num= round(Wildtype_reads_sum_num/Per_Othertype_num)#wildtye需要折叠多少行
adjust_umi_num=Wildtype_row_num*Per_Othertype_num + Othertype_row_num*Per_Othertype_num #调整后的umi数量
row_num=Othertype_row_num+Wildtype_row_num#矩阵总行数
col_num=Per_Othertype_num*2+2 #矩阵总列数，第一列为ID，第二列为对应基因型当前总的reads数
matrix_num=row_num*col_num
Per_Othertype_begin=Per_Othertype_num+2
compute_begin=Per_Othertype_num+3
weip_other_begin=Wildtype_row_num-1
wei_other_begin=Wildtype_row_num+1
#fix_error_rate=Identity*Error

m <- matrix(data = rep(0, matrix_num), nrow = row_num, ncol = col_num)#A
colnames(m) <- c(1:col_num)
rownames(m) <- c(1:row_num)
for (i in 2:4) {
  assign(paste0("m", i), m)
}
m[, 1] <- 1:row_num #A
m[, 2] <- Per_Othertype_num
m[, 3:Per_Othertype_begin] <- 1
m1 <- m
l <- 0
ll <- 0
m5 <- matrix(data = rep(0, matrix_num), nrow = row_num, ncol = col_num)
m5[, 1] <- 1:row_num
colnames(m5) <- c(1:col_num)
rownames(m5) <- c(1:row_num)
print("ok")
##############################################################################################
#Step2 Enzyme function
##############################################################################################
double <- function(x) {2 * as.numeric(x)}
# Multiplication factors for cycles > 10
if(Len==25){
  multiplication_factors <- c(2^1,2^1,2^1,2^1,2^1,2^1,2^1,2^1,2^1,2^1    #cycle 10-20
                              , 2^0.41,2^0.41,2^0.41,2^0.41,2^0.41                                   #cycle 20-25
                              , 2^0.408,2^0.408,2^0.408,2^0.408,2^0.408                              #cycle 25-30
                              , 2^0.06,2^0.06,2^0.06,2^0.06,2^0.06                                   #cycle 30-35     
                              , 2^(-0.07) , 2^(-0.07) , 2^(-0.07) , 2^(-0.07) , 2^(-0.07)            #cycle 35-40
  )
}else if(Len==250){
  multiplication_factors <- c(2^1,2^1,2^1,2^1,2^1,2^1,2^1,2^1,2^1,2^1    #cycle 10-20
                              , 2^0.524,2^0.524,2^0.524,2^0.524,2^0.524                                   #cycle 20-25
                              , 2^0.166,2^0.166,2^0.166,2^0.166,2^0.166                                   #cycle 25-30
                              , 2^0.002,2^0.002,2^0.002,2^0.002,2^0.002                                                  #cycle 30-35
                              , 2^0.13,2^0.13,2^0.13,2^0.13,2^0.13                                   #cycle 35-40
  )
}else if(Len==750){
  multiplication_factors <- c(2^0.848 , 2^0.848 , 2^0.848 , 2^0.848 , 2^0.848#cycle 10-15
                              , 2^0.864 , 2^0.864 , 2^0.864 , 2^0.864 , 2^0.864                          #cycle 15-20
                              , 2^0.402 , 2^0.402 , 2^0.402 , 2^0.402 , 2^0.402                               #cycle 20-25
                              , 2^0.26 , 2^0.26 , 2^0.26 , 2^0.26 , 2^0.26                          #cycle 25-30
                              , 2^(-0.024) , 2^(-0.024) , 2^(-0.024) , 2^(-0.024) , 2^(-0.024)           #cycle 30-35
                              , 2^0.12 , 2^0.12 , 2^0.12 , 2^0.12 , 2^0.12                          #cycle 35-40
  )
}else{
  multiplication_factors <- c(2^0.848 , 2^0.848 , 2^0.848 , 2^0.848 , 2^0.848#cycle 10-15
                              , 2^0.864 , 2^0.864 , 2^0.864 , 2^0.864 , 2^0.864                          #cycle 15-20
                              , 2^0.402 , 2^0.402 , 2^0.402 , 2^0.402 , 2^0.402                               #cycle 20-25
                              , 2^0.26 , 2^0.26 , 2^0.26 , 2^0.26 , 2^0.26                          #cycle 25-30
                              , 2^(-0.024) , 2^(-0.024) , 2^(-0.024) , 2^(-0.024) , 2^(-0.024)           #cycle 30-35
                              , 2^0.12 , 2^0.12 , 2^0.12 , 2^0.12 , 2^0.12                          #cycle 35-40
  )
}

##############################################################################################
#Step3 Cycle
##############################################################################################
cycle <- Cycle  # Set the number of cycles
for (c in 1:cycle) {
  m <- m1
  kk <- kk + 1 
  a <- sum(as.numeric(m[, 2]))
  #  n1 <- rbinom(1, a, fix_error_rate) 
  n1 <- rbinom(1, a, Error)
  wei_row=Othertype_row_num+1
  wei <- matrix(data = rep(0, wei_row), nrow = wei_row, ncol = 1)
  total_weight <- sum(as.numeric(m[1:Wildtype_row_num, 2]))
  wei[1, 1] <- total_weight
  wei[2:wei_row, 1] <- as.numeric(m[wei_other_begin:row_num, 2])
  
  if (c <= 10) {
    m1[, 2:col_num] <- apply(m1[, 2:col_num], 2, double)
  } else {
    factor_index <- c - 10
    if (factor_index <= length(multiplication_factors)) {
      factor <- multiplication_factors[factor_index]
      m1[, 2:col_num] <- round(apply(m1[, 2:col_num], 2, function(x) as.numeric(x) * factor))
    }
  }
  
  if (n1 == 0) {
    for (i in 2:4) {
      assign(paste0("m", i), m5)
    }
    
    k <- 0
    l <- l + 1
    table[1, l] <- k
    table[2, l] <- 0
    table[3, l] <- adjust_umi_num
  } 
  else 
  {
    
    ak <- sample(c(1:wei_row), size = n1, replace = TRUE, prob = wei)
    ss <- as.data.frame(table(ak))
    ss <- as.matrix(ss)
    
    for (i in 1:nrow(ss)) {
      p <- as.numeric(ss[i, 1])
      np <- as.numeric(ss[i, 2])
      if (p==1) { 
        wei_wt <- matrix(data=rep(0,Wildtype_row_num),nrow=Wildtype_row_num,ncol=1)
        for (i in 1:Wildtype_row_num){
          wei_wt[i, 1]<- as.numeric(m[i, 2])
        }
        ak_wt <- sample(c(1:Wildtype_row_num),size = np,replace = TRUE, prob = wei_wt)
        ss_wt <- as.data.frame(table(ak_wt))
        ss_wt <- as.matrix(ss_wt) 
        for (ii in 1:nrow(ss_wt)){
          p_wt <- as.numeric(ss_wt[ii,1]) 
          np_wt <- as.numeric(ss_wt[ii,2]) 
          weip_wt <- as.numeric(m[p_wt, 3:col_num])
          akp_wt <- sample(c(3:col_num), size = np_wt, replace = TRUE, prob = weip_wt)
          ssp_wt <- as.data.frame(table(akp_wt))
          ssp_wt <- as.matrix(ssp_wt)
          m2[p_wt, as.numeric(ssp_wt[, 1])] <- -as.numeric(ssp_wt[, 2])
          m2[p_wt, 2] <- -as.numeric(np_wt)
        }
      }else{
        weip <- as.numeric(m[p+weip_other_begin, 3:col_num])
        akp <- sample(c(3:col_num), size = np, replace = TRUE, prob = weip)
        ssp <- as.data.frame(table(akp)) 
        ssp <- as.matrix(ssp) 
        m2[p+weip_other_begin, as.numeric(ssp[, 1])] <- -as.numeric(ssp[, 2]) 
        m2[p+weip_other_begin, 2] <- -as.numeric(np) 
      }
    }
    ak2 <- sample(c(1:wei_row), size = n1, replace = TRUE, prob = wei)
    ss2 <- as.data.frame(table(ak2))
    ss2 <- as.matrix(ss2)
    for (i in 1:nrow(ss2)) {
      p <- as.numeric(ss2[i, 1])
      np <- as.numeric(ss2[i, 2]) 
      if (p==1) { 
        wei_wt <- matrix(data=rep(0,Wildtype_row_num),nrow=Wildtype_row_num,ncol=1)
        for (i in 1:Wildtype_row_num){
          wei_wt[i, 1]<- as.numeric(m[i, 2])
        }          
        ak_wt <- sample(c(1:Wildtype_row_num),size = np,replace = TRUE, prob = wei_wt)
        ss_wt <- as.data.frame(table(ak_wt))
        ss_wt <- as.matrix(ss_wt) 
        for (ii in 1:nrow(ss_wt)){
          p_wt <- as.numeric(ss_wt[ii,1])
          np_wt <- as.numeric(ss_wt[ii,2]) 
          weip_wt <- as.numeric(m[p_wt, 3:col_num])
          akp_wt <- sample(c(3:col_num), size = np_wt, replace = TRUE, prob = weip_wt)
          ssp_wt <- as.data.frame(table(akp_wt))
          ssp_wt <- as.matrix(ssp_wt)
          m3[p_wt, as.numeric(ssp_wt[, 1])] <- as.numeric(m3[p_wt, as.numeric(ssp_wt[, 1])]) - as.numeric(ssp_wt[, 2])
          for (jj in 1:nrow(ssp_wt)) {
            if (3 <= as.numeric(ssp_wt[jj, 1]) & as.numeric(ssp_wt[jj, 1]) <= Per_Othertype_begin) {
              ssp_wt[jj, 1] <- as.numeric(ssp_wt[jj, 1]) + Per_Othertype_num
            }
          } 
          m3[p_wt, as.numeric(ssp_wt[, 1])] <- as.numeric(m3[p_wt, as.numeric(ssp_wt[, 1])]) + as.numeric(ssp_wt[, 2])
        }
      }else{
        weip <- as.numeric(m[p+weip_other_begin, 3:col_num]) 
        akp <- sample(c(3:col_num), size = np, replace = TRUE, prob = weip) 
        ssp <- as.data.frame(table(akp)) 
        ssp <- as.matrix(ssp) 
        m3[p+weip_other_begin, as.numeric(ssp[, 1])] <- as.numeric(m3[p+weip_other_begin, as.numeric(ssp[, 1])]) - as.numeric(ssp[, 2]) 
        
        for (j in 1:nrow(ssp)) {
          if (3 <= as.numeric(ssp[j, 1]) & as.numeric(ssp[j, 1]) <= Per_Othertype_begin) {
            ssp[j, 1] <- as.numeric(ssp[j, 1]) + Per_Othertype_num
          }
        }
        
        m3[p+weip_other_begin, as.numeric(ssp[, 1])] <- as.numeric(m3[p+weip_other_begin, as.numeric(ssp[, 1])]) + as.numeric(ssp[, 2])
      }
    }
    
    m4 <- apply(m3, 2, as.numeric) + apply(m2, 2, as.numeric)
    m4[, 1] <- 1:row_num
    m1 <- apply(m1, 2, as.numeric) + apply(m4, 2, as.numeric)
    m1[, 1] <- 1:row_num
    for (i in 2:4) {
      assign(paste0("m", i), m5)
    }
    
    if (Throughput==1){
      if (Mode==1){    
        sub_m1 <- m1[, 3:ncol(m1)]
        nrow_m1 <- nrow(sub_m1)
        ncol_m1 <- ncol(sub_m1)
        weights <- as.vector(sub_m1)
        indices <- seq_along(weights)
        samples <- sample(indices, size = Throughput_reads, replace = TRUE, prob = weights)
        count_matrix <- matrix(0, nrow = nrow_m1, ncol = ncol_m1)
        sample_table <- table(samples)
        count_matrix[as.integer(names(sample_table))] <- as.integer(sample_table)
        
        TOO_count=0
        tOM_count=0
        major_umi_num=0
        minor_umi_num=0
        if (c %in% c(20, 30, 40)){
          for (i in 1:nrow(count_matrix)) {
            for (j in 1:ncol(count_matrix)) {
              if (count_matrix[i, j] != 0 && count_matrix[i, j] > 2) { if ((j + Per_Othertype_num <= ncol(count_matrix) && count_matrix[i, j + Per_Othertype_num] == 0)|| (j - Per_Othertype_num > 0 && count_matrix[i, j - Per_Othertype_num] == 0))
              {TOO_count <- TOO_count+1}}}}
          
          for (i in 1:nrow(count_matrix)) {
            for (j in 1:ncol(count_matrix)) {
              if (count_matrix[i, j] != 0) {
                if ((j + Per_Othertype_num <= ncol(count_matrix) &&count_matrix[i, j + Per_Othertype_num] != 0 &&(count_matrix[i, j] + count_matrix[i, j + Per_Othertype_num] > 2)) ||
                    (j - Per_Othertype_num > 0 &&count_matrix[i, j - Per_Othertype_num] != 0 &&(count_matrix[i, j] + count_matrix[i, j - Per_Othertype_num] > 2))
                )
                {tOM_count <- tOM_count+1}}}}  
          
          #计算major和minor的数量
          for (i in 1:nrow(count_matrix)) {
            for (j in 1:ncol(count_matrix)) {
              if (count_matrix[i, j] != 0) {
                
                # Check right side
                right_j <- j + Per_Othertype_num
                if (right_j <= ncol(count_matrix)) {
                  if (count_matrix[i, right_j] < count_matrix[i, j] && count_matrix[i, right_j] + count_matrix[i, j] > 2) {
                    major_umi_num <- major_umi_num + count_matrix[i, j]
                    minor_umi_num <- minor_umi_num + count_matrix[i, right_j]
                  }
                }
                
                # Check left side
                left_j <- j - Per_Othertype_num
                if (left_j > 0) {
                  if (count_matrix[i, left_j] < count_matrix[i, j] && count_matrix[i, left_j] + count_matrix[i, j] > 2) {
                    major_umi_num <- major_umi_num + count_matrix[i, j]
                    minor_umi_num <- minor_umi_num + count_matrix[i, left_j]
                  }
                }
                
              }
            }
          }
          
          
          TOM_count<- tOM_count/2
          tho_ratio <- TOM_count / TOO_count
          #    major_umi_num_ratio <- major_umi_num/(major_umi_num+minor_umi_num)
          minor_umi_num_ratio <- minor_umi_num/(major_umi_num+minor_umi_num)
          ll <- ll + 1
          mode_table[1, ll] <- tho_ratio
          mode_table[2, ll] <- TOM_count
          mode_table[3, ll] <- TOO_count
          mode_table[4, ll] <- minor_umi_num_ratio
          mode_table[5, ll] <- major_umi_num
          mode_table[6, ll] <- minor_umi_num
        }}else if(Mode==0){
          #    print("Mode:Throughput")
          sub_m1 <- m1[, 3:ncol(m1)]
          nrow_m1 <- nrow(sub_m1)
          ncol_m1 <- ncol(sub_m1)
          weights <- as.vector(sub_m1)
          indices <- seq_along(weights)
          samples <- sample(indices, size = Throughput_reads, replace = TRUE, prob = weights)
          count_matrix <- matrix(0, nrow = nrow_m1, ncol = ncol_m1)
          sample_table <- table(samples)
          count_matrix[as.integer(names(sample_table))] <- as.integer(sample_table)
          TOO_count=0
          tOM_count=0
          major_umi_num=0
          minor_umi_num=0
          for (i in 1:nrow(count_matrix)) {# 遍历每一列（只到第 1000 列，因为超出 1000 列后没有对应的 1000 列之后的单元格）
            for (j in 1:ncol(count_matrix)) {# 检查当前位置和1000列以后的对应位置是否都不为0
              if (count_matrix[i, j] != 0 && count_matrix[i, j] > 2) { if ((j + Per_Othertype_num <= ncol(count_matrix) && count_matrix[i, j + Per_Othertype_num] == 0)|| (j - Per_Othertype_num > 0 && count_matrix[i, j - Per_Othertype_num] == 0))
              {TOO_count <- TOO_count+1}}}}
          
          for (i in 1:nrow(count_matrix)) {# 遍历每一列（只到第 1000 列，因为超出 1000 列后没有对应的 1000 列之后的单元格）
            for (j in 1:ncol(count_matrix)) {# 检查当前位置和1000列以后的对应位置是否都不为0
              if (count_matrix[i, j] != 0) {
                if ((j + Per_Othertype_num <= ncol(count_matrix) &&count_matrix[i, j + Per_Othertype_num] != 0 &&(count_matrix[i, j] + count_matrix[i, j + Per_Othertype_num] > 2)) ||
                    (j - Per_Othertype_num > 0 &&
                     count_matrix[i, j - Per_Othertype_num] != 0 &&
                     (count_matrix[i, j] + count_matrix[i, j - Per_Othertype_num] > 2))
                )
                {tOM_count <- tOM_count+1}}}}  
          
          for (i in 1:nrow(count_matrix)) {
            for (j in 1:ncol(count_matrix)) {
              if (count_matrix[i, j] != 0) {
                
                # Check right side
                right_j <- j + Per_Othertype_num
                if (right_j <= ncol(count_matrix)) {
                  if (count_matrix[i, right_j] < count_matrix[i, j] && count_matrix[i, right_j] + count_matrix[i, j] > 2) {
                    major_umi_num <- major_umi_num + count_matrix[i, j]
                    minor_umi_num <- minor_umi_num + count_matrix[i, right_j]
                  }
                }
                
                # Check left side
                left_j <- j - Per_Othertype_num
                if (left_j > 0) {
                  if (count_matrix[i, left_j] < count_matrix[i, j] && count_matrix[i, left_j] + count_matrix[i, j] > 2) {
                    major_umi_num <- major_umi_num + count_matrix[i, j]
                    minor_umi_num <- minor_umi_num + count_matrix[i, left_j]
                  }
                }
                
              }
            }
          }
          
          TOM_count<- tOM_count/2
          tho_ratio <- TOM_count / TOO_count
          minor_umi_num_ratio <- minor_umi_num/(major_umi_num+minor_umi_num)
          
          l <- l + 1
          table[1, l] <- tho_ratio
          table[2, l] <- TOM_count
          table[3, l] <- TOO_count
          table[4, l] <- minor_umi_num_ratio
          table[5, l] <- major_umi_num
          table[6, l] <- minor_umi_num
        }}else if (Throughput==0)  {
          if (Mode==1){
            #    print("Mode:stand&mode")
            if (c %in% c(20, 30, 40)){
              h <- 0
              major_umi_num=0
              minor_umi_num=0
              for (i in 1:nrow(m1)) {
                for (j in compute_begin:col_num) {
                  if (m1[i, j] == 0) {
                    h <- h + 1
                  }
                }
              }
              sub_m1 <- m1[, 3:ncol(m1)]
              for (i in 1:nrow(sub_m1)) {
                for (j in 1:ncol(sub_m1)) {
                  if (sub_m1[i, j] != 0) {
                    
                    # Check right side
                    right_j <- j + Per_Othertype_num
                    if (right_j <= ncol(sub_m1)) {
                      if (sub_m1[i, right_j] < sub_m1[i, j] && sub_m1[i, right_j] + sub_m1[i, j] > 2) {
                        major_umi_num <- major_umi_num + sub_m1[i, j]
                        minor_umi_num <- minor_umi_num + sub_m1[i, right_j]
                      }
                    }
                    
                    # Check left side
                    left_j <- j - Per_Othertype_num
                    if (left_j > 0) {
                      if (sub_m1[i, left_j] < sub_m1[i, j] && sub_m1[i, left_j] + sub_m1[i, j] > 2) {
                        major_umi_num <- major_umi_num + sub_m1[i, j]
                        minor_umi_num <- minor_umi_num + sub_m1[i, left_j]
                      }
                    }
                    
                  }
                }
              }    
              
              
              k <- (adjust_umi_num - h) / h
              kkk=adjust_umi_num - h
              minor_umi_num_ratio <- minor_umi_num/(major_umi_num+minor_umi_num)
              ll <- ll + 1
              mode_table[1, ll] <- k
              mode_table[2, ll] <- kkk
              mode_table[3, ll] <- h
              mode_table[4, ll] <- minor_umi_num_ratio
              mode_table[5, ll] <- major_umi_num
              mode_table[6, ll] <- minor_umi_num    
            }}else if (Mode==0){
              #    print("Mode:stand")
              h <- 0
              major_umi_num=0
              minor_umi_num=0
              for (i in 1:nrow(m1)) {
                for (j in compute_begin:col_num) {
                  if (m1[i, j] == 0) {
                    h <- h + 1
                  }
                }
              }
              sub_m1 <- m1[, 3:ncol(m1)]
              for (i in 1:nrow(sub_m1)) {
                for (j in 1:ncol(sub_m1)) {
                  if (sub_m1[i, j] != 0) {      
                    # Check right side
                    right_j <- j + Per_Othertype_num
                    if (right_j <= ncol(sub_m1)) {
                      if (sub_m1[i, right_j] < sub_m1[i, j] && sub_m1[i, right_j] + sub_m1[i, j] > 2) {
                        major_umi_num <- major_umi_num + sub_m1[i, j]
                        minor_umi_num <- minor_umi_num + sub_m1[i, right_j]
                      }
                    }
                    
                    # Check left side
                    left_j <- j - Per_Othertype_num
                    if (left_j > 0) {
                      if (sub_m1[i, left_j] < sub_m1[i, j] && sub_m1[i, left_j] + sub_m1[i, j] > 2) {
                        major_umi_num <- major_umi_num + sub_m1[i, j]
                        minor_umi_num <- minor_umi_num + sub_m1[i, left_j]
                      }
                    }
                    
                  }
                }
              }    
              k <- (adjust_umi_num - h) / h
              kkk=adjust_umi_num - h
              minor_umi_num_ratio <- minor_umi_num/(major_umi_num+minor_umi_num)
              l <- l + 1
              table[1, l] <- k
              table[2, l] <- kkk
              table[3, l] <- h
              table[4, l] <- minor_umi_num_ratio
              table[5, l] <- major_umi_num
              table[6, l] <- minor_umi_num    
            }}
    
  }
}

##################################################################################
#Write result
#################################################################################
if (Mode==1){
  my_table <- mode_table
}else if (Mode==0){
  my_table <- table}
file_path <- Output
write.table(my_table, file = file_path, row.names = FALSE, quote = FALSE, sep = "\t")
print("Simulation finished!")



