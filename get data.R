k <- 12

T_used<- 96

# only pred_h = 0
pred_h<-0

library(forecast)
library(readr)


# setwd("C:/Users/subing/OneDrive - The University of Hong Kong - Connect/桌面/monthly data all")


#Short_term_interest_rates <- read_csv("Short-term interest rates.csv")
# 
#Long_term_interest_rates <- read_csv("Long-term interest rates.csv")
#CPI_total <- read_csv("CPI_total.csv")
#total_production <- read_csv("total production.csv")
# Construction <- read_csv("Construction.csv")
# PPI_manuf <- read_csv("PPI_manuf.csv")
# PPI_total <- read_csv("PPI_total.csv")
# CPI_food <- read_csv("CPI_food.csv")
# CPI_other <- read_csv("CPI_other.csv")

PPI_domestic <- read_csv("PPI_domestic.csv")

CPI_energy <- read_csv("CPI_energy.csv")
Manufacturing <- read_csv("Manufacturing.csv")


PPI_domestic <- na.omit(PPI_domestic)


colname<-Reduce(intersect, list(CPI_energy$Location, PPI_domestic$Location, #PPI_total$Location,
                                Manufacturing$Location))


# colname<-setdiff(colname, c("Estonia"))

library('dse')

fun_data<-function(data,difflog = TRUE){
  data<-as.data.frame(data)
  
  data<-data[data[,1] %in% colname,]
  
  # remove location
  data<- data[,-1]
  
  data<- as.matrix(data)
  data<-apply(data,2, as.numeric)
  
  for(i in 1:dim(data)[1]){
    
    data_i <-as.numeric(data[i,])
    if(difflog){
      data_i_diff<-diff(log(data_i))
    }else{
      data_i_diff<-diff((data_i))
    }
    
    if(is.na(mean(data_i_diff,na.rm=T))){
      data_i_diff[is.na(data_i_diff)] <- mean(data,na.rm=T)
    }else{
      data_i_diff[is.na(data_i_diff)] <- mean(data_i_diff,na.rm=T)
    }
    
    # remove 0
    data[i,-1] <- data_i_diff
    
  }
  
  data<- data[,-1]
  
  # 10 years data
  return(data[,((12 * (k -1) + 1) + 5 ):( (12 * (k -1 + 10)) + 6)])
  #return(data[,((12 * (k -1) + 1) - 5 ):( (12 * (k -1 + 10))  )])
  
}

# total has 24 years data 

# Share_prices_<-fun_data(Share_prices)

# CPI_total<-fun_data(CPI_total)

# Short_term_interest_rates<-fun_data(Short_term_interest_rates, FALSE)
# 
# Long_term_interest_rates<-fun_data(Long_term_interest_rates, FALSE)
# total_production<- fun_data(total_production)

# PPI_manuf <- fun_data(PPI_manuf)

# CPI_total <-  fun_data(CPI_total)

# Construction <- fun_data(Construction)
# CPI_food <-fun_data(CPI_food)
#CPI_other<-fun_data(CPI_other)

# PPI_total <-  fun_data(PPI_total)

PPI_domestic <- fun_data(PPI_domestic)
CPI_energy <-fun_data(CPI_energy)
Manufacturing<-fun_data(Manufacturing)



ymat<- PPI_domestic[,-1]

T_len <- dim(ymat)[2] + 1

xmat1<- CPI_energy[,- T_len] - apply(CPI_energy[,- T_len], 1, mean)

xmat4<- Manufacturing[,- T_len]  - apply(Manufacturing[,- T_len], 1, mean)



# xmat1<- CPI_food[,- T_len] - apply(CPI_food[,- T_len], 1, mean)
# xmat2<-  CPI_other[,- T_len]  - apply(CPI_other[,- T_len], 1, mean)
#xmat3<-  Short_term_interest_rates[,- T_len]  -  apply(Short_term_interest_rates[,- T_len], 1, mean)
#Long_term_interest_rates[,- T_len]  -  apply(Long_term_interest_rates[,- T_len], 1, mean)