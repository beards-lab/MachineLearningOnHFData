rm(list = ls())
library(survival)
library(readxl)
# 从Excel文件加载数据，假设时间列名为 "Time"，事件指示符列名为 "Event"
# 以下是假定的列名，您需替换为您数据的具体列名
data <- read_excel("coxregression.xlsx")

# 创建一个Surv对象，它结合了生存时间和事件发生情况
surv_object <- Surv(time = data$Time, event = data$Event)


