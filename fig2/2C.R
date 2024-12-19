# 安装必要的包
if (!requireNamespace("karyoploteR", quietly = TRUE)) {
  install.packages("karyoploteR")
}
library(karyoploteR)

# 加载数据
mutation_data <- read.csv("data.2C.csv")  # 确保文件格式正确

# 定义颜色映射
mutation_colors <- c(
  "Probably Deleterious" = "red",
  "Possibly Deleterious" = "purple",
  "Benign" = "green"
)

# 构建基因组图谱
kp <- plotKaryotype(genome="hg19", plot.type=1)  # 根据需求更换基因组版本

# 添加突变数据到染色体图谱，形状为三角形
for (mutation_type in unique(mutation_data$Type)) {
  type_data <- subset(mutation_data, Type == mutation_type)
  # 使用自定义符号绘制三角形
  kpPolygon(kp,
            chr = type_data$Chromosome,
            x0 = type_data$Start - 1,  # 三角形左边顶点
            x1 = type_data$Start + 1,  # 三角形右边顶点
            y0 = rep(0, nrow(type_data)),  # 底边 y 值
            y1 = rep(0.5, nrow(type_data)),  # 顶部 y 值
            col = mutation_colors[mutation_type],
            border = NA)
}

# 添加图例
legend("topright", legend = names(mutation_colors), col = mutation_colors, pch=17)  # pch=17 是三角形标记
