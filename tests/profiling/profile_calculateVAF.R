library(profvis)
library(SVseek)

profvis({
    calculateVAF(data_fusion$gr, data_fusion$breakpoint, "qiaseq")
})
