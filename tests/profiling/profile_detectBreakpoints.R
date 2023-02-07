library(profvis)
library(SVseek)

profvis({
    detectBreakpoints(data_fusion$gr, data_fusion$genes)
})
