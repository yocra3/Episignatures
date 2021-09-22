library(ggplot2)
 a <- read.table("cot.txt")
 a$step <- 1:4
 a$run <- rep(seq_len(nrow(a)/4) , each = 4)
   a$config <- rep(seq_len(nrow(a)/20) , each = 20)

ggplot(subset(a, step > 1), aes(x = step, y = log(V8), group = run, col = factor(config))) + geom_line() + geom_point()
ggplot(subset(a, step > 1), aes(x = step, y = V11, group = run, col = factor(config))) + geom_line() + geom_point()
filter(a, step == 4) %>% group_by(config) %>% summarize(se = sd(V8), m = mean(V8), min = min(V8), max = max(V8))
