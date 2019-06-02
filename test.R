library(ggplot2)

# Set up the vectors
true_gambling <- c("tSafe", "tGamble")
2rediction <- c("pSafe","pGamble")

# Create the data frame
df <- expand.grid(true_gambling, prediction)
        #true: S, G  pred:
df$value <- c(10,20, # S
              30,40) # G

#Plot the Data
g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = c("green", "red", "red", "green")) +
    theme(legend.position="none") + xlab("") + ylab("") + ggtitle("subjet 1")
g + scale_size_continuous(range=c(10,30)) + geom_text(aes(label = value))

