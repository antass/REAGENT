lnames <- load("~/Meta_Analysis/kyoto/aracne/data/airway_ARACNE.rda")
air.network[1:10, 1:10]

# Number of connections
sum(air.network>0)

# Number of nodes in network
nrow(air.network[apply(air.network, 1, sum)>0, ])

# Average strength of connections
summary(as.numeric(air.network))

library(Rgraphviz, lib.loc="/home/aniat/R/x86_64-unknown-linux-gnu-library/")
svg("~/Meta_Analysis/kyoto/aracne/Airway_aracne_full_network.svg")
plot( as( air.network ,"graphNEL") )
dev.off()


