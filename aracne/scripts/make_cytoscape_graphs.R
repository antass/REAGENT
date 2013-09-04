library (RCytoscape)

# create a vector of 25 random 0's and 1's
values = round (runif (25))

# arrange these into a 5 by 5 matrix
m = matrix (values, nrow=5, ncol=5, byrow=TRUE)

# give the matrix simple integer row and column names.  These become the node names
rownames (m) = 1:5
colnames (m) = 1:5

# if you want to remove self-loops, enable this line
#diag (m) = 0   # remove any self-loops

# now create a Bioconductor graph of subclass 'graphAM'
g = new ('graphAM', adjMat=m, edgemode='directed')

# these next few lines are useful if you run this script repetetively, since RCy requires that all window names are unique,
# these lines allows you to detect and then delete any pre-existing window with the same name

if (!exists ('cy'))   # get access to application-level Cytoscape operations.
  cy = CytoscapeConnection ()

window.name = 'adjacency matrix graph simple demo'

if (window.name %in% as.character (getWindowList (cy)))
  deleteWindow (cy, window.name)

# now create a window, send the graph, render it, lay it out

cw = CytoscapeWindow (window.name, g)
displayGraph (cw)
redraw (cw)
layout (cw, 'jgraph-spring')


### ARACNE
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_ARACNE.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_cancer_FS_ARACNE.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_premalignancy_FS_ARACNE.rda")  # joint.aracne.p.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/allegro_normal_FS_ARACNE.rda")  # all.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/airway_normal_FS_ARACNE.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/greentea_normal_FS_ARACNE.rda")  # joint.aracne.p.fs.networklname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_ARACNE.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/lam_normal_FS_ARACNE.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_FS_ARACNE.rda")  # joint.aracne.p.fs.networklname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_ARACNE.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/airway_cancer_FS_ARACNE.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/greentea_premalignancy_FS_ARACNE.rda")  # joint.aracne.p.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/lam_premalignancy_FS_ARACNE.rda")  # joint.aracne.p.fs.network
### CLR
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_CLR.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_cancer_FS_CLR.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_premalignancy_FS_CLR.rda")  # joint.aracne.p.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/allegro_normal_FS_CLR.rda")  # all.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/airway_normal_FS_CLR.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/greentea_normal_FS_CLR.rda")  # joint.aracne.p.fs.networklname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_CLR.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/lam_normal_FS_CLR.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_FS_CLR.rda")  # joint.aracne.p.fs.networklname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_CLR.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/airway_cancer_FS_CLR.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/greentea_premalignancy_FS_CLR.rda")  # joint.aracne.p.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/lam_premalignancy_FS_CLR.rda")  # joint.aracne.p.fs.network

### ARACNE
c.net = new ('graphAM', adjMat=joint.aracne.c.fs.network, edgemode='undirected')
n.net = new ('graphAM', adjMat=joint.aracne.n.fs.network, edgemode='undirected')
p.net = new ('graphAM', adjMat=joint.aracne.p.fs.network, edgemode='undirected')
air.c.net = new ('graphAM', adjMat=air.aracne.c.fs.network, edgemode='undirected')
air.n.net = new ('graphAM', adjMat=air.aracne.n.fs.network, edgemode='undirected')
grt.p.net = new ('graphAM', adjMat=grt.aracne.p.fs.network, edgemode='undirected')
all.c.net = new ('graphAM', adjMat=all.aracne.c.fs.network, edgemode='undirected')
all.n.net = new ('graphAM', adjMat=all.aracne.n.fs.network, edgemode='undirected')
lam.p.net = new ('graphAM', adjMat=lam.aracne.p.fs.network, edgemode='undirected')
grt.n.net = new ('graphAM', adjMat=grt.aracne.n.fs.network, edgemode='undirected')
lam.n.net = new ('graphAM', adjMat=lam.aracne.n.fs.network, edgemode='undirected')
### CLR
c.net.clr = new ('graphAM', adjMat=joint.clr.c.fs.network, edgemode='undirected')
n.net.clr = new ('graphAM', adjMat=joint.clr.n.fs.network, edgemode='undirected')
p.net.clr = new ('graphAM', adjMat=joint.clr.p.fs.network, edgemode='undirected')
air.c.net.clr = new ('graphAM', adjMat=air.clr.c.fs.network, edgemode='undirected')
air.n.net.clr = new ('graphAM', adjMat=air.clr.n.fs.network, edgemode='undirected')
grt.p.net.clr = new ('graphAM', adjMat=grt.clr.p.fs.network, edgemode='undirected')
all.c.net.clr = new ('graphAM', adjMat=all.clr.c.fs.network, edgemode='undirected')
all.n.net.clr = new ('graphAM', adjMat=all.clr.n.fs.network, edgemode='undirected')
lam.p.net.clr = new ('graphAM', adjMat=lam.clr.p.fs.network, edgemode='undirected')
grt.n.net.clr = new ('graphAM', adjMat=grt.clr.n.fs.network, edgemode='undirected')
lam.n.net.clr = new ('graphAM', adjMat=lam.clr.n.fs.network, edgemode='undirected')
 

if (!exists ('cy'))   # get access to application-level Cytoscape operations.
  cy = CytoscapeConnection ()

### ARACNE
window.name1 = 'N'; net1 <- n.net #'Normal ARACNE network'
window.name2 = 'C'; net2 <- c.net #Cancer ARACNE network'
window.name3 = 'P'; net3 <- p.net #'Premalignancy ARACNE network'
window.name4 = 'N air'; net4 <- air.n.net#Normal Airway ARACNE network'
window.name5 = 'C air'; net5 <- air.c.net#Cancer Airway ARACNE network'
window.name6 = 'P grt'; net6 <- grt.p.net#Premalignancy Green Tea ARACNE network'
window.name7 = 'N all'; net7 <- all.n.net#Normal Airway ARACNE network'
window.name8 = 'C all'; net8 <- all.c.net#Cancer Airway ARACNE network'
window.name9 = 'P lam'; net9 <- lam.p.net#Premalignancy Green Tea ARACNE network'
window.name10 = 'N grt'; net10 <- grt.n.net#Normal Airway ARACNE network'
window.name11 = 'N lam'; net11 <- lam.n.net#Normal Airway ARACNE network'
### CLR
window.name12 = 'N CLR'; net.clr12 <- n.net.clr #'Normal ARACNE net.clrwork'
window.name13= 'C CLR'; net.clr13 <- c.net.clr #Cancer ARACNE net.clrwork'
window.name14 = 'P CLR'; net.clr14 <- p.net.clr #'Premalignancy ARACNE net.clrwork'
window.name15 = 'N air CLR'; net.clr15 <- air.n.net.clr#Normal Airway ARACNE net.clrwork'
window.name16 = 'C air CLR'; net.clr16 <- air.c.net.clr#Cancer Airway ARACNE net.clrwork'
window.name17 = 'P grt CLR'; net.clr17 <- grt.p.net.clr#Premalignancy Green Tea ARACNE net.clrwork'
window.name18 = 'N all CLR'; net.clr18 <- all.n.net.clr#Normal Airway ARACNE net.clrwork'
window.name19 = 'C all CLR'; net.clr19 <- all.c.net.clr#Cancer Airway ARACNE net.clrwork'
window.name20 = 'P lam CLR'; net.clr20 <- lam.p.net.clr#Premalignancy Green Tea ARACNE net.clrwork'
window.name21 = 'N grt CLR'; net.clr21 <- grt.n.net.clr#Normal Airway ARACNE net.clrwork'
window.name22 = 'N lam CLR'; net.clr22 <- lam.n.net.clr#Normal Airway ARACNE net.clrwork'

### 




# specify which network to plot
# id <- 3
for (id in c(1)) {
  window.name = get(paste0("window.name", id))
  net <- get(paste0("net", id))
  
  # Activate RCytoscape RPC plugin
  if (window.name %in% as.character (getWindowList (cy)))
    deleteWindow (cy, window.name)
  
  cw = CytoscapeWindow (window.name, net)
  displayGraph (cw)
  redraw (cw)
#   layout (cw, 'jgraph-spring')
}

for (id in 12) {
  window.name = get(paste0("window.name", id))
  net <- get(paste0("net.clr", id))
  
  # Activate RCytoscape RPC plugin
  if (window.name %in% as.character (getWindowList (cy)))
    deleteWindow (cy, window.name)
  
  cw = CytoscapeWindow (window.name, net)
  displayGraph (cw)
  redraw (cw)
  #   layout (cw, 'jgraph-spring')
}