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


###
load("~/ma/kyoto/aracne/data/airway_cancer_symbol_FS2_ARACNE.rda")
load("~/ma/kyoto/aracne/data/airway_normal_symbol_FS2_ARACNE.rda")
load("~/ma/kyoto/aracne/data/allegro_cancer_symbol_FS2_ARACNE.rda")
load("~/ma/kyoto/aracne/data/allegro_normal_symbol_FS2_ARACNE.rda")
load("~/ma/kyoto/aracne/data/airway_symbol_FS2_ARACNE.rda")
load("~/ma/kyoto/aracne/data/allegro_symbol_FS2_ARACNE.rda")
load("~/ma/kyoto/aracne/data/airway_symbol_ARACNE_noNULL.rda")
load("~/ma/kyoto/aracne/data/allegro_symbol_ARACNE_noNULL.rda")


air.gc = new ('graphAM', adjMat=air.cancer.aracne.fs2.network, edgemode='undirected', )
air.gn = new ('graphAM', adjMat=air.normal.aracne.fs2.network, edgemode='undirected')
all.gc = new ('graphAM', adjMat=all.cancer.aracne.fs2.network, edgemode='undirected')
all.gn = new ('graphAM', adjMat=all.normal.aracne.fs2.network, edgemode='undirected')
air.gffs = new ('graphAM', adjMat=air.aracne.fs2.network, edgemode='undirected')
all.gffs = new ('graphAM', adjMat=all.aracne.fs2.network, edgemode='undirected')
air.gf = new ('graphAM', adjMat=air.network, edgemode='undirected')
all.gf = new ('graphAM', adjMat=all.network, edgemode='undirected')

if (!exists ('cy'))   # get access to application-level Cytoscape operations.
  cy = CytoscapeConnection ()

window.name1 = 'Airway ARACNE cancer network'
net1 <- air.gc

window.name2 = 'Airway ARACNE normal network'
net2 <- air.gn

window.name3 = 'Airway ARACNE full network'
net3 <- air.gffs

window.name4 = 'Airway ARACNE original network'
net4 <- air.gf

window.name5 = 'Allegro ARACNE cancer network'
net5 <- all.gc

window.name6 = 'Allegro ARACNE normal network'
net6 <- all.gn

window.name7 = 'Allegro ARACNE full network'
net7 <- all.gffs

window.name8 = 'Allegro ARACNE original network'
net8 <- all.gf

window.name = window.name4
net <- net4

if (window.name %in% as.character (getWindowList (cy)))
  deleteWindow (cy, window.name)

cw = CytoscapeWindow (window.name, net)
displayGraph (cw)
redraw (cw)
layout (cw, 'jgraph-spring')