
#Load the data and library
setwd("C:/Users/SAGAR/OneDrive/Desktop/Sagar/NCSU/MEM/DSC 595")

load("C:/Users/SAGAR/OneDrive/Desktop/Sagar/NCSU/MEM/DSC 595/wiring.rda")

library(igraph)

#Details of the rdpos graph object
rdpos

#Plot the graph 
par(mar=c(2,2,2,2))
plot(rdpos,
     main = "ROETHLISBERGER & DICKSON BANK WIRING ROOM VISUALIZATION")

#The names of the vertices (nodes) in the rdpos
V(rdpos)$name

#Define Fruchterman-Reingold layout
fr_layout <- layout_with_fr(rdpos)

# Define basic colors for the three groups
colors <- setNames(c("red", "green", "blue"), c("I", "S", "W"))

# Assign colors to nodes based on their type (I, S, W)
node_colors <- sapply(V(rdpos)$name, function(name) colors[substr(name, 1, 1)])

#Define labels as the name of the nodes
labels = V(rdpos)$name

#Define edge curvature
curved <- 0.2

# Set plotting margin
par(mar = c(2,2,2,2))

# Plot the graph with node labels using defined parameters
plot(rdpos,
     layout = fr_layout,
     vertex.label = labels,
     vertex.color = node_colors,
     edge.curved = curved,
     main = "ROETHLISBERGER & DICKSON BANK WIRING ROOM VISUALIZATION"
)

# Add a legend
legend("bottomright",
       legend = c("Inspector", "Solderer", "Wireman/Assembler"),
       col = c("red", "green", "blue"),
       pt.cex = 1.5,
       pch = 21, 
       bty = "n",
       title = "Categories")

#Interactive Graph for ROETHLISBERGER & DICKSON BANK WIRING ROOM Network:

library(visNetwork)

V(rdpos)$name

# Create a nodes data frame
nodes <- data.frame(id = 1:vcount(rdpos), label = V(rdpos)$name)

# Create an edges data frame
# Extract edges from the igraph object and map names to IDs
edges <- get.data.frame(rdpos, what = "edges")
edges$from <- match(edges$from, nodes$label)
edges$to <- match(edges$to, nodes$label)

edges

# Create the network visualization
visNetwork(nodes, edges) %>%
  visNodes(color = "PURPLE", shape= "diamond") %>%
  visEdges(smooth = TRUE) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

##Compute centrality scores:
#1.Degree centrality
deg <- degree(rdpos)
deg
#2. Closeness centrality
close <- closeness(rdpos)
close
#3. Betweenness centrality
bet <- betweenness(rdpos)
bet
#4. Eigenvector centrality
eigens_vec<- eigen_centrality(rdpos, scale = FALSE)
eigens_vec$vector

# Centralization with isolates
centr_eigen(rdpos, directed=FALSE)$centralization

# Remove isolates from rdpos
rdpos_without_isolates <- delete.vertices(rdpos, V(rdpos)[degree(rdpos) == 0])

# Centralization without isolates
centr_eigen(rdpos_without_isolates, directed=FALSE)$centralization


# Define the colors for each group
colors <- c("I" = "purple", "S" = "blue", "W" = "yellow")

# Assign colors to nodes based on their type (I, S, W)
node_colors <- ifelse(grepl("^I", V(rdpos_without_isolates)$name), "purple",
                      ifelse(grepl("^S", V(rdpos_without_isolates)$name), "blue", "yellow"))

#Plot without 
plot(rdpos_without_isolates,
     vertex.color = node_colors,
     main = "ROETHLISBERGER & DICKSON BANK WIRING ROOM VISUALIZATION"
)

#Finding Largest Cliques:
largest_cliques(rdpos_without_isolates)

#Data Frame of Largest Cliques:
cliques <- as.data.frame(largest_cliques(rdpos_without_isolates))
cliques

#Computing Coreness:
coreness <- graph.coreness(rdpos_without_isolates)
coreness

#Calculates the modularity:
# Create a numerical membership vector based on the first letter of node names
# Assign '1' for inspectors (I), '2' for solderers (S), and '3' for wiremen/assemblers (W)
membership <- sapply(V(rdpos_without_isolates)$name, function(name) {
  if (substr(name, 1, 1) == "I") {
    return(1)
  } else if (substr(name, 1, 1) == "S") {
    return(2)
  } else if (substr(name, 1, 1) == "W") {
    return(3)
  } else {
    return(NA) # For nodes that do not start with I, S, or W
  }
})

# Ensure membership is a numeric vector
membership <- as.numeric(membership)
membership

mod <- modularity(rdpos_without_isolates, membership)
mod

#Calculate assortativity coefficient for 'name'
assortativity_nominal(
  rdpos_without_isolates, 
  as.integer(as.factor(V(rdpos_without_isolates)$name)))

# Ensure the necessary libraries are loaded
library(igraph)
library(ggraph)
library(ggplot2)

#Performs the Louvain community detection on the rdpos network understand:
louvain_comm <- cluster_louvain(rdpos_without_isolates)
louvain_comm

# Get community memberships from the Louvain community detection result
communities <- membership(louvain_comm)
communities

# Assign these memberships to the vertices of your graph
V(rdpos_without_isolates)$louvain_community <- communities
V(rdpos_without_isolates)$louvain_community

vertex.attributes(rdpos_without_isolates)

# Plot the network
set.seed(123)

ggraph(rdpos_without_isolates, layout = "fr") +
  geom_edge_link(color = "grey", alpha = 0.7) +
  geom_node_point(aes(color = as.factor(louvain_community)), size = 5) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +  # Replace 'name' with the actual attribute for labels
  labs(title = "ROETHLISBERGER & DICKSON BANK WIRING ROOM Network") +
  theme_minimal()


# Load necessary libraries
library(igraph)
library(blockmodeling)

# Set a seed value to ensure reproducibility
set.seed(54321)

# Get adjacency matrix from the igraph object 'rdpos_without_isolates'
mat <- as.matrix(get.adjacency(rdpos_without_isolates))

# Estimate blockmodel partitions with k = 6 using the optRandomParC command
class6 <- optRandomParC(M=mat, k=6, rep=10, approach="ss", blocks="com")
class6

# Retrieve the best partition from the blockmodeling result
best_partition <- class6$best$best1$clu
best_partition

# Assign the block designations to the igraph object 'drug_connect'
V(rdpos_without_isolates)$block <- best_partition
V(rdpos_without_isolates)$block

# Plot the graph with vertices colored by their block designation
par(mar=c(1,1,1,1), mfrow=c(1,1))
plot(rdpos_without_isolates,
     vertex.label=NA,  # Hide vertex labels for a cleaner plot
     vertex.size=5,  # Set vertex size
     edge.arrow.size=.5,  # Set edge arrow size
     vertex.color=V(rdpos_without_isolates)$block,  # Color vertices by block
     main="Rdpos Friendship Graph with Blockmodel Partitions k=6")  # Title for the plot

