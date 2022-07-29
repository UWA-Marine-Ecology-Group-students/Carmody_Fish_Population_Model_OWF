library(sfnetworks)
library(tidyr)
library(dplyr)
library(sf)
library(TSP)
library(igraph)


#### CREATE A SF POINT USING THE CENTROIDS OF THE POLYGONS AND WORK OUT THE CLOSEST POINTS TO EACH ONE ####
points_sf <- st_as_sf(points, coords = c("X", "Y"), agr="identity", crs=4283)
points_sp <- st_cast(points_sf$geometry, "POINT", crs=4283) 

## Figure out which six points are the closest - won't all have six close ones but can start from there
dist.mat <- st_distance(points_sp) # Great Circle distance since in lat/lon
# Number within 1.5km: Subtract 1 to exclude the point itself
# num.5km <- apply(dist.mat, 1, function(x) {
#   sum(x < 0.05) - 1
# })

# nn.dist <- apply(dist.mat, 1, function(x) {
#   return(sort(x, partial = 2)[2])
# })

closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[2] })
second.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[3] })
third.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[4] })
fourth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[5] })
fifth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[6] })
sixth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[7] })
seventh.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[8] })
eighth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[9] })
ninth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[10] })
tenth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[11] })

neighbours <- points_sf %>% 
  dplyr::select(ID) %>% 
  mutate(`closest` = closest) %>% 
  mutate(second = second.closest) %>% 
  mutate(third = third.closest) %>% 
  mutate(fourth = fourth.closest) %>% 
  mutate(fifth = fifth.closest) %>% 
  mutate(sixth = sixth.closest) %>% 
  mutate(seventh = seventh.closest) %>% 
  mutate(eighth = eighth.closest) %>% 
  mutate(ninth = ninth.closest) %>% 
  mutate(tenth = tenth.closest)

closest <- as.data.frame(closest) %>%
  rename(ID = "closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y")) 
closest <- st_cast(st_geometry(closest), "POINT") 

second.closest <- as.data.frame(second.closest) %>%
  rename(ID = "second.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
second.closest <- st_cast(st_geometry(second.closest), "POINT")

third.closest <- as.data.frame(third.closest) %>%
  rename(ID = "third.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
third.closest <- st_cast(st_geometry(third.closest), "POINT")

fourth.closest <- as.data.frame(fourth.closest) %>%
  rename(ID = "fourth.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
fourth.closest <- st_cast(st_geometry(fourth.closest), "POINT")

fifth.closest <- as.data.frame(fifth.closest) %>%
  rename(ID = "fifth.closest") %>% 
  inner_join(., points, by="ID") %>%  
  st_as_sf(., coords = c("X", "Y"))
fifth.closest <- st_cast(st_geometry(fifth.closest), "POINT")

sixth.closest <- as.data.frame(sixth.closest) %>%
  rename(ID = "sixth.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
sixth.closest <- st_cast(st_geometry(sixth.closest), "POINT")

seventh.closest <- as.data.frame(seventh.closest) %>%
  rename(ID = "seventh.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
seventh.closest <- st_cast(st_geometry(seventh.closest), "POINT")

eighth.closest <- as.data.frame(eighth.closest) %>%
  rename(ID = "eighth.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
eighth.closest <- st_cast(st_geometry(eighth.closest), "POINT")

ninth.closest <- as.data.frame(ninth.closest) %>%
  rename(ID = "ninth.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
ninth.closest <- st_cast(st_geometry(ninth.closest), "POINT")

tenth.closest <- as.data.frame(tenth.closest) %>%
  rename(ID = "tenth.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
tenth.closest <- st_cast(st_geometry(tenth.closest), "POINT")

#### CONNECT THE POINTS TO THEIR NEIGHBOURS TO FORM A NETWORK ####
n <- nrow(points)

# Closest
linestrings.closest <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.closest <- st_multilinestring(do.call("rbind", linestrings.closest))

# Second closest
linestrings.second <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], second.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.second <- st_multilinestring(do.call("rbind", linestrings.second))

# Third closest
linestrings.third <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], third.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.third<- st_multilinestring(do.call("rbind", linestrings.third))

# Fourth closest
linestrings.fourth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], fourth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.fourth <- st_multilinestring(do.call("rbind", linestrings.fourth))

# Fifth closest
linestrings.fifth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], fifth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.fifth <- st_multilinestring(do.call("rbind", linestrings.fifth))

# Sixth closest
linestrings.sixth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], sixth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.sixth <- st_multilinestring(do.call("rbind", linestrings.sixth)) # You have to turn this into a multilinestring when you put rbind it 

# Seventh closest
linestrings.seventh <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], seventh.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.seventh <- st_multilinestring(do.call("rbind", linestrings.seventh)) # You have to turn this into a multilinestring when you put rbind it 

# Eighth closest
linestrings.eighth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], eighth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.eighth <- st_multilinestring(do.call("rbind", linestrings.eighth)) # You have to turn this into a multilinestring when you put rbind it 

# Ninth closest
linestrings.ninth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], ninth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.ninth <- st_multilinestring(do.call("rbind", linestrings.ninth)) # You have to turn this into a multilinestring when you put rbind it 

# tenth closest
linestrings.tenth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], tenth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.tenth <- st_multilinestring(do.call("rbind", linestrings.tenth)) # You have to turn this into a multilinestring when you put rbind it 

# Linestring neighbours
linestrings.neighbours <- lapply(X = 1:1538, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], points_sp[x+1]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.neighbours <- st_multilinestring(do.call("rbind", linestrings.neighbours))
neighbour.lines <- st_cast(st_sfc(multilinestring.neighbours),"LINESTRING")
neighbour.lines <- st_as_sf(neighbour.lines) %>% 
  mutate(length = st_length(x)*111)
lines.to.remove <- which(neighbour.lines$length>20)
neighbour.lines <- neighbour.lines[-c(lines.to.remove), ]
neighbour.lines <- st_cast(st_geometry(neighbour.lines), "MULTILINESTRING")
st_crs(neighbour.lines) <- 4283

## Join together
close.lines <- st_union(multilinestring.closest, neighbour.lines)
                  
lines <- st_union(multilinestring.closest, multilinestring.second, multilinestring.third, 
                  multilinestring.fourth, multilinestring.fifth, multilinestring.sixth)
lines <- st_cast(st_geometry(lines), "LINESTRING")
lines <- st_as_sf(lines)
neighbour.lines <- st_cast(st_geometry(neighbour.lines), "LINESTRING")
lines <- rbind(neighbour.lines, lines)

lines <- st_cast(st_geometry(lines), "LINESTRING")

st_crs(lines) <- 4283
st_crs(points_sp) <- 4283
st_crs(points_sf) <- 4283

#### SET UP THE SF NETWORK AND CREATE A DISTANCE MATRIX ####
network <- as_sfnetwork(lines, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

## Add our point IDs to the nodes of the network so we can check everything is working right
nodes <- activate(network, "nodes")

# network <- st_join(nodes, points_sf, join = st_intersects)
# 
# network <- activate(network, "nodes") %>% 
#   arrange(., ID) %>% 
#   dplyr::distinct(., .keep_all = TRUE) %>%  # Remove any duplicate IDs where nodes
#   filter(!is.na(ID))

## Calculate the distances from each point to every other point on the network
net <- activate(network, "nodes")
cost_matrix <- (st_network_cost(net, from=points_sf))/1000

dim(cost_matrix) # Check that the dimensions match up to how many points you think you should have in the network

test <- st_network_paths(network, from=1, to=1514, algorithm=c("bellman-ford"))

ids <- as.data.frame(st_nearest_feature(net)) # These are the IDs of the nearest feature to every point
row.names(cost_matrix) <- ids # Set our row names to the point IDs
colnames(cost_matrix) <- ids # Set our column names to the point IDS

names <- as.data.frame(colnames(cost_matrix))

test <- as.data.frame(ids)

# Shit that doesn't work
## Check to see if any points are missing and if so add them to the network
nodes <- st_as_sf(network, "nodes") %>% 
  st_drop_geometry() %>% 
  unlist() %>% 
  unname()
ids <- points$ID
missing <- setdiff(ids, nodes)
add.nodes <- neighbours[as.numeric(missing), ]
add.nodes <- add.nodes[,1:2]

# Format nodes correctly if any are missing
add.nodes <- st_cast(st_geometry(add.nodes), "POINT")
add.nodes <- st_as_sf(st_sfc(add.nodes, crs = 4283))
add.nodes <- st_intersection(add.nodes, points_sp)
add.nodes <- add.nodes %>% 
  filter(ID %in% c(missing))

# Add new nodes
network <- st_network_blend(network, add.nodes) # There should be more nodes than before, there might now be duplicates as well because of the way the network is structured

edges <- st_as_sf(network, "edges") %>% 
  st_drop_geometry() %>% 
  unlist() %>% 
  unname()
