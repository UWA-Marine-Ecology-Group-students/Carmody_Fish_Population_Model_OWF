library(sfnetworks)
library(dplyr)
library(sf)
library(TSP)


#### CREATE A SF POINT USING THE CENTROIDS OF THE POLYGONS AND WORK OUT THE CLOSEST POINTS TO EACH ONE ####
points$ID <- as.integer(points$ID)
points2 <- st_as_sf(points, coords = c("X", "Y"))
points2 <- st_cast(st_geometry(points2), "POINT")


## Figure out which six points are the closest - won't all have six close ones but can start from there
dist.mat <- st_distance(points2) # Great Circle distance since in lat/lon
# Number within 1.5km: Subtract 1 to exclude the point itself
num.5km <- apply(dist.mat, 1, function(x) {
  sum(x < 0.05) - 1
})
rm(num.3000)

nn.dist <- apply(dist.mat, 1, function(x) {
  return(sort(x, partial = 2)[2])
})

closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[2] })
second.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[3] })
third.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[4] })
fourth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[5] })
fifth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[6] })
sixth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[7] })

neighbours <- water %>% 
  dplyr::select(ID) %>% 
  mutate(`closest` = closest) %>% 
  mutate(second = second.closest) %>% 
  mutate(third = third.closest) %>% 
  mutate(fourth = fourth.closest) %>% 
  mutate(fifth = fifth.closest) %>% 
  mutate(sixth = sixth.closest)

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

#### CONNECT THE POINTS TO THEIR NEIGHBOURS TO FORM A NETWORK ####
n <- nrow(points)

# Closest
linestrings.closest <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.closest <- st_multilinestring(do.call("rbind", linestrings.closest))

# Second closest
linestrings.second <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], second.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.second <- st_multilinestring(do.call("rbind", linestrings.second))

# Third closest
linestrings.third <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], third.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.third<- st_multilinestring(do.call("rbind", linestrings.third))

# Fourth closest
linestrings.fourth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], fourth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.fourth <- st_multilinestring(do.call("rbind", linestrings.fourth))

# Fifth closest
linestrings.fifth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], fifth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.fifth <- st_multilinestring(do.call("rbind", linestrings.fifth))

# Sixth closest
linestrings.sixth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], sixth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.sixth <- st_multilinestring(do.call("rbind", linestrings.sixth))

connected <- st_combine(c(multilinestring.closest, multilinestring.second, multilinestring.third, multilinestring.fourth, multilinestring.fifth, multilinestring.sixth))
connected <- st_cast(connected, "LINESTRING") # Needs to be a line string rather than multiline for the next step

#### SET UP THE SF NETWORK AND CREATE A DISTANCE MATRIX ####
network <- as_sfnetwork(connected, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

## Calculate the distances from each point to every other point on the network
net <- activate(network, "nodes")
cost_matrix < st_network_cost(net)
dim(cost_matrix) # Check that the dimensions match up to how many points you think you should have in the network









