# Pathfinding
An attempt at creating a new pathfinding algorithm, which would essentially perform a two-way meet in the middle A* search, except that instead of aiming for the source or destination node, each frontier would aim for the closest node on the opposite frontier. Inspired by the art of creating [lichtenberg figures in wood](http://imgur.com/5IJ9VJo).
# Outcomes
The algorithm ended up being much slower than A* (not surprisingly), sometimes finds a less than optimal path, and encounters program issues if a path is not available. 
