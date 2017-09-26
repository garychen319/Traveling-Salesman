Traveling Salesman Problem

Files Included: written.pdf, Graph.java, Edge.java, Vertex.java, Display.java, README.txt

------------------------------------------------------------------------------------------

1.generateRandomVertices
For this function, I made n vertices with random position 1-100 as the x,y values. I then stored the vertices and edges using addVertex and addEdge.

2. nearestNeightborTsp
I kept a list of traveledEdges, which is the return of the function, and this is the list of edges traversed by the nearest neighbor algorithm. Another list unvisitedNodes keeps track of the nodes yet to be seen. At first, this list is an exact copy of vertexNames but decrements as the program goes on. In this function, I have a pivot which marks the current vertex I'm searching, and then a while loop which goes through all other unvisited nodes and finds the cheapest edge from the pivot to an univsited node. When the cheapest edge is found, I add it to traveledEdges and also add the distance to dist (the total distnace). TraveledEdges is returned by the function.

3. bruteForceTsp()
Just like nearestneighbor, I first declare traveledEdges, unvisitedNodes which share the same delcaration as in the previous function. A list of list<vertex> permutationsList stores different lists of the different arrangement of verticies. These arrangements are obtained by calling a helper function below computePermutations. computePermutations is a recursive function which goes through all unvisitedNodes (ends when unvisitedNodes is empty) and creates permutation variations stored in perm, a List<vertex>, and eventually adds perm to permutationsList. It does this by recursion and in each iteration it will decrement the univisited and add another number to perm, until unvisited has nothing left in it.

Once all the permutations are obtained, the rest of bruteForceTsp goes through each list of verticies in permutationsList, calcuates which one has the smallest distance by adding up the values of the edges that conncet the verticies. The permutation that yields the smallest distance is remembered and the combination of edges that makes up the smallest distance is saved to traveledEdges, which is returned by the function.

For more specific details of the function please refer to the line by line comments in the graph.java file.
