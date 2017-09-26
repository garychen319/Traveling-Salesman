import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;


public class Graph {

  // Keep a fast index to nodes in the map
  private Map<Integer, Vertex> vertexNames;

  /**
   * Construct an empty Graph with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Graph() {
    vertexNames = new HashMap<>();
  }

  /**
   * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the graph
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the graph
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the graph
   */
  public Collection<Vertex> getVertices() {
    return vertexNames.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(String name) {
    return vertexNames.get(name);
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(int nameU, int nameV, Double cost) {
    if (!vertexNames.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertexNames.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertexNames.get(nameU);
    Vertex targetVertex = vertexNames.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param name
   *          (String) name of vertex u
   * @param name2
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(int name, int name2, double cost) {
    addEdge(name, name2, cost);
    addEdge(name2, name, cost);
  }


  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
    return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
  }

  /**
   * Computes euclidean distance between two vertices as described by their
   * coordinates
   * 
   * @param u
   *          (Vertex) vertex u
   * @param v
   *          (Vertex) vertex v
   * @return (double) distance between two vertices
   */
  public double computeEuclideanDistance(Vertex u, Vertex v) {
    return computeEuclideanDistance(u.x, u.y, v.x, v.y);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
    for (Vertex u : getVertices())
      for (Edge uv : u.adjacentEdges) {
        Vertex v = uv.target;
        uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
      }
  }



  // STUDENT CODE STARTS HERE

  public void generateRandomVertices(int n) {
    vertexNames = new HashMap<>(); // reset the vertex hashmap
    
    Random rand = new Random();
    
    for (int i=0; i<n; i++){
    	Vertex temp = new Vertex(i, rand.nextInt(101), rand.nextInt(101)); //make vertex at these random locations, 1-100 inclusive
    	addVertex(temp);
    }
    //go through all possible combinations of vertex pairs
    for(Vertex node1 : vertexNames.values()){
    	for(Vertex node2 : vertexNames.values()){
    		if(node1 != node2){ //if they are not the same vertex
    			double edgeDist = computeEuclideanDistance(node1, node2); //calculate distance between them
    			addUndirectedEdge(node1.name, node2.name, edgeDist); //add this edge
    		}
    	}
    }
    computeAllEuclideanDistances(); // compute distances
  }
  
  public List<Edge> nearestNeighborTsp() {
	List<Edge> traveledEdges = new ArrayList<>(); //The list of edges that is traversed, returned by function
	List<Vertex> unvisitedNodes = new ArrayList<>(); //The list of vertices that hasn't been traveled to yet.
    for(Vertex node:vertexNames.values()){ //copies values in the vertexNames hash table into unvisitedNodes
    	unvisitedNodes.add(node);
    }
    
    Vertex start = unvisitedNodes.get(0); //start at the 0th node of unvisited list
    unvisitedNodes.remove(0); //remove it from list (it's now visited)
    Vertex pivot = start; //pivot = node algorithm is currently searching
    
    while(!unvisitedNodes.isEmpty()){ //while there are still nodes to visit
      Edge cheapestEdge = null;
      double cheapestVal = 999; //larger than 100 root 2
      Vertex next = null;
      
      if(pivot.adjacentEdges.size() != 0){
	      for(Edge e:pivot.adjacentEdges){ //search through all adjacent edges
	    	if(e.distance < cheapestVal && unvisitedNodes.contains(e.target)){ //look for the lowest cost edge that links to an unvisited target
	    		cheapestEdge = e; //remember the edge
	    		cheapestVal = e.distance; //mark as current lowest distance
	            next = e.target; //remember the target of the edge
	        }
	      }
      }
      pivot = next; //move pivot along to e.target
      unvisitedNodes.remove(pivot); //mark is as visited
      traveledEdges.add(cheapestEdge);//add it to the list of edges traversed
    }
 
	  //Add final edge to complete cycle
	  double dist = computeEuclideanDistance(pivot,start); //distance from final node to start node
	  Edge finalEdge = new Edge(pivot, start, dist); //Create an edge between final node and start node
	  traveledEdges.add(finalEdge); //Add this edge to the traveled path (to complete the cycle)
	  return traveledEdges;
  }

  
  
  public List<Edge> bruteForceTsp() {
	  List<Edge> traveledEdges = new ArrayList<>(); //The list of edges that is traversed, returned by function
	    List<Vertex> unvisitedNodes = new ArrayList<>(); //The list of vertices that hasn't been traveled to yet.
	    for(Vertex node:vertexNames.values()){ //copies values in the vertexNames hash table into unvisitedNodes
	        unvisitedNodes.add(node);
	      }
	    List<List<Vertex>> permutationsList = new ArrayList<>();
	    List<Vertex> emptyList = new ArrayList<>();
	    computePermutations(permutationsList, unvisitedNodes, emptyList);
	    
	    double finalDist = 9999999; //initialize to something bigger than the path distance
	    for(List<Vertex> permutation : permutationsList){
	      
		  //start by adding the final edge from end to start, the edge to complete the cycle
		  Vertex start = permutation.get(0);
		  Vertex end = permutation.get(permutation.size()-1);
		  double dist = computeEuclideanDistance(end,start); //distance from final node to start node
		  Edge finalEdge = new Edge(end, start, dist); //Create an edge between final node and start node
		  
		  Double tempDist = dist; //stores temporary distance, first add distance of final edge
		  List<Edge> tempEdges = new LinkedList<Edge>(); //stores temporary path
		  tempEdges.add(finalEdge); //add the final edge to this list first
		
		  //iterate through all vertices in the permutation and add their edges to tempEdges, distances to tempDist
		  for(int i = 0; i < permutation.size() - 1; i++){
	    	Vertex v = permutation.get(i);
	    	Vertex x = permutation.get(i+1);
	    	
	    	double d = computeEuclideanDistance(v,x);
	    	tempDist += d; //add the distance of the edge to the distance list
	    	
	    	Edge e = new Edge(v, x, d);
	    	tempEdges.add(e); //add the edge to the temporary edge list
	      }

	      if(tempDist < finalDist){ //if found a path smaller than previous smallest path
	    	traveledEdges = tempEdges; //replace old edges list
	    	finalDist = tempDist; //replace old smallest distance
	      }
	    }
	    return traveledEdges;
	  }

	  private void computePermutations(List<List<Vertex>> permList, List<Vertex> unvisited, List<Vertex> perm){
		/*permList is the overall list of permutations, unvisited is unvisited nodes which decreases throughout recursion
		perm is permutation sequence that increases throughout recursion*/
		  
	    if(unvisited.size() == 0){ //all nodes visited, base case, done.
	    	permList.add(perm); //add permutation calculated in this cycle to big permutation list of lists.
	    }

      for(int i = 0; i < unvisited.size(); i++){ //loop through each of the remaining vertices and permute
    	//initialize 2 new lists to avoid overwriting, need the original for other recursive calls
        List<Vertex> perm1 = new ArrayList<>(perm);
        List<Vertex> unvisited1 = new ArrayList<>(unvisited);
        
        Vertex v = unvisited.get(i); //take another node from unvisited
        perm1.add(v); //add a vertex to the permutation
        unvisited1.remove(v); //remove from unvisited
        computePermutations(permList, unvisited1, perm1); //recurse with new unvisited and perm (unvisited size--, perm size ++)
	    }
	  }
  
  // STUDENT CODE ENDS HERE



  /**
   * Prints out the adjacency list of the graph for debugging
   */
  public void printAdjacencyList() {
    for (int u : vertexNames.keySet()) {
      StringBuilder sb = new StringBuilder();
      sb.append(u);
      sb.append(" -> [ ");
      for (Edge e : vertexNames.get(u).adjacentEdges) {
        sb.append(e.target.name);
        sb.append("(");
        sb.append(e.distance);
        sb.append(") ");
      }
      sb.append("]");
      System.out.println(sb.toString());
    }
  }
}
