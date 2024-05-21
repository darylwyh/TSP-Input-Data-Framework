/**
 * An implementation of the traveling salesman problem in Java using dynamic programming to improve
 * the time complexity from O(n!) to O(n^2 * 2^n).
 *
 * <p>Time Complexity: O(n^2 * 2^n) Space Complexity: O(n * 2^n)
 *
 * @author William Fiset, william.alexandre.fiset@gmail.com
 */
// package com.williamfiset.algorithms.graphtheory;
package williamfiset; 

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class TspDynamicProgrammingIterative {

  private final int N, start;
  private final double[][] distance; // TODO: distance matrix for ISP data 
  private List<Integer> tour = new ArrayList<>(); // solution tour 
  private double minTourCost = Double.POSITIVE_INFINITY; 
  private boolean ranSolver = false;

  public TspDynamicProgrammingIterative(double[][] distance) {
    this(0, distance); // start with 0 as start node 
  }

  public TspDynamicProgrammingIterative(int start, double[][] distance) {
    N = distance.length; // row length => distance[0].length 

    if (N <= 2) throw new IllegalStateException("N <= 2 not yet supported.");
    if (N != distance[0].length) throw new IllegalStateException("Matrix must be square (n x n)");
    if (start < 0 || start >= N) throw new IllegalArgumentException("Invalid start node.");
    if (N > 32)
      throw new IllegalArgumentException(
          "Matrix too large! A matrix that size for the DP TSP problem with a time complexity of"
              + "O(n^2*2^n) requires way too much computation for any modern home computer to handle");

    this.start = start; // init start node 
    this.distance = distance; // init matrix 
  }

  // Returns the optimal tour for the traveling salesman problem.
  public List<Integer> getTour() {
    if (!ranSolver) solve();
    return tour;
  }

  // Returns the minimal tour cost.
  public double getTourCost() {
    if (!ranSolver) solve();
    return minTourCost;
  }

  // Solves the traveling salesman problem and caches solution.
  public void solve() {

    if (ranSolver) return;
    // exp: N = 6, END_STATE = 00111111, all cities visited 
    final int END_STATE = (1 << N) - 1; 
    Double[][] memo = new Double[N][1 << N];
 
    // Add all outgoing edges from the starting node to memo table.
    // Stores all distances in memoization table. As base case. 
    for (int end = 0; end < N; end++) {
      if (end == start) continue; // don't add an edge from a node to itself
      // both the starting and the destination nodes visited
      memo[end][(1 << start) | (1 << end)] = distance[start][end];
    }

    // minimum number of cities required for such tour cycle is 3
    for (int r = 3; r <= N; r++) {
      // iterates over all combinations of cities with r cities visited out of N
      // subset represents each combination of cities
      for (int subset : combinations(r, N)) {
        if (notIn(start, subset)) continue; // ensure starting city in subset
        // iterates over all possible cities to move to from the starting city
        for (int next = 0; next < N; next++) {
          if (next == start || notIn(next, subset)) continue;
          int subsetWithoutNext = subset ^ (1 << next);
          double minDist = Double.POSITIVE_INFINITY;
          for (int end = 0; end < N; end++) {
            if (end == start || end == next || notIn(end, subset)) continue;
            double newDistance = memo[end][subsetWithoutNext] + distance[end][next];
            if (newDistance < minDist) {
              minDist = newDistance;
            }
          }
          memo[next][subset] = minDist;
        }
      }
    }

    // Connect tour back to starting node and minimize cost.
    for (int i = 0; i < N; i++) {
      if (i == start) continue;
      // memo[..]: cost of completing tour from city i to starting city
      double tourCost = memo[i][END_STATE] + distance[i][start];
      if (tourCost < minTourCost) {
        minTourCost = tourCost;
      }
    }

    int lastIndex = start;
    int state = END_STATE;
    tour.add(start);

    // Reconstruct TSP path from memo table.
    for (int i = 1; i < N; i++) {

      int bestIndex = -1; // index of best city to visit next
      double bestDist = Double.POSITIVE_INFINITY;
      for (int j = 0; j < N; j++) {
        if (j == start || notIn(j, state)) continue;
        double newDist = memo[j][state] + distance[j][lastIndex];
        if (newDist < bestDist) {
          bestIndex = j;
          bestDist = newDist;
        }
      }

      tour.add(bestIndex);
      // mark the chosen city as visited by removing it from current state
      state = state ^ (1 << bestIndex); 
      lastIndex = bestIndex;
    }

    tour.add(start);
    Collections.reverse(tour);

    ranSolver = true; // bool to tell tsp is solved 
  }

  private static boolean notIn(int elem, int subset) {
    return ((1 << elem) & subset) == 0;
  }

  // This method generates all bit sets of size n where r bits
  // are set to one. The result is returned as a list of integer masks.
  public static List<Integer> combinations(int r, int n) {
    List<Integer> subsets = new ArrayList<>();
    combinationsHelper(0, 0, r, n, subsets);
    return subsets;
  }

  // To find all the combinations of size r we need to recurse until we have
  // selected r elements (aka r = 0), otherwise if r != 0 then we still need to select
  // an element which is found after the position of our last selected element
  private static void combinationsHelper(int set, int at, int r, int n, List<Integer> subsets) {

    // Return early if there are more elements left to select than what is available.
    int elementsLeftToPick = n - at;
    if (elementsLeftToPick < r) return;

    // We selected 'r' elements so we found a valid subset!
    if (r == 0) {
      subsets.add(set);
    } else {
      for (int i = at; i < n; i++) {
        // Try including this element
        set ^= (1 << i);

        combinationsHelper(set, i + 1, r - 1, n, subsets);

        // Backtrack and try the instance where we did not include this element
        set ^= (1 << i);
      }
    }
  }

  public static void main(String[] args) {
    // Create adjacency matrix
    int n = 6;
    double[][] distanceMatrix = new double[n][n];
    for (double[] row : distanceMatrix) java.util.Arrays.fill(row, 10000);
    distanceMatrix[5][0] = 10;
    distanceMatrix[1][5] = 12;
    distanceMatrix[4][1] = 2;
    distanceMatrix[2][4] = 4;
    distanceMatrix[3][2] = 6;
    distanceMatrix[0][3] = 8;

    int startNode = 0;
    TspDynamicProgrammingIterative solver =
        new TspDynamicProgrammingIterative(startNode, distanceMatrix);

    // Prints: [0, 3, 2, 4, 1, 5, 0]
    System.out.println("Tour: " + solver.getTour());

    // Print: 42.0
    System.out.println("Tour cost: " + solver.getTourCost());
  }
}