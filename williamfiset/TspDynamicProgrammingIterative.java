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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.HashMap;
import java.util.HashSet;

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

    if (N <= 2)
      throw new IllegalStateException("N <= 2 not yet supported.");
    if (N != distance[0].length)
      throw new IllegalStateException("Matrix must be square (n x n)");
    if (start < 0 || start >= N)
      throw new IllegalArgumentException("Invalid start node.");
    if (N > 32)
      throw new IllegalArgumentException(
          "Matrix too large! A matrix that size for the DP TSP problem with a time complexity of"
              + "O(n^2*2^n) requires way too much computation for any modern home computer to handle");

    this.start = start; // init start node
    this.distance = distance; // init matrix
  }

  // Returns the optimal tour for the traveling salesman problem.
  public List<Integer> getTour() {
    if (!ranSolver)
      solve();
    return tour;
  }

  // Returns the minimal tour cost.
  public double getTourCost() {
    if (!ranSolver)
      solve();
    return minTourCost;
  }

  // Solves the traveling salesman problem and caches solution.
  public void solve() {

    if (ranSolver)
      return;
    // exp: N = 6, END_STATE = 00111111, all cities visited
    final int END_STATE = (1 << N) - 1;
    Double[][] memo = new Double[N][1 << N];

    // Add all outgoing edges from the starting node to memo table.
    // Stores all distances in memoization table. As base case.
    for (int end = 0; end < N; end++) {
      if (end == start)
        continue; // don't add an edge from a node to itself
      // both the starting and the destination nodes visited
      memo[end][(1 << start) | (1 << end)] = distance[start][end];
    }

    // minimum number of cities required for such tour cycle is 3
    for (int r = 3; r <= N; r++) {
      // iterates over all combinations of cities with r cities visited out of N
      // subset represents each combination of cities
      for (int subset : combinations(r, N)) {
        if (notIn(start, subset))
          continue; // ensure starting city in subset
        // iterates over all possible cities to move to from the starting city
        for (int next = 0; next < N; next++) {
          if (next == start || notIn(next, subset))
            continue;
          int subsetWithoutNext = subset ^ (1 << next);
          double minDist = Double.POSITIVE_INFINITY;
          for (int end = 0; end < N; end++) {
            if (end == start || end == next || notIn(end, subset))
              continue;
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
      if (i == start)
        continue;
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
        if (j == start || notIn(j, state))
          continue;
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
  // selected r elements (aka r = 0), otherwise if r != 0 then we still need to
  // select
  // an element which is found after the position of our last selected element
  private static void combinationsHelper(int set, int at, int r, int n, List<Integer> subsets) {

    // Return early if there are more elements left to select than what is
    // available.
    int elementsLeftToPick = n - at;
    if (elementsLeftToPick < r)
      return;

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

  // calc Euclidean Distance between two points
  private static double calculateEuclidDistance(double[] point1, double[] point2) {
    double dx = point1[0] - point2[0];
    double dy = point1[1] - point2[1];
    return Math.sqrt(dx * dx + dy * dy);
  }

  // Function to read data from file
  private static double[][] readCoordinatesFromFile(String cityFile, String coordinatesFile) {
    // HashMap<Name of City, Coor of City>
    HashMap<String, double[]> cityCoordinates = new HashMap<>();
    try (BufferedReader br = new BufferedReader(new FileReader(coordinatesFile))) {
      String line;
      while ((line = br.readLine()) != null) {
        String[] parts = line.trim().split(",");

        String cityName = parts[0].trim();
        double latitude = Double.parseDouble(parts[1].trim());
        double longitude = Double.parseDouble(parts[2].trim());
        cityCoordinates.put(cityName, new double[] { latitude, longitude });
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
    // Read from weights.intra file
    try (BufferedReader br = new BufferedReader(new FileReader(cityFile))) {
      String line;
      Set<String> uniqueCoordinates = new HashSet<>(); // Set to store unique coordinates
      List<double[]> uniqueCoordinatesList = new ArrayList<>(); // List to store unique coordinates in order
      while ((line = br.readLine()) != null && uniqueCoordinates.size() < 20) {
        String[] parts = line.trim().split("\\s+");
        String city1 = parts[0].replaceAll("\\+", " ").split(",")[0].trim();
        String city2 = parts[1].replaceAll("\\+", " ").split(",")[0].trim();

        if (cityCoordinates.containsKey(city1) && cityCoordinates.containsKey(city2)) {
          double[] coord1 = cityCoordinates.get(city1);
          double[] coord2 = cityCoordinates.get(city2);
          String coordinateString1 = Arrays.toString(coord1);
          String coordinateString2 = Arrays.toString(coord2);
          if (!uniqueCoordinates.contains(coordinateString1)) {
            uniqueCoordinates.add(coordinateString1);
            uniqueCoordinatesList.add(coord1);
          }
          if (!uniqueCoordinates.contains(coordinateString2)) {
            uniqueCoordinates.add(coordinateString2);
            uniqueCoordinatesList.add(coord2);
          }
        }
      }
      // Convert List to array
      double[][] coordinates = new double[uniqueCoordinatesList.size()][2];
      for (int i = 0; i < uniqueCoordinatesList.size(); i++) {
        coordinates[i] = uniqueCoordinatesList.get(i);
      }
      return coordinates;
    } catch (IOException e) {
      e.printStackTrace();
      return null;
    }
  }

  private static double[][] normalizeCoordinates(double[][] coordinates) {
    // Find minimum and maximum values of latitude and longitude
    double minLat = Double.MAX_VALUE, maxLat = Double.MIN_VALUE; // Latitude
    double minLon = Double.MAX_VALUE, maxLon = Double.MIN_VALUE; // Longitude
    for (double[] coord : coordinates) {
      minLat = Math.min(minLat, coord[0]);
      maxLat = Math.max(maxLat, coord[0]);
      minLon = Math.min(minLon, coord[1]);
      maxLon = Math.max(maxLon, coord[1]);
    }

    // Normalize each coordinate
    double[][] normalizedCoordinates = new double[coordinates.length][2];
    for (int i = 0; i < coordinates.length; i++) {
      double lat = coordinates[i][0];
      double lon = coordinates[i][1];
      double normalizedLat = (lat - minLat) / (maxLat - minLat);
      double normalizedLon = (lon - minLon) / (maxLon - minLon);
      normalizedCoordinates[i][0] = normalizedLat;
      normalizedCoordinates[i][1] = normalizedLon;
    }

    return normalizedCoordinates;
  }

  private static double[][] equirectangularProjection(double[][] coordinates) {
    // Find minimum and maximum latitude and longitude
    double minLat = Double.MAX_VALUE, maxLat = Double.MIN_VALUE;
    double minLon = Double.MAX_VALUE, maxLon = Double.MIN_VALUE;
    for (double[] coord : coordinates) {
      minLat = Math.min(minLat, coord[0]);
      maxLat = Math.max(maxLat, coord[0]);
      minLon = Math.min(minLon, coord[1]);
      maxLon = Math.max(maxLon, coord[1]);
    }

    // Project coordinates and scale to fit within (0, 1)
    double[][] projectedCoordinates = new double[coordinates.length][2];
    for (int i = 0; i < coordinates.length; i++) {
      double lat = coordinates[i][0];
      double lon = coordinates[i][1];

      // Project latitude and longitude
      double x = (lon - minLon) / (maxLon - minLon);
      double y = (lat - minLat) / (maxLat - minLat);

      // Adjust scaling to fit within (0, 1)
      double epsilon = 1e-6; // Small positive value
      x = epsilon + (1 - 2 * epsilon) * x;
      y = epsilon + (1 - 2 * epsilon) * y;

      projectedCoordinates[i][0] = x;
      projectedCoordinates[i][1] = y;
    }
    return projectedCoordinates;
  }

  public static void main(String[] args) {

    // Read coordinates from file
    double[][] coordinates = readCoordinatesFromFile("ISP data\\weights.intra", "ISP data\\city_coordinates.txt");
    // normalize
    coordinates = equirectangularProjection(coordinates);

    System.out.println("coors: ");
    for (int j = 0; j < coordinates[0].length; j++) {
      for (int i = 0; i < coordinates.length; i++) {
        System.out.println(coordinates[i][j]);
      }
    }
    // Calculate the number of points
    int n = coordinates.length;
    System.out.printf("n is %d ", n);

    // Calculate Euclidean distance matrix
    double[][] distanceMatrix = new double[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i == j) {
          distanceMatrix[i][j] = 0; // Distance from a point to itself is 0
        } else {
          distanceMatrix[i][j] = calculateEuclidDistance(coordinates[i], coordinates[j]);
        }
      }
    }
    int startNode = 0;
    TspDynamicProgrammingIterative solver = new TspDynamicProgrammingIterative(startNode, distanceMatrix);

    System.out.println("Tour: " + solver.getTour());
    System.out.println("Tour cost: " + solver.getTourCost());

   // Write coordinates to a text file
   String fileName = "points20_second.txt";
   try (BufferedWriter writer = new BufferedWriter(new FileWriter(fileName))) {
       for (int i = 0; i < coordinates.length; i++) {
           writer.write(String.format("%.6f %.6f", coordinates[i][0], coordinates[i][1]));
           writer.newLine(); // Move to the next line after writing each coordinate pair
       }
       System.out.println("Coordinates have been written to " + fileName);
   } catch (IOException e) {
       e.printStackTrace();
   }
  }
}

/*
 * 
 * 
 * // calc distance matrix
 * private static double[][] calculateDistanceMatrix(double[][] coordinates) {
 * int n = coordinates.length;
 * double[][] distanceMatrix = new double[n][n];
 * for (int i = 0; i < n; i++) {
 * for (int j = 0; j < n; j++) {
 * if (i == j) {
 * distanceMatrix[i][j] = 0; // Distance from a point to itself is 0
 * } else {
 * distanceMatrix[i][j] = calculateEuclidDistance(coordinates[i],
 * coordinates[j]);
 * }
 * }
 * }
 * return distanceMatrix;
 * }
 * 
 * // Create adjacency matrix
 * int n = 6;
 * double[][] distanceMatrix = new double[n][n];
 * for (double[] row : distanceMatrix)
 * java.util.Arrays.fill(row, 10000);
 * distanceMatrix[5][0] = 10;
 * distanceMatrix[1][5] = 12;
 * distanceMatrix[4][1] = 2;
 * distanceMatrix[2][4] = 4;
 * distanceMatrix[3][2] = 6;
 * distanceMatrix[0][3] = 8;
 * 
 * 
 * private static double[][] normalizeCoordinates(double[][] coordinates) {
 * // Find minimum and maximum values of latitude and longitude
 * double minLat = Double.MAX_VALUE, maxLat = Double.MIN_VALUE; // Latitude
 * double minLon = Double.MAX_VALUE, maxLon = Double.MIN_VALUE; // Longtitude
 * for (double[] coord : coordinates) {
 * minLat = Math.min(minLat, coord[0]);
 * maxLat = Math.max(maxLat, coord[0]);
 * minLon = Math.min(minLon, coord[1]);
 * maxLon = Math.max(maxLon, coord[1]);
 * }
 * 
 * // Normalize each coordinate
 * double[][] normalizedCoordinates = new double[coordinates.length][2];
 * for (int i = 0; i < coordinates.length; i++) {
 * double lat = coordinates[i][0];
 * double lon = coordinates[i][1];
 * double normalizedLat = (lat - minLat) / (maxLat - minLat);
 * double normalizedLon = (lon - minLon) / (maxLon - minLon);
 * normalizedCoordinates[i][0] = normalizedLat;
 * normalizedCoordinates[i][1] = normalizedLon;
 * }
 * return normalizedCoordinates;
 * }
 */
