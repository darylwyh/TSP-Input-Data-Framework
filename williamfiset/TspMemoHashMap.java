package williamfiset;

import com.sun.management.OperatingSystemMXBean;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.*;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.*;

public class TspMemoHashMap {

    private final int N, start;
    private final double[][] distance;
    private List<Integer> tour = new ArrayList<>();
    private double minTourCost = Double.POSITIVE_INFINITY;
    private boolean ranSolver = false;

    public TspMemoHashMap(double[][] distance) {
        this(0, distance);
    }

    public TspMemoHashMap(int start, double[][] distance) {
        N = distance.length;

        if (N <= 2) {
            throw new IllegalStateException("N <= 2 not yet supported.");
        }
        if (N != distance[0].length) {
            throw new IllegalArgumentException("Matrix must be square (n x n)");
        }
        if (start < 0 || start >= N) {
            throw new IllegalArgumentException("Invalid start node.");
        }

        this.start = start;
        this.distance = distance;
    }

    // Returns the optimal tour for the traveling salesman problem.
    public List<Integer> getTour() {
        if (!ranSolver) {
            solve();
        }
        return tour;
    }

    // Returns the minimal tour cost.
    public double getTourCost() {
        if (!ranSolver) {
            solve();
        }
        return minTourCost;
    }

    private void solve() {

        if (ranSolver) {
            return;
        }

        final int END_STATE = (1 << N) - 1;
        Map<Integer, Double> memo = new HashMap<>();

        // Add all outgoing edges from the starting node to memo table.
        for (int end = 0; end < N; end++) {
            if (end == start) {
                continue;
            }
            memo.put((1 << start) | (1 << end) | (end << 16), distance[start][end]);
        }

        for (int r = 3; r <= N; r++) {
            for (int subset : combinations(r, N)) {
                if (notIn(start, subset)) {
                    continue;
                }
                for (int next = 0; next < N; next++) {
                    if (next == start || notIn(next, subset)) {
                        continue;
                    }
                    int subsetWithoutNext = subset ^ (1 << next);
                    double minDist = Double.POSITIVE_INFINITY;
                    for (int end = 0; end < N; end++) {
                        if (end == start || end == next || notIn(end, subset)) {
                            continue;
                        }
                        double newDistance = memo.getOrDefault((subsetWithoutNext | (end << 16)), Double.POSITIVE_INFINITY) + distance[end][next];
                        if (newDistance < minDist) {
                            minDist = newDistance;
                        }
                    }
                    memo.put((subset | (next << 16)), minDist);
                }
            }
        }

        // Connect tour back to starting node and minimize cost.
        for (int i = 0; i < N; i++) {
            if (i == start) {
                continue;
            }
            double tourCost = memo.getOrDefault((END_STATE | (i << 16)), Double.POSITIVE_INFINITY) + distance[i][start];
            if (tourCost < minTourCost) {
                minTourCost = tourCost;
            }
        }

        // Reconstruct TSP path.
        int lastIndex = start;
        int state = END_STATE;
        tour.add(start);

        for (int i = 1; i < N; i++) {

            int bestIndex = -1;
            double bestDist = Double.POSITIVE_INFINITY;
            for (int j = 0; j < N; j++) {
                if (j == start || notIn(j, state)) {
                    continue;
                }
                double newDist = memo.getOrDefault((state | (j << 16)), Double.POSITIVE_INFINITY) + distance[j][lastIndex];
                if (newDist < bestDist) {
                    bestIndex = j;
                    bestDist = newDist;
                }
            }

            tour.add(bestIndex);
            state = state ^ (1 << bestIndex);
            lastIndex = bestIndex;
        }

        tour.add(start);
        ranSolver = true;
    }

    private static boolean notIn(int elem, int subset) {
        return ((1 << elem) & subset) == 0;
    }

    // This method generates all bit sets of size n where r bits are set to one.
    public static List<Integer> combinations(int r, int n) {
        List<Integer> subsets = new ArrayList<>();
        combinations(0, 0, r, n, subsets);
        return subsets;
    }

    // To find all the combinations of size r we need to recurse until we have
    // selected r elements (i.e r = 0), otherwise if r != 0 then we still need to select
    // an element which is found after the position of our last selected element
    private static void combinations(int set, int at, int r, int n, List<Integer> subsets) {
        int elementsRemaining = n - at;
        if (elementsRemaining < r) {
            return;
        }

        if (r == 0) {
            subsets.add(set);
        } else {
            for (int i = at; i < n; i++) {
                set |= (1 << i);
                combinations(set, i + 1, r - 1, n, subsets);
                set &= ~(1 << i);
            }
        }
    }

    // calc Euclidean Distance between two points
    private static double calculateEuclidDistance(double[] point1, double[] point2) {
        double dx = point1[0] - point2[0];
        double dy = point1[1] - point2[1];
        return Math.sqrt(dx * dx + dy * dy);
    }

    // Function to read data from file, for RocketFuel 
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
                cityCoordinates.put(cityName, new double[]{latitude, longitude});
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        // Read from weights.intra file
        try (BufferedReader br = new BufferedReader(new FileReader(cityFile))) {
            int pointsSize = 20;
            String line;
            Set<String> uniqueCoordinates = new HashSet<>(); // Set to store unique coordinates
            List<double[]> uniqueCoordinatesList = new ArrayList<>(); // List to store unique coordinates in order
            while ((line = br.readLine()) != null && uniqueCoordinates.size() < pointsSize) {
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

    // Reads coordinates from a .tsp file
    public static double[][] readtsp(String filename, int numPoints) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(filename));
        List<double[]> coordinatesList = new ArrayList<>();
        String line;
        boolean readNodes = false;
        while ((line = br.readLine()) != null) {
            if (line.trim().equals("NODE_COORD_SECTION")) {
                readNodes = true;
                continue;
            }
            if (readNodes) {
                if (line.trim().equals("EOF")) {
                    break;
                }
                String[] parts = line.trim().split("\\s+");
                double x = Double.parseDouble(parts[1]);
                double y = Double.parseDouble(parts[2]);
                coordinatesList.add(new double[]{x, y});
                if (coordinatesList.size() == numPoints) {
                    break;
                }
            }
        }
        br.close();
        return coordinatesList.toArray(new double[0][0]);
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

    private static void printHeapUsage() {
        MemoryMXBean memoryMXBean = ManagementFactory.getMemoryMXBean();
        MemoryUsage heapMemoryUsage = memoryMXBean.getHeapMemoryUsage();

        long usedHeapMemory = heapMemoryUsage.getUsed();
        long maxHeapMemory = heapMemoryUsage.getMax();

        System.out.println("Used heap memory: " + usedHeapMemory / (1024 * 1024) + " MB");
        System.out.println("Max heap memory: " + maxHeapMemory / (1024 * 1024) + " MB");
    }

    private static void printCpuUsage() {
        OperatingSystemMXBean osBean = ManagementFactory.getPlatformMXBean(OperatingSystemMXBean.class);
        double cpuLoad = osBean.getProcessCpuLoad() * 100;
        System.out.println("CPU Load: " + cpuLoad + " %");
    }

    public static void main(String[] args) {
        // Measure start time
        long startTime = System.currentTimeMillis();

        // Read coordinates from file
        double[][] coordinates = readCoordinatesFromFile("ISP data\\weights.intra", "ISP data\\city_coordinates.txt");

        // normalize
        coordinates = equirectangularProjection(coordinates);
        //coordinates = MDSProjection.mdsProjection(coordinates);

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

        // Write coordinates to a text file to pointsSize
        String fileName = "points20_RocketFuel.txt";
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fileName))) {
            for (int i = 0; i < coordinates.length; i++) {
                writer.write(String.format("%.6f %.6f", coordinates[i][0], coordinates[i][1]));
                writer.newLine(); // Move to the next line after writing each coordinate pair
            }
            System.out.println("Coordinates have been written to " + fileName);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Measure end time
        long endTime = System.currentTimeMillis();
        long elapsedTime = endTime - startTime;

        // Print elapsed time
        System.out.println("Elapsed time: " + elapsedTime + " milliseconds");

        // Print final heap usage  // Print CPU usage
        printHeapUsage();
        printCpuUsage();

    }
}
