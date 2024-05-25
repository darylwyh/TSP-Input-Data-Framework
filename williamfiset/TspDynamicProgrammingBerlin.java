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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List; // Make sure to import the correct MXBean for CPU load

public class TspDynamicProgrammingBerlin {

    private final int N, start;
    private final double[][] distanceMatrix;
    private List<Integer> tour = new ArrayList<>();
    private double minTourCost = Double.POSITIVE_INFINITY;
    private boolean ranSolver = false;

    public TspDynamicProgrammingBerlin(double[][] coordinates) {
        this(0, coordinates);
    }

    public TspDynamicProgrammingBerlin(int start, double[][] coordinates) {
        N = coordinates.length;

        if (N <= 2) {
            throw new IllegalStateException("N <= 2 not yet supported.");
        }
        if (start < 0 || start >= N) {
            throw new IllegalArgumentException("Invalid start node.");
        }

        // Create distance matrix based on pairwise distances between coordinates
        distanceMatrix = new double[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                distanceMatrix[i][j] = distance(coordinates[i], coordinates[j]);
            }
        }

        this.start = start;
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

    // Solves the traveling salesman problem and caches solution.
    public void solve() {

        if (ranSolver) {
            return;
        }

        final int END_STATE = (1 << N) - 1;
        Double[][] memo = new Double[N][1 << N];

        // Add all outgoing edges from the starting node to memo table.
        for (int end = 0; end < N; end++) {
            if (end == start) {
                continue;
            }
            memo[end][(1 << start) | (1 << end)] = distanceMatrix[start][end];
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
                        double newDistance = memo[end][subsetWithoutNext] + distanceMatrix[end][next];
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
            if (i == start) {
                continue;
            }
            double tourCost = memo[i][END_STATE] + distanceMatrix[i][start];
            if (tourCost < minTourCost) {
                minTourCost = tourCost;
            }
        }

        int lastIndex = start;
        int state = END_STATE;
        tour.add(start);

        // Reconstruct TSP path from memo table.
        for (int i = 1; i < N; i++) {

            int bestIndex = -1;
            double bestDist = Double.POSITIVE_INFINITY;
            for (int j = 0; j < N; j++) {
                if (j == start || notIn(j, state)) {
                    continue;
                }
                double newDist = memo[j][state] + distanceMatrix[j][lastIndex];
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
        Collections.reverse(tour);

        ranSolver = true;
    }

    private double distance(double[] coord1, double[] coord2) {
        // Euclidean distance formula
        double xi = coord1[0];
        double yi = coord1[1];
        double xj = coord2[0];
        double yj = coord2[1];
        return Math.sqrt(Math.pow(xi - xj, 2) + Math.pow(yi - yj, 2));
    }

    private static boolean notIn(int elem, int subset) {
        return ((1 << elem) & subset) == 0;
    }

    // This method generates all bit sets of size n where r bits
    // are set to one. The result is returned as a list of integer masks.
    public static List<Integer> combinations(int r, int n) {
        List<Integer> subsets = new ArrayList<>();
        combinations(0, 0, r, n, subsets);
        return subsets;
    }

    // To find all the combinations of size r we need to recurse until we have
    // selected r elements (aka r = 0), otherwise if r != 0 then we still need to
    // select
    // an element which is found after the position of our last selected element
    private static void combinations(int set, int at, int r, int n, List<Integer> subsets) {

        // Return early if there are more elements left to select than what is
        // available.
        int elementsLeftToPick = n - at;
        if (elementsLeftToPick < r) {
            return;
        }

        // We selected 'r' elements so we found a valid subset!
        if (r == 0) {
            subsets.add(set);
        } else {
            for (int i = at; i < n; i++) {
                // Try including this element
                set ^= (1 << i);

                combinations(set, i + 1, r - 1, n, subsets);

                // Backtrack and try the instance where we did not include this element
                set ^= (1 << i);
            }
        }
    }

    // Reads coordinates from a .tsp file
    public static double[][] readCoordinates(String filename, int numPoints) throws IOException {
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
        try {
            // Measure start time
            long startTime = System.currentTimeMillis();
            // Print initial heap usage
            printHeapUsage();

            int n = 23; // Choose up to 25 nodes
            String filename = "ISP data\\berlin52.tsp";
            //String filename = "ISP data\\ch71009.tsp";
            double[][] coordinates = readCoordinates(filename, n);
            coordinates = equirectangularProjection(coordinates);

            //String fileName = "points24_china.txt";
            String writefileName = "points23_berlin.txt";
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(writefileName))) {
                for (int i = 0; i < coordinates.length; i++) {
                    writer.write(String.format("%.6f %.6f", coordinates[i][0], coordinates[i][1]));
                    writer.newLine(); // Move to the next line after writing each coordinate pair
                }
                System.out.println("Coordinates have been written to " + writefileName);
            } catch (IOException e) {
                e.printStackTrace();
            }

            int startNode = 0;

            TspDynamicProgrammingModified solver = new TspDynamicProgrammingModified(startNode, coordinates);

            // Print tour and tour cost
            System.out.println("Tour: " + solver.getTour());
            System.out.println("Tour cost: " + solver.getTourCost());

            // Measure end time
            long endTime = System.currentTimeMillis();
            long elapsedTime = endTime - startTime;

            // Print elapsed time
            System.out.println("Elapsed time: " + elapsedTime + " milliseconds");

            // Print final heap usage  // Print CPU usage
            printHeapUsage(); printCpuUsage();

            

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
