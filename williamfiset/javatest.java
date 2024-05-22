package williamfiset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class javatest {

    private static double calculateEuclidDistance(double[] point1, double[] point2) {
        double dx = point1[0] - point2[0];
        double dy = point1[1] - point2[1];
        return Math.sqrt(dx * dx + dy * dy);
    }

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
            int count = 0;
            double[][] coordinates = new double[20][2];
            while ((line = br.readLine()) != null && count < 20) {
                String[] parts = line.trim().split("\\s+");
                String city1 = parts[0].replaceAll("\\+", " ").split(",")[0].trim();
                String city2 = parts[1].replaceAll("\\+", " ").split(",")[0].trim();
                
                
                if (cityCoordinates.containsKey(city1) && cityCoordinates.containsKey(city2)) {
                    coordinates[count] = cityCoordinates.get(city1);
                    count++;
                    if (count >= 20)
                        break;
                    coordinates[count] = cityCoordinates.get(city2);
                    count++;
                }
            }
            return coordinates;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static void main(String[] args) {

        // Read coordinates from file
        double[][] coordinates = readCoordinatesFromFile("ISP data\\weights.intra", "ISP data\\city_coordinates.txt");

        // Calculate the number of points
        int n = coordinates.length;

        System.out.println("Coordinates: \n");
        for (int j = 0; j < coordinates[0].length; j++) {
            for (int i = 0; i < coordinates.length; i++) {
                System.out.println(coordinates[i][j]);
            }
            System.out.println('\n');
        }

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

        // Prints: [0, 3, 2, 4, 1, 5, 0]
        System.out.println("Tour: " + solver.getTour());

        // Print: 42.0
        System.out.println("Tour cost: " + solver.getTourCost());
    }
}

/*
 * String city1 = parts[0].replace("+", " "); // "Townsville, Austrailia4282"
 * String[] city1parts = city1.split(","); // Split the input string based on
 * the comma
 * String firstPart = city1parts[0].trim();
 * city1 = firstPart;
 * 
 * String city2 = parts[1].replace("+", " ");
 * String[] city2parts = city1.split(","); // Split the input string based on
 * the comma
 * String firstPart2 = city2parts[0].trim();
 * city2 = firstPart2;
 */