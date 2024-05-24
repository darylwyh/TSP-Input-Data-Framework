package DataNormalization;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import williamfiset.TspDynamicProgrammingModified;

public class norm {
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
                if (line.trim().equals("EOF"))
                    break;
                String[] parts = line.trim().split("\\s+");
                double x = Double.parseDouble(parts[1]);
                double y = Double.parseDouble(parts[2]);
                coordinatesList.add(new double[] { x, y });
                if (coordinatesList.size() == numPoints)
                    break;
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

    public static void main(String[] args) {
        try {
            int n = 50;  
            // String filename = "ISP data\\berlin52.tsp";
            String filename = "ISP data\\ch71009.tsp";
            double[][] coordinates = readCoordinates(filename, n);
            coordinates = equirectangularProjection(coordinates);

            String fileName = "points50_china.txt";
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(fileName))) {
                for (int i = 0; i < coordinates.length; i++) {
                    writer.write(String.format("%.6f %.6f", coordinates[i][0], coordinates[i][1]));
                    writer.newLine(); // Move to the next line after writing each coordinate pair
                }
                System.out.println("Coordinates have been written to " + fileName);
            } catch (IOException e) {
                e.printStackTrace();
            }
 
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
