package williamfiset; 

public class javatest {
    public static void main(String[] args) {
        int n = 6;
        double[][] distanceMatrix = new double[n][n];
        for (double[] row : distanceMatrix)
            java.util.Arrays.fill(row, 10000);
        distanceMatrix[5][0] = 10;
        distanceMatrix[1][5] = 12;
        distanceMatrix[4][1] = 2;
        distanceMatrix[2][4] = 4;
        distanceMatrix[3][2] = 6;
        distanceMatrix[0][3] = 8;
 
        System.out.println(distanceMatrix.length);
    }
}