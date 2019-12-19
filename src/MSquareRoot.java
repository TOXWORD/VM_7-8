import java.io.File;
import java.io.IOException;
import java.util.Scanner;

public class MSquareRoot {
    private double[][] matrix;
    private double[] vectorB;
    private double[] answers;
    private byte rows;
    private byte cols;
    private double determinant;
    private double[] convertedVectorB;
    private double[][] symmetricMatrix;
    private double[][] L;
    private double[][] U;

    public MSquareRoot(double[][] sharedMatrix) {
        rows = 5;
        cols = 5;
        determinant = 1;
        matrix = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.matrix[i][j] = sharedMatrix[i][j];
            }
        }
        vectorB = new double[rows];
        for (int i = 0; i < rows; i++) {
            this.vectorB[i] = sharedMatrix[i][sharedMatrix[0].length - 1];
        }
        answers = new double[rows];
        convertedVectorB = new double[rows];
        symmetricMatrix = new double[rows][cols];
        L = new double[rows][cols];
        U = new double[rows][cols];
    }

    public void findSymmetric() throws ArrayIndexOutOfBoundsException {
        double[][] transposedMatrix = Matrix.transposeMatrix(matrix);
        symmetricMatrix = Matrix.matrixCompositionMatrix(transposedMatrix, matrix);
        convertedVectorB = Vector.matrixCompositionVector(transposedMatrix, vectorB);
    }

    public void findLU() throws ArrayIndexOutOfBoundsException {
        double sum = 0;
        for (int i = 0; i < symmetricMatrix.length; i++) {
            for (int k = 0; k < i; k++) {
                sum += U[k][i] * U[k][i];
            }
            U[i][i] = Math.sqrt(symmetricMatrix[i][i] - sum);
            for (int j = i + 1; j < symmetricMatrix.length; j++) {
                sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += U[k][i] * U[k][j];
                }
                U[i][j] = (symmetricMatrix[i][j] - sum) / U[i][i];
            }
            sum = 0;
        }
        L = Matrix.transposeMatrix(U);
        for (int i = 0; i < L.length; i++) {
            determinant *= L[i][i];
        }
    }

    public void findSolution() throws ArrayIndexOutOfBoundsException {
        this.findSymmetric();
        this.findLU();
        double[] vectorY = {0, 0, 0, 0, 0};
        double sum = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < i; j++) {
                sum += vectorY[j] * L[i][j];
            }
            vectorY[i] = (convertedVectorB[i] - sum) / L[i][i];
            sum = 0;
        }
        for (int i = rows - 1; i >= 0; i--) {
            for (int j = rows - 1; j >= i + 1; j--) {
                sum += answers[j] * U[i][j];
            }
            answers[i] = (vectorY[i] - sum) / U[i][i];
            sum = 0;
        }
    }

    public double[] getAnswers() {
        return answers;
    }
}