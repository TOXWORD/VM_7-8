import java.util.Arrays;

public class Krylov {
    private int length;
    private double[][] matrix;
    private double[][] symmetricMatrix;
    private double[][] sharedMatrix;
    private double[] polynomial;
    private double[] residual;
    private double[] eigenvector;
    private double eigenvalue;
    private double[] selfVector;

    public Krylov() {
        length = 5;
        matrix = new double[][]{{0.4974, 0.0000, -0.1299, 0.0914, 0.1523},
                {-0.0305, 0.3248, 0.0000, -0.0619, 0.0203},
                {0.0102, -0.0914, 0.5887, 0.0112, 0.0355},
                {0.0305, 0.0000, -0.0741, 0.5887, 0.0000},
                {0.0203, -0.0305, 0.1472, -0.0122, 0.4263}};
        symmetricMatrix = Matrix.matrixCompositionMatrix(Matrix.transposeMatrix(matrix), matrix);
        sharedMatrix = new double[length][length + 1];
        polynomial = new double[length + 1];
        residual = new double[length];
        eigenvector = new double[length];
        eigenvalue = 0;
    }

    public void findSolution() {
        convertToEquationMatrix();
        findPolynomial();
        setSelfVector();
        outEigenVector();
        findResidual();
    }

    private void findEigenvector() {
        /*Looking for coefficients for vectorsY*/
        eigenvector = new double[length];
        double[] coefficients = new double[length];
        coefficients[0] = 1;
        for (int i = 1; i < coefficients.length; i++) {
            coefficients[i] = coefficients[i - 1] * eigenvalue + polynomial[i];
        }
        /*Rolling around*/
        double[][] vectorsY = new double[length][length];
        for (int i = 0; i < vectorsY.length; i++) {
            for (int j = 0; j < vectorsY[0].length; j++) {
                vectorsY[j][i] = sharedMatrix[i][j];
            }
        }
        /*Looking for eigenvector*/
        for (int i = 0; i < eigenvector.length; i++) {
            eigenvector = Vector.vectorPlusVector(Vector.vectorCompositionNumber(vectorsY[i], coefficients[i]), eigenvector);
        }

        eigenvector = Vector.orthonormalize(eigenvector);
    }

    private void findPolynomial() {
        /*Find coefficients*/
        MSquareRoot msr = new MSquareRoot(sharedMatrix);
        msr.findSolution();
        double[] coefficients = msr.getAnswers();
        /*Fill polynomial*/
        polynomial[0] = 1;
        for (int i = 1; i < polynomial.length; i++) {
            polynomial[i] = -coefficients[i - 1];
        }
    }

    private void findResidual() {
        residual = Vector.matrixCompositionVector(symmetricMatrix, eigenvector);
        residual = Vector.vectorMinusVector(residual, Vector.vectorCompositionNumber(eigenvector, eigenvalue));
    }

    private void convertToEquationMatrix() {
        double[] vector = new double[length];
        vector[0] = 1;
        for (int i = length - 1; i >= 0; i--) {
            putVectorToEquationMatrix(vector, i);
            vector = Vector.matrixCompositionVector(symmetricMatrix, vector);
        }
        putVectorToEquationMatrix(vector, length);
    }

    private void putVectorToEquationMatrix(double[] vector, int col) {
        for (int i = 0; i < length; i++) {
            sharedMatrix[i][col] = vector[i];
        }
    }

    public void outputPolynomial() {
        System.out.println("Polynomial: ");
        for (int i = 0; i < polynomial.length; i++) {
            System.out.printf("%.5f ", polynomial[i]);
        }
    }

    private void setSelfVector() {
        selfVector = new double[]{0.27249, 0.09105, 0.46520, 0.36970, 0.12232};
    }

    public void outEigenVector() {
        System.out.println("Self vectors:");
        for (int i = 0; i < length; i++) {
            eigenvalue = selfVector[i];
            findEigenvector();
            System.out.printf("%.5f", eigenvalue);
            System.out.println();
            System.out.println(Arrays.toString(eigenvector));
            findResidual();
            System.out.println("Residual norm: " + Vector.norm(residual));
            double res = 0;
            for (int j = 0; j < polynomial.length - 1; j++) {
                res += polynomial[j] * Math.pow(eigenvalue, 5 - j);
            }
            System.out.println("Polynomial residual: ");
            System.out.println(res);
            System.out.println();
        }
    }

}