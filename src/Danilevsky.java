import java.util.Arrays;

public class Danilevsky {
    private double[][] matrix;
    private double[][] frobenius;
    private double[][] similarity;
    private double[][] multiply;
    private double[][] inverseMultiply;
    private double[] eigenvector;
    private double[] selfVector;
    private double eigenvalue;
    private double[] polynomial;
    private double[] residual;
    private int length;

    public Danilevsky() {
        length = 5;
        eigenvalue = 0;
        matrix = new double[][]{{0.4974, 0.0000, -0.1299, 0.0914, 0.1523},
                {-0.0305, 0.3248, 0.0000, -0.0619, 0.0203},
                {0.0102, -0.0914, 0.5887, 0.0112, 0.0355},
                {0.0305, 0.0000, -0.0741, 0.5887, 0.0000},
                {0.0203, -0.0305, 0.1472, -0.0122, 0.4263}};
        matrix = Matrix.matrixCompositionMatrix(Matrix.transposeMatrix(matrix), matrix);
        frobenius = Matrix.copyMatrix(matrix);
        eigenvector = new double[length];
        similarity = Matrix.copyMatrix(Matrix.IDENTITY_MATRIX);
        polynomial = new double[length + 1];
        selfVector = new double[length];
    }


    public void findSolution() {
        convertToFrobenius();
        findPolynomial();
        setSelfVector();
        outEigenVector();
        findResidual();
        outputPolynomial();
    }

    private void findEigenvector() {
        eigenvector[length - 1] = 1;
        for (int i = length - 2; i >= 0; i--) {
            eigenvector[i] = eigenvector[i + 1] * eigenvalue;
        }
        eigenvector = Vector.matrixCompositionVector(similarity, eigenvector);
        eigenvector = Vector.orthonormalize(eigenvector);

    }

    private void convertToFrobenius() {
        for (int i = length - 1; i >= 1; i--) {
            /*Annul matrices*/
            inverseMultiply = Matrix.copyMatrix(Matrix.IDENTITY_MATRIX);
            multiply = Matrix.copyMatrix(Matrix.IDENTITY_MATRIX);
            /*Find elementary matrices*/
            inverseMultiply = Matrix.changeRow(inverseMultiply, frobenius[i], i - 1);
            this.convertMultiply(frobenius[i], i - 1, i - 1);
            /*Similarity matrix*/
            similarity = Matrix.matrixCompositionMatrix(similarity, multiply);
            /*Multiplying*/
            frobenius = Matrix.matrixCompositionMatrix(inverseMultiply, frobenius);
            frobenius = Matrix.matrixCompositionMatrix(frobenius, multiply);
        }
    }

    private void convertMultiply(double[] row, int index, int divider) {
        for (int i = 0; i < length; i++) {
            if (i != index) {
                multiply[index][i] = -row[i] / row[divider];
            } else {
                multiply[index][i] = 1 / row[divider];
            }
        }
    }

    public void printF() {

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                System.out.print(frobenius[i][j] + " ");
            }
            System.out.println();
        }
    }

    private void findPolynomial() {
        polynomial[polynomial.length - 1] = 1;
        for (int i = 0; i < polynomial.length - 1; i++) {
            polynomial[i] = -frobenius[0][length - i - 1];
        }
    }

    public void outputPolynomial() {

        for (int i = 0; i < polynomial.length; i++) {
            System.out.printf("%.5f ", polynomial[polynomial.length - 1 - i]);
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
            for (int j = 1; j < polynomial.length; j++) {
                res += polynomial[j] * Math.pow(eigenvalue, j);
            }
            System.out.println("Polynomial residual: ");
            System.out.println(res);
            System.out.println();
        }
    }

    private void findResidual() {
        residual = Vector.matrixCompositionVector(matrix, eigenvector);
        residual = Vector.vectorMinusVector(residual, Vector.vectorCompositionNumber(eigenvector, eigenvalue));
    }

}