public class Matrix {
    public static final double[][] IDENTITY_MATRIX = {
            {1, 0, 0, 0, 0},
            {0, 1, 0, 0, 0},
            {0, 0, 1, 0, 0},
            {0, 0, 0, 1, 0},
            {0, 0, 0, 0, 1}
    };

    public static double[][] changeRow(double[][] mtr, double[] newRow, int index) throws ArrayIndexOutOfBoundsException {
        mtr[index] = newRow;
        return mtr;
    }

    public static double norm(double[][] mtr) throws ArrayIndexOutOfBoundsException {
        double sum1 = 0;
        double sum2 = 0;
        for (int i = 0; i < mtr.length; i++) {
            for (int j = 0; j < mtr[0].length; j++) {
                sum1 += Math.abs(mtr[i][j]);
            }
            if (sum1 > sum2) {
                sum2 = sum1;
                sum1 = 0;
            }
        }
        return sum2;
    }

    public static void outputMatrix(double[][] mtr) throws ArrayIndexOutOfBoundsException {
        for (int i = 0; i < mtr.length; i++) {
            for (int j = 0; j < mtr[0].length; j++) {
                System.out.printf("%7.4f ", mtr[i][j]);
            }
            System.out.println();
        }
    }

    public static double[][] copyMatrix(double[][] mtr) throws ArrayIndexOutOfBoundsException {
        double[][] tmp = new double[mtr.length][mtr[0].length];
        for (int i = 0; i < mtr.length; i++) {
            for (int j = 0; j < mtr[0].length; j++) {
                tmp[i][j] = mtr[i][j];
            }
        }
        return tmp;
    }

    public static double[][] transposeMatrix(double[][] mtr) throws ArrayIndexOutOfBoundsException {
        double[][] tmp = new double[mtr[0].length][mtr.length];
        for (int i = 0; i < mtr.length; i++) {
            for (int j = 0; j < mtr[0].length; j++) {
                tmp[j][i] = mtr[i][j];
            }
        }
        return tmp;
    }

    public static double[][] matrixCompositionMatrix(double[][] a, double[][] b) throws ArrayIndexOutOfBoundsException {
        if (a[0].length == b.length) {
            double[][] result = new double[a.length][b[0].length];
            for (int i = 0; i < a.length; i++)
                for (int j = 0; j < b[0].length; j++) {
                    result[i][j] = 0;
                    for (int k = 0; k < b.length; k++) {
                        result[i][j] += a[i][k] * b[k][j];
                    }
                }
            return result;
        } else {
            throw new ArrayIndexOutOfBoundsException("Incorrect length of matrix.");
        }
    }

    public static double[][] matrixMinusMatrix(double[][] mtr1, double[][] mtr2) throws ArrayIndexOutOfBoundsException {
        if (mtr1.length == mtr2.length && mtr1[0].length == mtr2[0].length) {
            double[][] mtr = new double[mtr1.length][mtr1.length];
            for (int i = 0; i < mtr1.length; i++) {
                for (int j = 0; j < mtr1.length; j++) {
                    mtr[i][j] = mtr1[i][j] - mtr2[i][j];
                }
            }
            return mtr;
        } else {
            throw new ArrayIndexOutOfBoundsException("Incorrect length of matrix.");
        }
    }

    public static double[][] matrixPlusMatrix(double[][] mtr1, double[][] mtr2) throws ArrayIndexOutOfBoundsException {
        if (mtr1.length == mtr2.length && mtr1[0].length == mtr2[0].length) {
            double[][] mtr = new double[mtr1.length][mtr1.length];
            for (int i = 0; i < mtr1.length; i++) {
                for (int j = 0; j < mtr1.length; j++) {
                    mtr[i][j] = mtr1[i][j] + mtr2[i][j];
                }
            }
            return mtr;
        } else {
            throw new ArrayIndexOutOfBoundsException("Incorrect length of matrix.");
        }
    }
}