public class Vector {
    public static double[] orthonormalize(double[] vector) throws ArrayIndexOutOfBoundsException {
        double vectorLength = 0;
        for (int i = 0; i < vector.length; i++) {
            vectorLength += (vector[i] * vector[i]);
        }
        vectorLength = Math.sqrt(vectorLength);
        double[] tmp = Vector.vectorCompositionNumber(vector, 1 / vectorLength);
        return tmp;
    }

    public static double norm(double[] vector) throws ArrayIndexOutOfBoundsException {
        double norm = Math.abs(vector[0]);
        for (int i = 1; i < vector.length; i++) {
            if (Math.abs(vector[i]) > norm) {
                norm = Math.abs(vector[i]);
            }
        }
        return norm;
    }

    public static double[] vectorCompositionNumber(double[] vector, double coeff) throws ArrayIndexOutOfBoundsException {
        double tmp[] = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            tmp[i] = vector[i] * coeff;
        }
        return tmp;
    }

    public static double scalarComposition(double[] vector1, double[] vector2) throws ArrayIndexOutOfBoundsException {
        if (vector1.length == vector2.length) {
            double tmp = 0;
            for (int i = 0; i < vector1.length; i++) {
                tmp += (vector1[i] * vector2[i]);
            }
            return tmp;
        } else {
            throw new ArrayIndexOutOfBoundsException("Incorrect length of vector.");
        }
    }

    public static void outputVector(double[] vector) throws ArrayIndexOutOfBoundsException {
        if (vector[0] >= 0) {
            System.out.printf("(%8.6f; ", vector[0]);
        } else {
            System.out.printf("(%9.6f; ", vector[0]);
        }
        for (int i = 1; i < vector.length - 1; i++) {
            System.out.printf("%9.6f; ", vector[i]);
        }
        System.out.printf("%9.6f)%n", vector[vector.length - 1]);
    }

    public static double[] copyVector(double[] from) throws ArrayIndexOutOfBoundsException {
        double[] vector = new double[from.length];
        for (int i = 0; i < from.length; i++) {
            vector[i] = from[i];
        }
        return vector;
    }

    public static double[] matrixCompositionVector(double[][] mtr, double[] vec) throws ArrayIndexOutOfBoundsException {
        if (mtr.length == vec.length) {
            double[] tmp = new double[mtr.length];
            for (int i = 0; i < mtr.length; i++) {
                tmp[i] = 0;
                for (int j = 0; j < mtr.length; j++) {
                    tmp[i] += mtr[i][j] * vec[j];
                }
            }
            return tmp;
        } else {
            throw new ArrayIndexOutOfBoundsException("Incorrect length of matrix/vector.");
        }
    }

    public static double[] vectorMinusVector(double[] first, double[] second) throws ArrayIndexOutOfBoundsException {
        if (first.length == second.length) {
            double[] tmp = new double[first.length];
            for (int i = 0; i < first.length; i++) {
                tmp[i] = first[i] - second[i];
            }
            return tmp;
        } else {
            throw new ArrayIndexOutOfBoundsException("Incorrect length of vector.");
        }
    }

    public static double[] vectorPlusVector(double[] first, double[] second) throws ArrayIndexOutOfBoundsException {
        if (first.length == second.length) {
            double[] tmp = new double[first.length];
            for (int i = 0; i < first.length; i++) {
                tmp[i] = first[i] + second[i];
            }
            return tmp;
        } else {
            throw new ArrayIndexOutOfBoundsException("Incorrect length of vector.");
        }
    }
}