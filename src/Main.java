public class Main {
    public static void main(String[] args) {

        System.out.println("Krylov: ");
        Krylov fb = new Krylov();
        fb.findSolution();
        fb.outputPolynomial();
        System.out.println();
        System.out.println();

        System.out.println("Danilevsky: ");
        Danilevsky d = new Danilevsky();
        d.findSolution();
//            Vector.outputVector(d.getEigenvector());
//            d.outputResidual();
//            d.outputPolynomial();
        //     System.out.println(Vector.norm(d.getResidual()));

    }
}