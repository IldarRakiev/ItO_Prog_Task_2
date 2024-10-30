package org.example;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.ejml.simple.SimpleMatrix;
public class InteriorPointMethod {
    public static SimpleMatrix interiorPointMethod(SimpleMatrix A, SimpleMatrix b,
                                                   SimpleMatrix c, SimpleMatrix x0, double alpha, double epsilon) {
        SimpleMatrix x = x0.copy();
        SimpleMatrix x_temp = x0.copy();
        SimpleMatrix x_new;
        boolean solved = false;
        while (!solved) {
            try {
                SimpleMatrix D = SimpleMatrix.diag(x.getDDRM().getData());
                SimpleMatrix A_temp = A.mult(D);
                SimpleMatrix Dc = D.mult(c);
                SimpleMatrix A_temp_t = A_temp.transpose();
                SimpleMatrix AAt = A_temp.mult(A_temp_t);
                SimpleMatrix AAtInv = AAt.invert();
                Page 2
                SimpleMatrix P =
                        SimpleMatrix.identity(x.numRows()).minus(A_temp_t.mult(AAtInv).mult(A_temp));
                SimpleMatrix Cp = P.mult(Dc);
                double v = 0;
                for (int i = 0; i < Cp.numRows(); i++) {
                    double val = Cp.get(i);
                    if (val < 0) {
                        v = Math.max(v, -val);
                    }
                }
                if (v < 1e-6) {
                    return x;
                }
                for (int i = 0; i < x_temp.numRows(); i++) {
                    x_temp.set(i, 1 + alpha * Cp.get(i) / v);
                }
                x_new = D.mult(x_temp);
                solved = true;
                for (int i = 0; i < x.numRows(); i++) {
                    if (Math.abs(x_new.get(i) - x.get(i)) > epsilon) {
                        solved = false;
                        break;
                    }
                }
                x = x_new;
            } catch (Exception e) {
                System.out.println("The method is not applicable!");
                solved = true;
                return null;
            }
        }
        return x;
    }
    public static BigDecimal round(double value){
        BigDecimal bd = BigDecimal.valueOf(value);
        bd = bd.setScale(5, RoundingMode.HALF_UP);
        return bd;
    }
    public static void printVector(SimpleMatrix vector){
        System.out.print("[");
        for (int i = 0; i < vector.numRows() - 1; i++) {
            System.out.print(round(vector.get(i, 0)) + ", ");
        }
        System.out.println(round(vector.get(vector.numRows() - 1)) + "]");
    }
    public static void tester() {
        double[][][] aDataList = {
                {{4, 6}},
                {{1, 1, 2}, {2, 1, 1}, {3, 4, 5}},
                {{1, 2, 1, 0}, {3, 1, 0, 1}},
                {{1, 2, 1 ,0, 0}, {2, 3, 0, 1, 0}, {1, 1, 0, 0, 1}},
                {{1, 2, 3}},
                {{4, 6, 1}}
        };
        double[][] bDataList = {
                {30}, {9, 7, 26}, {4, 3}, {5, 6, 7}, {6}, {31}
        };
        double[][] cDataList = {
                {3, 4}, {9, 10, 16}, {4, 3, 0, 0}, {3, 4, 0, 0, 0}, {2, 5, 7}, {3, 3}
        };
        double[][] x0DataList = {
                {3, 3}, {1, 2, 3}, {0.5, 1, 1.5, 0.5}, {2, 0.5, 2, 0.5, 4.5}, {1, 1,
                1}, {1, 1}
        };
        double[] epsilonList = {
                1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6
        };
        double alpha1 = 0.5;
        double alpha2 = 0.9;
        int precision = 6;
        for (int i = 0; i < aDataList.length; i++){
            SimpleMatrix A = new SimpleMatrix(aDataList[i]);
            SimpleMatrix b = new SimpleMatrix(bDataList[i].length, 1, true,
                    bDataList[i]);
            SimpleMatrix c = new SimpleMatrix(cDataList[i].length, 1, true,
                    cDataList[i]);
            SimpleMatrix x0 = new SimpleMatrix(x0DataList[i].length, 1, true,
                    x0DataList[i]);
            System.out.println("Test " + (i+1) + ":");
            SimpleMatrix result1 = interiorPointMethod(A, b, c, x0, alpha1,
                    epsilonList[i]);
            SimpleMatrix result2 = null;
            if (result1 != null && !result1.hasUncountable()){
                result2 = interiorPointMethod(A, b, c, x0, alpha2, epsilonList[i]);
            }
            if (result1 == null) {
            } else if (result1.hasUncountable()) {
                System.out.println("The problem does not have solution!");
            } else {
                System.out.println("A vector of decision variables");
                System.out.print(" = 0.5: ");
                printVector(result1);
                System.out.print(" = 0.9: ");
                printVector(result2);
                System.out.print("by simplex method: ");
                Object[] result = simplexMethod(cDataList[i], aDataList[i],
                        bDataList[i], precision, epsilonList[i]);
                if (result[0] instanceof String) {
                    System.out.println(result[0]);
                } else {
                    double[] solution = (double[]) result[0];
                    System.out.println(Arrays.toString(solution));
                }
                System.out.println();
                System.out.println("Maximum value of the objective function");
                System.out.print(" = 0.5: ");
                System.out.println(round(c.dot(result1)));
                System.out.print(" = 0.9: ");
                System.out.println(round(c.dot(result2)));
            }
            System.out.println("------------------------------");
        }
    }
    public static void main(String[] args) {
        tester();
    }
    public static Object[] initializeTableau(double[] C, double[][] A, double[] b) {
        int nConstraints = A.length;
        int nVariables = C.length;
        List<List<Double>> tableau = new ArrayList<>();
        for (int i = 0; i < nConstraints; i++) {
            List<Double> row = new ArrayList<>();
            for (double value : A[i]) {
                row.add(value);
            }
            for (int j = 0; j < nConstraints; j++) {
                row.add(0.0);
            }
            row.add(b[i]);
            row.set(nVariables + i, 1.0);
            tableau.add(row);
        }
        List<Double> costRow = new ArrayList<>();
        for (double ci : C) {
            costRow.add(-ci);
        }
        for (int j = 0; j < nConstraints + 1; j++) {
            costRow.add(0.0);
        }
        tableau.add(costRow);
        List<Integer> basis = new ArrayList<>();
        for (int i = 0; i < nConstraints; i++) {
            basis.add(nVariables + i);
        }
        return new Object[]{tableau, basis};
    }
    public static Integer pivotColumn(List<List<Double>> tableau) {
        List<Double> lastRow = tableau.get(tableau.size() - 1);
        double minValue = Double.POSITIVE_INFINITY;
        Integer colIndex = null;
        for (int i = 0; i < lastRow.size() - 1; i++) {
            if (lastRow.get(i) < minValue) {
                minValue = lastRow.get(i);
                colIndex = i;
            }
        }
        if (minValue >= 0) {
            return null;
        }
        return colIndex;
    }
    public static Integer pivotRow(List<List<Double>> tableau, int col) {
        int nConstraints = tableau.size() - 1;
        List<Double> rhs = new ArrayList<>();
        List<Double> lhs = new ArrayList<>();
        List<Double> ratios = new ArrayList<>();
        for (int i = 0; i < nConstraints; i++) {
            rhs.add(tableau.get(i).get(tableau.get(i).size() - 1));
            lhs.add(tableau.get(i).get(col));
        }
        boolean unbounded = true;
        for (int i = 0; i < nConstraints; i++) {
            if (lhs.get(i) > 0) {
                ratios.add(rhs.get(i) / lhs.get(i));
                unbounded = false;
            } else {
                ratios.add(Double.POSITIVE_INFINITY);
            }
        }
        if (unbounded) {
            return null;
        }
        double minRatio = Double.POSITIVE_INFINITY;
        Integer rowIndex = null;
        for (int i = 0; i < ratios.size(); i++) {
            if (ratios.get(i) < minRatio) {
                minRatio = ratios.get(i);
                rowIndex = i;
            }
        }
        if (minRatio == Double.POSITIVE_INFINITY) {
            return null;
        }
        return rowIndex;
    }
    public static void pivot(List<List<Double>> tableau, int row, int col) {
        double pivotValue = tableau.get(row).get(col);
        for (int j = 0; j < tableau.get(row).size(); j++) {
            tableau.get(row).set(j, tableau.get(row).get(j) / pivotValue);
        }
        for (int i = 0; i < tableau.size(); i++) {
            if (i != row) {
                double factor = tableau.get(i).get(col);
                for (int j = 0; j < tableau.get(i).size(); j++) {
                    tableau.get(i).set(j, tableau.get(i).get(j) - factor *
                            tableau.get(row).get(j));
                }
            }
        }
    }
    public static Object[] simplexMethod(double[] C, double[][] A, double[] b, int
            precision, double epsilon) {
        Object[] initResult = initializeTableau(C, A, b);
        List<List<Double>> tableau = (List<List<Double>>) initResult[0];
        List<Integer> basis = (List<Integer>) initResult[1];
        double prevZ = Double.NEGATIVE_INFINITY;
        while (true) {
            Integer col = pivotColumn(tableau);
            if (col == null) {
                double[] solution = new double[C.length];
                for (int i = 0; i < basis.size(); i++) {
                    if (basis.get(i) < C.length) {
                        solution[basis.get(i)] =
                        tableau.get(i).get(tableau.get(i).size() - 1);
                    }
                }
                solution = Arrays.stream(solution)
                        .map(val -> Math.round(val * Math.pow(10, precision)) /
                                Math.pow(10, precision))
                        .toArray();
                double optimalValue = Math.round(-tableau.get(tableau.size() -
                        1).get(tableau.get(0).size() - 1) * Math.pow(10, precision)) /
                        Math.pow(10, precision);
                return new Object[]{solution, optimalValue};
            }
            Integer row = pivotRow(tableau, col);
            if (row == null) {
                return new Object[]{"Unbounded solution"};
            }
            basis.set(row, col);
            pivot(tableau, row, col);
            double currentZ = -tableau.get(tableau.size() -
                    1).get(tableau.get(0).size() - 1);
            if (Math.abs(currentZ - prevZ) < epsilon) {
                if (tableau.get(tableau.size() - 1).stream().anyMatch(value -> value <
                        0)) {
                    return new Object[]{"The method is not applicable!"};
                }
            }
            prevZ = currentZ;
        }
    }
}