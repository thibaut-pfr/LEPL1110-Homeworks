#include <stdio.h>
#include <math.h>
#include "glfem.h"

double integrate(double x[3], double y[3], double (*f) (double, double)) {

    double xi[3] = {1./6, 1./6, 2./3};
    double psi[3] = {1./6, 2./3, 1./6};

    double xLoc;
    double yLoc;
    double I = 0;

    double jacobian_matrix_determinant = fabs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]));
    //printf("Determinant of the transformation matrix: %f\n", jacobian_matrix_determinant);
    
    for (int i = 0; i < 3; i++) {
        xLoc = x[0] * (1 - xi[i] - psi[i]) + x[1] * xi[i] + x[2] * psi[i];
        yLoc = y[0] * (1 - xi[i] - psi[i]) + y[1] * xi[i] + y[2] * psi[i];
        //printf("x: %f, y: %f\n", xLoc, yLoc);
        I += f(xLoc,yLoc) * 1./6;
    }

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
        
    return I * jacobian_matrix_determinant;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n) {

    if (n == 0) {
    return integrate(x, y, f);
    }

    else {
        double xi[6] = {0.0, 0.5, 0.0, 0.5, 1.0, 0.0};
        double psi[6] = {0.0, 0.0, 0.5, 0.5, 0.0, 1.0};
        int triangles[4][3] = {{0,1,2},{2,3,5},{1,2,3},{1,3,4}};

        double newxLoc[3];
        double newyLoc[3];
        double I = 0;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                newxLoc[j] = x[0] * (1 - xi[triangles[i][j]] - psi[triangles[i][j]]) + x[1] * xi[triangles[i][j]] + x[2] * psi[triangles[i][j]];
                newyLoc[j] = y[0] * (1 - xi[triangles[i][j]] - psi[triangles[i][j]]) + y[1] * xi[triangles[i][j]] + y[2] * psi[triangles[i][j]];
            }
            I += integrateRecursive(newxLoc, newyLoc, f, n - 1);
        }
        return I;
    }
}