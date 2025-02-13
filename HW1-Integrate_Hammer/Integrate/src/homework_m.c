#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];
    double jacob = 0.0;

    double w[3] = {1.0/6.0, 1.0/6.0, 1.0/6.0};
    double xi[3] = {1.0/6.0,1.0/6.0,2.0/3.0};
    double eta[3] = {1.0/6.0,2.0/3.0,1.0/6.0};

    for(int i = 0; i < 3; i++){
        xLoc[i] = (1-xi[i]-eta[i])*x[0] + xi[i]*x[1] + eta[i]*x[2];
        yLoc[i] = (1-xi[i]-eta[i])*y[0] + xi[i]*y[1] + eta[i]*y[2];
    }

    // Jacobien de mise à l'échelle : rapport des surfaces avant et après standardisation
    jacob = fabs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]));
    for(int i = 0; i < 3; i++) I += w[i]*jacob*f(xLoc[i],yLoc[i]);
    

// Pour dessiner l'element, les sommets du triangle :-)
// Decommenter la ligne pour dessiner aussi les points d'integration
    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    
    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    double I = 0;

    // Integrate recursively the triangle
    if(n == 0){
        I = integrate(x,y,f);
    }else{
        double x0[3] = {x[0],(x[0]+x[1])/2,(x[0]+x[2])/2};
        double y0[3] = {y[0],(y[0]+y[1])/2,(y[0]+y[2])/2};
        I += integrateRecursive(x0,y0,f,n-1);

        double x1[3] = {x[1],(x[0]+x[1])/2,(x[1]+x[2])/2};
        double y1[3] = {y[1],(y[0]+y[1])/2,(y[1]+y[2])/2};
        I += integrateRecursive(x1,y1,f,n-1);

        double x2[3] = {x[2],(x[0]+x[2])/2,(x[1]+x[2])/2};
        double y2[3] = {y[2],(y[0]+y[2])/2,(y[1]+y[2])/2};
        I += integrateRecursive(x2,y2,f,n-1);

        double x3[3] = {(x[1]+x[2])/2,(x[0]+x[1])/2,(x[0]+x[2])/2};
        double y3[3] = {(y[1]+y[2])/2,(y[0]+y[1])/2,(y[0]+y[2])/2};
        I += integrateRecursive(x3,y3,f,n-1);
    }   
     
    return I;
}
