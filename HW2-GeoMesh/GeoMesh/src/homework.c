#include "fem.h"
#include <math.h>


double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;


    double dNotch = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - r0;
    double dHole = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1)) - r1;



    if (dNotch < 0) h = h0;
    else if (dHole < 0) h = h1;  
    else if (dNotch < d0){
        double x = dNotch;
        double x0 = 0;
        double x1 = d0;
        double y0 = h0;
        double y1 = h;
        double y0_prime = 0;
        double y1_prime = 0;

        double t = (x - x0) / (x1 - x0);
        double h00 = (1 + 2 * t) * (1 - t) * (1 - t);
        double h10 = t * (1 - t) * (1 - t);
        double h01 = t * t * (3 - 2 * t);
        double h11 = t * t * (t - 1);
        h =  h00 * y0 + h10 * (x1 - x0) * y0_prime + h01 * y1 + h11 * (x1 - x0) * y1_prime;
    }
    else if (dHole < d1){
        double x = dHole;
        double x0 = 0;
        double x1 = d1;
        double y0 = h1;
        double y1 = h;
        double y0_prime = 0;
        double y1_prime = 0;

        double t = (x - x0) / (x1 - x0);
        double h00 = (1 + 2 * t) * (1 - t) * (1 - t);
        double h10 = t * (1 - t) * (1 - t);
        double h01 = t * t * (3 - 2 * t);
        double h11 = t * t * (t - 1);
        h =  h00 * y0 + h10 * (x1 - x0) * y0_prime + h01 * y1 + h11 * (x1 - x0) * y1_prime;
    }

    return h;
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(x0, y0, 0, w, h, -1, 0,&ierr);   
    ErrorGmsh(ierr);
    int idNotch = gmshModelOccAddDisk(x0, y0, 0, r0, r0, -1,NULL,0,NULL,0,&ierr); 
    ErrorGmsh(ierr);
    int idHole  = gmshModelOccAddDisk(x1, y1, 0, r1, r1, -1,NULL,0,NULL,0,&ierr);    
    ErrorGmsh(ierr);
    
    int plate[] = {2,idPlate};
    int notch[] = {2,idNotch};
    int hole[]  = {2,idHole};
    gmshModelOccCut(plate,2,notch,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate,2,hole,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);
 
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);
//
    
}