

#include"fem.h"


#ifndef NORENUMBER 

double *GLarray;

int femIntCompare(const void *a, const void *b)
{
    return (GLarray[*(int*)a] - GLarray[*(int*)b]);
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;
    int * map = malloc(sizeof(int)*theMesh->nodes->nNodes);
    
    for(i = 0; i < theMesh->nodes->nNodes; i++) 
        map[i] = i;

    switch (renumType) {
        case FEM_NO :
            break;        
        case FEM_XNUM : 
            GLarray = theMesh->nodes->X;
            qsort(map, theMesh->nodes->nNodes, sizeof(int), femIntCompare);
        case FEM_YNUM : 
            GLarray = theMesh->nodes->Y;
            qsort(map, theMesh->nodes->nNodes, sizeof(int), femIntCompare);
            break;
        default : Error("Unexpected renumbering option");}

    for(i = 0; i < theMesh->nodes->nNodes; i++) 
        theMesh->nodes->number[map[i]] = i;
    
    free(map);
}


#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    int nLocalNode = theMesh->nLocalNode;
    int map[4]; int myMin, myMax, myBand = 0;
    for(int i = 0; i < theMesh->nElem; ++i) {
        int *elem = &(theMesh->elem[nLocalNode*i]);
        for(int j = 0; j < nLocalNode; j++) {
            map[j] = theMesh->nodes->number[elem[j]];
        }
        myMax = map[0];
        myMin = map[0];
        for(int j = 1; j < nLocalNode; ++j) {
            myMax = (myMax < map[j]) ? map[j] : myMax;
            myMin = (myMin > map[j]) ? map[j] : myMin;
        }
        myBand = (myBand < myMax - myMin) ? myMax - myMin : myBand;
    }
    return(++myBand);
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    for(int i = 0; i < nLoc; i++) {
        int iGlo = map[i];
        myBandSystem->B[iGlo] += Bloc[i];
        for(int j = 0; j < nLoc; j++) {
            int jGlo = map[j];
            if(jGlo >= iGlo) myBandSystem->A[iGlo][jGlo] += Aloc[i*nLoc+j];
        }
    }
}


#endif
#ifndef NOBANDELIMINATE


double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    for (k = 0; k < size; k++) {
        if(fabs(A[k][k]) <= 1e-8) Error("femBandSystemEliminate : null pivot");
        jend = k+band < size ? k+band : size;
        for (i = k+1; i < jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i; j < jend; j++) 
                A[i][j] -= factor * A[k][j];
            B[i] -= factor * B[k];
        }
    }

    // backward substitution
    for(int i = size-1; i >= 0; i--) {
        double sum = 0.0;
        jend = i+band < size ? i+band : size;
        for(int j = i+1; j < jend; j++) 
            sum += A[i][j] * B[j];
        B[i] = (B[i] - sum) / A[i][i];
    }


    return(myBand->B);
}


#endif

