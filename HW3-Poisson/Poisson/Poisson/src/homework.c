#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{   
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges;

    int nBoundary = 0;

    // thEdges->elem contains the nodes of the edges on the boundary
    // We need to identify the nodes that are on the boundary
    // We will do this by counting the number of times each unique node appears in theEdges->elem
    // A node that appears twice should be counted only one time
    // First we create a list of the unique nodes
    int nNodes = theEdges->nodes->nNodes;
    int *uniqueNodes = calloc(nNodes, sizeof(int));

    for (int i = 0; i < theEdges->nElem; i++) {
        for (int j = 0; j < 2; j++) {
            uniqueNodes[theEdges->elem[2*i+j]]++;
        }
    }

    // Now we count the number of nodes that appear at least once
    for (int i = 0; i < nNodes; i++) {
        if (uniqueNodes[i] > 0) {
            nBoundary++;
        }
    }
    
    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary*sizeof(int));
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");
    int k = 0;
    for (int i = 0; i < nBoundary; i++) {
        while (uniqueNodes[k] == 0) {
            k++;
        }
        theBoundary->elem[i] = k;
        k++; 
    }
    free(uniqueNodes);
}
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{

    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    geoMeshFree(theProblem->geo);
    free(theProblem);
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    
    int nLocal = theMesh->nLocalNode;
    for (int j = 0; j < nLocal; j++) {
        int iNode = theMesh->elem[iElem*nLocal+j];
        x[j] = theMesh->nodes->X[iNode];
        y[j] = theMesh->nodes->Y[iNode];
        map[j] = iNode; 
    }

}

# endif
# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem)
{

    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo,"Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
 
    if (theSpace->n > 4) Error("Unexpected discrete space size !");  
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,map[4];
    int nLocal = theMesh->nLocalNode;

    double **A = theSystem->A;
    double *B = theSystem->B;

    // Linear system assembly
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        // Local assembly
        femPoissonLocal(theProblem,iElem,map,x,y);
        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double weight = theRule->weight[iInteg];
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdsi = 0.0;
            double dydsi = 0.0;
            double dxdet = 0.0;
            double dydet = 0.0;

            for(i = 0; i < nLocal; i++) {
                dxdsi += x[i]*dphidxsi[i];
                dydsi += y[i]*dphidxsi[i];
                dxdet += x[i]*dphideta[i];
                dydet += y[i]*dphideta[i];
            }


            double detJac = dxdsi*dydet - dydsi*dxdet;
            if (detJac <= 0.0){ // permutation des noeuds pour les remettre dans le bon ordre
                int node = theMesh->elem[iElem*nLocal];
                theMesh->elem[iElem*nLocal] = theMesh->elem[iElem*nLocal+2];
                theMesh->elem[iElem*nLocal+2] = node;
            };

            for(i = 0; i < nLocal; i++) {
                dphidx[i] = (dydet*dphidxsi[i]-dydsi*dphideta[i])/detJac;
                dphidy[i] = (-dxdet*dphidxsi[i]+dxdsi*dphideta[i])/detJac;
            }

            for(i = 0; i <theSpace->n; i++) {
                B[map[i]] += phi[i]*detJac*weight;
                for(j = 0; j < theSpace->n; j++) {
                    A[map[i]][map[j]] += (dphidx[i]*dphidx[j]+dphidy[i]*dphidy[j])*detJac*weight;
                }
            }
        }
    }

    // Constraints
    for(iEdge = 0; iEdge < theBoundary->nElem; iEdge++) {
        int iNode = theBoundary->elem[iEdge];
        femFullSystemConstrain(theSystem,iNode,0.0);
    }

    // Solve
    femFullSystemEliminate(theSystem);
}

# endif



