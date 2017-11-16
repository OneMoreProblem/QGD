#include "extendedFaceStencil.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <HashTable.H>

#include "emptyFvPatch.H"
#include "coupledFvPatch.H"
#include "wedgeFvPatch.H"
#include "symmetryFvPatch.H"
#include "symmetryPlaneFvPatch.H"

void Foam::extendedFaceStencil::faceScalarGPar(const List3<scalar>& procVfValues, const surfaceScalarField& sF, surfaceVectorField& gradIF)
{

    surfaceScalarField tField = sF; 

    faceScalarDer(procVfValues,sF.boundaryField(),0,tField);
    gradIF.primitiveFieldRef().replace(0, tField);
    faceScalarDer(procVfValues,sF.boundaryField(),1,tField);
    gradIF.primitiveFieldRef().replace(1, tField);
    faceScalarDer(procVfValues,sF.boundaryField(),2,tField);
    gradIF.primitiveFieldRef().replace(2, tField);

};

void Foam::extendedFaceStencil::faceScalarGPar(const tmp<List3<scalar>>& tprocVfValues, const tmp<surfaceScalarField>& tsF, tmp<surfaceVectorField>& tmpgradIF)
{
}

