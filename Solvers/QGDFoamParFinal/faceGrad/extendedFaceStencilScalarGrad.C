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

//- Calculate gradient of volume scalar function on the faces
//
// \param iF         Internal scalar field.
//                   Allowable values: constant reference to the volScalarField.
//
// \return           Gradient of iF (vector field) which was computed on the faces of mesh.
tmp<surfaceVectorField> Foam::extendedFaceStencil::faceScalarGrad(const volScalarField& iF)
{
 //   Info << "inF = " << iF << endl;
    surfaceScalarField sF = linearInterpolate(iF);   
    
    tmp<surfaceVectorField> tgradIF(0*fvc::snGrad(iF)  * mesh_.Sf() / mesh_.magSf());
    surfaceVectorField& gradIF = tgradIF.ref();
    
    gradIF = faceScalarGrad(iF,sF);

    return tgradIF;
};
     
tmp<surfaceVectorField> Foam::extendedFaceStencil::faceScalarGrad(const tmp<volScalarField>& tiF)
{
    return faceScalarGrad(tiF());
}



tmp<surfaceVectorField> Foam::extendedFaceStencil::faceScalarGrad(const volScalarField& iF, const surfaceScalarField& sF)
{

    tmp<surfaceVectorField> tgradIF(0*fvc::snGrad(iF)  * mesh_.Sf() / mesh_.magSf());
    surfaceVectorField& gradIF = tgradIF.ref();
    //scalarField tField = sF;
    surfaceScalarField tField = sF*0;

//    const faceList& faces = mesh_.faces();

    faceScalarDer(iF.primitiveField(),sF.primitiveField(),0,tField);
    gradIF.primitiveFieldRef().replace(0, tField);
    faceScalarDer(iF.primitiveField(),sF.primitiveField(),1,tField);
    gradIF.primitiveFieldRef().replace(1, tField);
    faceScalarDer(iF.primitiveField(),sF.primitiveField(),2,tField);
    gradIF.primitiveFieldRef().replace(2, tField);

    
    //update boundary field
    forAll(mesh_.boundaryMesh(), ipatch)
    {
        bool notConstrain = true;
        const fvPatch& fvp = mesh_.boundary()[ipatch];
        if
        (
            isA<emptyFvPatch>(fvp) ||
            isA<wedgeFvPatch>(fvp) ||
            isA<coupledFvPatch>(fvp) ||
            isA<symmetryFvPatch>(fvp) ||
            isA<symmetryPlaneFvPatch>(fvp)
//            fvp.coupled()
        )
        {
            notConstrain = false;
        }

        if (notConstrain)
        {
            gradIF.boundaryFieldRef()[ipatch] = iF.boundaryField()[ipatch].snGrad()*mesh_.Sf().boundaryField()[ipatch] / mesh_.magSf().boundaryField()[ipatch];
        }
    }

    if(!Pstream::parRun())
    {
        return tgradIF;
    }
    
    /*
     *
     * Update processor patches for parallel case
     *
     */
    //allocate storage for near-patch field


    List3<scalar> procVfValues(nProcPatches_); //array of values from neighb. processors
    formVfValues(iF,procVfValues);

//    Pout << "procVfValues = " << procVfValues << endl;
//    Pout << "procGdf_ = " << procGdf_ << endl;
//    Pout << "procWf2_ = " << procWf2_ << endl;

    //Step 3. Calculate gradient at faces on processor patches
//    Info << "iF = "   << iF           << endl;
//    Info << "PVfV = " << procVfValues << endl;
    faceScalarDer(procVfValues,sF,0,tField);
    gradIF.boundaryFieldRef().replace(0, tField.boundaryFieldRef());
    faceScalarDer(procVfValues,sF,1,tField);
    gradIF.boundaryFieldRef().replace(1, tField.boundaryFieldRef());
    faceScalarDer(procVfValues,sF,2,tField);
    gradIF.boundaryFieldRef().replace(2, tField.boundaryFieldRef());

//    Info <<"PgradIF = " << gradIF << endl;

    return tgradIF;


};

tmp<surfaceVectorField> Foam::extendedFaceStencil::faceScalarGrad(const tmp<volScalarField>& tiF, const tmp<surfaceScalarField>& tsf)
{
    return faceScalarGrad(tiF());
}

//
//END-OF-FILE
//


