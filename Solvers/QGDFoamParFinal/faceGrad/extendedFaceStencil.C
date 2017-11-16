#include "extendedFaceStencil.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include "wedgeFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "symmetryFvPatch.H"
#include "emptyFvPatch.H"
#include <HashTable.H>

// constructors
Foam::extendedFaceStencil::extendedFaceStencil(const IOobject& io, const bool isTime)
:
    regIOobject(io, isTime),
    mesh_(refCast<const fvMesh>(io.db()))
{
    findNeighbours();
    calculateWeights();
}


Foam::extendedFaceStencil::extendedFaceStencil(const regIOobject& rio)
:
    regIOobject(rio),
    mesh_(refCast<const fvMesh>(rio.db()))
{
    findNeighbours();
    calculateWeights();
}

Foam::extendedFaceStencil::extendedFaceStencil(const regIOobject& rio, bool registerCopy)
:
    regIOobject(rio, registerCopy),
    mesh_(refCast<const fvMesh>(rio.db()))
{
    findNeighbours();
    calculateWeights();
}

Foam::extendedFaceStencil::~extendedFaceStencil()
{
}

void Foam::extendedFaceStencil::rename(const word& newName)
{
};

bool Foam::extendedFaceStencil::readData(Istream&)
{
  return false;
};

bool Foam::extendedFaceStencil::read()
{
  return false;
};

bool Foam::extendedFaceStencil::modified()  const
{
  return false;
};

bool Foam::extendedFaceStencil::readIfModified()
{
  return false;
};

bool Foam::extendedFaceStencil::writeData(Ostream&) const
{
  return false;
};

bool Foam::extendedFaceStencil::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType
)  const
{
  return false;
};

bool Foam::extendedFaceStencil::write()  const
{
  return false;
};

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
tmp<surfaceTensorField> Foam::extendedFaceStencil::faceVectorGrad(const volVectorField& iVF, const surfaceVectorField& sVF)
{
//    surfaceVectorField gradComp0col = faceScalarGrad(iVF.component(0));
//    surfaceVectorField gradComp1col = faceScalarGrad(iVF.component(1));
//    surfaceVectorField gradComp2col = faceScalarGrad(iVF.component(2));
    
    tmp<surfaceTensorField> tgradIVF(0*fvc::snGrad(iVF) * mesh_.Sf() / mesh_.magSf());
    surfaceScalarField tField = sVF.component(0)*0;
    surfaceTensorField& gradIVF = tgradIVF.ref();

    faceScalarDer(iVF.primitiveField().component(0),sVF.primitiveField().component(0),0,tField);
    gradIVF.primitiveFieldRef().replace(0,tField);

    faceScalarDer(iVF.primitiveField().component(0),sVF.primitiveField().component(0),1,tField);
    gradIVF.primitiveFieldRef().replace(3,tField);

    faceScalarDer(iVF.primitiveField().component(0),sVF.primitiveField().component(0),2,tField);
    gradIVF.primitiveFieldRef().replace(6,tField);

    faceScalarDer(iVF.primitiveField().component(1),sVF.primitiveField().component(1),0,tField);
    gradIVF.primitiveFieldRef().replace(1,tField);

    faceScalarDer(iVF.primitiveField().component(1),sVF.primitiveField().component(1),1,tField);
    gradIVF.primitiveFieldRef().replace(4,tField);

    faceScalarDer(iVF.primitiveField().component(1),sVF.primitiveField().component(1),2,tField);
    gradIVF.primitiveFieldRef().replace(7,tField);

    faceScalarDer(iVF.primitiveField().component(2),sVF.primitiveField().component(2),0,tField);
    gradIVF.primitiveFieldRef().replace(2,tField);

    faceScalarDer(iVF.primitiveField().component(2),sVF.primitiveField().component(2),1,tField);
    gradIVF.primitiveFieldRef().replace(5,tField);

    faceScalarDer(iVF.primitiveField().component(2),sVF.primitiveField().component(2),2,tField);
    gradIVF.primitiveFieldRef().replace(8,tField);


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
        )
        {
            notConstrain = false;
        }

        if (notConstrain)
        {
            gradIVF.boundaryFieldRef()[ipatch] = mesh_.Sf().boundaryField()[ipatch]/mesh_.magSf().boundaryField()[ipatch]*iVF.boundaryField()[ipatch].snGrad();
        }
    }



    //set internal field
//    gradIVF.primitiveFieldRef().replace(0, gradComp0col.primitiveField().component(0));
//    gradIVF.primitiveFieldRef().replace(1, gradComp1col.primitiveField().component(0));
//    gradIVF.primitiveFieldRef().replace(2, gradComp2col.primitiveField().component(0));
//    
//    gradIVF.primitiveFieldRef().replace(3, gradComp0col.primitiveField().component(1));
//    gradIVF.primitiveFieldRef().replace(4, gradComp1col.primitiveField().component(1));
//    gradIVF.primitiveFieldRef().replace(5, gradComp2col.primitiveField().component(1));
//    
//    gradIVF.primitiveFieldRef().replace(6, gradComp0col.primitiveField().component(2));
//    gradIVF.primitiveFieldRef().replace(7, gradComp1col.primitiveField().component(2));
//    gradIVF.primitiveFieldRef().replace(8, gradComp2col.primitiveField().component(2));
/*    
    //set external fields
    forAll(mesh_.boundaryMesh(), patchi)
    {
        forAll(mesh_.boundary()[patchi], facei)
        {
            gradIVF.boundaryFieldRef()[patchi][facei].component(0) = 
                gradComp0col.boundaryField()[patchi][facei].component(0);
            gradIVF.boundaryFieldRef()[patchi][facei].component(1) = 
                gradComp1col.boundaryField()[patchi][facei].component(0);
            gradIVF.boundaryFieldRef()[patchi][facei].component(2) = 
                gradComp2col.boundaryField()[patchi][facei].component(0);

            gradIVF.boundaryFieldRef()[patchi][facei].component(3) = 
                gradComp0col.boundaryField()[patchi][facei].component(1);
            gradIVF.boundaryFieldRef()[patchi][facei].component(4) = 
               gradComp1col.boundaryField()[patchi][facei].component(1);
            gradIVF.boundaryFieldRef()[patchi][facei].component(5) = 
                    gradComp2col.boundaryField()[patchi][facei].component(1);

            gradIVF.boundaryFieldRef()[patchi][facei].component(6) = 
                gradComp0col.boundaryField()[patchi][facei].component(2);
            gradIVF.boundaryFieldRef()[patchi][facei].component(7) = 
                gradComp1col.boundaryField()[patchi][facei].component(2);
            gradIVF.boundaryFieldRef()[patchi][facei].component(8) = 
                gradComp2col.boundaryField()[patchi][facei].component(2);
        }
    }
*/
    if(!Pstream::parRun())
    {
        return tgradIVF;
    }

    List<List3<scalar>> procVfValues(nProcPatches_);

    formVfValues(iVF,procVfValues);

    faceScalarDer(procVfValues[0],sVF.component(0),0,tField);
    gradIVF.boundaryFieldRef().replace(0,tField.boundaryField());

    faceScalarDer(procVfValues[0],sVF.component(0),1,tField);
    gradIVF.boundaryFieldRef().replace(3,tField.boundaryField());

    faceScalarDer(procVfValues[0],sVF.component(0),2,tField);
    gradIVF.boundaryFieldRef().replace(6,tField.boundaryField());

    faceScalarDer(procVfValues[1],sVF.component(1),0,tField);
    gradIVF.boundaryFieldRef().replace(1,tField.boundaryField());

    faceScalarDer(procVfValues[1],sVF.component(1),1,tField);
    gradIVF.boundaryFieldRef().replace(4,tField.boundaryField());

    faceScalarDer(procVfValues[1],sVF.component(1),2,tField);
    gradIVF.boundaryFieldRef().replace(7,tField.boundaryField());

    faceScalarDer(procVfValues[2],sVF.component(2),0,tField);
    gradIVF.boundaryFieldRef().replace(2,tField.boundaryField());

    faceScalarDer(procVfValues[2],sVF.component(2),1,tField);
    gradIVF.boundaryFieldRef().replace(5,tField.boundaryField());

    faceScalarDer(procVfValues[2],sVF.component(2),2,tField);
    gradIVF.boundaryFieldRef().replace(8,tField.boundaryField());

    return tgradIVF;
    
};

tmp<surfaceTensorField> Foam::extendedFaceStencil::faceVectorGrad(const tmp<volVectorField>& tiVF, const tmp<surfaceVectorField>& tsVF)
{
    return faceVectorGrad(tiVF(),tsVF());
}


//- Calculate divergence of volume vector field on the faces.
//
// \param iVF        Internal vector field.
//                   Allowable values: constant reference to the volVectorField.
//
// \return           Divergence of iVF (scalar field) which was computed on the faces of mesh.
tmp<surfaceScalarField> Foam::extendedFaceStencil::faceVectorDiv(const volVectorField& iVF)
{

    surfaceVectorField sVF = linearInterpolate(iVF);
    surfaceScalarField tField = sVF.component(0)*0;
    tmp<surfaceScalarField> tdivIVF(0*fvc::snGrad(iVF) & mesh_.Sf() / mesh_.magSf());
    surfaceScalarField& divIVF = tdivIVF.ref();

    faceScalarDer(iVF.primitiveField().component(0),sVF.primitiveField().component(0),0,tField);
    divIVF.primitiveFieldRef() = tField;

    faceScalarDer(iVF.primitiveField().component(1),sVF.primitiveField().component(1),1,tField);
    divIVF.primitiveFieldRef() += tField;

    faceScalarDer(iVF.primitiveField().component(2),sVF.primitiveField().component(2),2,tField);
    divIVF.primitiveFieldRef() += tField;



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
        )
        {
            notConstrain = false;
        }

        if (notConstrain)
        {
            divIVF.boundaryFieldRef()[ipatch] = mesh_.Sf().boundaryField()[ipatch]/mesh_.magSf().boundaryField()[ipatch]&iVF.boundaryField()[ipatch].snGrad();
        }
    }

   
/*
    surfaceVectorField gradComp0 = faceScalarGrad(iVF.component(0));
    surfaceVectorField gradComp1 = faceScalarGrad(iVF.component(1));
    surfaceVectorField gradComp2 = faceScalarGrad(iVF.component(2));

    tmp<surfaceScalarField> tdivIVF(0*fvc::snGrad(iVF) & mesh_.Sf() / mesh_.magSf());
    surfaceScalarField& divIVF = tdivIVF.ref();
    
    divIVF.primitiveFieldRef() = gradComp0.primitiveField().component(0)
                               + gradComp1.primitiveField().component(1)
                               + gradComp2.primitiveField().component(2);
    
    forAll(mesh_.boundary(), patchi)
    {
        divIVF.boundaryFieldRef()[patchi] = 
            gradComp0.boundaryField()[patchi].component(0)
            +
            gradComp1.boundaryField()[patchi].component(1)
            +
            gradComp2.boundaryField()[patchi].component(2);
    }
*/


    if(!Pstream::parRun())
    {
        return tdivIVF;
    }

    List<List3<scalar>> procVfValues(nProcPatches_); 
    formVfValues(iVF,procVfValues);

    
    faceScalarDer(procVfValues[0],sVF.component(0),0,tField);
    divIVF.boundaryFieldRef() = tField.boundaryField();

    faceScalarDer(procVfValues[1],sVF.component(1),1,tField);
    divIVF.boundaryFieldRef() += tField.boundaryField();

    faceScalarDer(procVfValues[2],sVF.component(2),2,tField);
    divIVF.boundaryFieldRef() += tField.boundaryField();

    
    return tdivIVF;
};

tmp<surfaceScalarField> Foam::extendedFaceStencil::faceVectorDiv(const tmp<volVectorField>& tiVF)
{
    return faceVectorDiv(tiVF());
}

//- Calculate divergence of volume tensor field on the faces.
//
// \param iTF        Internal tensor field.
//                   Allowable values: constant reference to the volTensorField.
//
// \return           Divergence of iTF (vector field) which was computed on the faces of mesh.
tmp<surfaceVectorField> Foam::extendedFaceStencil::faceTensorDiv(const volTensorField& iTF)
{

    surfaceTensorField sTF = linearInterpolate(iTF);
    surfaceScalarField tField = sTF.component(0)*0;
    tmp<surfaceVectorField> tdivITF(0*fvc::snGrad(iTF.component(0)) * mesh_.Sf() / mesh_.magSf());
    surfaceVectorField& divITF = tdivITF.ref();
    surfaceScalarField divComp = tField;

    faceScalarDer(iTF.primitiveField().component(0),sTF.primitiveField().component(0),0,tField);
    divComp = tField;
    faceScalarDer(iTF.primitiveField().component(1),sTF.primitiveField().component(1),1,tField);
    divComp += tField;
    faceScalarDer(iTF.primitiveField().component(2),sTF.primitiveField().component(2),2,tField);
    divComp += tField;
    divITF.primitiveFieldRef().replace(0,divComp);

    faceScalarDer(iTF.primitiveField().component(3),sTF.primitiveField().component(3),0,tField);
    divComp = tField;
    faceScalarDer(iTF.primitiveField().component(4),sTF.primitiveField().component(4),1,tField);
    divComp += tField;
    faceScalarDer(iTF.primitiveField().component(5),sTF.primitiveField().component(5),2,tField);
    divComp += tField;
    divITF.primitiveFieldRef().replace(1,divComp);

    faceScalarDer(iTF.primitiveField().component(6),sTF.primitiveField().component(6),0,tField);
    divComp = tField;
    faceScalarDer(iTF.primitiveField().component(7),sTF.primitiveField().component(7),1,tField);
    divComp += tField;
    faceScalarDer(iTF.primitiveField().component(8),sTF.primitiveField().component(8),2,tField);
    divComp += tField;
    divITF.primitiveFieldRef().replace(2,divComp);

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
        )
        {
            notConstrain = false;
        }

        if (notConstrain)
        {
            divITF.boundaryFieldRef()[ipatch] = mesh_.Sf().boundaryField()[ipatch]/mesh_.magSf().boundaryField()[ipatch]&iTF.boundaryField()[ipatch].snGrad();
        }
    }


    if (!Pstream::parRun()) 
    {
        return tdivITF;
    }


    List<List3<scalar>> procVfValues(nProcPatches_);
    formVfValues(iTF,procVfValues);

    faceScalarDer(procVfValues[0],sTF.component(0),0,tField);
    divComp = tField;
    faceScalarDer(procVfValues[1],sTF.component(1),1,tField);
    divComp += tField;
    faceScalarDer(procVfValues[2],sTF.component(2),2,tField);
    divComp += tField;
    divITF.boundaryFieldRef().replace(0,divComp.boundaryField());

    faceScalarDer(procVfValues[3],sTF.component(3),0,tField);
    divComp = tField;
    faceScalarDer(procVfValues[4],sTF.component(4),1,tField);
    divComp += tField;
    faceScalarDer(procVfValues[5],sTF.component(5),2,tField);
    divComp += tField;
    divITF.boundaryFieldRef().replace(1,divComp.boundaryField());

    faceScalarDer(procVfValues[6],sTF.component(6),0,tField);
    divComp = tField;
    faceScalarDer(procVfValues[7],sTF.component(7),1,tField);
    divComp += tField;
    faceScalarDer(procVfValues[8],sTF.component(8),2,tField);
    divComp += tField;
    divITF.boundaryFieldRef().replace(2,divComp.boundaryField());



/*    
    tmp<surfaceVectorField> gradComp0 (faceScalarGrad(iTF.component(0)));
    tmp<surfaceVectorField> gradComp1 (faceScalarGrad(iTF.component(1)));
    tmp<surfaceVectorField> gradComp2 (faceScalarGrad(iTF.component(2)));
    
    tmp<surfaceVectorField> gradComp3 (faceScalarGrad(iTF.component(3)));
    tmp<surfaceVectorField> gradComp4 (faceScalarGrad(iTF.component(4)));
    tmp<surfaceVectorField> gradComp5 (faceScalarGrad(iTF.component(5)));

    tmp<surfaceVectorField> gradComp6 (faceScalarGrad(iTF.component(6)));
    tmp<surfaceVectorField> gradComp7 (faceScalarGrad(iTF.component(7)));
    tmp<surfaceVectorField> gradComp8 (faceScalarGrad(iTF.component(8)));

    tmp<surfaceScalarField> divComp0 (gradComp0().component(0) + gradComp3().component(1) + gradComp6().component(2));
    tmp<surfaceScalarField> divComp1 (gradComp1().component(0) + gradComp4().component(1) + gradComp7().component(2));
    tmp<surfaceScalarField> divComp2 (gradComp2().component(0) + gradComp5().component(1) + gradComp8().component(2));

    tmp<surfaceVectorField> tdivITF(0*fvc::snGrad(iTF.component(0)) * mesh_.Sf() / mesh_.magSf());
    surfaceVectorField& divITF = tdivITF.ref();
    
    divITF.primitiveFieldRef().replace(0, divComp0().primitiveField());
    divITF.primitiveFieldRef().replace(1, divComp1().primitiveField());
    divITF.primitiveFieldRef().replace(2, divComp2().primitiveField());
    
    forAll(mesh_.boundary(), patchi)
    {
        forAll(mesh_.boundary()[patchi], facei)
        {
            divITF.boundaryFieldRef()[patchi][facei].component(0) = 
                divComp0().boundaryField()[patchi][facei];
            divITF.boundaryFieldRef()[patchi][facei].component(1) = 
                divComp1().boundaryField()[patchi][facei];
            divITF.boundaryFieldRef()[patchi][facei].component(2) = 
                divComp2().boundaryField()[patchi][facei];
        }
    }
*/    
    return tdivITF;
}

tmp<surfaceVectorField> Foam::extendedFaceStencil::faceTensorDiv(const tmp<volTensorField>& tiTF)
{
    return faceTensorDiv(tiTF());
}

//
//END-OF-FILE
//


