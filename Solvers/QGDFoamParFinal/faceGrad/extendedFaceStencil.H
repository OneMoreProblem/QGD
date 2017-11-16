#ifndef extendedFaceStencil_H
#define extendedFaceStencil_H

#include "regIOobject.H"
#include "labelList.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"
#include "vector.H"
#include "List.H"

#include "processorFvPatch.H"

namespace Foam
{


class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class extendedFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedFaceStencil
:
    public regIOobject
{

private:

    template <class T>
    using List2 = List<List<T> >;

    template <class T>
    using List3 = List2<List<T> >;
    
    template <class T>
    using labelHashTable = HashTable<T,label>;
    
    template <class T>
    using Triple = FixedList<T,3>;


    
//    template <class T>
//    using tripleLabelHashTable = HashTable<T,Triple<T> >

private:

        // List of neighbour cells for each face
        labelListList neighbourCells_; 
        
private:

        // Mesh
        const fvMesh& mesh_;
        
        //Data for internal field
        List2<vector> GdfAll_;
        List2<scalar> wf2All_;
        
        //Data for patches
        /* For faces connecting 2 processors  */
        List<label> procPairs_; // if >= 0, connects processor through patch
        List<label> neigProcs_; // No of neighbouring processor
        labelHashTable<label> idProcPatchPairs_; //hash for processors, connected through the patch
        label nProcPatches_;
        
        /* Least square method weights at processor patches */
        List3<vector> procGdf_;
        List3<scalar> procWf2_;
        
        /* Addressing information */
        List3<label> myProcPatchCells_;
        List2<label> ownEnd_;
        List2<label> neiStart_;
        List2<label> neiEnd_;
        List2<label> corStart_;
        List2<label> corEnd_;
        
        /* Global information about domain decomposition */
        //labelHashTable <List<label> > pointProcs_; //list of processors containing point (global id)
        labelHashTable <List<label> > pointCells_; //list of cells containing cell (global id)
        labelHashTable <label> cellProc_;
        
        /* Addressing information for corner neighbouring processors */
        List2<Triple<label> > corAddr_;
        List2<label> corCellIds_;
        labelHashTable<label> corProcIds_;
        
public:

// Constructors

        //- Construct from IOobject. Optional flag for if IOobject is the
        //  top level regIOobject.
        extendedFaceStencil(const IOobject&, const bool isTime = false);

        //- Construct as copy
        extendedFaceStencil(const regIOobject&);
        
        //- Construct as copy, transferring registry registration to copy
        //  if registerCopy is true
        extendedFaceStencil(const regIOobject&, bool registerCopy);

        //- Destructor
        ~extendedFaceStencil();

        void rename(const word& newName);

        bool readData(Istream&);

        bool read();

        bool modified() const;

        bool readIfModified();

        bool writeData(Ostream&) const;

        bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType
        ) const;

        bool write() const;

        void findNeighbours();

        void calculateWeights();
        

        void faceScalarDer(const Field<scalar>& iF,const Field<scalar>& sF,int com, surfaceScalarField& rField);
        void faceScalarDer(const tmp<Field<scalar>>& tiF,const tmp<Field<scalar>>& tsF, int com, tmp<surfaceScalarField>& trField);

        void faceScalarDer
        (
            const List3<scalar>& procVfValues,
            const surfaceScalarField& sF,
            int derComp,
            surfaceScalarField& rField
        );

        void faceScalarDer
        (
            const tmp<List3<scalar>>& tprocVfValues,
            const tmp<surfaceScalarField>& tsF,
            int derComp,
            tmp<surfaceScalarField>& trField
        );

        tmp<surfaceVectorField> faceScalarGrad(const volScalarField& iF);
        tmp<surfaceVectorField> faceScalarGrad(const tmp<volScalarField>& tiF);

        tmp<surfaceVectorField> faceScalarGrad(const volScalarField& iF, const surfaceScalarField& sF);
        tmp<surfaceVectorField> faceScalarGrad(const tmp<volScalarField>& tiF, const tmp<surfaceScalarField>& tsF);

        tmp<surfaceTensorField> faceVectorGrad(const volVectorField& iVF, const surfaceVectorField& sVF);
        tmp<surfaceTensorField> faceVectorGrad(const tmp<volVectorField>& tiVF,const tmp<surfaceVectorField>& tsVF);

        tmp<surfaceScalarField> faceVectorDiv(const volVectorField& iVF);
        tmp<surfaceScalarField> faceVectorDiv(const tmp<volVectorField>& tiVF);

        tmp<surfaceVectorField> faceTensorDiv(const volTensorField& iTF);
        tmp<surfaceVectorField> faceTensorDiv(const tmp<volTensorField>& iTF);
        
        #include "extendedFaceStencilFormVfValues.H"
};

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
