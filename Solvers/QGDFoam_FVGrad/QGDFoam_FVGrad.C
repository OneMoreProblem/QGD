/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    QGDFoam_FVGrad

Description
    Quasi-Gasdynamic solver.
    Now only for orthogonal mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOdictionary paramsQGDDict
    (
        IOobject
        (
            "paramsQGD",            // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );
    
    surfaceScalarField alphaField
    (
        IOobject
        (
            "alphaQGD",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // --- Directed interpolation of primitive fields onto faces

//      rho
        surfaceScalarField rhof ("rhof",fvc::interpolate(rho));

//      U
        surfaceVectorField Uf ("Uf", fvc::interpolate(U));
//      p
        surfaceScalarField pf = fvc::interpolate(p);
        surfaceVectorField gradPf ("gradPf",fvc::interpolate(fvc::grad(p)));

        volScalarField rPsi("rPsi", 1.0/psi);

        volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv() * rPsi));

// for QGD
//      c
        surfaceScalarField cf ("cf", fvc::interpolate(c));
//      h
         surfaceScalarField hQGD = 1.0 / mesh.surfaceInterpolation::deltaCoeffs();
//      tau
        surfaceScalarField tauQGDf("tauQGD", alphaField * hQGD / (cf ));
// Fluxes for first equation

        surfaceVectorField rhoW1("rhoW1", tauQGDf * fvc::interpolate(fvc::div((rho * (U * U)))));

        surfaceVectorField rhoW2 = tauQGDf * gradPf;

        surfaceVectorField jm("jm", Uf*rhof - rhoW1 - rhoW2);
        surfaceScalarField phiJm("phiJm", jm & mesh.Sf());
// End fluxes for first equation

// Fluxes for second equation
        surfaceVectorField phiP
        (
          pf*mesh.Sf()
        );

        surfaceVectorField phiUW
        (
          "phiUW",
          tauQGDf * Uf *
          (
            rhof * ((fvc::snGrad(U) & mesh.Sf()/mesh.magSf() ) * (Uf & mesh.Sf() ) )
            +
            ( gradPf & mesh.Sf() )
          )
         );

        surfaceScalarField gammaf ("gammaf",fvc::interpolate(thermo.Cp() / thermo.Cv()));

        surfaceVectorField phiR
        (
          "phiR",
          tauQGDf * Uf * ( gradPf & mesh.Sf() )
          +
          tauQGDf * mesh.Sf() *
          (
            (pf * (fvc::snGrad(U) & mesh.Sf()/mesh.magSf()) * gammaf)
          )
        );

        surfaceTensorField gradUf
        (
          "gradUf",
          fvc::interpolate(fvc::grad((U)))
        );

        surfaceScalarField divUf
        (
          "divUf",
          fvc::interpolate(fvc::div(U))
        );

        tensor unit(1, 0, 0, 0, 1, 0, 0, 0, 1);

        surfaceTensorField Pif
        (
          "Pif",
          tauQGDf * 
          (
//            Uf * (rhof * (Uf & gradUf) + gradPf)
            ((fvc::interpolate(rho*(U*U)) & gradUf) + Uf*gradPf)
            +
            unit * ( (Uf & gradPf) + (gammaf * pf * divUf) )
          )
        );

        surfaceVectorField phiJmU("phiJmU", (jm * Uf) & mesh.Sf());
        surfaceVectorField phiPixx("phiPixx", phiUW + phiR);
        surfaceVectorField phiPi("phiPi", Pif & mesh.Sf());
// End fluxes for second equation

// Fluxes for third equation
        surfaceScalarField ef("ef", fvc::interpolate(e));
        surfaceScalarField rhoEf = fvc::interpolate(rhoE);
        surfaceScalarField Hf("Hf", fvc::interpolate((rhoE + p)/rho));
// term for q
        volScalarField gammam1 = thermo.Cp() / thermo.Cv() - 1.0;
        surfaceScalarField gammam1f = fvc::interpolate(gammam1);

        surfaceVectorField qf
        (
          "qf",
            - tauQGDf * ( ( (fvc::interpolate(rho*(U*U))) & fvc::interpolate(fvc::grad(e))) + (pf * rhof * Uf * (Uf & fvc::interpolate(fvc::grad(1/rho)))))
        );

        surfaceScalarField phiQ
        (
          "phiQ",
          qf & mesh.Sf()
        );

        surfaceScalarField phiPiU
        (
          "phiPiU",
          (Pif & Uf) & mesh.Sf()
        );

        surfaceScalarField phiJmH("phiJmH", phiJm * Hf);
// End for third equation
// ******************************************************************* //

        #include "readTimeControls.H"
        #include "centralCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Solve density
        solve
        (
            fvm::ddt(rho)
            +
            fvc::div(phiJm)
        );

//  *********************** //

        // --- Solve momentum
        solve
        (
            fvm::ddt(rhoU)
            + 
            fvc::div(phiJmU)
            +
            fvc::div(phiP)
            -
            fvc::div(phiPi)
         );

//      Correct velocity
        U.ref() =
            rhoU()
           /rho();
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

// ************************************************************************************** //

       //--- Solve energy
       solve
       (
           fvm::ddt(rhoE)
         + fvc::div(phiJmH)
         + fvc::div(phiQ)
         - fvc::div(phiPiU)
       );
       
       if (runTime.outputTime())
       {
           rhoE.write();
           phiJmH.write();
//	   phiJm.write();
//	   Hf.write();
//	   HC.write();
       }
       
        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        rhoE.boundaryFieldRef() == rho.boundaryField()*
        (e.boundaryField() + 0.5*magSqr(U.boundaryField()));
       
       thermo.correct();

       p.ref() =
           rho()
          /psi();
       p.correctBoundaryConditions();
       rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
// *********************************************************************** //
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	
	//	Info<< "Fields: " << T << rho << p << U << endl; 
	
	Info<< "max/min T = "<< max(T)<< "/" << min(T) << endl;
	Info<< "max/min p = "<< max(p)<< "/" << min(p) << endl;
	Info<< "max/min rho = "<< max(rho)<< "/" << min(rho) << endl;
	Info<< "max/min U = "<< max(U)<< "/" << min(U) << endl;
    }



    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
