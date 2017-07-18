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
    QGDFoam

Description
    Quasi-Gasdynamic solver.
    Now only for orthogonal mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // --- Directed interpolation of primitive fields onto faces

//      rho
        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));
        surfaceScalarField rhof ("rhof",fvc::interpolate(rho));

//      U
        surfaceVectorField Uf ("Uf", fvc::interpolate(U));

//      rho * U
        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));
        surfaceVectorField rhoUf ("rhoUf",fvc::interpolate(rhoU));

//      p
        surfaceScalarField pf = fvc::interpolate(p);
        surfaceVectorField gradPf ("gradPf", fvc::snGrad(p) * mesh.Sf() / mesh.magSf());

        volScalarField rPsi("rPsi", 1.0/psi);

        surfaceScalarField phivf("phivf", Uf & mesh.Sf());

	volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv() * rPsi));

//        c.write();

//      for QGD
//      c
        surfaceScalarField cf ("cf", fvc::interpolate(c));
//      h
        surfaceScalarField hQGD = 1.0 / mesh.surfaceInterpolation::deltaCoeffs();
//      tau
        surfaceScalarField tauQGDf("tauQGD", 0.5 * hQGD / (cf ));
//        tauQGDf.write();
//      end for QGD

// For first equation


        phi = phivf*rhof;

        surfaceScalarField phiRhoW1 = tauQGDf * (fvc::snGrad(U & U * rho)*mesh.Sf()/mesh.magSf()) & mesh.Sf();

        surfaceScalarField phiRhoW2 = fvc::snGrad(p) * mesh.magSf();
        phiRhoW2 *= tauQGDf;

        surfaceScalarField phiRhoW("phiRhoW", phiRhoW1 + phiRhoW2);

        surfaceScalarField phiJm("phiJm", phi - phiRhoW);
//        phi.write();
// End for first equation

// For second equation
        surfaceVectorField phiRhoUU
        (
          "phiRhoUU",
        // What is correct???
        //phivf*rhoUf                           // first var
          phivf*fvc::interpolate(rho)*Uf        // second var
        );

        surfaceVectorField phiP
        (
          pf*mesh.Sf()
        );


        surfaceVectorField phiRhoUW
        (
    	    "phiRhoUW",
            phivf * (
            ( fvc::snGrad(U & U * rho)*mesh.Sf()/mesh.magSf() )
//            ( fvc::snGrad(U) & Uf * rhof *mesh.Sf()/mesh.magSf() ) it was wrong !!!
            + gradPf
            ) * tauQGDf
        );

// Find Pixx*Sf= (uwPi + rPi)*Sf=phiUW + phiR

        surfaceVectorField phiUW
        (
          "phiUW",
          tauQGDf * Uf *
          (
            rhof * ((fvc::snGrad(U) & mesh.Sf()/mesh.magSf() )* (Uf & mesh.Sf() ) )
            +
            ( gradPf & mesh.Sf() )
          )
        );

        surfaceScalarField gammaf ("gammaf",fvc::interpolate(thermo.Cp() / thermo.Cv()));
//        gammaf.write();

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
        
        surfaceVectorField phiJmU("phiJmU", phiRhoUU - phiRhoUW);
        surfaceVectorField phiPixx("phiPixx", phiUW + phiR);
//        phiJmU.write();
//        phiPixx.write();
//        tauQGDf.write();
// End for second equation

// For third equation
       surfaceScalarField ef("ef", fvc::interpolate(e));
       
       surfaceScalarField rhoEf = fvc::interpolate(rhoE);
//       rhoE.write();
//       p.write();
//       rho.write();
//       Hf.write();
       surfaceScalarField Hf("Hf", fvc::interpolate((rhoE + p)/rho));
       surfaceScalarField Hf1("Hf1", 0.5*magSqr(Uf) + ef + pf/rhof);

       surfaceScalarField phiEp
       (
           "phiEp",
           ( Uf & mesh.Sf() ) * (rhof*Hf)
       );

        surfaceVectorField divPhiUf = fvc::snGrad(U & U * rho) * mesh.Sf() / mesh.magSf();

        phiRhoW1 = divPhiUf & mesh.Sf() * tauQGDf;

        phiRhoW2 = gradPf & mesh.Sf() * tauQGDf;

        phiRhoW = phiRhoW1 + phiRhoW2;

     surfaceScalarField phiWE
     (
         "phiWE",
         phiRhoW*(Hf)
     );

// term for q
       volScalarField gammam1 = thermo.Cp() / thermo.Cv() - 1.0;
       surfaceScalarField gammam1f = fvc::interpolate(gammam1);

       surfaceScalarField qNS
       (
         "qNS",
         - fvc::interpolate(thermo.Cp()) * (tauQGDf * pf * fvc::snGrad(T) * mesh.magSf())
       );

       surfaceScalarField phiQ
       (
         "phiQ",
         - (tauQGDf * rhof * (Uf & Uf) * ( fvc::snGrad(p/rho)*mesh.Sf()/mesh.magSf() & mesh.Sf() ) )/gammam1f
         - (pf * (Uf &  Uf) * ((fvc::snGrad(1.0/rho)*mesh.Sf()/mesh.magSf() ) & mesh.Sf() ) * tauQGDf * rhof)
//         - fvc::interpolate(thermo.Cp()) * (tauQGDf * pf * fvc::snGrad(T) * mesh.magSf())
       );
       
//       qNS.write();
//       phiQ.write();

       surfaceScalarField Pixx
       (
         "Pixx",
         tauQGDf *
         (
           rhof * (Uf & (fvc::snGrad(U)*mesh.Sf()/mesh.magSf()) & Uf)
           +
           ( Uf & gradPf )
         )
         +
         tauQGDf *
         (
           (Uf & gradPf)
         )
       );

       surfaceVectorField ortX = mesh.Sf()/mesh.magSf();
       forAll(ortX, faceI)
       {
          ortX[faceI]=vector(1,0,0);
       }

       surfaceScalarField phiPixxU
       (
         "phiPixxU",
         Pixx * (Uf & mesh.Sf())
         +
         (pf * tauQGDf * gammaf * fvc::snGrad(U) & Uf * mesh.magSf())
//          ( pf * fvc::interpolate(fvc::div(Uf & mesh.Sf())) * gammaf )
       );

//       surfaceScalarField phiJmH("phiJmH", phiEp - phiWE);
	surfaceScalarField phiJmH("phiJmH", phiJm * Hf);
//	phiJm.write();
//	Hf.write();
// End for third equation
// ******************************************************************* //

        #include "centralCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (runTime.outputTime())
        {
    	    psi.write();
    	    volScalarField cpV ("cp", thermo.Cp());
    	    cpV.write();
    	    volScalarField cvV ("cv", thermo.Cv());
            cvV.write();
        }

//        surfaceVectorField phiU = Uf * phi;
//
//        volVectorField divPhiU = fvc::div(phiU);
//        surfaceVectorField divPhiUf = fvc::interpolate(divPhiU);
//
//        surfaceVectorField rhoW("rhoW", (divPhiUf + gradPf) * tauQGDf);
        

        // --- Solve density
        solve(
        fvm::ddt(rho)
        + fvc::div(phiJm)
        );

//  *********************** //


//	phiRhoUU.write();
//	phiRhoUW.write();
//	rhoU.write();
//	phiJmU.write();
//	phiP.write();
//	phiPixx.write();

        // --- Solve momentum
        solve(fvm::ddt(rhoU)
         + fvc::div(phiJmU)
         + fvc::div(phiP)
         - fvc::div(phiPixx)
         );

//	rhoU.write();

//      Correct velocity
        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

// ************************************************************************************** //

//	phiJmH.write();
//	phiQ.write();
//	phiPixxU.write();

       //--- Solve energy
       solve
       (
           fvm::ddt(rhoE)
         + fvc::div(phiJmH)
         + fvc::div(phiQ)
         - fvc::div(phiPixxU)
       );

//       rhoE.write();
       e = rhoE/rho - 0.5*magSqr(U);

       e.correctBoundaryConditions();

       thermo.correct();

//       rhoE.boundaryField() =
//           rho.boundaryField()*
//           (
//               e.boundaryField() + 0.5*magSqr(U.boundaryField())
//           );

       p.dimensionedInternalField() =
           rho.dimensionedInternalField()
          /psi.dimensionedInternalField();
       p.correctBoundaryConditions();
       rho.boundaryField() = psi.boundaryField()*p.boundaryField();

//       turbulence->correct();


// *********************************************************************** //
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
