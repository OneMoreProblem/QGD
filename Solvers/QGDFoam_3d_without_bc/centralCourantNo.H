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

Global
    centralCourantNo

Description
    Calculates the mean and maximum wave speed based Courant Numbers.

\*---------------------------------------------------------------------------*/

if(runTime.controlDict().lookupOrDefault<bool>("adjustTimeStep", false))
{
//    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
//    surfaceScalarField dTField ("dTField", hQGD / (mag(Uf) + mag(cf)));
//
//    Info << "maxCo = " << maxCo << endl;
//
//    runTime.setDeltaT(maxCo * min(dTField.internalField()));
//
//    CoNum = 0.3;//max(amaxSfbyDelta/mesh.magSf()).value()*runTime.deltaTValue();
//
//    //meanCoNum = 0.3;

    surfaceScalarField Cof
    (
        "Cof",
        (mag(Uf) + cf)*runTime.deltaT()/hQGD
    );
    
    CoNum = max(Cof).value();
}

Info << "Courant Number = " << CoNum << endl;
Info << "Delta T: " << runTime.deltaTValue() << endl;

// ************************************************************************* //
