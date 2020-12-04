/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "incompressibleInterPhaseTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mixture>
Foam::incompressibleInterPhaseTransportModel<Mixture>::
incompressibleInterPhaseTransportModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& rhoPhi,
    const Mixture& mixture
)
:
    varRho_(false),
    phi_(phi),
    rhoPhi_(rhoPhi)
{
    {
        IOdictionary turbulenceProperties
        (
            IOobject
            (
                turbulenceModel::propertiesName,
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        const bool varRho
        (
            turbulenceProperties.getOrDefault<bool>("varRho", false)
        );

        varRho_ = varRho;
    }

    if (varRho_)
    {
        rhoIncTurbulence_ =
        (
            incompressible::phaseIncompressibleTurbulenceModel::New
            (
                rho,
                U,
                rhoPhi,
                phi,
                mixture
            )
        );
    }
    else
    {
        incTurbulence_ = incompressible::turbulenceModel::New
        (
            U,
            phi,
            mixture
        );

        incTurbulence_->validate();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Mixture>
Foam::tmp<Foam::fvVectorMatrix>
Foam::incompressibleInterPhaseTransportModel<Mixture>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    if (varRho_)
    {
       return rhoIncTurbulence_->divDevRhoReff(U);
    }
    else
    {
        return incTurbulence_->divDevRhoReff(rho, U);
    }
}

template<class Mixture>
void Foam::incompressibleInterPhaseTransportModel<Mixture>::correct()
{
    if (varRho_)
    {
        rhoIncTurbulence_->correct();
    }
    else
    {
        incTurbulence_->correct();
    }
}


// ************************************************************************* //
