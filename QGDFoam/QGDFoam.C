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
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

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
        surfaceScalarField rho_int ("rho_int",fvc::interpolate(rho));

//      U
        surfaceVectorField U_int = fvc::interpolate(U);

//      rho * U
        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));
        surfaceVectorField rhoU_int = fvc::interpolate(rhoU);

//      p
        surfaceScalarField p_int = fvc::interpolate(p);
        volVectorField gradP("gradP", fvc::grad(p));
        surfaceVectorField gradP_int ("gradPint", fvc::interpolate(gradP));

// from rhoCentralFOAM
        volScalarField rPsi("rPsi", 1.0/psi);
        surfaceScalarField rPsi_int = fvc::interpolate(rPsi);
        surfaceScalarField psi_int = fvc::interpolate(psi);

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

//      flux throught the faces --- ???
        surfaceScalarField phiv_int("phiv_int", U_int & mesh.Sf());

        volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name())*mesh.magSf()
        );
        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name())*mesh.magSf()
        );
        
        

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
// end from rhoCentralFOAM

//      for QGD
//      c
        surfaceScalarField c_int ("c_int", fvc::interpolate(c));
        c.write();
        c_int.write();
//      h
        surfaceScalarField hQGD = 1.0 / mesh.surfaceInterpolation::deltaCoeffs();
//      tau
        surfaceScalarField tauQGD_int("tauQGD", 0.5 * hQGD / c_int);
//      end for QGD

// ******************************************************************* //

        #include "centralCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

// 1. Solve the density equation:
//   \[
//      \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho U + \rho W),\\
//      where \rho W = \tau \nabla \cdot (\rho U \times U + I p) = \tau (\nabla \cdot (\rho U \times U) + \grad p),
//            \tau = \alpha \frac{h}{c}, \quad \alpha = 0.5.
//   \]
//
// - Calculate the mass flux throught faces:
//   \[
//      (\rho U \times U) \cdot S_f.
//   \]
// - Find the divergence of mass flux and interpolate it to face centers.
// - Find the pressure gradient field.
// - Find the $\tau$ (tau_QGD).
// - Find the $\rho W$ as sum of the divergence of mass flux and pressure gradient.
// - Multiply $\rho W$ to $S_f$ (scalar product) and to $\tau$.
// - Solve density equation and find density.

        phi = phiv_int*rho_int;

        surfaceVectorField phiU = U_int * phi;

        volVectorField divPhiU = fvc::div(phiU);
        surfaceVectorField divPhiU_int = fvc::interpolate(divPhiU);

        surfaceScalarField rhoW1 = divPhiU_int & mesh.Sf();
        rhoW1 *= tauQGD_int;

        surfaceScalarField rhoW2 = gradP_int & mesh.Sf();
        rhoW2 *= tauQGD_int;

        surfaceScalarField rhoW("rhoW", rhoW1 + rhoW2);
        
        surfaceVectorField rhoWW("rhoWW", (divPhiU_int + gradP_int) * tauQGD_int);
        

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi)
        - fvc::div(rhoW)
        );

rhoWW.write();
tauQGD_int.write();
gradP_int.write();
rho_int.write();

//  *********************** //

        surfaceVectorField phiUp
        (
            phiv_int*rhoU_int
          + p_int*mesh.Sf()
        );

        surfaceVectorField phiW
        (
            phiv_int * (divPhiU_int + gradP_int) * tauQGD_int
        );

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));
        solve(fvm::ddt(rhoU) + fvc::div(phiUp) 
         -fvc::div(phiW)
         );

//      Correct velocity
//        U.dimensionedInternalField() =
//            rhoU.dimensionedInternalField()
//           /rho.dimensionedInternalField();
//        U.correctBoundaryConditions();
//        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
// ********************** //

//        surfaceScalarField sigmaDotU
//        (
//            "sigmaDotU",
//            (
//                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
//              + (mesh.Sf() & fvc::interpolate(tauMC))
//            )
//            & (a_pos*U_pos + a_neg*U_neg)
//        );

       surfaceScalarField e_int = fvc::interpolate(e);

       surfaceScalarField phiEp
       (
           "phiEp",
           phiv_int*(rho_int*(e_int + 0.5*magSqr(U_int)) + p_int)
       );

     surfaceScalarField phiWE
     (
         "phiWE",
         rhoW*(e_int + 0.5*magSqr(U_int) + (p_int / rho_int))
     );

/*       volScalarField gammam1 = thermo.Cp() / thermo.Cv() - 1.0;
       surfaceScalarField gammam1_int = fvc::interpolate(gammam1);

       surfaceScalarField Rq
       (
           "Rq",
           (
             ( U_int
             &
             fvc::interpolate(fvc::grad(p/rho))
             )
             /
             gammam1_int
             + 
             ( p_int * ( U_int & fvc::interpolate(fvc::grad(1.0/rho)) ) )
           )
           * tauQGD_int * rho_int
       );

       surfaceScalarField phiURq = phiv_int * Rq;

       surfaceScalarField kappaQGD = tauQGD_int * p_int /gammam1_int;
       
       surfaceVectorField qNS = fvc::interpolate(fvc::grad(p/rho));
       */

//       solve
//       (
//           fvm::ddt(rhoE)
//         + fvc::div(phiEp)
//         - fvc::div(phiWE)
//       );

//       e.correctBoundaryConditions();
       
       thermo.correct();
       
//       rhoE.boundaryField() =
//           rho.boundaryField()*
//           (
//               e.boundaryField() + 0.5*magSqr(U.boundaryField())
//           );

//       p.dimensionedInternalField() =
//           rho.dimensionedInternalField()
//          /psi.dimensionedInternalField();
//       p.correctBoundaryConditions();
//       rho.boundaryField() = psi.boundaryField()*p.boundaryField();

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
