// EnergyPlus, Copyright (c) 1996-2016, The Board of Trustees of the University of Illinois and
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights
// reserved.
//
// If you have questions about your rights to use or distribute this software, please contact
// Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without Lawrence Berkeley National Laboratory's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the
// features, functionality or performance of the source code ("Enhancements") to anyone; however,
// if you choose to make your Enhancements available either publicly, or directly to Lawrence
// Berkeley National Laboratory, without imposing a separate written license agreement for such
// Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free
// perpetual license to install, use, modify, prepare derivative works, incorporate into other
// computer software, distribute, and sublicense such enhancements or derivative works thereof,
// in binary and source code form.

// EnergyPlus Headers
#include <EnergyPlus.hh>
#include <General.hh>
#include <DataGlobals.hh>
#include <DataBranchAirLoopPlant.hh>
#include <UtilityRoutines.hh>
#include <Curves.hh>
#include "Eigen/Dense"

namespace EnergyPlus {

namespace Curves {

Real64 calculateMoodyFrictionFactor(Real64 const ReynoldsNumber, Real64 const RoughnessRatio)
{

	// FUNCTION INFORMATION:
	//       AUTHOR         Edwin Lee
	//       DATE WRITTEN   August 2009
	//       MODIFIED       na
	//       RE-ENGINEERED  na

	// PURPOSE OF THIS FUNCTION:
	// This will evaluate the moody friction factor based on Reynolds number and roughness ratio

	// METHODOLOGY EMPLOYED:
	// General empirical correlations for friction factor based on Moody Chart data

	// REFERENCES:
	// Haaland, SE (1983). "Simple and Explicit Formulas for the Friction Factor in Turbulent Flow".
	//   Trans. ASIVIE, J. of Fluids Engineering 103: 89-90.

	// Return value
	Real64 CalculateMoodyFrictionFactor;

	// Locals
	// FUNCTION ARGUMENT DEFINITIONS:

	// FUNCTION PARAMETER DEFINITIONS:
	// na

	// INTERFACE BLOCK SPECIFICATIONS:
	// na

	// DERIVED TYPE DEFINITIONS:
	// na

	// FUNCTION LOCAL VARIABLE DECLARATIONS:
	static bool FrictionFactorErrorHasOccurred(false);

	//Check for no flow before calculating values
	if (ReynoldsNumber == 0.0) {
		return 0.0;
	}

	//Check for no roughness also here
	if (RoughnessRatio == 0.0) {
		return 0.0;
	}

	//Calculate the friction factor
	Real64 Term1 = std::pow(RoughnessRatio / 3.7, 1.11);
	Real64 Term2 = 6.9 / ReynoldsNumber;
	Real64 Term3 = -1.8 * std::log10(Term1 + Term2);
	if (Term3 != 0.0) {
		CalculateMoodyFrictionFactor = std::pow(Term3, -2.0);
	} else {
		if (!FrictionFactorErrorHasOccurred) {
			std::string RR = General::RoundSigDigits(RoughnessRatio, 7);
			std::string Re = General::RoundSigDigits(ReynoldsNumber, 1);
			ShowSevereError("Plant Pressure System: Error in moody friction factor calculation");
			ShowContinueError("Current Conditions: Roughness Ratio=" + RR + "; Reynolds Number=" + Re);
			ShowContinueError("These conditions resulted in an unhandled numeric issue.");
			ShowContinueError("Please contact EnergyPlus support/development team to raise an alert about this issue");
			ShowContinueError("This issue will occur only one time.  The friction factor has been reset to 0.04 for calculations");
			FrictionFactorErrorHasOccurred = true;
		}
		CalculateMoodyFrictionFactor = 0.04;
	}

	return CalculateMoodyFrictionFactor;

}

Real64 PlantPressure::value(Real64 massFlow, Real64 density, Real64 viscosity, Real64, Real64)
{

	// FUNCTION INFORMATION:
	//       AUTHOR         Edwin Lee
	//       DATE WRITTEN   August 2009
	//       MODIFIED       na
	//       RE-ENGINEERED  na

	// PURPOSE OF THIS FUNCTION:
	// This will evaluate the pressure drop for components which use pressure information

	// METHODOLOGY EMPLOYED:
	// Friction factor pressure drop equation:
	// DP = [f*(L/D) + K] * (rho * V^2) / 2

	// REFERENCES:
	// na

	// Using/Aliasing

	// Return value
	Real64 PressureCurveValue;

	// Locals
	// FUNCTION ARGUMENT DEFINITIONS:

	// FUNCTION PARAMETER DEFINITIONS:
	// na

	// INTERFACE BLOCK SPECIFICATIONS:
	// na

	// DERIVED TYPE DEFINITIONS:
	// na

	// FUNCTION LOCAL VARIABLE DECLARATIONS:
	//Real64 Diameter;
	//Real64 MinorLossCoeff;
	//Real64 Length;
	//Real64 Roughness;
	//bool IsConstFPresent;
	//Real64 ConstantF;
	Real64 frictionFactor;
	//Real64 CrossSectArea;
	//Real64 Velocity;
	//Real64 ReynoldsNumber;
	//Real64 RoughnessRatio;

	//Retrieve data from structure
	//Diameter = PressureCurve(PressureCurveIndex).EquivDiameter;
	//MinorLossCoeff = PressureCurve(PressureCurveIndex).MinorLossCoeff;
	//Length = PressureCurve(PressureCurveIndex).EquivLength;
	//Roughness = PressureCurve(PressureCurveIndex).EquivRoughness;
	//IsConstFPresent = PressureCurve(PressureCurveIndex).ConstantFPresent;
	//ConstantF = PressureCurve(PressureCurveIndex).ConstantF;

	//Intermediate calculations
	Real64 crossSectArea = (DataGlobals::Pi / 4.0) * pow_2(equivDiameter);
	Real64 velocity = massFlow / (density * crossSectArea);
	Real64 ReynoldsNumber = density * equivDiameter * velocity / viscosity; //assuming mu here
	Real64 roughnessRatio = equivRoughness / equivDiameter;

	//If we don't have any flow then exit out
	if (massFlow < DataBranchAirLoopPlant::MassFlowTolerance) {
		curveInput1 = massFlow;
		curveInput2 = density;
		curveInput3 = velocity;
		curveOutput = 0.0;
		return 0.0;
	}

	//Calculate the friction factor
	if (constantFPresent) { //use the constant value
		frictionFactor = constantF;
	} else { // must calculate f
		frictionFactor = calculateMoodyFrictionFactor(ReynoldsNumber, roughnessRatio);
	}

	//Pressure drop calculation
	PressureCurveValue = (frictionFactor * (equivLength / equivDiameter) + minorLossCoeff) * (density * pow_2(velocity))*0.5;

	if (EMSOverrideOn) {
		PressureCurveValue = EMSOverrideCurveValue;
	}

	curveInput1 = massFlow;
	curveInput2 = density;
	curveInput3 = velocity;
	curveOutput = PressureCurveValue;

	return PressureCurveValue;
}

static unsigned intervalByBisection(Real64 v, const std::vector<Real64> &x, unsigned i0, unsigned i1)
{
	unsigned delta = i1 - i0;
	if (delta == 1) {
		return i0;
	}
	unsigned mid = (unsigned)0.5*delta;
	if (v < x[mid]) {
		return intervalByBisection(v, x, i0, mid);
	} else if (v == x[mid]) {
		return mid;
	}
	return intervalByBisection(v, x, mid, i1);
}

Real64 Table1D::compute(Real64 v1, Real64 v2, Real64 v3, Real64 v4, Real64 v5) const
{
	unsigned interval = intervalByBisection(v1, x1, 0, x1.size() - 1);
	return 0.0;
}

Polynomial* Polynomial::leastSquaresFit(unsigned order, const std::vector<Real64> &x, const std::vector<double> &y)
{
	if (order < 2) {
		// Error message
		return nullptr;
	}
	unsigned n = std::min(x.size(), y.size());
	if (order < n) {
		// Error message
		return nullptr;
	}
	Eigen::VectorXd b(n);
	Eigen::MatrixXd A(n, order);
	for (unsigned i = 0; i < n; i++) {
		b[i] = y[i];
		A(i, 0) = 1.0;
		for (unsigned j = 1; j < order; j++) {
			A(i, j) = A(i, j - 1)*x[i];
		}
	}
	Eigen::VectorXd c = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	std::vector<Real64> coeffs(order);
	for (unsigned j = 0; j < order; j++) {
		coeffs[j] = c[j];
	}
	Polynomial *poly =  new Polynomial(coeffs);
	poly->m_fromTabularData = true;
	return poly;
}

} // Curves

} // EnergyPlus

