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

#ifndef Curves_hh_INCLUDED
#define Curves_hh_INCLUDED

// ObjexxFCL Headers
#include <ObjexxFCL/Array1D.hh>
#include <ObjexxFCL/Array1S.hh>
#include <ObjexxFCL/Array2D.hh>
#include <ObjexxFCL/Array2S.hh>
#include <ObjexxFCL/Array5D.hh>
#include <ObjexxFCL/Optional.hh>

// EnergyPlus Headers
#include <EnergyPlus.hh>
#include <DataGlobals.hh>

namespace EnergyPlus {

namespace Curves {

	struct Curve
	{
		enum class Type {
			Linear, BiLinear, Quadratic, BiQuadratic, Cubic, QuadraticLinear, BiCubic, TriQuadratic,
			Exponent, Quartic, FuncPressDrop, MultiVariableLookup, FanPressureRise, ExponentialSkewNormal, Sigmoid,
			RectangularHyperbola1, RectangularHyperbola2, ExponentialDecay, DoubleExponentialDecay, QuadLinear,
			CubicLinear, ChillerPartLoadCustom
		};
		// Members
		std::string name; // Curve Name
		//int ObjectType; // Curve object type (e.g., integer for Curve:Linear above)
		//int CurveType; // Curve type (see parameter definitions above)
		//int InterpolationType; // table interpolation method
		//int DataFormat; // format of tabular data
		//int TableIndex; // Index to tablular data (0 if a standard curve object)
		//int TableVariables; // Number of independent variables (0 if a standard curve object)
		//int NumIVLowErrorIndex; // Index to table object error message for too few IV's
		//int NumIVHighErrorIndex; // Index to table object error message for too many IV's
		//int X1SortOrder; // sort order for table data for X1
		//int X2SortOrder; // sort order for table data for X2
		//Real64 Coeff1; // constant coefficient
		//Real64 Coeff2; // linear coeff (1st independent variable)
		//Real64 Coeff3; // quadratic coeff (1st independent variable)
		//Real64 Coeff4; // linear coeff (2nd ind var) or cubic coeff
		//Real64 Coeff5; // quadratic coeff (2nd independent variable)
		//Real64 Coeff6; // cross coeff (1st & 2nd ind var)
		//Real64 Coeff7; // cubic coeff for bicubic (1st ind var)
		//Real64 Coeff8; // cubic coeff for bicubic (2nd ind var)
		//Real64 Coeff9; // cross coeff for bicubic (1st quadratic & 2nd linear)
		//Real64 Coeff10; // cross coeff for bicubic (1st linear & 2nd quadratic)
		//Real64 Coeff11; // cross coeff
		//Real64 Coeff12; // cross coeff
		//Real64 Var1Max; // maximum of 1st independent variable
		//Real64 Var1Min; // minimum of 1st independent variable
		//Real64 Var2Max; // maximum of 2nd independent variable
		//Real64 Var2Min; // minimum of 2nd independent variable
		//Real64 Var3Max; // maximum of 3rd independent variable
		//Real64 Var3Min; // minimum of 3rd independent variable
		//Real64 Var4Max; // maximum of 4th independent variable
		//Real64 Var4Min; // minimum of 4th independent variable
		//Real64 Var5Max; // maximum of 5th independent variable
		//Real64 Var5Min; // minimum of 5th independent variable
		Real64 curveMin; // minimum value of curve output
		Real64 curveMax; // maximum value of curve output
		bool curveMinPresent; // If TRUE, then cap minimum curve output
		bool curveMaxPresent; // if TRUE, then cap maximum curve output
		//Array1D< TriQuadraticCurveDataStruct > Tri2ndOrder; // structure for triquadratic curve data
		bool EMSOverrideOn; // if TRUE, then EMS is calling to override curve value
		Real64 EMSOverrideCurveValue; // Value of curve result EMS is directing to use
		// report variables
		Real64 CurveOutput; // curve output or result
		Real64 CurveInput1; // curve input #1 (e.g., x or X1 variable)
		Real64 CurveInput2; // curve input #1 (e.g., y or X2 variable)
		Real64 CurveInput3; // curve input #1 (e.g., z or X3 variable)
		Real64 CurveInput4; // curve input #1 (e.g., X4 variable)
		Real64 CurveInput5; // curve input #1 (e.g., X5 variable)

		// Default Constructor
		Curve() :
			curveMin(0.0),
			curveMax(0.0),
			curveMinPresent(false),
			curveMaxPresent(false),
			EMSOverrideOn(false),
			EMSOverrideCurveValue(0.0),
			CurveOutput(0.0),
			CurveInput1(0.0),
			CurveInput2(0.0),
			CurveInput3(0.0),
			CurveInput4(0.0),
			CurveInput5(0.0)
		{}
		virtual ~Curve(){}
		virtual Type type() const = 0;
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const = 0;
		virtual Real64 value(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0 ) const
		{
			Real64 result = compute(v1, v2, v3, v4, v5);
			if (curveMinPresent) result = max(result, curveMin);
			if (curveMaxPresent) result = min(result, curveMax);
			return result;
		}

	};

	struct Linear : public Curve
	{
		Real64 coeff1; // constant coefficient
		Real64 coeff2; // linear coeff (1st independent variable)
		Real64 var1Max; // maximum of 1st independent variable
		Real64 var1Min; // minimum of 1st independent variable
		Linear() :
			Curve(),
			coeff1(0.0),
			coeff2(0.0),
			var1Max(0.0),
			var1Min(0.0)
		{}
		virtual ~Linear(){}
		virtual Type type() const
		{
			return Curve::Type::Linear;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
		{
			return coeff1 + v1 * coeff2;
		}
		virtual Real64 value(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
		{
			v1 = max(min(v1, var1Max), var1Min);
			Real64 result = compute(v1, v2, v3, v4, v5);
			if(curveMinPresent) result = max(result, curveMin);
			if(curveMaxPresent) result = min(result, curveMax);
			return result;
		}
	};

	struct Quadratic : public Linear
	{
		Real64 coeff3; // quadratic coeff (1st independent variable)
		Quadratic() :
			Linear(),
			coeff3(0.0)
		{}
		virtual ~Quadratic();
		virtual Type type() const
		{
			return Curve::Type::Quadratic;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
		{
			return coeff1 + v1 * (coeff2 + v1 * coeff3);
		}
	};

	struct Exponent : public Quadratic
	{
		Exponent() :
			Quadratic()
		{}
		virtual ~Exponent();
		virtual Type type() const
		{
			return Curve::Type::Exponent;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
		{
			return coeff1 + coeff2 * std::pow(v1, coeff3);
		}
	};

	struct RectangularHyperbola1 : public Quadratic
	{
		RectangularHyperbola1() :
			Quadratic()
		{}
		virtual ~RectangularHyperbola1();
		virtual Type type() const
		{
			return Curve::Type::RectangularHyperbola1;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
		{
			Real64 curveValueNumer = coeff1 * v1;
			Real64 curveValueDenom = coeff2 + v1;
			return (curveValueNumer / curveValueDenom) + coeff3;
		}
	};

	struct RectangularHyperbola2 : public Quadratic
	{
		RectangularHyperbola2() :
			Quadratic()
		{}
		virtual ~RectangularHyperbola2();
		virtual Type type() const
		{
			return Curve::Type::RectangularHyperbola2;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
		{
			Real64 curveValueNumer = coeff1 * v1;
			Real64 curveValueDenom = coeff2 + v1;
			return ( curveValueNumer / curveValueDenom ) + ( coeff3 * v1 );
		}
	};

	struct ExponentialDecay : public Quadratic
	{
		ExponentialDecay() :
			Quadratic()
		{}
		virtual ~ExponentialDecay();
		virtual Type type() const
		{
			return Curve::Type::ExponentialDecay;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
		{
			return coeff1 + coeff2 * std::exp(coeff3 * v1);
		}
	};

	struct Cubic : public Quadratic
	{
		Real64 coeff4; // linear coeff (2nd ind var) or cubic coeff
		Cubic() :
			Quadratic(),
			coeff4(0.0)
		{}
		virtual ~Cubic();
		virtual Type type() const
		{
			return Curve::Type::Cubic;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
		{
			return coeff1 + v1 * (coeff2 + v1 * (coeff3 + v1 * coeff4));
		}
	};

	struct FanPressureRise : public Cubic
	{
		FanPressureRise() :
			Cubic()
		{}
		virtual ~FanPressureRise();
		virtual Type type() const
		{
			return Curve::Type::FanPressureRise;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
		{
			return v1 * (coeff1 * v1 + coeff2 + coeff3 * std::sqrt(v2)) + coeff4 * v2;
		}
	};

	struct ExponentialSkewNormal : public Cubic
	{
		ExponentialSkewNormal() :
			Cubic()
		{}
		virtual ~ExponentialSkewNormal();
		virtual Type type() const
		{
			return Curve::Type::ExponentialSkewNormal;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
		{
			static Real64 const sqrt_2_inv(1.0 / std::sqrt(2.0));
			Real64 coeffZ1 = (v1 - coeff1) / coeff2;
			Real64 coeffZ2 = (coeff4 * v1 * std::exp(coeff3 * v1) - coeff1) / coeff2;
			Real64 coeffZ3 = -coeff1 / coeff2;
			//    CurveValueNumer = EXP(-0.5d0 * CoeffZ1**2) * (1.0d0 + SIGN(1.0d0,CoeffZ2) * ErfFunction(ABS(CoeffZ2)/SQRT(2.0d0)))
			//    CurveValueDenom = EXP(-0.5d0 * CoeffZ3**2) * (1.0d0 + SIGN(1.0d0,CoeffZ3) * ErfFunction(ABS(CoeffZ3)/SQRT(2.0d0)))
			Real64 curveValueNumer = std::exp(-0.5 * (coeffZ1 * coeffZ1)) * (1.0 + sign(1.0, coeffZ2) * std::erf(std::abs(coeffZ2) * sqrt_2_inv));
			Real64 curveValueDenom = std::exp(-0.5 * (coeffZ3 * coeffZ3)) * (1.0 + sign(1.0, coeffZ3) * std::erf(std::abs(coeffZ3) * sqrt_2_inv));
			return curveValueNumer / curveValueDenom;
		}
	};

	struct Quartic : public Cubic
	{
		Real64 coeff5; // quadratic coeff (2nd independent variable)
		Quartic() :
			Cubic(),
			coeff5(0.0)
		{}
		virtual ~Quartic();
		virtual Type type() const
		{
			return Curve::Type::Quartic;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
		{
			return coeff1 + v1 * (coeff2 + v1 * (coeff3 + v1 * (coeff4 + v1 * coeff5)));
		}
	};

	struct Sigmoid : public Quartic
	{
		Sigmoid() :
			Quartic()
		{}
		virtual ~Sigmoid();
		virtual Type type() const
		{
			return Curve::Type::Sigmoid;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
		{
			Real64 curveValueExp = std::exp((coeff3 - v1) / coeff4);
			return coeff1 + coeff2 / std::pow(1.0 + curveValueExp, coeff5);
		}
	};

	struct DoubleExponentialDecay : public Quartic
	{
		DoubleExponentialDecay() :
			Quartic()
		{}
		virtual ~DoubleExponentialDecay();
		virtual Type type() const
		{
			return Curve::Type::DoubleExponentialDecay;
		}
		virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
		{
			return coeff1 + coeff2 * std::exp(coeff3 * v1) + coeff4 * std::exp(coeff5 * v1);
		}
	};

} // Curves

} // EnergyPlus

#endif
