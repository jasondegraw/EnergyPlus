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
		Exponent, Quartic, PlantPressure, MultiVariableLookup, FanPressureRise, ExponentialSkewNormal, Sigmoid,
		RectangularHyperbola1, RectangularHyperbola2, ExponentialDecay, DoubleExponentialDecay, QuadLinear,
		CubicLinear, ChillerPartLoadWithLift, Table1D, Polynomial
	};
	enum class Interpolation {
		Linear, Quadratic, Cubic, Quartic
	};
	// Members
	std::string name; // curve Name
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
	bool curveMinPresent; // if TRUE, then cap minimum curve output
	bool curveMaxPresent; // if TRUE, then cap maximum curve output
	//Array1D< TriQuadraticCurveDataStruct > Tri2ndOrder; // structure for triquadratic curve data
	bool EMSOverrideOn; // if TRUE, then EMS is calling to override curve value
	Real64 EMSOverrideCurveValue; // value of curve result EMS is directing to use
	// report variables
	Real64 curveOutput; // curve output or result
	Real64 curveInput1; // curve input #1 (e.g., x or X1 variable)
	Real64 curveInput2; // curve input #1 (e.g., y or X2 variable)
	Real64 curveInput3; // curve input #1 (e.g., z or X3 variable)
	Real64 curveInput4; // curve input #1 (e.g., X4 variable)
	Real64 curveInput5; // curve input #1 (e.g., X5 variable)

	// Default Constructor
	Curve() :
		curveMin(0.0),
		curveMax(0.0),
		curveMinPresent(false),
		curveMaxPresent(false),
		EMSOverrideOn(false),
		EMSOverrideCurveValue(0.0),
		curveOutput(0.0),
		curveInput1(0.0),
		curveInput2(0.0),
		curveInput3(0.0),
		curveInput4(0.0),
		curveInput5(0.0)
	{}
	Curve(const std::string &name) :
		name(name),
		curveMin(0.0),
		curveMax(0.0),
		curveMinPresent(false),
		curveMaxPresent(false),
		EMSOverrideOn(false),
		EMSOverrideCurveValue(0.0),
		curveOutput(0.0),
		curveInput1(0.0),
		curveInput2(0.0),
		curveInput3(0.0),
		curveInput4(0.0),
		curveInput5(0.0)
	{}
	virtual ~Curve(){}
	virtual Type type() const = 0;
	virtual Real64 value(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) = 0;

};

struct Curve1D : public Curve
{
	Real64 var1Max; // maximum of 1st independent variable
	Real64 var1Min; // minimum of 1st independent variable
	Curve1D() :
		Curve(),
		var1Max(0.0),
		var1Min(0.0)
	{}
	virtual ~Curve1D(){}
	virtual Type type() const = 0;
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const = 0;
	virtual Real64 value(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
	{
		curveInput1 = v1;
		v1 = max(min(v1, var1Max), var1Min);
		Real64 result = compute(v1, v2, v3, v4, v5);
		if (curveMinPresent) {
			result = max(result, curveMin);
		}
		if (curveMaxPresent) {
			result = min(result, curveMax);
		}
		curveOutput = result;
		return result;
	}
};

struct Polynomial : public Curve1D
{
	std::vector<Real64> coeffs; // all the coefficients
	Polynomial(const std::vector<Real64> &coefficients = std::vector<Real64>()) :
		Curve1D(),
		coeffs(coefficients)
	{}
	virtual ~Polynomial(){}
	virtual Type type() const
	{
		switch (coeffs.size()) {
		case 2:
			return Type::Linear;
		case 3:
			return Type::Quadratic;
		case 4:
			return Type::Cubic;
		case 5:
			return Type::Quartic;
		default:
			return Type::Polynomial;
		}
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		Real64 val = 1.0;
		Real64 result = 0.0;
		for (auto coeff : coeffs) {
			result += val*coeff;
			val *= v1;
		}
		return result;
	}
};

struct Exponent : public Curve1D
{
	Real64 coeff1; // constant coefficient
	Real64 coeff2; // linear coeff (1st independent variable)
	Real64 coeff3; // quadratic coeff (1st independent variable)
	Exponent() :
		Curve1D(),
		coeff1(0.0),
		coeff2(0.0),
		coeff3(0.0)
	{}
	virtual ~Exponent(){}
	virtual Type type() const
	{
		return Curve::Type::Exponent;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return coeff1 + coeff2 * std::pow(v1, coeff3);
	}
};

struct RectangularHyperbola1 : public Exponent
{
	RectangularHyperbola1() :
		Exponent()
	{}
	virtual ~RectangularHyperbola1(){}
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

struct RectangularHyperbola2 : public Exponent
{
	RectangularHyperbola2() :
		Exponent()
	{}
	virtual ~RectangularHyperbola2(){}
	virtual Type type() const
	{
		return Curve::Type::RectangularHyperbola2;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		Real64 curveValueNumer = coeff1 * v1;
		Real64 curveValueDenom = coeff2 + v1;
		return ( curveValueNumer / curveValueDenom ) + ( coeff3 * v1 );
	}
};

struct ExponentialDecay : public Exponent
{
	ExponentialDecay() :
		Exponent()
	{}
	virtual ~ExponentialDecay(){}
	virtual Type type() const
	{
		return Curve::Type::ExponentialDecay;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return coeff1 + coeff2 * std::exp(coeff3 * v1);
	}
};

struct ExponentialSkewNormal : public Exponent
{
	Real64 coeff4; // linear coeff (2nd ind var) or cubic coeff
	ExponentialSkewNormal() :
		Exponent(),
		coeff4(0.0)
	{}
	virtual ~ExponentialSkewNormal(){}
	virtual Type type() const
	{
		return Curve::Type::ExponentialSkewNormal;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
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

struct Sigmoid : public ExponentialSkewNormal
{
	Real64 coeff5; // quadratic coeff (2nd independent variable)
	Sigmoid() :
		ExponentialSkewNormal(),
		coeff5(0.0)
	{}
	virtual ~Sigmoid(){}
	virtual Type type() const
	{
		return Curve::Type::Sigmoid;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		Real64 curveValueExp = std::exp((coeff3 - v1) / coeff4);
		return coeff1 + coeff2 / std::pow(1.0 + curveValueExp, coeff5);
	}
};

struct DoubleExponentialDecay : public Sigmoid
{
	DoubleExponentialDecay() :
		Sigmoid()
	{}
	virtual ~DoubleExponentialDecay(){}
	virtual Type type() const
	{
		return Curve::Type::DoubleExponentialDecay;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return coeff1 + coeff2 * std::exp(coeff3 * v1) + coeff4 * std::exp(coeff5 * v1);
	}
};

struct Curve2D : public Curve1D
{
	Real64 var2Max; // maximum of 2nd independent variable
	Real64 var2Min; // minimum of 2nd independent variable
	Curve2D() :
		Curve1D(),
		var2Max(0.0),
		var2Min(0.0)
	{}
	virtual ~Curve2D(){}
	virtual Type type() const = 0;
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const = 0;
	virtual Real64 value(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
	{
		curveInput1 = v1;
		curveInput2 = v2;
		v1 = max(min(v1, var1Max), var1Min);
		v2 = max(min(v2, var2Max), var2Min);
		Real64 result = compute(v1, v2, v3, v4, v5);
		if (curveMinPresent) {
			result = max(result, curveMin);
		}
		if (curveMaxPresent) {
			result = min(result, curveMax);
		}
		curveOutput = result;
		return result;
	}
};

struct FanPressureRise : public Curve2D
{
	Real64 coeff1; // constant coefficient
	Real64 coeff2; // linear coeff (1st independent variable)
	Real64 coeff3; // quadratic coeff (1st independent variable)
	Real64 coeff4; // linear coeff (2nd ind var) or cubic coeff
	FanPressureRise() :
		Curve2D(),
		coeff1(0.0),
		coeff2(0.0),
		coeff3(0.0),
		coeff4(0.0)
	{}
	virtual ~FanPressureRise(){}
	virtual Type type() const
	{
		return Curve::Type::FanPressureRise;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return v1 * (coeff1 * v1 + coeff2 + coeff3 * std::sqrt(v2)) + coeff4 * v2;
	}
};

struct BiQuadratic : public Curve2D
{
	Real64 coeff1; // constant coefficient
	Real64 coeff2; // linear coeff (1st independent variable)
	Real64 coeff3; // quadratic coeff (1st independent variable)
	Real64 coeff4; // linear coeff (2nd ind var) or cubic coeff
	Real64 coeff5; // quadratic coeff (2nd independent variable)
	Real64 coeff6; // cross coeff (1st & 2nd ind var)
	BiQuadratic() :
		Curve2D(),
		coeff1(0.0),
		coeff2(0.0),
		coeff3(0.0),
		coeff4(0.0),
		coeff5(0.0),
		coeff6(0.0)
	{}
	virtual ~BiQuadratic(){}
	virtual Type type() const
	{
		return Curve::Type::BiQuadratic;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return coeff1 + v1 * (coeff2 + v1 * coeff3) + v2 * (coeff4 + v2 * coeff5) + v1 * v2 * coeff6;
	}
};

struct BiCubic : public BiQuadratic
{
	Real64 coeff7; // cubic coeff for bicubic (1st ind var)
	Real64 coeff8; // cubic coeff for bicubic (2nd ind var)
	Real64 coeff9; // cross coeff for bicubic (1st quadratic & 2nd linear)
	Real64 coeff10; // cross coeff for bicubic (1st linear & 2nd quadratic)
	BiCubic() :
		BiQuadratic(),
		coeff7(0.0),
		coeff9(0.0),
		coeff10(0.0)
	{}
	virtual ~BiCubic(){}
	virtual Type type() const
	{
		return Curve::Type::BiCubic;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return coeff1 + v1 * coeff2 + v1 * v1 * coeff3 + v2 * coeff4 + v2 * v2 * coeff5 + v1 * v2 * coeff6
			+ v1 * v1 * v1 * coeff7 + v2 * v2 * v2 * coeff8 + v1 * v1 * v2 * coeff9 + v1 * v2 * v2 * coeff10;
	}
};

struct QuadraticLinear : public Curve2D
{
	Real64 coeff1; // constant coefficient
	Real64 coeff2; // linear coeff (1st independent variable)
	Real64 coeff3; // quadratic coeff (1st independent variable)
	Real64 coeff4; // linear coeff (2nd ind var) or cubic coeff
	Real64 coeff5; // quadratic coeff (2nd independent variable)
	Real64 coeff6; // cross coeff (1st & 2nd ind var)
	QuadraticLinear() :
		Curve2D(),
		coeff1(0.0),
		coeff2(0.0),
		coeff3(0.0),
		coeff4(0.0),
		coeff5(0.0),
		coeff6(0.0)
	{}
	virtual ~QuadraticLinear(){}
	virtual Type type() const
	{
		return Curve::Type::QuadraticLinear;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return (coeff1 + v1 * (coeff2 + v1 * coeff3)) + (coeff4 + v1 * (coeff5 + v1 * coeff6)) * v2;
	}
};

struct CubicLinear : public QuadraticLinear
{
	CubicLinear() :
		QuadraticLinear()
	{}
	virtual ~CubicLinear(){}
	virtual Type type() const
	{
		return Curve::Type::CubicLinear;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return (coeff1 + v1 * (coeff2 + v1 * (coeff3 + v1 * coeff4))) + (coeff5 + v1 * coeff6) * v2;
	}
};

struct Curve3D : public Curve2D
{
	Real64 var3Max; // maximum of 3rd independent variable
	Real64 var3Min; // minimum of 3rd independent variable
	Curve3D() :
		Curve2D(),
		var3Max(0.0),
		var3Min(0.0)
	{}
	virtual ~Curve3D(){}
	virtual Type type() const = 0;
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const = 0;
	virtual Real64 value(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
	{
		curveInput1 = v1;
		curveInput2 = v2;
		curveInput3 = v3;
		v1 = max(min(v1, var1Max), var1Min);
		v2 = max(min(v2, var2Max), var2Min);
		v3 = max(min(v3, var3Max), var3Min);
		Real64 result = compute(v1, v2, v3, v4, v5);
		if (curveMinPresent) {
			result = max(result, curveMin);
		}
		if (curveMaxPresent) {
			result = min(result, curveMax);
		}
		curveOutput = result;
		return result;
	}
};

struct TriQuadratic : public Curve3D
{
	Real64 coeffA0;
	Real64 coeffA1;
	Real64 coeffA2;
	Real64 coeffA3;
	Real64 coeffA4;
	Real64 coeffA5;
	Real64 coeffA6;
	Real64 coeffA7;
	Real64 coeffA8;
	Real64 coeffA9;
	Real64 coeffA10;
	Real64 coeffA11;
	Real64 coeffA12;
	Real64 coeffA13;
	Real64 coeffA14;
	Real64 coeffA15;
	Real64 coeffA16;
	Real64 coeffA17;
	Real64 coeffA18;
	Real64 coeffA19;
	Real64 coeffA20;
	Real64 coeffA21;
	Real64 coeffA22;
	Real64 coeffA23;
	Real64 coeffA24;
	Real64 coeffA25;
	Real64 coeffA26;
	TriQuadratic() :
		Curve3D(),
		coeffA0(0.0),
		coeffA1(0.0),
		coeffA2(0.0),
		coeffA3(0.0),
		coeffA4(0.0),
		coeffA5(0.0),
		coeffA6(0.0),
		coeffA7(0.0),
		coeffA8(0.0),
		coeffA9(0.0),
		coeffA10(0.0),
		coeffA11(0.0),
		coeffA12(0.0),
		coeffA13(0.0),
		coeffA14(0.0),
		coeffA15(0.0),
		coeffA16(0.0),
		coeffA17(0.0),
		coeffA18(0.0),
		coeffA19(0.0),
		coeffA20(0.0),
		coeffA21(0.0),
		coeffA22(0.0),
		coeffA23(0.0),
		coeffA24(0.0),
		coeffA25(0.0),
		coeffA26(0.0)
	{}
	virtual ~TriQuadratic(){}
	virtual Type type() const
	{
		return Curve::Type::TriQuadratic;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		Real64 const v1s(v1 * v1);
		Real64 const v2s(v2 * v2);
		Real64 const v3s(v3 * v3);
		return coeffA0 + coeffA1 * v1s + coeffA2 * v1 + coeffA3 * v2s + coeffA4 * v2 + coeffA5 * v3s + coeffA6 * v3
			+ coeffA7 * v1s * v2s + coeffA8 * v1 * v2 + coeffA9 * v1 * v2s + coeffA10 * v1s * v2 + coeffA11 * v1s * v3s
			+ coeffA12 * v1 * v3 + coeffA13 * v1 * v3s + coeffA14 * v1s * v3 + coeffA15 * v2s * v3s + coeffA16 * v2 * v3
			+ coeffA17 * v2 * v3s + coeffA18 * v2s * v3 + coeffA19 * v1s * v2s * v3s + coeffA20 * v1s * v2s * v3
			+ coeffA21 * v1s * v2 * v3s + coeffA22 * v1 * v2s * v3s + coeffA23 * v1s * v2 * v3 + coeffA24 * v1 * v2s * v3
			+ coeffA25 * v1 * v2 * v3s + coeffA26 * v1 * v2 * v3;
	}
};

struct ChillerPartLoadWithLift : public Curve3D
{
	Real64 coeff1;
	Real64 coeff2;
	Real64 coeff3;
	Real64 coeff4;
	Real64 coeff5;
	Real64 coeff6;
	Real64 coeff7;
	Real64 coeff8;
	Real64 coeff9;
	Real64 coeff10;
	Real64 coeff11;
	Real64 coeff12;
	ChillerPartLoadWithLift() :
		Curve3D(),
		coeff1(0.0),
		coeff2(0.0),
		coeff3(0.0),
		coeff4(0.0),
		coeff5(0.0),
		coeff6(0.0),
		coeff7(0.0),
		coeff8(0.0),
		coeff9(0.0),
		coeff10(0.0),
		coeff11(0.0),
		coeff12(0.0)
	{}
	virtual ~ChillerPartLoadWithLift(){}
	virtual Type type() const
	{
		return Curve::Type::ChillerPartLoadWithLift;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		Real64 const v1s(v1 * v1);
		Real64 const v2s(v2 * v2);
		Real64 const v3s(v3 * v3);
		return coeff1 + coeff2*v1 + coeff3*v1*v1 + coeff4*v2 + coeff5*v2*v2 + coeff6*v1*v2 + coeff7*v1*v1*v1
			+ coeff8*v2*v2*v2 + coeff9*v1*v1*v2 + coeff10*v1*v2*v2 + coeff11*v1*v1*v2*v2 + coeff12*v3*v2*v2*v2;
	}
};

struct Curve4D : public Curve3D
{
	Real64 var4Max; // maximum of 4th independent variable
	Real64 var4Min; // minimum of 4th independent variable
	Curve4D() :
		Curve3D(),
		var4Max(0.0),
		var4Min(0.0)
	{}
	virtual ~Curve4D(){}
	virtual Type type() const = 0;
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const = 0;
	virtual Real64 value(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0)
	{
		curveInput1 = v1;
		curveInput2 = v2;
		curveInput3 = v3;
		curveInput4 = v4;
		v1 = max(min(v1, var1Max), var1Min);
		v2 = max(min(v2, var2Max), var2Min);
		v3 = max(min(v3, var3Max), var3Min);
		v4 = max(min(v4, var4Max), var4Min);
		Real64 result = compute(v1, v2, v3, v4, v5);
		if (curveMinPresent) {
			result = max(result, curveMin);
		}
		if (curveMaxPresent) {
			result = min(result, curveMax);
		}
		curveOutput = result;
		return result;
	}
};

struct QuadLinear : public Curve4D
{
	Real64 coeff1; // constant coefficient
	Real64 coeff2; // linear coeff (1st independent variable)
	Real64 coeff3; // quadratic coeff (1st independent variable)
	Real64 coeff4; // linear coeff (2nd ind var) or cubic coeff
	Real64 coeff5; // quadratic coeff (2nd independent variable)
	QuadLinear() :
		Curve4D(),
		coeff1(0.0),
		coeff2(0.0),
		coeff3(0.0),
		coeff4(0.0),
		coeff5(0.0)
	{}
	virtual ~QuadLinear(){}
	virtual Type type() const
	{
		return Curve::Type::QuadLinear;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const
	{
		return coeff1 + v1 * coeff2 + v2 * coeff3 + v3 * coeff4 + v4 * coeff5;
	}
};

struct PlantPressure : public Curve
{
	// Members
	Real64 equivDiameter; // - An effective diameter for calculation of Re & e/D [m]
	Real64 minorLossCoeff; // - K factor                                          [-]
	Real64 equivLength; // - An effective length to apply friction calculation [m]
	Real64 equivRoughness; // - An effective roughness (e) to calculate e/D       [m]
	bool constantFPresent; // - Signal for if a constant value of f was entered
	Real64 constantF; // - Constant value of f (if applicable)               [-]
	//Real64 CurveInput1; // - MassFlow                                         [kg/s]
	//Real64 CurveInput2; // - Density                                          [kg/m3]
	//Real64 CurveInput3; // - Velocity                                         [m/s]

	// Default Constructor
	PlantPressure() :
		Curve(),
		equivDiameter(0.0),
		minorLossCoeff(0.0),
		equivLength(0.0),
		equivRoughness(0.0),
		constantFPresent(false),
		constantF(0.0)
	{}
	virtual ~PlantPressure(){}
	virtual Type type() const
	{
		return Curve::Type::PlantPressure;
	}
	virtual Real64 value(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0);
};

struct Table1D : public Curve1D
{
	std::vector<Real64> x1; // independent variable
	std::vector<Real64> y; // dependent variable
	const Interpolation interpolation; // type of interpolation
	Table1D(Interpolation interp = Interpolation::Linear) :
		Curve1D(),
		interpolation(interp)
	{}
	virtual ~Table1D(){}
	virtual Type type() const
	{
		return Curve::Type::Table1D;
	}
	virtual Real64 compute(Real64 v1, Real64 v2 = 0, Real64 v3 = 0, Real64 v4 = 0, Real64 v5 = 0) const;

};

} // Curves

} // EnergyPlus

#endif
