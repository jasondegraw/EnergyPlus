// EnergyPlus Headers
#include <DataAirflowNetwork.hh>
#include <DataPrecisionGlobals.hh>

#include <ObjexxFCL/gio.hh>
#include <Psychrometrics.hh>
#include <AirflowNetworkSolver.hh>

namespace EnergyPlus {

namespace DataAirflowNetwork {

// MODULE INFORMATION:
//       AUTHOR         Lixing Gu, Don Shirey, and Muthusamy V. Swami
//       DATE WRITTEN   Aug. 2003
//       MODIFIED       na
//       RE-ENGINEERED  na

// PURPOSE OF THIS MODULE:
// This module should contain the information that is needed to simulate
// performance of air distribution system, including pressure, temperature
// and moisture levels at each node, and airflow and sensible and latent energy losses
// at each element

// Using/Aliasing
using namespace DataPrecisionGlobals;

// Data
// module should be available to other modules and routines.  Thus,
// all variables in this module must be PUBLIC.

// MODULE PARAMETER DEFINITIONS:
int const CompTypeNum_DOP(1); // Detailed large opening component
int const CompTypeNum_SOP(2); // Simple opening component
int const CompTypeNum_SCR(3); // Surface crack component
int const CompTypeNum_SEL(4); // Surface effective leakage ratio component
int const CompTypeNum_PLR(5); // Distribution system crack component
int const CompTypeNum_DWC(6); // Distribution system duct component
int const CompTypeNum_CVF(7); // Distribution system constant volume fan component
int const CompTypeNum_FAN(8); // Distribution system detailed fan component
int const CompTypeNum_MRR(9); // Distribution system multiple curve fit power law resistant flow component
int const CompTypeNum_DMP(10); // Distribution system damper component
int const CompTypeNum_ELR(11); // Distribution system effective leakage ratio component
int const CompTypeNum_CPD(12); // Distribution system constant pressure drop component
int const CompTypeNum_COI(13); // Distribution system coil component
int const CompTypeNum_TMU(14); // Distribution system terminal unit component
int const CompTypeNum_EXF(15); // Zone exhaust fan
int const CompTypeNum_HEX(16); // Distribution system heat exchanger
int const CompTypeNum_HOP(17); // Horizontal opening component
int const CompTypeNum_RVD(18); // Reheat VAV terminal damper

// EPlus component Type
int const EPlusTypeNum_SCN(1); // Supply connection
int const EPlusTypeNum_RCN(2); // Return connection
int const EPlusTypeNum_RHT(3); // Reheat terminal
int const EPlusTypeNum_FAN(4); // Fan
int const EPlusTypeNum_COI(5); // Heating or cooling coil
int const EPlusTypeNum_HEX(6); // Heat ecxchanger
int const EPlusTypeNum_RVD(7); // Reheat VAV terminal damper

// EPlus node type
int const EPlusTypeNum_ZIN(1); // Zone inlet node
int const EPlusTypeNum_ZOU(2); // Zone outlet node
int const EPlusTypeNum_SPL(3); // Splitter node
int const EPlusTypeNum_MIX(4); // Mixer node
int const EPlusTypeNum_OAN(5); // Outside air system node
int const EPlusTypeNum_EXT(6); // OA system inlet node
int const EPlusTypeNum_FIN(7); // Fan Inlet node
int const EPlusTypeNum_FOU(8); // Fan Outlet Node
int const EPlusTypeNum_COU(9); // Coil Outlet Node
int const EPlusTypeNum_HXO(10); // Heat exchanger Outlet Node
int const EPlusTypeNum_DIN(11); // Damper Inlet node
int const EPlusTypeNum_DOU(12); // Damper Outlet Node
int const EPlusTypeNum_SPI(13); // Splitter inlet Node
int const EPlusTypeNum_SPO(14); // Splitter Outlet Node

int const iWPCCntr_Input(1);
int const iWPCCntr_SurfAvg(2);

// DERIVED TYPE DEFINITIONS:

// MODULE VARIABLE DECLARATIONS:
// Node simulation variable in air distribution system
// Link simulation variable in air distribution system
// Sensible and latent exchange variable in air distribution system

int SimulateAirflowNetwork(1);
// Vent Control  DistSys Control  Flag    Description
//  NONE           NONE           0      No AirflowNetwork and SIMPLE
//  SIMPLE         NONE           1      Simple calculations only
//  MULTIZONE      NONE           2      Perform multizone calculations only
//  NONE           DISTSYS        3      Perform distribution system durin system on time only
//  SIMPLE         DISTSYS        4      Perform distribution system durin system on time and simple calculations during off time
//  MULTIZONE      DISTSYS        5      Perform distribution system durin system on time and multizone calculations during off time

int const AirflowNetworkControlSimple(1); // Simple calculations only
int const AirflowNetworkControlMultizone(2); // Perform multizone calculations only
int const AirflowNetworkControlSimpleADS(4); // Perform distribution system durin system
// on time and simple calculations during off time
int const AirflowNetworkControlMultiADS(5); // Perform distribution system durin system on time
// and multizone calculations during off time

FArray1D_bool AirflowNetworkZoneFlag;

int NumOfNodesMultiZone(0); // Number of nodes for multizone calculation
int NumOfNodesDistribution(0); // Number of nodes for distribution system calculation
int NumOfLinksMultiZone(0); // Number of links for multizone calculation
int NumOfLinksDistribution(0); // Number of links for distribution system calculation

int AirflowNetworkNumOfNodes(0); // Number of nodes for AirflowNetwork calculation
// = NumOfNodesMultiZone+NumOfNodesDistribution
int AirflowNetworkNumOfComps(0); // Number of components for AirflowNetwork calculation
int AirflowNetworkNumOfLinks(0); // Number of links for AirflowNetwork calculation
// = NumOfLinksMultiZone+NumOfLinksDistribution
// RoomAirManager use
int AirflowNetworkNumOfSurfaces(0); // The number of surfaces for multizone calculation
int AirflowNetworkNumOfZones(0); // The number of zones for multizone calculation

bool RollBackFlag(false); // Roll back flag when system time steo down shifting
FArray1D< Real64 > ANZT; // Local zone air temperature for roll back use
FArray1D< Real64 > ANZW; // Local zone air humidity ratio for roll back use
FArray1D< Real64 > ANCO; // Local zone air CO2 for roll back use
FArray1D< Real64 > ANGC; // Local zone air generic contaminant for roll back use
int AirflowNetworkNumOfExhFan(0); // Number of zone exhaust fans
FArray1D_bool AirflowNetworkZoneExhaustFan; // Logical to use zone exhaust fans
bool AirflowNetworkFanActivated(false); // Supply fan activation flag
bool AirflowNetworkUnitarySystem(false); // set to TRUE for unitary systems (to make answers equal, will remove eventually)
// Multispeed HP only
int MultiSpeedHPIndicator(0); // Indicator for multispeed heat pump use
// Addiitonal airflow needed for an VAV fan to compensate the leakage losses and supply pathway pressure losses [kg/s]
Real64 VAVTerminalRatio(0.0); // The terminal flow ratio when a supply VAV fan reach its max flow rate
bool VAVSystem(false); // This flag is used to represent a VAV system

//     NOTICE
//     Copyright © 1996-2014 The Board of Trustees of the University of Illinois
//     and The Regents of the University of California through Ernest Orlando Lawrence
//     Berkeley National Laboratory.  All rights reserved.
//     Portions of the EnergyPlus software package have been developed and copyrighted
//     by other individuals, companies and institutions.  These portions have been
//     incorporated into the EnergyPlus software package under license.   For a complete
//     list of contributors, see "Notice" located in main.cc.
//     NOTICE: The U.S. Government is granted for itself and others acting on its
//     behalf a paid-up, nonexclusive, irrevocable, worldwide license in this data to
//     reproduce, prepare derivative works, and perform publicly and display publicly.
//     Beginning five (5) years after permission to assert copyright is granted,
//     subject to two possible five year renewals, the U.S. Government is granted for
//     itself and others acting on its behalf a paid-up, non-exclusive, irrevocable
//     worldwide license in this data to reproduce, prepare derivative works,
//     distribute copies to the public, perform publicly and display publicly, and to
//     permit others to do so.
//     TRADEMARKS: EnergyPlus is a trademark of the US Department of Energy.

// Object Data
FArray1D< AirflowNetworkNodeSimuData > AirflowNetworkNodeSimu;
FArray1D< AirflowNetworkLinkSimuData > AirflowNetworkLinkSimu;
FArray1D< AirflowNetworkExchangeProp > AirflowNetworkExchangeData;
FArray1D< AirflowNetworkExchangeProp > AirflowNetworkMultiExchangeData;
FArray1D< AirflowNetworkLinkReportData > AirflowNetworkLinkReport;
FArray1D< AirflowNetworkNodeReportData > AirflowNetworkNodeReport;
FArray1D< AirflowNetworkLinkReportData > AirflowNetworkLinkReport1;
AirflowNetworkSimuProp AirflowNetworkSimu("", "NoMultizoneOrDistribution", "Input", 0, "", "", "", 500, 0, 1.0e-5, 1.0e-5, -0.5, 500.0, 0.0, 1.0, 0, 1.0e-4, 0, 0, 0, 0, "ZeroNodePressures", false); // unique object name | AirflowNetwork control | Wind pressure coefficient input control | Integer equivalent for WPCCntr field | CP Array name at WPCCntr = "INPUT" | Building type | Height Selection | Maximum number of iteration | Initialization flag | Relative airflow convergence | Absolute airflow convergence | Convergence acceleration limit | Maximum pressure change in an element [Pa] | Azimuth Angle of Long Axis of Building | Ratio of Building Width Along Short Axis to Width Along Long Axis | Number of wind directions | Minimum pressure difference | Exterior large opening error count during HVAC system operation | Exterior large opening error index during HVAC system operation | Large opening error count at Open factor > 1.0 | Large opening error error index at Open factor > 1.0 | Initialization flag type
FArray1D< AirflowNetworkNodeProp > AirflowNetworkNodeData;
FArray1D< AirflowNetworkCompProp > AirflowNetworkCompData;
FArray1D< AirflowNetworkLinkageProp > AirflowNetworkLinkageData;
FArray1D< MultizoneZoneProp > MultizoneZoneData;
FArray1D< MultizoneSurfaceProp > MultizoneSurfaceData;
FArray1D< MultizoneCompDetOpeningProp > MultizoneCompDetOpeningData;
FArray1D< MultizoneCompSimpleOpeningProp > MultizoneCompSimpleOpeningData;
FArray1D< MultizoneCompHorOpeningProp > MultizoneCompHorOpeningData;
FArray1D< MultizoneSurfaceCrackStdCndns > MultizoneSurfaceStdConditionsCrackData;
FArray1D< MultizoneSurfaceCrackProp > MultizoneSurfaceCrackData;
FArray1D< MultizoneSurfaceELAProp > MultizoneSurfaceELAData;
FArray1D< MultizoneExternalNodeProp > MultizoneExternalNodeData;
FArray1D< MultizoneCPArrayProp > MultizoneCPArrayData;
FArray1D< MultizoneCPArrayProp > MultizoneCPArrayDataSingleSided;
FArray1D< MultizoneCPValueProp > MultizoneCPValueData;
FArray1D< MultizoneCPValueProp > MultizoneCPValueDataTemp; // temporary CP values
FArray1D< MultizoneCPValueProp > MultizoneCPValueDataTempUnMod; // temporary CPValues, without modifcation factor
FArray1D< DeltaCpProp > DeltaCp;
FArray1D< DeltaCpProp > EPDeltaCP;
FArray1D< MultizoneCompExhaustFanProp > MultizoneCompExhaustFanData;
FArray1D< DisSysNodeProp > DisSysNodeData;
FArray1D< DisSysCompLeakProp > DisSysCompLeakData;
FArray1D< DisSysCompELRProp > DisSysCompELRData;
FArray1D< DisSysCompDuctProp > DisSysCompDuctData;
FArray1D< DisSysCompDamperProp > DisSysCompDamperData;
FArray1D< DisSysCompCVFProp > DisSysCompCVFData;
FArray1D< DisSysCompDetFanProp > DisSysCompDetFanData;
FArray1D< DisSysCompCoilProp > DisSysCompCoilData;
FArray1D< DisSysCompHXProp > DisSysCompHXData;
FArray1D< DisSysCompTermUnitProp > DisSysCompTermUnitData;
FArray1D< DisSysCompCPDProp > DisSysCompCPDData;
FArray1D< AiflowNetworkReportProp > AirflowNetworkReportData;

} // DataAirflowNetwork

namespace AirflowNetwork {

using namespace DataAirflowNetwork;
using namespace AirflowNetworkSolver;
using namespace Psychrometrics;

int SurfaceCrack::calcAfe( // Returns number of flows, either 1 or 2
    bool laminarInit, // Initialization flag. If true, use laminar relationship
    Real64 pressureDrop, // Total pressure drop across a component (P1 - P2) [Pa]
    int const i, // Linkage number
    int const n, // Node 1 number
    int const M, // Node 2 number
    FArray1A< Real64 > &F, // Airflow through the component [kg/s]
    FArray1A< Real64 > &DF // Partial derivative:  DF/DP
    )
{
    F.dim(2);
    DF.dim(2);

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    Real64 CDM;
    Real64 FL;
    Real64 FT;
    Real64 RhozNorm;
    Real64 VisczNorm;
    Real64 expn;
    Real64 Ctl;
    Real64 coef;
    Real64 Corr;
    Real64 VisAve;
    Real64 Tave;
    Real64 RhoCor;

    // Formats
    static gio::Fmt Format_901("(A5,I3,6X,4E16.7)");

    // FLOW:
    // Crack standard condition from given inputs
    Corr = MultizoneSurfaceData(i).Factor;
    RhozNorm = PsyRhoAirFnPbTdbW(StandardP, StandardT, StandardW);
    VisczNorm = 1.71432e-5 + 4.828e-8 * StandardT;

    expn = FlowExpo;
    VisAve = 0.5 * (VISCZ(n) + VISCZ(M));
    Tave = 0.5 * (TZ(n) + TZ(M));
    if(pressureDrop >= 0.0) {
        coef = FlowCoef / SQRTDZ(n) * Corr;
    } else {
        coef = FlowCoef / SQRTDZ(M) * Corr;
    }

    int NF = 1;
    if(laminarInit) {
        // Initialization by linear relation.
        if(pressureDrop >= 0.0) {
            RhoCor = (TZ(n) + KelvinConv) / (Tave + KelvinConv);
            Ctl = std::pow(RhozNorm / RHOZ(n) / RhoCor, expn - 1.0) * std::pow(VisczNorm / VisAve, 2.0 * expn - 1.0);
            DF(1) = coef * RHOZ(n) / VISCZ(n) * Ctl;
        } else {
            RhoCor = (TZ(M) + KelvinConv) / (Tave + KelvinConv);
            Ctl = std::pow(RhozNorm / RHOZ(M) / RhoCor, expn - 1.0) * std::pow(VisczNorm / VisAve, 2.0 * expn - 1.0);
            DF(1) = coef * RHOZ(M) / VISCZ(M) * Ctl;
        }
        F(1) = -DF(1) * pressureDrop;
    } else {
        // Standard calculation.
        if(pressureDrop >= 0.0) {
            // Flow in positive direction.
            // Laminar flow.
            RhoCor = (TZ(n) + KelvinConv) / (Tave + KelvinConv);
            Ctl = std::pow(RhozNorm / RHOZ(n) / RhoCor, expn - 1.0) * std::pow(VisczNorm / VisAve, 2.0 * expn - 1.0);
            CDM = coef * RHOZ(n) / VISCZ(n) * Ctl;
            FL = CDM * pressureDrop;
            // Turbulent flow.
            if(expn == 0.5) {
                FT = coef * SQRTDZ(n) * std::sqrt(pressureDrop) * Ctl;
            } else {
                FT = coef * SQRTDZ(n) * std::pow(pressureDrop, expn) * Ctl;
            }
        } else {
            // Flow in negative direction.
            // Laminar flow.
            RhoCor = (TZ(M) + KelvinConv) / (Tave + KelvinConv);
            Ctl = std::pow(RhozNorm / RHOZ(M) / RhoCor, expn - 1.0) * std::pow(VisczNorm / VisAve, 2.0 * expn - 1.0);
            CDM = coef * RHOZ(M) / VISCZ(M) * Ctl;
            FL = CDM * pressureDrop;
            // Turbulent flow.
            if(expn == 0.5) {
                FT = -coef * SQRTDZ(M) * std::sqrt(-pressureDrop) * Ctl;
            } else {
                FT = -coef * SQRTDZ(M) * std::pow(-pressureDrop, expn) * Ctl;
            }
        }
        // Select laminar or turbulent flow.
        if(LIST >= 4) gio::write(Unit21, Format_901) << " scr: " << i << pressureDrop << FL << FT;
        if(std::abs(FL) <= std::abs(FT)) {
            F(1) = FL;
            DF(1) = CDM;
        } else {
            F(1) = FT;
            DF(1) = FT * expn / pressureDrop;
        }
    }
    return NF;
}

int DetailedOpening::calcAfe(
    bool laminarInit, // Initialization flag.If true, use laminar relationship
    Real64 pressureDrop, // Total pressure drop across a component (P1 - P2) [Pa]
    int IL, // Linkage number
    int n, // Node 1 number
    int m, // Node 2 number
    FArray1A< Real64 > &F, // Airflow through the component [kg/s]
    FArray1A< Real64 > &DF // Partial derivative:  DF/DP
    )
{

    // SUBROUTINE INFORMATION:
    //       AUTHOR         Lixing Gu
    //       DATE WRITTEN   Oct. 2005
    //       MODIFIED       na
    //       RE-ENGINEERED  This subroutine is revised based on a vertical large opening subroutine from COMIS

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine simulates airflow and pressure of a detailed large opening component.

    // METHODOLOGY EMPLOYED:
    // Purpose:  This routine calculates the massflow and its derivative
    //       through a large opening in both flow directions. As input
    //       the density profiles RhoProfF/T are required aswell as the
    //       effective pressure difference profile DpProfNew, which is the
    //       sum of the stack pressure difference profile DpProf and the
    //       difference of the actual pressures at reference height. The
    //       profiles are calculated in the routine PresProfile.
    //       The massflow and its derivative are calculated for each
    //       interval representing a step of the pressure difference
    //       profile. The total flow and derivative are obtained by
    //       summation over the whole opening.
    //       The calculation is split into different cases representing
    //       different situations of the opening:
    //       - closed opening (opening factor = 0): summation of top and
    //         bottom crack (crack length = lwmax) plus "integration" over
    //         a vertically distributed crack of length (2*lhmax+lextra).
    //       - type 1: normal rectangular opening: "integration" over NrInt
    //         openings with width actlw and height actlh/NrInt
    //       - type 2: horizontally pivoted window: flow direction assumed
    //         strictly perpendicular to the plane of the opening
    //         -> "integration" over normal rectangular openings at top
    //         and bottom of LO plus a rectangular opening in series with two
    //         triangular openings in the middle of the LO (most general
    //         situation). The geometry is defined by the input parameters
    //         actlw(=lwmax), actlh, axisheight, opening angle.
    //       Assuming the massflow perpendicular to the opening plane in all
    //       cases the ownheightfactor has no influence on the massflow.

    // REFERENCES:
    // Helmut E. Feustel and Alison Rayner-Hooson, "COMIS Fundamentals," LBL-28560,
    // Lawrence Berkeley National Laboratory, Berkeley, CA, May 1990

    // USE STATEMENTS:
    using DataGlobals::PiOvr2;

    // Argument array dimensioning
    F.dim(2);
    DF.dim(2);

    // Locals
    // SUBROUTINE ARGUMENT DEFINITIONS:

    // SUBROUTINE PARAMETER DEFINITIONS:
    Real64 const RealMax(0.1e+37);
    Real64 const RealMin(1e-37);
    static Real64 const sqrt_1_2(std::sqrt(1.2));

    // INTERFACE BLOCK SPECIFICATIONS
    // na

    // DERIVED TYPE DEFINITIONS
    // na

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:

    static Real64 const sqrt_2(std::sqrt(2.0));

    Real64 Width;
    Real64 Height;

    Real64 fma12; // massflow in direction "from-to" [kg/s]
    Real64 fma21; // massflow in direction "to-from" [kg/s]
    Real64 dp1fma12; // derivative d fma12 / d Dp [kg/s/Pa]
    Real64 dp1fma21; // derivative d fma21 / d Dp [kg/s/Pa]
    FArray1D< Real64 > DpProfNew(NrInt + 2); // Differential pressure profile for Large Openings, taking into account fixed
    // pressures and actual zone pressures at reference height
    Real64 Fact; // Actual opening factor
    Real64 DifLim; // Limit for the pressure difference where laminarization takes place [Pa]
    Real64 Cfact;
    Real64 FvVeloc;

    Real64 ActLh;
    Real64 ActLw;
    Real64 Lextra;
    Real64 Axishght;
    Real64 ActCD;
    Real64 Cs;
    Real64 expn;
    Real64 Type;
    Real64 Interval;
    Real64 fmasum;
    Real64 dfmasum;
    Real64 Prefact;
    FArray1D< Real64 > EvalHghts(NrInt + 2);
    Real64 h2;
    Real64 h4;
    Real64 alpha;
    Real64 rholink;
    Real64 c1;
    Real64 c2;
    Real64 DpZeroOffset;
    Real64 area;
    Real64 WFact;
    Real64 HFact;
    int i;
    int Loc;
    int iNum;

    // FLOW:
    // Get component properties
    DifLim = 1.0e-4;
    //CompNum = AirflowNetworkCompData(j).TypeNum;
    Width = MultizoneSurfaceData(IL).Width;
    Height = MultizoneSurfaceData(IL).Height;
    Fact = MultizoneSurfaceData(IL).OpenFactor;
    Loc = (AirflowNetworkLinkageData(IL).DetOpenNum - 1) * (NrInt + 2);
    iNum = NumFac;
    ActCD = 0.0;

    if(iNum == 2) {
        if(Fact <= OpenFac2) {
            WFact = WidthFac1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (WidthFac2 - WidthFac1);
            HFact = HeightFac1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (HeightFac2 - HeightFac1);
            Cfact = DischCoeff1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (DischCoeff2 - DischCoeff1);
        } else {
            ShowFatalError("Open Factor is above the maximum input range for opening factors in AirflowNetwork:MultiZone:Component:DetailedOpening = " + Name);
        }
    }

    if(iNum == 3) {
        if(Fact <= OpenFac2) {
            WFact = WidthFac1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (WidthFac2 - WidthFac1);
            HFact = HeightFac1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (HeightFac2 - HeightFac1);
            Cfact = DischCoeff1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (DischCoeff2 - DischCoeff1);
        } else if(Fact <= OpenFac3) {
            WFact = WidthFac2 + (Fact - OpenFac2) / (OpenFac3 - OpenFac2) * (WidthFac3 - WidthFac2);
            HFact = HeightFac2 + (Fact - OpenFac2) / (OpenFac3 - OpenFac2) * (HeightFac3 - HeightFac2);
            Cfact = DischCoeff2 + (Fact - OpenFac2) / (OpenFac3 - OpenFac2) * (DischCoeff3 - DischCoeff2);
        } else {
            ShowFatalError("Open Factor is above the maximum input range for opening factors in AirflowNetwork:MultiZone:Component:DetailedOpening = " + Name);
        }
    }

    if(iNum == 4) {
        if(Fact <= OpenFac2) {
            WFact = WidthFac1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (WidthFac2 - WidthFac1);
            HFact = HeightFac1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (HeightFac2 - HeightFac1);
            Cfact = DischCoeff1 + (Fact - OpenFac1) / (OpenFac2 - OpenFac1) * (DischCoeff2 - DischCoeff1);
        } else if(Fact <= OpenFac3) {
            WFact = WidthFac2 + (Fact - OpenFac2) / (OpenFac3 - OpenFac2) * (WidthFac3 - WidthFac2);
            HFact = HeightFac2 + (Fact - OpenFac2) / (OpenFac3 - OpenFac2) * (HeightFac3 - HeightFac2);
            Cfact = DischCoeff2 + (Fact - OpenFac2) / (OpenFac3 - OpenFac2) * (DischCoeff3 - DischCoeff2);
        } else if(Fact <= OpenFac4) {
            WFact = WidthFac3 + (Fact - OpenFac3) / (OpenFac4 - OpenFac3) * (WidthFac4 - WidthFac3);
            HFact = HeightFac3 + (Fact - OpenFac3) / (OpenFac4 - OpenFac3) * (HeightFac4 - HeightFac3);
            Cfact = DischCoeff3 + (Fact - OpenFac3) / (OpenFac4 - OpenFac3) * (DischCoeff4 - DischCoeff3);
        } else {
            ShowFatalError("Open Factor is above the maximum input range for opening factors in AirflowNetwork:MultiZone:Component:DetailedOpening = " + Name);
        }
    }

    // calculate DpProfNew
    for(i = 1; i <= NrInt + 2; ++i) {
        DpProfNew(i) = pressureDrop + DpProf(Loc + i) - DpL(1, IL);
    }

    // Get opening data based on the opening factor
    if(Fact == 0) {
        ActLw = MultizoneSurfaceData(IL).Width;
        ActLh = MultizoneSurfaceData(IL).Height;
        Cfact = 0.0;
    } else {
        ActLw = MultizoneSurfaceData(IL).Width * WFact;
        ActLh = MultizoneSurfaceData(IL).Height * HFact;
        ActCD = Cfact;
    }

    Cs = FlowCoef;
    expn = FlowExpo;
    Type = LVOType;
    if(Type == 1) {
        Lextra = LVOValue;
        Axishght = 0.0;
    } else if(Type == 2) {
        Lextra = 0.0;
        Axishght = LVOValue;
        ActLw = MultizoneSurfaceData(IL).Width;
        ActLh = MultizoneSurfaceData(IL).Height;
    }

    // Add window multiplier with window close
    if(MultizoneSurfaceData(IL).Multiplier > 1.0) Cs *= MultizoneSurfaceData(IL).Multiplier;
    // Add window multiplier with window open
    if(Fact > 0.0) {
        if(MultizoneSurfaceData(IL).Multiplier > 1.0) ActLw *= MultizoneSurfaceData(IL).Multiplier;
    }

    // Add recurring warnings
    if(Fact > 0.0) {
        if(ActLw == 0.0) {
            ++WidthErrCount;
            if(WidthErrCount < 2) {
                ShowWarningError("The actual width of the AirflowNetwork:MultiZone:Component:DetailedOpening of " + Name + " is 0.");
                ShowContinueError("The actual width is set to 1.0E-6 m.");
                ShowContinueErrorTimeStamp("Occurrence info:");
            } else {
                ShowRecurringWarningErrorAtEnd("The actual width of the AirflowNetwork:MultiZone:Component:DetailedOpening of " + Name + " is 0 error continues.", WidthErrIndex, ActLw, ActLw);
            }
            ActLw = 1.0e-6;
        }
        if(ActLh == 0.0) {
            ++HeightErrCount;
            if(HeightErrCount < 2) {
                ShowWarningError("The actual height of the AirflowNetwork:MultiZone:Component:DetailedOpening of " + Name + " is 0.");
                ShowContinueError("The actual height is set to 1.0E-6 m.");
                ShowContinueErrorTimeStamp("Occurrence info:");
            } else {
                ShowRecurringWarningErrorAtEnd("The actual width of the AirflowNetwork:MultiZone:Component:DetailedOpening of " + Name + " is 0 error continues.", HeightErrIndex, ActLh, ActLh);
            }
            ActLh = 1.0e-6;
        }
    }
    // Initialization:
    int NF = 1;
    Interval = ActLh / NrInt;
    fma12 = 0.0;
    fma21 = 0.0;
    dp1fma12 = 0.0;
    dp1fma21 = 0.0;

    // Closed LO
    if(Cfact == 0) {
        DpZeroOffset = DifLim;
        // bottom crack
        if(DpProfNew(1) > 0) {
            if(std::abs(DpProfNew(1)) <= DpZeroOffset) {
                dfmasum = Cs * ActLw * std::pow(DpZeroOffset, expn) / DpZeroOffset;
                fmasum = DpProfNew(1) * dfmasum;
            } else {
                fmasum = Cs * ActLw * std::pow(DpProfNew(1), expn);
                dfmasum = fmasum * expn / DpProfNew(1);
            }
            fma12 += fmasum;
            dp1fma12 += dfmasum;
        } else {
            if(std::abs(DpProfNew(1)) <= DpZeroOffset) {
                dfmasum = -Cs * ActLw * std::pow(DpZeroOffset, expn) / DpZeroOffset;
                fmasum = DpProfNew(1) * dfmasum;
            } else {
                fmasum = Cs * ActLw * std::pow(-DpProfNew(1), expn);
                dfmasum = fmasum * expn / DpProfNew(1);
            }
            fma21 += fmasum;
            dp1fma21 += dfmasum;
        }
        // top crack
        if(DpProfNew(NrInt + 2) > 0) {
            if(std::abs(DpProfNew(NrInt + 2)) <= DpZeroOffset) {
                dfmasum = Cs * ActLw * std::pow(DpZeroOffset, expn) / DpZeroOffset;
                fmasum = DpProfNew(NrInt + 2) * dfmasum;
            } else {
                fmasum = Cs * ActLw * std::pow(DpProfNew(NrInt + 2), expn);
                dfmasum = fmasum * expn / DpProfNew(NrInt + 2);
            }
            fma12 += fmasum;
            dp1fma12 += dfmasum;
        } else {
            if(std::abs(DpProfNew(NrInt + 2)) <= DpZeroOffset) {
                dfmasum = -Cs * ActLw * std::pow(DpZeroOffset, expn) / DpZeroOffset;
                fmasum = DpProfNew(NrInt + 2) * dfmasum;
            } else {
                fmasum = Cs * ActLw * std::pow(-DpProfNew(NrInt + 2), expn);
                dfmasum = fmasum * expn / DpProfNew(NrInt + 2);
            }
            fma21 += fmasum;
            dp1fma21 += dfmasum;
        }
        // side and extra cracks
        Prefact = Interval * (2 + Lextra / ActLh) * Cs;
        for(i = 2; i <= NrInt + 1; ++i) {
            if(DpProfNew(i) > 0) {
                if(std::abs(DpProfNew(i)) <= DpZeroOffset) {
                    dfmasum = Prefact * std::pow(DpZeroOffset, expn) / DpZeroOffset;
                    fmasum = DpProfNew(i) * dfmasum;
                } else {
                    fmasum = Prefact * std::pow(DpProfNew(i), expn);
                    dfmasum = fmasum * expn / DpProfNew(i);
                }
                fma12 += fmasum;
                dp1fma12 += dfmasum;
            } else {
                if(std::abs(DpProfNew(i)) <= DpZeroOffset) {
                    dfmasum = -Prefact * std::pow(DpZeroOffset, expn) / DpZeroOffset;
                    fmasum = DpProfNew(i) * dfmasum;
                } else {
                    fmasum = Prefact * std::pow(-DpProfNew(i), expn);
                    dfmasum = fmasum * expn / DpProfNew(i);
                }
                fma21 += fmasum;
                dp1fma21 += dfmasum;
            }
        }
    }

    // Open LO, type 1
    if((Cfact != 0) && (Type == 1)) {
        DpZeroOffset = DifLim * 1e-3;
        Prefact = ActLw * ActCD * Interval * sqrt_2;
        for(i = 2; i <= NrInt + 1; ++i) {
            if(DpProfNew(i) > 0) {
                if(std::abs(DpProfNew(i)) <= DpZeroOffset) {
                    dfmasum = std::sqrt(RhoProfF(Loc + i) * DpZeroOffset) / DpZeroOffset;
                    fmasum = DpProfNew(i) * dfmasum;
                } else {
                    fmasum = std::sqrt(RhoProfF(Loc + i) * DpProfNew(i));
                    dfmasum = 0.5 * fmasum / DpProfNew(i);
                }
                fma12 += fmasum;
                dp1fma12 += dfmasum;
            } else {
                if(std::abs(DpProfNew(i)) <= DpZeroOffset) {
                    dfmasum = -std::sqrt(RhoProfT(Loc + i) * DpZeroOffset) / DpZeroOffset;
                    fmasum = DpProfNew(i) * dfmasum;
                } else {
                    fmasum = std::sqrt(-RhoProfT(Loc + i) * DpProfNew(i));
                    dfmasum = 0.5 * fmasum / DpProfNew(i);
                }
                fma21 += fmasum;
                dp1fma21 += dfmasum;
            }
        }

        fma12 *= Prefact;
        fma21 *= Prefact;
        dp1fma12 *= Prefact;
        dp1fma21 *= Prefact;

    }

    // Open LO, type 2
    if((Cfact != 0) && (Type == 2)) {
        // Initialization
        DpZeroOffset = DifLim * 1e-3;
        // New definition for opening factors for LVO type 2: opening angle = 90 degrees --> opening factor = 1.0
        // should be PIOvr2 in below?
        alpha = Fact * PiOvr2;
        Real64 const cos_alpha(std::cos(alpha));
        Real64 const tan_alpha(std::tan(alpha));
        h2 = Axishght * (1.0 - cos_alpha);
        h4 = Axishght + (ActLh - Axishght) * cos_alpha;
        EvalHghts(1) = 0.0;
        EvalHghts(NrInt + 2) = ActLh;
        // New definition for opening factors for LVO type 2: pening angle = 90 degrees --> opening factor = 1.0
        if(Fact == 1.0) {
            h2 = Axishght;
            h4 = Axishght;
        }

        for(i = 2; i <= NrInt + 1; ++i) {
            EvalHghts(i) = Interval * (i - 1.5);
        }

        // Calculation of massflow and its derivative
        for(i = 2; i <= NrInt + 1; ++i) {
            if(DpProfNew(i) > 0) {
                rholink = RhoProfF(Loc + i);
            } else {
                rholink = RhoProfT(Loc + i);
            }

            if((EvalHghts(i) <= h2) || (EvalHghts(i) >= h4)) {
                if(std::abs(DpProfNew(i)) <= DpZeroOffset) {
                    dfmasum = ActCD * ActLw * Interval * std::sqrt(2.0 * rholink * DpZeroOffset) / DpZeroOffset * sign(1, DpProfNew(i));
                    fmasum = DpProfNew(i) * dfmasum;
                } else {
                    fmasum = ActCD * ActLw * Interval * std::sqrt(2.0 * rholink * std::abs(DpProfNew(i)));
                    dfmasum = 0.5 * fmasum / DpProfNew(i);
                }
            } else {
                // triangular opening at the side of LO
                c1 = ActCD * ActLw * Interval * std::sqrt(2.0 * rholink);
                c2 = 2 * ActCD * std::abs(Axishght - EvalHghts(i)) * tan_alpha * Interval * std::sqrt(2.0 * rholink);
                if((c1 != 0) && (c2 != 0)) {
                    if(std::abs(DpProfNew(i)) <= DpZeroOffset) {
                        dfmasum = std::sqrt(DpZeroOffset / (1 / c1 / c1 + 1 / c2 / c2)) / DpZeroOffset * sign(1, DpProfNew(i));
                        fmasum = DpProfNew(i) * dfmasum;
                    } else {
                        fmasum = std::sqrt(std::abs(DpProfNew(i)) / (1 / c1 / c1 + 1 / c2 / c2));
                        dfmasum = 0.5 * fmasum / DpProfNew(i);
                    }
                } else {
                    fmasum = 0.0;
                    dfmasum = 0.0;
                }
            }

            if(DpProfNew(i) > 0) {
                fma12 += fmasum;
                dp1fma12 += dfmasum;
            } else {
                fma21 += fmasum;
                dp1fma21 += dfmasum;
            }
        }

    }

    // Calculate some velocity in the large opening
    area = ActLh * ActLw * ActCD;
    if(area > (Cs + RealMin)) {
        if(area > RealMin) {
            FvVeloc = (fma21 + fma12) / area;
        } else {
            FvVeloc = 0.0;
        }
    } else {
        // here the average velocity over the full area, may blow half in half out.
        // velocity= Fva/Nett area=Fma/Rho/(Cm/( (2**N)* SQRT(1.2) ) )
        if(Cs > 0.0) {
            // get the average Rho for this closed window
            for(i = 2; i <= NrInt + 1; ++i) {
                rholink = 0.0;
                if(DpProfNew(i) > 0) {
                    rholink = RhoProfF(Loc + i);
                } else {
                    rholink = RhoProfT(Loc + i);
                }
                rholink /= NrInt;
                rholink = 1.2;
            }
            FvVeloc = (fma21 + fma12) * std::pow(2.0, expn) * sqrt_1_2 / (rholink * Cs);
        } else {
            FvVeloc = 0.0;
        }
    }

    // Output mass flow rates and associated derivatives
    F(1) = fma12 - fma21;
    DF(1) = dp1fma12 - dp1fma21;
    F(2) = 0.0;
    if(fma12 != 0.0 && fma21 != 0.0) {
        F(2) = fma21;
    }
    DF(2) = 0.0;
    return NF;
}

    }

} // EnergyPlus
