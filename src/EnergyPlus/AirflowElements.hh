#ifndef AIRFLOWELEMENTS_HH
#define AIRFLOWELEMENTS_HH

// ObjexxFCL Headers
#include <ObjexxFCL/Array1D.hh>

// EnergyPlus Headers
#include <EnergyPlus.hh>
#include <DataGlobals.hh>

namespace EnergyPlus {

namespace AirflowNetwork {

	// Using/Aliasing

	// Data
	// module should be available to other modules and routines.  Thus,
	// all variables in this module must be PUBLIC.

	// MODULE PARAMETER DEFINITIONS:
	extern int const CompTypeNum_DOP; // Detailed large opening component
	extern int const CompTypeNum_SOP; // Simple opening component
	extern int const CompTypeNum_SCR; // Surface crack component
	extern int const CompTypeNum_SEL; // Surface effective leakage ratio component
	extern int const CompTypeNum_PLR; // Distribution system crack component
	extern int const CompTypeNum_DWC; // Distribution system duct component
	extern int const CompTypeNum_CVF; // Distribution system constant volume fan component
	extern int const CompTypeNum_FAN; // Distribution system detailed fan component
	extern int const CompTypeNum_MRR; // Distribution system multiple curve fit power law resistant flow component
	extern int const CompTypeNum_DMP; // Distribution system damper component
	extern int const CompTypeNum_ELR; // Distribution system effective leakage ratio component
	extern int const CompTypeNum_CPD; // Distribution system constant pressure drop component
	extern int const CompTypeNum_COI; // Distribution system coil component
	extern int const CompTypeNum_TMU; // Distribution system terminal unit component
	extern int const CompTypeNum_EXF; // Zone exhaust fan
	extern int const CompTypeNum_HEX; // Distribution system heat exchanger
	extern int const CompTypeNum_HOP; // Horizontal opening component
	extern int const CompTypeNum_RVD; // Reheat VAV terminal damper

	// EPlus component Type
	extern int const EPlusTypeNum_SCN; // Supply connection
	extern int const EPlusTypeNum_RCN; // Return connection
	extern int const EPlusTypeNum_RHT; // Reheat terminal
	extern int const EPlusTypeNum_FAN; // Fan
	extern int const EPlusTypeNum_COI; // Heating or cooling coil
	extern int const EPlusTypeNum_HEX; // Heat ecxchanger
	extern int const EPlusTypeNum_RVD; // Reheat VAV terminal damper

	// EPlus node type
	extern int const EPlusTypeNum_ZIN; // Zone inlet node
	extern int const EPlusTypeNum_ZOU; // Zone outlet node
	extern int const EPlusTypeNum_SPL; // Splitter node
	extern int const EPlusTypeNum_MIX; // Mixer node
	extern int const EPlusTypeNum_OAN; // Outside air system node
	extern int const EPlusTypeNum_EXT; // OA system inlet node
	extern int const EPlusTypeNum_FIN; // Fan Inlet node
	extern int const EPlusTypeNum_FOU; // Fan Outlet Node
	extern int const EPlusTypeNum_COU; // Coil Outlet Node
	extern int const EPlusTypeNum_HXO; // Heat exchanger Outlet Node
	extern int const EPlusTypeNum_DIN; // Damper Inlet node
	extern int const EPlusTypeNum_DOU; // Damper Outlet Node
	extern int const EPlusTypeNum_SPI; // Splitter inlet Node
	extern int const EPlusTypeNum_SPO; // Splitter Outlet Node

	extern int const iWPCCntr_Input;
	extern int const iWPCCntr_SurfAvg;

	// DERIVED TYPE DEFINITIONS:

	// MODULE VARIABLE DECLARATIONS:
	// Node simulation variable in air distribution system
	// Link simulation variable in air distribution system
	// Sensible and latent exchange variable in air distribution system

	extern int SimulateAirflowNetwork;
	// Vent Control  DistSys Control  Flag    Description
	//  NONE           NONE           0      No AirflowNetwork and SIMPLE
	//  SIMPLE         NONE           1      Simple calculations only
	//  MULTIZONE      NONE           2      Perform multizone calculations only
	//  NONE           DISTSYS        3      Perform distribution system durin system on time only
	//  SIMPLE         DISTSYS        4      Perform distribution system durin system on time and simple calculations during off time
	//  MULTIZONE      DISTSYS        5      Perform distribution system durin system on time and multizone calculations during off time

	extern int const AirflowNetworkControlSimple; // Simple calculations only
	extern int const AirflowNetworkControlMultizone; // Perform multizone calculations only
	extern int const AirflowNetworkControlSimpleADS; // Perform distribution system durin system
	// on time and simple calculations during off time
	extern int const AirflowNetworkControlMultiADS; // Perform distribution system durin system on time
	// and multizone calculations during off time

	extern Array1D_bool AirflowNetworkZoneFlag;

	extern int NumOfNodesMultiZone; // Number of nodes for multizone calculation
	extern int NumOfNodesDistribution; // Number of nodes for distribution system calculation
	extern int NumOfLinksMultiZone; // Number of links for multizone calculation
	extern int NumOfLinksDistribution; // Number of links for distribution system calculation

	extern int AirflowNetworkNumOfNodes; // Number of nodes for AirflowNetwork calculation
	// = NumOfNodesMultiZone+NumOfNodesDistribution
	extern int AirflowNetworkNumOfComps; // Number of components for AirflowNetwork calculation
	extern int AirflowNetworkNumOfLinks; // Number of links for AirflowNetwork calculation
	// = NumOfLinksMultiZone+NumOfLinksDistribution
	// RoomAirManager use
	extern int AirflowNetworkNumOfSurfaces; // The number of surfaces for multizone calculation
	extern int AirflowNetworkNumOfZones; // The number of zones for multizone calculation

	extern bool RollBackFlag; // Roll back flag when system time steo down shifting
	extern Array1D< Real64 > ANZT; // Local zone air temperature for roll back use
	extern Array1D< Real64 > ANZW; // Local zone air humidity ratio for roll back use
	extern Array1D< Real64 > ANCO; // Local zone air CO2 for roll back use
	extern Array1D< Real64 > ANGC; // Local zone air generic contaminant for roll back use
	extern int AirflowNetworkNumOfExhFan; // Number of zone exhaust fans
	extern Array1D_bool AirflowNetworkZoneExhaustFan; // Logical to use zone exhaust fans
	extern bool AirflowNetworkFanActivated; // Supply fan activation flag
	extern bool AirflowNetworkUnitarySystem; // set to TRUE for unitary systems (to make answers equal, will remove eventually)
	// Multispeed HP only
	extern int MultiSpeedHPIndicator; // Indicator for multispeed heat pump use
	// Addiitonal airflow needed for an VAV fan to compensate the leakage losses and supply pathway pressure losses [kg/s]
	extern Real64 VAVTerminalRatio; // The terminal flow ratio when a supply VAV fan reach its max flow rate
	extern bool VAVSystem; // This flag is used to represent a VAV system

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

	// Types

	struct DetailedOpening // Large detailed opening component
	{
		// Members
		std::string Name; // Name of large detailed opening component
		Real64 FlowCoef; // Air Mass Flow Coefficient When Window or Door Is Closed
		Real64 FlowExpo; // Air Mass Flow exponent When Window or Door Is Closed
		std::string TypeName; // Name of Large vertical opening type
		int LVOType; // Large vertical opening type number
		Real64 LVOValue; // Extra crack length for LVO type 1 with multiple openable parts,
		// or Height of pivoting axis for LVO type 2
		int NumFac; // Number of Opening Factor Values
		Real64 OpenFac1; // Opening factor #1
		Real64 DischCoeff1; // Discharge coefficient for opening factor #1
		Real64 WidthFac1; // Width factor for for Opening factor #1
		Real64 HeightFac1; // Height factor for opening factor #1
		Real64 StartHFac1; // Start height factor for opening factor #1
		Real64 OpenFac2; // Opening factor #2
		Real64 DischCoeff2; // Discharge coefficient for opening factor #2
		Real64 WidthFac2; // Width factor for for Opening factor #2
		Real64 HeightFac2; // Height factor for opening factor #2
		Real64 StartHFac2; // Start height factor for opening factor #2
		Real64 OpenFac3; // Opening factor #3
		Real64 DischCoeff3; // Discharge coefficient for opening factor #3
		Real64 WidthFac3; // Width factor for for Opening factor #3
		Real64 HeightFac3; // Height factor for opening factor #3
		Real64 StartHFac3; // Start height factor for opening factor #3
		Real64 OpenFac4; // Opening factor #4
		Real64 DischCoeff4; // Discharge coefficient for opening factor #4
		Real64 WidthFac4; // Width factor for for Opening factor #4
		Real64 HeightFac4; // Height factor for opening factor #4
		Real64 StartHFac4; // Start height factor for opening factor #4
		Real64 OpenFactor; // Opening factor
		int WidthErrCount; // Width error count
		int WidthErrIndex; // Width error index
		int HeightErrCount; // Height error count
		int HeightErrIndex; // Height error index

		// Default Constructor
		DetailedOpening() :
			FlowCoef( 0.0 ),
			FlowExpo( 0.0 ),
			TypeName( "NONPIVOTED" ),
			LVOType( 0 ),
			LVOValue( 0.0 ),
			NumFac( 0 ),
			OpenFac1( 0.0 ),
			DischCoeff1( 0.0 ),
			WidthFac1( 0.0 ),
			HeightFac1( 0.0 ),
			StartHFac1( 0.0 ),
			OpenFac2( 0.0 ),
			DischCoeff2( 0.0 ),
			WidthFac2( 0.0 ),
			HeightFac2( 0.0 ),
			StartHFac2( 0.0 ),
			OpenFac3( 0.0 ),
			DischCoeff3( 0.0 ),
			WidthFac3( 0.0 ),
			HeightFac3( 0.0 ),
			StartHFac3( 0.0 ),
			OpenFac4( 0.0 ),
			DischCoeff4( 0.0 ),
			WidthFac4( 0.0 ),
			HeightFac4( 0.0 ),
			StartHFac4( 0.0 ),
			OpenFactor( 0.0 ),
			WidthErrCount( 0 ),
			WidthErrIndex( 0 ),
			HeightErrCount( 0 ),
			HeightErrIndex( 0 )
		{}

		// Member Constructor
		DetailedOpening(
			std::string const & Name, // Name of large detailed opening component
			Real64 const FlowCoef, // Air Mass Flow Coefficient When Window or Door Is Closed
			Real64 const FlowExpo, // Air Mass Flow exponent When Window or Door Is Closed
			std::string const & TypeName, // Name of Large vertical opening type
			int const LVOType, // Large vertical opening type number
			Real64 const LVOValue, // Extra crack length for LVO type 1 with multiple openable parts,
			int const NumFac, // Number of Opening Factor Values
			Real64 const OpenFac1, // Opening factor #1
			Real64 const DischCoeff1, // Discharge coefficient for opening factor #1
			Real64 const WidthFac1, // Width factor for for Opening factor #1
			Real64 const HeightFac1, // Height factor for opening factor #1
			Real64 const StartHFac1, // Start height factor for opening factor #1
			Real64 const OpenFac2, // Opening factor #2
			Real64 const DischCoeff2, // Discharge coefficient for opening factor #2
			Real64 const WidthFac2, // Width factor for for Opening factor #2
			Real64 const HeightFac2, // Height factor for opening factor #2
			Real64 const StartHFac2, // Start height factor for opening factor #2
			Real64 const OpenFac3, // Opening factor #3
			Real64 const DischCoeff3, // Discharge coefficient for opening factor #3
			Real64 const WidthFac3, // Width factor for for Opening factor #3
			Real64 const HeightFac3, // Height factor for opening factor #3
			Real64 const StartHFac3, // Start height factor for opening factor #3
			Real64 const OpenFac4, // Opening factor #4
			Real64 const DischCoeff4, // Discharge coefficient for opening factor #4
			Real64 const WidthFac4, // Width factor for for Opening factor #4
			Real64 const HeightFac4, // Height factor for opening factor #4
			Real64 const StartHFac4, // Start height factor for opening factor #4
			Real64 const OpenFactor, // Opening factor
			int const WidthErrCount, // Width error count
			int const WidthErrIndex, // Width error index
			int const HeightErrCount, // Height error count
			int const HeightErrIndex // Height error index
		) :
			Name( Name ),
			FlowCoef( FlowCoef ),
			FlowExpo( FlowExpo ),
			TypeName( TypeName ),
			LVOType( LVOType ),
			LVOValue( LVOValue ),
			NumFac( NumFac ),
			OpenFac1( OpenFac1 ),
			DischCoeff1( DischCoeff1 ),
			WidthFac1( WidthFac1 ),
			HeightFac1( HeightFac1 ),
			StartHFac1( StartHFac1 ),
			OpenFac2( OpenFac2 ),
			DischCoeff2( DischCoeff2 ),
			WidthFac2( WidthFac2 ),
			HeightFac2( HeightFac2 ),
			StartHFac2( StartHFac2 ),
			OpenFac3( OpenFac3 ),
			DischCoeff3( DischCoeff3 ),
			WidthFac3( WidthFac3 ),
			HeightFac3( HeightFac3 ),
			StartHFac3( StartHFac3 ),
			OpenFac4( OpenFac4 ),
			DischCoeff4( DischCoeff4 ),
			WidthFac4( WidthFac4 ),
			HeightFac4( HeightFac4 ),
			StartHFac4( StartHFac4 ),
			OpenFactor( OpenFactor ),
			WidthErrCount( WidthErrCount ),
			WidthErrIndex( WidthErrIndex ),
			HeightErrCount( HeightErrCount ),
			HeightErrIndex( HeightErrIndex )
		{}

	};

	struct SimpleOpening // Large simple opening component
	{
		// Members
		std::string Name; // Name of large simple opening component
		Real64 FlowCoef; // Air Mass Flow Coefficient When Window or Door Is Closed
		Real64 FlowExpo; // Air Mass Flow exponent When Window or Door Is Closed
		Real64 MinRhoDiff; // Minimum density difference for two-way flow
		Real64 DischCoeff; // Discharge coefficient at full opening
		Real64 OpenFactor; // Opening factor

		// Default Constructor
		SimpleOpening() :
			FlowCoef( 0.0 ),
			FlowExpo( 0.0 ),
			MinRhoDiff( 0.0 ),
			DischCoeff( 0.0 ),
			OpenFactor( 0.0 )
		{}

		// Member Constructor
		SimpleOpening(
			std::string const & Name, // Name of large simple opening component
			Real64 const FlowCoef, // Air Mass Flow Coefficient When Window or Door Is Closed
			Real64 const FlowExpo, // Air Mass Flow exponent When Window or Door Is Closed
			Real64 const MinRhoDiff, // Minimum density difference for two-way flow
			Real64 const DischCoeff, // Discharge coefficient at full opening
			Real64 const OpenFactor // Opening factor
		) :
			Name( Name ),
			FlowCoef( FlowCoef ),
			FlowExpo( FlowExpo ),
			MinRhoDiff( MinRhoDiff ),
			DischCoeff( DischCoeff ),
			OpenFactor( OpenFactor )
		{}

	};

	struct HorizontalOpening // Large horizontal opening component
	{
		// Members
		std::string Name; // Name of large horizontal opening component
		Real64 FlowCoef; // Air Mass Flow Coefficient When Window or Door Is Closed
		Real64 FlowExpo; // Air Mass Flow exponent When Window or Door Is Closed
		Real64 Slope; // Sloping plane angle
		Real64 DischCoeff; // Discharge coefficient at full opening

		// Default Constructor
		HorizontalOpening() :
			FlowCoef( 0.0 ),
			FlowExpo( 0.0 ),
			Slope( 0.0 ),
			DischCoeff( 0.0 )
		{}

		// Member Constructor
		HorizontalOpening(
			std::string const & Name, // Name of large horizontal opening component
			Real64 const FlowCoef, // Air Mass Flow Coefficient When Window or Door Is Closed
			Real64 const FlowExpo, // Air Mass Flow exponent When Window or Door Is Closed
			Real64 const Slope, // Sloping plane angle
			Real64 const DischCoeff // Discharge coefficient at full opening
		) :
			Name( Name ),
			FlowCoef( FlowCoef ),
			FlowExpo( FlowExpo ),
			Slope( Slope ),
			DischCoeff( DischCoeff )
		{}

	};

	struct SurfaceCrackStandardConditions // Surface crack standard conditions
	{
		// Members
		std::string Name; // Name of standard conditions component
		Real64 StandardT; // Standard temperature for crack data
		Real64 StandardP; // Standard borometric pressure for crack data
		Real64 StandardW; // Standard humidity ratio for crack data

		// Default Constructor
		SurfaceCrackStandardConditions() :
			StandardT( 0.0 ),
			StandardP( 0.0 ),
			StandardW( 0.0 )
		{}

		// Member Constructor
		MultizoneSurfaceCrackStdCndns(
			std::string const & Name, // Name of standard conditions component
			Real64 const StandardT, // Standard temperature for crack data
			Real64 const StandardP, // Standard borometric pressure for crack data
			Real64 const StandardW // Standard humidity ratio for crack data
		) :
			Name( Name ),
			StandardT( StandardT ),
			StandardP( StandardP ),
			StandardW( StandardW )
		{}

	};

	struct SurfaceCrack // Surface crack component
	{
		// Members
		std::string Name; // Name of crack component
		std::string ExternalNodeNames; // Name of external node.Not requird for internal surface
		Real64 FlowCoef; // Air Mass Flow Coefficient When Window or Door Is Closed
		Real64 FlowExpo; // Air Mass Flow exponent When Window or Door Is Closed
		Real64 StandardT; // Standard temperature for crack data
		Real64 StandardP; // Standard borometric pressure for crack data
		Real64 StandardW; // Standard humidity ratio for crack data

		// Default Constructor
		SurfaceCrack() :
			FlowCoef( 0.0 ),
			FlowExpo( 0.0 ),
			StandardT( 0.0 ),
			StandardP( 0.0 ),
			StandardW( 0.0 )
		{}

		// Member Constructor
		SurfaceCrack(
			std::string const & Name, // Name of crack component
			std::string const & ExternalNodeNames, // Name of external node.Not requird for internal surface
			Real64 const FlowCoef, // Air Mass Flow Coefficient When Window or Door Is Closed
			Real64 const FlowExpo, // Air Mass Flow exponent When Window or Door Is Closed
			Real64 const StandardT, // Standard temperature for crack data
			Real64 const StandardP, // Standard borometric pressure for crack data
			Real64 const StandardW // Standard humidity ratio for crack data
		) :
			Name( Name ),
			ExternalNodeNames( ExternalNodeNames ),
			FlowCoef( FlowCoef ),
			FlowExpo( FlowExpo ),
			StandardT( StandardT ),
			StandardP( StandardP ),
			StandardW( StandardW )
		{}

	};

	struct SurfaceEffectiveLeakageArea // Surface effective leakage area component
	{
		// Members
		std::string Name; // Name of effective leakage area component
		Real64 ELA; // Effective leakage area
		Real64 DischCoeff; // Discharge coefficient
		Real64 RefDeltaP; // Reference pressure difference
		Real64 FlowExpo; // Air Mass Flow exponent When Window or Door Is Closed
		Real64 TestDeltaP; // Testing pressure difference
		Real64 TestDisCoef; // Testing Discharge coefficient

		// Default Constructor
    SurfaceEffectiveLeakageArea() :
			ELA( 0.0 ),
			DischCoeff( 0.0 ),
			RefDeltaP( 0.0 ),
			FlowExpo( 0.0 ),
			TestDeltaP( 0.0 ),
			TestDisCoef( 0.0 )
		{}

		// Member Constructor
    SurfaceEffectiveLeakageArea(
			std::string const & Name, // Name of effective leakage area component
			Real64 const ELA, // Effective leakage area
			Real64 const DischCoeff, // Discharge coefficient
			Real64 const RefDeltaP, // Reference pressure difference
			Real64 const FlowExpo, // Air Mass Flow exponent When Window or Door Is Closed
			Real64 const TestDeltaP, // Testing pressure difference
			Real64 const TestDisCoef // Testing Discharge coefficient
		) :
			Name( Name ),
			ELA( ELA ),
			DischCoeff( DischCoeff ),
			RefDeltaP( RefDeltaP ),
			FlowExpo( FlowExpo ),
			TestDeltaP( TestDeltaP ),
			TestDisCoef( TestDisCoef )
		{}

	};

	struct ExhaustFan // Zone exhaust fan component
	{
		// Members
		std::string Name; // Name of exhaust fan component
		Real64 FlowRate; // mass flow rate
		int SchedPtr; // Schedule pointer
		Real64 FlowCoef; // Air Mass Flow Coefficient
		Real64 FlowExpo; // Air Mass Flow exponent
		Real64 StandardT; // Standard temperature for crack data
		Real64 StandardP; // Standard borometric pressure for crack data
		Real64 StandardW; // Standard humidity ratio for crack data
		int InletNode; // Inlet node number
		int OutletNode; // Outlet node number
		int EPlusZoneNum; // Zone number

		// Default Constructor
		ExhaustFan() :
			FlowRate( 0.0 ),
			SchedPtr( 0 ),
			FlowCoef( 0.0 ),
			FlowExpo( 0.0 ),
			StandardT( 0.0 ),
			StandardP( 0.0 ),
			StandardW( 0.0 ),
			InletNode( 0 ),
			OutletNode( 0 ),
			EPlusZoneNum( 0 )
		{}

		// Member Constructor
		ExhaustFan(
			std::string const & Name, // Name of exhaust fan component
			Real64 const FlowRate, // mass flow rate
			int const SchedPtr, // Schedule pointer
			Real64 const FlowCoef, // Air Mass Flow Coefficient
			Real64 const FlowExpo, // Air Mass Flow exponent
			Real64 const StandardT, // Standard temperature for crack data
			Real64 const StandardP, // Standard borometric pressure for crack data
			Real64 const StandardW, // Standard humidity ratio for crack data
			int const InletNode, // Inlet node number
			int const OutletNode, // Outlet node number
			int const EPlusZoneNum // Zone number
		) :
			Name( Name ),
			FlowRate( FlowRate ),
			SchedPtr( SchedPtr ),
			FlowCoef( FlowCoef ),
			FlowExpo( FlowExpo ),
			StandardT( StandardT ),
			StandardP( StandardP ),
			StandardW( StandardW ),
			InletNode( InletNode ),
			OutletNode( OutletNode ),
			EPlusZoneNum( EPlusZoneNum )
		{}

	};

	struct DuctLeak // duct leak component
	{
		// Members
		std::string Name; // Name of component leak
		Real64 FlowCoef; // Air Mass Flow Coefficient
		Real64 FlowExpo; // Air Mass Flow exponent

		// Default Constructor
		DuctLeak() :
			FlowCoef( 0.0 ),
			FlowExpo( 0.0 )
		{}

		// Member Constructor
		DuctLeak(
			std::string const & Name, // Name of component leak
			Real64 const FlowCoef, // Air Mass Flow Coefficient
			Real64 const FlowExpo // Air Mass Flow exponent
		) :
			Name( Name ),
			FlowCoef( FlowCoef ),
			FlowExpo( FlowExpo )
		{}

	};

	struct EffectiveLeakageRatio // effective leakage ratio component
	{
		// Members
		std::string Name; // Name of component leak
		Real64 ELR; // Value of effective leakage ratio
		Real64 FlowRate; // Maximum airflow rate
		Real64 RefPres; // Reference pressure difference
		Real64 FlowExpo; // Air Mass Flow exponent

		// Default Constructor
    EffectiveLeakageRatio() :
			ELR( 0.0 ),
			FlowRate( 0.0 ),
			RefPres( 0.0 ),
			FlowExpo( 0.0 )
		{}

		// Member Constructor
    EffectiveLeakageRatio(
			std::string const & Name, // Name of component leak
			Real64 const ELR, // Value of effective leakage ratio
			Real64 const FlowRate, // Maximum airflow rate
			Real64 const RefPres, // Reference pressure difference
			Real64 const FlowExpo // Air Mass Flow exponent
		) :
			Name( Name ),
			ELR( ELR ),
			FlowRate( FlowRate ),
			RefPres( RefPres ),
			FlowExpo( FlowExpo )
		{}

	};

	struct Duct // Duct component
	{
		// Members
		std::string Name; // Name of duct component
		Real64 L; // Duct length [m]
		Real64 D; // Hydrolic diameter [m]
		Real64 A; // Cross section area [m2]
		Real64 Rough; // Surface roughness [m]
		Real64 TurDynCoef; // Turbulent dynamic loss coefficient
		Real64 UThermal; // Overall heat transmittance [W/m2.K]
		Real64 UMoisture; // Overall moisture transmittance [kg/m2]
		Real64 MThermal; // Thermal capacity [J/K]
		Real64 MMoisture; // Mositure capacity [kg]
		Real64 LamDynCoef; // Laminar dynamic loss coefficient
		Real64 LamFriCoef; // Laminar friction loss coefficient
		Real64 InitLamCoef; // Coefficient of linear initialization
		Real64 RelRough; // e/D: relative roughness,
		Real64 RelL; // L/D: relative length,
		Real64 g; // 1/sqrt(Darcy friction factor),
		Real64 A1; // 1.14 - 0.868589*ln(e/D),

		// Default Constructor
		Duct() :
			L( 0.0 ),
			D( 0.0 ),
			A( 0.0 ),
			Rough( 0.0 ),
			TurDynCoef( 0.0 ),
			UThermal( 0.0 ),
			UMoisture( 0.0 ),
			MThermal( 0.0 ),
			MMoisture( 0.0 ),
			LamDynCoef( 0.0 ),
			LamFriCoef( 0.0 ),
			InitLamCoef( 0.0 ),
			RelRough( 0.0 ),
			RelL( 0.0 ),
			g( 0.0 ),
			A1( 0.0 )
		{}

		// Member Constructor
		Duct(
			std::string const & Name, // Name of duct component
			Real64 const L, // Duct length [m]
			Real64 const D, // Hydrolic diameter [m]
			Real64 const A, // Cross section area [m2]
			Real64 const Rough, // Surface roughness [m]
			Real64 const TurDynCoef, // Turbulent dynamic loss coefficient
			Real64 const UThermal, // Overall heat transmittance [W/m2.K]
			Real64 const UMoisture, // Overall moisture transmittance [kg/m2]
			Real64 const MThermal, // Thermal capacity [J/K]
			Real64 const MMoisture, // Mositure capacity [kg]
			Real64 const LamDynCoef, // Laminar dynamic loss coefficient
			Real64 const LamFriCoef, // Laminar friction loss coefficient
			Real64 const InitLamCoef, // Coefficient of linear initialization
			Real64 const RelRough, // e/D: relative roughness,
			Real64 const RelL, // L/D: relative length,
			Real64 const g, // 1/sqrt(Darcy friction factor),
			Real64 const A1 // 1.14 - 0.868589*ln(e/D),
		) :
			Name( Name ),
			L( L ),
			D( D ),
			A( A ),
			Rough( Rough ),
			TurDynCoef( TurDynCoef ),
			UThermal( UThermal ),
			UMoisture( UMoisture ),
			MThermal( MThermal ),
			MMoisture( MMoisture ),
			LamDynCoef( LamDynCoef ),
			LamFriCoef( LamFriCoef ),
			InitLamCoef( InitLamCoef ),
			RelRough( RelRough ),
			RelL( RelL ),
			g( g ),
			A1( A1 )
		{}

	};

	struct Damper // Damper component
	{
		// Members
		std::string Name; // Name of damper component
		Real64 LTP; // Value for laminar turbulent transition
		Real64 LamFlow; // Laminar flow coefficient
		Real64 TurFlow; // Turbulent flow coefficient
		Real64 FlowExpo; // Air Mass Flow exponent
		Real64 FlowMin; // Minimum control air mass rate
		Real64 FlowMax; // Maximum control air mass rate
		Real64 A0; // First polynomial coefficient of the control variable (constant coefficient)
		Real64 A1; // Second polynomial coefficient of the control variable (linear coefficient)
		Real64 A2; // Third polynomial coefficient of the control variable (quadratic coefficient)
		Real64 A3; // Fourth polynomial coefficient of the control variable (cubic coefficient)

		// Default Constructor
		Damper() :
			LTP( 0.0 ),
			LamFlow( 0.0 ),
			TurFlow( 0.0 ),
			FlowExpo( 0.0 ),
			FlowMin( 0.0 ),
			FlowMax( 0.0 ),
			A0( 0.0 ),
			A1( 0.0 ),
			A2( 0.0 ),
			A3( 0.0 )
		{}

		// Member Constructor
		Damper(
			std::string const & Name, // Name of damper component
			Real64 const LTP, // Value for laminar turbulent transition
			Real64 const LamFlow, // Laminar flow coefficient
			Real64 const TurFlow, // Turbulent flow coefficient
			Real64 const FlowExpo, // Air Mass Flow exponent
			Real64 const FlowMin, // Minimum control air mass rate
			Real64 const FlowMax, // Maximum control air mass rate
			Real64 const A0, // First polynomial coefficient of the control variable (constant coefficient)
			Real64 const A1, // Second polynomial coefficient of the control variable (linear coefficient)
			Real64 const A2, // Third polynomial coefficient of the control variable (quadratic coefficient)
			Real64 const A3 // Fourth polynomial coefficient of the control variable (cubic coefficient)
		) :
			Name( Name ),
			LTP( LTP ),
			LamFlow( LamFlow ),
			TurFlow( TurFlow ),
			FlowExpo( FlowExpo ),
			FlowMin( FlowMin ),
			FlowMax( FlowMax ),
			A0( A0 ),
			A1( A1 ),
			A2( A2 ),
			A3( A3 )
		{}

	};

	struct ConstantVolumeFan // Constant volume fan component
	{
		// Members
		std::string Name; // Name of detailed fan component
		Real64 FlowRate; // Air volume flow rate
		Real64 Ctrl; // Control ratio
		int FanTypeNum; // Fan type: Constant volume or ONOFF
		int FanIndex; // Fan index
		int InletNode; // Inlet node number
		int OutletNode; // Outlet node number
		Real64 MaxAirMassFlowRate; // Max Specified MAss Flow Rate of Damper [kg/s]

		// Default Constructor
    ConstantVolumeFan() :
			FlowRate( 0.0 ),
			Ctrl( 0.0 ),
			FanTypeNum( 0 ),
			FanIndex( 0 ),
			InletNode( 0 ),
			OutletNode( 0 ),
			MaxAirMassFlowRate( 0.0 )
		{}

		// Member Constructor
    ConstantVolumeFan(
			std::string const & Name, // Name of detailed fan component
			Real64 const FlowRate, // Air volume flow rate
			Real64 const Ctrl, // Control ratio
			int const FanTypeNum, // Fan type: Constant volume or ONOFF
			int const FanIndex, // Fan index
			int const InletNode, // Inlet node number
			int const OutletNode, // Outlet node number
			Real64 const MaxAirMassFlowRate // Max Specified MAss Flow Rate of Damper [kg/s]
		) :
			Name( Name ),
			FlowRate( FlowRate ),
			Ctrl( Ctrl ),
			FanTypeNum( FanTypeNum ),
			FanIndex( FanIndex ),
			InletNode( InletNode ),
			OutletNode( OutletNode ),
			MaxAirMassFlowRate( MaxAirMassFlowRate )
		{}

	};

	struct DetailedFan // Detailed fan component
	{
		// Members
		std::string Name; // Name of constant volume fan component
		Real64 FlowCoef; // Coefficient for linear initialization
		Real64 FlowExpo; // Turbulent flow coefficient
		Real64 RhoAir; // Reference air density
		Real64 Qfree; // Free delivery flow at P=0
		Real64 Pshut; // Shutoff pressure at Q=0
		Real64 TranRat; // Flow coefficient at laminar/turbulent transition
		int n; // Number of ranges for fan performance curve
		Array1D< Real64 > Coeff; // Coefficients of fan performance curve.
		//Each range has a min flow rate and 4 coeffieincts

		// Default Constructor
    DetailedFan() :
			FlowCoef( 0.0 ),
			FlowExpo( 0.0 ),
			RhoAir( 0.0 ),
			Qfree( 0.0 ),
			Pshut( 0.0 ),
			TranRat( 0.0 )
		{}

		// Member Constructor
    DetailedFan(
			std::string const & Name, // Name of constant volume fan component
			Real64 const FlowCoef, // Coefficient for linear initialization
			Real64 const FlowExpo, // Turbulent flow coefficient
			Real64 const RhoAir, // Reference air density
			Real64 const Qfree, // Free delivery flow at P=0
			Real64 const Pshut, // Shutoff pressure at Q=0
			Real64 const TranRat, // Flow coefficient at laminar/turbulent transition
			int const n, // Number of ranges for fan performance curve
			Array1< Real64 > const & Coeff // Coefficients of fan performance curve.
		) :
			Name( Name ),
			FlowCoef( FlowCoef ),
			FlowExpo( FlowExpo ),
			RhoAir( RhoAir ),
			Qfree( Qfree ),
			Pshut( Pshut ),
			TranRat( TranRat ),
			n( n ),
			Coeff( Coeff )
		{}

	};

	struct Coil // Coil component
	{
		// Members
		std::string Name; // Name of coil component
		std::string EPlusType; // EnergyPlus coil type
		Real64 L; // Air path length
		Real64 D; // Air path hydraulic diameter

		// Default Constructor
		Coil() :
			L( 0.0 ),
			D( 0.0 )
		{}

		// Member Constructor
		Coil(
			std::string const & Name, // Name of coil component
			std::string const & EPlusType, // EnergyPlus coil type
			Real64 const L, // Air path length
			Real64 const D // Air path hydraulic diameter
		) :
			Name( Name ),
			EPlusType( EPlusType ),
			L( L ),
			D( D )
		{}

	};

	struct DisSysCompHXProp // Coil component
	{
		// Members
		std::string Name; // Name of coil component
		std::string EPlusType; // EnergyPlus coil type
		Real64 L; // Air path length
		Real64 D; // Air path hydraulic diameter
		bool CoilParentExists; // Is a coil component

		// Default Constructor
		DisSysCompHXProp() :
			L( 0.0 ),
			D( 0.0 ),
			CoilParentExists( false )
		{}

		// Member Constructor
		DisSysCompHXProp(
			std::string const & Name, // Name of coil component
			std::string const & EPlusType, // EnergyPlus coil type
			Real64 const L, // Air path length
			Real64 const D, // Air path hydraulic diameter
			bool const CoilParentExists // Is a coil component
		) :
			Name( Name ),
			EPlusType( EPlusType ),
			L( L ),
			D( D ),
			CoilParentExists( CoilParentExists )
		{}

	};

	struct TerminalUnit // Terminal unit component
	{
		// Members
		std::string Name; // Name of coil component
		std::string EPlusType; // EnergyPlus coil type
		Real64 L; // Air path length
		Real64 D; // Air path hydraulic diameter
		int DamperInletNode; // Damper inlet node number
		int DamperOutletNode; // Damper outlet node number

		// Default Constructor
    TerminalUnit() :
			L( 0.0 ),
			D( 0.0 ),
			DamperInletNode( 0 ),
			DamperOutletNode( 0 )
		{}

		// Member Constructor
    TerminalUnit(
			std::string const & Name, // Name of coil component
			std::string const & EPlusType, // EnergyPlus coil type
			Real64 const L, // Air path length
			Real64 const D, // Air path hydraulic diameter
			int const DamperInletNode, // Damper inlet node number
			int const DamperOutletNode // Damper outlet node number
		) :
			Name( Name ),
			EPlusType( EPlusType ),
			L( L ),
			D( D ),
			DamperInletNode( DamperInletNode ),
			DamperOutletNode( DamperOutletNode )
		{}

	};

	struct ConstantPressureDrop // Constant pressure drop component
	{
		// Members
		std::string Name; // Name of constant pressure drop component
		Real64 A; // cross section area
		Real64 DP; // Pressure difference across the component

		// Default Constructor
    ConstantPressureDrop() :
			A( 0.0 ),
			DP( 0.0 )
		{}

		// Member Constructor
    ConstantPressureDrop(
			std::string const & Name, // Name of constant pressure drop component
			Real64 const A, // cross section area
			Real64 const DP // Pressure difference across the component
		) :
			Name( Name ),
			A( A ),
			DP( DP )
		{}

	};

} // AirflowNetwork

} // EnergyPlus

#endif
