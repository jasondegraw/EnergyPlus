// EnergyPlus, Copyright (c) 1996-2020, The Board of Trustees of the University of Illinois,
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy), Oak Ridge
// National Laboratory, managed by UT-Battelle, Alliance for Sustainable Energy, LLC, and other
// contributors. All rights reserved.
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
//     similar designation, without the U.S. Department of Energy's prior written consent.
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

#include "AutosizingFixture.hh"
#include <gtest/gtest.h>

#include <EnergyPlus/Autosizing/HeatingAirflowUASizing.hh>
#include <EnergyPlus/DataEnvironment.hh>
#include <EnergyPlus/DataHVACGlobals.hh>
#include <EnergyPlus/DataSizing.hh>
#include <EnergyPlus/WaterCoils.hh>
#include <ObjexxFCL/Array1D.hh>

namespace EnergyPlus {

TEST_F(AutoSizingFixture, HeatingAirflowUASizingGauntlet)
{
    // this global state is what would be set up by E+ currently
    DataEnvironment::StdRhoAir = 1.2;
    // there is definitely a better way to do this...
    Array1D<EnergyPlus::DataSizing::TermUnitSizingData> tmpTermUnitData;
    Array1D<EnergyPlus::DataSizing::ZoneSizingData> tmpFinalZoneSizing;
    Array1D<EnergyPlus::DataSizing::ZoneEqSizingData> tmpZoneEqSizing;
    Array1D<EnergyPlus::DataSizing::SystemSizingData> tmpFinalSysSizing;
    Array1D<EnergyPlus::DataSizing::SystemSizingInputData> tmpSysSizingInput;
    Array1D<DataAirLoop::OutsideAirSysProps> tmpOutsideAirSys;
    Array1D<DataSizing::ZoneEqSizingData> tmpOASysEqSizing;
    std::vector<AirLoopHVACDOAS::AirLoopDOAS> tmpAirloopDOAS;

    CommonFlags baseFlags;
    baseFlags.compType = DataHVACGlobals::cAllCoilTypes(DataHVACGlobals::Coil_HeatingWater);
    baseFlags.compName = "MyWaterCoil";
    baseFlags.curZoneEqNum = 1;

    // create the sizer and set up the flags to specify the sizing configuration
    HeatingAirflowUASizer sizer;
    HeatingAirflowUASizerFlags flags;

    // ZONE EQUIPMENT TESTING
    baseFlags.curTermUnitSizingNum = 1;
    // this isn't an actual object field (internal data used for sizing) or reported to eio, only for unit testing
    flags.sizingString = "Heating Coil Airflow For UA";

    // reset eio stream
    has_eio_output(true);

    // Test #1 - Zone Equipment, no autosizing
    baseFlags.termUnitSingDuct = true;
    Real64 inputValue = 5;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    AutoSizingResultType result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_FALSE(sizer.wasAutoSized);
    EXPECT_NEAR(5.0, sizer.autoSizedValue, 0.01); // hard-sized value
    sizer.autoSizedValue = 0.0;                   // reset for next test

    baseFlags.printWarningFlag = true; // this field isn't reported to eio, only for unit testing
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_FALSE(sizer.wasAutoSized);
    EXPECT_NEAR(5.0, sizer.autoSizedValue, 0.01); // hard-sized value
    sizer.autoSizedValue = 0.0;                   // reset for next test
    baseFlags.printWarningFlag = false;

    std::string eiooutput =
        std::string({"! <Component Sizing Information>, Component Type, Component Name, Input Field Description, Value\n"
                     " Component Sizing Information, Coil:Heating:Water, MyWaterCoil, User-Specified Heating Coil Airflow For UA, 5.00000\n"});

    EXPECT_TRUE(compare_eio_stream(eiooutput, true));

    // now allocate sizing arrays for testing autosized field
    tmpTermUnitData.allocate(1);
    tmpTermUnitData(1).AirVolFlow = 5;
    tmpFinalZoneSizing.allocate(1);
    tmpZoneEqSizing.allocate(1);

    baseFlags.zoneSizingRunDone = true;

    // Test 2 - Zone Equipment, Single Duct TU
    baseFlags.termUnitSingDuct = true;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    // do sizing
    sizer.zoneSizingInput.allocate(1);
    sizer.zoneSizingInput(1).ZoneNum = 1;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    EXPECT_NEAR(5.0, tmpTermUnitData(1).AirVolFlow, 0.01);
    EXPECT_NEAR(1.2, DataEnvironment::StdRhoAir, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 3 - Zone Equipment, Powered Induction TU
    baseFlags.termUnitSingDuct = false;
    baseFlags.termUnitPIU = true;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 4 - Zone Equipment, Induction TU
    baseFlags.termUnitPIU = false;
    baseFlags.termUnitIU = true;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 5 - Zone Equipment, Zone Eq Fan Coil
    baseFlags.termUnitIU = false;
    baseFlags.zoneEqFanCoil = true;
    tmpTermUnitData(1).AirVolFlow = 0.0;
    tmpFinalZoneSizing(1).DesHeatVolFlow = 5.0;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 6 - Zone Equipment, Other Equipment
    baseFlags.zoneEqFanCoil = false;
    baseFlags.otherEqType = true;
    tmpFinalZoneSizing(1).DesHeatVolFlow = 0.0;
    tmpZoneEqSizing(1).AirVolFlow = 5.0;
    tmpZoneEqSizing(1).SystemAirFlow = true;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 7 - Zone Equipment, Other Equipment
    tmpZoneEqSizing(1).AirVolFlow = 0.0;
    tmpZoneEqSizing(1).HeatingAirVolFlow = 5.0;
    tmpZoneEqSizing(1).SystemAirFlow = false;
    tmpZoneEqSizing(1).HeatingAirFlow = true;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 8 - Zone Equipment, Other Equipment
    tmpZoneEqSizing(1).HeatingAirVolFlow = 0.0;
    tmpZoneEqSizing(1).HeatingAirFlow = false;
    tmpFinalZoneSizing(1).DesHeatMassFlow = 5.0;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(5.0, sizer.autoSizedValue, 0.01); // uses a mass flow rate for sizing
    sizer.autoSizedValue = 0.0;                   // reset for next test

    // reset eio stream
    has_eio_output(true);
    eiooutput = "";

    // AIRLOOP EQUIPMENT TESTING - CurDuctType not set, no reporting
    // Test 9 - Airloop Equipment
    baseFlags.curZoneEqNum = 0;
    baseFlags.numZoneSizingInput = 0;
    baseFlags.curTermUnitSizingNum = 0;
    baseFlags.otherEqType = false;
    tmpZoneEqSizing.deallocate();
    tmpZoneEqSizing.deallocate();
    tmpFinalZoneSizing.deallocate();

    baseFlags.curSysNum = 1;
    baseFlags.curSysNum = 1;
    baseFlags.numPrimaryAirSys = 1;
    baseFlags.numSysSizInput = 1;
    baseFlags.sysSizingRunDone = false;
    // start with a hard-sized value as the user input, no system sizing arrays
    inputValue = 5.0;
    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_FALSE(sizer.wasAutoSized);
    EXPECT_NEAR(5.0, sizer.autoSizedValue, 0.01); // hard-sized value
    sizer.autoSizedValue = 0.0;                   // reset for next test
    EXPECT_TRUE(compare_eio_stream(eiooutput, true));

    // Test 10 - Airloop Equipment - CurDuctType not set
    baseFlags.curSysNum = 1;
    baseFlags.numPrimaryAirSys = 1;
    baseFlags.numSysSizInput = 1;
    baseFlags.sysSizingRunDone = true;
    tmpFinalSysSizing.allocate(1);
    tmpSysSizingInput.allocate(1);
    tmpSysSizingInput(1).AirLoopNum = 1;

    tmpFinalSysSizing(1).DesMainVolFlow = 5.0; // CurDuctType not set
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    baseFlags.printWarningFlag = true; // this field isn't reported to eio, only for unit testing

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // <Component Sizing Information> header already reported above (and flag set false). Only coil sizing information reported here.
    eiooutput = std::string({" Component Sizing Information, Coil:Heating:Water, MyWaterCoil, Design Size Heating Coil Airflow For UA, 6.00000\n"});

    EXPECT_TRUE(compare_eio_stream(eiooutput, true));

    // Test 11 - Airloop Equipment - CurDuctType = Main, SysAirMinFlowRat = 0
    DataSizing::CurDuctType = DataHVACGlobals::Main;
    tmpFinalSysSizing(1).DesMainVolFlow = 5.0;
    tmpFinalSysSizing(1).SysAirMinFlowRat = 0.0;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    baseFlags.printWarningFlag = false; // this field isn't reported to eio, only for unit testing
    baseFlags.curDuctType = DataSizing::CurDuctType;

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 12 - Airloop Equipment - CurDuctType = Main, SysAirMinFlowRat = 0.5
    tmpFinalSysSizing(1).SysAirMinFlowRat = 0.5;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(3.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 13 - Airloop Equipment - CurDuctType = Cooling, SysAirMinFlowRat = 0
    DataSizing::CurDuctType = DataHVACGlobals::Cooling;
    tmpFinalSysSizing(1).DesMainVolFlow = 0.0;
    tmpFinalSysSizing(1).DesCoolVolFlow = 5.0;
    tmpFinalSysSizing(1).SysAirMinFlowRat = 0.0;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    baseFlags.curDuctType = DataSizing::CurDuctType;

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 14 - Airloop Equipment - CurDuctType = Cooling, SysAirMinFlowRat = 0.5
    tmpFinalSysSizing(1).SysAirMinFlowRat = 0.5;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(3.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 15 - Airloop Equipment - CurDuctType = Heating, SysAirMinFlowRat doesn't matter
    DataSizing::CurDuctType = DataHVACGlobals::Heating;
    tmpFinalSysSizing(1).DesCoolVolFlow = 0.0;
    tmpFinalSysSizing(1).DesHeatVolFlow = 5.0;
    tmpFinalSysSizing(1).SysAirMinFlowRat = 0.5;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;
    baseFlags.curDuctType = DataSizing::CurDuctType;

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // OUTDOOR AIR SYSTEM EQUIPMENT TESTING
    // Test 16 - Outdoor Air System Equipment, no DOAS air loop
    tmpFinalSysSizing(1).DesHeatVolFlow = 0.0;
    tmpFinalSysSizing(1).DesOutAirVolFlow = 5.0;
    tmpOASysEqSizing.allocate(1);
    tmpOutsideAirSys.allocate(1);
    baseFlags.curOASysNum = 1;
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(6.0, sizer.autoSizedValue, 0.01);
    sizer.autoSizedValue = 0.0; // reset for next test

    // Test 17 - Outdoor Air System Equipment with DOAS system
    tmpFinalSysSizing(1).DesOutAirVolFlow = 0.0;
    tmpOutsideAirSys(1).AirLoopDOASNum = 0;
    AirLoopHVACDOAS::AirLoopDOAS thisDOAS;
    thisDOAS.SizingMassFlow = 5.0;
    tmpAirloopDOAS.push_back(thisDOAS);
    // start with an auto-sized value as the user input
    inputValue = EnergyPlus::DataSizing::AutoSize;

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_TRUE(sizer.wasAutoSized);
    EXPECT_NEAR(5.0, sizer.autoSizedValue, 0.01); // uses a mass flow rate for sizing
    sizer.autoSizedValue = 0.0;                   // reset for next test

    // reset eio stream
    has_eio_output(true);

    // Test 17 - Outdoor Air System Equipment with DOAS system, hard-sized air flow rate
    // start with an auto-sized value as the user input
    inputValue = 5.0;
    tmpAirloopDOAS[0].SizingMassFlow = 3.0;
    baseFlags.printWarningFlag = true; // this field isn't reported to eio, only for unit testing

    // do sizing
    sizer.wasAutoSized = false;
    sizer.setParameters(baseFlags,
                        flags,
                        tmpTermUnitData,
                        tmpFinalZoneSizing,
                        tmpZoneEqSizing,
                        tmpSysSizingInput,
                        tmpFinalSysSizing,
                        tmpOutsideAirSys,
                        tmpOASysEqSizing,
                        tmpAirloopDOAS);
    result = sizer.size(inputValue);
    EXPECT_EQ(AutoSizingResultType::NoError, result);
    EXPECT_FALSE(sizer.wasAutoSized);
    EXPECT_NEAR(5.0, sizer.autoSizedValue, 0.01); // hard-sized value
    sizer.autoSizedValue = 0.0;                   // reset for next test

    // <Component Sizing Information> header already reported above (and flag set false). Only coil sizing information reported here.
    eiooutput =
        std::string({" Component Sizing Information, Coil:Heating:Water, MyWaterCoil, Design Size Heating Coil Airflow For UA, 3.00000\n"
                     " Component Sizing Information, Coil:Heating:Water, MyWaterCoil, User-Specified Heating Coil Airflow For UA, 5.00000\n"});
    EXPECT_TRUE(compare_eio_stream(eiooutput, true));
}

} // namespace EnergyPlus
