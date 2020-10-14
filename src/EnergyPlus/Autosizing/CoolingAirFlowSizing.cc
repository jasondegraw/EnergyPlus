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

#include <EnergyPlus/Autosizing/CoolingAirFlowSizing.hh>
#include <EnergyPlus/Data/EnergyPlusData.hh>
#include <EnergyPlus/DataEnvironment.hh>
#include <EnergyPlus/DataHVACGlobals.hh>
#include <EnergyPlus/General.hh>
#include <EnergyPlus/OutputReportPredefined.hh>
#include <EnergyPlus/UtilityRoutines.hh>
#include <EnergyPlus/WeatherManager.hh>

namespace EnergyPlus {

Real64 CoolingAirFlowSizer::size(EnergyPlusData &state, Real64 _originalValue, bool &errorsFound)
{
    if (!this->checkInitialized(errorsFound)) {
        return 0.0;
    }
    this->preSize(_originalValue);
    std::string DDNameFanPeak = "";
    std::string dateTimeFanPeak = "";

    if (this->dataEMSOverrideON) {
        this->autoSizedValue = this->dataEMSOverride;
    } else if (this->dataConstantUsedForSizing > 0 && this->dataFractionUsedForSizing > 0) {
        this->autoSizedValue = this->dataConstantUsedForSizing * this->dataFractionUsedForSizing;
    } else {
        if (this->curZoneEqNum > 0) {
            if (!this->wasAutoSized && !this->sizingDesRunThisZone) {
                this->autoSizedValue = _originalValue;
                if (UtilityRoutines::SameString(this->compType, "Coil:Cooling:DX:TwoStageWithHumidityControlMode")) {
                    this->autoSizedValue /= (1.0 - this->dataBypassFrac); // back out bypass fraction applied in GetInput
                    this->originalValue /= (1.0 - this->dataBypassFrac);  // back out bypass fraction applied in GetInput
                }
            } else if (this->zoneEqSizing(this->curZoneEqNum).DesignSizeFromParent) {
                this->autoSizedValue = this->zoneEqSizing(this->curZoneEqNum).AirVolFlow;
            } else {
                auto const SELECT_CASE_var(this->zoneEqSizing(this->curZoneEqNum).SizingMethod(DataHVACGlobals::CoolingAirflowSizing));
                if ((SELECT_CASE_var == DataSizing::SupplyAirFlowRate) || (SELECT_CASE_var == DataSizing::None) ||
                    (SELECT_CASE_var == DataSizing::FlowPerFloorArea)) {
                    if (this->zoneEqSizing(this->curZoneEqNum).SystemAirFlow) {
                        this->autoSizedValue = max(this->zoneEqSizing(this->curZoneEqNum).AirVolFlow,
                                                   this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow,
                                                   this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow);
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else {
                        if (DataSizing::ZoneCoolingOnlyFan) {
                            this->autoSizedValue = this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow;
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (DataSizing::ZoneHeatingOnlyFan) {
                            this->autoSizedValue = this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow;
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        } else if (this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow) {
                            this->autoSizedValue = this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow;
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                            this->autoSizedValue = this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow;
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                            this->autoSizedValue = max(this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow,
                                                       this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow);
                            if (this->autoSizedValue == this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow) {
                                if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                    this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                    DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                    dateTimeFanPeak =
                                        General::TrimSigDigits(
                                            state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                        "/" +
                                        General::TrimSigDigits(
                                            state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                        " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                                }
                            } else if (this->autoSizedValue == this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow) {
                                if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                    this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                    DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                    dateTimeFanPeak =
                                        General::TrimSigDigits(
                                            state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                        "/" +
                                        General::TrimSigDigits(
                                            state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                        " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                                }
                            }
                        } else {
                            this->autoSizedValue = max(this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow,
                                                       this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow);
                            if (this->autoSizedValue == this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow) {
                                if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                    this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                    DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                    dateTimeFanPeak =
                                        General::TrimSigDigits(
                                            state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                        "/" +
                                        General::TrimSigDigits(
                                            state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                        " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                                }
                            } else if (this->autoSizedValue == this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow) {
                                if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                    this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                    DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                    dateTimeFanPeak =
                                        General::TrimSigDigits(
                                            state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                        "/" +
                                        General::TrimSigDigits(
                                            state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                        " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                                }
                            }
                        }
                    }
                } else if (SELECT_CASE_var == DataSizing::FractionOfAutosizedCoolingAirflow) {
                    if (DataSizing::ZoneCoolingOnlyFan) {
                        this->autoSizedValue = this->dataFracOfAutosizedCoolingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow;
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else if (DataSizing::ZoneHeatingOnlyFan) {
                        this->autoSizedValue = this->dataFracOfAutosizedHeatingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow;
                        if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow) {
                        this->autoSizedValue = this->dataFracOfAutosizedCoolingAirflow * this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow;
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue = this->dataFracOfAutosizedHeatingAirflow * this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow;
                        if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue =
                            max(this->dataFracOfAutosizedCoolingAirflow * this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow,
                                this->dataFracOfAutosizedHeatingAirflow * this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow);
                        if (this->autoSizedValue ==
                            this->dataFracOfAutosizedCoolingAirflow * this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow) {
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->autoSizedValue ==
                                   this->dataFracOfAutosizedHeatingAirflow * this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow) {
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        }
                    } else {
                        this->autoSizedValue =
                            max(this->dataFracOfAutosizedCoolingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow,
                                this->dataFracOfAutosizedHeatingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow);
                        if (this->autoSizedValue ==
                            this->dataFracOfAutosizedCoolingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow) {
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->autoSizedValue ==
                                   this->dataFracOfAutosizedHeatingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow) {
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        }
                    }
                } else if (SELECT_CASE_var == DataSizing::FractionOfAutosizedHeatingAirflow) {
                    if (DataSizing::ZoneCoolingOnlyFan) {
                        this->autoSizedValue = this->dataFracOfAutosizedCoolingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow;
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else if (DataSizing::ZoneHeatingOnlyFan) {
                        this->autoSizedValue = this->dataFracOfAutosizedHeatingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow;
                        if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow) {
                        this->autoSizedValue = this->dataFracOfAutosizedCoolingAirflow * this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow;
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue = this->dataFracOfAutosizedHeatingAirflow * this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow;
                        if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue =
                            max(this->dataFracOfAutosizedCoolingAirflow * this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow,
                                this->dataFracOfAutosizedHeatingAirflow * this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow);
                        if (this->autoSizedValue ==
                            this->dataFracOfAutosizedCoolingAirflow * this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow) {
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->autoSizedValue ==
                                   this->dataFracOfAutosizedHeatingAirflow * this->zoneEqSizing(this->curZoneEqNum).HeatingAirVolFlow) {
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        }
                    } else {
                        this->autoSizedValue =
                            max(this->dataFracOfAutosizedCoolingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow,
                                this->dataFracOfAutosizedHeatingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow);
                        if (this->autoSizedValue ==
                            this->dataFracOfAutosizedCoolingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow) {
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->autoSizedValue ==
                                   this->dataFracOfAutosizedHeatingAirflow * this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow) {
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        }
                    }
                } else if (SELECT_CASE_var == DataSizing::FlowPerCoolingCapacity) {
                    if (DataSizing::ZoneCoolingOnlyFan) {
                        this->autoSizedValue = this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity;
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else if (DataSizing::ZoneHeatingOnlyFan) {
                        this->autoSizedValue = this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity;
                        if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow) {
                        this->autoSizedValue = this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity;
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue = this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity;
                        if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue = max(this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity,
                                                   this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity);
                        if (this->autoSizedValue == this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity) {
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->autoSizedValue == this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity) {
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        }
                    } else {
                        this->autoSizedValue = max(this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity,
                                                   this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity);
                        if (this->autoSizedValue == this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity) {
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->autoSizedValue == this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity) {
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        }
                    }
                } else if (SELECT_CASE_var == DataSizing::FlowPerHeatingCapacity) {
                    if (DataSizing::ZoneCoolingOnlyFan) {
                        this->autoSizedValue = this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity;
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else if (DataSizing::ZoneHeatingOnlyFan) {
                        this->autoSizedValue = this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity;
                        if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow) {
                        this->autoSizedValue = this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity;
                        if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && !this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue = this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity;
                        if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                            this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                            DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                            dateTimeFanPeak =
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) + "/" +
                                General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                        }
                    } else if (this->zoneEqSizing(this->curZoneEqNum).HeatingAirFlow && this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue = max(this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity,
                                                   this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity);
                        if (this->autoSizedValue == this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity) {
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->autoSizedValue == this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity) {
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        }
                    } else {
                        this->autoSizedValue = max(this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity,
                                                   this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity);
                        if (this->autoSizedValue == this->dataFlowPerCoolingCapacity * this->dataAutosizedCoolingCapacity) {
                            if (this->finalZoneSizing(this->curZoneEqNum).CoolDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).CoolDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).CoolDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtCoolMax);
                            }
                        } else if (this->autoSizedValue == this->dataFlowPerHeatingCapacity * this->dataAutosizedHeatingCapacity) {
                            if (this->finalZoneSizing(this->curZoneEqNum).HeatDDNum > 0 &&
                                this->finalZoneSizing(this->curZoneEqNum).HeatDDNum <= DataEnvironment::TotDesDays) {
                                DDNameFanPeak = state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Title;
                                dateTimeFanPeak =
                                    General::TrimSigDigits(state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).Month) +
                                    "/" +
                                    General::TrimSigDigits(
                                        state.dataWeatherManager->DesDayInput(this->finalZoneSizing(this->curZoneEqNum).HeatDDNum).DayOfMonth) +
                                    " " + coilSelectionReportObj->getTimeText(this->finalZoneSizing(this->curZoneEqNum).TimeStepNumAtHeatMax);
                            }
                        }
                    }
                } else {
                    if (DataSizing::ZoneCoolingOnlyFan) {
                        this->autoSizedValue = this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow;
                    } else if (this->termUnitIU && (this->curTermUnitSizingNum > 0)) {
                        this->autoSizedValue = this->termUnitSizing(this->curTermUnitSizingNum).AirVolFlow;
                    } else if (this->zoneEqFanCoil) {
                        this->autoSizedValue = this->zoneEqSizing(this->curZoneEqNum).AirVolFlow;
                    } else if (this->zoneHeatingOnlyFan) {
                        this->autoSizedValue = this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow;
                    } else if (this->zoneEqSizing(this->curZoneEqNum).CoolingAirFlow) {
                        this->autoSizedValue = this->zoneEqSizing(this->curZoneEqNum).CoolingAirVolFlow;
                    } else {
                        this->autoSizedValue =
                            max(this->finalZoneSizing(this->curZoneEqNum).DesCoolVolFlow, this->finalZoneSizing(this->curZoneEqNum).DesHeatVolFlow);
                    }
                }
            }
        } else if (this->curSysNum > 0) {
            if (!this->wasAutoSized && !this->sizingDesRunThisAirSys) {
                this->autoSizedValue = _originalValue;
                if (UtilityRoutines::SameString(this->compType, "Coil:Cooling:DX:TwoStageWithHumidityControlMode")) {
                    this->autoSizedValue /= (1.0 - this->dataBypassFrac); // back out bypass fraction applied in GetInput
                    this->originalValue /= (1.0 - this->dataBypassFrac);  // back out bypass fraction applied in GetInput
                }
            } else {
                if (this->curOASysNum > 0) {
                    if (this->oaSysEqSizing(this->curOASysNum).AirFlow) {
                        // Parent object sets flow rate
                        this->autoSizedValue = this->oaSysEqSizing(this->curOASysNum).AirVolFlow;
                    } else if (this->oaSysEqSizing(this->curOASysNum).CoolingAirFlow) {
                        // Parent object sets flow rate
                        this->autoSizedValue = this->oaSysEqSizing(this->curOASysNum).CoolingAirVolFlow;
                    } else if (outsideAirSys(this->curOASysNum).AirLoopDOASNum > -1) {
                        this->autoSizedValue = this->airloopDOAS[outsideAirSys(this->curOASysNum).AirLoopDOASNum].SizingMassFlow /
                                               DataEnvironment::StdRhoAir;
                    } else {
                        this->autoSizedValue = this->finalSysSizing(this->curSysNum).DesOutAirVolFlow;
                    }
                } else if (this->dataAirFlowUsedForSizing > 0.0) {
                    this->autoSizedValue = this->dataAirFlowUsedForSizing;
                } else {
                    if (this->unitarySysEqSizing(this->curSysNum).AirFlow) {
                        this->autoSizedValue = this->unitarySysEqSizing(this->curSysNum).AirVolFlow;
                    } else if (this->unitarySysEqSizing(this->curSysNum).CoolingAirFlow) {
                        this->autoSizedValue = this->unitarySysEqSizing(this->curSysNum).CoolingAirVolFlow;
                    } else {
                        if (this->curDuctType == DataHVACGlobals::Main) {
                            this->autoSizedValue = this->finalSysSizing(this->curSysNum).DesMainVolFlow;
                        } else if (this->curDuctType == DataHVACGlobals::Cooling) {
                            this->autoSizedValue = this->finalSysSizing(this->curSysNum).DesCoolVolFlow;
                        } else if (this->curDuctType == DataHVACGlobals::Heating) {
                            this->autoSizedValue = this->finalSysSizing(this->curSysNum).DesHeatVolFlow;
                        } else if (this->curDuctType == DataHVACGlobals::Other) {
                            this->autoSizedValue = this->finalSysSizing(this->curSysNum).DesMainVolFlow;
                        } else {
                            this->autoSizedValue = this->finalSysSizing(this->curSysNum).DesMainVolFlow;
                        }
                    }
                }
            }
        } else if (this->dataNonZoneNonAirloopValue > 0) {
            this->autoSizedValue = this->dataNonZoneNonAirloopValue;
        } else if (!this->wasAutoSized) {
            this->autoSizedValue = this->originalValue;
        } else {
            std::string msg = this->callingRoutine + ' ' + this->compType + ' ' + this->compName + ", Developer Error: Component sizing incomplete.";
            ShowSevereError(msg);
            this->addErrorMessage(msg);
            msg = "SizingString = " + this->sizingString + ", SizingResult = " + General::TrimSigDigits(this->autoSizedValue, 1);
            ShowContinueError(msg);
            this->addErrorMessage(msg);
            errorsFound = true;
        }

        // override sizing string
        if (this->overrideSizeString) {
            if (UtilityRoutines::SameString(this->compType, "ZoneHVAC:FourPipeFanCoil")) {
                this->sizingString = "Maximum Supply Air Flow Rate [m3/s]";
                if (this->isEpJSON) this->sizingString = "maximum_supply_air_flow_rate [m3/s]";
            } else if (this->coilType_Num == DataHVACGlobals::CoilDX_CoolingTwoSpeed) {
                if (this->dataDXSpeedNum == 1) { // mode 1 is high speed in DXCoils loop
                    if (this->isEpJSON) {
                        this->sizingString = "high_speed_rated_air_flow_rate [m3/s]";
                    } else {
                        this->sizingString = "High Speed Rated Air Flow Rate [m3/s]";
                    }
                } else if (this->dataDXSpeedNum == 2) {
                    if (this->isEpJSON) {
                        this->sizingString = "low_speed_rated_air_flow_rate [m3/s]";
                    } else {
                        this->sizingString = "Low Speed Rated Air Flow Rate [m3/s]";
                    }
                }
            } else if (this->isEpJSON) {
                this->sizingString = "cooling_supply_air_flow_rate [m3/s]";
            }
        }
        if (this->dataScalableSizingON) {
            if (this->zoneAirFlowSizMethod == DataSizing::SupplyAirFlowRate || this->zoneAirFlowSizMethod == DataSizing::None) {
                this->sizingStringScalable = "(scaled by flow / zone) ";
            } else if (this->zoneAirFlowSizMethod == DataSizing::FlowPerFloorArea) {
                this->sizingStringScalable = "(scaled by flow / area) ";
            } else if (this->zoneAirFlowSizMethod == DataSizing::FractionOfAutosizedCoolingAirflow ||
                       this->zoneAirFlowSizMethod == DataSizing::FractionOfAutosizedHeatingAirflow) {
                this->sizingStringScalable = "(scaled by fractional multiplier) ";
            } else if (this->zoneAirFlowSizMethod == DataSizing::FlowPerCoolingCapacity ||
                       this->zoneAirFlowSizMethod == DataSizing::FlowPerHeatingCapacity) {
                this->sizingStringScalable = "(scaled by flow / capacity) ";
            }
        }
    }
    this->select2StgDXHumCtrlSizerOutput(errorsFound);

    if (this->isCoilReportObject) {
        // SizingResult is airflow in m3/s
        coilSelectionReportObj->setCoilAirFlow(this->compName, this->compType, this->autoSizedValue, this->wasAutoSized);
    }
    if (this->isFanReportObject) {
        //  fill fan peak day and time here
        if (this->dataScalableSizingON) {
            DDNameFanPeak = "Scaled size, not from any peak";
            dateTimeFanPeak = "Scaled size, not from any peak";
        }
        OutputReportPredefined::PreDefTableEntry(OutputReportPredefined::pdchFanDesDay, this->compName, DDNameFanPeak);
        OutputReportPredefined::PreDefTableEntry(OutputReportPredefined::pdchFanPkTime, this->compName, dateTimeFanPeak);
    }
    return this->autoSizedValue;
}

void CoolingAirFlowSizer::clearState()
{
    BaseSizerWithScalableInputs::clearState();
}

} // namespace EnergyPlus
