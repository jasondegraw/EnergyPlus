// EnergyPlus, Copyright (c) 1996-2018, The Board of Trustees of the University of Illinois,
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

#ifndef ErrorTracking_hh_INCLUDED
#define ErrorTracking_hh_INCLUDED

#include <array>
#include <string>
#include <vector>
#include <iostream>

// EnergyPlus Headers
#define NO_USING_OBJEXXFCL
#include <EnergyPlus.hh>
#include <UtilityRoutines.hh>
#undef NO_USING_OBJEXXFCL
#include <ErrorCodes.hh>

namespace EnergyPlus {

namespace ErrorTracking {

    enum class EventType
    {
        Info,
        Debug,
        Warning,
        Error, // An event has occurred that will eventually require termination
        Fatal
    };

    struct LoggedEvent
    {
        LoggedEvent(const EventType type,
                    const std::string &label,
                    const std::string &showFormat,
                    const std::string &storeFormat)
            : type(type), label(label), showFormat(showFormat), storeFormat(storeFormat), count(0)
        {
        }

        const EventType type;                      // Type of event
        const std::string label;                   // A string label for the event
        const std::string showFormat;             // Format format for display to the user
        const std::string storeFormat;            // Format format for storage formats
        int count;                                 // The number of times this event has been logged

    };

    struct TrackedEvent : public LoggedEvent
    {
        enum class Tracking
        {
            None,
            Min,
            Max,
            Sum,
            MinMax,
            MinSum,
            MaxSum,
            MinMaxSum
        };
        TrackedEvent(const EventType type,
                     const std::string &label,
                     const std::string &summaryFormat,
                     const std::string &showFormat,
                     const std::string &storeFormat,
                     const std::vector<Tracking> &tracking = {},
                     const std::array<std::string, 3> labels = {{{}, {}, {}}})
            : LoggedEvent(type, label, showFormat, storeFormat), summaryFormat(summaryFormat), tracking(tracking), min_value(std::numeric_limits<double>::infinity()),
              min_label(labels[0]), max_value(-std::numeric_limits<double>::infinity()), max_label(labels[1]), sum_value(0), sum_label(labels[2])
        {
        }

        template <typename... Args> void update(const Args &... args)
        {
            if (internal_update(tracking.begin(), args...)) {
                ++count;
            }
        }

        const std::string summaryFormat;
        const std::vector<Tracking> tracking;
        double min_value;
        const std::string min_label;
        double max_value;
        const std::string max_label;
        double sum_value;
        const std::string sum_label;

    private:
        bool internal_update(std::vector<Tracking>::const_iterator &iter)
        {
            return false;
        }

        template <typename... Args> bool internal_update(std::vector<Tracking>::const_iterator &iter, int val, const Args &... args)
        {
            return internal_update(iter, static_cast<double> val, args...);
        }

        template <typename... Args> bool internal_update(std::vector<Tracking>::const_iterator &iter, const std::string &val, const Args &... args)
        {
            if (iter != types.end()) {
                ++iter;
                return internal_update(iter, args...);
            }
            return true;
        }

        template <typename... Args> bool internal_update(std::vector<Tracking>::const_iterator &iter, double val, const Args &... args)
        {
            if (iter != types.end()) {
                switch (*iter) {
                case Type::Min:
                    min_value = std::min(min_value, val);
                    break;
                case Type::Max:
                    max_value = std::max(max_value, val);
                    break;
                case Type::Sum:
                    sum_value += val;
                    break;
                case Type::MinMax:
                    min_value = std::min(min_value, val);
                    max_value = std::max(max_value, val);
                    break;
                case Type::MinSum:
                    min_value = std::min(min_value, val);
                    sum_value += val;
                    break;
                case Type::MaxSum:
                    max_value = std::max(max_value, val);
                    sum_value += val;
                    break;
                case Type::MinMaxSum:
                    min_value = std::min(min_value, val);
                    max_value = std::max(max_value, val);
                    sum_value += val;
                    break;
                default:
                    // None is do nothing
                    break;
                }
                ++iter;
                return internal_update(iter, args...);
            }
            return true;
        }

        template <typename... Args> bool internal_update(std::vector<Tracking>::const_iterator &iter, double val)
        {
            if (iter != types.end()) {
                switch (*iter) {
                case Type::Min:
                    min_value = std::min(min_value, val);
                    break;
                case Type::Max:
                    max_value = std::max(max_value, val);
                    break;
                case Type::Sum:
                    sum_value += val;
                    break;
                case Type::MinMax:
                    min_value = std::min(min_value, val);
                    max_value = std::max(max_value, val);
                    break;
                case Type::MinSum:
                    min_value = std::min(min_value, val);
                    sum_value += val;
                    break;
                case Type::MaxSum:
                    max_value = std::max(max_value, val);
                    sum_value += val;
                    break;
                case Type::MinMaxSum:
                    min_value = std::min(min_value, val);
                    max_value = std::max(max_value, val);
                    sum_value += val;
                    break;
                default:
                    // None is do nothing
                    break;
                }
            }
            return true;
        }
    };

    struct RecurringError : public TrackedEvent
    {
        enum class TrackingType
        {
            None,
            Minimum,
            Maximum,
            Sum
        };
        RecurringError(const TrackedEvent &event, int showlimit = std::numeric_limits<int>::infinity())
            : TrackedEvent(event.type, event.label, event.summaryFormat, event.showFormat, event.storeFormat, event.tracking), showlimit(showlimit)
        {
        }
        RecurringError(const EventType type,
                       const std::string &label,
                       const std::string &summaryFormat,
                       const std::string &showFormat,
                       const std::string &storeFormat,
                       const std::vector<Tracking> &types = {},
                       const std::array<std::string, 3> labels = {{{}, {}, {}}},
                       int showlimit = 0)
            : TrackedEvent(type, label, summaryFormat, showFormat, storeFormat, types, labels), showlimit(showlimit)
        {
        }

        template <typename... Args> static RecurringError from_code(RecurringCode code, int showlimit, const Args &... args)
        {
            size_t index = (size_t)code;
            if (index >= RECURRING_COUNT) {
                return RecurringError(EventType::Error, "DEV0000", fmt::sprintf("Internal Error: Index %d out of range", index));
            }
            std::string oneLiner = fmt::vsprintf(events[index].oneLiner, fmt::make_printf_args(args...));
            return RecurringError(events[index].type,
                                  events[index].label,
                                  oneLiner,
                                  events[index].extended,
                                  events[index].parameters,
                                  events[index].types,
                                  {events[index].min_label, events[index].max_label, events[index].sum_label},
                                  showlimit);
        }

        const int showlimit;
        static std::array<TrackedEvent, RECURRING_COUNT> events;
    };

    typedef RecurringError *RecurringHandle;

    class Tracker
    {
    public:
        Tracker() : fatalExceptions(0), errorExceptions(0), warningExceptions(0), fatalBadIndex(0), errorBadIndex(0), warningBadIndex(0)
        {
        }

        template <typename... Args> void warning(WarningCode code, const char *file, int line, const Args &... args)
        {
            size_t index = (size_t)code;
            if (index >= WARNING_COUNT) {
                write(fmt::sprintf(" **  Error  ** [DEV0000] Internal Error: Warning index %d out of range", index));
                ++warningBadIndex;
                return;
            }
            auto &object = warnings[index];
            try {
                write(fmt::vsprintf(object.format, fmt::make_printf_args(args...)));
            } catch (...) {
                write(fmt::sprintf(
                    " ** Warning ** [%s] Exception thrown during warning processing at line %d in %s", errorcodes::warningcodes[index], line, file));
                ++warning_exceptions;
            }
        }

        template <typename... Args> void error(ErrorCode code, const char *file, int line, const Args &... args)
        {
            size_t index = (size_t)code;
            if (index >= ERROR_COUNT) {
                write(fmt::sprintf(" **  Error  ** [DEV0000] Internal Error: Error index %d out of range", index));
                ++errorBadIndex;
                return;
            }

            auto &object = errors[index];
            try {
                write(fmt::vsprintf(object.format, fmt::make_printf_args(args...)));
            } catch (...) {
                write(fmt::sprintf(
                    " **  Error  ** [%s] Exception thrown during error processing at line %d in %s", errorcodes::errorcodes[index], line, file));
                ++error_exceptions;
            }
        }

        /*
        template <typename... Args> RecurringHandle create_recurring_error(RecurringCode code, int showlimit, const Args &... args)
        {
            //stored.emplace_back(RecurringError::from_code(code, showlimit, args...));
            return &(stored.back());
        }

        template <typename... Args> void show_recurring(RecurringHandle &handle, const char *file, int line, const Args &... args)
        {
            try {
                handle->update(args...);
                write(fmt::vsprintf(handle->format, fmt::make_printf_args(args...)));
            } catch (...) {
                write(fmt::sprintf(" **  Error  ** [%s] Exception thrown during error processing at line %d in %s", handle->label, line, file));
            }
        }
        */

        template <typename... Args> void fatal(FatalCode code, const char *file, int line, const Args &... args)
        {
            size_t index = (size_t)code;
            if (index >= FATAL_COUNT) {
                ++fatalBadIndex;
                //ShowFatalError(fmt::sprintf(" **  Error  ** [DEV0000] Internal Error: Error index %d out of range", index));
                ShowFatalError(fmt::sprintf("Internal Error: Error index %d out of range", index));
            }

            auto &object = fatals[index];
            ++object.count;
            std::string message;
            try {
                message = fmt::vsprintf(object.showFormat, fmt::make_printf_args(args...));
            } catch (...) {
                ++fatalExceptions;
                ShowFatalError(
                    fmt::sprintf(
                    //" **  Error  ** [%s] Exception thrown during fatal processing at line %d in %s", object.label, line, file));
                    "Exception thrown during fatal processing at line %d in %s",
                    //object.label,
                    line,
                    file));
            }
            ShowFatalError(message);
        }

        /*void summarize_recurring()
        {
            std::cout << "\n===== Recurring Error Summary =====\nThe following recurring error messages occurred.\n";
            for (auto &s : stored) {
                std::cout << s.summary << std::endl;
                std::cout << fmt::sprintf(" **   ~~~   ** This error occurred %d total times;\n", s.count);
                if (!s.max_label.empty()) {
                    std::cout << fmt::sprintf(" **   ~~~   ** Maximum value of \"%s\": %f;\n", s.max_label, s.max_value);
                }
                if (!s.min_label.empty()) {
                    std::cout << fmt::sprintf(" **   ~~~   ** Minimum value of \"%s\": %f;\n", s.min_label, s.min_value);
                }
                if (!s.sum_label.empty()) {
                    std::cout << fmt::sprintf(" **   ~~~   ** Minimum value of \"%s\": %f;\n", s.sum_label, s.min_value);
                }
            }
        }*/

        std::string diagnostics()
        {
            std::string result("\n===== Error Diagnostic Report =====\nThe following issues occurred during error processing.");
            // write(fmt::sprintf(" %d exceptions during fatal processing", fatal_exceptions));
            // write(fmt::sprintf(" %d exceptions during error processing", error_exceptions));
            // write(fmt::sprintf(" %d exceptions during warning processing", warning_exceptions));
        }

        void write(const std::string &message)
        {
            std::cout << message << std::endl;
        }

        // std::string last_error_message;

        // std::vector<RecurringError> stored;
        int fatalExceptions;
        int errorExceptions;
        int warningExceptions;
        int fatalBadIndex;
        int errorBadIndex;
        int recurringBadIndex;
        int warningBadIndex;

        // static std::array<LoggingEvent, 1> warnings;
        // static std::array<LoggingEvent, 1> errors;
        static std::array<LoggedEvent, FATAL_COUNT> fatals;
    };


} // namespace ErrorTracking


extern ErrorTracking::Tracker tracker;

} // namespace EnergyPlus

#endif
