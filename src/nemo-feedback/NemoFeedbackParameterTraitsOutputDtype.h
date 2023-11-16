/*
 * (C) Crown Copyright 2023 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/ParameterTraits.h"

namespace nemo_feedback {

/// Enum type for obs variable data types
enum class OutputDtype {
    Float,
    Double
};

/// Helps with the conversion of OutputDtype values to/from strings.
struct OutputDtypeParameterTraitsHelper {
  typedef OutputDtype EnumType;
  static constexpr char enumTypeName[] = "OutputDtype";
  static constexpr util::NamedEnumerator<EnumType> namedValues[] = {
    { EnumType::Float, "float" },
    { EnumType::Double, "double" }
  };
};

}  // namespace nemo_feedback

namespace oops {

/// Specialization of ParameterTraits for OutputDtype.
template <>
struct ParameterTraits<nemo_feedback::OutputDtype> :
    public EnumParameterTraits<nemo_feedback::OutputDtypeParameterTraitsHelper>
{};

}  // namespace oops
