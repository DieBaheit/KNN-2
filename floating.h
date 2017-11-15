#ifndef ANN_FLOATING_H
#define ANN_FLOATING_H

#include <cmath> // supplies abs, fabs and max.
#include <cstdint> // supplies fixed width integer types.
#include <limits> // supplies numeric_limits (including epsilon)
#include <type_traits> // supplies information on types.

namespace ann
{
  /// Tolerable distance in Units in the Last Place for comparing floating
  /// point numbers.
  constexpr std::uint_fast8_t ULP_TOLERANCE = 10;
  /// Product of this and corresponding epsilon yields tolerance for comparing
  /// floating point numbers.
  constexpr std::uint_fast8_t ABSOLUTE_TOLERANCE_FACTOR = 10 * ULP_TOLERANCE;


  /// Supplies typedefs to integer types with the same bit width as T. Types
  /// resolve to void, if no integer type of matching width exists.
  template <class T> struct SameWidthAs
  {
    typedef // signed integer
      typename
      std::conditional< sizeof(T) == 1, std::int8_t,
        typename
        std::conditional< sizeof(T) == 2, std::int16_t,
          typename
          std::conditional< sizeof(T) == 4, std::int32_t,
            typename
            std::conditional< sizeof(T) == 8, std::int64_t, bool >::type >::type >::type >::type
      Int;

    typedef // unsigned integer
      typename
      std::conditional< sizeof(T) == 1, std::uint8_t, 
        typename
        std::conditional< sizeof(T) == 2, std::uint16_t, 
          typename
          std::conditional< sizeof(T) == 4, std::uint32_t, 
            typename
            std::conditional< sizeof(T) == 8, std::uint64_t, bool >::type >::type >::type >::type
      UInt;
  };


  /// Holds information on any floating point type T; undefined for other types.
  /// Only use the single template argument version of this class.
  template <class T, class Enable = void> struct FloatingFacts;
  /// Specialization for floating point types.
  template <class T>
  struct
  FloatingFacts
  < T,
    typename std::enable_if< std::is_floating_point< T >::value >::type
  >
  {
    /// Unsigned integral type used to specify bit width of type T.
    typedef std::uint_fast16_t BitWidth;
    /// Shortcut to additional information on type T.
    typedef std::numeric_limits<T> Traits;

    /// Bit width of type T.
    static constexpr BitWidth BITS = sizeof(T) * 8;
    /// Bit width of the mantissa of type T.
    static constexpr BitWidth MANTISSA_BITS = Traits::digits - 1;
    /// Bit width of the exponent of type T.
    static constexpr BitWidth EXPONENT_BITS = BITS - MANTISSA_BITS - 1;

    /// True, iff this type is compatible with ULP based relative comparison.
    static constexpr bool can_exploit_ulp
      = std::numeric_limits< T >::is_iec559 &&
          !std::is_void< typename SameWidthAs< T >::Int >::value;
  };



  /// Compares floating point values in a rather reasonable way. This is
  /// probably as good as a generic comparison gets, but will still fail in
  /// many situations. You should always consider not only the magnitude of
  /// the compared values (what this function does), but also pay heed to the
  /// history of each value to know what to consider equal.
  /// This version utilizes ULP based relative comparison. It is only available
  /// for floating point types of applicable format and bit width.
  /// Thanks to: http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
  template <typename T>
  typename std::enable_if
           < std::is_floating_point< T >::value && FloatingFacts< T >::can_exploit_ulp,
             bool
           >::type
  equal
  (
    T xA,
    T yA,
    T absoluteToleranceA
      = ABSOLUTE_TOLERANCE_FACTOR * std::numeric_limits< T >::epsilon(),
    std::uint_fast8_t ulpToleranceA = ULP_TOLERANCE
  )
  {
    // Have an integal type of appropriate width
    typedef typename SameWidthAs<T>::Int SameWidthInt;
    // Absolute comparison needed when x and y are close to 0. Also covers ==.
    if (std::fabs(xA - yA) <= absoluteToleranceA)  return true;
    // Signs need to match (except for 0, which is covered above).
    if ((xA < 0) != (yA < 0))  return false;
    // ULP based relative comparison
    return ulpToleranceA
      >= std::abs
         ( reinterpret_cast<SameWidthInt&>(xA) -
           reinterpret_cast<SameWidthInt&>(yA)
         );
  }


  /// This version of equality comparison is used for types, that can not be
  /// compared with the ULP trick. It's not inferior to the ULP version, but
  /// behaves slightly different and looks less cool.
  template <typename T>
  typename std::enable_if
           < std::is_floating_point< T >::value && !FloatingFacts< T >::can_exploit_ulp,
             bool
           >::type
  equal
  (
    T xA,
    T yA,
    T absoluteToleranceA
      = ABSOLUTE_TOLERANCE_FACTOR * std::numeric_limits<T>::epsilon(),
    T relativeToleranceA = ULP_TOLERANCE * std::numeric_limits<T>::epsilon()
  )
  {
    T differenceL = std::fabs(xA - yA);
    // Absolute comparison needed when x and y are close to 0. Also covers ==.
    if (differenceL <= absoluteToleranceA)  return true;
    // The difference could still be relatively small.
    if (differenceL <= relativeToleranceA * std::max(std::fabs(xA), std::fabs(yA)))
      return true;
    // Nope, those differ.
    return false;
  }


  /// Binary equality predicate used for non floating point types. This is not
  /// meant to replace ==, but to offer floating point comparison of templated
  /// types without 'damaging' comparison of other types.
  /// By the way: Please notice the nice overloading based on type properties.
  template <typename T>
  typename std::enable_if<!std::is_floating_point<T>::value, bool>::type
  equal
  (
    T xA,
    T yA
  )
  {
    return xA == yA;
  }

  /// This is a shorcut to function equal with (not always sane) default parameters.
  template <typename T> bool eq (T xA, T yA) { return equal(xA, yA); }

} // namespace ann

#endif //ndef ANN_FLOATING_H
