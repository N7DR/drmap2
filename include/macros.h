// $Id: macros.h 13 2018-12-15 00:34:28Z n7dr $

// Released under the GNU Public License, version 2
//   see: https://www.gnu.org/licenses/gpl-2.0.html

// Principal author: N7DR

// Copyright owners:
//    N7DR

#ifndef MACROS_H
#define MACROS_H

/*! \file   macros.h

    Macros and templates.
*/

#include <algorithm>
#include <chrono>
#include <deque>
#include <experimental/functional>  // for not_fn
#include <future>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <ranges>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cmath>

// syntactic suger for time-related use
using centiseconds =  std::chrono::duration<long, std::centi>;           ///< hundredths of a second
using deciseconds  =  std::chrono::duration<long, std::deci>;            ///< tenths of a second

namespace SR  = std::ranges;
namespace SRV = std::ranges::views;

/// Syntactic sugar for read/write access
// See: "C++ Move Semantics", p. 81; 5.1.3 Using Move Semantics to Solve the Dilemma

#if (!defined(READ_AND_WRITE))

#define READ_AND_WRITE_RET(y)                                       \
/*! Read access to _##y */                                      \
  inline const decltype(_##y)& y(void) const& { return _##y; }    \
  inline decltype(_##y) y(void) && { return std::move(_##y); }    \
\
/*! Write access to _##y */                                     \
  inline auto y(const decltype(_##y)& n) -> decltype(*this)& { _##y = n; return (*this); }

#define READ_AND_WRITE(y)                                       \
/*! Read access to _##y */                                      \
  inline const decltype(_##y)& y(void) const& { return _##y; }    \
  inline decltype(_##y) y(void) && { return std::move(_##y); }    \
\
/*! Write access to _##y */                                     \
  inline void y(const decltype(_##y)& n) { _##y = n; }

#endif    // !READ_AND_WRITE

/// Syntactic sugar for read/write access with thread locking
#if (!defined(SAFE_READ_AND_WRITE))

#define SAFE_READ_AND_WRITE(y, z)                                           \
/*! Read access to _##y */                                                  \
  inline const decltype(_##y)& y(void) const& { SAFELOCK(z); return _##y; }    \
  inline decltype(_##y) y(void) && { SAFELOCK(z); return std::move(_##y); }    \
/*! Write access to _##y */                                                 \
  inline void y(const decltype(_##y)& n) { SAFELOCK(z); _##y = n; }

#endif    // !SAFE_READ_AND_WRITE

/// Read and write with mutex that is part of the object
#if (!defined(SAFE_READ_AND_WRITE_WITH_INTERNAL_MUTEX))

#define SAFE_READ_AND_WRITE_WITH_INTERNAL_MUTEX(y, z)                                           \
/*! Read access to _##y */                                                  \
  inline const decltype(_##y)& y(void) const& { SAFELOCK(z); return _##y; }    \
  inline decltype(_##y) y(void) && { SAFELOCK(z); return std::move(_##y); }    \
/*! Write access to _##y */                                                 \
  inline void y(const decltype(_##y)& n) { SAFELOCK(z); _##y = n; }

#endif    // !SAFE_READ_AND_WRITE_WITH_INTERNAL_MUTEX

/// Syntactic sugar for read-only access
#if (!defined(READ))

#define READ(y)                                                 \
/*! Read-only access to _##y */                                 \
inline const decltype(_##y)& y(void) const& { return _##y; }    \
inline decltype(_##y) y(void) && { return std::move(_##y); } 

#endif    // !READ

/// Syntactic sugar for read-only access with thread locking
#if (!defined(SAFEREAD))

#define SAFEREAD(y, z)                                                      \
/*! Read-only access to _##y */                                             \
  inline const decltype(_##y)& y(void) const& { SAFELOCK(z); return _##y; }    \
  inline decltype(_##y) y(void) && { SAFELOCK(z); return std::move(_##y); } 
#endif    // !SAFEREAD

// alternative name for SAFEREAD
#if (!defined(SAFE_READ))

/// read with a mutex
#define SAFE_READ(y, z)                                                      \
/*! Read-only access to _##y */                                             \
  inline const decltype(_##y)& y(void) const& { SAFELOCK(z); return _##y; }    \
  inline decltype(_##y) y(void) && { SAFELOCK(z); return std::move(_##y); }

#endif    // !SAFE_READ

/// read with mutex that is part of the object
#if (!defined(SAFEREAD_WITH_INTERNAL_MUTEX))

#define SAFEREAD_WITH_INTERNAL_MUTEX(y, z)                                                      \
/*! Read-only access to _##y */                                             \
  inline const decltype(_##y)& y(void) const& { SAFELOCK(z); return _##y; }    \
  inline decltype(_##y) y(void) && { SAFELOCK(z); return std::move(_##y); }

#endif    // !SAFEREAD_WITH_INTERNAL_MUTEX

// alternative name for SAFEREAD_WITH_INTERNAL_MUTEX
#if (!defined(SAFE_READ_WITH_INTERNAL_MUTEX))

/// read with an internal mutex
#define SAFE_READ_WITH_INTERNAL_MUTEX(y, z)                                                      \
/*! Read-only access to _##y */                                             \
  inline const decltype(_##y)& y(void) const& { SAFELOCK(z); return _##y; }    \
  inline decltype(_##y) y(void) && { SAFELOCK(z); return std::move(_##y); }

#endif    // !SAFE_READ_WITH_INTERNAL_MUTEX

/// Syntactic sugar for write access
#if (!defined(WRITE))

#define WRITE(y)                                       \
/*! Write access to _##y */                                     \
  inline void y(const decltype(_##y)& n) { _##y = n; }

// R/W for strings that should be accessible via string_view
#define READ_AND_WRITE_STR(y)                                       \
/*! Read access to _##y */                                      \
  [[nodiscard]] inline const decltype(_##y)& y(void) const& { return _##y; }    \
  [[nodiscard]] inline decltype(_##y) y(void) && { return std::move(_##y); }    \
\
/*! Write access to _##y */                                     \
  inline void y(const decltype(_##y)& n) { _##y = n; } \
  inline void y(const std::string_view n_sv) { _##y = n_sv; }

#endif    // !READ_AND_WRITE

/// Error classes are all similar, so define a macro for them
#if (!defined(ERROR_CLASS))

#define ERROR_CLASS(z) \
\
class z : public x_error \
{ \
protected: \
\
public: \
\
/*! \brief      Construct from error code and reason \
    \param  n   error code \
    \param  s   reason \
*/ \
  inline z(const int n, const std::string& s = (std::string)"") : \
    x_error(n, s) \
  { } \
}

#endif      // !ERROR_CLASS

// a fairly large number of useful type-related functions

template<typename T>
using base_type = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

// see: https://stackoverflow.com/questions/51032671/idiomatic-way-to-write-concept-that-says-that-type-is-a-stdvector

template <class, template<class...> class>
inline constexpr bool is_specialization = false;

template <template<class...> class T, class... Args>
inline constexpr bool is_specialization<T<Args...>, T> = true;

// arrays are different; https://stackoverflow.com/questions/62940931/c20-concepts-and-templates-fixe-size-and-resizable-array
template<class A>
inline constexpr bool is_stdarray = false;

template<class T, std::size_t I>
inline constexpr bool is_stdarray<std::array<T, I>> = true;

// basic types + string
template <class T> concept is_int         = std::is_same_v<T, int>;
template <class T> concept is_uint        = std::is_same_v<T, unsigned int>;
template <class T> concept is_string      = std::is_same_v<T, std::string>;
template <class T> concept is_string_view = std::is_same_v<T, std::string_view>;

// standard containers
template <class T> concept is_array              = is_stdarray<T>;
template <class T> concept is_deque              = is_specialization<T, std::deque>;
template <class T> concept is_list               = is_specialization<T, std::list>;
template <class T> concept is_map                = is_specialization<T, std::map>;
template <class T> concept is_multimap           = is_specialization<T, std::multimap>;
template <class T> concept is_multiset           = is_specialization<T, std::multiset>;
template <class T> concept is_queue              = is_specialization<T, std::queue>;
template <class T> concept is_set                = is_specialization<T, std::set>;
template <class T> concept is_unordered_map      = is_specialization<T, std::unordered_map>;
template <class T> concept is_unordered_multimap = is_specialization<T, std::unordered_multimap>;
template <class T> concept is_unordered_multiset = is_specialization<T, std::unordered_multiset>;
template <class T> concept is_unordered_set      = is_specialization<T, std::unordered_set>;
template <class T> concept is_vector             = is_specialization<T, std::vector>;

// combination concepts
template <class T> concept is_mum   = is_map<T>      or is_unordered_map<T>;
template <class T> concept is_mmumm = is_multimap<T> or is_unordered_multimap<T>;
template <class T> concept is_sus   = is_set<T>      or is_unordered_set<T>;
template <class T> concept is_ssuss = is_multiset<T> or is_unordered_multiset<T>;
template <class T> concept is_ssv   = is_string<T>   or is_string_view<T>;

// combinations of combinations
template <class T> concept ANYSET        = is_sus<T> or is_ssuss<T>;
template <class T> concept is_mum_or_sus = is_mum<T> or is_sus<T>;

// heterogeneous lookup for strings
// https://schneide.blog/2024/10/23/heterogeneous-lookup-in-unordered-c-containers/
struct stringly_hash
{ using is_transparent = void;

  [[nodiscard]] inline size_t operator()(const char* rhs) const
    { return std::hash<std::string_view>{}(rhs); }

  [[nodiscard]] inline size_t operator()(const std::string_view rhs) const
    { return std::hash<std::string_view>{}(rhs); }

  [[nodiscard]] inline size_t operator()(const std::string& rhs) const
    { return std::hash<std::string>{}(rhs); }
};

// make the versions with heterogeneous lookup easily available

// heterogenous lookup for unordered maps with string keys
template <typename ValueType>
using UNORDERED_STRING_MAP = std::unordered_map<std::string, ValueType, stringly_hash, std::equal_to<>>;

// heterogenous lookup for unordered multimaps with string keys
template <typename ValueType>
using UNORDERED_STRING_MULTIMAP = std::unordered_multimap<std::string, ValueType, stringly_hash, std::equal_to<>>;

// heterogenous lookup for unordered sets of strings
using UNORDERED_STRING_SET = std::unordered_set<std::string, stringly_hash, std::equal_to<>>;

// heterogenous lookup for unordered multisets of strings
using UNORDERED_STRING_MULTISET = std::unordered_multiset<std::string, stringly_hash, std::equal_to<>>;

// heterogenous lookup for ordered maps with string keys
template <typename ValueType>
using STRING_MAP = std::map<std::string, ValueType, std::less<>>;

// heterogenous lookup for ordered maps with string keys, with custom ordering
template <typename ValueType, typename Comparison>
using CUSTOM_STRING_MAP = std::map<std::string, ValueType, Comparison>;

// heterogenous lookup for multimaps with string keys
template <typename ValueType>
using STRING_MULTIMAP = std::multimap<std::string, ValueType, std::less<>>;

// heterogenous lookup for ordered sets of strings
using STRING_SET = std::set<std::string, std::less<>>;

// heterogenous lookup for ordered sets of strings with custom ordering
template <typename Comparison>
using CUSTOM_STRING_SET = std::set<std::string, Comparison>;

// heterogenous lookup for ordered multisets of strings
using STRING_MULTISET = std::multiset<std::string, std::less<>>;

// classes for tuples... it seems like there should be a way to do this with TMP,
// but the level-breaking caused by the need to control textual names seems to make
// this impossible without resorting to #defines. So since I don't immediately see
// a way to avoid #defines completely while keeping the usual access syntax
// (i.e., obj.name()), we might as well use a sledgehammer and do everything with #defines.

// we could access with something like obj.at<name>, but that would mean a different access
// style for this kind of object as compared to ordinary classes using the READ and READ_AND_WRITE
// macros

/// tuple class (1) -- complete overkill
#define WRAPPER_1(nm, a0, a1)                          \
                                                       \
class nm : public std::tuple < a0 >                    \
{                                                      \
protected:                                             \
                                                       \
public:                                                \
                                                       \
  explicit nm( const a0 & X )                                           \
    { std::get<0>(*this) = X;                          \
    }                                                  \
                                                       \
  inline a0 a1(void) const                             \
    { return std::get<0>(*this); }                     \
                                                       \
  inline void a1(const a0 & var)                       \
    { std::get<0>(*this) = var; }                      \
}

/// tuple class (2)
#define WRAPPER_2(nm, a0, a1, b0, b1)                  \
                                                       \
class nm : public std::tuple < a0, b0 >                \
{                                                      \
protected:                                             \
                                                       \
public:                                                \
                                                       \
  nm( a0 X, b0 Y)                                      \
    { std::get<0>(*this) = X;                          \
      std::get<1>(*this) = Y;                          \
    }                                                  \
                                                       \
  inline a0 a1(void) const                             \
    { return std::get<0>(*this); }                     \
                                                       \
  inline void a1(const a0 & var)                       \
    { std::get<0>(*this) = var; }                      \
                                                       \
  inline b0 b1(void) const                             \
    { return std::get<1>(*this); }                     \
                                                       \
  inline void b1(const b0 & var)                       \
    { std::get<1>(*this) = var; }                      \
}

/// tuple class (2)
#define WRAPPER_2_NC(nm, a0, a1, b0, b1)               \
                                                       \
class nm : public std::tuple < a0, b0 >                \
{                                                      \
protected:                                             \
                                                       \
public:                                                \
                                                       \
  nm( a0 X, b0 Y)                                      \
    { std::get<0>(*this) = X;                          \
      std::get<1>(*this) = Y;                          \
    }                                                  \
                                                       \
  inline a0 a1(void) const                             \
    { return std::get<0>(*this); }                     \
                                                       \
  inline void a1(a0 var)                               \
    { std::get<0>(*this) = var; }                      \
                                                       \
  inline b0 b1(void) const                             \
    { return std::get<1>(*this); }                     \
                                                       \
  inline void b1(b0 var)                               \
    { std::get<1>(*this) = var; }                      \
}

/// tuple class (3)
#define WRAPPER_3(nm, a0, a1, b0, b1, c0, c1)          \
                                                       \
class nm : public std::tuple < a0, b0, c0 >            \
{                                                      \
protected:                                             \
                                                       \
public:                                                \
                                                       \
  nm( a0 X, b0 Y, c0 Z)                                \
    { std::get<0>(*this) = X;                          \
      std::get<1>(*this) = Y;                          \
      std::get<2>(*this) = Z;                          \
    }                                                  \
                                                       \
  nm( void ) { }                                       \
                                                       \
  inline a0 a1(void) const                             \
    { return std::get<0>(*this); }                     \
                                                       \
  inline void a1(const a0 & var)                       \
    { std::get<0>(*this) = var; }                      \
                                                       \
  inline b0 b1(void) const                             \
    { return std::get<1>(*this); }                     \
                                                       \
  inline void b1(const b0 & var)                       \
    { std::get<1>(*this) = var; }                      \
                                                       \
  inline c0 c1(void) const                             \
    { return std::get<2>(*this); }                     \
                                                       \
  inline void c1(const c0 & var)                       \
    { std::get<2>(*this) = var; }                      \
};                                                     \
                                                       \
inline std::ostream& operator<<(std::ostream& ost, const nm& type)  \
{ ost << #nm << ": " << std::endl                                   \
      << #a1 << ": " << type.a1() << std::endl                      \
      << #b1 << ": " << type.b1() << std::endl                      \
      << #c1 << ": " << type.c1();                                  \
                                                                    \
  return ost;                                                       \
}

/// tuple class (3)
#define WRAPPER_3_NC(nm, a0, a1, b0, b1, c0, c1)        \
                                                        \
class nm : public std::tuple < a0, b0, c0 >             \
{                                                       \
protected:                                              \
                                                        \
public:                                                 \
                                                        \
  nm( a0 X, b0 Y, c0 Z)                                 \
    { std::get<0>(*this) = X;                           \
      std::get<1>(*this) = Y;                           \
      std::get<2>(*this) = Z;                           \
    }                                                   \
                                                        \
  inline a0 a1(void) const                              \
    { return std::get<0>(*this); }                      \
                                                        \
  inline void a1(a0 var)                                \
    { std::get<0>(*this) = var; }                       \
                                                        \
  inline b0 b1(void) const                              \
    { return std::get<1>(*this); }                      \
                                                        \
  inline void b1(b0 var)                                \
    { std::get<1>(*this) = var; }                       \
                                                        \
  inline c0 c1(void) const                              \
    { return std::get<2>(*this); }                      \
                                                        \
  inline void c1(c0 var)                                \
    { std::get<2>(*this) = var; }                       \
}

/// tuple class (3)
#define WRAPPER_3_SERIALIZE(nm, a0, a1, b0, b1, c0, c1)                     \
                                                                            \
class nm : public std::tuple < a0, b0, c0 >                                 \
{                                                                           \
protected:                                                                  \
                                                                            \
public:                                                                     \
                                                                            \
  nm( a0 X, b0 Y, c0 Z)                                                     \
    { std::get<0>(*this) = X;                                               \
      std::get<1>(*this) = Y;                                               \
      std::get<2>(*this) = Z;                                               \
    }                                                                       \
                                                                            \
  nm(void) = default;                                                       \
                                                                            \
  inline a0 a1(void) const                                                  \
    { return std::get<0>(*this); }                                          \
                                                                            \
  inline void a1(const a0 & var)                                            \
    { std::get<0>(*this) = var; }                                           \
                                                                            \
  inline b0 b1(void) const                                                  \
    { return std::get<1>(*this); }                                          \
                                                                            \
  inline void b1(const b0 & var)                                            \
    { std::get<1>(*this) = var; }                                           \
                                                                            \
  inline c0 c1(void) const                                                  \
    { return std::get<2>(*this); }                                          \
                                                                            \
  inline void c1(const c0 & var)                                            \
    { std::get<2>(*this) = var; }                                           \
                                                                            \
  template<typename Archive>                                                \
  void serialize(Archive& ar, const unsigned version)                       \
    { ar & std::get<0>(*this)                                               \
         & std::get<1>(*this)                                               \
         & std::get<2>(*this);                                              \
    }                                                                       \
}

/// tuple class (4)
#define WRAPPER_4_NC(nm, a0, a1, b0, b1, c0, c1, d0, d1)                    \
                                                                            \
class nm : public std::tuple < a0, b0, c0, d0 >                             \
{                                                                           \
protected:                                                                  \
                                                                            \
public:                                                                     \
                                                                            \
  nm( a0 X, b0 Y, c0 Z, d0 A)                                               \
    { std::get<0>(*this) = X;                                               \
      std::get<1>(*this) = Y;                                               \
      std::get<2>(*this) = Z;                                               \
      std::get<3>(*this) = A;                                               \
}                                                                           \
                                                                            \
  inline a0 a1(void) const                                                  \
    { return std::get<0>(*this); }                                          \
                                                                            \
  inline void a1(a0 var)                                                    \
    { std::get<0>(*this) = var; }                                           \
                                                                            \
  inline b0 b1(void) const                                                  \
    { return std::get<1>(*this); }                                          \
                                                                            \
  inline void b1(b0 var)                                                    \
    { std::get<1>(*this) = var; }                                           \
                                                                            \
  inline c0 c1(void) const                                                  \
    { return std::get<2>(*this); }                                          \
                                                                            \
  inline void c1(c0 var)                                                    \
    { std::get<2>(*this) = var; }                                           \
                                                                            \
  inline d0 d1(void) const                                                  \
    { return std::get<3>(*this); }                                          \
                                                                            \
  inline void d1(d0 var)                                                    \
    { std::get<3>(*this) = var; }                                           \
}

/// tuple class (4)
#define WRAPPER_4_SERIALIZE(nm, a0, a1, b0, b1, c0, c1, d0, d1)             \
                                                                            \
class nm : public std::tuple < a0, b0, c0, d0  >                            \
{                                                                           \
protected:                                                                  \
                                                                            \
public:                                                                     \
                                                                            \
  nm( a0 X, b0 Y, c0 Z, d0 A)                                               \
    { std::get<0>(*this) = X;                                               \
      std::get<1>(*this) = Y;                                               \
      std::get<2>(*this) = Z;                                               \
      std::get<3>(*this) = A;                                               \
    }                                                                       \
                                                                            \
  nm(void) = default;                                                       \
                                                                            \
  inline a0 a1(void) const                                                  \
    { return std::get<0>(*this); }                                          \
                                                                            \
  inline void a1(const a0 & var)                                            \
    { std::get<0>(*this) = var; }                                           \
                                                                            \
  inline b0 b1(void) const                                                  \
    { return std::get<1>(*this); }                                          \
                                                                            \
  inline void b1(const b0 & var)                                            \
    { std::get<1>(*this) = var; }                                           \
                                                                            \
  inline c0 c1(void) const                                                  \
    { return std::get<2>(*this); }                                          \
                                                                            \
  inline void c1(const c0 & var)                                            \
    { std::get<2>(*this) = var; }                                           \
                                                                            \
  inline d0 d1(void) const                                                  \
    { return std::get<3>(*this); }                                          \
                                                                            \
  inline void d1(const d0 & var)                                            \
    { std::get<3>(*this) = var; }                                           \
                                                                            \
  template<typename Archive>                                                \
  void serialize(Archive& ar, const unsigned version)                       \
    { ar & std::get<0>(*this)                                               \
         & std::get<1>(*this)                                               \
         & std::get<2>(*this)                                               \
         & std::get<3>(*this);                                              \
    }                                                                       \
};                                                                          \
                                                                            \
inline std::ostream& operator<<(std::ostream& ost, const nm& type)          \
{ ost << #nm << ": " << std::endl                                           \
      << #a1 << ": " << type.a1() << std::endl                              \
      << #b1 << ": " << type.b1() << std::endl                              \
      << #c1 << ": " << type.c1() << std::endl                              \
      << #d1 << ": " << type.d1();                                          \
                                                                            \
  return ost;                                                               \
}

/// tuple class (5)
#define WRAPPER_5_NC(nm, a0, a1, b0, b1, c0, c1, d0, d1, e0, e1)            \
                                                                            \
class nm : public std::tuple < a0, b0, c0, d0, e0 >                         \
{                                                                           \
protected:                                                                  \
                                                                            \
public:                                                                     \
                                                                            \
  nm( a0 X, b0 Y, c0 Z, d0 A, e0 B)                                         \
    { std::get<0>(*this) = X;                                               \
      std::get<1>(*this) = Y;                                               \
      std::get<2>(*this) = Z;                                               \
      std::get<3>(*this) = A;                                               \
      std::get<4>(*this) = B;                                               \
    }                                                                       \
                                                                            \
  nm(void) = default;                                                       \
                                                                            \
  inline a0 a1(void) const                                                  \
    { return std::get<0>(*this); }                                          \
                                                                            \
  inline void a1(a0 var)                                                    \
    { std::get<0>(*this) = var; }                                           \
                                                                            \
  inline b0 b1(void) const                                                  \
    { return std::get<1>(*this); }                                          \
                                                                            \
  inline void b1(b0 var)                                                    \
    { std::get<1>(*this) = var; }                                           \
                                                                            \
  inline c0 c1(void) const                                                  \
    { return std::get<2>(*this); }                                          \
                                                                            \
  inline void c1(c0 var)                                                    \
    { std::get<2>(*this) = var; }                                           \
                                                                            \
  inline d0 d1(void) const                                                  \
    { return std::get<3>(*this); }                                          \
                                                                            \
  inline void d1(d0 var)                                                    \
    { std::get<3>(*this) = var; }                                           \
                                                                            \
  inline e0 e1(void) const                                                  \
    { return std::get<4>(*this); }                                          \
                                                                            \
  inline void e1(e0 var)                                                    \
    { std::get<4>(*this) = var; }                                           \
}

/// tuple class (6)
#define WRAPPER_6_NC(nm, a0, a1, b0, b1, c0, c1, d0, d1, e0, e1, f0, f1)    \
                                                                            \
class nm : public std::tuple < a0, b0, c0, d0, e0, f0 >                     \
{                                                                           \
protected:                                                                  \
                                                                            \
public:                                                                     \
                                                                            \
  nm( a0 X, b0 Y, c0 Z, d0 A, e0 B, f0 C)                                   \
    { std::get<0>(*this) = X;                                               \
      std::get<1>(*this) = Y;                                               \
      std::get<2>(*this) = Z;                                               \
      std::get<3>(*this) = A;                                               \
      std::get<4>(*this) = B;                                               \
      std::get<5>(*this) = C;                                               \
}                                                                           \
                                                                            \
  inline a0 a1(void) const                                                  \
    { return std::get<0>(*this); }                                          \
                                                                            \
  inline void a1(a0 var)                                                    \
    { std::get<0>(*this) = var; }                                           \
                                                                            \
  inline b0 b1(void) const                                                  \
    { return std::get<1>(*this); }                                          \
                                                                            \
  inline void b1(b0 var)                                                    \
    { std::get<1>(*this) = var; }                                           \
                                                                            \
  inline c0 c1(void) const                                                  \
    { return std::get<2>(*this); }                                          \
                                                                            \
  inline void c1(c0 var)                                                    \
    { std::get<2>(*this) = var; }                                           \
                                                                            \
  inline d0 d1(void) const                                                  \
    { return std::get<3>(*this); }                                          \
                                                                            \
  inline void d1(d0 var)                                                    \
    { std::get<3>(*this) = var; }                                           \
                                                                            \
  inline e0 e1(void) const                                                  \
    { return std::get<4>(*this); }                                          \
                                                                            \
  inline void e1(e0 var)                                                    \
    { std::get<4>(*this) = var; }                                           \
                                                                            \
  inline f0 f1(void) const                                                  \
    { return std::get<5>(*this); }                                          \
                                                                            \
  inline void f1(f0 var)                                                    \
    { std::get<5>(*this) = var; }                                           \
}

/// tuple class (7)
#define WRAPPER_7_NC(nm, a0, a1, b0, b1, c0, c1, d0, d1, e0, e1, f0, f1, g0, g1)    \
                                                                                    \
class nm : public std::tuple < a0, b0, c0, d0, e0, f0, g0 >                         \
{                                                                                   \
protected:                                                                          \
                                                                                    \
public:                                                                             \
                                                                                    \
  nm( a0 X, b0 Y, c0 Z, d0 A, e0 B, f0 C, g0 D)                                     \
    { std::get<0>(*this) = X;                                                       \
      std::get<1>(*this) = Y;                                                       \
      std::get<2>(*this) = Z;                                                       \
      std::get<3>(*this) = A;                                                       \
      std::get<4>(*this) = B;                                                       \
      std::get<5>(*this) = C;                                                       \
      std::get<6>(*this) = D;                                                       \
}                                                                                   \
                                                                                    \
  inline a0 a1(void) const                                                          \
    { return std::get<0>(*this); }                                                  \
                                                                                    \
  inline void a1(a0 var)                                                            \
    { std::get<0>(*this) = var; }                                                   \
                                                                                    \
  inline b0 b1(void) const                                                          \
    { return std::get<1>(*this); }                                                  \
                                                                                    \
  inline void b1(b0 var)                                                            \
    { std::get<1>(*this) = var; }                                                   \
                                                                                    \
  inline c0 c1(void) const                                                          \
    { return std::get<2>(*this); }                                                  \
                                                                                    \
  inline void c1(c0 var)                                                            \
    { std::get<2>(*this) = var; }                                                   \
                                                                                    \
  inline d0 d1(void) const                                                          \
    { return std::get<3>(*this); }                                                  \
                                                                                    \
  inline void d1(d0 var)                                                            \
    { std::get<3>(*this) = var; }                                                   \
                                                                                    \
  inline e0 e1(void) const                                                          \
    { return std::get<4>(*this); }                                                  \
                                                                                    \
  inline void e1(e0 var)                                                            \
    { std::get<4>(*this) = var; }                                                   \
                                                                                    \
  inline f0 f1(void) const                                                          \
    { return std::get<5>(*this); }                                                  \
                                                                                    \
  inline void f1(f0 var)                                                            \
    { std::get<5>(*this) = var; }                                                   \
                                                                                    \
  inline g0 g1(void) const                                                          \
    { return std::get<6>(*this); }                                                  \
                                                                                    \
  inline void g1(g0 var)                                                            \
    { std::get<6>(*this) = var; }                                                   \
}

/// tuple class (8)
#define WRAPPER_8_NC(nm, a0, a1, b0, b1, c0, c1, d0, d1, e0, e1, f0, f1, g0, g1, h0, h1)    \
                                                                                            \
class nm : public std::tuple < a0, b0, c0, d0, e0, f0, g0, h0 >                             \
{                                                                                           \
protected:                                                                                  \
                                                                                            \
public:                                                                                     \
                                                                                            \
  nm( a0 X, b0 Y, c0 Z, d0 A, e0 B, f0 C, g0 D, h0 E)                                       \
    { std::get<0>(*this) = X;                                                               \
      std::get<1>(*this) = Y;                                                               \
      std::get<2>(*this) = Z;                                                               \
      std::get<3>(*this) = A;                                                               \
      std::get<4>(*this) = B;                                                               \
      std::get<5>(*this) = C;                                                               \
      std::get<6>(*this) = D;                                                               \
      std::get<7>(*this) = E;                                                               \
}                                                                                           \
                                                                                            \
  inline a0 a1(void) const                                                                  \
    { return std::get<0>(*this); }                                                          \
                                                                                            \
  inline void a1(a0 var)                                                                    \
    { std::get<0>(*this) = var; }                                                           \
                                                                                            \
  inline b0 b1(void) const                                                                  \
    { return std::get<1>(*this); }                                                          \
                                                                                            \
  inline void b1(b0 var)                                                                    \
    { std::get<1>(*this) = var; }                                                           \
                                                                                            \
  inline c0 c1(void) const                                                                  \
    { return std::get<2>(*this); }                                                          \
                                                                                            \
  inline void c1(c0 var)                                                                    \
    { std::get<2>(*this) = var; }                                                           \
                                                                                            \
  inline d0 d1(void) const                                                                  \
    { return std::get<3>(*this); }                                                          \
                                                                                            \
  inline void d1(d0 var)                                                                    \
    { std::get<3>(*this) = var; }                                                           \
                                                                                            \
  inline e0 e1(void) const                                                                  \
    { return std::get<4>(*this); }                                                          \
                                                                                            \
  inline void e1(e0 var)                                                                    \
    { std::get<4>(*this) = var; }                                                           \
                                                                                            \
  inline f0 f1(void) const                                                                  \
    { return std::get<5>(*this); }                                                          \
                                                                                            \
  inline void f1(f0 var)                                                                    \
    { std::get<5>(*this) = var; }                                                           \
                                                                                            \
  inline g0 g1(void) const                                                                  \
    { return std::get<6>(*this); }                                                          \
                                                                                            \
  inline void g1(g0 var)                                                                    \
    { std::get<6>(*this) = var; }                                                           \
                                                                                            \
  inline h0 h1(void) const                                                                  \
    { return std::get<7>(*this); }                                                          \
                                                                                            \
  inline void h1(h0 var)                                                                    \
    { std::get<7>(*this) = var; }                                                           \
}

/*! \brief      Is an object a member of a set, unordered_set, map or unordered map?
    \param  s   set, unordered_set, map or unordered map to be tested
    \param  v   object to be tested for membership
    \return     Whether <i>v</i> is a member of <i>s</i>
*/
template <class T, class U>
inline bool operator>(const T& s, const U& v)
  requires (is_sus<T> or is_mum<T>) and (std::is_same_v<typename T::key_type, U>)
{ return s.find(v) != s.cend(); }

/*! \brief      Is an object a key of a map or unordered map; if so return the value, otherwise return a provided default
    \param  m   map or unordered map to be searched
    \param  k   target key
    \param  d   default
    \return     if <i>k</i> is a member of <i>m</i>, the corresponding value, otherwise the default
*/
template <class C, class K>
  requires (is_mum<C>) and (std::is_constructible_v<typename C::key_type, K>)
typename C::mapped_type MUM_VALUE(const C& m, const K& k, const typename C::mapped_type& d = typename C::mapped_type())
{ const auto cit { m.find(k) };

  return ( (cit == m.cend()) ? d : cit->second );
}

/*! \brief                      Invert a mapping from map<T, set<T> > to map<T, set<T> >, where final keys are the elements of the original set
    \param  original_mapping    original mapping
    \return                     inverted mapping
*/
template <typename M>  // M = map<T, set<T> >
auto INVERT_MAPPING(const M& original_mapping) -> std::map<typename M::key_type, typename M::key_type>
{ using RT = std::invoke_result_t<decltype(INVERT_MAPPING<M>), const M&>;
  
  RT rv;

  for (auto cit { original_mapping.cbegin() }; cit != original_mapping.cend(); ++cit)
  { for (const auto& p : cit->second)
      rv += { p, cit->first };
  }

  return rv;
}

// convenient syntactic sugar for some STL functions

/*! \brief          Write a <i>map<key, value></i> object to an output stream
    \param  ost     output stream
    \param  mp      object to write
    \return         the output stream
*/
//std::ostream& operator<<(std::ostream& ost, const std::map<auto, auto>& mp)
//{ for (auto cit { mp.cbegin() }; cit != mp.cend(); ++cit)
//    ost << "map[" << cit->first << "]: " << cit->second << std::endl;
//
//  return ost;
//}

template <class T>
  requires is_map<T>
std::ostream& operator<<(std::ostream& ost, const T& mp)
{ for (auto cit = mp.cbegin(); cit != mp.cend(); ++cit)
    ost << "map[" << cit->first << "]: " << cit->second << std::endl;

  return ost;
}

/*! \brief          Apply a function to all in a container
    \param  first   container
    \param  fn      function
    \return         <i>fn</i>
*/

// see https://en.cppreference.com/w/cpp/algorithm/ranges/for_each
template< std::input_iterator I, std::sentinel_for<I> S, class Proj = std::identity, std::indirectly_unary_invocable<std::projected<I, Proj>> Fun >
constexpr std::ranges::for_each_result<I, Fun> FOR_ALL( I first, S last, Fun f, Proj proj = {} ) { return std::ranges::for_each(first, last, f, proj); }; 

template< std::ranges::input_range R, class Proj = std::identity, std::indirectly_unary_invocable<std::projected<std::ranges::iterator_t<R>, Proj>> Fun >
constexpr std::ranges::for_each_result<std::ranges::borrowed_iterator_t<R>, Fun> FOR_ALL( R&& r, Fun f, Proj proj = {} ) { return std::ranges::for_each(std::forward<R>(r), f, proj); }; 

// https://en.cppreference.com/w/cpp/algorithm/ranges/all_any_none_of
template< std::input_iterator I, std::sentinel_for<I> S, class Proj = std::identity, std::indirect_unary_predicate<std::projected<I, Proj>> Pred >
constexpr bool ANY_OF( I first, S last, Pred pred, Proj proj = {} ) { return std::ranges::any_of(first, last, pred, proj); };

template< std::ranges::input_range R, class Proj = std::identity, std::indirect_unary_predicate<std::projected<std::ranges::iterator_t<R>,Proj>> Pred >
constexpr bool ANY_OF( R&& r, Pred pred, Proj proj = {} ) { return std::ranges::any_of(std::forward<R>(r), pred, proj); };

template< std::input_iterator I, std::sentinel_for<I> S, class Proj = std::identity, std::indirect_unary_predicate<std::projected<I, Proj>> Pred >
constexpr bool ALL_OF( I first, S last, Pred pred, Proj proj = {} ) { return std::ranges::all_of(first, last, pred, proj); };

template< std::ranges::input_range R, class Proj = std::identity, std::indirect_unary_predicate<std::projected<std::ranges::iterator_t<R>,Proj>> Pred >
constexpr bool ALL_OF( R&& r, Pred pred, Proj proj = {} ) { return std::ranges::all_of(std::forward<R>(r), pred, proj); };

template< std::input_iterator I, std::sentinel_for<I> S, class Proj = std::identity, std::indirect_unary_predicate<std::projected<I, Proj>> Pred >
constexpr bool NONE_OF( I first, S last, Pred pred, Proj proj = {} )  { return std::ranges::none_of(first, last, pred, proj); };

template< std::ranges::input_range R, class Proj = std::identity, std::indirect_unary_predicate<std::projected<std::ranges::iterator_t<R>,Proj>> Pred >
constexpr bool NONE_OF( R&& r, Pred pred, Proj proj = {} ) { return std::ranges::none_of(std::forward<R>(r), pred, proj); };

// transform in place
// https://en.cppreference.com/w/cpp/algorithm/ranges/transform
template< std::ranges::input_range R, std::copy_constructible F, class Proj = std::identity >
constexpr auto TIP( R&& r, F op, Proj proj = {} ) -> std::ranges::in_out_result<std::ranges::borrowed_iterator_t<R>, decltype(begin(r))> 
  { return std::ranges::transform(r, r.begin(), op, proj); };

/*! \brief          Copy all in a container to another container
    \param  first   initial container
    \param  oi      iterator on final container
    \return         <i>oi</i>
*/
template<class Input, class OutputIterator>
inline OutputIterator COPY_ALL(Input& first, OutputIterator oi)
  { return (std::copy(first.begin(), first.end(), oi)); }

/*! \brief          Remove values in a container that match a predicate, and resize the container
    \param  first   container
    \param  pred    predicate to apply

    Does not work for maps
*/
//template <class Input, class UnaryPredicate>
//inline void REMOVE_IF_AND_RESIZE(Input& first, UnaryPredicate pred)
//  { first.erase(std::remove_if(first.begin(), first.end(), pred), first.end()); }

#if 0
/*! \brief          Remove map values that match a predicate, and resize the map
    \param  items   map
    \param  pred    predicate to apply
    
    I think that this is the same as erase_if() now
*/
template<typename M, typename PredicateT>
void REMOVE_IF_AND_RESIZE(M& items, const PredicateT& pred )
 requires (is_mum<M>)
{ for (auto it { items.begin() }; it != items.end(); )
  { if (pred(*it))
      it = items.erase(it);
    else
      ++it;
  }
}
#endif

/*! \brief          Keep only values in a container that match a predicate, and resize the container
    \param  first   container
    \param  pred    predicate to apply

    Does not work for maps
*/
template <class Input, class UnaryPredicate>
inline void KEEP_IF_AND_RESIZE(Input& first, UnaryPredicate pred)
  { first.erase(std::remove_if(first.begin(), first.end(), std::experimental::not_fn(pred)), first.end()); }

/*! \brief      Reverse the contents of a container
    \param  v   container
*/
template <class Input>
inline void REVERSE(Input& v)
  { std::reverse(v.begin(), v.end()); }

/*! \brief          Append from one container to another if a predicate is met
    \param  s       source container
    \param  d       destination container
    \param  pred    predicate
*/
template <class Input, class Output, typename PredicateT>
inline void APPEND_IF(const Input& s, Output& d, const PredicateT& pred)
  { std::copy_if(s.cbegin(), s.cend(), std::back_inserter(d), pred); }

/*! \brief          Create a container and append from another if a predicate is met
    \param  s       source container
    \param  pred    predicate
    \return         the new container

    Called as, for example, CREATE_AND_FILL<vector<string>>(in_vec, [](const string& s) { return (s == "BURBLE"s); } );
*/
template <typename Output, typename Input, typename PredicateT>
auto CREATE_AND_FILL(Input&& s, const PredicateT& pred) -> const Output
{ Output rv;

  std::ranges::copy_if(s, std::back_inserter(rv), pred);

  return rv;
}

/*! \brief          Find first value in a container that matches a predicate
    \param  v       container
    \param  pred    (boolean) predicate to apply
    \return         first value in <i>v</i> for which <i>pred</i> is true
*/
template <typename Input, typename UnaryPredicate>
inline auto FIND_IF(Input& v, UnaryPredicate pred) -> typename Input::iterator
  { return std::ranges::find_if(v, pred); }

/*! \brief          Find first value in a container that matches a predicate
    \param  v       container (const)
    \param  pred    (boolean) predicate to apply
    \return         first value in <i>v</i> for which <i>pred</i> is true
*/
template <typename Input, typename UnaryPredicate>
inline auto FIND_IF(const Input& v, UnaryPredicate pred) -> typename Input::const_iterator
  { return std::ranges::find_if(v, pred); }

/*! \brief              Bound a value within limits
    \param  val         value to bound
    \param  low_val     lower bound
    \param  high_val    upper bound
    \return             max(min(<i>val</i>, <i>max_val</i>), <i>min_val</i>)
*/
template <typename T>
inline T LIMIT(const T val, const T low_val, const T high_val)
  { return (val < low_val ? low_val : (val > high_val ? high_val : val)); }

/*! \brief              Bound a value within limits
    \param  val         value to bound
    \param  low_val     lower bound
    \param  high_val    upper bound
    \return             max(min(<i>val</i>, <i>max_val</i>), <i>min_val</i>)
*/
template <typename T, typename U, typename V>
inline T LIMIT(const T val, const U low_val, const V high_val)
  { return (val < static_cast<T>(low_val) ? static_cast<T>(low_val) : (val > static_cast<T>(high_val) ? static_cast<T>(high_val) : val)); }

/*! \brief              Return first true element in vector of pair<bool, value> elements
    \param  vec         the vector
    \param  def         value to return of no elements are true
    \return             either the first true element or, if none, <i>def</i>
*/
template <typename T>
T SELECT_FIRST_TRUE(std::vector<std::pair<bool, T>>& vec, const T& def)
{ for (const auto& element : vec)
  { if (element.first)
      return (element.second);
  }

  return def;
}

/*! \brief              Execute the same (void) function on multiple threads
    \param  n_threads   the number of threads on which to run the function <i>fn</i>
    \param  fn          the function to execute
    \param  args        arguments to <i>fn</i>

    This can be useful, for example, when reading large files, each line of which
    needs non-negligible processing and the order of processing is unimportant
*/
template <typename Function, typename... Args>
void EXECUTE_FUNCTION_MT(const unsigned int n_threads, Function&& fn, Args&&... args)
{ std::vector<std::future<void>> vec_futures;           // place to store the (void) futures

  for (unsigned int n { 0 }; n < n_threads; ++n)
    vec_futures.emplace_back(async(std::launch::async, fn, args...));

  for (auto& this_future : vec_futures)
    this_future.get();                                  // .get() blocks until the future is available
}

/*! \brief              Obtain a copy of the element with the minimum value in a container
    \param  container   the container
    \return             copy of element of <i>container</i> with the minimum value
*/
template <typename C>
inline auto MIN_ELEMENT(const C& container)
  { return ( *(std::min_element(begin(container), end(container))) ); }

/*! \brief              Obtain a copy of the element with the maximum value in a container
    \param  container   the container
    \return             copy of element of <i>container</i> with the maximum value
*/
template <typename C>
inline auto MAX_ELEMENT(const C& container)
  { return ( *(std::max_element(begin(container), end(container))) ); }

/*! \brief      Sort the contents of a container
    \param  v   container
    \param  f   sort function
*/
template <typename C, typename F = std::less<>>
inline void SORT(C& v, F f = F())
  { std::sort(v.begin(), v.end(), f); }

/*! \brief              A somewhat useful class to allow inverse mapping, as long as both keys and elements are unique
*/
template <typename T, typename U>
class invertible_unordered_map : public std::unordered_map<T, U>
{
protected:
  
public:

  const std::unordered_map<U, T> invert(void)
  { std::unordered_map<U, T> rv;
  
    FOR_ALL(*this, [&rv](const auto& forward_element) { rv.insert( { forward_element.second, forward_element.first } ); });
    
    return rv;
  }
};

// ------------------------------------------------------ container operations --------------------------------------

/*! \brief          Split a vector into chunks
    \param  invec   vector
    \param  n       number of chunks
    \return         vector of vectors
    
    Not robust against really silly parameter values
*/
template <typename T>
inline std::vector<T> SPLIT_VECTOR(const T& invec, const int n)
{ std::vector<T> rv(n);

  const size_t sz   { invec.size() };
  
  if (sz < static_cast<size_t>(n))       // just return the original vector wrapped in a vector if n is too small
  { rv += invec;
    return rv;
  }
  
  const size_t incr { std::max(sz / n, static_cast<size_t>(1)) };
  
  auto this_start { invec.cbegin() };
  auto this_end   { invec.cbegin() + incr }; 

  for (int index { 0 }; index < n; ++index)
  { rv[index] = T(this_start, this_end);
  
    this_start += incr;
    
    if (index == (n - 2))           // make sure that we go to the end on the last one
      this_end = invec.cend();
    else
      this_end += incr;
  } 

  return rv;
}

/*! \brief        Print the first portion of a vector
    \param  vec   vector
    \param  n     number of elements to print
*/
template <typename T>
  requires (is_vector<T>)
void PRINT_VECTOR(const T& vec, const int n)
{ const int limit { std::min(static_cast<int>(vec.size()), n) };
  
  for (int n { 0 }; n < limit; ++n)
    std::cout << "[" << n << "] : " << vec[n] << std::endl;
}

/*! \brief        Append one vector to another
    \param  dest  destination vector
    \param  src   source vector
    \return       <i>dest</i> with <i>src</i> appended
*/
template <typename V>
  requires (is_vector<V>)
void operator+=(V& dest, V&& src)
{ dest.reserve(dest.size() + src.size());
  dest.insert(dest.end(), src.begin(), src.end());
}

/*! \brief        Append one vector to another
    \param  dest  destination vector
    \param  src   source vector
    \return       <i>dest</i> with <i>src</i> appended
*/
template <typename V>
  requires (is_vector<V>)
inline void operator+=(V& dest, const V& src)
{ dest.reserve(dest.size() + src.size());
  dest.insert(dest.end(), src.begin(), src.end());
}

/*! \brief        Append one vector to another
    \param  dest  destination vector
    \param  src   source vector
    \return       <i>dest</i> with <i>src</i> appended
*/
template <typename V1, typename V2>
  requires is_vector<V1> and is_vector<V2> and std::convertible_to<typename V2::value_type, typename V1::value_type>
void operator+=(V1& dest, const V2& src)
{ dest.reserve(dest.size() + src.size());

  FOR_ALL(src, [&dest] (const auto& v2) { dest.push_back(static_cast<typename V1::value_type>(v2)); } );
}

/*! \brief        Append a vector to a single element
    \param  iv    initial value
    \param  src   source vector
    \return       vector containing <i>iv</i> with <i>src</i> appended
*/
template <typename T, typename V>
  requires is_vector<V> and std::convertible_to<typename V::value_type, T>
std::vector<T> operator+(T v, const V& src)
{ std::vector<T> rv;

  rv.reserve(1 + src.size());
  
  rv.push_back(v);
  FOR_ALL(src, [&rv] (const auto& v2) { rv.push_back(static_cast<typename V::value_type>(v2)); } );
  
  return rv;
}

/*! \brief        Concatenate two vectors
    \param  v1  first vector
    \param  v2  second vector
    \return     concatenation of <i>v1</i> and <i>v2</i>
*/
template <typename V>
  requires (is_vector<V>)
V operator+(const V& v1, V&& v2)
{ V rv(v1.size() + v2.size());

  rv = v1;
  rv += std::forward<V>(v2);
  
  return rv;
}

/*! \brief        Concatenate two vectors
    \param  v1  first vector
    \param  v2  second vector
    \return     concatenation of <i>v1</i> and <i>v2</i>
*/
template <typename V>
  requires (is_vector<V>)
V operator+(const V& v1, const V& v2)
{ V rv(v1.size() + v2.size());

  rv = v1;
  rv += v2;

  return rv;
}

/*! \brief              Add an element to a vector
    \param  dest        destination vector
    \param  element     element to append
    \return             <i>dest</i> with <i>element</i> appended
*/
template <typename V, typename E>
auto operator+(const V& v1, E&& element) -> V
  requires (is_vector<V> or std::derived_from<V, std::vector<typename V::value_type>>) and std::convertible_to<base_type<E>, typename V::value_type>
{ V rv(v1.size() + 1);

  rv = v1;
  rv.push_back(std::forward<E>(element));
  
  return rv; 
}

/*! \brief              Add an element to a vector
    \param  v1          destination vector
    \param  element     element to append
*/
template <typename V, typename E>
inline void operator+=(V& v1, E&& element)
  requires is_vector<V> and std::convertible_to<base_type<E>, typename V::value_type>
{ v1.push_back(std::forward<E>(element)); }

/*! \brief              Add an element to a vector
    \param  v1          destination vector
    \param  element     element to append
*/
template <typename V, typename E>
inline void operator+=(V& v1, const E& element)
  requires (is_vector<V> or std::derived_from<V, std::vector<typename V::value_type>>) and (std::convertible_to<base_type<E>, typename V::value_type>)
//  requires is_vector<V> and std::convertible_to<base_type<E>, typename V::value_type>
{ v1.push_back(element); }

/*! \brief              Append string_view element to a deque, list or vector of strings
    \param  c1          destination deque, list or vector
    \param  element     element to append
*/
template <typename C>
inline void operator+=(C& c1, std::string_view element)
  requires (is_deque<C> or is_list<C> or is_vector<C>) and (is_string<typename C::value_type>)
{ c1.emplace_back(std::string { element }); }

/*! \brief              Add an element to a set
    \param  s           initial set
    \param  element     element to append
    \return             <i>s</i> with <i>element</i> inserted
*/
template <typename S, typename E>
auto operator+(const S& s, const E& element) -> S
  requires is_set<S> and (std::is_same_v<typename S::value_type, base_type<E>>)
{ S rv { s };

  rv.insert(element);
  
  return rv; 
}

/*! \brief        insert one MUM/SUS into another
    \param  dest  destination MUM/SUS
    \param  src   source MUM/SUS
    \return       <i>dest</i> with <i>src</i> inserted
*/
template <typename MUMD, typename MUMS>
  requires ( ( (is_mum<MUMD>) and (is_mum<MUMS>) and
             (std::is_same_v<typename MUMD::key_type, typename MUMS::key_type>) and (std::is_same_v<typename MUMD::mapped_type, typename MUMS::mapped_type>) ) or
             ( (is_sus<MUMD>) and (is_sus<MUMS>) and (std::is_same_v<typename MUMD::value_type, typename MUMS::value_type>) ) )
inline void operator+=(MUMD& dest, const MUMS& src)
  { dest.insert(src.cbegin(), src.cend()); }

/*! \brief              Add an element to a MUM
    \param  m           destination MUM
    \param  element     element to insert
*/
template <typename C, typename K, typename V>
inline void operator+=(C& mum, std::pair<K, V>&& element)
  requires ( (is_mum<C> and (std::is_same_v<typename C::key_type, K>) and (std::is_same_v<typename C::mapped_type, V>)) or
             (is_mmumm<C> and (std::is_same_v<typename C::key_type, K>) and (std::is_same_v<typename C::mapped_type, V>))
           )
  { mum.insert(std::move(element)); }

/*! \brief              Add an element to a MUM
    \param  m           destination MUM
    \param  element     element to insert
*/
/* The perceived "wisdom" is that requires clauses are clearer than SFINAE,. In this case, that is clearly untrue.
    Compare this to the above function, where I purposefully used requires clauses to achieve the same result
*/
//template <typename C, typename B>
//inline void operator+=(C& mum, const std::pair<typename C::key_type, typename C::mapped_type>& il)
//  requires (is_mum<C> or is_mmumm<C>) or (std::is_base_of(B, C) and is_mum<B>)
//  { mum.insert(il); }


template <typename C>
inline void operator+=(C& mum, const std::pair<typename C::key_type, typename C::mapped_type>& il)
  requires (is_mum<C> or is_mmumm<C>)
  { mum.insert(il); }
  
/*! \brief              Add an element to a set, unordered set, multiset or unordered multiset
    \param  m           destination set or unordered set
    \param  element     element to insert
*/
template <typename C, typename T>
inline void operator+=(C& sus, T&& element)
  requires (is_sus<C> or is_multiset<C> or is_unordered_multiset<C>) and (std::is_same_v<typename C::value_type, base_type<T>>)
  { sus.insert(std::forward<T>(element)); }

template <typename C>
inline void operator+=(C& sus, const typename C::value_type& element)
  requires (is_sus<C> or is_multiset<C> or is_unordered_multiset<C>)
  { sus.insert(element); }

/// add all the elements of a [vector] container to a set or unordered set
template <typename C, typename T>
inline void operator+=(C& sus, T&& container)
  requires is_sus<C> and is_vector<T> and (std::is_same_v<base_type<typename C::key_type>, base_type<typename T::value_type> >)
  { std::for_each( container.begin(), container.end(), [&sus](auto&& element) { sus += std::forward<typename T::value_type>(element); } ); }

/*! \brief              Remove an element from a set, map, unordered set or unordered map
    \param  sus         destination set or unordered set
    \param  element     element to remove, or an iterator into <i>sus</i>
*/
template <typename C, typename T>
inline void operator-=(C& sus, const T& element)
  requires (is_sus<C> or is_mum<C>) and (std::is_same_v<typename C::key_type, base_type<T>> or std::is_same_v<base_type<T>, typename C::iterator>)
  { sus.erase(element); }

/*! \brief      Find the most common value in a vector
    \param  v   the target vector
    \return     the most common element in <i>v</i> and the number of times that it appears
*/
template <typename V>
  requires is_vector<V>
auto MOST_COMMON(const V& v) -> std::pair<typename V::value_type, int>
{ using RT = std::invoke_result_t<decltype(MOST_COMMON<V>), const V&>;
  
  std::unordered_map<typename V::value_type, int> um;

  FOR_ALL(v, [&um](const typename V::value_type& val) { um[val]++; });
  
  RT rv;
 
  FOR_ALL(um, [&rv](const RT& kcount) { if (kcount.second > rv.second) rv = kcount; } );

  return rv;
}

/*! \brief      Is an object a member of a set, an unordered_set, or a key in a map or unordered map?
    \param  s   set or unordered_set  to be tested
    \param  v   object to be tested for membership
    \return     Whether <i>v</i> is a member/key of <i>s</i>
*/
template <class T, class U>
inline bool contains(const T& s, const U& v)
  requires (is_sus<T> or is_mum<T>) and (std::is_same_v<typename T::key_type, U>)
  { return s.contains(v); }

/*! \brief          Is an object a member of a vector?
    \param  vec     vector  to be tested
    \param  val     object to be tested for membership
    \return         Whether <i>val</i> is a member of <i>vec</i>
*/  
template <class T, class U>
inline bool contains(const T& vec, const U& val)
  requires (is_vector<T> or std::derived_from<T, std::vector<typename T::value_type>>) and std::is_same_v<typename T::value_type, U> 
  { return std::ranges::find(vec, val) != vec.cend(); }
  
// define a hash function for pairs
// http://stackoverflow.com/questions/13485979/hash-function-of-unordered-set/13486174#13486174
// http://www.cplusplus.com/reference/functional/hash/
// https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
namespace std
{ template <typename T, typename U>
  struct hash< std::pair<T, U> >

  { using result_type = size_t;

    result_type operator()( const std::pair<T, U>& k ) const
    { result_type res { 17 };
            
      res = res * 31 + hash<T>()( k.first );
      res = res * 31 + hash<U>()( k.second );
      
      return res;
    }
  };
}

/*! \brief          Create a vector that is a transformed version of another vector
    \param  vec     original vector
    \param  op      transformation operation to be applied to each member of <i>vec</i>
    \return         <i>vec</i> with the operation <i>op</i> applied to each member
*/
template<typename T, typename OP>
auto VTV(const std::vector<T>& vec, OP op) -> std::vector<std::invoke_result_t<decltype(op), T>>
{ auto r { vec | std::views::transform(op) };

  std::vector<std::invoke_result_t<decltype(op), T>> rv;
  
  if constexpr(is_vector<decltype(rv)>)
    rv.reserve(vec.size());

  std::ranges::move(r, std::back_inserter(rv));
 
  return rv;
}

/*! \brief          Sum all the elements in a vector, with an optional initial value
    \param  vec     original vector
    \param  init    initial value of the output
    \return         the sum of <i>init</i> and all the elements in <i>vec</i>
*/
template<typename T>
T ACCUMULATE(const std::vector<T>& vec, const T init = 0)
{ T sum { init };
      
  return std::accumulate(vec.cbegin(), vec.cend(), sum);
}

/*! \brief          Sum all the elements in a vector
    \param  vec     vector to sum
    \return         the sum of all the elements in <i>vec</i>
*/
template<typename T>
inline T SUM(const std::vector<T>& vec)
  { return ACCUMULATE(vec); }

/*! \brief          Convert a vector of vectors to a single vector
    \param  vv      vector of vectors
    \return         the contents of <i>vv</i> in a single vector
    
    Returns an empty vector if all the inner vectors in <i>vv</i> are not the same length
*/
template <typename T>
std::vector<T> FLATTEN(std::vector<std::vector<T>>&& vv)
{ std::vector<T> rv;

  if (vv.empty())
    return rv;
    
  const size_t len0 { vv[0].size() };
  
  if (len0 == 0)
    return rv;
    
  for (size_t n { 1 }; n < vv.size(); ++n)
    if (vv[n].size() != len0)
      return rv;

// serially append vectors to rv  
  for (auto&& v : vv)
    rv += std::move(v);
    
  return rv;
}

/*! \brief          Build a two-dimensional vector from a one-dimensional one
    \param  v       the original one-dimensional vector
    \param  n_cols  the number of columns in the output vector
    \return         <i>v</i>, rearranged into a vector (rows) of vectors (columns)
*/
template <typename T>
std::vector<std::vector<T>> BUILD_VV(const std::vector<T>& v, const size_t n_cols)
{ std::vector<std::vector<T>> rv;

  if ( (v.size() % n_cols) != 0)
    return rv;
  
  const size_t n_rows { v.size() / n_cols };
  
  size_t offset { 0 };
  
  for (size_t transfer { 1 }; transfer <= n_rows; ++transfer)
  { rv += std::vector<T>(v.begin() + offset, v.begin() + offset + n_cols); // does this invalidate v?
    offset += n_cols;
  }
  
  return rv;
}

/*! \brief       Given a single index into an X * Y vector, return the X and Y indices
    \param  SX   the length of a row (i.e., the number of columns)
    \param  idx  the index into a vector
    \return      the column index (i.e., x) and row index (i.e., y) equivalent to <i>idx</i>
    
    Assumes that the vector is laid out in the order: Y=0: X[0]...x[N-1]; Y=1: X[0]...x[N-1]; etc.
*/
inline std::pair<size_t, size_t> INDICES(const size_t SX, const size_t idx)
  { return { (idx % SX), (idx / SX) }; }

/*! \brief      Given a pair of indices into an X * Y vector laid out as a single vector, return the equivalent single index
    \param  SX  the length of a row (i.e., the number of columns)
    \param  X   the x value
    \param  Y   the y value
    \return     the equivalent index into <i>idx</i>
    
    Assumes that the vector is laid out in the order: Y=0: X[0]...x[N-1]; Y=1: X[0]...x[N-1]; etc.
*/
inline size_t INDEX(const size_t SX, const size_t X, const size_t Y)
  { return (Y * SX) + X; }

/*! \brief          Calculate the means of columns in a two-dimensional vector laid out as a single vector
    \param  invec   the input vector
    \param  n_cols  the number of columns in the vector
    \return         a vector of length <i>n_cols</i> each element of which is the sum of the equivalent columns in <i>invec</i>
    
    Assumes that the vector is laid out in the order: Y=0: X[0]...x[N-1]; Y=1: X[0]...x[N-1]; etc.
*/
template <typename T>
std::vector<T> COLUMNS_MEAN(const std::vector<T>& invec, const size_t n_cols)
{ std::vector<T> rv;

  rv.reserve(n_cols);

  const size_t n_rows { invec.size() / n_cols };

//      Y             X
  std::vector<std::vector<T>> v2(n_rows, std::vector<T>(n_cols));

// create a real two-dimensional vector from invec  
  for (size_t idx { 0 }; idx < invec.size(); ++idx)
  { const auto [ x, y ] { INDICES(n_cols, idx) };

    v2[y][x] = invec[idx];
  }

// now calculate the means and append them to tv  
  for (size_t n_col { 0 }; n_col < n_cols; ++n_col)
  { auto sum { 0 };
  
    for (size_t n_row { 0 }; n_row < n_rows; ++n_row)
      sum += v2[n_row][n_col];
      
    const auto mean { sum / n_rows };
    
    rv += std::move(mean);
  }
  
  return rv;
}

/*! \brief          Calculate the means of columns in a two-dimensional vector
    \param  invv    the input vector
    \return         a vector, each element of which is the sum of the equivalent columns in <i>invec</i>
*/
template <typename T>
std::vector<T> COLUMNS_MEAN(const std::vector<std::vector<T>>& invv)
{ std::vector<T> rv;

  if (invv.empty())
    return rv;

  const size_t n_rows { invv.size() };
  const size_t n_cols { invv[0].size() };
  
  rv.reserve(n_cols);

  for (size_t n_col { 0 }; n_col < n_cols; ++n_col)
  { T sum { 0 };
  
    for (size_t n_row { 0 }; n_row < n_rows; ++n_row)
      sum += invv[n_row][n_col];
      
    const T mean { sum / n_rows };
    
    rv += std::move(mean);
  }
  
  return rv;
}

/*! \brief              Return rows of a two-dimensional vector
    \param  vv          the input vector
    \param  start_row   the number of the starting row to be returned
    \param  n_rows      the number of rows to be returned
    \return             rows numbered <i>start_row</i> to <i>start_row</i> + <i>n_rows</i> - 1 of <i>vv</i>
    
    It is not an error to try to read past the end of <i>vv</i>
*/
template <typename T>
std::vector<std::vector<T>> GET_ROWS(const std::vector<std::vector<T>>& vv, const size_t start_row, const size_t n_rows)
{ std::vector<std::vector<T>> rv;

  for (size_t counter { 1 }; ((counter <= n_rows) and (start_row + counter - 1) < vv.size()); ++counter)  // not an error to try to read past the end; simply ignore if n_rows is too large
    rv += vv[start_row + counter - 1];
    
  return rv;
}

/*! \brief              Calculate the mean of part of a vector
    \param  invec       the input vector
    \param  start_idx   initial index into <i>invec</i>
    \param  nv          number of values to be used 
    \return             the mean of the velues in the target portion of <i>invec</i>
    
    Ignores any attempts to use values beyond the end of the vector
*/
template <typename T>
T V_MEAN(const std::vector<T>& invec, const size_t start_idx, const size_t nv)
{ T sum { 0 };

  size_t n_processed { 0 };

  for (size_t idx { start_idx }; (n_processed < nv) and (idx < invec.size()); ++idx)
  { sum += invec[idx];
    n_processed++;
  }

  return sum / n_processed;
}

/*! \brief      Return the fractional part of a value
    \param  v   value
    \return     the fractional part of <i>v</i>
    
    The sign of the result is (according to the definition of modf) the same as the sign of <i>v</i>
*/
template <typename T>
T frac_part(const T v)
{ double int_part;           // to be thrown away

  return modf(v, &int_part);
}

/*! \brief      Convert a range to a vector
    \param  r   range
    \return     <i>r</i> as a vector
*/
#if 0
template <std::ranges::range R>
auto RANGE_VECTOR(R&& r)
{ auto r_common { r | std::views::common };

  return std::vector(r_common.begin(), r_common.end());
}
#endif

/*! \brief      Convert a range to a list
    \param  r   range
    \return     <i>r</i> as a list
*/
template <std::ranges::range R>
auto RANGE_LIST(R&& r)
{ auto r_common { std::forward<R>(r) | std::views::common };

  return std::list(r_common.begin(), r_common.end());
}

/*! \brief      Convert a range to a set
    \param  r   range
    \return     <i>r</i> as a set
*/
template <std::ranges::range R>
auto RANGE_SET(R&& r)
{ auto r_common { std::forward<R>(r) | std::views::common };

  return std::set(r_common.begin(), r_common.end());
}

#if 0
template <std::ranges::range R>
  requires is_string<typename R::value_type>
auto RANGE_SET(R&& r)
{ auto r_common { std::forward<R>(r) | std::views::common };

  return STRING_SET(r_common.begin(), r_common.end());
}
#endif

/*! \brief      Convert a range to an unordered set
    \param  r   range
    \return     <i>r</i> as an unordered set
*/
template <std::ranges::range R>
auto RANGE_USET(R&& r)
{ auto r_common { std::forward<R>(r) | std::views::common };

  return std::unordered_set(r_common.begin(), r_common.end());
}

/*! \brief      Convert a range to a particular container type
    \param  r   range
    \return     <i>r</i> as a particular container
*/
template <typename RT, std::ranges::range R>
auto RANGE_CONTAINER(R&& r) -> RT
{ auto r_common { std::forward<R>(r) | std::views::common };

  return RT(r_common.begin(), r_common.end());
}

/*! \brief          Convert a value to a value in the principal domain
    \param  value   value to shift to the principal range
    \param  minv    minimum value in the principal range
    \param  maxv    maximum value in the principal range
    \return         <i>value</i> shifted to the range [minv, maxv)
*/
template <typename T, typename U, typename V>
  requires ( (std::convertible_to<U, T> and std::convertible_to<V, T>) and (!is_int<T>) )
auto to_principal_domain(T value, U minv, V maxv) -> T
{ const T min_T { static_cast<T>(minv) };
  const T max_T { static_cast<T>(maxv) };

  if (value < min_T)
   { const T width          { max_T - min_T };
     const T initial_skew   { (min_T - value) };
     const T skew_in_widths { initial_skew / width };
     
     const T rv { value + (static_cast<int>(skew_in_widths) + 1) * width };  // the cast is faster than using floor()
    
     if (rv < min_T or rv >= max_T)
     { std::cerr << "value < min_T ERROR in to_principal_domain(" << value << ", " << minv << ", " << maxv << ") : " << rv << std::endl;
       std::cerr << "width = " << width << "; initial_skew = " << initial_skew << "; skew_in_widths = " << skew_in_widths << std::endl;
       throw std::exception();
       exit(-1);
     }
    
     return rv;
  }   

  if (value >= max_T)
  { const T width          { max_T - min_T };
    const T skew_in_widths { (value - max_T) / width };
    const T rv { minv + (skew_in_widths - static_cast<int>(skew_in_widths)) * width };  // the cast is faster than using floor()
    
    if (rv < min_T or rv >= max_T)
    { std::cerr << "value >= max_T ERROR in to_principal_domain(" << value << ", " << minv << ", " << maxv << ") : " << rv << std::endl;
      throw std::exception();
      exit(-1);
    }
   
    return rv;
  }
  
  return value;
}

#if 0
// for integer types
template <typename T, typename U, typename V>
  requires ( (std::is_convertible_v<U, T> and std::is_convertible_v<V, T>) and (is_int<T>) )
auto to_principal_domain(T value, U minv, V maxv) -> T
{ std::cout << "int to_principal_domain: " << value << ", " << minv << ", " << maxv << std::endl;
  
  const T min_T { static_cast<T>(minv) };
  const T max_T { static_cast<T>(maxv) };

  if (value < min_T)
   { const T width          { max_T - min_T };
     const T initial_skew   { (min_T - value) };
     const double skew_in_widths { static_cast<double>(initial_skew) / width };
     
     const T rv { value + (static_cast<int>(skew_in_widths) + 1) * width };  // the cast is faster than using floor()
     
     if (rv < min_T or rv >= max_T)
     { std::cerr << "value < min_T ERROR in to_principal_domain(" << value << ", " << minv << ", " << maxv << ") : " << rv << std::endl;
       throw std::exception();
       exit(-1);
     }
     
     std::cout << "small value; to_principal_domain returning " << rv << std::endl;
    
     return rv;
  }   

  if (value > max_T)        // NB >, not >= for ints
  { const T      width          { max_T - min_T };
    const double skew_in_widths { static_cast<double>((value - max_T)) / width };
    const T      rv             { minv + static_cast<T>((skew_in_widths - static_cast<int>(skew_in_widths)) * width + 0.5) };  // the cast is faster than using floor()
    
    if (rv < min_T or rv >= max_T)
    { std::cerr << "value >= max_T ERROR in to_principal_domain(" << value << ", " << minv << ", " << maxv << ") : " << rv << std::endl;
      throw std::exception();
      exit(-1);
    }
   
    std::cout << "large value; to_principal_domain returning " << rv << std::endl;
   
    return rv;
  }
  
  std::cout << "to_principal_domain returning " << value << std::endl;
  
  return value;
}
#endif

/*! \brief          Convert a value to exponent and a mantissa, base 10
    \param  val     val to convert
    \return         exponent and mantissa
*/
template <typename T>
std::pair<T, T> decompose10(const T val)
{ const auto e { floor(log10(fabs(val))) };
  const auto s { (val) / pow(10, e) };
  
  return { e, s };
}

// accumulate in place
template <typename InputIt, typename T, typename OP>
auto ACCIP(InputIt first, InputIt last, T init, OP op) -> T
{
    for (; first != last; ++first)
    { init += (op(*first));
    }
    return init;
}

// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
 template <typename T> inline constexpr
 int signum(T x, std::false_type is_signed) {
     return T(0) < x;
 }

 template <typename T> inline constexpr
 int signum(T x, std::true_type is_signed) {
     return (T(0) < x) - (x < T(0));
 }

 template <typename T> inline constexpr
 int signum(T x) {
     return signum(x, std::is_signed<T>());
 }

/*! \brief      Convert a value from radians to an integral number of decidegrees
    \param  r   angle in radians
    \return     <i>r</i>, converted to decidegrees
*/
template <typename I>
auto TOD10(const std::floating_point auto v) -> I
{ using FP = decltype(v);

  const bool is_negative { v < 0 };
  const FP   pv          { is_negative ? -v : v };    // +ve value
  const FP   r10         { pv * (180 / 3.14159265358979) * 10 };
  const FP   frac        { r10 - floor(r10) };
  
  if (fabs(frac - 0.5) < 0.01)      // close to an edge value; round to even
  { const int fl { static_cast<int>(floor(r10)) };
    const int last_digit { fl % 10 };
    
    switch (last_digit)
    { case 0 :
      case 2 :
      case 4 :
      case 6 :
      case 8 :
        return (is_negative ? -fl : fl);
        
      case 1 :
      case 3 :
      case 5 :
      case 7 :
      case 9 :
        return (is_negative ? -(fl + 1) : (fl + 1));
        
      default :
      { std::cerr << "ERROR in TOD10" << std::endl;
        exit(-1);
      }
    }
  }
  
  return static_cast<I>(round(v * (180 / 3.14159265358979) * 10)); 
}

/*! \brief      Round upwards to the next higher sensible number
    \param  x   number to round
    \return     <i>x</i> rounded upward to what should be a sensible number
*/
template <typename T>
T auto_round(const T x)
{ const double lg { log10(x) };

  if (lg == static_cast<int>(lg))
    return x;

  const int q    { static_cast<int>(lg) };
  const int fact { static_cast<int>(x / pow(10, q)) + 1 };

  return static_cast<T>(fact * pow(10, q));
}

/*! \brief      Calculate and return the pairwise sum of two vectors
    \param      v1  first vector
    \param      v2  second vector
    \return     vector with each element equal to the sum of the equivalent elements in <i>v1</i> and <i>v2</i>
*/
template <typename T>
std::vector<T> PAIRWISE_SUM(const std::vector<T>& v1, const std::vector<T>& v2)
{ const size_t common_size { std::min(v1.size(), v2.size()) };

  std::vector<T> rv;
  rv.reserve(common_size);
  
  for (size_t idx { 0 }; idx < common_size; ++idx)
    rv += v1[idx] + v2[idx];
    
  return rv;
}

/*! \brief      Return a pair with members switched
    \param      tu  original pair
    \return     <i>tu</i> with the two elements swapped
*/
template <typename T, typename U>
inline auto INVERT_PAIR(const std::pair<T, U>& tu) -> std::pair<U, T>
  { return std::pair { tu.second, tu.first }; }

// temporary, because .at() does not support heterogeneous lookup
template <typename M, typename SV>
  requires is_mum<M> and is_string<typename M::key_type> and is_string_view<SV>
typename M::mapped_type AT(const M& m, const SV sv)
{ const auto it { m.find(sv) };

  if (it == m.end())
    throw std::exception();

  return it->second;
}

template <typename M, typename S>
  requires is_mum<M> and is_string<typename M::key_type> and is_string<S>
inline typename M::mapped_type AT(const M& m, const S& s)
  { return m.at(s); }

#endif    // MACROS_H
