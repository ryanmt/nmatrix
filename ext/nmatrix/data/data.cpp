/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2013, Ruby Science Foundation
// NMatrix is Copyright (c) 2013, Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == data.cpp
//
// Functions and data for dealing the data types.

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Project Includes
 */

#include "types.h"
#include "data.h"

/*
 * Global Variables
 */

namespace nm {
  const char* const EWOP_OPS[nm::NUM_EWOPS] = {
    "+",
    "-",
    "*",
    "/",
    "**",
    "%",
    "==",
    "!=",
    "<",
    ">",
    "<=",
    ">="
  };

  const std::string EWOP_NAMES[nm::NUM_EWOPS] = {
    "add",
    "sub",
    "mul",
    "div",
    "pow",
    "mod",
    "eqeq",
    "neq",
    "lt",
    "gt",
    "leq",
    "geq"
  };


  template <typename Type>
  Complex<Type>::Complex(const RubyObject& other) {
    switch(TYPE(other.rval)) {
    case T_COMPLEX:
      r = NUM2DBL(rb_funcall(other.rval, rb_intern("real"), 0));
      i = NUM2DBL(rb_funcall(other.rval, rb_intern("imag"), 0));
      break;
    case T_FLOAT:
    case T_RATIONAL:
    case T_FIXNUM:
    case T_BIGNUM:
      r = NUM2DBL(other.rval);
      i = 0.0;
      break;
    default:
      rb_raise(rb_eTypeError, "not sure how to convert this type of VALUE to a complex");
    }
  }


  template <typename Type>
  Rational<Type>::Rational(const RubyObject& other) {
    switch (TYPE(other.rval)) {
    case T_RATIONAL:
      n = NUM2LONG(rb_funcall(other.rval, rb_intern("numerator"), 0));
      d = NUM2LONG(rb_funcall(other.rval, rb_intern("denominator"), 0));
      break;
    case T_FIXNUM:
    case T_BIGNUM:
      n = NUM2LONG(other.rval);
      d = 1;
      break;
    case T_COMPLEX:
    case T_FLOAT:
      rb_raise(rb_eTypeError, "cannot convert float to a rational");
      break;
    default:
      rb_raise(rb_eTypeError, "not sure how to convert this type of VALUE to a rational");
    }
  }


} // end of namespace nm
extern "C" {

const char* const DTYPE_NAMES[nm::NUM_DTYPES] = {
	"byte",
	"int8",
	"int16",
	"int32",
	"int64",
	"float32",
	"float64",
	"complex64",
	"complex128",
	"rational32",
	"rational64",
	"rational128",
	"object"
};


const size_t DTYPE_SIZES[nm::NUM_DTYPES] = {
	sizeof(uint8_t),
	sizeof(int8_t),
	sizeof(int16_t),
	sizeof(int32_t),
	sizeof(int64_t),
	sizeof(float32_t),
	sizeof(float64_t),
	sizeof(nm::Complex64),
	sizeof(nm::Complex128),
	sizeof(nm::Rational32),
	sizeof(nm::Rational64),
	sizeof(nm::Rational128),
	sizeof(nm::RubyObject)
};


const nm::dtype_t Upcast[nm::NUM_DTYPES][nm::NUM_DTYPES] = {
  { nm::BYTE, nm::INT16, nm::INT16, nm::INT32, nm::INT64, nm::FLOAT32, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::RATIONAL32, nm::RATIONAL64, nm::RATIONAL128, nm::RUBYOBJ},
  { nm::INT16, nm::INT8, nm::INT16, nm::INT32, nm::INT64, nm::FLOAT32, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::RATIONAL32, nm::RATIONAL64, nm::RATIONAL128, nm::RUBYOBJ},
  { nm::INT16, nm::INT16, nm::INT16, nm::INT32, nm::INT64, nm::FLOAT32, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::RATIONAL32, nm::RATIONAL64, nm::RATIONAL128, nm::RUBYOBJ},
  { nm::INT32, nm::INT32, nm::INT32, nm::INT32, nm::INT64, nm::FLOAT32, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::RATIONAL32, nm::RATIONAL64, nm::RATIONAL128, nm::RUBYOBJ},
  { nm::INT64, nm::INT64, nm::INT64, nm::INT64, nm::INT64, nm::FLOAT32, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::RATIONAL32, nm::RATIONAL64, nm::RATIONAL128, nm::RUBYOBJ},
  { nm::FLOAT32, nm::FLOAT32, nm::FLOAT32, nm::FLOAT32, nm::FLOAT32, nm::FLOAT32, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::FLOAT64, nm::FLOAT64, nm::FLOAT64, nm::RUBYOBJ},
  { nm::FLOAT64, nm::FLOAT64, nm::FLOAT64, nm::FLOAT64, nm::FLOAT64, nm::FLOAT64, nm::FLOAT64, nm::COMPLEX128, nm::COMPLEX128, nm::FLOAT64, nm::FLOAT64, nm::FLOAT64, nm::RUBYOBJ},
  { nm::COMPLEX64, nm::COMPLEX64, nm::COMPLEX64, nm::COMPLEX64, nm::COMPLEX64, nm::COMPLEX64, nm::COMPLEX128, nm::COMPLEX64, nm::COMPLEX128, nm::COMPLEX64, nm::COMPLEX64, nm::COMPLEX64, nm::RUBYOBJ},
  { nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::COMPLEX128, nm::RUBYOBJ},
  { nm::RATIONAL32, nm::RATIONAL32, nm::RATIONAL32, nm::RATIONAL32, nm::RATIONAL32, nm::FLOAT64, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::RATIONAL32, nm::RATIONAL64, nm::RATIONAL128, nm::RUBYOBJ},
  { nm::RATIONAL64, nm::RATIONAL64, nm::RATIONAL64, nm::RATIONAL64, nm::RATIONAL64, nm::FLOAT64, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::RATIONAL64, nm::RATIONAL64, nm::RATIONAL128, nm::RUBYOBJ},
  { nm::RATIONAL128, nm::RATIONAL128, nm::RATIONAL128, nm::RATIONAL128, nm::RATIONAL128, nm::FLOAT64, nm::FLOAT64, nm::COMPLEX64, nm::COMPLEX128, nm::RATIONAL128, nm::RATIONAL128, nm::RATIONAL128, nm::RUBYOBJ},
  { nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ, nm::RUBYOBJ}
};


/*
 * Forward Declarations
 */

/*
 * Functions
 */

/*
 * Converts a RubyObject
 */
void rubyval_to_cval(VALUE val, nm::dtype_t dtype, void* loc) {
  using namespace nm;
	switch (dtype) {
		case BYTE:
			*reinterpret_cast<uint8_t*>(loc)			= static_cast<uint8_t>(RubyObject(val));
			break;

		case INT8:
			*reinterpret_cast<int8_t*>(loc)				= static_cast<int8_t>(RubyObject(val));
			break;

		case INT16:
			*reinterpret_cast<int16_t*>(loc)			= static_cast<int16_t>(RubyObject(val));
			break;

		case INT32:
			*reinterpret_cast<int32_t*>(loc)			= static_cast<int32_t>(RubyObject(val));
			break;

		case INT64:
			*reinterpret_cast<int64_t*>(loc)			= static_cast<int64_t>(RubyObject(val));
			break;

		case FLOAT32:
			*reinterpret_cast<float32_t*>(loc)		= static_cast<float32_t>(RubyObject(val));
			break;

		case FLOAT64:
			*reinterpret_cast<float64_t*>(loc)		= static_cast<float64_t>(RubyObject(val));
			break;

		case COMPLEX64:
			*reinterpret_cast<Complex64*>(loc)		= RubyObject(val).to<Complex64>();
			break;

		case COMPLEX128:
			*reinterpret_cast<Complex128*>(loc)		= RubyObject(val).to<Complex64>();
			break;

		case RATIONAL32:
			*reinterpret_cast<Rational32*>(loc)		= RubyObject(val).to<Rational32>();
			break;

		case RATIONAL64:
			*reinterpret_cast<Rational64*>(loc)		= RubyObject(val).to<Rational64>();
			break;

		case RATIONAL128:
			*reinterpret_cast<Rational128*>(loc)	= RubyObject(val).to<Rational128>();
			break;

		case RUBYOBJ:
		  *reinterpret_cast<VALUE*>(loc)        = val;
			//rb_raise(rb_eTypeError, "Attempting a bad conversion from a Ruby value.");
			break;

	  default:
	    rb_raise(rb_eTypeError, "Attempting a bad conversion.");
	    break;
	}
}

/*
 * Create a RubyObject from a regular C value (given a dtype). Does not return a VALUE! To get a VALUE, you need to
 * look at the rval property of what this function returns.
 */
nm::RubyObject rubyobj_from_cval(void* val, nm::dtype_t dtype) {
  using namespace nm;
	switch (dtype) {
		case BYTE:
			return RubyObject(*reinterpret_cast<uint8_t*>(val));

		case INT8:
			return RubyObject(*reinterpret_cast<int8_t*>(val));

		case INT16:
			return RubyObject(*reinterpret_cast<int16_t*>(val));

		case INT32:
			return RubyObject(*reinterpret_cast<int32_t*>(val));

		case INT64:
			return RubyObject(*reinterpret_cast<int64_t*>(val));

		case FLOAT32:
			return RubyObject(*reinterpret_cast<float32_t*>(val));

		case FLOAT64:
			return RubyObject(*reinterpret_cast<float64_t*>(val));

		case COMPLEX64:
			return RubyObject(*reinterpret_cast<Complex64*>(val));

		case COMPLEX128:
			return RubyObject(*reinterpret_cast<Complex128*>(val));
			
		case RATIONAL32:
			return RubyObject(*reinterpret_cast<Rational32*>(val));
			
		case RATIONAL64:
			return RubyObject(*reinterpret_cast<Rational64*>(val));
			
		case RATIONAL128:
			return RubyObject(*reinterpret_cast<Rational128*>(val));

	  default:
	    throw;
	    rb_raise(nm_eDataTypeError, "Conversion to RubyObject requested from unknown/invalid data type (did you try to convert from a VALUE?)");
	}
	return Qnil;
}



/*
 * Allocate and return a piece of data of the correct dtype, converted from a
 * given RubyObject.
 */
void* rubyobj_to_cval(VALUE val, nm::dtype_t dtype) {
  size_t size =  DTYPE_SIZES[dtype];
  void* ret_val = ALLOC_N(char, size);

  rubyval_to_cval(val, dtype, ret_val);

  return ret_val;
}


void nm_init_data() {
    nm::RubyObject obj(INT2FIX(1));
    nm::Rational32 x(obj);
    nm::Rational64 y(obj);
    nm::Rational128 z(obj);
    nm::Complex64 a(obj);
    nm::Complex128 b(obj);
}


} // end of extern "C" block
