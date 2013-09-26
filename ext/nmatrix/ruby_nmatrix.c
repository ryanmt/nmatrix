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
// == ruby_nmatrix.c
//
// Ruby-facing NMatrix C functions. Not compiled directly -- included
// into nmatrix.cpp.
//

/*
 * Forward Declarations
 */

static VALUE nm_init(int argc, VALUE* argv, VALUE nm);
static VALUE nm_init_copy(VALUE copy, VALUE original);
static VALUE nm_init_transposed(VALUE self);
static VALUE nm_read(int argc, VALUE* argv, VALUE self);
static VALUE nm_write(int argc, VALUE* argv, VALUE self);
static VALUE nm_init_yale_from_old_yale(VALUE shape, VALUE dtype, VALUE ia, VALUE ja, VALUE a, VALUE from_dtype, VALUE nm);
static VALUE nm_alloc(VALUE klass);
static VALUE nm_dtype(VALUE self);
static VALUE nm_stype(VALUE self);
static VALUE nm_default_value(VALUE self);
static size_t effective_dim(STORAGE* s);
static VALUE nm_effective_dim(VALUE self);
static VALUE nm_dim(VALUE self);
static VALUE nm_offset(VALUE self);
static VALUE nm_shape(VALUE self);
static VALUE nm_supershape(VALUE self);
static VALUE nm_capacity(VALUE self);
static VALUE nm_each_with_indices(VALUE nmatrix);
static VALUE nm_each_stored_with_indices(VALUE nmatrix);
static VALUE nm_each_ordered_stored_with_indices(VALUE nmatrix);

static SLICE* get_slice(size_t dim, int argc, VALUE* arg, size_t* shape);
static VALUE nm_xslice(int argc, VALUE* argv, void* (*slice_func)(const STORAGE*, SLICE*), void (*delete_func)(NMATRIX*), VALUE self);
static VALUE nm_mset(int argc, VALUE* argv, VALUE self);
static VALUE nm_mget(int argc, VALUE* argv, VALUE self);
static VALUE nm_mref(int argc, VALUE* argv, VALUE self);
static VALUE nm_is_ref(VALUE self);

static VALUE is_symmetric(VALUE self, bool hermitian);

static VALUE nm_guess_dtype(VALUE self, VALUE v);
static VALUE nm_min_dtype(VALUE self, VALUE v);

/*
 * Macro defines an element-wise accessor function for some operation.
 *
 * This is only responsible for the Ruby accessor! You still have to write the actual functions, obviously.
 */
#define DEF_ELEMENTWISE_RUBY_ACCESSOR(oper, name)                 \
static VALUE nm_ew_##name(VALUE left_val, VALUE right_val) {  \
  return elementwise_op(nm::EW_##oper, left_val, right_val);  \
}

/*
 * Macro declares a corresponding accessor function prototype for some element-wise operation.
 */
#define DECL_ELEMENTWISE_RUBY_ACCESSOR(name)    static VALUE nm_ew_##name(VALUE left_val, VALUE right_val);

DECL_ELEMENTWISE_RUBY_ACCESSOR(add)
DECL_ELEMENTWISE_RUBY_ACCESSOR(subtract)
DECL_ELEMENTWISE_RUBY_ACCESSOR(multiply)
DECL_ELEMENTWISE_RUBY_ACCESSOR(divide)
DECL_ELEMENTWISE_RUBY_ACCESSOR(power)
DECL_ELEMENTWISE_RUBY_ACCESSOR(mod)
DECL_ELEMENTWISE_RUBY_ACCESSOR(eqeq)
DECL_ELEMENTWISE_RUBY_ACCESSOR(neq)
DECL_ELEMENTWISE_RUBY_ACCESSOR(lt)
DECL_ELEMENTWISE_RUBY_ACCESSOR(gt)
DECL_ELEMENTWISE_RUBY_ACCESSOR(leq)
DECL_ELEMENTWISE_RUBY_ACCESSOR(geq)

static VALUE elementwise_op(nm::ewop_t op, VALUE left_val, VALUE right_val);

static VALUE nm_symmetric(VALUE self);
static VALUE nm_hermitian(VALUE self);

static VALUE nm_eqeq(VALUE left, VALUE right);

static VALUE matrix_multiply_scalar(NMATRIX* left, VALUE scalar);
static VALUE matrix_multiply(NMATRIX* left, NMATRIX* right);
static VALUE nm_multiply(VALUE left_v, VALUE right_v);
static VALUE nm_det_exact(VALUE self);
static VALUE nm_complex_conjugate_bang(VALUE self);

static nm::dtype_t	interpret_dtype(int argc, VALUE* argv, nm::stype_t stype);
static void*		interpret_initial_value(VALUE arg, nm::dtype_t dtype);
static size_t*	interpret_shape(VALUE arg, size_t* dim);
static nm::stype_t	interpret_stype(VALUE arg);

/* Singleton methods */
static VALUE nm_upcast(VALUE self, VALUE t1, VALUE t2);


#ifdef BENCHMARK
static double get_time(void);
#endif

///////////////////
// Ruby Bindings //
///////////////////

void Init_nmatrix() {


	///////////////////////
	// Class Definitions //
	///////////////////////

	cNMatrix = rb_define_class("NMatrix", rb_cObject);
	//cNVector = rb_define_class("NVector", cNMatrix);

	// Special exceptions

	/*
	 * Exception raised when there's a problem with data.
	 */
	nm_eDataTypeError    = rb_define_class("DataTypeError", rb_eStandardError);

	/*
	 * Exception raised when something goes wrong with the storage of a matrix.
	 */
	nm_eStorageTypeError = rb_define_class("StorageTypeError", rb_eStandardError);

	///////////////////
	// Class Methods //
	///////////////////

	rb_define_alloc_func(cNMatrix, nm_alloc);

	///////////////////////
  // Singleton Methods //
  ///////////////////////

	rb_define_singleton_method(cNMatrix, "upcast", (METHOD)nm_upcast, 2); /* in ext/nmatrix/nmatrix.cpp */
	rb_define_singleton_method(cNMatrix, "guess_dtype", (METHOD)nm_guess_dtype, 1);
	rb_define_singleton_method(cNMatrix, "min_dtype", (METHOD)nm_min_dtype, 1);

	//////////////////////
	// Instance Methods //
	//////////////////////

	rb_define_method(cNMatrix, "initialize", (METHOD)nm_init, -1);
	rb_define_method(cNMatrix, "initialize_copy", (METHOD)nm_init_copy, 1);
	rb_define_singleton_method(cNMatrix, "read", (METHOD)nm_read, -1);

	rb_define_method(cNMatrix, "write", (METHOD)nm_write, -1);

	// Technically, the following function is a copy constructor.
	rb_define_method(cNMatrix, "transpose", (METHOD)nm_init_transposed, 0);

	rb_define_method(cNMatrix, "dtype", (METHOD)nm_dtype, 0);
	rb_define_method(cNMatrix, "stype", (METHOD)nm_stype, 0);
	rb_define_method(cNMatrix, "cast_full",  (METHOD)nm_cast, 3);
	rb_define_method(cNMatrix, "default_value", (METHOD)nm_default_value, 0);
	rb_define_protected_method(cNMatrix, "__list_default_value__", (METHOD)nm_list_default_value, 0);
	rb_define_protected_method(cNMatrix, "__yale_default_value__", (METHOD)nm_yale_default_value, 0);

	rb_define_method(cNMatrix, "[]", (METHOD)nm_mref, -1);
	rb_define_method(cNMatrix, "slice", (METHOD)nm_mget, -1);
	rb_define_method(cNMatrix, "[]=", (METHOD)nm_mset, -1);
	rb_define_method(cNMatrix, "is_ref?", (METHOD)nm_is_ref, 0);
	rb_define_method(cNMatrix, "dimensions", (METHOD)nm_dim, 0);
	rb_define_method(cNMatrix, "effective_dimensions", (METHOD)nm_effective_dim, 0);

	rb_define_protected_method(cNMatrix, "__list_to_hash__", (METHOD)nm_to_hash, 0); // handles list and dense, which are n-dimensional

	rb_define_method(cNMatrix, "shape", (METHOD)nm_shape, 0);
	rb_define_method(cNMatrix, "supershape", (METHOD)nm_supershape, 0);
	rb_define_method(cNMatrix, "offset", (METHOD)nm_offset, 0);
	rb_define_method(cNMatrix, "det_exact", (METHOD)nm_det_exact, 0);
	rb_define_method(cNMatrix, "complex_conjugate!", (METHOD)nm_complex_conjugate_bang, 0);

	rb_define_protected_method(cNMatrix, "__dense_each__", (METHOD)nm_dense_each, 0);
	rb_define_protected_method(cNMatrix, "__dense_map__", (METHOD)nm_dense_map, 0);
	rb_define_protected_method(cNMatrix, "__dense_map_pair__", (METHOD)nm_dense_map_pair, 1);
	rb_define_method(cNMatrix, "each_with_indices", (METHOD)nm_each_with_indices, 0);
	rb_define_method(cNMatrix, "each_stored_with_indices", (METHOD)nm_each_stored_with_indices, 0);
	rb_define_method(cNMatrix, "each_ordered_stored_with_indices", (METHOD)nm_each_ordered_stored_with_indices, 0);
	rb_define_protected_method(cNMatrix, "__list_map_merged_stored__", (METHOD)nm_list_map_merged_stored, 2);
	rb_define_protected_method(cNMatrix, "__yale_map_merged_stored__", (METHOD)nm_yale_map_merged_stored, 2);
	rb_define_protected_method(cNMatrix, "__yale_map_stored__", (METHOD)nm_yale_map_stored, 0);
	rb_define_protected_method(cNMatrix, "__yale_stored_diagonal_each_with_indices__", (METHOD)nm_yale_stored_diagonal_each_with_indices, 0);
	rb_define_protected_method(cNMatrix, "__yale_stored_nondiagonal_each_with_indices__", (METHOD)nm_yale_stored_nondiagonal_each_with_indices, 0);

	rb_define_method(cNMatrix, "==",	  (METHOD)nm_eqeq,				1);

	rb_define_method(cNMatrix, "+",			(METHOD)nm_ew_add,			1);
	rb_define_method(cNMatrix, "-",			(METHOD)nm_ew_subtract,	1);
  rb_define_method(cNMatrix, "*",			(METHOD)nm_ew_multiply,	1);
	rb_define_method(cNMatrix, "/",			(METHOD)nm_ew_divide,		1);
  rb_define_method(cNMatrix, "**",    (METHOD)nm_ew_power,    1);
  rb_define_method(cNMatrix, "%",     (METHOD)nm_ew_mod,      1);

	rb_define_method(cNMatrix, "=~", (METHOD)nm_ew_eqeq, 1);
	rb_define_method(cNMatrix, "!~", (METHOD)nm_ew_neq, 1);
	rb_define_method(cNMatrix, "<=", (METHOD)nm_ew_leq, 1);
	rb_define_method(cNMatrix, ">=", (METHOD)nm_ew_geq, 1);
	rb_define_method(cNMatrix, "<", (METHOD)nm_ew_lt, 1);
	rb_define_method(cNMatrix, ">", (METHOD)nm_ew_gt, 1);

	/////////////////////////////
	// Helper Instance Methods //
	/////////////////////////////
	rb_define_protected_method(cNMatrix, "__yale_vector_set__", (METHOD)nm_vector_set, -1);

	/////////////////////////
	// Matrix Math Methods //
	/////////////////////////
	rb_define_method(cNMatrix, "dot",		(METHOD)nm_multiply,		1);

	rb_define_method(cNMatrix, "symmetric?", (METHOD)nm_symmetric, 0);
	rb_define_method(cNMatrix, "hermitian?", (METHOD)nm_hermitian, 0);

	rb_define_method(cNMatrix, "capacity", (METHOD)nm_capacity, 0);

	/////////////
	// Aliases //
	/////////////

	rb_define_alias(cNMatrix, "dim", "dimensions");
	rb_define_alias(cNMatrix, "effective_dim", "effective_dimensions");
	rb_define_alias(cNMatrix, "equal?", "eql?");

	///////////////////////
	// Symbol Generation //
	///////////////////////

	nm_init_ruby_constants();

	//////////////////////////
	// YaleFunctions module //
	//////////////////////////

	nm_init_yale_functions();

	/////////////////
	// BLAS module //
	/////////////////

	nm_math_init_blas();

	///////////////
	// IO module //
	///////////////
	nm_init_io();

	/////////////////////////////////////////////////
	// Force compilation of necessary constructors //
	/////////////////////////////////////////////////
	nm_init_data();
}


//////////////////
// Ruby Methods //
//////////////////


/*
 * Slice constructor.
 */
static SLICE* alloc_slice(size_t dim) {
  SLICE* slice = ALLOC(SLICE);
  slice->coords = ALLOC_N(size_t, dim);
  slice->lengths = ALLOC_N(size_t, dim);
  return slice;
}


/*
 * Slice destructor.
 */
static void free_slice(SLICE* slice) {
  xfree(slice->coords);
  xfree(slice->lengths);
  xfree(slice);
}


/*
 * Allocator.
 */
static VALUE nm_alloc(VALUE klass) {
  NMATRIX* mat = ALLOC(NMATRIX);
  mat->storage = NULL;

  // DO NOT MARK This STRUCT. It has no storage allocated, and no stype, so mark will do an invalid something.
  return Data_Wrap_Struct(klass, NULL, nm_delete, mat);
}

/*
 * Find the capacity of an NMatrix. The capacity only differs from the size for
 * Yale matrices, which occasionally allocate more space than they need. For
 * list and dense, capacity gives the number of elements in the matrix.
 *
 * If you call this on a slice, it may behave unpredictably. Most likely it'll
 * just return the original matrix's capacity.
 */
static VALUE nm_capacity(VALUE self) {
  VALUE cap;

  switch(NM_STYPE(self)) {
  case nm::YALE_STORE:
    cap = UINT2NUM(reinterpret_cast<YALE_STORAGE*>(NM_STORAGE_YALE(self)->src)->capacity);
    break;

  case nm::DENSE_STORE:
    cap = UINT2NUM(nm_storage_count_max_elements( NM_STORAGE_DENSE(self) ));
    break;

  case nm::LIST_STORE:
    cap = UINT2NUM(nm_list_storage_count_elements( NM_STORAGE_LIST(self) ));
    break;

  default:
    rb_raise(nm_eStorageTypeError, "unrecognized stype in nm_capacity()");
  }

  return cap;
}


/*
 * Mark function.
 */
void nm_mark(NMATRIX* mat) {
  STYPE_MARK_TABLE(mark)
  mark[mat->stype](mat->storage);
}


/*
 * Destructor.
 */
void nm_delete(NMATRIX* mat) {
  static void (*ttable[nm::NUM_STYPES])(STORAGE*) = {
    nm_dense_storage_delete,
    nm_list_storage_delete,
    nm_yale_storage_delete
  };
  ttable[mat->stype](mat->storage);

  xfree(mat);
}

/*
 * Slicing destructor.
 */
void nm_delete_ref(NMATRIX* mat) {
  static void (*ttable[nm::NUM_STYPES])(STORAGE*) = {
    nm_dense_storage_delete_ref,
    nm_list_storage_delete_ref,
    nm_yale_storage_delete_ref
  };
  ttable[mat->stype](mat->storage);

  xfree(mat);
}

/*
 * Register the addresses of an array of VALUEs with the gc to avoid collection
 * while using them internally.
 */
void nm_register_values(VALUE* values, size_t n) {
  if (values) {
    for (size_t i = n; i-- > 0;) {
      rb_gc_register_address(values + i);
    }
  }
}

/*
 * Unregister the addresses of an array of VALUEs with the gc to allow normal
 * garbage collection to occur again.
 */
void nm_unregister_values(VALUE* values, size_t n) {
  if (values) {
    for (size_t i = n; i-- > 0;) {
      rb_gc_unregister_address(values + i);
    }
  }
}

/*
 * call-seq:
 *     dtype -> Symbol
 *
 * Get the data type (dtype) of a matrix, e.g., :byte, :int8, :int16, :int32,
 * :int64, :float32, :float64, :complex64, :complex128, :rational32,
 * :rational64, :rational128, or :object (the last is a Ruby object).
 */
static VALUE nm_dtype(VALUE self) {
  ID dtype = rb_intern(DTYPE_NAMES[NM_DTYPE(self)]);
  return ID2SYM(dtype);
}


/*
 * call-seq:
 *     upcast(first_dtype, second_dtype) -> Symbol
 *
 * Given a binary operation between types t1 and t2, what type will be returned?
 *
 * This is a singleton method on NMatrix, e.g., NMatrix.upcast(:int32, :int64)
 */
static VALUE nm_upcast(VALUE self, VALUE t1, VALUE t2) {

  nm::dtype_t d1    = nm_dtype_from_rbsymbol(t1),
              d2    = nm_dtype_from_rbsymbol(t2);

  return ID2SYM(rb_intern( DTYPE_NAMES[ Upcast[d1][d2] ] ));
}


/*
 * call-seq:
       default_value -> ...
 *
 * Get the default value for the matrix. For dense, this is undefined and will return Qnil. For list, it is user-defined.
 * For yale, it's going to be some variation on zero, but may be Qfalse or Qnil.
 */
static VALUE nm_default_value(VALUE self) {
  switch(NM_STYPE(self)) {
  case nm::YALE_STORE:
    return nm_yale_default_value(self);
  case nm::LIST_STORE:
    return nm_list_default_value(self);
  case nm::DENSE_STORE:
  default:
    return Qnil;
  }
}


/*
 * call-seq:
 *     each_with_indices -> Enumerator
 *
 * Iterate over all entries of any matrix in standard storage order (as with #each), and include the indices.
 */
static VALUE nm_each_with_indices(VALUE nmatrix) {
  volatile VALUE nm = nmatrix;

  switch(NM_STYPE(nm)) {
  case nm::YALE_STORE:
    return nm_yale_each_with_indices(nm);
  case nm::DENSE_STORE:
    return nm_dense_each_with_indices(nm);
  case nm::LIST_STORE:
    return nm_list_each_with_indices(nm, false);
  default:
    rb_raise(nm_eDataTypeError, "Not a proper storage type");
  }
}

/*
 * call-seq:
 *     each_stored_with_indices -> Enumerator
 *
 * Iterate over the stored entries of any matrix. For dense and yale, this iterates over non-zero
 * entries; for list, this iterates over non-default entries. Yields dim+1 values for each entry:
 * i, j, ..., and the entry itself.
 */
static VALUE nm_each_stored_with_indices(VALUE nmatrix) {
  volatile VALUE nm = nmatrix;

  switch(NM_STYPE(nm)) {
  case nm::YALE_STORE:
    return nm_yale_each_stored_with_indices(nm);
  case nm::DENSE_STORE:
    return nm_dense_each_with_indices(nm);
  case nm::LIST_STORE:
    return nm_list_each_with_indices(nm, true);
  default:
    rb_raise(nm_eDataTypeError, "Not a proper storage type");
  }
}


/*
 * call-seq:
 *     each_ordered_stored_with_indices -> Enumerator
 *
 * Very similar to #each_stored_with_indices. The key difference is that it enforces matrix ordering rather
 * than storage ordering, which only matters if your matrix is Yale.
 */
static VALUE nm_each_ordered_stored_with_indices(VALUE nmatrix) {
  volatile VALUE nm = nmatrix;

  switch(NM_STYPE(nm)) {
  case nm::YALE_STORE:
    return nm_yale_each_ordered_stored_with_indices(nm);
  case nm::DENSE_STORE:
    return nm_dense_each_with_indices(nm);
  case nm::LIST_STORE:
    return nm_list_each_with_indices(nm, true);
  default:
    rb_raise(nm_eDataTypeError, "Not a proper storage type");
  }
}


/*
 * Equality operator. Returns a single true or false value indicating whether
 * the matrices are equivalent.
 *
 * For elementwise, use =~ instead.
 *
 * This method will raise an exception if dimensions do not match.
 */
static VALUE nm_eqeq(VALUE left, VALUE right) {
  NMATRIX *l, *r;

  CheckNMatrixType(left);
  CheckNMatrixType(right);

  UnwrapNMatrix(left, l);
  UnwrapNMatrix(right, r);

  if (l->stype != r->stype)
    rb_raise(rb_eNotImpError, "comparison between different matrix stypes not yet implemented");

  bool result = false;

  switch(l->stype) {
  case nm::DENSE_STORE:
    result = nm_dense_storage_eqeq(l->storage, r->storage);
    break;
  case nm::LIST_STORE:
    result = nm_list_storage_eqeq(l->storage, r->storage);
    break;
  case nm::YALE_STORE:
    result = nm_yale_storage_eqeq(l->storage, r->storage);
    break;
  }

  return result ? Qtrue : Qfalse;
}

DEF_ELEMENTWISE_RUBY_ACCESSOR(ADD, add)
DEF_ELEMENTWISE_RUBY_ACCESSOR(SUB, subtract)
DEF_ELEMENTWISE_RUBY_ACCESSOR(MUL, multiply)
DEF_ELEMENTWISE_RUBY_ACCESSOR(DIV, divide)
DEF_ELEMENTWISE_RUBY_ACCESSOR(POW, power)
DEF_ELEMENTWISE_RUBY_ACCESSOR(MOD, mod)
DEF_ELEMENTWISE_RUBY_ACCESSOR(EQEQ, eqeq)
DEF_ELEMENTWISE_RUBY_ACCESSOR(NEQ, neq)
DEF_ELEMENTWISE_RUBY_ACCESSOR(LEQ, leq)
DEF_ELEMENTWISE_RUBY_ACCESSOR(GEQ, geq)
DEF_ELEMENTWISE_RUBY_ACCESSOR(LT, lt)
DEF_ELEMENTWISE_RUBY_ACCESSOR(GT, gt)

/*
 * call-seq:
 *     hermitian? -> Boolean
 *
 * Is this matrix hermitian?
 *
 * Definition: http://en.wikipedia.org/wiki/Hermitian_matrix
 *
 * For non-complex matrices, this function should return the same result as symmetric?.
 */
static VALUE nm_hermitian(VALUE self) {
  return is_symmetric(self, true);
}


/*
 * call-seq:
 *     complex_conjugate -> NMatrix
 *
 * Transform the matrix (in-place) to its complex conjugate. Only works on complex matrices.
 *
 * FIXME: For non-complex matrices, someone needs to implement a non-in-place complex conjugate (which doesn't use a bang).
 * Bang should imply that no copy is being made, even temporarily.
 */
static VALUE nm_complex_conjugate_bang(VALUE self) {
  NMATRIX* m;
  void* elem;
  size_t size, p;

  UnwrapNMatrix(self, m);

  if (m->stype == nm::DENSE_STORE) {

    size = nm_storage_count_max_elements(NM_STORAGE(self));
    elem = NM_STORAGE_DENSE(self)->elements;

  } else if (m->stype == nm::YALE_STORE) {

    size = nm_yale_storage_get_size(NM_STORAGE_YALE(self));
    elem = NM_STORAGE_YALE(self)->a;

  } else {
    rb_raise(rb_eNotImpError, "please cast to yale or dense (complex) first");
  }

  // Walk through and negate the imaginary component
  if (NM_DTYPE(self) == nm::COMPLEX64) {

    for (p = 0; p < size; ++p) {
      reinterpret_cast<nm::Complex64*>(elem)[p].i = -reinterpret_cast<nm::Complex64*>(elem)[p].i;
    }

  } else if (NM_DTYPE(self) == nm::COMPLEX128) {

    for (p = 0; p < size; ++p) {
      reinterpret_cast<nm::Complex128*>(elem)[p].i = -reinterpret_cast<nm::Complex128*>(elem)[p].i;
    }

  } else {
    rb_raise(nm_eDataTypeError, "can only calculate in-place complex conjugate on matrices of type :complex64 or :complex128");
  }

  return self;
}

/*
 * Helper function for creating a matrix. You have to create the storage and pass it in, but you don't
 * need to worry about deleting it.
 */
NMATRIX* nm_create(nm::stype_t stype, STORAGE* storage) {
  NMATRIX* mat = ALLOC(NMATRIX);

  mat->stype   = stype;
  mat->storage = storage;

  return mat;
}

/*
 * @see nm_init
 */
static VALUE nm_init_new_version(int argc, VALUE* argv, VALUE self) {
  VALUE shape_ary, initial_ary, hash;
  //VALUE shape_ary, default_val, capacity, initial_ary, dtype_sym, stype_sym;
  // Mandatory args: shape, dtype, stype
  // FIXME: This is the one line of code standing between Ruby 1.9.2 and 1.9.3.
#ifndef OLD_RB_SCAN_ARGS // Ruby 1.9.3 and higher
  rb_scan_args(argc, argv, "11:", &shape_ary, &initial_ary, &hash); // &stype_sym, &dtype_sym, &default_val, &capacity);
#else // Ruby 1.9.2 and lower
  if (argc == 3)
    rb_scan_args(argc, argv, "12", &shape_ary, &initial_ary, &hash);
  else if (argc == 2) {
    VALUE unknown_arg;
    rb_scan_args(argc, argv, "11", &shape_ary, &unknown_arg);
    if (!NIL_P(unknown_arg) && TYPE(unknown_arg) == T_HASH) {
      hash        = unknown_arg;
      initial_ary = Qnil;
    } else {
      initial_ary  = unknown_arg;
      hash        = Qnil;
    }
  }
#endif

  // Get the shape.
  size_t  dim;
  size_t* shape = interpret_shape(shape_ary, &dim);
  void*   init;
  void*   v = NULL;
  size_t  v_size = 0;

  nm::stype_t stype = nm::DENSE_STORE;
  nm::dtype_t dtype = nm::RUBYOBJ;
  VALUE dtype_sym = Qnil, stype_sym = Qnil, default_val_num = Qnil, capacity_num = Qnil;
  size_t capacity = 0;
  if (!NIL_P(hash)) {
    dtype_sym       = rb_hash_aref(hash, ID2SYM(nm_rb_dtype));
    stype_sym       = rb_hash_aref(hash, ID2SYM(nm_rb_stype));
    capacity_num    = rb_hash_aref(hash, ID2SYM(nm_rb_capacity));
    default_val_num = rb_hash_aref(hash, ID2SYM(nm_rb_default));
  }

  //     stype ||= :dense
  stype = !NIL_P(stype_sym) ? nm_stype_from_rbsymbol(stype_sym) : nm::DENSE_STORE;

  //     dtype ||= h[:dtype] || guess_dtype(initial_ary) || :object
  if (NIL_P(initial_ary) && NIL_P(dtype_sym))
    dtype = nm::RUBYOBJ;
  else if (NIL_P(dtype_sym))
    dtype = nm_dtype_guess(initial_ary);
  else
    dtype = nm_dtype_from_rbsymbol(dtype_sym);

  //   if stype != :dense
  //     if initial_ary.nil?
  //       init = h[:default] || 0
  //     elsif initial_ary.is_a?(Array)
  //       init = initial_ary.size > 1 ? (h[:default] || 0) : initial_ary[0]
  //     else
  //       init = initial_ary # not an array, just a value
  //     end
  //   end
  if (stype != nm::DENSE_STORE) {
    if (!NIL_P(default_val_num))
      init = rubyobj_to_cval(default_val_num, dtype);
    else if (NIL_P(initial_ary))
      init = NULL;
    else if (TYPE(initial_ary) == T_ARRAY)
      init = RARRAY_LEN(initial_ary) == 1 ? rubyobj_to_cval(rb_ary_entry(initial_ary, 0), dtype) : NULL;
    else
      init = rubyobj_to_cval(initial_ary, dtype);
  }

  // capacity = h[:capacity] || 0
  if (stype == nm::YALE_STORE) {
    if (!NIL_P(capacity_num)) capacity = FIX2INT(capacity_num);
  }

  if (!NIL_P(initial_ary)) {
    v = interpret_initial_value(initial_ary, dtype);

    if (TYPE(initial_ary) == T_ARRAY) 	v_size = RARRAY_LEN(initial_ary);
    else                                v_size = 1;
  }

  // :object matrices MUST be initialized.
  else if (stype == nm::DENSE_STORE && dtype == nm::RUBYOBJ) {
    // Pretend [nil] was passed for RUBYOBJ.
    v          = ALLOC(VALUE);
    *(VALUE*)v = Qnil;

    v_size = 1;

  }

	NMATRIX* nmatrix;
  UnwrapNMatrix(self, nmatrix);

  nmatrix->stype = stype;

  switch (stype) {
  	case nm::DENSE_STORE:
  		nmatrix->storage = (STORAGE*)nm_dense_storage_create(dtype, shape, dim, v, v_size);
  		break;

  	case nm::LIST_STORE:
  		nmatrix->storage = (STORAGE*)nm_list_storage_create(dtype, shape, dim, init);
  		break;

  	case nm::YALE_STORE:
  		nmatrix->storage = (STORAGE*)nm_yale_storage_create(dtype, shape, dim, capacity);
  		nm_yale_storage_init((YALE_STORAGE*)(nmatrix->storage), init);
  		break;
  }

  // If we're not creating a dense, and an initial array was provided, use that and multi-slice-set
  // to set the contents of the matrix right now.
  if (stype != nm::DENSE_STORE && v_size > 1) {
    VALUE* slice_argv = ALLOCA_N(VALUE, dim);
    size_t* tmp_shape = ALLOC_N(size_t, dim);
    for (size_t m = 0; m < dim; ++m) {
      slice_argv[m] = ID2SYM(nm_rb_mul); // :* -- full range
      tmp_shape[m]  = shape[m];
    }

    SLICE* slice = get_slice(dim, dim, slice_argv, shape);
    // Create a temporary dense matrix and use it to do a slice assignment on self.
    NMATRIX* tmp          = nm_create(nm::DENSE_STORE, (STORAGE*)nm_dense_storage_create(dtype, tmp_shape, dim, v, v_size));
    volatile VALUE rb_tmp = Data_Wrap_Struct(CLASS_OF(self), nm_mark, nm_delete, tmp);
    if (stype == nm::YALE_STORE)  nm_yale_storage_set(self, slice, rb_tmp);
    else                          nm_list_storage_set(self, slice, rb_tmp);

    free_slice(slice);

    // We need to free v if it's not the same size as tmp -- because tmp will have made a copy instead.
    if (nm_storage_count_max_elements(tmp->storage) != v_size)
      xfree(v);

    // nm_delete(tmp); // This seems to enrage the garbage collector (because rb_tmp is still available). It'd be better if we could force it to free immediately, but no sweat.
  }

  return self;
}

/*
 * call-seq:
 *     new(shape) -> NMatrix
 *     new(shape, initial_value) -> NMatrix
 *     new(shape, initial_array) -> NMatrix
 *     new(shape, initial_value, options) -> NMatrix
 *     new(shape, initial_array, options) -> NMatrix
 *
 * Create a new NMatrix.
 *
 * The only mandatory argument is shape, which may be a positive integer or an array of positive integers.
 *
 * It is recommended that you supply an initialization value or array of values. Without one, Yale and List matrices will
 * be initialized to 0; and dense matrices will be undefined.
 *
 * Additional options may be provided using keyword arguments. The keywords are +:dtype, +:stype+, +:capacity+, and
 * +:default+. Only Yale uses a capacity argument, which is used to reserve the initial size of its storage vectors.
 * List and Yale both accept a default value (which itself defaults to 0). This default is taken from the initial value
 * if such a value is given; it is more likely to be required when an initial array is provided.
 *
 * The storage type, or stype, is used to specify whether we want a +:dense+, +:list+, or +:yale+ matrix; dense is the
 * default.
 *
 * The data type, or dtype, can be one of: :byte, :int8, :int16, :int32, :int64, :float32, :float64, :complex64,
 * :complex128, :rational128, or :object. The constructor will attempt to guess it from the initial value/array/default
 * provided, if any. Otherwise, the default is :object, which stores any type of Ruby object.
 *
 * In addition to the above, there is a legacy constructor from the alpha version. To use that version, you must be
 * providing exactly four arguments. It is now deprecated.
 *
 * There is one additional constructor for advanced users, which takes seven arguments and is only for creating Yale
 * matrices with known IA, JA, and A arrays. This is used primarily internally for IO, e.g., reading Matlab matrices,
 * which are stored in old Yale (not our Yale) format. But be careful; there are no overflow warnings. All of these
 * constructors are defined for power-users. Everyone else should probably resort to the shortcut functions defined in
 * shortcuts.rb.
 */
static VALUE nm_init(int argc, VALUE* argv, VALUE nm) {

  if (argc <= 3) { // Call the new constructor unless all four arguments are given (or the 7-arg version is given)
  	return nm_init_new_version(argc, argv, nm);
  }

  /* First, determine stype (dense by default) */
  nm::stype_t stype;
  size_t  offset = 0;

  if (!SYMBOL_P(argv[0]) && TYPE(argv[0]) != T_STRING) {
    stype = nm::DENSE_STORE;

  } else {
    // 0: String or Symbol
    stype  = interpret_stype(argv[0]);
    offset = 1;
  }

  // If there are 7 arguments and Yale, refer to a different init function with fewer sanity checks.
  if (argc == 7) {
  	if (stype == nm::YALE_STORE) {
			return nm_init_yale_from_old_yale(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], nm);

		} else {
			rb_raise(rb_eArgError, "Expected 2-4 arguments (or 7 for internal Yale creation)");
		}
  }

	// 1: Array or Fixnum
	size_t dim;
  size_t* shape = interpret_shape(argv[offset], &dim);

  // 2-3: dtype
  nm::dtype_t dtype = interpret_dtype(argc-1-offset, argv+offset+1, stype);

  size_t init_cap = 0, init_val_len = 0;
  void* init_val  = NULL;
  if (!SYMBOL_P(argv[1+offset]) || TYPE(argv[1+offset]) == T_ARRAY) {
  	// Initial value provided (could also be initial capacity, if yale).

    if (stype == nm::YALE_STORE && NM_RUBYVAL_IS_NUMERIC(argv[1+offset])) {
      init_cap = FIX2UINT(argv[1+offset]);

    } else {
    	// 4: initial value / dtype
      init_val = interpret_initial_value(argv[1+offset], dtype);

      if (TYPE(argv[1+offset]) == T_ARRAY) 	init_val_len = RARRAY_LEN(argv[1+offset]);
      else                                  init_val_len = 1;
    }

  } else {
  	// DType is RUBYOBJ.

    if (stype == nm::DENSE_STORE) {
    	/*
    	 * No need to initialize dense with any kind of default value unless it's
    	 * an RUBYOBJ matrix.
    	 */
      if (dtype == nm::RUBYOBJ) {
      	// Pretend [nil] was passed for RUBYOBJ.
      	init_val = ALLOC(VALUE);
        *(VALUE*)init_val = Qnil;

        init_val_len = 1;

      } else {
      	init_val = NULL;
      }
    } else if (stype == nm::LIST_STORE) {
    	init_val = ALLOC_N(char, DTYPE_SIZES[dtype]);
      std::memset(init_val, 0, DTYPE_SIZES[dtype]);
    }
  }

  // TODO: Update to allow an array as the initial value.
	NMATRIX* nmatrix;
  UnwrapNMatrix(nm, nmatrix);

  nmatrix->stype = stype;

  switch (stype) {
  	case nm::DENSE_STORE:
  		nmatrix->storage = (STORAGE*)nm_dense_storage_create(dtype, shape, dim, init_val, init_val_len);
  		break;

  	case nm::LIST_STORE:
  		nmatrix->storage = (STORAGE*)nm_list_storage_create(dtype, shape, dim, init_val);
  		break;

  	case nm::YALE_STORE:
  		nmatrix->storage = (STORAGE*)nm_yale_storage_create(dtype, shape, dim, init_cap);
  		nm_yale_storage_init((YALE_STORAGE*)(nmatrix->storage), NULL);
  		break;
  }

  return nm;
}


/*
 * Helper for nm_cast which uses the C types instead of the Ruby objects. Called by nm_cast.
 */
NMATRIX* nm_cast_with_ctype_args(NMATRIX* self, nm::stype_t new_stype, nm::dtype_t new_dtype, void* init_ptr) {
  NMATRIX* lhs = ALLOC(NMATRIX);
  lhs->stype   = new_stype;

  // Copy the storage
  CAST_TABLE(cast_copy);
  lhs->storage = cast_copy[lhs->stype][self->stype](self->storage, new_dtype, init_ptr);

  return lhs;
}


/*
 * call-seq:
 *     cast_full(stype) -> NMatrix
 *     cast_full(stype, dtype, sparse_basis) -> NMatrix
 *
 * Copy constructor for changing dtypes and stypes.
 */
VALUE nm_cast(VALUE self, VALUE new_stype_symbol, VALUE new_dtype_symbol, VALUE init) {
  nm::dtype_t new_dtype = nm_dtype_from_rbsymbol(new_dtype_symbol);
  nm::stype_t new_stype = nm_stype_from_rbsymbol(new_stype_symbol);

  CheckNMatrixType(self);
  NMATRIX *rhs;

  UnwrapNMatrix( self, rhs );

  void* init_ptr = ALLOCA_N(char, DTYPE_SIZES[new_dtype]);
  rubyval_to_cval(init, new_dtype, init_ptr);

  return Data_Wrap_Struct(CLASS_OF(self), nm_mark, nm_delete, nm_cast_with_ctype_args(rhs, new_stype, new_dtype, init_ptr));
}

/*
 * Copy constructor for transposing.
 */
static VALUE nm_init_transposed(VALUE self) {
  static STORAGE* (*storage_copy_transposed[nm::NUM_STYPES])(const STORAGE* rhs_base) = {
    nm_dense_storage_copy_transposed,
    nm_list_storage_copy_transposed,
    nm_yale_storage_copy_transposed
  };

  NMATRIX* lhs = nm_create( NM_STYPE(self),
                            storage_copy_transposed[NM_STYPE(self)]( NM_STORAGE(self) )
                          );

  return Data_Wrap_Struct(CLASS_OF(self), nm_mark, nm_delete, lhs);
}

/*
 * Copy constructor for no change of dtype or stype (used for #initialize_copy hook).
 */
static VALUE nm_init_copy(VALUE copy, VALUE original) {
  NMATRIX *lhs, *rhs;

  CheckNMatrixType(original);

  if (copy == original) return copy;

  UnwrapNMatrix( original, rhs );
  UnwrapNMatrix( copy,     lhs );

  lhs->stype = rhs->stype;

  // Copy the storage
  CAST_TABLE(ttable);
  lhs->storage = ttable[lhs->stype][rhs->stype](rhs->storage, rhs->storage->dtype, NULL);

  return copy;
}

/*
 * Get major, minor, and release components of NMatrix::VERSION. Store in function parameters.
 */
static void get_version_info(uint16_t& major, uint16_t& minor, uint16_t& release) {
  // Get VERSION and split it on periods. Result is an Array.
  VALUE version = rb_funcall(rb_const_get(cNMatrix, rb_intern("VERSION")), rb_intern("split"), 1, rb_str_new_cstr("."));
  VALUE* ary    = RARRAY_PTR(version); // major, minor, and release

  // Convert each to an integer
  VALUE  maj    = rb_funcall(ary[0], rb_intern("to_i"), 0);
  VALUE  min    = rb_funcall(ary[1], rb_intern("to_i"), 0);
  VALUE  rel    = rb_funcall(ary[2], rb_intern("to_i"), 0);

  major   = static_cast<uint16_t>(nm::RubyObject(maj));
  minor   = static_cast<uint16_t>(nm::RubyObject(min));
  release = static_cast<uint16_t>(nm::RubyObject(rel));
}


/*
 * Interpret the NMatrix::write symmetry argument (which should be nil or a symbol). Return a symm_t (enum).
 */
static nm::symm_t interpret_symm(VALUE symm) {
  if (symm == Qnil) return nm::NONSYMM;

  ID rb_symm = rb_intern("symmetric"),
     rb_skew = rb_intern("skew"),
     rb_herm = rb_intern("hermitian");
     // nm_rb_upper, nm_rb_lower already set

  ID symm_id = rb_to_id(symm);

  if (symm_id == rb_symm)            return nm::SYMM;
  else if (symm_id == rb_skew)       return nm::SKEW;
  else if (symm_id == rb_herm)       return nm::HERM;
  else if (symm_id == nm_rb_upper)   return nm::UPPER;
  else if (symm_id == nm_rb_lower)   return nm::LOWER;
  else                            rb_raise(rb_eArgError, "unrecognized symmetry argument");

  return nm::NONSYMM;
}



void read_padded_shape(std::ifstream& f, size_t dim, size_t* shape) {
  size_t bytes_read = 0;

  // Read shape
  for (size_t i = 0; i < dim; ++i) {
    size_t s;
    f.read(reinterpret_cast<char*>(&s), sizeof(size_t));
    shape[i] = s;

    std::cerr << "shape(" << i << ")=" << s << std::endl;

    bytes_read += sizeof(size_t);
  }

  // Ignore padding
  f.ignore(bytes_read % 8);
}


void write_padded_shape(std::ofstream& f, size_t dim, size_t* shape) {
  size_t bytes_written = 0;

  // Write shape
  for (size_t i = 0; i < dim; ++i) {
    size_t s = shape[i];
    std::cerr << "writing shape(" << i << ")=" << shape[i] << std::endl;
    f.write(reinterpret_cast<const char*>(&s), sizeof(size_t));

    bytes_written += sizeof(size_t);
  }

  // Pad with zeros
  size_t zero = 0;
  while (bytes_written % 8) {
    f.write(reinterpret_cast<const char*>(&zero), sizeof(size_t));

    bytes_written += sizeof(IType);
  }
}


void read_padded_yale_elements(std::ifstream& f, YALE_STORAGE* storage, size_t length, nm::symm_t symm, nm::dtype_t dtype) {
  NAMED_DTYPE_TEMPLATE_TABLE_NO_ROBJ(ttable, nm::read_padded_yale_elements, void, std::ifstream&, YALE_STORAGE*, size_t, nm::symm_t)

  ttable[dtype](f, storage, length, symm);
}


void write_padded_yale_elements(std::ofstream& f, YALE_STORAGE* storage, size_t length, nm::symm_t symm, nm::dtype_t dtype) {
  NAMED_DTYPE_TEMPLATE_TABLE_NO_ROBJ(ttable, nm::write_padded_yale_elements, void, std::ofstream& f, YALE_STORAGE*, size_t, nm::symm_t)

  ttable[dtype](f, storage, length, symm);
}


void read_padded_dense_elements(std::ifstream& f, DENSE_STORAGE* storage, nm::symm_t symm, nm::dtype_t dtype) {
  NAMED_DTYPE_TEMPLATE_TABLE_NO_ROBJ(ttable, nm::read_padded_dense_elements, void, std::ifstream&, DENSE_STORAGE*, nm::symm_t)

  ttable[dtype](f, storage, symm);
}


void write_padded_dense_elements(std::ofstream& f, DENSE_STORAGE* storage, nm::symm_t symm, nm::dtype_t dtype) {
  NAMED_DTYPE_TEMPLATE_TABLE_NO_ROBJ(ttable, nm::write_padded_dense_elements, void, std::ofstream& f, DENSE_STORAGE*, nm::symm_t)

  ttable[dtype](f, storage, symm);
}


/*
 * Helper function to get exceptions in the module Errno (e.g., ENOENT). Example:
 *
 *     rb_raise(rb_get_errno_exc("ENOENT"), RSTRING_PTR(filename));
 */
static VALUE rb_get_errno_exc(const char* which) {
  return rb_const_get(rb_const_get(rb_cObject, rb_intern("Errno")), rb_intern(which));
}



/*
 * Binary file writer for NMatrix standard format. file should be a path, which we aren't going to
 * check very carefully (in other words, this function should generally be called from a Ruby
 * helper method). Function also takes a symmetry argument, which allows us to specify that we only want to
 * save the upper triangular portion of the matrix (or if the matrix is a lower triangular matrix, only
 * the lower triangular portion). nil means regular storage.
 */
static VALUE nm_write(int argc, VALUE* argv, VALUE self) {
  using std::ofstream;

  std::cerr << "a" << std::endl;

  if (argc < 1 || argc > 2) {
    rb_raise(rb_eArgError, "Expected one or two arguments");
  }
  VALUE file = argv[0],
        symm = argc == 1 ? Qnil : argv[1];

  NMATRIX* nmatrix;
  UnwrapNMatrix( self, nmatrix );

  nm::symm_t symm_ = interpret_symm(symm);

  if (nmatrix->storage->dtype == nm::RUBYOBJ) {
    rb_raise(rb_eNotImpError, "Ruby Object writing is not implemented yet");
  }

  // Get the dtype, stype, itype, and symm and ensure they're the correct number of bytes.
  uint8_t st = static_cast<uint8_t>(nmatrix->stype),
          dt = static_cast<uint8_t>(nmatrix->storage->dtype),
          sm = static_cast<uint8_t>(symm_);
  uint16_t dim = nmatrix->storage->dim;

  //FIXME: Cast the matrix to the smallest possible index type. Write that in the place of IType.

  // Check arguments before starting to write.
  if (nmatrix->stype == nm::LIST_STORE) rb_raise(nm_eStorageTypeError, "cannot save list matrix; cast to yale or dense first");
  if (symm_ != nm::NONSYMM) {
    if (dim != 2) rb_raise(rb_eArgError, "symmetry/triangularity not defined for a non-2D matrix");
    if (nmatrix->storage->shape[0] != nmatrix->storage->shape[1])
      rb_raise(rb_eArgError, "symmetry/triangularity not defined for a non-square matrix");
    if (symm_ == nm::HERM &&
          dt != static_cast<uint8_t>(nm::COMPLEX64) && dt != static_cast<uint8_t>(nm::COMPLEX128) && dt != static_cast<uint8_t>(nm::RUBYOBJ))
      rb_raise(rb_eArgError, "cannot save a non-complex matrix as hermitian");
  }

  ofstream f(RSTRING_PTR(file), std::ios::out | std::ios::binary);

  // Get the NMatrix version information.
  uint16_t major, minor, release, null16 = 0;
  get_version_info(major, minor, release);

  std::cerr << "b" << std::endl;


  // WRITE FIRST 64-BIT BLOCK
  f.write(reinterpret_cast<const char*>(&major),   sizeof(uint16_t));
  f.write(reinterpret_cast<const char*>(&minor),   sizeof(uint16_t));
  f.write(reinterpret_cast<const char*>(&release), sizeof(uint16_t));
  f.write(reinterpret_cast<const char*>(&null16),  sizeof(uint16_t));

  uint8_t ZERO = 0;
  // WRITE SECOND 64-BIT BLOCK
  f.write(reinterpret_cast<const char*>(&dt), sizeof(uint8_t));
  f.write(reinterpret_cast<const char*>(&st), sizeof(uint8_t));
  f.write(reinterpret_cast<const char*>(&ZERO),sizeof(uint8_t));
  f.write(reinterpret_cast<const char*>(&sm), sizeof(uint8_t));
  f.write(reinterpret_cast<const char*>(&null16), sizeof(uint16_t));
  f.write(reinterpret_cast<const char*>(&dim), sizeof(uint16_t));

  std::cerr << "c" << std::endl;
  // Write shape (in 64-bit blocks)
  write_padded_shape(f, nmatrix->storage->dim, nmatrix->storage->shape);

  std::cerr << "d" << std::endl;
  if (nmatrix->stype == nm::DENSE_STORE) {
    std::cerr << "e" << std::endl;
    write_padded_dense_elements(f, reinterpret_cast<DENSE_STORAGE*>(nmatrix->storage), symm_, nmatrix->storage->dtype);
    std::cerr << "f" << std::endl;
  } else if (nmatrix->stype == nm::YALE_STORE) {
    YALE_STORAGE* s = reinterpret_cast<YALE_STORAGE*>(nmatrix->storage);
    uint32_t ndnz   = s->ndnz,
             length = nm_yale_storage_get_size(s);
    std::cerr << "g" << std::endl;
    f.write(reinterpret_cast<const char*>(&ndnz),   sizeof(uint32_t));
    f.write(reinterpret_cast<const char*>(&length), sizeof(uint32_t));

    std::cerr << "h" << std::endl;
    write_padded_yale_elements(f, s, length, symm_, s->dtype);
    std::cerr << "i" << std::endl;
  }

  f.close();

  return Qtrue;
}


/*
 * Binary file reader for NMatrix standard format. file should be a path, which we aren't going to
 * check very carefully (in other words, this function should generally be called from a Ruby
 * helper method).
 *
 * Note that currently, this function will by default refuse to read files that are newer than
 * your version of NMatrix. To force an override, set the second argument to anything other than nil.
 *
 * Returns an NMatrix Ruby object.
 */
static VALUE nm_read(int argc, VALUE* argv, VALUE self) {
  using std::ifstream;

  VALUE file, force_;

  // Read the arguments
  rb_scan_args(argc, argv, "11", &file, &force_);
  bool force   = (force_ != Qnil && force_ != Qfalse);


  if (!RB_FILE_EXISTS(file)) { // FIXME: Errno::ENOENT
    rb_raise(rb_get_errno_exc("ENOENT"), "%s", RSTRING_PTR(file));
  }

  // Open a file stream
  ifstream f(RSTRING_PTR(file), std::ios::in | std::ios::binary);

  uint16_t major, minor, release;
  get_version_info(major, minor, release); // compare to NMatrix version

  uint16_t fmajor, fminor, frelease, null16;

  // READ FIRST 64-BIT BLOCK
  f.read(reinterpret_cast<char*>(&fmajor),   sizeof(uint16_t));
  f.read(reinterpret_cast<char*>(&fminor),   sizeof(uint16_t));
  f.read(reinterpret_cast<char*>(&frelease), sizeof(uint16_t));
  f.read(reinterpret_cast<char*>(&null16),   sizeof(uint16_t));

  int ver  = major * 10000 + minor * 100 + release,
      fver = fmajor * 10000 + fminor * 100 + release;
  if (fver > ver && force == false) {
    rb_raise(rb_eIOError, "File was created in newer version of NMatrix than current");
  }
  if (null16 != 0) fprintf(stderr, "Warning: Expected zero padding was not zero\n");

  uint8_t dt, st, it, sm;
  uint16_t dim;

  // READ SECOND 64-BIT BLOCK
  f.read(reinterpret_cast<char*>(&dt), sizeof(uint8_t));
  f.read(reinterpret_cast<char*>(&st), sizeof(uint8_t));
  f.read(reinterpret_cast<char*>(&it), sizeof(uint8_t)); // FIXME: should tell how few bytes indices are stored as
  f.read(reinterpret_cast<char*>(&sm), sizeof(uint8_t));
  f.read(reinterpret_cast<char*>(&null16), sizeof(uint16_t));
  f.read(reinterpret_cast<char*>(&dim), sizeof(uint16_t));
  std::cerr << "dim = " << dim << std::endl;

  if (null16 != 0) fprintf(stderr, "Warning: Expected zero padding was not zero\n");
  nm::stype_t stype = static_cast<nm::stype_t>(st);
  nm::dtype_t dtype = static_cast<nm::dtype_t>(dt);
  nm::symm_t  symm  = static_cast<nm::symm_t>(sm);
  //nm::itype_t itype = static_cast<nm::itype_t>(it);

  // READ NEXT FEW 64-BIT BLOCKS
  size_t* shape = ALLOC_N(size_t, dim);
  read_padded_shape(f, dim, shape);

  STORAGE* s;
  if (stype == nm::DENSE_STORE) {
    s = nm_dense_storage_create(dtype, shape, dim, NULL, 0);

    read_padded_dense_elements(f, reinterpret_cast<DENSE_STORAGE*>(s), symm, dtype);

  } else if (stype == nm::YALE_STORE) {
    uint32_t ndnz, length;

    // READ YALE-SPECIFIC 64-BIT BLOCK
    f.read(reinterpret_cast<char*>(&ndnz),     sizeof(uint32_t));
    f.read(reinterpret_cast<char*>(&length),   sizeof(uint32_t));

    s = nm_yale_storage_create(dtype, shape, dim, length); // set length as init capacity

    read_padded_yale_elements(f, reinterpret_cast<YALE_STORAGE*>(s), length, symm, dtype);
  } else {
    rb_raise(nm_eStorageTypeError, "please convert to yale or dense before saving");
  }

  NMATRIX* nm = nm_create(stype, s);

  // Return the appropriate matrix object (Ruby VALUE)
  // FIXME: This should probably return CLASS_OF(self) instead of cNMatrix, but I don't know how that works for
  // FIXME: class methods.
  switch(stype) {
  case nm::DENSE_STORE:
  case nm::YALE_STORE:
    return Data_Wrap_Struct(cNMatrix, nm_mark, nm_delete, nm);
  default: // this case never occurs (due to earlier rb_raise)
    return Qnil;
  }

}



/*
 * Create a new NMatrix helper for handling internal ia, ja, and a arguments.
 *
 * This constructor is only called by Ruby code, so we can skip most of the
 * checks.
 */
static VALUE nm_init_yale_from_old_yale(VALUE shape, VALUE dtype, VALUE ia, VALUE ja, VALUE a, VALUE from_dtype, VALUE nm) {
  size_t dim     = 2;
  size_t* shape_  = interpret_shape(shape, &dim);
  nm::dtype_t dtype_  = nm_dtype_from_rbsymbol(dtype);
  char *ia_       = RSTRING_PTR(ia),
       *ja_       = RSTRING_PTR(ja),
       *a_        = RSTRING_PTR(a);
  nm::dtype_t from_dtype_ = nm_dtype_from_rbsymbol(from_dtype);
  NMATRIX* nmatrix;

  UnwrapNMatrix( nm, nmatrix );

  nmatrix->stype   = nm::YALE_STORE;
  nmatrix->storage = (STORAGE*)nm_yale_storage_create_from_old_yale(dtype_, shape_, ia_, ja_, a_, from_dtype_);

  return nm;
}

/*
 * Check to determine whether matrix is a reference to another matrix.
 */
static VALUE nm_is_ref(VALUE self) {
  if (NM_SRC(self) == NM_STORAGE(self)) return Qfalse;
  else                                  return Qtrue;
}

/*
 * call-seq:
 *     slice -> ...
 *
 * Access the contents of an NMatrix at given coordinates, using copying.
 *
 *     n.slice(3,3)  # => 5.0
 *     n.slice(0..1,0..1) #=> matrix [2,2]
 *
 */
static VALUE nm_mget(int argc, VALUE* argv, VALUE self) {
  static void* (*ttable[nm::NUM_STYPES])(const STORAGE*, SLICE*) = {
    nm_dense_storage_get,
    nm_list_storage_get,
    nm_yale_storage_get
  };
  return nm_xslice(argc, argv, ttable[NM_STYPE(self)], nm_delete, self);
}

/*
 * call-seq:
 *     matrix[indices] -> ...
 *
 * Access the contents of an NMatrix at given coordinates by reference.
 *
 *     n[3,3]  # => 5.0
 *     n[0..1,0..1] #=> matrix [2,2]
 *
 */
static VALUE nm_mref(int argc, VALUE* argv, VALUE self) {
  static void* (*ttable[nm::NUM_STYPES])(const STORAGE*, SLICE*) = {
    nm_dense_storage_ref,
    nm_list_storage_ref,
    nm_yale_storage_ref
  };
  return nm_xslice(argc, argv, ttable[NM_STYPE(self)], nm_delete_ref, self);
}

/*
 * Modify the contents of an NMatrix in the given cell
 *
 *     n[3,3] = 5.0
 *
 * Also returns the new contents, so you can chain:
 *
 *     n[3,3] = n[2,3] = 5.0
 */
static VALUE nm_mset(int argc, VALUE* argv, VALUE self) {
  size_t dim = NM_DIM(self); // last arg is the value

  if ((size_t)(argc) > NM_DIM(self)+1) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for %u)", argc, effective_dim(NM_STORAGE(self))+1);
  } else {
    SLICE* slice = get_slice(dim, argc-1, argv, NM_STORAGE(self)->shape);

    static void (*ttable[nm::NUM_STYPES])(VALUE, SLICE*, VALUE) = {
      nm_dense_storage_set,
      nm_list_storage_set,
      nm_yale_storage_set
    };

    ttable[NM_STYPE(self)](self, slice, argv[argc-1]);

    free_slice(slice);

    return argv[argc-1];
  }
  return Qnil;
}

/*
 * Matrix multiply (dot product): against another matrix or a vector.
 *
 * For elementwise, use * instead.
 *
 * The two matrices must be of the same stype (for now). If dtype differs, an upcast will occur.
 */
static VALUE nm_multiply(VALUE left_v, VALUE right_v) {
  NMATRIX *left, *right;

  UnwrapNMatrix( left_v, left );

  if (NM_RUBYVAL_IS_NUMERIC(right_v))
    return matrix_multiply_scalar(left, right_v);

  else if (TYPE(right_v) == T_ARRAY)
    rb_raise(rb_eNotImpError, "please convert array to nx1 or 1xn NMatrix first");

  else { // both are matrices (probably)
    CheckNMatrixType(right_v);
    UnwrapNMatrix( right_v, right );

    if (left->storage->shape[1] != right->storage->shape[0])
      rb_raise(rb_eArgError, "incompatible dimensions");

    if (left->stype != right->stype)
      rb_raise(rb_eNotImpError, "matrices must have same stype");

    return matrix_multiply(left, right);

  }

  return Qnil;
}


/*
 * call-seq:
 *     dim -> Integer
 *
 * Get the number of dimensions of a matrix.
 *
 * In other words, if you set your matrix to be 3x4, the dim is 2. If the
 * matrix was initialized as 3x4x3, the dim is 3.
 *
 * Use #effective_dim to get the dimension of an NMatrix which acts as a vector (e.g., a column or row).
 */
static VALUE nm_dim(VALUE self) {
  return INT2FIX(NM_STORAGE(self)->dim);
}

/*
 * call-seq:
 *     shape -> Array
 *
 * Get the shape (dimensions) of a matrix.
 */
static VALUE nm_shape(VALUE self) {
  STORAGE* s   = NM_STORAGE(self);

  // Copy elements into a VALUE array and then use those to create a Ruby array with rb_ary_new4.
  VALUE* shape = ALLOCA_N(VALUE, s->dim);
  for (size_t index = 0; index < s->dim; ++index)
    shape[index] = INT2FIX(s->shape[index]);

  return rb_ary_new4(s->dim, shape);
}


/*
 * call-seq:
 *     offset -> Array
 *
 * Get the offset (slice position) of a matrix. Typically all zeros, unless you have a reference slice.
 */
static VALUE nm_offset(VALUE self) {
  STORAGE* s   = NM_STORAGE(self);

  // Copy elements into a VALUE array and then use those to create a Ruby array with rb_ary_new4.
  VALUE* offset = ALLOCA_N(VALUE, s->dim);
  for (size_t index = 0; index < s->dim; ++index)
    offset[index] = INT2FIX(s->offset[index]);

  return rb_ary_new4(s->dim, offset);
}


/*
 * call-seq:
 *     supershape -> Array
 *
 * Get the shape of a slice's parent.
 */
static VALUE nm_supershape(VALUE self) {

  STORAGE* s   = NM_STORAGE(self);
  if (s->src == s) return nm_shape(self); // easy case (not a slice)
  else s = s->src;

  VALUE* shape = ALLOCA_N(VALUE, s->dim);
  for (size_t index = 0; index < s->dim; ++index)
    shape[index] = INT2FIX(s->shape[index]);

  return rb_ary_new4(s->dim, shape);
}

/*
 * call-seq:
 *     stype -> Symbol
 *
 * Get the storage type (stype) of a matrix, e.g., :yale, :dense, or :list.
 */
static VALUE nm_stype(VALUE self) {
  ID stype = rb_intern(STYPE_NAMES[NM_STYPE(self)]);
  return ID2SYM(stype);
}

/*
 * call-seq:
 *     symmetric? -> Boolean
 *
 * Is this matrix symmetric?
 */
static VALUE nm_symmetric(VALUE self) {
  return is_symmetric(self, false);
}


/*
 * Gets the dimension of a matrix which might be a vector (have one or more shape components of size 1).
 */
static size_t effective_dim(STORAGE* s) {
  size_t d = 0;
  for (size_t i = 0; i < s->dim; ++i) {
    if (s->shape[i] != 1) d++;
  }
  return d;
}


/*
 * call-seq:
 *     effective_dim -> Fixnum
 *
 * Returns the number of dimensions that don't have length 1. Guaranteed to be less than or equal to #dim.
 */
static VALUE nm_effective_dim(VALUE self) {
  return INT2FIX(effective_dim(NM_STORAGE(self)));
}


/*
 * Get a slice of an NMatrix.
 */
static VALUE nm_xslice(int argc, VALUE* argv, void* (*slice_func)(const STORAGE*, SLICE*), void (*delete_func)(NMATRIX*), VALUE self) {
  VALUE result = Qnil;
  STORAGE* s = NM_STORAGE(self);

  if (NM_DIM(self) < (size_t)(argc)) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for %u)", argc, effective_dim(s));
  } else {
    SLICE* slice = get_slice(NM_DIM(self), argc, argv, s->shape);

    if (slice->single) {
      static void* (*ttable[nm::NUM_STYPES])(const STORAGE*, SLICE*) = {
        nm_dense_storage_ref,
        nm_list_storage_ref,
        nm_yale_storage_ref
      };

      if (NM_DTYPE(self) == nm::RUBYOBJ)  result = *reinterpret_cast<VALUE*>( ttable[NM_STYPE(self)](s, slice) );
      else                                result = rubyobj_from_cval( ttable[NM_STYPE(self)](s, slice), NM_DTYPE(self) ).rval;

    } else {

      NMATRIX* mat  = ALLOC(NMATRIX);
      mat->stype    = NM_STYPE(self);
      mat->storage  = (STORAGE*)((*slice_func)( s, slice ));

      result        = Data_Wrap_Struct(CLASS_OF(self), nm_mark, delete_func, mat);
    }

    free_slice(slice);
  }

  return result;
}

//////////////////////
// Helper Functions //
//////////////////////

static VALUE elementwise_op(nm::ewop_t op, VALUE left_val, VALUE right_val) {

	NMATRIX* left;
	NMATRIX* result;

	CheckNMatrixType(left_val);
	UnwrapNMatrix(left_val, left);

  if (TYPE(right_val) != T_DATA || (RDATA(right_val)->dfree != (RUBY_DATA_FUNC)nm_delete && RDATA(right_val)->dfree != (RUBY_DATA_FUNC)nm_delete_ref)) {
    // This is a matrix-scalar element-wise operation.
    std::string sym;
    switch(left->stype) {
    case nm::DENSE_STORE:
      sym = "__dense_scalar_" + nm::EWOP_NAMES[op] + "__";
      break;
    case nm::YALE_STORE:
      sym = "__yale_scalar_" + nm::EWOP_NAMES[op] + "__";
      break;
    case nm::LIST_STORE:
      sym = "__list_scalar_" + nm::EWOP_NAMES[op] + "__";
      break;
    default:
      rb_raise(rb_eNotImpError, "unknown storage type requested scalar element-wise operation");
    }
    return rb_funcall(left_val, rb_intern(sym.c_str()), 1, right_val);

  } else {

    // Check that the left- and right-hand sides have the same dimensionality.
    if (NM_DIM(left_val) != NM_DIM(right_val)) {
      rb_raise(rb_eArgError, "The left- and right-hand sides of the operation must have the same dimensionality.");
    }

    // Check that the left- and right-hand sides have the same shape.
    if (memcmp(&NM_SHAPE(left_val, 0), &NM_SHAPE(right_val, 0), sizeof(size_t) * NM_DIM(left_val)) != 0) {
      rb_raise(rb_eArgError, "The left- and right-hand sides of the operation must have the same shape.");
    }

    NMATRIX* right;
    UnwrapNMatrix(right_val, right);

    if (left->stype == right->stype) {
      std::string sym;

      switch(left->stype) {
      case nm::DENSE_STORE:
        sym = "__dense_elementwise_" + nm::EWOP_NAMES[op] + "__";
        break;
      case nm::YALE_STORE:
        sym = "__yale_elementwise_" + nm::EWOP_NAMES[op] + "__";
        break;
      case nm::LIST_STORE:
        sym = "__list_elementwise_" + nm::EWOP_NAMES[op] + "__";
        break;
      default:
        rb_raise(rb_eNotImpError, "unknown storage type requested element-wise operation");
      }
      return rb_funcall(left_val, rb_intern(sym.c_str()), 1, right_val);

    } else {
      rb_raise(rb_eArgError, "Element-wise operations are not currently supported between matrices with differing stypes.");
    }
  }

	return Data_Wrap_Struct(CLASS_OF(left_val), nm_mark, nm_delete, result);
}

/*
 * Check to determine whether matrix is a reference to another matrix.
 */
bool is_ref(const NMATRIX* matrix) {
  return matrix->storage->src != matrix->storage;
}

/*
 * Helper function for nm_symmetric and nm_hermitian.
 */
static VALUE is_symmetric(VALUE self, bool hermitian) {
  NMATRIX* m;
  UnwrapNMatrix(self, m);

  if (m->storage->shape[0] == m->storage->shape[1] and m->storage->dim == 2) {
		if (NM_STYPE(self) == nm::DENSE_STORE) {
      if (hermitian) {
        nm_dense_storage_is_hermitian((DENSE_STORAGE*)(m->storage), m->storage->shape[0]);

      } else {
      	nm_dense_storage_is_symmetric((DENSE_STORAGE*)(m->storage), m->storage->shape[0]);
      }

    } else {
      // TODO: Implement, at the very least, yale_is_symmetric. Model it after yale/transp.template.c.
      rb_raise(rb_eNotImpError, "symmetric? and hermitian? only implemented for dense currently");
    }

  }

  return Qfalse;
}

///////////////////////
// Utility Functions //
///////////////////////

/*
 * Guess the dtype given a Ruby VALUE and return it as a symbol.
 *
 * Not to be confused with nm_dtype_guess, which returns an nm::dtype_t. (This calls that.)
 */
static VALUE nm_guess_dtype(VALUE self, VALUE v) {
  return ID2SYM(rb_intern(DTYPE_NAMES[nm_dtype_guess(v)]));
}

/*
 * Get the minimum allowable dtype for a Ruby VALUE and return it as a symbol.
 */
static VALUE nm_min_dtype(VALUE self, VALUE v) {
  return ID2SYM(rb_intern(DTYPE_NAMES[nm_dtype_min(v)]));
}

/*
 * Helper for nm_dtype_min(), handling integers.
 */
nm::dtype_t nm_dtype_min_fixnum(int64_t v) {
  if (v >= 0 && v <= UCHAR_MAX) return nm::BYTE;
  else {
    v = std::abs(v);
    if (v <= CHAR_MAX) return nm::INT8;
    else if (v <= SHRT_MAX) return nm::INT16;
    else if (v <= INT_MAX) return nm::INT32;
    else return nm::INT64;
  }
}

/*
 * Helper for nm_dtype_min(), handling rationals.
 */
nm::dtype_t nm_dtype_min_rational(VALUE vv) {
  nm::Rational128* v = ALLOCA_N(nm::Rational128, 1);
  rubyval_to_cval(vv, nm::RATIONAL128, v);

  int64_t i = std::max(std::abs(v->n), v->d);
  if (i <= SHRT_MAX) return nm::INT16;
  else if (i <= INT_MAX) return nm::INT32;
  else return nm::INT64;
}

/*
 * Return the minimum dtype required to store a given value.
 *
 * This is kind of arbitrary. For Float, it always returns :float32 for example, since in some cases neither :float64
 * not :float32 are sufficient.
 *
 * This function is used in upcasting for scalar math. We want to ensure that :int8 + 1 does not return an :int64, basically.
 *
 * FIXME: Eventually, this function should actually look at the value stored in Fixnums (for example), so that it knows
 * whether to return :int64 or :int32.
 */
nm::dtype_t nm_dtype_min(VALUE v) {

  switch(TYPE(v)) {
  case T_FIXNUM:
    return nm_dtype_min_fixnum(FIX2LONG(v));
  case T_BIGNUM:
    return nm::INT64;
  case T_FLOAT:
    return nm::FLOAT32;
  case T_COMPLEX:
    return nm::COMPLEX64;
  case T_RATIONAL:
    return nm_dtype_min_rational(v);
  case T_STRING:
    return RSTRING_LEN(v) == 1 ? nm::BYTE : nm::RUBYOBJ;
  case T_TRUE:
  case T_FALSE:
  case T_NIL:
  default:
    return nm::RUBYOBJ;
  }
}


/*
 * Guess the data type given a value.
 *
 * TODO: Probably needs some work for Bignum.
 */
nm::dtype_t nm_dtype_guess(VALUE v) {
  switch(TYPE(v)) {
  case T_TRUE:
  case T_FALSE:
  case T_NIL:
    return nm::RUBYOBJ;
  case T_STRING:
    return RSTRING_LEN(v) == 1 ? nm::BYTE : nm::RUBYOBJ;

#if SIZEOF_INT == 8
  case T_FIXNUM:
    return nm::INT64;

  case T_RATIONAL:
    return nm::RATIONAL128;

#else
# if SIZEOF_INT == 4
  case T_FIXNUM:
    return nm::INT32;

  case T_RATIONAL:
    return nm::RATIONAL64;

#else
  case T_FIXNUM:
    return nm::INT16;

  case T_RATIONAL:
    return nm::RATIONAL32;
# endif
#endif

  case T_BIGNUM:
    return nm::INT64;

#if SIZEOF_FLOAT == 4
  case T_COMPLEX:
    return nm::COMPLEX128;

  case T_FLOAT:
    return nm::FLOAT64;

#else
# if SIZEOF_FLOAT == 2
  case T_COMPLEX:
    return nm::COMPLEX64;

  case T_FLOAT:
    return nm::FLOAT32;
# endif
#endif

  case T_ARRAY:
  	/*
  	 * May be passed for dense -- for now, just look at the first element.
     *
  	 * TODO: Look at entire array for most specific type.
  	 */

    return nm_dtype_guess(RARRAY_PTR(v)[0]);

  default:
    RB_P(v);
    rb_raise(rb_eArgError, "Unable to guess a data type from provided parameters; data type must be specified manually.");
  }
}



/*
 * Allocate and return a SLICE object, which will contain the appropriate coordinate and length information for
 * accessing some part of a matrix.
 */
static SLICE* get_slice(size_t dim, int argc, VALUE* arg, size_t* shape) {
  VALUE beg, end;
  int excl;

  SLICE* slice = alloc_slice(dim);
  slice->single = true;

  // r is the shape position; t is the slice position. They may differ when we're dealing with a
  // matrix where the effective dimension is less than the dimension (e.g., a vector).
  for (size_t r = 0, t = 0; r < dim; ++r) {
    VALUE v = t == argc ? Qnil : arg[t];

    // if the current shape indicates a vector and fewer args were supplied than necessary, just use 0
    if (argc - t + r < dim && shape[r] == 1) {
      slice->coords[r]  = 0;
      slice->lengths[r] = 1;

    } else if (FIXNUM_P(v)) { // this used CLASS_OF before, which is inefficient for fixnum

      slice->coords[r]  = FIX2UINT(v);
      slice->lengths[r] = 1;
      t++;

    } else if (SYMBOL_P(v) && rb_to_id(v) == nm_rb_mul) { // :* means the whole possible range

      slice->coords[r]  = 0;
      slice->lengths[r] = shape[r];
      slice->single     = false;

    } else if (TYPE(arg[t]) == T_HASH) { // 3:5 notation (inclusive)
      VALUE begin_end   = rb_funcall(v, rb_intern("shift"), 0); // rb_hash_shift
      slice->coords[r]  = FIX2UINT(rb_ary_entry(begin_end, 0));
      slice->lengths[r] = FIX2UINT(rb_ary_entry(begin_end, 1)) - slice->coords[r];

      if (RHASH_EMPTY_P(v)) t++; // go on to the next

      slice->single = false;

    } else if (CLASS_OF(v) == rb_cRange) {
      rb_range_values(arg[t], &beg, &end, &excl);
      slice->coords[r]  = FIX2UINT(beg);
      // Exclude last element for a...b range
      slice->lengths[r] = FIX2UINT(end) - slice->coords[r] + (excl ? 0 : 1);

      slice->single     = false;

      t++;

    } else {
      rb_raise(rb_eArgError, "expected Fixnum, Range, or Hash for slice component instead of %s", rb_obj_classname(v));
    }

    if (slice->coords[r] > shape[r] || slice->coords[r] + slice->lengths[r] > shape[r])
      rb_raise(rb_eRangeError, "slice is larger than matrix in dimension %u (slice component %u)", r, t);
  }

  return slice;
}

#ifdef BENCHMARK
/*
 * A simple function used when benchmarking NMatrix.
 */
static double get_time(void) {
  struct timeval t;
  struct timezone tzp;

  gettimeofday(&t, &tzp);

  return t.tv_sec + t.tv_usec*1e-6;
}
#endif

/*
 * The argv parameter will be either 1 or 2 elements.  If 1, could be either
 * initial or dtype.  If 2, is initial and dtype. This function returns the
 * dtype.
 */
static nm::dtype_t interpret_dtype(int argc, VALUE* argv, nm::stype_t stype) {
  int offset;

  switch (argc) {
  	case 1:
  		offset = 0;
  		break;

  	case 2:
  		offset = 1;
  		break;

  	default:
  		rb_raise(rb_eArgError, "Need an initial value or a dtype.");
  		break;
  }

  if (SYMBOL_P(argv[offset])) {
  	return nm_dtype_from_rbsymbol(argv[offset]);

  } else if (TYPE(argv[offset]) == T_STRING) {
  	return nm_dtype_from_rbstring(StringValue(argv[offset]));

  } else if (stype == nm::YALE_STORE) {
  	rb_raise(rb_eArgError, "Yale storage class requires a dtype.");

  } else {
  	return nm_dtype_guess(argv[0]);
  }
}

/*
 * Convert an Ruby value or an array of Ruby values into initial C values.
 */
static void* interpret_initial_value(VALUE arg, nm::dtype_t dtype) {
  unsigned int index;
  void* init_val;

  if (TYPE(arg) == T_ARRAY) {
  	// Array
    init_val = ALLOC_N(char, DTYPE_SIZES[dtype] * RARRAY_LEN(arg));
    NM_CHECK_ALLOC(init_val);
    for (index = 0; index < RARRAY_LEN(arg); ++index) {
    	rubyval_to_cval(RARRAY_PTR(arg)[index], dtype, (char*)init_val + (index * DTYPE_SIZES[dtype]));
    }

  } else {
  	// Single value
    init_val = rubyobj_to_cval(arg, dtype);
  }

  return init_val;
}

/*
 * Convert the shape argument, which may be either a Ruby value or an array of
 * Ruby values, into C values.  The second argument is where the dimensionality
 * of the matrix will be stored.  The function itself returns a pointer to the
 * array describing the shape, which must be freed manually.
 */
static size_t* interpret_shape(VALUE arg, size_t* dim) {
  size_t* shape;

  if (TYPE(arg) == T_ARRAY) {
    *dim = RARRAY_LEN(arg);
    shape = ALLOC_N(size_t, *dim);

    for (size_t index = 0; index < *dim; ++index) {
      shape[index] = FIX2UINT( RARRAY_PTR(arg)[index] );
    }

  } else if (FIXNUM_P(arg)) {
    *dim = 2;
    shape = ALLOC_N(size_t, *dim);

    shape[0] = FIX2UINT(arg);
    shape[1] = FIX2UINT(arg);

  } else {
    rb_raise(rb_eArgError, "Expected an array of numbers or a single Fixnum for matrix shape");
  }

  return shape;
}

/*
 * Convert a Ruby symbol or string into an storage type.
 */
static nm::stype_t interpret_stype(VALUE arg) {
  if (SYMBOL_P(arg)) {
  	return nm_stype_from_rbsymbol(arg);

  } else if (TYPE(arg) == T_STRING) {
  	return nm_stype_from_rbstring(StringValue(arg));

  } else {
  	rb_raise(rb_eArgError, "Expected storage type");
  }
}

//////////////////
// Math Helpers //
//////////////////

STORAGE* matrix_storage_cast_alloc(NMATRIX* matrix, nm::dtype_t new_dtype) {
  if (matrix->storage->dtype == new_dtype && !is_ref(matrix))
    return matrix->storage;

  CAST_TABLE(cast_copy_storage);
  return cast_copy_storage[matrix->stype][matrix->stype](matrix->storage, new_dtype, NULL);
}

STORAGE_PAIR binary_storage_cast_alloc(NMATRIX* left_matrix, NMATRIX* right_matrix) {
  STORAGE_PAIR casted;
  nm::dtype_t new_dtype = Upcast[left_matrix->storage->dtype][right_matrix->storage->dtype];

  casted.left  = matrix_storage_cast_alloc(left_matrix, new_dtype);
  casted.right = matrix_storage_cast_alloc(right_matrix, new_dtype);

  return casted;
}

static VALUE matrix_multiply_scalar(NMATRIX* left, VALUE scalar) {
  rb_raise(rb_eNotImpError, "matrix-scalar multiplication not implemented yet");
  return Qnil;
}

static VALUE matrix_multiply(NMATRIX* left, NMATRIX* right) {
  ///TODO: multiplication for non-dense and/or non-decimal matrices

  // Make sure both of our matrices are of the correct type.
  STORAGE_PAIR casted = binary_storage_cast_alloc(left, right);

  size_t*  resulting_shape   = ALLOC_N(size_t, 2);
  resulting_shape[0] = left->storage->shape[0];
  resulting_shape[1] = right->storage->shape[1];

  // Sometimes we only need to use matrix-vector multiplication (e.g., GEMM versus GEMV). Find out.
  bool vector = false;
  if (resulting_shape[1] == 1) vector = true;

  static STORAGE* (*storage_matrix_multiply[nm::NUM_STYPES])(const STORAGE_PAIR&, size_t*, bool) = {
    nm_dense_storage_matrix_multiply,
    nm_list_storage_matrix_multiply,
    nm_yale_storage_matrix_multiply
  };

  STORAGE* resulting_storage = storage_matrix_multiply[left->stype](casted, resulting_shape, vector);
  NMATRIX* result = nm_create(left->stype, resulting_storage);

  // Free any casted-storage we created for the multiplication.
  // TODO: Can we make the Ruby GC take care of this stuff now that we're using it?
  // If we did that, we night not have to re-create these every time, right? Or wrong? Need to do
  // more research.
  static void (*free_storage[nm::NUM_STYPES])(STORAGE*) = {
    nm_dense_storage_delete,
    nm_list_storage_delete,
    nm_yale_storage_delete
  };

  if (left->storage != casted.left)   free_storage[result->stype](casted.left);
  if (right->storage != casted.right) free_storage[result->stype](casted.right);

  if (result) return Data_Wrap_Struct(cNMatrix, nm_mark, nm_delete, result);
  return Qnil; // Only if we try to multiply list matrices should we return Qnil.
}

/*
 * Calculate the exact determinant of a dense matrix.
 *
 * Returns nil for dense matrices which are not square or number of dimensions other than 2.
 *
 * Note: Currently only implemented for 2x2 and 3x3 matrices.
 */
static VALUE nm_det_exact(VALUE self) {
  if (NM_STYPE(self) != nm::DENSE_STORE) rb_raise(nm_eStorageTypeError, "can only calculate exact determinant for dense matrices");

  if (NM_DIM(self) != 2 || NM_SHAPE0(self) != NM_SHAPE1(self)) return Qnil;

  // Calculate the determinant and then assign it to the return value
  void* result = ALLOCA_N(char, DTYPE_SIZES[NM_DTYPE(self)]);
  nm_math_det_exact(NM_SHAPE0(self), NM_STORAGE_DENSE(self)->elements, NM_SHAPE0(self), NM_DTYPE(self), result);

  return rubyobj_from_cval(result, NM_DTYPE(self)).rval;
}

/////////////////
// Exposed API //
/////////////////

/*
 * Create a dense matrix. Used by the NMatrix GSL fork. Unlike nm_create, this one copies all of the
 * arrays and such passed in -- so you don't have to allocate and pass a new shape object for every
 * matrix you want to create, for example. Same goes for elements.
 *
 * Returns a properly-wrapped Ruby object as a VALUE.
 *
 * *** Note that this function is for API only. Please do not use it internally.
 *
 * TODO: Add a column-major option for libraries that use column-major matrices.
 */
VALUE rb_nmatrix_dense_create(nm::dtype_t dtype, size_t* shape, size_t dim, void* elements, size_t length) {
  NMATRIX* nm;
  size_t nm_dim;
  size_t* shape_copy;

  // Do not allow a dim of 1. Treat it as a column or row matrix.
  if (dim == 1) {
    nm_dim				= 2;
    shape_copy		= ALLOC_N(size_t, nm_dim);
    shape_copy[0]	= shape[0];
    shape_copy[1]	= 1;

  } else {
    nm_dim			= dim;
    shape_copy	= ALLOC_N(size_t, nm_dim);
    memcpy(shape_copy, shape, sizeof(size_t)*nm_dim);
  }

  // Copy elements
  void* elements_copy = ALLOC_N(char, DTYPE_SIZES[dtype]*length);
  memcpy(elements_copy, elements, DTYPE_SIZES[dtype]*length);

  // allocate and create the matrix and its storage
  nm = nm_create(nm::DENSE_STORE, nm_dense_storage_create(dtype, shape_copy, dim, elements_copy, length));

  // tell Ruby about the matrix and its storage, particularly how to garbage collect it.
  return Data_Wrap_Struct(cNMatrix, nm_mark, nm_delete, nm);
}

/*
 * Create a dense vector. Used by the NMatrix GSL fork.
 *
 * Basically just a convenience wrapper for rb_nmatrix_dense_create().
 *
 * Returns a properly-wrapped Ruby NMatrix object as a VALUE. Included for backwards compatibility
 * for when NMatrix had an NVector class.
 */
VALUE rb_nvector_dense_create(nm::dtype_t dtype, void* elements, size_t length) {
  size_t dim = 1, shape = length;
  return rb_nmatrix_dense_create(dtype, &shape, dim, elements, length);
}
