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
// == math.h
//
// Header file for math functions, interfacing with BLAS, etc.
//
// For instructions on adding CBLAS and CLAPACK functions, see the
// beginning of math.cpp.
//
// Some of these functions are from ATLAS. Here is the license for
// ATLAS:
//
/*
 *             Automatically Tuned Linear Algebra Software v3.8.4
 *                    (C) Copyright 1999 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef MATH_H
#define MATH_H

/*
 * Standard Includes
 */

extern "C" { // These need to be in an extern "C" block or you'll get all kinds of undefined symbol errors.
  #include <cblas.h>

  #ifdef HAVE_CLAPACK_H
    #include <clapack.h>
  #endif
}

#include <algorithm> // std::min, std::max
#include <limits> // std::numeric_limits

/*
 * Project Includes
 */

/*
 * Macros
 */
#define REAL_RECURSE_LIMIT 4

/*
 * Data
 */


extern "C" {
  /*
   * C accessors.
   */
  void nm_math_det_exact(const int M, const void* elements, const int lda, nm::dtype_t dtype, void* result);
  void nm_math_transpose_generic(const size_t M, const size_t N, const void* A, const int lda, void* B, const int ldb, size_t element_size);
  void nm_math_init_blas(void);

}


namespace nm {
  namespace math {

/*
 * Types
 */


/*
 * Functions
 */


template <typename DType>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const DType* alpha, const DType* A, const int lda, const DType* beta, DType* C, const int ldc) {
  rb_raise(rb_eNotImpError, "syrk not yet implemented for non-BLAS dtypes");
}

template <typename DType>
inline void herk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const DType* alpha, const DType* A, const int lda, const DType* beta, DType* C, const int ldc) {
  rb_raise(rb_eNotImpError, "herk not yet implemented for non-BLAS dtypes");
}

template <>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const float* alpha, const float* A, const int lda, const float* beta, float* C, const int ldc) {
  cblas_ssyrk(Order, Uplo, Trans, N, K, *alpha, A, lda, *beta, C, ldc);
}

template <>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const double* alpha, const double* A, const int lda, const double* beta, double* C, const int ldc) {
  cblas_dsyrk(Order, Uplo, Trans, N, K, *alpha, A, lda, *beta, C, ldc);
}

template <>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const Complex64* alpha, const Complex64* A, const int lda, const Complex64* beta, Complex64* C, const int ldc) {
  cblas_csyrk(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template <>
inline void syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const Complex128* alpha, const Complex128* A, const int lda, const Complex128* beta, Complex128* C, const int ldc) {
  cblas_zsyrk(Order, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


template <>
inline void herk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const Complex64* alpha, const Complex64* A, const int lda, const Complex64* beta, Complex64* C, const int ldc) {
  cblas_cherk(Order, Uplo, Trans, N, K, alpha->r, A, lda, beta->r, C, ldc);
}

template <>
inline void herk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N,
                 const int K, const Complex128* alpha, const Complex128* A, const int lda, const Complex128* beta, Complex128* C, const int ldc) {
  cblas_zherk(Order, Uplo, Trans, N, K, alpha->r, A, lda, beta->r, C, ldc);
}


template <typename DType>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const DType* alpha,
                 const DType* A, const int lda, DType* B, const int ldb) {
  rb_raise(rb_eNotImpError, "trmm not yet implemented for non-BLAS dtypes");
}

template <>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const float* alpha,
                 const float* A, const int lda, float* B, const int ldb) {
  cblas_strmm(order, side, uplo, ta, diag, m, n, *alpha, A, lda, B, ldb);
}

template <>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const double* alpha,
                 const double* A, const int lda, double* B, const int ldb) {
  cblas_dtrmm(order, side, uplo, ta, diag, m, n, *alpha, A, lda, B, ldb);
}

template <>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const Complex64* alpha,
                 const Complex64* A, const int lda, Complex64* B, const int ldb) {
  cblas_ctrmm(order, side, uplo, ta, diag, m, n, alpha, A, lda, B, ldb);
}

template <>
inline void trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const Complex128* alpha,
                 const Complex128* A, const int lda, Complex128* B, const int ldb) {
  cblas_ztrmm(order, side, uplo, ta, diag, m, n, alpha, A, lda, B, ldb);
}



// Yale: numeric matrix multiply c=a*b
template <typename DType>
inline void numbmm(const unsigned int n, const unsigned int m, const unsigned int l, const IType* ia, const IType* ja, const DType* a, const bool diaga,
            const IType* ib, const IType* jb, const DType* b, const bool diagb, IType* ic, IType* jc, DType* c, const bool diagc) {
  const unsigned int max_lmn = std::max(std::max(m, n), l);
  IType next[max_lmn];
  DType sums[max_lmn];

  DType v;

  IType head, length, temp, ndnz = 0;
  IType minmn = std::min(m,n);
  IType minlm = std::min(l,m);

  for (IType idx = 0; idx < max_lmn; ++idx) { // initialize scratch arrays
    next[idx] = std::numeric_limits<IType>::max();
    sums[idx] = 0;
  }

  for (IType i = 0; i < n; ++i) { // walk down the rows
    head = std::numeric_limits<IType>::max()-1; // head gets assigned as whichever column of B's row j we last visited
    length = 0;

    for (IType jj = ia[i]; jj <= ia[i+1]; ++jj) { // walk through entries in each row
      IType j;

      if (jj == ia[i+1]) { // if we're in the last entry for this row:
        if (!diaga || i >= minmn) continue;
        j   = i;      // if it's a new Yale matrix, and last entry, get the diagonal position (j) and entry (ajj)
        v   = a[i];
      } else {
        j   = ja[jj]; // if it's not the last entry for this row, get the column (j) and entry (ajj)
        v   = a[jj];
      }

      for (IType kk = ib[j]; kk <= ib[j+1]; ++kk) {

        IType k;

        if (kk == ib[j+1]) { // Get the column id for that entry
          if (!diagb || j >= minlm) continue;
          k  = j;
          sums[k] += v*b[k];
        } else {
          k  = jb[kk];
          sums[k] += v*b[kk];
        }

        if (next[k] == std::numeric_limits<IType>::max()) {
          next[k] = head;
          head    = k;
          ++length;
        }
      } // end of kk loop
    } // end of jj loop

    for (IType jj = 0; jj < length; ++jj) {
      if (sums[head] != 0) {
        if (diagc && head == i) {
          c[head] = sums[head];
        } else {
          jc[n+1+ndnz] = head;
          c[n+1+ndnz]  = sums[head];
          ++ndnz;
        }
      }

      temp = head;
      head = next[head];

      next[temp] = std::numeric_limits<IType>::max();
      sums[temp] = 0;
    }

    ic[i+1] = n+1+ndnz;
  }
} /* numbmm_ */


/*
template <typename DType, typename IType>
inline void new_yale_matrix_multiply(const unsigned int m, const IType* ija, const DType* a, const IType* ijb, const DType* b, YALE_STORAGE* c_storage) {
  unsigned int n = c_storage->shape[0],
               l = c_storage->shape[1];

  // Create a working vector of dimension max(m,l,n) and initial value IType::max():
  std::vector<IType> mask(std::max(std::max(m,l),n), std::numeric_limits<IType>::max());

  for (IType i = 0; i < n; ++i) { // A.rows.each_index do |i|

    IType j, k;
    size_t ndnz;

    for (IType jj = ija[i]; jj <= ija[i+1]; ++jj) { // walk through column pointers for row i of A
      j = (jj == ija[i+1]) ? i : ija[jj];   // Get the current column index (handle diagonals last)

      if (j >= m) {
        if (j == ija[jj]) rb_raise(rb_eIndexError, "ija array for left-hand matrix contains an out-of-bounds column index %u at position %u", jj, j);
        else              break;
      }

      for (IType kk = ijb[j]; kk <= ijb[j+1]; ++kk) { // walk through column pointers for row j of B
        if (j >= m) continue; // first of all, does B *have* a row j?
        k = (kk == ijb[j+1]) ? j : ijb[kk];   // Get the current column index (handle diagonals last)

        if (k >= l) {
          if (k == ijb[kk]) rb_raise(rb_eIndexError, "ija array for right-hand matrix contains an out-of-bounds column index %u at position %u", kk, k);
          else              break;
        }

        if (mask[k] == )
      }

    }
  }
}
*/

// Yale: Symbolic matrix multiply c=a*b
inline size_t symbmm(const unsigned int n, const unsigned int m, const unsigned int l, const IType* ia, const IType* ja, const bool diaga,
            const IType* ib, const IType* jb, const bool diagb, IType* ic, const bool diagc) {
  unsigned int max_lmn = std::max(std::max(m,n), l);
  IType mask[max_lmn];  // INDEX in the SMMP paper.
  IType j, k; /* Local variables */
  size_t ndnz = n;

  for (IType idx = 0; idx < max_lmn; ++idx)
    mask[idx] = std::numeric_limits<IType>::max();

  if (ic) { // Only write to ic if it's supplied; otherwise, we're just counting.
    if (diagc)  ic[0] = n+1;
    else        ic[0] = 0;
  }

  IType minmn = std::min(m,n);
  IType minlm = std::min(l,m);

  for (IType i = 0; i < n; ++i) { // MAIN LOOP: through rows

    for (IType jj = ia[i]; jj <= ia[i+1]; ++jj) { // merge row lists, walking through columns in each row

      // j <- column index given by JA[jj], or handle diagonal.
      if (jj == ia[i+1]) { // Don't really do it the last time -- just handle diagonals in a new yale matrix.
        if (!diaga || i >= minmn) continue;
        j = i;
      } else j = ja[jj];

      for (IType kk = ib[j]; kk <= ib[j+1]; ++kk) { // Now walk through columns K of row J in matrix B.
        if (kk == ib[j+1]) {
          if (!diagb || j >= minlm) continue;
          k = j;
        } else k = jb[kk];

        if (mask[k] != i) {
          mask[k] = i;
          ++ndnz;
        }
      }
    }

    if (diagc && mask[i] == std::numeric_limits<IType>::max()) --ndnz;

    if (ic) ic[i+1] = ndnz;
  }

  return ndnz;
} /* symbmm_ */


// In-place quicksort (from Wikipedia) -- called by smmp_sort_columns, below. All functions are inclusive of left, right.
namespace smmp_sort {
  const size_t THRESHOLD = 4;  // switch to insertion sort for 4 elements or fewer

  template <typename DType>
  void print_array(DType* vals, IType* array, IType left, IType right) {
    for (IType i = left; i <= right; ++i) {
      std::cerr << array[i] << ":" << vals[i] << "  ";
    }
    std::cerr << std::endl;
  }

  template <typename DType>
  IType partition(DType* vals, IType* array, IType left, IType right, IType pivot) {
    IType pivotJ = array[pivot];
    DType pivotV = vals[pivot];

    // Swap pivot and right
    array[pivot] = array[right];
    vals[pivot]  = vals[right];
    array[right] = pivotJ;
    vals[right]  = pivotV;

    IType store = left;
    for (IType idx = left; idx < right; ++idx) {
      if (array[idx] <= pivotJ) {
        // Swap i and store
        std::swap(array[idx], array[store]);
        std::swap(vals[idx],  vals[store]);
        ++store;
      }
    }

    std::swap(array[store], array[right]);
    std::swap(vals[store],  vals[right]);

    return store;
  }

  // Recommended to use the median of left, right, and mid for the pivot.
  template <typename I>
  inline I median(I a, I b, I c) {
    if (a < b) {
      if (b < c) return b; // a b c
      if (a < c) return c; // a c b
                 return a; // c a b

    } else { // a > b
      if (a < c) return a; // b a c
      if (b < c) return c; // b c a
                 return b; // c b a
    }
  }


  // Insertion sort is more efficient than quicksort for small N
  template <typename DType>
  void insertion_sort(DType* vals, IType* array, IType left, IType right) {
    for (IType idx = left; idx <= right; ++idx) {
      IType col_to_insert = array[idx];
      DType val_to_insert = vals[idx];

      IType hole_pos = idx;
      for (; hole_pos > left && col_to_insert < array[hole_pos-1]; --hole_pos) {
        array[hole_pos] = array[hole_pos - 1];  // shift the larger column index up
        vals[hole_pos]  = vals[hole_pos - 1];   // value goes along with it
      }

      array[hole_pos] = col_to_insert;
      vals[hole_pos]  = val_to_insert;
    }
  }


  template <typename DType>
  void quicksort(DType* vals, IType* array, IType left, IType right) {

    if (left < right) {
      if (right - left < THRESHOLD) {
        insertion_sort(vals, array, left, right);
      } else {
        // choose any pivot such that left < pivot < right
        IType pivot = median<IType>(left, right, (IType)(((unsigned long)left + (unsigned long)right) / 2));
        pivot = partition(vals, array, left, right, pivot);

        // recursively sort elements smaller than the pivot
        quicksort<DType>(vals, array, left, pivot-1);

        // recursively sort elements at least as big as the pivot
        quicksort<DType>(vals, array, pivot+1, right);
      }
    }
  }


}; // end of namespace smmp_sort


/*
 * For use following symbmm and numbmm. Sorts the matrix entries in each row according to the column index.
 * This utilizes quicksort, which is an in-place unstable sort (since there are no duplicate entries, we don't care
 * about stability).
 *
 * TODO: It might be worthwhile to do a test for free memory, and if available, use an unstable sort that isn't in-place.
 *
 * TODO: It's actually probably possible to write an even faster sort, since symbmm/numbmm are not producing a random
 * ordering. If someone is doing a lot of Yale matrix multiplication, it might benefit them to consider even insertion
 * sort.
 */
template <typename DType>
inline void smmp_sort_columns(const size_t n, const IType* ia, IType* ja, DType* a) {
  for (size_t i = 0; i < n; ++i) {
    if (ia[i+1] - ia[i] < 2) continue; // no need to sort rows containing only one or two elements.
    else if (ia[i+1] - ia[i] <= smmp_sort::THRESHOLD) {
      smmp_sort::insertion_sort<DType>(a, ja, ia[i], ia[i+1]-1); // faster for small rows
    } else {
      smmp_sort::quicksort<DType>(a, ja, ia[i], ia[i+1]-1);      // faster for large rows (and may call insertion_sort as well)
    }
  }
}


/*
 * From ATLAS 3.8.0:
 *
 * Computes one of two LU factorizations based on the setting of the Order
 * parameter, as follows:
 * ----------------------------------------------------------------------------
 *                       Order == CblasColMajor
 * Column-major factorization of form
 *   A = P * L * U
 * where P is a row-permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapazoidal if M > N), and U is upper triangular
 * (upper trapazoidal if M < N).
 *
 * ----------------------------------------------------------------------------
 *                       Order == CblasRowMajor
 * Row-major factorization of form
 *   A = P * L * U
 * where P is a column-permutation matrix, L is lower triangular (lower
 * trapazoidal if M > N), and U is upper triangular with unit diagonals (upper
 * trapazoidal if M < N).
 *
 * ============================================================================
 * Let IERR be the return value of the function:
 *    If IERR == 0, successful exit.
 *    If (IERR < 0) the -IERR argument had an illegal value
 *    If (IERR > 0 && Order == CblasColMajor)
 *       U(i-1,i-1) is exactly zero.  The factorization has been completed,
 *       but the factor U is exactly singular, and division by zero will
 *       occur if it is used to solve a system of equations.
 *    If (IERR > 0 && Order == CblasRowMajor)
 *       L(i-1,i-1) is exactly zero.  The factorization has been completed,
 *       but the factor L is exactly singular, and division by zero will
 *       occur if it is used to solve a system of equations.
 */
template <typename DType>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, DType* A, const int lda) {
#ifdef HAVE_CLAPACK_H
  rb_raise(rb_eNotImpError, "not yet implemented for non-BLAS dtypes");
#else
  rb_raise(rb_eNotImpError, "only LAPACK version implemented thus far");
#endif
  return 0;
}

#ifdef HAVE_CLAPACK_H
template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, float* A, const int lda) {
  return clapack_spotrf(order, uplo, N, A, lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, double* A, const int lda) {
  return clapack_dpotrf(order, uplo, N, A, lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex64* A, const int lda) {
  return clapack_cpotrf(order, uplo, N, reinterpret_cast<void*>(A), lda);
}

template <>
inline int potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex128* A, const int lda) {
  return clapack_zpotrf(order, uplo, N, reinterpret_cast<void*>(A), lda);
}
#endif



// Copies an upper row-major array from U, zeroing U; U is unit, so diagonal is not copied.
//
// From ATLAS 3.8.0.
template <typename DType>
static inline void trcpzeroU(const int M, const int N, DType* U, const int ldu, DType* C, const int ldc) {

  for (int i = 0; i != M; ++i) {
    for (int j = i+1; j < N; ++j) {
      C[j] = U[j];
      U[j] = 0;
    }

    C += ldc;
    U += ldu;
  }
}


/*
 * Un-comment the following lines when we figure out how to calculate NB for each of the ATLAS-derived
 * functions. This is probably really complicated.
 *
 * Also needed: ATL_MulByNB, ATL_DivByNB (both defined in the build process for ATLAS), and ATL_mmMU.
 *
 */

/*

template <bool RowMajor, bool Upper, typename DType>
static int trtri_4(const enum CBLAS_DIAG Diag, DType* A, const int lda) {

  if (RowMajor) {
    DType *pA0 = A, *pA1 = A+lda, *pA2 = A+2*lda, *pA3 = A+3*lda;
    DType tmp;
    if (Upper) {
      DType A01 = pA0[1], A02 = pA0[2], A03 = pA0[3],
                          A12 = pA1[2], A13 = pA1[3],
                                        A23 = pA2[3];

      if (Diag == CblasNonUnit) {
        pA0->inverse();
        (pA1+1)->inverse();
        (pA2+2)->inverse();
        (pA3+3)->inverse();

        pA0[1] = -A01 * pA1[1] * pA0[0];
        pA1[2] = -A12 * pA2[2] * pA1[1];
        pA2[3] = -A23 * pA3[3] * pA2[2];

        pA0[2] = -(A01 * pA1[2] + A02 * pA2[2]) * pA0[0];
        pA1[3] = -(A12 * pA2[3] + A13 * pA3[3]) * pA1[1];

        pA0[3] = -(A01 * pA1[3] + A02 * pA2[3] + A03 * pA3[3]) * pA0[0];

      } else {

        pA0[1] = -A01;
        pA1[2] = -A12;
        pA2[3] = -A23;

        pA0[2] = -(A01 * pA1[2] + A02);
        pA1[3] = -(A12 * pA2[3] + A13);

        pA0[3] = -(A01 * pA1[3] + A02 * pA2[3] + A03);
      }

    } else { // Lower
      DType A10 = pA1[0],
            A20 = pA2[0], A21 = pA2[1],
            A30 = PA3[0], A31 = pA3[1], A32 = pA3[2];
      DType *B10 = pA1,
            *B20 = pA2,
            *B30 = pA3,
            *B21 = pA2+1,
            *B31 = pA3+1,
            *B32 = pA3+2;


      if (Diag == CblasNonUnit) {
        pA0->inverse();
        (pA1+1)->inverse();
        (pA2+2)->inverse();
        (pA3+3)->inverse();

        *B10 = -A10 * pA0[0] * pA1[1];
        *B21 = -A21 * pA1[1] * pA2[2];
        *B32 = -A32 * pA2[2] * pA3[3];
        *B20 = -(A20 * pA0[0] + A21 * (*B10)) * pA2[2];
        *B31 = -(A31 * pA1[1] + A32 * (*B21)) * pA3[3];
        *B30 = -(A30 * pA0[0] + A31 * (*B10) + A32 * (*B20)) * pA3;
      } else {
        *B10 = -A10;
        *B21 = -A21;
        *B32 = -A32;
        *B20 = -(A20 + A21 * (*B10));
        *B31 = -(A31 + A32 * (*B21));
        *B30 = -(A30 + A31 * (*B10) + A32 * (*B20));
      }
    }

  } else {
    rb_raise(rb_eNotImpError, "only row-major implemented at this time");
  }

  return 0;

}


template <bool RowMajor, bool Upper, typename DType>
static int trtri_3(const enum CBLAS_DIAG Diag, DType* A, const int lda) {

  if (RowMajor) {

    DType tmp;

    if (Upper) {
      DType A01 = pA0[1], A02 = pA0[2], A03 = pA0[3],
                          A12 = pA1[2], A13 = pA1[3];

      DType *B01 = pA0 + 1,
            *B02 = pA0 + 2,
            *B12 = pA1 + 2;

      if (Diag == CblasNonUnit) {
        pA0->inverse();
        (pA1+1)->inverse();
        (pA2+2)->inverse();

        *B01 = -A01 * pA1[1] * pA0[0];
        *B12 = -A12 * pA2[2] * pA1[1];
        *B02 = -(A01 * (*B12) + A02 * pA2[2]) * pA0[0];
      } else {
        *B01 = -A01;
        *B12 = -A12;
        *B02 = -(A01 * (*B12) + A02);
      }

    } else { // Lower
      DType *pA0=A, *pA1=A+lda, *pA2=A+2*lda;
      DType A10=pA1[0],
            A20=pA2[0], A21=pA2[1];

      DType *B10 = pA1,
            *B20 = pA2;
            *B21 = pA2+1;

      if (Diag == CblasNonUnit) {
        pA0->inverse();
        (pA1+1)->inverse();
        (pA2+2)->inverse();
        *B10 = -A10 * pA0[0] * pA1[1];
        *B21 = -A21 * pA1[1] * pA2[2];
        *B20 = -(A20 * pA0[0] + A21 * (*B10)) * pA2[2];
      } else {
        *B10 = -A10;
        *B21 = -A21;
        *B20 = -(A20 + A21 * (*B10));
      }
    }


  } else {
    rb_raise(rb_eNotImpError, "only row-major implemented at this time");
  }

  return 0;

}

template <bool RowMajor, bool Upper, bool Real, typename DType>
static void trtri(const enum CBLAS_DIAG Diag, const int N, DType* A, const int lda) {
  DType *Age, *Atr;
  DType tmp;
  int Nleft, Nright;

  int ierr = 0;

  static const DType ONE = 1;
  static const DType MONE -1;
  static const DType NONE = -1;

  if (RowMajor) {

    // FIXME: Use REAL_RECURSE_LIMIT here for float32 and float64 (instead of 1)
    if ((Real && N > REAL_RECURSE_LIMIT) || (N > 1)) {
      Nleft = N >> 1;
#ifdef NB
      if (Nleft > NB) NLeft = ATL_MulByNB(ATL_DivByNB(Nleft));
#endif

      Nright = N - Nleft;

      if (Upper) {
        Age = A + Nleft;
        Atr = A + (Nleft * (lda+1));

        nm::math::trsm<DType>(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, Diag,
                              Nleft, Nright, ONE, Atr, lda, Age, lda);

        nm::math::trsm<DType>(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, Diag,
                              Nleft, Nright, MONE, A, lda, Age, lda);

      } else { // Lower
        Age = A + ((Nleft*lda));
        Atr = A + (Nleft * (lda+1));

        nm::math::trsm<DType>(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, Diag,
                              Nright, Nleft, ONE, A, lda, Age, lda);
        nm::math::trsm<DType>(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, Diag,
                              Nright, Nleft, MONE, Atr, lda, Age, lda);
      }

      ierr = trtri<RowMajor,Upper,Real,DType>(Diag, Nleft, A, lda);
      if (ierr) return ierr;

      ierr = trtri<RowMajor,Upper,Real,DType>(Diag, Nright, Atr, lda);
      if (ierr) return ierr + Nleft;

    } else {
      if (Real) {
        if (N == 4) {
          return trtri_4<RowMajor,Upper,Real,DType>(Diag, A, lda);
        } else if (N == 3) {
          return trtri_3<RowMajor,Upper,Real,DType>(Diag, A, lda);
        } else if (N == 2) {
          if (Diag == CblasNonUnit) {
            A->inverse();
            (A+(lda+1))->inverse();

            if (Upper) {
              *(A+1)     *=   *A;         // TRI_MUL
              *(A+1)     *=   *(A+lda+1); // TRI_MUL
            } else {
              *(A+lda)   *=   *A;         // TRI_MUL
              *(A+lda)   *=   *(A+lda+1); // TRI_MUL
            }
          }

          if (Upper) *(A+1)   = -*(A+1);      // TRI_NEG
          else       *(A+lda) = -*(A+lda);    // TRI_NEG
        } else if (Diag == CblasNonUnit) A->inverse();
      } else { // not real
        if (Diag == CblasNonUnit) A->inverse();
      }
    }

  } else {
    rb_raise(rb_eNotImpError, "only row-major implemented at this time");
  }

  return ierr;
}


template <bool RowMajor, bool Real, typename DType>
int getri(const int N, DType* A, const int lda, const int* ipiv, DType* wrk, const int lwrk) {

  if (!RowMajor) rb_raise(rb_eNotImpError, "only row-major implemented at this time");

  int jb, nb, I, ndown, iret;

  const DType ONE = 1, NONE = -1;

  int iret = trtri<RowMajor,false,Real,DType>(CblasNonUnit, N, A, lda);
  if (!iret && N > 1) {
    jb = lwrk / N;
    if (jb >= NB) nb = ATL_MulByNB(ATL_DivByNB(jb));
    else if (jb >= ATL_mmMU) nb = (jb/ATL_mmMU)*ATL_mmMU;
    else nb = jb;
    if (!nb) return -6; // need at least 1 row of workspace

    // only first iteration will have partial block, unroll it

    jb = N - (N/nb) * nb;
    if (!jb) jb = nb;
    I = N - jb;
    A += lda * I;
    trcpzeroU<DType>(jb, jb, A+I, lda, wrk, jb);
    nm::math::trsm<DType>(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
                          jb, N, ONE, wrk, jb, A, lda);

    if (I) {
      do {
        I -= nb;
        A -= nb * lda;
        ndown = N-I;
        trcpzeroU<DType>(nb, ndown, A+I, lda, wrk, ndown);
        nm::math::gemm<DType>(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
                              nb, N, ONE, wrk, ndown, A, lda);
      } while (I);
    }

    // Apply row interchanges

    for (I = N - 2; I >= 0; --I) {
      jb = ipiv[I];
      if (jb != I) nm::math::swap<DType>(N, A+I*lda, 1, A+jb*lda, 1);
    }
  }

  return iret;
}
*/



template <bool is_complex, typename DType>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, DType* A, const int lda) {

  int Nleft, Nright;
  const DType ONE = 1;
  DType *G, *U0 = A, *U1;

  if (N > 1) {
    Nleft = N >> 1;
    #ifdef NB
      if (Nleft > NB) Nleft = ATL_MulByNB(ATL_DivByNB(Nleft));
    #endif

    Nright = N - Nleft;

    // FIXME: There's a simpler way to write this next block, but I'm way too tired to work it out right now.
    if (uplo == CblasUpper) {
      if (order == CblasRowMajor) {
        G = A + Nleft;
        U1 = G + Nleft * lda;
      } else {
        G = A + Nleft * lda;
        U1 = G + Nleft;
      }
    } else {
      if (order == CblasRowMajor) {
        G = A + Nleft * lda;
        U1 = G + Nleft;
      } else {
        G = A + Nleft;
        U1 = G + Nleft * lda;
      }
    }

    lauum<is_complex, DType>(order, uplo, Nleft, U0, lda);

    if (is_complex) {

      nm::math::herk<DType>(order, uplo,
                            uplo == CblasLower ? CblasConjTrans : CblasNoTrans,
                            Nleft, Nright, &ONE, G, lda, &ONE, U0, lda);

      nm::math::trmm<DType>(order, CblasLeft, uplo, CblasConjTrans, CblasNonUnit, Nright, Nleft, &ONE, U1, lda, G, lda);
    } else {
      nm::math::syrk<DType>(order, uplo,
                            uplo == CblasLower ? CblasTrans : CblasNoTrans,
                            Nleft, Nright, &ONE, G, lda, &ONE, U0, lda);

      nm::math::trmm<DType>(order, CblasLeft, uplo, CblasTrans, CblasNonUnit, Nright, Nleft, &ONE, U1, lda, G, lda);
    }
    lauum<is_complex, DType>(order, uplo, Nright, U1, lda);

  } else {
    *A = *A * *A;
  }
}


#ifdef HAVE_CLAPACK_H
template <bool is_complex>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, float* A, const int lda) {
  clapack_slauum(order, uplo, N, A, lda);
}

template <bool is_complex>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, double* A, const int lda) {
  clapack_dlauum(order, uplo, N, A, lda);
}

template <bool is_complex>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex64* A, const int lda) {
  clapack_clauum(order, uplo, N, A, lda);
}

template <bool is_complex>
inline void lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int N, Complex128* A, const int lda) {
  clapack_zlauum(order, uplo, N, A, lda);
}
#endif


/*
* Function signature conversion for calling LAPACK's lauum functions as directly as possible.
*
* For documentation: http://www.netlib.org/lapack/double/dlauum.f
*
* This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
*/
template <bool is_complex, typename DType>
inline int clapack_lauum(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, void* a, const int lda) {
  if (n < 0)              rb_raise(rb_eArgError, "n cannot be less than zero, is set to %d", n);
  if (lda < n || lda < 1) rb_raise(rb_eArgError, "lda must be >= max(n,1); lda=%d, n=%d\n", lda, n);

  if (uplo == CblasUpper) lauum<is_complex, DType>(order, uplo, n, reinterpret_cast<DType*>(a), lda);
  else                    lauum<is_complex, DType>(order, uplo, n, reinterpret_cast<DType*>(a), lda);

  return 0;
}




/*
 * Macro for declaring LAPACK specializations of the getrf function.
 *
 * type is the DType; call is the specific function to call; cast_as is what the DType* should be
 * cast to in order to pass it to LAPACK.
 */
#define LAPACK_GETRF(type, call, cast_as)                                     \
template <>                                                                   \
inline int getrf(const enum CBLAS_ORDER Order, const int M, const int N, type * A, const int lda, int* ipiv) { \
  int info = call(Order, M, N, reinterpret_cast<cast_as *>(A), lda, ipiv);    \
  if (!info) return info;                                                     \
  else {                                                                      \
    rb_raise(rb_eArgError, "getrf: problem with argument %d\n", info);        \
    return info;                                                              \
  }                                                                           \
}

/* Specialize for ATLAS types */
/*LAPACK_GETRF(float,      clapack_sgetrf, float)
LAPACK_GETRF(double,     clapack_dgetrf, double)
LAPACK_GETRF(Complex64,  clapack_cgetrf, void)
LAPACK_GETRF(Complex128, clapack_zgetrf, void)
*/



/*
* Function signature conversion for calling LAPACK's potrf functions as directly as possible.
*
* For documentation: http://www.netlib.org/lapack/double/dpotrf.f
*
* This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
*/
template <typename DType>
inline int clapack_potrf(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, void* a, const int lda) {
  return potrf<DType>(order, uplo, n, reinterpret_cast<DType*>(a), lda);
}



template <typename DType>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, DType* a, const int lda) {
  rb_raise(rb_eNotImpError, "potri not yet implemented for non-BLAS dtypes");
  return 0;
}


#ifdef HAVE_CLAPACK_H
template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, float* a, const int lda) {
  return clapack_spotri(order, uplo, n, a, lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, double* a, const int lda) {
  return clapack_dpotri(order, uplo, n, a, lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, Complex64* a, const int lda) {
  return clapack_cpotri(order, uplo, n, reinterpret_cast<void*>(a), lda);
}

template <>
inline int potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, Complex128* a, const int lda) {
  return clapack_zpotri(order, uplo, n, reinterpret_cast<void*>(a), lda);
}
#endif


/*
 * Function signature conversion for calling LAPACK's potri functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/lapack/double/dpotri.f
 *
 * This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
 */
template <typename DType>
inline int clapack_potri(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const int n, void* a, const int lda) {
  return potri<DType>(order, uplo, n, reinterpret_cast<DType*>(a), lda);
}




}} // end namespace nm::math


#endif // MATH_H
