#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "types.h"

#include <initializer_list>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <cassert>


template<typename DType, T_DIM dim> class Vec;
template<typename DType, T_DIM dim_row, T_DIM dim_col> class Mat;

template<typename DType, T_DIM dim>
class Vec {
    void _init( unsigned int d ) {
        if (d) printf("[ERROR] To handle variadic initializers at compile time only.");
	}
	template< class _DType , class ... _DTypes > void _init( T_DIM d , _DType v , _DTypes ... values ) {
		coords[d] = (DType)v;
		if( d+1<dim ) _init( d+1 , values... );
	}
    DType coords[dim] = {};

public:
	// DType coords[dim];
	Vec() { memset( coords , 0 , sizeof(DType)*dim ); }
	Vec( const Vec & v ){ memcpy( coords , v.coords , sizeof(DType)*dim ); }
	template< class ... _DTypes > Vec( _DTypes ... values ) { static_assert( sizeof...(values)==dim || sizeof...(values)==0 , "[ERROR] Invalid number of coefficients" ); _init( 0 , values... ); }
	template< class _DType > Vec( const Vec< _DType , dim >& v ){ for( T_DIM d=0 ; d<dim ; d++ ) coords[d] = (DType) v.coords[d]; }
	inline       DType& operator[] ( T_IDX i )       { assert(0 <= i && i < dim); return coords[i]; }
	inline const DType& operator[] ( T_IDX i ) const { assert(0 <= i && i < dim); return coords[i]; }
	inline       DType& operator() ( T_IDX i )       { assert(0 <= i && i < dim); return coords[i]; }
	inline const DType& operator() ( T_IDX i ) const { assert(0 <= i && i < dim); return coords[i]; }
	inline Vec  operator - ( void ) const { Vec q ; for( T_DIM d=0 ; d<dim ; d++ ) q.coords[d] = - coords[d] ; return q; }

	inline Vec& operator += (Vec<DType ,dim> v) { for( T_DIM d=0 ; d<dim ; d++ ) coords[d] += v.coords[d] ; return *this; }
	inline Vec  operator +  (Vec<DType ,dim> v) const { Vec q ; for( T_DIM d=0 ; d<dim ; d++ ) q.coords[d] = coords[d] + v.coords[d] ; return q; }
	inline Vec& operator -= (Vec<DType ,dim> v) { return (*this)+=(-v); }
	inline Vec  operator -  (Vec<DType ,dim> v) const { return (*this)+ (-v); }
	inline Vec& operator *= (DType r ) { for( T_DIM d=0 ; d<dim ; d++ ) coords[d] *= r ; return *this; }
	inline Vec  operator *  (DType r ) const { Vec q ; for( T_DIM d=0 ; d<dim ; d++ ) q.coords[d] = coords[d] * r ; return q; }
	inline Vec& operator /= (DType r ) { for( T_DIM d=0 ; d<dim ; d++ ) coords[d] /= r ; return *this; }
	inline Vec  operator /  (DType r ) const { Vec q ; for( T_DIM d=0 ; d<dim ; d++ ) q.coords[d] = coords[d] / r ; return q; }
	inline Vec& operator *= (Vec<DType ,dim> v) { for( T_DIM d=0 ; d<dim ; d++ ) coords[d] *= v.coords[d] ; return *this; }
	inline Vec  operator *  (Vec<DType ,dim> v) const { Vec q ; for( T_DIM d=0 ; d<dim ; d++ ) q.coords[d] = coords[d] * v.coords[d] ; return q; }
	inline Vec& operator /= (Vec<DType ,dim> v) { for( T_DIM d=0 ; d<dim ; d++ ) coords[d] /= v.coords[d] ; return *this; }
	inline Vec  operator /  (Vec<DType ,dim> v) const { Vec q ; for( T_DIM d=0 ; d<dim ; d++ ) q.coords[d] = coords[d] / v.coords[d] ; return q; }
	inline bool operator == (const Vec<DType ,dim>& v) const { for( T_DIM d=0 ; d<dim ; d++ ) if (coords[d] != v.coords[d]) { return false; } return true; }
	inline bool operator != (const Vec<DType ,dim>& v) const { return !(*this == v); }

	inline Vec normalized() const { DType norm = std::sqrt(dot( *this , *this )); if (norm > DType(0)) return *this / norm; return *this;}
	
	inline DType dot( const Vec& v2 ) const { return Vec::dot(*this, v2); }
	inline Vec cross( const Vec& v2 ) const { return Vec::cross(*this, v2); }
	inline DType square_norm() const { return dot( *this , *this ); }
	
	inline static DType dot( const Vec& v1 , const Vec& v2 ) { DType dot = 0 ; for( T_DIM d=0 ; d<dim ; d++ ) dot += v1.coords[d] * v2.coords[d] ; return dot; }
	inline static Vec cross( const Vec& v1 , const Vec& v2 );
	inline static DType square_norm(const Vec& v){ return dot( v , v ); }

	friend std::ostream &operator << ( std::ostream &os , const Vec &p ) {
		os << "(";
		for( T_DIM d=0 ; d<dim ; d++ ) {
			if( d ) os << ", ";
			os << p[d];
		}
		return os << ")";
	}
};

template<typename DType, T_DIM dim>
Vec<DType, dim> Vec<DType, dim>::cross(const Vec& v, const Vec& w) {
	static_assert( dim == 3, "[ERROR] this function supports 3d points only" );

    return Vec<DType, dim>(
        v[1] * w[2] - v[2] * w[1],
        v[2] * w[0] - v[0] * w[2],
        v[0] * w[1] - v[1] * w[0]
    );
}



/******************************
 * 2D matrices
******************************/

template<typename DType, T_DIM dim_row, T_DIM dim_col>
class Mat {
    DType entry[dim_row * dim_col] = {};

public:
	Mat() { memset(entry, 0, sizeof(DType)*dim_row*dim_col); }
	Mat(const Mat & mat){ memcpy(entry , mat.entry, sizeof(DType)*dim_row*dim_col); }
	Mat(std::initializer_list<DType> init_list){
        memset(entry, 0, sizeof(DType)*dim_row*dim_col);
        T_DIM d=0; auto iter=init_list.begin();
        for (; d<dim_row*dim_col && iter<init_list.end(); d++, iter++) entry[d] = *iter;
    }
	template<class _DType> Mat(const Mat<_DType, dim_row, dim_col>& mat){ for(T_DIM d=0; d<dim_row*dim_col; d++) entry[d] = (DType) mat.entry[d]; }
	inline       DType& operator() (T_IDX i, T_IDX j)       { assert(0 <= i && i < dim_row && 0 <= j && j < dim_col); return entry[i*dim_col + j]; }
	inline const DType& operator() (T_IDX i, T_IDX j) const { assert(0 <= i && i < dim_row && 0 <= j && j < dim_col); return entry[i*dim_col + j]; }
    inline T_DIM rows() const { return dim_row; }
    inline T_DIM cols() const { return dim_col; }
	inline Mat  operator - ( void ) const { Mat ret ; for( T_DIM d=0 ; d<dim_row*dim_col ; d++ ) ret.entry[d] = -entry[d] ; return ret; }
	inline Mat& operator += (const Mat& mat)       {          for(T_DIM d=0; d<dim_row*dim_col; d++)     entry[d] += mat.entry[d] ; return *this; }
	inline Mat  operator +  (const Mat& mat) const { Mat ret; for(T_DIM d=0; d<dim_row*dim_col; d++) ret.entry[d]  = entry[d] + mat.entry[d]; return ret; }
	inline Mat& operator -= (const Mat& mat)       { return (*this) += (-mat); }
	inline Mat  operator -  (const Mat& v)   const { return (*this) + (-v); }
	inline Mat& operator *= (const DType& r)       {          for(T_DIM d=0; d<dim_row*dim_col; d++)     entry[d] *= r ; return *this; }
	inline Mat  operator *  (const DType& r) const { Mat ret; for(T_DIM d=0; d<dim_row*dim_col; d++) ret.entry[d]  = entry[d] * r; return ret; }
	inline Mat& operator /= (const DType& r)       {          for(T_DIM d=0; d<dim_row*dim_col; d++)     entry[d] /= r ; return *this; }
	inline Mat  operator /  (const DType& r) const { Mat ret; for(T_DIM d=0; d<dim_row*dim_col; d++) ret.entry[d]  = entry[d] / r; return ret; }
	inline Mat& operator *= (const Mat& mat)       {          for(T_DIM d=0; d<dim_row*dim_col; d++)     entry[d] *= mat.entry[d]; return *this; }
	inline Mat  operator *  (const Mat& mat) const { Mat ret; for(T_DIM d=0; d<dim_row*dim_col; d++) ret.entry[d]  = entry[d] * mat.entry[d]; return ret; }
	inline Mat& operator /= (const Mat& mat)       {          for(T_DIM d=0; d<dim_row*dim_col; d++)     entry[d] /= mat.entry[d] ; return *this; }
	inline Mat  operator /  (const Mat& mat) const { Mat ret; for(T_DIM d=0; d<dim_row*dim_col; d++) ret.entry[d]  = entry[d] / mat.entry[d]; return ret; }
	inline bool operator == (const Mat& mat) const {          for(T_DIM d=0; d<dim_row*dim_col; d++) if (entry[d] != mat.entry[d]) { return false; } return true; }
	inline bool operator != (const Mat& mat) const { return !(*this == mat); }

	Mat<DType, dim_col, dim_row> transpose() const;
    template<T_DIM dim_col2> Mat<DType, dim_row, dim_col2> matmul(const Mat<DType, dim_col, dim_col2>& mat) const;
    Vec<DType, dim_row> matmul(const Vec<DType, dim_col>& vec) const;
    DType det() const;

private:
    DType det_1_by_1() const { return (*this)(0,0); }
    DType det_2_by_2() const { return (*this)(0,0) * (*this)(1,1) - (*this)(0,1) * (*this)(1,0); }
    DType det_3_by_3() const;

	
	friend std::ostream& operator << (std::ostream& os , const Mat& mat) {
		os << "Mat([";
		for(T_DIM r=0; r<dim_row; r++) {
            if (r > 0) os << "     ";
            os << "[";
            for(T_DIM c=0; c<dim_col; c++) {
                os << mat(r, c);
                if (c < dim_col-1) os << ", ";
            }
            os << "]";
            if (r < dim_row-1) os << "\n";
		}
		return os << "])";
	}
};

template<typename DType, T_DIM dim_row, T_DIM dim_col>
Mat<DType, dim_col, dim_row> Mat<DType, dim_row, dim_col>::transpose() const {
	Mat<DType, dim_col, dim_row> ret;
	for (T_IDX i = 0; i < dim_row; i++) {
		for (T_IDX j = 0; j < dim_col; j++) {
			ret(j,i) = (*this)(i,j);
		}
	}
	return ret;
}


template<typename DType, T_DIM dim_row, T_DIM dim_col> template<T_DIM dim_col2>
Mat<DType, dim_row, dim_col2> Mat<DType, dim_row, dim_col>::matmul(const Mat<DType, dim_col, dim_col2>& mat) const {
    Mat<DType, dim_row, dim_col2> ret;
    for (T_DIM r = 0; r < dim_row; r++)
        for (T_DIM c = 0; c < dim_col2; c++)
            for (T_DIM k = 0; k < dim_col; k++) {
                ret(r, c) += (*this)(r, k) * mat(k, c);
    }
    return ret;
}

template<typename DType, T_DIM dim_row, T_DIM dim_col>
Vec<DType, dim_row> Mat<DType, dim_row, dim_col>::matmul(const Vec<DType, dim_col>& vec) const {
    Vec<DType, dim_row> ret;
    for (T_DIM d = 0; d < dim_row; d++) {
        for (T_DIM k = 0; k < dim_col; k++)
            ret[d] += (*this)(d, k) * vec[k];
    }
    return ret;
}


template<typename DType, T_DIM dim_row, T_DIM dim_col>
DType Mat<DType, dim_row, dim_col>::det() const {
	/// @note supports only less than 3-by-3 matrices
	static_assert(dim_row == dim_col);
	static_assert(dim_row <= 3);
	if (dim_row == dim_col && dim_row == 1) return det_1_by_1();
	if (dim_row == dim_col && dim_row == 2) return det_2_by_2();
	if (dim_row == dim_col && dim_row == 3) return det_3_by_3();
}

template<typename DType, T_DIM dim_row, T_DIM dim_col>
DType Mat<DType, dim_row, dim_col>::det_3_by_3() const {
	return (*this)(0,0) * ((*this)(1,1) * (*this)(2,2) - (*this)(1,2) * (*this)(2,1))
 		 - (*this)(0,1) * ((*this)(1,0) * (*this)(2,2) - (*this)(1,2) * (*this)(2,0))
 		 + (*this)(0,2) * ((*this)(1,0) * (*this)(2,1) - (*this)(1,1) * (*this)(2,0));
}


#endif