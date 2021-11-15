#ifndef GAIA_SMALL_VECTORS_HPP
#define GAIA_SMALL_VECTORS_HPP

#include <ostream>
#include <cmath>
#include <omp.h>

namespace gaia{

  /*
    This class code is based on the vector class
    in exafmm (https://github.com/exafmm/exafmm/blob/master/include/vec.h)
  */
  template<int N, typename T>
  class vec {
  private:
    T* data;
  public:
    vec(){data = new T[N];}
    vec(const T &v){ data = new T[N]; for (int i=0; i<N; i++) data[i] = v;}
    vec(const vec &v) { data = new T[N]; for(int i=0; i<N; i++) data[i] = v[i];}
    ~vec(){delete[] data;}
    const vec &operator=(const T v) {for (int i=0; i<N; i++) data[i] = v; return *this;}
    const vec &operator+=(const T v) {for (int i=0; i<N; i++) data[i] += v;return *this;}
    const vec &operator-=(const T v) {for (int i=0; i<N; i++) data[i] -= v;return *this;}
    const vec &operator*=(const T v) {for (int i=0; i<N; i++) data[i] *= v;return *this;}
    const vec &operator/=(const T v) {for (int i=0; i<N; i++) data[i] /= v;return *this;}
    const vec &operator>=(const T v) {for (int i=0; i<N; i++) data[i] >= v;return *this;}
    const vec &operator<=(const T v) {for (int i=0; i<N; i++) data[i] <= v;return *this;}
    const vec &operator&=(const T v) {for (int i=0; i<N; i++) data[i] &= v;return *this;}
    const vec &operator|=(const T v) {for (int i=0; i<N; i++) data[i] |= v;return *this;}
    const vec &operator=(const vec & v) {for (int i=0; i<N; i++) data[i] = v[i];return *this;}
    const vec &operator+=(const vec & v) {for (int i=0; i<N; i++) data[i] += v[i];return *this;}
    const vec &operator-=(const vec & v) {for (int i=0; i<N; i++) data[i] -= v[i];return *this;}
    const vec &operator*=(const vec & v) {for (int i=0; i<N; i++) data[i] *= v[i];return *this;}
    const vec &operator/=(const vec & v) {for (int i=0; i<N; i++) data[i] /= v[i];return *this;}
    const vec &operator>=(const vec & v) {for (int i=0; i<N; i++) data[i] >= v[i];return *this;}
    const vec &operator<=(const vec & v) {for (int i=0; i<N; i++) data[i] <= v[i];return *this;}
    const vec &operator&=(const vec & v) {for (int i=0; i<N; i++) data[i] &= v[i];return *this;}
    const vec &operator|=(const vec & v) {for (int i=0; i<N; i++) data[i] |= v[i];return *this;}
    vec operator+(const T v) const {return vec(*this) += v;}
    vec operator-(const T v) const {return vec(*this) -= v;}
    vec operator*(const T v) const {return vec(*this) *= v;}
    vec operator/(const T v) const {return vec(*this) /= v;}
    vec operator>(const T v) const {return vec(*this) >= v;}
    vec operator<(const T v) const {return vec(*this) <= v;}
    vec operator&(const T v) const {return vec(*this) &= v;}
    vec operator|(const T v) const {return vec(*this) |= v;}
    vec operator+(const vec & v) const {return vec(*this) += v;}
    vec operator-(const vec & v) const {return vec(*this) -= v;}
    vec operator*(const vec & v) const {return vec(*this) *= v;}
    vec operator/(const vec & v) const {return vec(*this) /= v;}
    vec operator>(const vec & v) const {return vec(*this) >= v;}
    vec operator<(const vec & v) const {return vec(*this) <= v;}
    vec operator&(const vec & v) const {return vec(*this) &= v;}
    vec operator|(const vec & v) const {return vec(*this) |= v;}
    vec operator-() const {vec temp;for (int i=0; i<N; i++) temp[i] = -data[i];return temp;}
    T &operator[](int i) {return data[i];}
    const T &operator[](int i) const {return data[i];}
    operator       T* ()       {return data;}
    operator const T* () const {return data;}
    friend T* K_(vec& v){return v.data;}
    friend T* K_(vec* v){return v->data;}
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {for (int i=0; i<N; i++) s << v[i] << ' ';return s;}
    friend T sum(const vec & v) {T temp = 0;for (int i=0; i<N; i++) temp += v[i];return temp;}
    friend T norm2(const vec & v) {T temp = 0;for (int i=0; i<N; i++) temp += v[i] * v[i];return sqrt(temp);}
    friend T norm22(const vec & v) {T temp = 0;for (int i=0; i<N; i++) temp += v[i] * v[i];return temp;}
    friend vec min(const vec & v, const vec & w) {vec temp;for (int i=0; i<N; i++) temp[i] = v[i] < w[i] ? v[i] : w[i];return temp;}
    friend vec max(const vec & v, const vec & w) {vec temp;for (int i=0; i<N; i++) temp[i] = v[i] > w[i] ? v[i] : w[i];return temp;}
    friend T min(const vec & v) {T temp = v[0];for (int i=1; i<N; i++) temp = temp < v[i] ? temp : v[i];return temp;}
    friend T max(const vec & v) {T temp = v[0];for (int i=1; i<N; i++) temp = temp > v[i] ? temp : v[i];return temp;}
    friend inline T norm2diff(vec& x, vec& y){
      T res = 0.;
      for(int ii=0; ii<N; ii++){res += (x[ii]-y[ii])*(x[ii]-y[ii]);}
      return sqrt(res);
    }
    friend inline T doprod(const vec<N,T>& x, const vec<N,T>& y){
      T res = 0.;
      for(int ii=0; ii<N; ii++){res += x[ii]*y[ii];}
      return res;
    }
    friend inline T scaleddoprod(const vec& x, const vec& y){return doprod(x,y)/(norm2(x)*norm2(y));}
    friend vec abs(const vec& v){vec temp; for(int i=0;i<N;i++){temp[i] = ((v[i] >= 0.) ? v[i] : - v[i]);}return temp; }
  };

} // GAIA

#endif
