#ifndef GAIA_VECTORS_HPP
#define GAIA_VECTORS_HPP

#include <ostream>
#include <fstream>
#include <string>
#include <cmath>
#include <omp.h>

namespace gaia{

  /*
    This class code is based on the vector class
    in exafmm (https://github.com/exafmm/exafmm/blob/master/include/vec.h)
  */
  template<typename T>
  class Vec {
  private:
  
    T* data; // Vector entries
    int N;   // Vector size
  
  public:
  
    // CONSTRUCTORS
    // Empty constructor
    Vec(){data = nullptr; N = -1;}
    // Standard constructor (dynamic allocation)
    Vec(int _N) : N(_N) {
      data = (T*)std::malloc(sizeof(T)*N);
    }
    // Parsing constructor (allow to consider any array of type T as a vector)
    Vec(int _N, T* _data) : N(_N) {data = _data;}
    // Initialize with given value
    Vec(int _N, const T &v) : N(_N) { data = new T[N]; for (int i=0; i<N; i++) data[i] = v;}
    // Copy an already declared vector
    Vec(int _N, const Vec &v) : N(_N) { data = new T[N]; for(int i=0; i<N; i++) data[i] = v[i];}
    ~Vec(){
      if(data){
	std::free(data);}N=0;
    }
    void allocate(int _N){
      N = _N;
      data = (T*)std::malloc(sizeof(T)*N);
    }
    void alignAllocate(int _N, int _align){
      N = _N;
      if(posix_memalign((void**)&data, _align, N*sizeof(T)) != 0){
	std::cerr << "Allocation error (alignAllocate in gaia_vector.hpp)" << std::endl;
	exit(1);
      }
    }
    void allocate(int _N, const T& val){
      N = _N;
      data = (T*)std::malloc(sizeof(T)*N);
      for (int i=0; i<N; i++) data[i] = val;
    }
    void free(){
      if(data){std::free(data);}
      N = 0;
    }
    void set(int _N, T* _data){N = _N; data = _data;}

    // BASIC OPERATORS
    const Vec &operator=(const T v) {
      for (int i=0; i<N; i++) data[i] = v;
      return *this;}
    const Vec &operator+=(const T v) {
      for (int i=0; i<N; i++) data[i] += v;
      return *this;}
    const Vec &operator-=(const T v) {
      for (int i=0; i<N; i++) data[i] -= v;
      return *this;}
    const Vec &operator*=(const T v) {
      for (int i=0; i<N; i++) data[i] *= v;
      return *this;}
    const Vec &operator/=(const T v) {
      for (int i=0; i<N; i++) data[i] /= v;
      return *this;}
    const Vec &operator=(const Vec & v) {N = v.N; if(data){delete[] data;} data = new T[N]; for (int i=0; i<N; i++) data[i] = v[i]; return *this;}
    const Vec &operator+=(const Vec & v) {for (int i=0; i<N; i++) data[i] += v[i];return *this;}
    const Vec &operator-=(const Vec & v) {for (int i=0; i<N; i++) data[i] -= v[i];return *this;}
    const Vec &operator*=(const Vec & v) {for (int i=0; i<N; i++) data[i] *= v[i];return *this;}
    const Vec &operator/=(const Vec & v) {for (int i=0; i<N; i++) data[i] /= v[i];return *this;}
    Vec operator+(const T v) const {return Vec(*this) += v;}
    Vec operator-(const T v) const {return Vec(*this) -= v;}
    Vec operator*(const T v) const {return Vec(*this) *= v;}
    Vec operator/(const T v) const {return Vec(*this) /= v;}
    Vec operator+(const Vec & v) const {return Vec(*this) += v;}
    Vec operator-(const Vec & v) const {return Vec(*this) -= v;}
    Vec operator*(const Vec & v) const {return Vec(*this) *= v;}
    Vec operator/(const Vec & v) const {return Vec(*this) /= v;}
    Vec operator-() const {Vec temp;for (int i=0; i<N; i++) temp[i] = -data[i];return temp;}
    friend std::ostream &operator<<(std::ostream & s, const Vec & v) {for (int i=0; i<v.N; i++) s << v[i] << ' ';return s;}

    // ACCESSORS
    T &operator[](int i) {return data[i];}
    const T &operator[](int i) const {return data[i];}
    operator       T* ()       {return data;}
    operator const T* () const {return data;}
    friend T* K_(Vec& v){return v.data;}
    friend T* K_(Vec* v){return v->data;}
    friend int Size(const Vec & v){return v.N;}
    friend const int* pSize(const Vec & v){return &v.N;}

    // SIMPLE USUAL FUNCTIONs (that will be discarded in a next version)
    friend T norm2(const Vec & v) {T temp = 0;for (int i=0; i<v.N; i++) temp += v[i] * v[i];return sqrt(temp);}
    friend Vec min(const Vec & v, const Vec & w) {Vec temp;for (int i=0; i<v.N; i++) temp[i] = v[i] < w[i] ? v[i] : w[i];return temp;}
    friend Vec max(const Vec & v, const Vec & w) {Vec temp;for (int i=0; i<v.N; i++) temp[i] = v[i] > w[i] ? v[i] : w[i];return temp;}
    friend T min(const Vec & v) {T temp = v[0];for (int i=1; i<v.N; i++) temp = temp < v[i] ? temp : v[i];return temp;}
    friend T max(const Vec & v) {T temp = v[0];for (int i=1; i<v.N; i++) temp = temp > v[i] ? temp : v[i];return temp;}
    friend inline T doprod(const Vec& x, const Vec& y){
      T res = 0.;
      for(int ii=0; ii<x.N; ii++){res += x[ii]*y[ii];}
      return res;
    }
    friend Vec abs(const Vec& v){Vec temp; for(int i=0;i<v.N;i++){temp[i] = ((v[i] >= 0.) ? v[i] : - v[i]);}return temp; }
    /** Norms **/
    friend double nrmInfRel(Vec& v, Vec& u){
      double tmpa = 0.;
      double tmpb = 0.;
      for(int i = 0; i < Size(v); i++){
	tmpa = std::max(std::abs(v[i]-u[i]),tmpa);
	tmpb = std::max(std::abs(v[i]     ),tmpb);
      }
      return (tmpa / tmpb);
    }
    friend double nrmInfRel(Vec& v, T* u){
      double tmpa = 0.;
      double tmpb = 0.;
      for(int i = 0; i < Size(v); i++){
	tmpa = std::max(std::abs(v[i]-u[i]),tmpa);
	tmpb = std::max(std::abs(v[i]     ),tmpb);
      }
      return (tmpa / tmpb);
    }
    friend double nrm2Rel(Vec& v, Vec& u){
      double tmpa = 0.;
      double tmpb = 0.;
      for(int i = 0; i < Size(v); i++){
	tmpa += std::abs(v[i]-u[i])*std::abs(v[i]-u[i]);
	tmpb += std::abs(v[i])*std::abs(v[i]);
      }
      return std::sqrt(tmpa / tmpb);
    }
    friend double nrm2Rel(Vec& v, T*& u){
      double tmpa = 0.;
      double tmpb = 0.;
      for(int i = 0; i < Size(v); i++){
	tmpa += std::abs(v[i]-u[i])*std::abs(v[i]-u[i]);
	tmpb += std::abs(v[i])*std::abs(v[i]);
      }
      return std::sqrt(tmpa / tmpb);
    }
    friend double nrm1Rel(Vec& v, Vec& u){
      double tmpa = 0.;
      double tmpb = 0.;
      for(int i = 0; i < Size(v); i++){
	tmpa += std::abs(v[i]-u[i]);
	tmpb += std::abs(v[i]);
      }
      return (tmpa / tmpb);
    }
    friend double nrm1Rel(Vec& v, T*& u){
      double tmpa = 0.;
      double tmpb = 0.;
      for(int i = 0; i < Size(v); i++){
	tmpa += std::abs(v[i]-u[i]);
	tmpb += std::abs(v[i]);
      }
      return (tmpa / tmpb);
    }
  };
  
  
} // gaia

#endif
