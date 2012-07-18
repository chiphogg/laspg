#ifndef G_OF_R_LATTICE_H_
#define G_OF_R_LATTICE_H_

#include <set>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

namespace g_of_r {

// I'm unsure whether to use float or double.  float is not really any faster,
// so it's strictly a storage question.  Do I want 7 significant figures or 15?
// This typedef makes it easy for me to change my mind.
//
// NOTE: when calculating g(r) for bulk Au out to 50+ nm, using 'float' causes
// noticeable noise, whereas using 'double' does not.
typedef float LatFloat;
//typedef double LatFloat;

namespace ublas = boost::numeric::ublas;

// The dimensionality of the vector space.  NOTE: simply changing this (say,
// from 3 to 2) will NOT produce a working 2D vector space!  Certain functions
// have parameter lists which assume a particular dimensionality!  I defined
// this to reduce the instances of arbitrary "3"'s popping up in the code.
const int kDim = 3;

const double root_two_pi = boost::math::constants::root_two_pi<double>();
const double pi = boost::math::constants::pi<double>();

// Anything more than kSigmaCutoff sigma's from the mean is considered to be 0.
const LatFloat kSigmaCutoff = 5;

/*******************************************************************************
* CLASS LatticeCell                                                            *
*******************************************************************************/

class LatticeCell : public ublas::bounded_vector<int, kDim> {
 public:
  LatticeCell();
  LatticeCell(int x, int y, int z);
  LatticeCell(const LatticeCell &v);
  bool operator<(const LatticeCell& v) const;
};

/*******************************************************************************
* CLASS LatticeVector                                                          *
*******************************************************************************/

class LatticeVector : public ublas::bounded_vector<LatFloat, kDim> {
 public:
  LatticeVector();
  LatticeVector(const LatticeCell &v);
  LatticeVector(LatFloat x, LatFloat y, LatFloat z);
  bool operator<(const LatticeVector& v) const;

 private:
  // This is used for floating point comparison's sake.  Any lattice vector
  // components which differ by less than kEpsilon (as a fraction of the unit
  // cell) are pretty well equivalent.
  static const LatFloat kEpsilon = 1e-8;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* types for atom displacements                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef std::pair<LatFloat, LatFloat> AtomDistance;
typedef std::pair<LatticeVector, LatFloat> AtomDisplacement;

// One AtomDistance is "smaller" than another when it contributes non-negligibly
// to a smaller r.  (Note the reversed comparison order, since our PQ wants the
// *smallest*, not the biggest.)
typedef struct {
  bool operator()(const AtomDistance &a, const AtomDistance &b) const {
    LatFloat a_smallest_r = (a.first - kSigmaCutoff * a.second);
    LatFloat b_smallest_r = (b.first - kSigmaCutoff * b.second);
    return a_smallest_r > b_smallest_r;
  }
} AtomDistanceComparator;

typedef std::priority_queue<AtomDistance, std::vector<AtomDistance>,
        AtomDistanceComparator> AtomDistancePQ;

/*******************************************************************************
* CLASS Metric                                                                 *
*******************************************************************************/

// Defines a 3-dimensional metric tensor.  Sample usage:
//    Metric g(1.0, 1.0, 2.0);
class Metric : public ublas::bounded_matrix<LatFloat, kDim, kDim> {
 public:
  Metric(LatFloat g11, LatFloat g22, LatFloat g33, 
         LatFloat g12 = 0.0, LatFloat g13 = 0.0, LatFloat g23 = 0.0)
  {
    (*this)(0, 0) = g11;
    (*this)(1, 1) = g22;
    (*this)(2, 2) = g33;
    (*this)(0, 1) = (*this)(1, 0) = g12;
    (*this)(0, 2) = (*this)(2, 0) = g13;
    (*this)(1, 2) = (*this)(2, 1) = g23;
  }

  // Computes the squared norm of v using this metric.  Sample usage:
  //    Metric g(1.0, 1.0, 2.0);
  //    LatticeVector v(1.0, 1.0, 1.0);
  //    LatFloat norm_squared = g.NormSquared(v); // 4.0
  LatFloat NormSquared(const LatticeVector &v) const;

  LatFloat CellVolume() const;
};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* SUBCLASS MetricCubic                                                         *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

class MetricCubic : public Metric {
 public:
  MetricCubic(LatFloat a): Metric(a * a, a * a, a * a) {}
};

/*******************************************************************************
* CLASS LatticeVectorSizeComparator                                                *
*******************************************************************************/

// Compares two LatticeVector objects, given a Metric, to see which vector is
// larger.  (This lets us make a priority queue of LatticeVector objects,
// smallest first.)
// Sample usage:
//    LatticeVector v1(3, 2, 2);
//    LatticeVector v2(2, 2, 2);
//    LatticeVectorSizeComparator::g = MetricCubic(1);
//    std::priority_queue<LatticeVector, std::vector<LatticeVector>,
//      LatticeVectorSizeComparator> pq;
//    pq.push(v1);
//    pq.push(v2);
//    assert(pq.top() == v2)
class LatticeVectorSizeComparator {
 public:
  static Metric g;
  bool operator() (const LatticeVector &v1, const LatticeVector &v2) const {
    return g.NormSquared(v1) > g.NormSquared(v2);
  }

};

/*******************************************************************************
* CLASS Atom                                                                   *
*******************************************************************************/

const int kNotAnAtom = -1;

class Atom {
 public:
  Atom(int a, LatticeVector v, LatFloat sigma):
      atomic_number_(a), location_(v), sigma_(sigma) {}
  Atom(): atomic_number_(kNotAnAtom), location_(), sigma_(0.0) {}
  Atom(const Atom &a) : atomic_number_(a.atomic_number_),
    location_(a.location_), sigma_(a.sigma_) {}

  int atomic_number() { return atomic_number_; }

  LatticeVector location() { return location_; }
  void set_location(const LatticeVector &v) { location_ = v; }

  LatFloat sigma() { return sigma_; }

 private:

  // The element of the periodic table.
  int atomic_number_;
  // The location within the unit cell.
  LatticeVector location_;
  // The standard deviation of this atom's position (assumed to be isotropic
  // Gaussian).
  LatFloat sigma_;
};

/*******************************************************************************
* CLASS UnitCell                                                               *
*******************************************************************************/

class UnitCell {
 public:
  UnitCell();

  // Adds an atom (of the type specified by "atomic_number") to this unit cell,
  // at lattice vector "v".
  void AddAtom(int atomic_number, LatticeVector v, LatFloat sigma);

  // Returns the number of atoms in the unit cell.
  int NumAtoms() const {
    return atoms_.size();
  }

  Atom operator()(unsigned i) const {
    return atoms_[i];
  }

  LatFloat NumberDensity(const Metric &g) const;

 private:
  std::vector<Atom> atoms_;
};

/*******************************************************************************
* CLASS Structure                                                              *
*******************************************************************************/

class Structure {
 public:
  // Returns the number of atoms in the unit cell.
  int NumAtoms() const {
    return atoms_.size();
  }

  // Retrieves the i'th atom in the structure (ordering is arbitrary).
  Atom operator()(unsigned i) const {
    return atoms_[i];
  }

  // Compute rho(r) for a given structure using the definition.
  ublas::vector<LatFloat> rho(const ublas::vector<LatFloat> &r, const Metric &g)
      const;

  // Compute S(Q) for a given structure using the Debye equation.
  ublas::vector<LatFloat> S(const ublas::vector<LatFloat> &q, const Metric &g)
      const;

 protected:
  std::vector<Atom> atoms_;

};

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* SUBCLASS Sphere                                                              *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

class Sphere : public Structure {
 public:
  Sphere(LatFloat radius, const UnitCell &uc, const Metric &g);
};

} // namespace g_of_r

#endif // G_OF_R_LATTICE_H_
