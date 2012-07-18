#include <cmath>
#include <ctime>
#include <iostream>
#include <queue>
#include <set>

#include "lattice.h"

namespace g_of_r {

LatticeVector::LatticeVector(LatFloat x, LatFloat y, LatFloat z) {
  (*this)(0) = x;
  (*this)(1) = y;
  (*this)(2) = z;
}

LatticeVector::LatticeVector() {
  this->clear();
}

LatticeVector::LatticeVector(const LatticeCell &v) {
  for (int i = 0; i < v.size(); ++i) {
    (*this)(i) = v(i);
  }
}

bool LatticeVector::operator<(const LatticeVector& v) const {
  for (unsigned i = 0; i < v.size(); ++i) {
    if (kEpsilon + (*this)(i) < v(i))
      return true;
    if (kEpsilon + v(i) < (*this)(i))
      return false;
  }
  return false;
}

LatticeCell::LatticeCell() {
  this->clear();
}

LatticeCell::LatticeCell(int x, int y, int z) {
  (*this)(0) = x;
  (*this)(1) = y;
  (*this)(2) = z;
}

bool LatticeCell::operator<(const LatticeCell& v) const {
  for (unsigned i = 0; i < v.size(); ++i) {
    if ((*this)(i) < v(i))
      return true;
    if ((*this)(i) > v(i))
      return false;
  }
  return false;
}

LatticeCell::LatticeCell(const LatticeCell &v) {
  for (int i = 0; i < v.size(); ++i) {
    (*this)(i) = static_cast<int>((v(i) > 0 ? 0.5 : -0.5) + v(i));
  }
}

LatFloat Metric::NormSquared(const LatticeVector &v) const {
  return inner_prod(v, prod(*this, v));
}

Metric LatticeVectorSizeComparator::g = MetricCubic(1);

UnitCell::UnitCell(): atoms_() {}

LatFloat Metric::CellVolume() const {
  const Metric &a = *this;
  double
      cofactor1 = a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)),
      cofactor2 = a(0, 1) * (a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0)),
      cofactor3 = a(0, 2) * (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0));
  return static_cast<LatFloat>(sqrt(cofactor1 - cofactor2 + cofactor3));
}

void UnitCell::AddAtom(int atomic_number, LatticeVector v, LatFloat sigma) {
  atoms_.push_back(Atom(atomic_number, v, sigma));
}

LatFloat UnitCell::NumberDensity(const Metric &g) const {
  return NumAtoms() / g.CellVolume();
}

ublas::vector<LatFloat> Structure::rho(const ublas::vector<LatFloat> &r,
                                   const Metric &g) const {
  const Structure &me = *this;  // It's convenient to avoid dereferencing
  ublas::vector<LatFloat> rho(r);
  rho.clear();
  unsigned N = NumAtoms();
  LatFloat factor = 1.0 / (2 * pi * N * root_two_pi);
  AtomDistancePQ distances;
  for (unsigned i = 0; i < (N - 1); ++i) {
    LatticeVector v_i(me(i).location());
    LatFloat sigma_i_sq = me(i).sigma() * me(i).sigma();
    for (unsigned j = i + 1; j < N; ++j) {
      LatticeVector v_j(me(j).location());
      v_j -= v_i;
      LatFloat r_ij = sqrt(g.NormSquared(v_j));
      LatFloat sigma = sqrt(sigma_i_sq + (me(j).sigma() * me(j).sigma()));
      distances.push(AtomDistance(r_ij, sigma));
    }
  }
  // Now, pop atoms off the priority queue to populate rho(r).
  // One particular interatomic distance (pair; first = dist., second = unc.):
  AtomDistance distance;
  // The smallest and largest r which a given atom affects:
  LatFloat atom_r_min, atom_r_max;
  // The smallest index we have to bother computing anything for:
  int i_min = 0;
  // Convenient variable names
  double r12, sigma12;
  while (!distances.empty()) {
    distance = distances.top();
    distances.pop();
    r12 = distance.first;
    sigma12 = distance.second;
    atom_r_min = r12 - kSigmaCutoff * sigma12;
    atom_r_max = r12 + kSigmaCutoff * sigma12;
    while (r[i_min] < atom_r_min && i_min < r.size())
      ++i_min;
    for (int i = i_min; i < r.size(); ++i) {
      if (r[i] > atom_r_max)
        break;
      // Precomputing 1/sqrt(2*pi*sigma^2) for each atomic distance could be nice!
      double dr = (r[i] - r12) / sigma12;
      rho[i] += factor * (r[i] / r12) * exp(-0.5 * dr * dr) / sigma12;
    }
  }
  return rho;
}

ublas::vector<LatFloat> Structure::S(const ublas::vector<LatFloat> &q,
                                   const Metric &g) const {
  //time_t time_curr, time_prev;
  unsigned N = NumAtoms();
  const Structure &me = *this;  // It's convenient to avoid dereferencing
  ublas::vector<LatFloat> s_vals(q.size());
  s_vals.clear();
  LatFloat factor = 2.0 / N;
  for (unsigned i = 0; i < (N - 1); ++i) {
    //time_prev = time(NULL);
    //std::cerr << "Processing atom " << i << "...";
    //std::cerr.flush();
    LatticeVector vi(me(i).location());
    LatFloat sigma_i_sq = me(i).sigma() * me(i).sigma();
    for (unsigned j = i + 1; j < N; ++j) {
      //std::cerr << me(i).location() << " <-> " << me(j).location() << std::endl;
      LatticeVector vj(me(j).location());
      vj -= vi;
      LatFloat r = sqrt(g.NormSquared(vj));
      LatFloat sigma_sq = sigma_i_sq + (me(j).sigma() * me(j).sigma());
      for (unsigned k = 0; k < q.size(); ++k) {
        s_vals[k] += factor * (sin(q[k] * r) / (q[k] * r)) *
            exp(-0.5 * q[k] * q[k] * sigma_sq);
      }
    }
    //time_curr = time(NULL);
    //std::cerr << "done!  It took " << (time_curr - time_prev) << "seconds."
    //    << std::endl;
  }
  for (unsigned k = 0; k < q.size(); ++k) {
    s_vals[k] += + 1;
  }
  return s_vals;
}

Sphere::Sphere(LatFloat radius, const UnitCell &uc, const Metric &g) {
  typedef std::set<LatticeCell> CellSet;
  CellSet discovered_cells;
  CellSet::iterator it;
  std::queue<LatticeCell> cells_to_process;
  LatticeCell v(0, 0, 0);

  // Start out with the central cell.  Then as long as the queue isn't empty, we
  // search the closest cell for atoms which are close enough.  If we find any,
  // add any previously-undiscovered neighbor cells to the queue.
  cells_to_process.push(v);
  discovered_cells.insert(v);
  while (!cells_to_process.empty()) {
    v = cells_to_process.front();
    cells_to_process.pop();
    // Search the unit cell for any atoms which might be within r_max, and add
    // them to the priority queue of distances.
    bool found_any = false;
    LatFloat u_mag;
    for (int a = 0; a < uc.NumAtoms(); ++a) {
      LatticeVector u(v);
      u += uc(a).location();
      u_mag = sqrt(g.NormSquared(u));
      if (u_mag < radius) {
        found_any = true;
        Atom a_new(uc(a));
        a_new.set_location(u);
        atoms_.push_back(a_new);
      }
    }
    if (!found_any)
      continue;
    // Add undiscovered neighboring lattice cells
    for (unsigned i = 0; i < v.size(); ++i) {
      for (int change = -1; change <= 1; change += 2) {
        LatticeCell v_new(v);
        v_new(i) += change;
        if (discovered_cells.find(v_new) == discovered_cells.end()) {
          cells_to_process.push(v_new);
          discovered_cells.insert(v_new);
        }
      }
    }
  }
}

} //namespace g_of_r
