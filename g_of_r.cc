#include <iostream>
#include <queue>
#include <vector>
#include <algorithm>

#include "lattice.h"

namespace ublas = boost::numeric::ublas;

using std::cout;
using std::cerr;
using std::endl;
using ublas::vector;
using ublas::matrix;

namespace g_of_r {

// Checks the unit cell to see if any atoms are closer than r_max to the origin.
// If they are, they get added to the priority queue of atomic distances.
//
// Args:
// INPUT:
//    dv:         Displacement vector to this unit cell.
//    uc:         Unit cell of the crystal.
//    sigma:      Uncertainty in the interatomic distance to each atom.
//    g:          Metric tensor for the crystal.
//    r_max:      Maximum distance we care about.
// OUTPUT:
//    distances:  Priority queue holding the interatomic distances.
bool AnyAtomsCloseEnough(const LatticeVector &dv, const UnitCell &uc, 
                         const std::vector<LatFloat> &sigma, const Metric &g,
                         LatFloat r_max,
                         AtomDistancePQ &distances) {
  bool found_any = false;
  LatFloat u_mag;
  for (int a = 0; a < uc.NumAtoms(); ++a) {
    LatticeVector u(dv);
    u += uc(a).location();
    u_mag = sqrt(g.NormSquared(u));
    if (u_mag - kSigmaCutoff * sigma[a] < r_max) {
      found_any = true;
      if (u_mag > sigma[a]) {  // Don't add the self-distance!
        distances.push(AtomDistance(u_mag, sigma[a]));
      }
    }
  }
  return found_any;
}

// Tabulates a list of atoms with non-negligible probability to be within r_max
// of v0.
//
// Args:
// INPUT:
//    v0:         The origin.
//    sigma_0:    The uncertainty in the origin.
//    r_max:      The maximum distance we care about.
//    uc:         The unit cell of the crystal.
//    g:          The metric tensor for the crystal.
// OUTPUT:
//    distances:  The priority queue holding the interatomic distances.
void ComputeAtomicDistances(LatticeVector v0, LatFloat sigma_0, LatFloat r_max,
                            const UnitCell &uc, const Metric &g,
                            AtomDistancePQ &distances) {
  typedef std::set<LatticeCell> CellSet;
  CellSet discovered_cells;
  CellSet::iterator it;
  std::queue<LatticeCell> cells_to_process;
  LatticeCell v; 

  // Precompute uncertainties in interatomic distances:
  std::vector<LatFloat> sigma(uc.NumAtoms(), sigma_0 * sigma_0);
  for (int a = 0; a < uc.NumAtoms(); ++a)
    sigma[a] = sqrt(sigma[a] + uc(a).sigma() * uc(a).sigma());

  // Start out with the central cell.  Then as long as the queue isn't empty, we
  // search the closest cell for atoms which are close enough.  If we find any,
  // add any previously-undiscovered neighbor cells to the queue.
  cells_to_process.push(v);
  discovered_cells.insert(v);
  while (!cells_to_process.empty()) {
    v = cells_to_process.front();
    cells_to_process.pop();
    LatticeVector dv(v);
    dv -= v0;
    // Search the unit cell for any atoms which might be within r_max, and add
    // them to the priority queue of distances.
    if (!AnyAtomsCloseEnough(dv, uc, sigma, g, r_max, distances))
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

// Calculates the partial rho(r) -- i.e., the contribution to rho(r) starting
// from a given atom.
//
// Args:
//    v0:  The position of the atom (the "origin")
//    sigma_0:  The standard deviation of the origin-atom's position
//    r:  The r-values where we evaluate g(r)
//    uc:  The unit cell of the crystal
//    g:  The metric tensor for the crystal
vector<LatFloat> RhoPartial(LatticeVector v0, LatFloat sigma_0, 
                          const vector<LatFloat> &r, const UnitCell &uc,
                          const Metric &g) {
  vector<LatFloat> rho(r.size());
  rho.clear();
  AtomDistancePQ distances;
  // the maximum distance we care about:
  LatFloat r_max = *(--r.end());
  cerr << "Computing atomic distances...";
  ComputeAtomicDistances(v0, sigma_0, r_max, uc, g, distances);
  cerr << "done! " << distances.size() << " atoms tallied." << endl;

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
      rho[i] += (r[i] / r12) * exp(-0.5 * dr * dr) / 
          (root_two_pi * sigma12);
    }
  }

  double norm_factor = 4.0 * pi;
  for (int i = 0; i < r.size(); ++i) {
    rho[i] /= (norm_factor * r[i] * r[i]);
  }
  return rho;
}

//vector<LatFloat> Rho(const vector<LatFloat> &r, const UnitCell &uc, const Metric &g, LatFloat r_max)
//{
//  vector<LatFloat> rho(r.size());
//  rho.clear();
//
//}

} // namespace g_of_r

using g_of_r::LatFloat;
using g_of_r::LatticeVector;
using g_of_r::pi;

ublas::vector<LatFloat> Sequence(LatFloat from, LatFloat to, int length) {
  ublas::vector<LatFloat> v(length);
  int num_intervals = length - 1;
  LatFloat interval_width = (to - from) / num_intervals;
  for (unsigned i = 0; i < length; i++)
    v(i) = from + i * interval_width;
  return v;
}

g_of_r::Metric Au_g() {
  LatFloat Au_const_nm = 0.40782;
  return g_of_r::MetricCubic(Au_const_nm);
}

g_of_r::UnitCell Au_uc(LatFloat sigma) {
  g_of_r::UnitCell uc;
  int atom_Au = 79;
  uc.AddAtom(atom_Au, LatticeVector(0.0, 0.0, 0.0), sigma);
  uc.AddAtom(atom_Au, LatticeVector(0.0, 0.5, 0.5), sigma);
  uc.AddAtom(atom_Au, LatticeVector(0.5, 0.0, 0.5), sigma);
  uc.AddAtom(atom_Au, LatticeVector(0.5, 0.5, 0.0), sigma);
  return uc;
}

vector<LatFloat> Au_G_of_r(vector<LatFloat> r, LatFloat sigma) {
  g_of_r::Metric g = Au_g();
  g_of_r::LatticeVector v0(0, 0, 0);
  g_of_r::UnitCell uc = Au_uc(sigma);
  double rho_0 = uc.NumberDensity(g);
  cerr << "rho_0: " << rho_0 << endl;
  vector<LatFloat> rho = g_of_r::RhoPartial(v0, sigma, r, uc, g), G=rho;
  for (int i = 0; i < rho.size(); ++i)
    G(i) = 4 * pi * r(i) * (rho(i) - rho_0);
  return G;
}

vector<LatFloat> g_from_G(const vector<LatFloat> &Gr, const vector<LatFloat> &r,
                          const g_of_r::Metric &g, const g_of_r::UnitCell &uc) {
  vector<LatFloat> gr(Gr);
  LatFloat factor = uc.NumberDensity(g) * 4 * pi;
  for (int i = 0; i < Gr.size(); ++i) {
    gr[i] = 1 + Gr[i] / (factor * r[i]);
  }
  return gr;
}

void compute_g_G_Au() {
  LatFloat dr = 0.001;    // nm
  LatFloat r_max =  7;    // nm
  LatFloat sigma = 0.01;  // nm
  ublas::vector<LatFloat> r(Sequence(dr, r_max, r_max / dr));
  g_of_r::Metric g = Au_g();
  g_of_r::UnitCell uc = Au_uc(sigma);
  cerr << "Computing G(r)...";
  cerr.flush();
  vector<LatFloat> Gr = Au_G_of_r(r, sigma);
  cerr << "done!" << endl << "Computing g(r)...";
  vector<LatFloat> gr = g_from_G(Gr, r, g, uc);
  cerr << "done!" << endl;

  cout << "r.nm" << "\t" << "G.of.r" << "\t" << "g.of.r" << endl;
  for (int i = 0; i < r.size(); ++i) {
    cout << r[i] << "\t" << Gr[i] << "\t" << gr[i] << endl;
  }
}

void PrintSOfQ() {
  LatFloat dq    =   0.10;  // nm^(-1)
  LatFloat q_max = 300.00;  // nm^(-1)
  ublas::vector<LatFloat> q(Sequence(dq, q_max, q_max / dq));

  LatFloat sigma = 0.01;  // nm
  LatFloat rad   = 1.00;  // nm

  g_of_r::Metric g = Au_g();
  int N = 10;
  LatFloat epsilon = 1e-7;
  ublas::matrix<LatFloat> S_of_q(q.size(), N);
  for (int i = 0; i < N; ++i) {
    g_of_r::UnitCell uc = Au_uc((sigma * i) / N + epsilon);
    g_of_r::Sphere au_2nm(rad, uc, g);
    ublas::vector<LatFloat> S(au_2nm.S(q, g));
    for (unsigned k = 0; k < q.size(); ++k) {
      S_of_q(k, i) = S(k);
    }
  }

  cout << "q.inv.nm";
  for (unsigned i = 0; i < N; ++i) {
    cout << "\t" << "S" << (i + 1);
  }
  for (unsigned k = 0; k < q.size(); ++k) {
    cout << endl << q(k);
    for (unsigned i = 0; i < N; ++i) {
      cout << "\t" << S_of_q(k, i);
    }
  }
  cout << endl;
}

void PrintSphereCoordinates() {
  g_of_r::Metric g = Au_g();
  g_of_r::UnitCell uc = Au_uc(0.01);
  LatFloat rad = 1.0;
  g_of_r::Sphere au_2nm(rad, uc, g);
  cout << "x\ty\tz" << endl;
  for (unsigned i = 0; i < au_2nm.NumAtoms(); ++i) {
    cout << au_2nm(i).location()(0);
    for (unsigned k = 1; k < g_of_r::kDim; ++k) {
      cout << "\t" << au_2nm(i).location()(k);
    }
    cout << endl;
  }
}

int main(int argc, char **argv) {
  // Define r-values where we evaluate rho
  LatFloat dr = 0.001;    // nm
  LatFloat r_max =  2.1;    // nm
  LatFloat sigma = 0.01;  // nm
  ublas::vector<LatFloat> r(Sequence(dr, r_max, r_max / dr));

  // Define the Au nanoparticle, d/2 = 1 nm
  g_of_r::Metric g = Au_g();
  g_of_r::UnitCell uc = Au_uc(0.01);
  LatFloat rad = 1.0;
  g_of_r::Sphere au_2nm(rad, uc, g);

  // Compute and output rho
  ublas::vector<LatFloat> rho(au_2nm.rho(r, g));
  cout << "r.nm\trho" << endl;
  for (unsigned i = 0; i < r.size(); ++i) {
    cout << r(i) << "\t" << rho(i) << endl;
  }
  return 0;
}
