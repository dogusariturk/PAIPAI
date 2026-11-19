// Regression test: interstitial sites must follow a variable-cell relaxation.
//
// A 6x6x6 bcc Nb cell with tetrahedral interstitial sites is "relaxed" by stretching the cell
// along x by about 4.6%, with every atom at an
// unchanged fractional position. Then:
//   1. updateCoordinatesFromContcar() must place every site at its stretched reference position,
//   2. reconcileInterstitialSitesFromContcar() must keep all interstitials on their sites
//      (the finiteT check), and
//   3. seedTrialCoordinatesFromContcar() must seed a newly occupied site at its stretched position.
// Before the fix, sites stayed at their unstretched Cartesian positions (0.87 A off at the far
// edge here), so checks 1 and 3 fail. Check 2 still passes on this ideal lattice, but in real
// relaxed structures local distortions add to the misplacement, relaxed interstitials land
// nearer to a neighbouring site, and the finiteT trial is rejected as a hop.
#include "structure.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

namespace {

const int N = 6;
const double A = 3.32;
const double L = A * N;
const double STRETCH = 20.83 / 19.92;
const int N_B = 8;
const double TOL = 1e-3;   // outputvasp writes 6 significant digits

const char* SAVE_FILE = "test_cell_rescale_SAVE";
const char* CONTCAR_FILE = "test_cell_rescale_CONTCAR";
const char* MAP_FILE = "test_cell_rescale_map.dat";

std::vector<std::vector<double>> lattice(const std::vector<std::vector<double>>& basis){
    std::vector<std::vector<double>> pts;
    for (const auto& b : basis)
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++)
                    pts.push_back({(i + b[0]) * A, (j + b[1]) * A, (k + b[2]) * A});
    return pts;
}

// Reference SAVE with the B atoms on the sites farthest along x, where the misplacement is largest.
std::vector<int> write_reference_save(){
    auto metals = lattice({{0, 0, 0}, {0.5, 0.5, 0.5}});
    auto sites = lattice({{0.5, 0.25, 0}, {0.5, 0.75, 0}, {0.25, 0.5, 0}, {0.75, 0.5, 0},
                          {0.5, 0, 0.25}, {0.5, 0, 0.75}, {0.25, 0, 0.5}, {0.75, 0, 0.5},
                          {0, 0.5, 0.25}, {0, 0.5, 0.75}, {0, 0.25, 0.5}, {0, 0.75, 0.5}});
    std::vector<int> order(sites.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b){ return sites[a][0] > sites[b][0]; });
    std::vector<int> occupied(order.begin(), order.begin() + N_B);
    std::sort(occupied.begin(), occupied.end());

    std::ofstream f(SAVE_FILE);
    f << "Nb + B\n1.0\n" << L << " 0 0\n0 " << L << " 0\n0 0 " << L << "\n";
    f << "Nb\n" << metals.size() << "\nB\n" << N_B << "\n" << sites.size() << "\nNo Shuffle\nCartesian\n";
    for (const auto& p : metals) f << p[0] << " " << p[1] << " " << p[2] << " 0\n";
    for (int s = 0; s < (int)sites.size(); s++){
        bool occ = std::binary_search(occupied.begin(), occupied.end(), s);
        f << sites[s][0] << " " << sites[s][1] << " " << sites[s][2] << " " << (occ ? 0 : -1) << "\n";
    }
    return occupied;
}

// Distance between a site position and the stretched reference site, with periodic images.
double stretched_error(const std::vector<double>& got, const std::vector<double>& ref){
    const double box[3] = {L * STRETCH, L, L};
    const double want[3] = {ref[0] * STRETCH, ref[1], ref[2]};
    double d2 = 0.0;
    for (int k = 0; k < 3; k++){
        double d = got[k] - want[k];
        d -= box[k] * std::round(d / box[k]);
        d2 += d * d;
    }
    return std::sqrt(d2);
}

bool report(const char* name, bool ok, const std::string& detail){
    std::cerr << (ok ? "PASS " : "FAIL ") << name << ": " << detail << "\n";
    return ok;
}

}  // namespace

int main(){
    std::ostringstream quiet;               // PAIPAI logs to stdout; keep the test output readable
    std::streambuf* keep = std::cout.rdbuf(quiet.rdbuf());

    std::vector<int> occupied = write_reference_save();
    Structure ref;
    if (!ref.readstruc(SAVE_FILE)){ std::cout.rdbuf(keep); std::cerr << "cannot read " << SAVE_FILE << "\n"; return 2; }
    ref.buildIntsiteMetalNeighborMap(2.5);
    ref.outputIntsiteMetalNeighborMap(MAP_FILE);

    // "Relaxed" structure: cell stretched along x, all fractional coordinates unchanged.
    Structure relaxed = ref;
    relaxed.cell_x1 *= STRETCH;
    for (auto& p : relaxed.pos) p[0] *= STRETCH;
    for (auto& p : relaxed.interstitial_pos) p[0] *= STRETCH;
    relaxed.outputvasp(CONTCAR_FILE);

    bool all_ok = true;

    Structure cage = ref;
    int updated = cage.updateCoordinatesFromContcar(CONTCAR_FILE, MAP_FILE);
    double max_err = 0.0;
    for (int s = 0; s < ref.num_interstitial; s++)
        max_err = std::max(max_err, stretched_error(cage.interstitial_pos[s], ref.interstitial_pos[s]));
    std::cout.rdbuf(keep);
    all_ok &= report("site positions follow the relaxed cell", updated && max_err < TOL,
                     "max site error " + std::to_string(max_err) + " A");

    std::cout.rdbuf(quiet.rdbuf());
    Structure check = ref;
    int n_reassigned = 0;
    int consistent = check.reconcileInterstitialSitesFromContcar(CONTCAR_FILE, MAP_FILE, 0, 1.5, "", n_reassigned);
    std::cout.rdbuf(keep);
    all_ok &= report("finiteT site check keeps every interstitial on its site", consistent == 1,
                     consistent ? "no hop" : "hop reported (see message above)");

    // Move one B to an empty site near the far edge and seed the trial from the relaxed CONTCAR.
    std::cout.rdbuf(quiet.rdbuf());
    int target = -1;
    for (int s = ref.num_interstitial - 1; s >= 0 && target < 0; s--)
        if (ref.interstitial_postype[s] < 0 && ref.interstitial_pos[s][0] > L - A) target = s;
    Structure trial = ref;
    trial.swapInterstitial(occupied[0], target);
    int seeded = trial.seedTrialCoordinatesFromContcar(ref, CONTCAR_FILE, MAP_FILE);
    double seed_err = stretched_error(trial.interstitial_pos[target], ref.interstitial_pos[target]);
    std::cout.rdbuf(keep);
    all_ok &= report("newly occupied site is seeded in the relaxed cell", seeded && seed_err < TOL,
                     "seed error " + std::to_string(seed_err) + " A");

    std::remove(SAVE_FILE);
    std::remove(CONTCAR_FILE);
    std::remove(MAP_FILE);
    return all_ok ? 0 : 1;
}
