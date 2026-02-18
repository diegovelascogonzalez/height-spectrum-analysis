// Height spectrum (v4, binned via k-index rings), 09/13/2025
//
// Build:
//   g++ heightspectrum_v4.cc -o heightspectrum -std=c++17 -lchemfiles -lfftw3 -lm
//
// Runtime note:
//   Ensure Chemfiles and FFTW are discoverable at runtime, e.g.:
//     export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
//
// What’s new vs v3:
//  - Optional detrending (per-leaflet mean removal) via REMOVE_MEAN toggle.
//  - Pipeline now bins by Fourier index radius r^2 = kx_eff^2 + ky^2 (r2c layout) per frame,
//    averages within each ring (bin), then computes time-averages & errors on the binned series.
//  - Outputs with suffix *_binned.dat. (Unbinned avg file is still written for continuity.)
//  - Fixed-blocking errors use M=10 consecutive blocks (as even as possible).

#include <chemfiles.hpp>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
#include <cmath>
#include <vector>
#include <tuple>
#include <fftw3.h>
#include <numeric>
#include <map>
#include <unordered_map>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <array>
#include <stdexcept>
#include <functional>
#include <limits>

using namespace chemfiles;

// ----------------------- User-configurable -----------------------
const std::string folder          = "/home/dlvg/Desktop/fixing_heightspec/60/systematic_analysis/oct_100percent/"; // <-- set me
const std::string trajectory_file = folder + "PN.dcd";                      // or PN.dcd
const std::string pdb_file        = folder + "PN.pdb";                      // or PN.pdb

// grid spacing in trajectory length units (Å, nm, etc.)
const double L_grid   = 16.44830556;    
const bool   centered = true;           // set true if coords are centered in [-L/2, L/2)

// Detrending (NEW): subtract per-leaflet spatial mean before FFT if true
const bool REMOVE_MEAN = false;         // toggle on/off detrending

// FP method guardrails (as in v2)
const double FP_REL_EPS     = 0.15;     // relative epsilon for plateau
const size_t FP_MIN_BLOCKS  = 16;       // minimum number of blocks at chosen level
const size_t FP_NEED_CONSEC = 2;        // need 2 consecutive small deltas (span of 3)

// Fixed-blocking (binned time series): use exactly M blocks (even split)
const size_t FIXED_BLOCKS_M = 10;

static constexpr double PI = 3.1415926535897932384626433832795;

// ----------------------- Utilities -----------------------
bool createDirectory(const std::string &path) {
    try {
        if (!std::filesystem::exists(path)) {
            std::filesystem::create_directory(path);
            std::cout << "Directory created: " << path << "\n";
        } else {
            std::cout << "Directory already exists: " << path << "\n";
        }
        return true;
    } catch (const std::exception &e) {
        std::cerr << "Error creating directory: " << e.what() << "\n";
        return false;
    }
}

void trajectory_info(chemfiles::Trajectory &traj, int &num_frames, int &num_atoms,
                     double &Lx, double &Ly, double &Lz) {
    using namespace chemfiles;
    size_t frame_count = 0;
    size_t atom_count  = 0;
    double total_x = 0.0, total_y = 0.0, total_z = 0.0;
    double total_area = 0.0;

    while (true) {
        try {
            Frame frame = traj.read();
            ++frame_count;
            if (frame_count == 1) {
                atom_count = frame.size();
            }
            auto lengths = frame.cell().lengths();
            total_x += lengths[0];
            total_y += lengths[1];
            total_z += lengths[2];
            total_area += lengths[0] * lengths[1];
        } catch (const Error &) { break; }
    }

    if (frame_count > 0) {
        num_frames = static_cast<int>(frame_count);
        num_atoms  = static_cast<int>(atom_count);

        Lx = total_x / frame_count;
        Ly = total_y / frame_count;
        Lz = total_z / frame_count;
        double mean_area = total_area / frame_count;

        std::cout << "Total frames: " << frame_count << "\n";
        std::cout << "Atoms/frame: " << atom_count << "\n";
        std::cout << "Mean box (length units): " << Lx << " × " << Ly << " × " << Lz << "\n";
        std::cout << "Mean XY area/frame: " << mean_area << "\n";

        std::ofstream area_out("A.dat");
        if (area_out) {
            area_out << std::setprecision(10) << std::fixed << mean_area << "\n";
            std::cout << "Saved mean area to A.dat\n";
        } else {
            std::cerr << "Failed to write A.dat\n";
        }
    } else {
        std::cerr << "No frames read from trajectory.\n";
    }
}

inline void check_boundaries(std::array<double,3> &position, const chemfiles::Vector3D &lengths) {
    for (size_t i = 0; i < 3; ++i) {
        position[i] = std::fmod(position[i], lengths[i]);
        if (position[i] < 0) position[i] += lengths[i];
    }
}

inline int wrap_index(int i, int N) { return (i % N + N) % N; }

// ------------------- v2-style weighted hole filling -------------------
void fill_empty_cells(
    std::vector<double> &z_grid,
    std::vector<int>    &z_counts,
    int Lx_n, int Ly_n,
    const std::function<std::size_t(int,int)> &idx)
{
    std::vector<double> z_copy = z_grid;
    std::vector<int>    n_copy = z_counts;

    for (int i = 0; i < Lx_n; ++i) {
        for (int j = 0; j < Ly_n; ++j) {
            std::size_t flat = idx(i, j);
            if (n_copy[flat] != 0) continue;

            int iu = wrap_index(i - 1, Lx_n);
            int id = wrap_index(i + 1, Lx_n);
            int jl = wrap_index(j - 1, Ly_n);
            int jr = wrap_index(j + 1, Ly_n);

            std::size_t up    = idx(iu, j);
            std::size_t down  = idx(id, j);
            std::size_t left  = idx(i, jl);
            std::size_t right = idx(i, jr);

            double z_sum = 0.0;
            int    w_sum = 0;

            for (std::size_t nb : {up, down, left, right}) {
                int w = n_copy[nb];
                if (w > 0) {
                    z_sum += z_copy[nb] * static_cast<double>(w);
                    w_sum += w;
                }
            }
            if (w_sum > 0) {
                z_grid[flat]   = z_sum / static_cast<double>(w_sum);
                z_counts[flat] = w_sum;
            }
        }
    }
}

// ------------------- FFT helpers -------------------
void fft2d_real_to_complex(
    const std::vector<double> &input,
    std::vector<std::complex<double>> &output,
    int Nx, int Ny)
{
    // FFTW r2c output size: Nx * (Ny/2 + 1)
    double *in = fftw_alloc_real(static_cast<int64_t>(Nx) * Ny);
    fftw_complex *out = fftw_alloc_complex(static_cast<int64_t>(Nx) * (Ny/2 + 1));

    std::copy(input.begin(), input.end(), in);
    fftw_plan plan = fftw_plan_dft_r2c_2d(Nx, Ny, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);

    output.resize(static_cast<size_t>(Nx) * (Ny/2 + 1));
    for (size_t i = 0; i < output.size(); ++i) {
        output[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

// Normalize to physical convention and average leaflets (v2 style)
void normalize_and_average_fft(
    const std::vector<std::complex<double>> &h_top,
    const std::vector<std::complex<double>> &h_bottom,
    std::vector<std::complex<double>> &hq,
    double Lx, double Ly,
    int Nx, int Ny)
{
    const int Ny_half = Ny / 2 + 1;
    const double norm = (Lx * Ly) / static_cast<double>(Nx * Ny); // dx*dy
    hq.resize(static_cast<size_t>(Nx) * Ny_half);
    for (int kx = 0; kx < Nx; ++kx) {
        for (int ky = 0; ky < Ny_half; ++ky) {
            size_t i = static_cast<size_t>(kx) * Ny_half + ky;
            hq[i] = 0.5 * (h_top[i] + h_bottom[i]) * norm;
        }
    }
}

// Phase shift for grid/cell-center alignment (v2 index-fraction form)
void apply_fft_phase_shift(
    std::vector<std::complex<double>> &hq,
    int Nx, int Ny)
{
    const int Ny_half = Ny / 2 + 1;
    for (int kx = 0; kx < Nx; ++kx) {
        int kx_eff = (kx <= Nx/2) ? kx : kx - Nx; // allow negative kx indices
        double phx = -PI * static_cast<double>(kx_eff) / static_cast<double>(Nx);
        for (int ky = 0; ky < Ny_half; ++ky) {
            double phy = -PI * static_cast<double>(ky) / static_cast<double>(Ny);
            double phase = phx + phy;
            hq[static_cast<size_t>(kx) * Ny_half + ky] *= std::polar(1.0, phase);
        }
    }
}

// ------------------- I/O for per-frame Fourier data -------------------


// ------------------- Spectrum & q-grid (per element) -------------------
void compute_power_spectrum(
    const std::vector<std::complex<double>> &hq,
    std::vector<double> &power,
    std::vector<double> &q_magnitudes,
    int Lx_n, int Ly_n,
    double Lx, double Ly)
{
    int Ny_half = Ly_n / 2 + 1;
    power.resize(static_cast<size_t>(Lx_n) * Ny_half);
    q_magnitudes.resize(static_cast<size_t>(Lx_n) * Ny_half);

    for (int kx = 0; kx < Lx_n; ++kx) {
        int kx_eff = (kx <= Lx_n / 2) ? kx : kx - Lx_n;
        double qx = (2.0 * PI / Lx) * static_cast<double>(kx_eff);
        for (int ky = 0; ky < Ny_half; ++ky) {
            double qy = (2.0 * PI / Ly) * static_cast<double>(ky);
            size_t idx = static_cast<size_t>(kx) * Ny_half + ky;
            const auto &h = hq[idx];
            power[idx] = std::norm(h);
            q_magnitudes[idx] = std::sqrt(qx*qx + qy*qy);
        }
    }
}

// ------------------- Leaflet assignment (from v1) -------------------
// Uses precomputed leaflet_ids: +1 = top, -1 = bottom
void assign_leaflet_heights_from_index(
    const chemfiles::Frame& frame,
    const std::vector<int>& leaflet_ids,  // +1 = top, -1 = bottom
    std::vector<double>& z_top,
    std::vector<int>& z_top_num,
    std::vector<double>& z_bottom,
    std::vector<int>& z_bottom_num,
    const std::function<std::size_t(int, int)>& idx,
    int Lx_n,
    int Ly_n,
    bool centered_coords
) {
    auto positions = frame.positions();
    auto cell = frame.cell();
    auto lengths = cell.lengths();

    for (size_t i = 0; i < positions.size(); ++i) {
        const auto& pos = positions[i];
        double x = centered_coords ? pos[0] + lengths[0] / 2.0 : pos[0];
        double y = centered_coords ? pos[1] + lengths[1] / 2.0 : pos[1];
        double z = centered_coords ? pos[2] + lengths[2] / 2.0 : pos[2];

        std::array<double, 3> wrapped = {x, y, z};
        check_boundaries(wrapped, lengths);

        int x_idx = std::min(static_cast<int>(std::floor(wrapped[0] * Lx_n / lengths[0])), Lx_n - 1);
        int y_idx = std::min(static_cast<int>(std::floor(wrapped[1] * Ly_n / lengths[1])), Ly_n - 1);
        std::size_t flat_idx = idx(x_idx, y_idx);

        if (leaflet_ids[i] == 1) {
            z_top[flat_idx] += z;
            z_top_num[flat_idx] += 1;
        } else if (leaflet_ids[i] == -1) {
            z_bottom[flat_idx] += z;
            z_bottom_num[flat_idx] += 1;
        }
    }

    for (size_t i = 0; i < z_top.size(); ++i) {
        if (z_top_num[i] > 0) z_top[i] /= z_top_num[i];
        if (z_bottom_num[i] > 0) z_bottom[i] /= z_bottom_num[i];
    }
}

// ------------------- FP error utilities -------------------
static std::string sigfig_fixed_label(double x, int sig = 6) {
    if (!(x > 0.0)) return std::string("0");
    int exp10 = static_cast<int>(std::floor(std::log10(x)));
    int decimals = std::max(0, sig - 1 - exp10);
    double scale = std::pow(10.0, sig - 1 - exp10);
    double rounded = std::round(x * scale) / scale;

    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(decimals) << rounded;

    std::string s = oss.str();
    if (s.find('.') != std::string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }
    if (s.empty()) s = "0";
    return s;
}

// Flyvbjerg–Petersen with guardrails + conservative fallback, with trace
std::vector<double> flyvbjerg_petersen_errors_with_trace(
    const std::vector<std::vector<double>> &data_series, // [frame][bin_index]
    const std::vector<double> &q_bins,                   // representative |q| per bin
    const std::string &trace_prefix = "blocking_trace_q")
{
    const size_t T = data_series.size();
    if (T == 0) return {};
    const size_t Q = data_series[0].size();
    if (q_bins.size() != Q) {
        throw std::runtime_error("flyvbjerg_petersen_errors_with_trace: q vector size mismatch");
    }

    std::vector<double> errs(Q, 0.0);
    createDirectory("./trace_error_fp_method");

    for (size_t q = 0; q < Q; ++q) {
        // Collect time series for this bin
        std::vector<double> x(T);
        for (size_t t = 0; t < T; ++t) x[t] = data_series[t][q];

        // Blocking ladder
        std::vector<double> blocked = x;
        std::vector<double> varmean;  // variance-of-the-mean at level k
        std::vector<size_t> Mks;      // number of blocks at level k

        while (blocked.size() > 1) {
            const size_t M = blocked.size();
            const double mean = std::accumulate(blocked.begin(), blocked.end(), 0.0) / static_cast<double>(M);

            double var = 0.0;
            for (double v : blocked) {
                const double d = v - mean;
                var += d * d;
            }
            var /= (M > 1 ? (M - 1) : 1);               // sample variance
            varmean.push_back(var / static_cast<double>(M)); // variance of overall mean
            Mks.push_back(M);

            // next level (pairwise averaging)
            std::vector<double> next;
            next.reserve(M / 2);
            for (size_t i = 0; i + 1 < M; i += 2) {
                next.push_back(0.5 * (blocked[i] + blocked[i+1]));
            }
            blocked.swap(next);
        }

        // Plateau detection with guardrails
        bool detected = false;
        double detected_plateau = 0.0;
        size_t detected_k = 0;

        for (size_t k = 1; k + FP_NEED_CONSEC < varmean.size(); ++k) {
            bool ok = true;
            for (size_t r = 0; r < FP_NEED_CONSEC; ++r) {
                const double v1 = varmean[k - 1 + r];
                const double v2 = varmean[k     + r];
                const double rel = std::abs(v2 - v1) / std::max({v1, v2, 1e-12});
                if (rel > FP_REL_EPS) { ok = false; break; }
            }
            if (ok && Mks[k] >= FP_MIN_BLOCKS) {
                double acc = 0.0;
                for (size_t r = 0; r <= FP_NEED_CONSEC; ++r) acc += varmean[k - 1 + r];
                detected_plateau = acc / static_cast<double>(FP_NEED_CONSEC + 1);
                detected_k = k;
                detected = true;
                break;
            }
        }

        // Reliable range (M_k >= FP_MIN_BLOCKS) and its maximum
        double reliable_max = 0.0;
        size_t reliable_k_of_max = 0;
        bool have_reliable = false;
        for (size_t k = 0; k < varmean.size(); ++k) {
            if (Mks[k] >= FP_MIN_BLOCKS) {
                if (!have_reliable || varmean[k] > reliable_max) {
                    reliable_max = varmean[k];
                    reliable_k_of_max = k;
                }
                have_reliable = true;
            }
        }

        // Choose final plateau (conservative)
        double chosen_plateau = 0.0;
        std::string status;

        if (detected) {
            if (have_reliable && detected_plateau < reliable_max) {
                chosen_plateau = reliable_max;
                status = "FOUND_CLIPPED_TO_RELIABLE_MAX";
            } else {
                chosen_plateau = detected_plateau;
                status = "FOUND_OK";
            }
        } else {
            if (have_reliable) {
                chosen_plateau = reliable_max;
                status = "NOT_FOUND_USED_RELIABLE_MAX";
            } else {
                chosen_plateau = (varmean.empty() ? 0.0 : varmean.back());
                status = "NOT_FOUND_NO_RELIABLE_USED_LAST";
            }
        }

        // Write trace file
        const double qabs = std::abs(q_bins[q]);
        const std::string qlabel6 = sigfig_fixed_label(qabs, 6);
        const std::string filepath =
            "./trace_error_fp_method/" + trace_prefix + std::to_string(q) + "_q" + qlabel6 + ".dat";

        std::ofstream out(filepath);
        out << std::setprecision(10);
        out << "# bin_index " << q
            << "  q_magnitude " << q_bins[q]
            << "  q_label_6sig " << qlabel6 << "\n";
        out << "# STATUS " << status
            << "  detected_plateau " << (detected ? detected_plateau : 0.0)
            << "  reliable_max " << (have_reliable ? reliable_max : 0.0)
            << "  reliable_k_of_max " << (have_reliable ? static_cast<long long>(reliable_k_of_max) : -1)
            << "  detected_k " << (detected ? static_cast<long long>(detected_k) : -1) << "\n";
        out << "# columns: level_k  var_of_mean  chosen_plateau  M_k\n";
        for (size_t k = 0; k < varmean.size(); ++k) {
            out << k << " " << varmean[k] << " " << chosen_plateau << " " << Mks[k] << "\n";
        }

        errs[q] = std::sqrt(chosen_plateau);
    }
    return errs;
}

// Fixed-blocking with exactly M blocks (even split) on binned series
std::vector<double> fixed_blocking_errors_M(
    const std::vector<std::vector<double>> &data_series, // [frame][bin_index]
    size_t desired_M,                                    // e.g., 10 blocks
    size_t &effective_M_out,                             // actual M used
    std::vector<size_t> &block_sizes_out)                // sizes of each block
{
    const size_t T = data_series.size();
    if (T == 0) { effective_M_out = 0; block_sizes_out.clear(); return {}; }
    const size_t Q = data_series[0].size();

    size_t M = std::min(desired_M, T);
    if (M < 2) {
        // Fallback: naive SEM per bin (unblocked)
        effective_M_out = 1;
        block_sizes_out.assign(1, T);
        std::vector<double> errs(Q, 0.0);
        for (size_t q = 0; q < Q; ++q) {
            double mean = 0.0;
            for (size_t t = 0; t < T; ++t) mean += data_series[t][q];
            mean /= static_cast<double>(T);
            double var = 0.0;
            for (size_t t = 0; t < T; ++t) {
                const double d = data_series[t][q] - mean;
                var += d*d;
            }
            var /= (T > 1 ? (T - 1) : 1);
            errs[q] = std::sqrt(var / static_cast<double>(T));
        }
        return errs;
    }

    // Distribute T samples into M blocks as evenly as possible
    const size_t base = T / M;
    const size_t rem  = T % M;
    std::vector<size_t> block_sizes(M, base);
    for (size_t i = 0; i < rem; ++i) block_sizes[i] += 1;

    while (!block_sizes.empty() && block_sizes.back() == 0) block_sizes.pop_back();
    M = block_sizes.size();
    if (M < 2) {
        effective_M_out = 1;
        block_sizes_out.assign(1, T);
        std::vector<double> errs(Q, 0.0);
        for (size_t q = 0; q < Q; ++q) {
            double mean = 0.0;
            for (size_t t = 0; t < T; ++t) mean += data_series[t][q];
            mean /= static_cast<double>(T);
            double var = 0.0;
            for (size_t t = 0; t < T; ++t) {
                const double d = data_series[t][q] - mean;
                var += d*d;
            }
            var /= (T > 1 ? (T - 1) : 1);
            errs[q] = std::sqrt(var / static_cast<double>(T));
        }
        return errs;
    }

    effective_M_out = M;
    block_sizes_out = block_sizes;

    std::vector<double> errs(Q, 0.0);

    // Compute block means per bin
    size_t offset = 0;
    std::vector<std::vector<double>> block_means(M, std::vector<double>(Q, 0.0));
    for (size_t m = 0; m < M; ++m) {
        const size_t B = block_sizes[m];
        for (size_t q = 0; q < Q; ++q) {
            double acc = 0.0;
            for (size_t t = 0; t < B; ++t) {
                acc += data_series[offset + t][q];
            }
            block_means[m][q] = acc / static_cast<double>(B);
        }
        offset += B;
    }

    // Sample variance of block means for each bin, then var_of_mean = var / M
    for (size_t q = 0; q < Q; ++q) {
        double mean = 0.0;
        for (size_t m = 0; m < M; ++m) mean += block_means[m][q];
        mean /= static_cast<double>(M);

        double var = 0.0;
        for (size_t m = 0; m < M; ++m) {
            const double d = block_means[m][q] - mean;
            var += d*d;
        }
        var /= (M > 1 ? (M - 1) : 1);
        errs[q] = std::sqrt(var / static_cast<double>(M));
    }

    return errs;
}

// ------------------- Binning helpers (rings in r2 = kx_eff^2+ky^2) -------------------
// Build mapping from (kx,ky) flat index to ring bin, using r2 integer.
// kx_eff = min(kx, Nx - kx) handles periodic wrap in x.
struct RingBins {
    // For each flat (kx,ky) index -> bin id [0..B-1]
    std::vector<int> flat_to_bin;
    // Unique r2 values per bin (monotonic increasing)
    std::vector<int> r2_vals;
    // Number of (kx,ky) members per bin
    std::vector<int> bin_degeneracy;
};

RingBins build_ring_bins(int Nx, int Ny) {
    const int Ny_half = Ny/2 + 1;
    const size_t total = static_cast<size_t>(Nx) * Ny_half;

    std::unordered_map<int, int> r2_to_bin; // r2 -> bin index
    r2_to_bin.reserve(total);

    // First collect all r2 present
    std::vector<int> r2_list;
    r2_list.reserve(total);

    for (int kx = 0; kx < Nx; ++kx) {
        int kx_eff = std::min(kx, Nx - kx); // non-negative
        for (int ky = 0; ky < Ny_half; ++ky) {
            int r2 = kx_eff*kx_eff + ky*ky;
            r2_list.push_back(r2);
        }
    }
    // Unique & sort r2
    std::sort(r2_list.begin(), r2_list.end());
    r2_list.erase(std::unique(r2_list.begin(), r2_list.end()), r2_list.end());

    // Map r2 -> bin id
    for (size_t i = 0; i < r2_list.size(); ++i) {
        r2_to_bin[r2_list[i]] = static_cast<int>(i);
    }

    // Fill flat_to_bin and count degeneracy
    std::vector<int> flat_to_bin(total, -1);
    std::vector<int> degeneracy(r2_list.size(), 0);

    for (int kx = 0; kx < Nx; ++kx) {
        int kx_eff = std::min(kx, Nx - kx);
        for (int ky = 0; ky < Ny_half; ++ky) {
            int r2 = kx_eff*kx_eff + ky*ky;
            int bin = r2_to_bin[r2];
            size_t flat = static_cast<size_t>(kx) * Ny_half + ky;
            flat_to_bin[flat] = bin;
            degeneracy[static_cast<size_t>(bin)]++;
        }
    }

    RingBins rb;
    rb.flat_to_bin = std::move(flat_to_bin);
    rb.r2_vals     = std::move(r2_list);
    rb.bin_degeneracy = std::move(degeneracy);
    return rb;
}

// --------------------------------------------------------------------------------------------
int main() {
    int num_frames = 0, num_atoms = 0;
    double Lx_mean = 0.0, Ly_mean = 0.0, Lz_mean = 0.0;


    try {
        // First pass for trajectory info
        chemfiles::Trajectory traj_info(trajectory_file);
        trajectory_info(traj_info, num_frames, num_atoms, Lx_mean, Ly_mean, Lz_mean);
        std::cout << "Grid spacing (length units): " << L_grid << "\n";
        std::cout << "Detrending (REMOVE_MEAN): " << (REMOVE_MEAN ? "ON" : "OFF") << "\n";

        // (Topology not required for index-based IDs; keep file open to verify availability)
        chemfiles::Trajectory pdb_traj(pdb_file);
        chemfiles::Frame pdb_frame = pdb_traj.read();
        (void)pdb_frame;

        // Grid sizes from mean box
        int Lx_n = static_cast<int>(std::ceil(Lx_mean / L_grid));
        int Ly_n = static_cast<int>(std::ceil(Ly_mean / L_grid));
        if (Lx_n <= 0 || Ly_n <= 0) {
            std::cerr << "Invalid grid sizes: Lx_n=" << Lx_n << ", Ly_n=" << Ly_n << "\n";
            return 1;
        }
        std::cout << "Lx_n = " << Lx_n << ", Ly_n = " << Ly_n << "\n";
        auto idx2d = [Ly_n](int i, int j) { return static_cast<std::size_t>(i) * Ly_n + j; };
        const int Ny_half = Ly_n/2 + 1;

        // Build ring bins on index space (r2 = kx_eff^2 + ky^2)
        RingBins bins = build_ring_bins(Lx_n, Ly_n);
        const size_t B = bins.r2_vals.size();
        std::cout << "Ring bins built: B = " << B << " unique r2 values\n";

        std::cout << "-------------------------------------------------\n";
        std::cout << "-------- Calculating height spectrum  -----------\n";
        std::cout << "-------------------------------------------------\n";

        // Second pass: extract coordinates
        chemfiles::Trajectory trajectory(trajectory_file);

        // Buffers: real-space grids (z) and r2c Fourier arrays
        std::vector<double> z_top        (static_cast<size_t>(Lx_n) * Ly_n, 0.0);
        std::vector<int>    z_top_num    (static_cast<size_t>(Lx_n) * Ly_n, 0);
        std::vector<double> z_bottom     (static_cast<size_t>(Lx_n) * Ly_n, 0.0);
        std::vector<int>    z_bottom_num (static_cast<size_t>(Lx_n) * Ly_n, 0);

        std::vector<std::complex<double>> h_top   (static_cast<size_t>(Lx_n) * Ny_half, {0.0,0.0});
        std::vector<std::complex<double>> h_bottom(static_cast<size_t>(Lx_n) * Ny_half, {0.0,0.0});
        std::vector<std::complex<double>> hq      (static_cast<size_t>(Lx_n) * Ny_half, {0.0,0.0});

        // Leaflet IDs (index-based split): first half top (+1), second half bottom (-1)
        // This matches the typical PN-only trajectories where atoms are ordered by leaflet.
        std::vector<int> leaflet_ids(static_cast<size_t>(num_atoms), 0);
        for (int i = 0; i < num_atoms; ++i) {
            leaflet_ids[static_cast<size_t>(i)] = (i < num_atoms/2) ? 1 : -1;
        }


        // For continuity: also compute an unbinned average spectrum (no errors here)
        std::vector<double> power_spectrum_sum(static_cast<size_t>(Lx_n) * Ny_half, 0.0);
        std::vector<double> first_frame_q; // q-values from the first usable frame

        // Binned time series (one vector per frame, length = number of rings B)
        std::vector<std::vector<double>> binned_power_frames;
        binned_power_frames.reserve(static_cast<size_t>(num_frames));

        // For each ring/bin we keep a representative q as the average of the per-frame mean q in that ring
        std::vector<double> q_bin_sum(B, 0.0);

        // Temporary per-frame arrays
        std::vector<double> power_tmp, q_tmp;
        size_t frame_index = 0;        // counts all frames read (including discarded)
        size_t used_frames = 0;        // counts frames actually used in averages
        size_t discarded_frames = 0;   // bad cell dimensions, etc.

        while (true) {
            try {
                auto frame = trajectory.read();
                frame_index++;

                const auto cell = frame.cell();
                const auto Lf   = cell.lengths();

                // Skip frames with invalid box (can happen in some trajectories)
                if (Lf[0] == 0.0 || Lf[1] == 0.0 || Lf[2] == 0.0) {
                    discarded_frames++;
                    continue;
                }

                // Clear per-frame accumulators
                std::fill(z_top.begin(),        z_top.end(),        0.0);
                std::fill(z_top_num.begin(),    z_top_num.end(),    0);
                std::fill(z_bottom.begin(),     z_bottom.end(),     0.0);
                std::fill(z_bottom_num.begin(), z_bottom_num.end(), 0);

                // 1) Map atoms to top/bottom leaflets (index-based split) and grid z(x,y)
                assign_leaflet_heights_from_index(
                    frame, leaflet_ids,
                    z_top, z_top_num,
                    z_bottom, z_bottom_num,
                    idx2d, Lx_n, Ly_n, centered
                );

                // 2) Fill empty grid cells by weighted neighbor averaging
                fill_empty_cells(z_top,    z_top_num,    Lx_n, Ly_n, idx2d);
                fill_empty_cells(z_bottom, z_bottom_num, Lx_n, Ly_n, idx2d);

                // 3) Optional detrending: remove spatial mean per leaflet (helps suppress q=0 drift)
                if (REMOVE_MEAN) {
                    auto demean = [](std::vector<double> &z) {
                        const double m = std::accumulate(z.begin(), z.end(), 0.0)
                                       / static_cast<double>(z.size());
                        for (double &v : z) v -= m;
                    };
                    demean(z_top);
                    demean(z_bottom);
                }

                // 4) FFT each leaflet (real -> complex, r2c layout)
                fft2d_real_to_complex(z_top,    h_top,    Lx_n, Ly_n);
                fft2d_real_to_complex(z_bottom, h_bottom, Lx_n, Ly_n);

                // 5) Midplane Fourier mode: h_q = (h_top + h_bottom)/2, normalized by box area
                normalize_and_average_fft(h_top, h_bottom, hq, Lf[0], Lf[1], Lx_n, Ly_n);

                // 6) Apply a phase factor so the real-space origin is centered consistently
                apply_fft_phase_shift(hq, Lx_n, Ly_n);

                // 7) Power spectrum for this frame: |h_q|^2 and corresponding q magnitudes
                compute_power_spectrum(hq, power_tmp, q_tmp, Lx_n, Ly_n, Lf[0], Lf[1]);

                // Unbinned running sum (for a legacy/continuity output)
                for (size_t i = 0; i < power_tmp.size(); ++i) {
                    power_spectrum_sum[i] += power_tmp[i];
                }
                if (used_frames == 0) {
                    first_frame_q = q_tmp; // use q-grid from first usable frame
                }

                // Ring-bin this frame (mean within each r2 ring)
                std::vector<double> binned_frame(B, 0.0);
                std::vector<double> qsum_frame(B, 0.0);
                std::vector<int>    cnt_frame(B, 0);
                for (size_t flat = 0; flat < power_tmp.size(); ++flat) {
                    int bin = bins.flat_to_bin[flat];
                    if (bin < 0) continue;
                    binned_frame[static_cast<size_t>(bin)] += power_tmp[flat];
                    qsum_frame[static_cast<size_t>(bin)]   += q_tmp[flat];
                    cnt_frame[static_cast<size_t>(bin)]    += 1;
                }
                for (size_t b = 0; b < B; ++b) {
                    if (cnt_frame[b] > 0) {
                        binned_frame[b] /= static_cast<double>(cnt_frame[b]);
                        q_bin_sum[b]    += qsum_frame[b] / static_cast<double>(cnt_frame[b]);
                    }
                }
                binned_power_frames.push_back(std::move(binned_frame));

                used_frames++;

            } catch (const chemfiles::Error &) {
                break; // EOF
            }
        }

                if (used_frames == 0) {
            std::cerr << "No usable frames after trajectory read.\n";
            return 1;
        }

        // Average unbinned spectrum (continuity output)
        {
            std::ofstream fout("height_spectrum_avg.dat");
            if (!fout) { std::cerr << "Failed to write height_spectrum_avg.dat\n"; }
            else {
                fout << std::fixed << std::setprecision(10);
                for (size_t i = 0; i < power_spectrum_sum.size(); ++i) {
                    double avg_power = power_spectrum_sum[i] / static_cast<double>(used_frames);
                    fout << (first_frame_q.empty()? 0.0 : first_frame_q[i]) << " " << avg_power << "\n";
                }
                std::cout << "Averaged (unbinned) height spectrum written to height_spectrum_avg.dat\n";
            }
        }

        // Final per-bin representative q: average across frames of per-frame mean q in each bin
        std::vector<double> q_bins(B, 0.0);
        for (size_t b = 0; b < B; ++b) q_bins[b] = q_bin_sum[b] / static_cast<double>(used_frames);

        // Average binned spectrum
        std::vector<double> avg_binned(B, 0.0);
        for (size_t f = 0; f < used_frames; ++f) {
            const auto &vec = binned_power_frames[f];
            for (size_t b = 0; b < B; ++b) avg_binned[b] += vec[b];
        }
        for (size_t b = 0; b < B; ++b) avg_binned[b] /= static_cast<double>(used_frames);

        // Write averaged binned spectrum (includes q=0 bin; keep for completeness)
        {
            std::ofstream fout("height_spectrum_avg_binned.dat");
            if (!fout) { std::cerr << "Failed to write height_spectrum_avg_binned.dat\n"; }
            else {
                fout << std::fixed << std::setprecision(10);
                for (size_t b = 0; b < B; ++b) {
                    fout << q_bins[b] << " " << avg_binned[b] << "\n";
                }
                std::cout << "Averaged (binned) height spectrum written to height_spectrum_avg_binned.dat\n";
            }
        }

        // FP errors on binned series (skip q=0 at write time)
        auto fp_errors_binned = flyvbjerg_petersen_errors_with_trace(binned_power_frames, q_bins);

        {
            std::ofstream fout_err("height_spectrum_avg_error_binned.dat");
            if (!fout_err) { std::cerr << "Failed to write height_spectrum_avg_error_binned.dat\n"; }
            else {
                fout_err << std::fixed << std::setprecision(10);
                for (size_t b = 0; b < B; ++b) {
                    if (q_bins[b] == 0.0) continue; // skip zero-mode ring
                    fout_err << q_bins[b] << " " << avg_binned[b] << " " << fp_errors_binned[b] << "\n";
                }
                std::cout << "Averaged (binned) spectrum with FP error written to height_spectrum_avg_error_binned.dat\n";
            }
        }

        // Fixed-blocking errors with M=10 blocks on binned series
        {
            size_t eff_M = 0;
            std::vector<size_t> block_sizes;
            auto errors_fixed = fixed_blocking_errors_M(binned_power_frames, FIXED_BLOCKS_M, eff_M, block_sizes);

            std::cout << "Fixed-blocking (binned): desired M=" << FIXED_BLOCKS_M
                      << ", effective M=" << eff_M << " (block sizes:";
            for (size_t bsz : block_sizes) std::cout << " " << bsz;
            std::cout << "), frames used=" << used_frames << "\n";

            std::ofstream fout_fix("height_spectrum_avg_error_fix_blocking_binned.dat");
            if (!fout_fix) {
                std::cerr << "Failed to write height_spectrum_avg_error_fix_blocking_binned.dat\n";
            } else {
                fout_fix << std::fixed << std::setprecision(10);
                for (size_t b = 0; b < B; ++b) {
                    if (q_bins[b] == 0.0) continue; // skip zero-mode ring
                    fout_fix << q_bins[b] << " " << avg_binned[b] << " " << errors_fixed[b] << "\n";
                }
                std::cout << "Averaged (binned) spectrum with fixed-blocking error written to height_spectrum_avg_error_fix_blocking_binned.dat\n";
            }
        }

    } catch (const chemfiles::Error &e) {
        std::cerr << "Chemfiles error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}