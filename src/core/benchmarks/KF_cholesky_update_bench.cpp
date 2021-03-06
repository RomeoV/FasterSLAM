#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "KF_cholesky_update.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

//
// DataLoader for Benchmark 1
//
auto data_loader(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H) {
#ifdef KF_YGLEE
    Vector2d x_init = {3.2403905331533212, -25.689432087069857};
    Matrix2d P_init = {0.20063369668655512, 0.018909593226709744, 0.018909593226709744, 0.011875705723671498};
#else
    double x_init[2] = {3.227460886446243, -25.613382543676146};
    double P_init[4] = {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843};
#endif
    std::copy(x_init, x_init+2, x);
    std::copy(P_init, P_init+4, P);
}

//
// DataLoader for Benchmark 2
//
void data_loader_unrolled4_avx(__m256d *x0x2,
                               __m256d *x1x3,
                               __m256d *P0,
                               __m256d *P1,
                               __m256d *P2,
                               __m256d *P3,
                               __m256d const v0v2,
                               __m256d const v1v3,
                               __m256d const R,
                               __m256d const H0,
                               __m256d const H1,
                               __m256d const H2,
                               __m256d const H3) 
{
    double x[4] __attribute__((aligned(32))) = {3.227460886446243, -25.613382543676146, 3.227460886446243, -25.613382543676146};
    double P[4] __attribute__((aligned(32))) = {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843};
    *x0x2 = _mm256_load_pd( x );
    *x1x3 = _mm256_load_pd( x );
    *P0   = _mm256_load_pd( P );
    *P1   = _mm256_load_pd( P );
    *P2   = _mm256_load_pd( P );
    *P3   = _mm256_load_pd( P );
}

template<typename func>
std::function<void (__m256d*,
                    __m256d*,
                    __m256d*,
                    __m256d*,
                    __m256d*,
                    __m256d*,
                    __m256d const,
                    __m256d const,
                    __m256d const,
                    __m256d const,
                    __m256d const,
                    __m256d const,
                    __m256d const)>
KF_cholesky_update_x4(func KF_cholesky_update_x1) {

    return [=] (__m256d *x0x2,
                __m256d *x1x3,
                __m256d *P0,
                __m256d *P1,
                __m256d *P2,
                __m256d *P3,
                __m256d const v0v2,
                __m256d const v1v3,
                __m256d const R,
                __m256d const H0,
                __m256d const H1,
                __m256d const H2,
                __m256d const H3) -> void {
        double x0x2_[4] __attribute__((aligned(32)));
        double x1x3_[4] __attribute__((aligned(32)));
        double v0v2_[4] __attribute__((aligned(32)));
        double v1v3_[4] __attribute__((aligned(32)));
        double   P0_[4] __attribute__((aligned(32)));
        double   P1_[4] __attribute__((aligned(32)));
        double   P2_[4] __attribute__((aligned(32)));
        double   P3_[4] __attribute__((aligned(32)));
        double    R_[4] __attribute__((aligned(32)));
        double   H0_[4] __attribute__((aligned(32)));
        double   H1_[4] __attribute__((aligned(32)));
        double   H2_[4] __attribute__((aligned(32)));
        double   H3_[4] __attribute__((aligned(32)));

        _mm256_store_pd(x0x2_, *x0x2);
        _mm256_store_pd(x1x3_, *x1x3);
        _mm256_store_pd(P0_, *P0);
        _mm256_store_pd(P1_, *P1);
        _mm256_store_pd(P2_, *P2);
        _mm256_store_pd(P3_, *P3);
        _mm256_store_pd(v0v2_, v0v2);
        _mm256_store_pd(v1v3_, v1v3);
        _mm256_store_pd(R_, R); 
        _mm256_store_pd(H0_, H0);
        _mm256_store_pd(H1_, H1);
        _mm256_store_pd(H2_, H2);
        _mm256_store_pd(H3_, H3);

        KF_cholesky_update_x1(x0x2_+0, P0_, v0v2_+0, R_, H0_);
        KF_cholesky_update_x1(x1x3_+0, P1_, v1v3_+0, R_, H1_);
        KF_cholesky_update_x1(x0x2_+2, P2_, v0v2_+2, R_, H2_);
        KF_cholesky_update_x1(x1x3_+2, P3_, v1v3_+2, R_, H3_);

        *x0x2 = _mm256_load_pd(x0x2_);
        *x1x3 = _mm256_load_pd(x1x3_);
        *P0   = _mm256_load_pd( P0_ );
        *P1   = _mm256_load_pd( P1_ );
        *P2   = _mm256_load_pd( P2_ );
        *P3   = _mm256_load_pd( P3_ );
    };
}

int main() {

    // ----------------------------------------- //
    // ----------------- TEST 1 ---------------- //
    // ----------------------------------------- //
#ifdef KF_YGLEE
    cVector2d v = {-0.017001037783700212, -0.010645013219889199};
    cMatrix2d R = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    cMatrix2d H = {0.073819203427568675, -0.99727164063023443, 0.038893076096335785, 0.0028789105989867826};    
#else
    double v[2] = {0.128762949830296, 0.019814250533567};
    double R[4] = {0.010000000000000, 0.0, 0.0, 0.000304617419787};
    double H[4] = {0.075431770172036, -0.997150965525639, 0.038902059641499, 0.002942835461779};
#endif
    Vector2d exact_x, x, x_v10, x_v11, x_rf, x_rf2;
    Matrix2d exact_P, P, P_v10, P_v11, P_rf, P_rf2;

    data_loader(exact_x, exact_P, v, R, H);
    KF_cholesky_update_base(exact_x, exact_P, v, R, H);
    
    data_loader(x_v10, P_v10, v, R, H);
    KF_cholesky_update_fused_ops(x_v10, P_v10, v, R, H);
    
    data_loader(x_v11, P_v11, v, R, H);
    KF_cholesky_update_fused_ops_avx(x_v11, P_v11, v, R, H);

#ifndef KF_YGLEE    
    data_loader(x_rf, P_rf, v, R, H);
    KF_cholesky_update_reduced_flops(x_rf, P_rf, v, R, H);
    
    data_loader(x_rf2, P_rf2, v, R, H);
    KF_cholesky_update_reduced_flops_avx(x_rf2, P_rf2, v, R, H);
#endif
     
    data_loader(x, P, v, R, H);
    KF_cholesky_update(x, P, v, R, H);

    // Check x is the same as in base function
    double error = 0.0;
    for (int i = 0; i < 2; i++) {
        error = fabs(        x[i] - exact_x[i] ); expect(that % error < 1e-10) << i;
        error = fabs(     x_v10[i] - exact_x[i] ); expect(that % error < 1e-10) << i + 10;
        error = fabs(     x_v11[i] - exact_x[i] ); expect(that % error < 1e-10) << i + 11;
    #ifndef KF_YGLEE
        error = fabs(     x_rf[i] - exact_x[i] ); expect(that % error < 1e-10) << i + 20;
        error = fabs(     x_rf2[i] - exact_x[i] ); expect(that % error < 1e-10) << i + 21;
    #endif
    }
    // Check P is the same as in base function
    for (int i = 0; i < 4; i++) {
        error = fabs(        P[i] - exact_P[i] ); expect(that % error < 1e-10) << i;
        error = fabs(     P_v10[i] - exact_P[i] ); expect(that % error < 1e-10) << i + 10;
        error = fabs(     P_v11[i] - exact_P[i] ); expect(that % error < 1e-10) << i + 11;
    #ifndef KF_YGLEE
        error = fabs(     P_rf[i] - exact_P[i] ); expect(that % error < 1e-10) << i + 20;
        error = fabs(     P_rf2[i] - exact_P[i] ); expect(that % error < 1e-10) << i + 21;
    #endif
    }
 
    // ----------------------------------------- //
    // -------------- BENCHMARK 1 -------------- //
    // ----------------------------------------- //

    Benchmark<decltype(&KF_cholesky_update)> bench("KF_cholesky_update benchmark (SINGLE CALL)"); 
    bench.data_loader = data_loader; 
    
#ifdef KF_YGLEE
    // Base
    bench.add_function(&KF_cholesky_update_base, "base", 0.0);
    bench.funcFlops[0] = KF_cholesky_update_base_flops(x, P, v, R, H);
    bench.funcBytes[0] = 8*KF_cholesky_update_base_memory(x, P, v, R, H);
    // Fused Ops
    bench.add_function(&KF_cholesky_update_fused_ops, "fused-ops-noAVX", 0.0);
    bench.funcFlops[1] = KF_cholesky_update_active_flops(x, P, v, R, H);
    bench.funcBytes[1] = 8*KF_cholesky_update_active_memory(x, P, v, R, H);
    // Fused Ops + AVX
    bench.add_function(&KF_cholesky_update_fused_ops_avx, "fused-ops-AVX", 0.0);
    bench.funcFlops[2] = KF_cholesky_update_active_flops(x, P, v, R, H);
    bench.funcBytes[2] = 8*KF_cholesky_update_active_memory(x, P, v, R, H);
    // Active = Fused Ops + AVX
    bench.add_function(&KF_cholesky_update, "active", 0.0);
    bench.funcFlops[3] = KF_cholesky_update_active_flops(x, P, v, R, H);
    bench.funcBytes[3] = 8*KF_cholesky_update_active_memory(x, P, v, R, H);
#else
    // Base
    bench.add_function(&KF_cholesky_update_base, "base", 0.0);
    bench.funcFlops[0] = KF_cholesky_update_base_flops(x, P, v, R, H);
    bench.funcBytes[0] = 8*KF_cholesky_update_base_memory(x, P, v, R, H);
    // Fused Ops
    bench.add_function(&KF_cholesky_update_fused_ops, "fused-ops-noAVX", 0.0);
    bench.funcFlops[1] = KF_cholesky_update_base_flops(x, P, v, R, H);
    bench.funcBytes[1] = 8*KF_cholesky_update_base_memory(x, P, v, R, H);
    // Fused Ops + AVX
    bench.add_function(&KF_cholesky_update_fused_ops_avx, "fused-ops-AVX", 0.0);
    bench.funcFlops[2] = KF_cholesky_update_base_flops(x, P, v, R, H);
    bench.funcBytes[2] = 8*KF_cholesky_update_base_memory(x, P, v, R, H);
    // Reduced Flops  
    bench.add_function(&KF_cholesky_update_reduced_flops, "reduced-flops-noAVX", 0.0);
    bench.funcFlops[3] = KF_cholesky_update_active_flops(x, P, v, R, H);
    bench.funcBytes[3] = 8*KF_cholesky_update_active_memory(x, P, v, R, H);
    // Reduced Flops + AVX
    bench.add_function(&KF_cholesky_update_reduced_flops_avx, "reduced-flops-AVX", 0.0);
    bench.funcFlops[4] = KF_cholesky_update_active_flops(x, P, v, R, H);
    bench.funcBytes[4] = 8*KF_cholesky_update_active_memory(x, P, v, R, H);
    // Active
    bench.add_function(&KF_cholesky_update, "active", 0.0);
    bench.funcFlops[5] = KF_cholesky_update_active_flops(x, P, v, R, H);
    bench.funcBytes[5] = 8*KF_cholesky_update_active_memory(x, P, v, R, H);
#endif

    bench.run_benchmark(x, P, v, R, H);
    
    //
    // END OF BENCHMARK 1
    //
#ifndef KF_YGLEE    
    // ----------------------------------------- //
    // ---------------- TEST 2 ----------------- //
    // ----------------------------------------- //

    double x_u4avx[4] __attribute__((aligned(32))) = {3.227460886446243, -25.613382543676146, 3.227460886446243, -25.613382543676146};
    double P_u4avx[4] __attribute__((aligned(32))) = {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843};
    double v_u4avx[4] __attribute__((aligned(32))) = {0.128762949830296, 0.019814250533567, 0.128762949830296, 0.019814250533567};
    double R_u4avx[4] __attribute__((aligned(32))) = {0.010000000000000, 0.0, 0.0, 0.000304617419787};
    double H_u4avx[4] __attribute__((aligned(32))) = {0.075431770172036, -0.997150965525639, 0.038902059641499, 0.002942835461779};
    
    __m256d x0x2 = _mm256_load_pd( x_u4avx );
    __m256d x1x3 = _mm256_load_pd( x_u4avx );
    __m256d P0   = _mm256_load_pd( P_u4avx );
    __m256d P1   = _mm256_load_pd( P_u4avx );
    __m256d P2   = _mm256_load_pd( P_u4avx );
    __m256d P3   = _mm256_load_pd( P_u4avx );
    __m256d v0v2 = _mm256_load_pd( v_u4avx );
    __m256d v1v3 = _mm256_load_pd( v_u4avx );
    __m256d RR   = _mm256_load_pd( R_u4avx );
    __m256d H0   = _mm256_load_pd( H_u4avx );
    __m256d H1   = _mm256_load_pd( H_u4avx );
    __m256d H2   = _mm256_load_pd( H_u4avx );
    __m256d H3   = _mm256_load_pd( H_u4avx );

    KF_cholesky_update_unrolled4_avx(&x0x2, &x1x3, &P0, &P1, &P2, &P3, v0v2, v1v3, RR, H0, H1, H2, H3);

    double x0x2_[4], x1x3_[4], P0_[4], P1_[4], P2_[4], P3_[4]; 
    _mm256_store_pd( x0x2_, x0x2 );
    _mm256_store_pd( x1x3_, x1x3 );
    _mm256_store_pd( P0_, P0 );
    _mm256_store_pd( P1_, P1 );
    _mm256_store_pd( P2_, P2 );
    _mm256_store_pd( P3_, P3 );

    double actual_x[4] = {3.470171202213126, -25.656742169761873, 3.470171202213126, -25.656742169761873};
    double actual_P[4] = {0.099441292170602, 0.008341518741166, 0.008341518741166, 0.005737523050281};
    for (int i = 0; i < 4; i++) {
        expect(fabs(x0x2_[i] - actual_x[i]) < 1e-10) << x0x2_[i] << " != " << actual_x[i];
        expect(fabs(x1x3_[i] - actual_x[i]) < 1e-10) << x1x3_[i] << " != " << actual_x[i];
    }
    for (int i = 0; i < 4; i++) {
        expect(fabs(P0[i] - actual_P[i]) < 1e-10) << P0[i] << " != " << actual_P[i];
        expect(fabs(P1[i] - actual_P[i]) < 1e-10) << P1[i] << " != " << actual_P[i];
        expect(fabs(P2[i] - actual_P[i]) < 1e-10) << P2[i] << " != " << actual_P[i];
        expect(fabs(P3[i] - actual_P[i]) < 1e-10) << P3[i] << " != " << actual_P[i];
    }
    
    // ----------------------------------------- //
    // -------------- BENCHMARK 2 -------------- //
    // ----------------------------------------- //
    auto KF_cholesky_update_base_x4              = KF_cholesky_update_x4(KF_cholesky_update_base);
    auto KF_cholesky_update_fused_ops_x4         = KF_cholesky_update_x4(KF_cholesky_update_fused_ops);
    auto KF_cholesky_update_fused_ops_avx_x4     = KF_cholesky_update_x4(KF_cholesky_update_fused_ops_avx);
    auto KF_cholesky_update_reduced_flops_x4     = KF_cholesky_update_x4(KF_cholesky_update_reduced_flops);
    auto KF_cholesky_update_reduced_flops_avx_x4 = KF_cholesky_update_x4(KF_cholesky_update_reduced_flops_avx);

    Benchmark<decltype(KF_cholesky_update_base_x4)> bench_u4avx("KF_cholesky_update_unrolled4_avx benchmark (UNROLLEDx4)");
    bench_u4avx.data_loader = data_loader_unrolled4_avx;

    bench_u4avx.add_function(KF_cholesky_update_base_x4, "base_x4", 0.0);
    bench_u4avx.funcFlops[0] = 4*KF_cholesky_update_base_flops(x, P, v, R, H); // From Bench-1
    bench_u4avx.funcBytes[0] = 4*8*KF_cholesky_update_base_memory(x, P, v, R, H);  // From Bench-1
    
    bench_u4avx.add_function(KF_cholesky_update_fused_ops_x4, "fused-ops-noAVX_x4", 0.0);
    bench_u4avx.funcFlops[1] = 4*KF_cholesky_update_base_flops(x, P, v, R, H); // From Bench-1
    bench_u4avx.funcBytes[1] = 4*8*KF_cholesky_update_base_memory(x, P, v, R, H);  // From Bench-1

    bench_u4avx.add_function(KF_cholesky_update_fused_ops_avx_x4, "fused-ops-AVX_x4", 0.0);
    bench_u4avx.funcFlops[2] = 4*KF_cholesky_update_base_flops(x, P, v, R, H); // From Bench-1
    bench_u4avx.funcBytes[2] = 4*8*KF_cholesky_update_base_memory(x, P, v, R, H);  // From Bench-1

    bench_u4avx.add_function(KF_cholesky_update_reduced_flops_x4, "reduced-flops-noAVX_x4", 0.0);
    bench_u4avx.funcFlops[3] = 4*KF_cholesky_update_active_flops(x, P, v, R, H); // From Bench-1
    bench_u4avx.funcBytes[3] = 4*8*KF_cholesky_update_active_memory(x, P, v, R, H);  // From Bench-1

    bench_u4avx.add_function(KF_cholesky_update_reduced_flops_avx_x4, "reduced-flops-AVX_x4", 0.0);
    bench_u4avx.funcFlops[4] = 4*KF_cholesky_update_active_flops(x, P, v, R, H); // From Bench-1
    bench_u4avx.funcBytes[4] = 4*8*KF_cholesky_update_active_memory(x, P, v, R, H);  // From Bench-1

    bench_u4avx.add_function(&KF_cholesky_update_unrolled4_avx, "unrolled4-AVX", 0.0);
    bench_u4avx.funcFlops[5] = 4*KF_cholesky_update_active_flops(x, P, v, R, H); // From Bench-1
    bench_u4avx.funcBytes[5] = 4*8*KF_cholesky_update_active_memory(x, P, v, R, H);  // From Bench-1

    bench_u4avx.run_benchmark(&x0x2, &x1x3, &P0, &P1, &P2, &P3, v0v2, v1v3, RR, H0, H1, H2, H3);
    
    //
    // END OF BENCHMARK 2
    //
#endif

    return 0;
}

