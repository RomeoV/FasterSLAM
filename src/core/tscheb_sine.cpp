#include "tscheb_sine.h"
#include "pi_to_pi.h"

// assumes angle is normalized!
void tscheb_dsines(double* alphas, size_t N, double* results) {

    double coeffs[] = {
        -0.10132118,          // x
         0.0066208798,        // x^3
        -0.00017350505,       // x^5
         0.0000025222919,     // x^7
        -0.000000023317787,   // x^9
         0.00000000013291342, // x^11

    };
    double pi_major = 3.1415927;
    double pi_minor = -0.00000008742278;

    for (size_t i = 0; i < N; i++) {
        double x2 = alphas[i]*alphas[i];
        double p11 = coeffs[5];
        double p9  = p11*x2 + coeffs[4];
        double p7  = p9*x2  + coeffs[3];
        double p5  = p7*x2  + coeffs[2];
        double p3  = p5*x2  + coeffs[1];
        double p1  = p3*x2  + coeffs[0];
        results[i] =  (alphas[i] - pi_major - pi_minor) *
            (alphas[i] + pi_major + pi_minor) * p1 * alphas[i];
    }
}

// assumes angle is normalized!
void tscheb_dsines_unrolled(double* alphas, size_t N, double* results) {

    double coeffs[] = {
        -0.10132118,          // x
         0.0066208798,        // x^3
        -0.00017350505,       // x^5
         0.0000025222919,     // x^7
        -0.000000023317787,   // x^9
         0.00000000013291342, // x^11

    };
    double pi_major = 3.1415927;
    double pi_minor = -0.00000008742278;

    size_t i;
    for (i = 0; i < N; i+=4) {
        double x0_2 = alphas[i]*alphas[i];
        double x1_2 = alphas[i]*alphas[i];
        double x2_2 = alphas[i]*alphas[i];
        double x3_2 = alphas[i]*alphas[i];

        double p0_11 = coeffs[5];
        double p1_11 = coeffs[5];
        double p2_11 = coeffs[5];
        double p3_11 = coeffs[5];

        double p0_9  = p0_11*x0_2 + coeffs[4];
        double p1_9  = p1_11*x1_2 + coeffs[4];
        double p2_9  = p2_11*x2_2 + coeffs[4];
        double p3_9  = p3_11*x3_2 + coeffs[4];

        double p0_7  = p0_9*x0_2  + coeffs[3];
        double p1_7  = p1_9*x1_2  + coeffs[3];
        double p2_7  = p2_9*x2_2  + coeffs[3];
        double p3_7  = p3_9*x3_2  + coeffs[3];

        double p0_5  = p0_7*x0_2  + coeffs[2];
        double p1_5  = p1_7*x1_2  + coeffs[2];
        double p2_5  = p2_7*x2_2  + coeffs[2];
        double p3_5  = p3_7*x3_2  + coeffs[2];

        double p0_3  = p0_5*x0_2  + coeffs[1];
        double p1_3  = p1_5*x1_2  + coeffs[1];
        double p2_3  = p2_5*x2_2  + coeffs[1];
        double p3_3  = p3_5*x3_2  + coeffs[1];

        double p0_1  = p0_3*x0_2  + coeffs[0];
        double p1_1  = p1_3*x1_2  + coeffs[0];
        double p2_1  = p2_3*x2_2  + coeffs[0];
        double p3_1  = p3_3*x3_2  + coeffs[0];

        results[i+0] =  (alphas[i+0] - pi_major - pi_minor) *
            (alphas[i+0] + pi_major + pi_minor) * p0_1 * alphas[i+0];
        results[i+1] =  (alphas[i+1] - pi_major - pi_minor) *
            (alphas[i+1] + pi_major + pi_minor) * p1_1 * alphas[i+1];
        results[i+2] =  (alphas[i+2] - pi_major - pi_minor) *
            (alphas[i+2] + pi_major + pi_minor) * p2_1 * alphas[i+2];
        results[i+3] =  (alphas[i+3] - pi_major - pi_minor) *
            (alphas[i+3] + pi_major + pi_minor) * p3_1 * alphas[i+3];
    }

    tscheb_dsines(alphas+i, N-i, results+i); // do the rest
}

double tscheb_dsine(double alpha, bool angle_is_normalized) {
    if (!angle_is_normalized) { alpha = pi_to_pi(alpha); }
    // Make sure alpha is in (-PI, PI]
        float coeffs[] = {
        -0.10132118,          // x
         0.0066208798,        // x^3
        -0.00017350505,       // x^5
         0.0000025222919,     // x^7
        -0.000000023317787,   // x^9
         0.00000000013291342, // x^11
    };
    double pi_major = 3.1415927;
    double pi_minor = -0.00000008742278;
    double alpha2 = alpha*alpha;
    double p11 = coeffs[5];
    double p9  = p11*alpha2 + coeffs[4];
    double p7  = p9*alpha2  + coeffs[3];
    double p5  = p7*alpha2  + coeffs[2];
    double p3  = p5*alpha2  + coeffs[1];
    double p1  = p3*alpha2  + coeffs[0];
    return (alpha - pi_major - pi_minor) *
    (alpha + pi_major + pi_minor) * p1 * alpha;
}

float tscheb_fsine(float alpha, bool angle_is_normalized) {
    if (! angle_is_normalized) { alpha = pi_to_pi(alpha); }

    float coeffs[] = {
        -0.10132118,          // x
         0.0066208798,        // x^3
        -0.00017350505,       // x^5
         0.0000025222919,     // x^7
        -0.000000023317787,   // x^9
         0.00000000013291342, // x^11
    };
    float pi_major = 3.1415927;
    float pi_minor = -0.00000008742278;
    float alpha2 = alpha*alpha;
    float p11 = coeffs[5];
    float p9  = p11*alpha2 + coeffs[4];
    float p7  = p9*alpha2  + coeffs[3];
    float p5  = p7*alpha2  + coeffs[2];
    float p3  = p5*alpha2  + coeffs[1];
    float p1  = p3*alpha2  + coeffs[0];
    return (alpha - pi_major - pi_minor) *
    (alpha + pi_major + pi_minor) * p1 * alpha;
}