#include "tscheb_sine.h"
#include "pi_to_pi.h"

// assumes angle is normalized
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

double tscheb_dsine(double alpha, bool angle_is_normalized) {
    if (!angle_is_normalized) { alpha = pi_to_pi_nongeneral(alpha); }
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