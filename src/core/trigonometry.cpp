#include "trigonometry.h"
#include "linalg.h"
const double step_size = M_PI / steps_pi;

const double num_steps = steps_pi / M_PI;

const __m256d step_size_vec = _mm256_set1_pd(step_size);
const __m256d num_steps_vec = _mm256_set1_pd(num_steps);
const __m256d minus_ones = _mm256_set1_pd((double)-1.0);
const __m256d pi_2_vec = _mm256_set1_pd((double)M_PI_2);

const __m256d offset_sin2 = _mm256_set1_pd((double)n_2pi * M_PI * num_steps);

double sin_table[N_sin];
double sin2_table[N_sin];

// One-sided sine table [- n_2pi * M_PI, + n_2pi * M_PI]
//double sin_table[N_sin]; //Here we only store positive sin for angles >=0


// Symmetric sine table [- n_2pi * M_PI, + n_2pi * M_PI]
//double sin2_table[N_sin];

void init_sin() {
    for (int i = 0; i < N_sin; i++) {
        sin_table[i] = sin(((double)i/steps_pi)*M_PI);
    }
}

void init_sin2() {
    //symmetric
    double low = - (double)n_2pi * M_PI;
    for (int i = 0; i < N_sin; i++) {
        sin2_table[i] = std::sin((double)low + ((double)i/steps_pi)*M_PI);
    }
}



double read_sin(double angle) {
    int idx = (int)(angle / M_PI * steps_pi);
    if (angle >= 0) {
        return sin_table[idx];
    } else {
        return -sin_table[-idx];
    }
}




//! ERROR in 0.x range (rad)!
float atan2_approximation1(float y, float x)
{
    //http://pubs.opengroup.org/onlinepubs/009695399/functions/atan2.html
    //Volkan SALMA

    const float ONEQTR_PI = M_PI / 4.0;
	const float THRQTR_PI = 3.0 * M_PI / 4.0;
	float r, angle;
	float abs_y = fabs(y) + 1e-10f;      // kludge to prevent 0/0 condition
	if ( x < 0.0f )
	{
		r = (x + abs_y) / (abs_y - x);
		angle = THRQTR_PI;
	}
	else
	{
		r = (x - abs_y) / (x + abs_y);
		angle = ONEQTR_PI;
	}
	angle += (0.1963f * r * r - 0.9817f) * r;
	if ( y < 0.0f )
		return( -angle );     // negate if in quad III or IV
	else
		return( angle );


}

//! Error in 0.00x range (appr)
float sqrt3(const float& n)
{
   static union {int i; float f;} u;
   u.i = 0x2035AD0C + (*(int*)&n >> 1);
   return n / u.f + u.f * 0.25f;
}

inline float rsqrtss_times_x( float* x)
{
   float out;
   __m128 in = _mm_load_ss( x );
   _mm_store_ss( &out, _mm_mul_ss( in, _mm_rsqrt_ss( in ) ) );
   // compiles to movss, movaps, rsqrtss, mulss, movss
   return out;
}

inline void rsqrtss_times_x_void( float* x, float* out)
{
   __m128 in = _mm_load_ss( x );
   _mm_store_ss( out, _mm_mul_ss( in, _mm_rsqrt_ss( in ) ) );
   // compiles to movss, movaps, rsqrtss, mulss, movss
}

__m256d read_sin_vec(__m256d angle) {
    angle = _mm256_mul_pd(num_steps_vec, angle);
    __m128i index = _mm256_cvtpd_epi32(angle);

    //Alternative:
    __m128i all_pos_index = _mm_sign_epi32(index, index);


    auto results = _mm256_i32gather_pd(sin_table, all_pos_index , 8);
    auto neg_results = _mm256_mul_pd(minus_ones, results);
    return _mm256_blendv_pd(results, neg_results,angle);
}

__m256d read_cos_vec(__m256d angle) {
    //This sub and the mul in read_sin_vec can be fused to an fmsub, but I keep it modular for now
    // In cos we can also use symmetry!
    angle = _mm256_add_pd(angle, pi_2_vec);
    return read_sin_vec(angle);
}

__m256d read_sin2_vec(__m256d angle) {
    angle = _mm256_fmadd_pd(num_steps_vec, angle, offset_sin2);
    __m128i index = _mm256_cvtpd_epi32(angle);
    return _mm256_i32gather_pd(sin2_table, index, 8);
}

__m256d read_cos2_vec(__m256d angle) {
    angle = _mm256_add_pd(angle, pi_2_vec);
    return read_sin2_vec(angle);
}

double read_cos(double angle) {
    return read_sin(angle + M_PI_2);
}


double Q_rsqrt( double number )
{
	long i;
	double x2, y;
	const double threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;                       // evil floating point bit level hacking
	i  = 0x5f3759df - ( i >> 1 );               // what the fuck? 
	y  = * ( double * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
//	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

	return y;
}