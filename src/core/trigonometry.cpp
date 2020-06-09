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
        sin2_table[i] = sin((double)low + ((double)i/steps_pi)*M_PI);
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
	//https://gist.github.com/volkansalma/2972237

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



__m256d atan2_approximation2(__m256d y, __m256d x){

	const float ONEQTR_PI = M_PI / 4.0;
	const float THRQTR_PI = 3.0 * M_PI / 4.0;

	__m256d constant = _mm256_set1_pd(1e-10f);
	// absolute value with mask
	//__m256d sign_bit = _mm256_set1_ps(-0.0f);
	//__m256d abs_y = _mm256_andnot_ps(sign_bit, y);
	// simple absolute value
	__m256d abs_y = _mm256_sqrt_pd( _mm256_mul_pd (y, y) );
	__m256d abs_y_simd = _mm256_add_pd(abs_y, constant);
	
	__m256d ONEQTR_PI_simd = _mm256_set_pd(ONEQTR_PI, ONEQTR_PI, ONEQTR_PI, ONEQTR_PI);
	__m256d THRQTR_PI_simd = _mm256_set_pd(THRQTR_PI, THRQTR_PI, THRQTR_PI, THRQTR_PI);

	__m256d zeros = _mm256_set1_pd(0.0f);
	__m256d vmask = _mm256_cmp_pd(x, zeros, _CMP_LT_OQ);

	__m256d angle = _mm256_blendv_pd(ONEQTR_PI_simd, THRQTR_PI_simd, vmask);
	__m256d r1 = _mm256_div_pd(_mm256_add_pd(x, abs_y_simd), _mm256_sub_pd(abs_y_simd, x));
	__m256d r2 = _mm256_div_pd(_mm256_sub_pd(x, abs_y_simd), _mm256_add_pd(abs_y_simd, x));

	__m256d sel_r = _mm256_blendv_pd(r2, r1, vmask);

	__m256d scalar = _mm256_set1_pd(0.1963f);
	__m256d scalar2 = _mm256_set1_pd(0.9817f);
	__m256d temp = _mm256_sub_pd(_mm256_mul_pd(scalar, _mm256_mul_pd(sel_r, sel_r)), scalar2);
	angle = _mm256_add_pd(angle, _mm256_mul_pd(temp, sel_r));

	__m256d vmask2 = _mm256_cmp_pd(y, zeros, _CMP_LT_OQ);
	__m256d angle_neg = _mm256_sub_pd(_mm256_set1_pd(0.0f), angle); 
	__m256d result = _mm256_blendv_pd(angle, angle_neg, vmask2);

	return result;
}

__m256d atan2_approximation3(__m256d y, __m256d x){
    // https://github.com/to-miz/sse_mathfun_extension/blob/master/sse_mathfun_extension.h?fbclid=IwAR1fAfFgZJldewkYYknL2E3JxamBC6KDYSZcY_kvQHqylpEFx6ApZLPm6EM
	/*v4sf x_eq_0 = _mm_cmpeq_ps( x, *(v4sf*)_ps_0 );
	v4sf x_gt_0 = _mm_cmpgt_ps( x, *(v4sf*)_ps_0 );
	v4sf x_le_0 = _mm_cmple_ps( x, *(v4sf*)_ps_0 );
	v4sf y_eq_0 = _mm_cmpeq_ps( y, *(v4sf*)_ps_0 );
	v4sf x_lt_0 = _mm_cmplt_ps( x, *(v4sf*)_ps_0 );
	v4sf y_lt_0 = _mm_cmplt_ps( y, *(v4sf*)_ps_0 );

	v4sf zero_mask = _mm_and_ps( x_eq_0, y_eq_0 );
	v4sf zero_mask_other_case = _mm_and_ps( y_eq_0, x_gt_0 );
	zero_mask = _mm_or_ps( zero_mask, zero_mask_other_case );

	v4sf pio2_mask = _mm_andnot_ps( y_eq_0, x_eq_0 );
	v4sf pio2_mask_sign = _mm_and_ps( y_lt_0, *(v4sf*)_ps_sign_mask );
	v4sf pio2_result = *(v4sf*)_ps_cephes_PIO2F;
	pio2_result = _mm_xor_ps( pio2_result, pio2_mask_sign );
	pio2_result = _mm_and_ps( pio2_mask, pio2_result );

	v4sf pi_mask = _mm_and_ps( y_eq_0, x_le_0 );
	v4sf pi = *(v4sf*)_ps_cephes_PIF;
	v4sf pi_result = _mm_and_ps( pi_mask, pi );

	v4sf swap_sign_mask_offset = _mm_and_ps( x_lt_0, y_lt_0 );
	swap_sign_mask_offset = _mm_and_ps( swap_sign_mask_offset, *(v4sf*)_ps_sign_mask );

	v4sf offset0 = _mm_setzero_ps();
	v4sf offset1 = *(v4sf*)_ps_cephes_PIF;
	offset1 = _mm_xor_ps( offset1, swap_sign_mask_offset );

	v4sf offset = _mm_andnot_ps( x_lt_0, offset0 );
	offset = _mm_and_ps( x_lt_0, offset1 );

	v4sf arg = _mm_div_ps( y, x );
	v4sf atan_result = atan_ps( arg );
	atan_result = _mm_add_ps( atan_result, offset );

	// select between zero_result, pio2_result and atan_result

	v4sf result = _mm_andnot_ps( zero_mask, pio2_result );
	atan_result = _mm_andnot_ps( pio2_mask, atan_result );
	atan_result = _mm_andnot_ps( pio2_mask, atan_result);
	result = _mm_or_ps( result, atan_result );
	result = _mm_or_ps( result, pi_result );*/
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

#ifdef __AVX2__
__m256d read_sin_vec(__m256d angle) {
    angle = _mm256_mul_pd(num_steps_vec, angle);
    __m128i index = _mm256_cvtpd_epi32(angle);

    //Alternative:
    __m128i all_pos_index = _mm_sign_epi32(index, index);


    auto results = _mm256_i32gather_pd(sin_table, all_pos_index , 8);
    auto neg_results = _mm256_mul_pd(minus_ones, results);
    return _mm256_blendv_pd(results, neg_results,angle);
}
#endif

#ifdef __AVX2__
__m256d read_cos_vec(__m256d angle) {
    //This sub and the mul in read_sin_vec can be fused to an fmsub, but I keep it modular for now
    // In cos we can also use symmetry!
    angle = _mm256_add_pd(angle, pi_2_vec);
    return read_sin_vec(angle);
}
#endif

#ifdef __AVX2__
__m256d read_sin2_vec(__m256d angle) {
#ifdef __FMA__
    angle = _mm256_fmadd_pd(num_steps_vec, angle, offset_sin2);
#else
    angle = _mm256_add_pd( _mm256_mul_pd(num_steps_vec, angle), offset_sin2);
#endif
    __m128i index = _mm256_cvtpd_epi32(angle);
    return _mm256_i32gather_pd(sin2_table, index, 8);
}
#endif

#ifdef __AVX2__
__m256d read_cos2_vec(__m256d angle) {
    angle = _mm256_add_pd(angle, pi_2_vec);
    return read_sin2_vec(angle);
}
#endif

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
