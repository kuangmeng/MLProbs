__constant const float LOG_ZERO = -2e20f;
__constant const float LOG_ONE = 0.0f;
__constant const float LOG_UNDERFLOW_THRESHOLD = LOG_UNDERFLOW;

#ifdef USE_TAYLOR_HMM
	#define fast_exp(x)		array_exp(x)
	#define log_exp_1(x)	array_log_exp_1(x)
#else
	#define fast_exp(x)		native_exp(x)
	#define log_exp_1(x)	native_log(native_exp(x) + 1)
#endif


#ifdef USE_TAYLOR_HMM

typedef double4 exp_vector_type;
typedef double exp_scalar_type;

__constant const exp_vector_type exp_05_0	= (exp_vector_type)(0.03254409303190190000, 0.16280432765779600000, 0.49929760485974900000, 0.99995149601363700000);
__constant const exp_vector_type exp_1_05	= (exp_vector_type)(0.01973899026052090000, 0.13822379685007000000, 0.48056651562365000000, 0.99326940370383500000);
__constant const exp_vector_type exp_2_1	= (exp_vector_type)(0.00940528203591384000, 0.09414963667859410000, 0.40825793595877300000, 0.93933625499130400000);
__constant const exp_vector_type exp_4_2	= (exp_vector_type)(0.00217245711583303000, 0.03484829428350620000, 0.22118199801337800000, 0.67049462206469500000);
__constant const exp_vector_type exp_8_4	= (exp_vector_type)(0.00012398771025456900, 0.00349155785951272000, 0.03727721426017900000, 0.17974997741536900000);
__constant const exp_vector_type exp_16_8	= (exp_vector_type)(0.00000051741713416603, 0.00002721456879608080, 0.00053418601865636800, 0.00464101989351936000);
__constant exp_scalar_type free_05_0	= 0.99999925508501600000;
__constant exp_scalar_type free_1_05	= 0.99906756856399500000;
__constant exp_scalar_type free_2_1		= 0.98369508190545300000;
__constant exp_scalar_type free_4_2		= 0.83556950223398500000;
__constant exp_scalar_type free_8_4		= 0.33249299994217400000;
__constant exp_scalar_type free_16_8	= 0.01507447981459420000;

float array_exp(float x)
{
	if (x > 0) return native_exp(x);
	if (x <= -16) return 0;
	
	exp_vector_type coeffs = (x > -2) 
		? ( (x > -0.5) ? exp_05_0 : ( (x > -1) ? exp_1_05 : exp_2_1) )
		: ( (x > -8) ? ( (x > -4) ? exp_4_2 : exp_8_4 ) : exp_16_8 );

	exp_scalar_type free = (x > -2) 
		? ( (x > -0.5) ? free_05_0 : ( (x > -1) ? free_1_05 : free_2_1) )
		: ( (x > -8) ? ( (x > -4) ? free_4_2 : free_8_4 ) : free_16_8 );

	return (((coeffs.s0 * x + coeffs.s1) * x + coeffs.s2) * x + coeffs.s3) * x + free;
}

__constant const float4 less_1 = (float4)(-0.009350833524763f, 0.130659527668286f, 0.498799810682272f, 0.693203116424741f);
__constant const float4 less_2 = (float4)(-0.014532321752540f, 0.139942324101744f, 0.495635523139337f, 0.692140569840976f);
__constant const float4 less_4 = (float4)(-0.004605031767994f, 0.063427417320019f, 0.695956496475118f, 0.514272634594009f);
__constant const float4 more_4 = (float4)(-0.000458661602210f, 0.009695946122598f, 0.930734667215156f, 0.168037164329057f);
	
float array_log_exp_1(float x) 
{
	float4 coeffs = (float4) (
		((x) <= 2.50f) 
		? (((x) <= 1.00f) ? less_1 : less_2) 
		: (((x) <= 4.50f) ? less_4 : more_4) );

	return mad(mad(mad(coeffs.s0, x, coeffs.s1), x, coeffs.s2), x, coeffs.s3);
}
#endif


/// Adds two log probabilities and stores in the first argument (local memory variant)
void logOfSum_local(local float *x, float y){
	float d = *x;
	float s;
	if (d < y) { s = d; }
	else { s = y; y = d; }
	d = y - s;
	*x = (s <= LOG_ZERO || d >= LOG_UNDERFLOW_THRESHOLD) ? y : log_exp_1(d) + s;
}

/// Adds two log probabilities and stores in the first argument (private memory variant)
void logOfSum_private(private float *x, float y){
	#define d *x
	float s;
	if (d < y) { s = d; }
	else { s = y; y = d; }
	d = y - s;
	*x = (s <= LOG_ZERO || d >= LOG_UNDERFLOW_THRESHOLD) ? y : log_exp_1(d) + s;
	#undef d
}

/// Adds two log probabilities and returns value
float logOfSum(float x, float y)
{
	#define d x
	float s;
	if (d < y) { s = d; }
	else { s = y; y = d; }
	d = y - s;
	return (s <= LOG_ZERO || d >= LOG_UNDERFLOW_THRESHOLD) ? y : log_exp_1(d) + s;
	#undef d
}





