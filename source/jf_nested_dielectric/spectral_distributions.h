/* 
 * JF Nested Dielectric by Jonah Friedman
 * 1.0.4
 * Copyright (c) 2014, Psyop Media Company, LLC and Jonah Friedman 
 * Open sourced under the 3 clause BSD license, see license.txt
 */

#include <math.h>

const int standard_observer_point_num = 42;
typedef struct CIE_standard_observer_data
{
	int num_indices;
	float reference_wavelengths[standard_observer_point_num] ;
	float X_curve[standard_observer_point_num] ;
	float Y_curve[standard_observer_point_num] ;
	float Z_curve[standard_observer_point_num] ;
	float start_wavelength ;
	float end_wavelength ;
} CIE_standard_observer_data;

const CIE_standard_observer_data CIE_data = 
{
	standard_observer_point_num,
	{370.0f, 380.0f, 390.0f, 400.0f, 410.0f, 420.0f, 430.0f, 440.0f, 450.0f, 460.0f, 470.0f, 480.0f, 490.0f, 500.0f, 510.0f, 520.0f, 530.0f, 540.0f, 550.0f, 560.0f, 570.0f, 580.0f, 590.0f, 600.0f, 610.0f, 620.0f, 630.0f, 640.0f, 650.0f, 660.0f, 670.0f, 680.0f, 690.0f, 700.0f, 710.0f, 720.0f, 730.0f, 740.0f, 750.0f, 760.0f, 770.0f, 780.0f},
	{0.0f, 0.0014f, 0.0042f, 0.0143f, 0.0435f, 0.1344f, 0.2839f, 0.3483f, 0.3362f, 0.2908f, 0.1954f, 0.0956f, 0.032f, 0.0049f, 0.0093f, 0.0633f, 0.1655f, 0.2904f, 0.4334f, 0.5945f, 0.7621f, 0.9163f, 1.0263f, 1.0622f, 1.0026f, 0.8544f, 0.6424f, 0.4479f, 0.2835f, 0.1649f, 0.0874f, 0.0468f, 0.0227f, 0.0114f, 0.0058f, 0.0029f, 0.0014f, 0.0007f, 0.0003f, 0.0002f, 0.0001f, 0.0f },
	{0.0f, 0.0f, 0.0001f, 0.0004f, 0.0012f, 0.004f, 0.0116f, 0.023f, 0.038f, 0.06f, 0.091f, 0.139f, 0.208f, 0.323f, 0.503f, 0.71f, 0.862f, 0.954f, 0.995f, 0.995f, 0.952f, 0.87f, 0.757f, 0.631f, 0.503f, 0.381f, 0.265f, 0.175f, 0.107f, 0.061f, 0.032f, 0.017f, 0.0082f, 0.0041f, 0.0021f, 0.001f, 0.0005f, 0.0003f, 0.0001f, 0.0001f, 0.0f, 0.0f},
	{0.0f, 0.0065f, 0.0201f, 0.0679f, 0.2074f, 0.6456f, 1.3856f, 1.7471f, 1.7721f, 1.6692f, 1.2876f, 0.813f, 0.4652f, 0.272f, 0.1582f, 0.0782f, 0.0422f, 0.0203f, 0.0087f, 0.0039f, 0.0021f, 0.0017f, 0.0011f, 0.0008f, 0.0003f, 0.0002f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
	370.0f,
	780.0f,
};

const int spectral_lookup_point_num = 128;
typedef struct spectral_wavelength_LUT { float values[ spectral_lookup_point_num ]; } spectral_wavelength_LUT;
typedef struct spectral_color_LUT { AtColor values[ spectral_lookup_point_num ]; } spectral_color_LUT;
typedef struct spectral_LUT
{
	spectral_wavelength_LUT wavelengths;
	spectral_color_LUT colors;
	int num_indices;
} spectral_LUT ;

typedef struct gamut_definition {
	AtMatrix matrix;
} gamut_definition ;


const char * enum_gamuts[] = 
{
	"sRGB",
	"AdobeRGB",
	"WideGamut",
	"ProPhotoRGB",
	"ACES",
	NULL
};

enum gamuts
{
	g_sRGB,
	g_AdobeRGB,
	g_WideGamut,
	g_ProPhotoRGB,
	g_ACES
};

// Untransposed matrices, so this is the matrix to convert sRGB to XYZ. The matrix is used in transposed form. 

const gamut_definition sRGB = 
{
	{
		3.2406f, 	-1.5372f, 	-0.4986f, 	0.0f,
		-0.9689f, 	1.8758f, 	0.0415f, 	0.0f,
		0.0557f, 	-0.204f, 	1.057f, 	0.0f,
		0.00000f, 	0.0000f, 	0.000f, 	1.0f
	}
};

const gamut_definition AdobeRGB = 
{
	{
		 2.041369f, 	-0.564946f, 	-0.344694f, 	0.0f,
		-0.969266f, 	 1.876011f, 	 0.041556f, 	0.0f,
		 0.013447f, 	 -0.11839f, 	  1.01541f, 	0.0f,
		      0.0f, 	      0.0f,  	      0.0f, 	1.0f
	}
};

const gamut_definition WideGamut = 
{
	{
		   1.6119f, 	-0.202822f, 	-0.302324f, 	0.0f,
		-0.509045f, 	  1.41188f, 	0.0660681f, 	0.0f,
		0.0260843f, 	-0.072347f, 	 0.962016f, 	0.0f,
		      0.0f, 	      0.0f, 	      0.0f, 	1.0f
	 }
};

const gamut_definition ProPhotoRGB = 
{
	{
		     1.39046f, 	 -0.264061f,	-0.0528021f,	0.0f,
		   -0.537655f, 	   1.48894f,	 0.0202733f,	0.0f,
		3.99943e-008f, -7.5953e-009f,	  0.918344f,	0.0f,
		         0.0f,          0.0f,	       0.0f,	1.0f
	}
};

const gamut_definition ACES = 
{
	{
		   1.175304441f,	-0.05598231868f,	-0.0128661219f, 	0.0f,
		 -0.3155550364f,	 1.25754023200f,	0.00899780455f, 	0.0f,
		-0.01940979856f,	 0.01898161422f,  	 0.9200991843f, 	0.0f,
		           0.0f, 	           0.0f,	          0.0f, 	1.0f
	}
};

const int distrib_point_num = 24;
typedef struct spectral_distrib { float values[ distrib_point_num ]; } spectral_distrib;

const char * enum_spectra[] = 
{
	"full",
	"daylight",
	"bluesky",
	"sunset",
	"LED",
	"flourescent",
	"sodium",
	NULL
};

enum spectra
{
	s_full,
	s_daylight,
	s_bluesky,
	s_sunset,
	s_LED,
	s_flourescent,
	s_sodium,
	s_test_spectrum,
};

const spectral_distrib full_spectrum = 
{
	{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f }
};
const spectral_distrib daylight_spectrum  = 
{
	{0.0f, 0.0f, 0.0f, 0.105724794711f, 0.514795077933f, 1.08665115291f, 1.76946007826f, 0.0416481841844f, 1.98951013066f, 1.80431517384f, 1.4753563571f, 2.18218125506f, 2.27387164259f, 3.02988831507f, 2.90409854811f, 3.16381911714f, 2.85982445729f, 2.20256856293f, 1.56304551484f, 1.66894343093f, 1.15405232944f, 0.183602213967f, 0.621613145473f, 0.0f}
};
const spectral_distrib bluesky_spectrum  = 
{
	{ 0.0f, 0.0f, 0.581244438946f, 1.04015058465f, 0.836351024333f, 2.7190010064f, 4.14577096723f, 4.75377359433f, 4.70900017045f, 4.09299445034f, 4.06596644656f, 4.32552899526f, 1.81264762271f, 1.13618612959f, 0.687774900221f, 0.433171337858f, 0.399991773025f, 0.431216834034f, 0.166830575769f, 0.271443193726f, 0.246500096022f, 0.386405723955f, 0.386405723955f, 0.386405723955f }
};
const spectral_distrib sunset_spectrum = 
{
	{ 0.0f, 0.0931405152941f, 0.0931405152941f, 0.200808738182f, 0.270306593172f, 0.394798166635f, 1.10765700895f, 1.23214861999f, 1.27902172049f, 0.529863153277f, 1.93457314145f, 2.21227481431f, 2.19156463352f, 1.76598318233f, 4.10407745505f, 1.13906280562f, 4.1281669099f, 2.49399724575f, 4.65709755023f, 2.74298031754f, 2.18448539022f, 3.60452364743f, 3.29324993618f, 3.58141160579f }
};
const spectral_distrib LED_spectrum = 
{
	{ 0.0f, 0.375236211642f, 0.375236211642f, 0.69467822458f, 3.38446533473f, 0.663320371622f, 7.12646340456f, 4.73564977113f, 1.36917642849f, 1.58630467417f, 1.75247935507f, 5.16742653807f, 1.99316574821f, 1.99316574821f, 1.99316574821f, 1.99316574821f, 1.76384134639f, 1.50442096537f, 6.62517484146f, 6.38364862305f, 0.52784190245f, 0.286315975595f, 0.0f, 0.0f }
};
const spectral_distrib flourescent_earthbulb_spectrum = 
{
	{ 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 5.30869678966f, 0.0f, 2.26210501546f, 0.0f, 0.0f, 7.21745470099f, 0.0f, 0.600013773541f, 0.600013773541f, 1.28209667816f, 0.6116218031f, 5.53848670045f, 0.0f, 0.186742162669f, 0.186742162669f, 0.186742162669f, 0.0f, 0.0f, 0.0f }
};
const spectral_distrib sodium_spectrum  = 
{
	{ 0.0f, 0.0f, 0.0f, 0.0f, 0.0810187457709f, 0.0810187457709f, 0.0810187457709f, 0.0540247784866f, 0.0540247784866f, 0.0540247784866f, 0.0540247784866f, 0.148282970918f, 0.202152090537f, 0.0f, 5.24820173173f, 0.0f, 0.0f, 0.0f, 0.0f, 0.72166964264f, 0.160129629301f, 0.0532031411811f, 0.0933312939392f, 0.0f}
};
const spectral_distrib test_spectrum = 
{
	{1.0f, 1.0f, 1.0f, 1.0f, 10.0f, 10.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 10.0f, 10.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 10.0f, 1.0f, 10.0f, 1.0f, 1.0f }
};


const spectral_distrib * getSpectralPreset(int spectral_preset)
{
	if (spectral_preset == s_full) return &full_spectrum;
	else if (spectral_preset == s_daylight) return &daylight_spectrum;
	else if (spectral_preset == s_bluesky) return &bluesky_spectrum;
	else if (spectral_preset == s_sunset) return &sunset_spectrum;
	else if (spectral_preset == s_LED) return &LED_spectrum;
	else if (spectral_preset == s_flourescent) return &flourescent_earthbulb_spectrum;
	else if (spectral_preset == s_sodium) return &sodium_spectrum;
	else if (spectral_preset == s_test_spectrum) return &test_spectrum;
	else return &test_spectrum;	
}




const int subdiv_distrib_point_num = 512;
typedef struct subdiv_spectral_distrib { float values[subdiv_distrib_point_num];} subdiv_spectraL_distrib;







float jfnd_roundf(float value)
{
  return floor(value + 0.5f);
}

float linear_interpolate( float float1, float float2, float weight)
{
	return (float1 * (1.0f - weight) ) + (float2 * weight );
}

float interpolated_distribution_value( const spectral_distrib *distrib, float index )
{
	const float index_fraction = index - floor(index);
	const int index_whole = (int) jfnd_roundf( index - index_fraction );

	const float lower_value = distrib->values[index_whole];
	float upper_value;
	if (index_whole >= distrib_point_num - 1)
	{
		upper_value = distrib->values[index_whole + 1];
	}
	else
	{
		upper_value = lower_value;
	}

	//const float returnvalue = (lower_value * (1.0f - index_fraction) ) + (upper_value * (index_fraction) );
	const float returnvalue = linear_interpolate(lower_value, upper_value, index_fraction);

	return returnvalue;
}


int spectral_LUT_index_from_sample(float sample)
{
	return (int) jfnd_roundf( sample * (float) (spectral_lookup_point_num - 1) );
}

float spectral_LUT_findex_from_sample(float sample)
{
	return sample * (float) (spectral_lookup_point_num - 1) ;
}


spectral_wavelength_LUT build_spectral_wavelength_LUT( const spectral_distrib *distrib )
{
	float total_area = 0;

	subdiv_spectral_distrib subdiv_distrib ;
	subdiv_spectral_distrib subdiv_distrib_iterative ;
	float subdiv_accumulator = 0;
	subdiv_distrib.values[0] = 0.0f;
	subdiv_distrib_iterative.values[0] = 0.0f;
	for( int i = 1; i < subdiv_distrib_point_num; i++ )
	{
		// to do: make this linear interpolation into some sort of a cubic interpolation
		const float index = ( float (i) / float (subdiv_distrib_point_num - 1) ) * (float) (distrib_point_num - 1);
		subdiv_distrib.values[i] = interpolated_distribution_value( distrib, index );
		subdiv_accumulator += subdiv_distrib.values[i];
		subdiv_distrib_iterative.values[i] = subdiv_accumulator;
	}

	spectral_wavelength_LUT wavelength_LUT;

	const float spacing = (float) subdiv_accumulator / (float) (spectral_lookup_point_num - 1) ;

	for( int lut_index = 0; lut_index < spectral_lookup_point_num; lut_index++ )
	{
		const float search_value = spacing * (float) lut_index; // search value is the ideal perfect value to find in the subdiv array

		for( int subdiv_index = 0; subdiv_index < subdiv_distrib_point_num; subdiv_index++ )
		{			
			const float distance = std::abs( search_value - subdiv_distrib_iterative.values[subdiv_index] ); // difference between this index and search value
			const float distance_prev = (subdiv_index != 0) ? std::abs(search_value - subdiv_distrib_iterative.values[subdiv_index - 1]): distance + 10.0f; // difference between previous index and search value
			const float distance_next = (subdiv_index != subdiv_distrib_point_num - 1) ? std::abs(search_value - subdiv_distrib_iterative.values[subdiv_index + 1]) : distance + 10.0f; // difference between next index and search value

			if (distance < distance_next) //next index is worse and not better
			{
				wavelength_LUT.values[lut_index] = linear_interpolate(
					CIE_data.start_wavelength, 
					CIE_data.end_wavelength, 
					(float) subdiv_index / (float) subdiv_distrib_point_num );
				break;
			}
		}
	}
	return wavelength_LUT;
}

float find_wavelength_index( float wavelength )
{
	const float start_wavelength = CIE_data.start_wavelength;
	const float end_wavelength = CIE_data.end_wavelength;

	const float indexFraction = (wavelength - start_wavelength) / (end_wavelength - start_wavelength) ;
	const float index = indexFraction * (float) (standard_observer_point_num - 1);
	if (wavelength > CIE_data.end_wavelength || wavelength < CIE_data.start_wavelength)
		return 5;
	else
		return index;
}

AtColor wavelength_to_color( float wavelength, const gamut_definition * color_gamut, AtColor whiteColor = AI_RGB_WHITE )
{
	float index = find_wavelength_index( wavelength );

	const float index_fraction = fmod(index, 1.0f);
	const int index_whole = int(floor( index ));
	const AtVector xyz_color_1 = {CIE_data.X_curve[index_whole],CIE_data.Y_curve[index_whole], CIE_data.Z_curve[index_whole]};
	const AtVector xyz_color_2 = {CIE_data.X_curve[index_whole + 1],CIE_data.Y_curve[index_whole + 1], CIE_data.Z_curve[index_whole + 1]};;

	AtVector whitepoint = {whiteColor.r, whiteColor.g, whiteColor.b};
	AtVector standardWhite = AI_V3_ONE;
	AiM4VectorByMatrixMult(&standardWhite, color_gamut->matrix, &standardWhite);
	AiM4VectorByMatrixMult(&whitepoint, color_gamut->matrix, &whitepoint);
	AtVector whiteshift = whitepoint/standardWhite;	

	const AtVector xyz_color = AiV3Lerp( index_fraction, xyz_color_1, xyz_color_2 ) * whiteshift;

	AtVector outColor;
	AtColor return_color;	

	AiM4VectorByMatrixTMult(&outColor, color_gamut->matrix, &xyz_color);

	return_color.r = outColor.x;
	return_color.g = outColor.y;
	return_color.b = outColor.z;

	return return_color;
}




/*
 * ---------------------------------------------------
 * - Functions for the shader proper 
 * ---------------------------------------------------
 */

spectral_LUT build_nonuniform_spectral_LUT( int spectral_preset, int gamut_preset, float saturation, bool clamp_negative_colors, AtColor xyz_whiteshift )
{
	const spectral_distrib *distrib;
	distrib = getSpectralPreset(spectral_preset);


	const gamut_definition * gamut_def;
	if (gamut_preset == g_sRGB) gamut_def = &sRGB;
	else if (gamut_preset == g_AdobeRGB) gamut_def = &AdobeRGB;
	else if (gamut_preset == g_WideGamut) gamut_def = &WideGamut;
	else if (gamut_preset == g_ProPhotoRGB) gamut_def = &ProPhotoRGB;
	else if (gamut_preset == g_ACES) gamut_def = &ACES;
	else  gamut_def = &sRGB;

	spectral_LUT out_LUT;
	out_LUT.wavelengths = build_spectral_wavelength_LUT( distrib );
	AtColor colorAccumulator = AI_RGB_BLACK;

	for( int i = 0; i < spectral_lookup_point_num; i++ )
	{
		out_LUT.colors.values[i] = wavelength_to_color( out_LUT.wavelengths.values[i], gamut_def, xyz_whiteshift );
		out_LUT.colors.values[i] = AiColorLerp(saturation, AI_RGB_WHITE, out_LUT.colors.values[i]);
		if (clamp_negative_colors)
		{
			out_LUT.colors.values[i] = AiColorClamp( out_LUT.colors.values[i], 0.0, 10000.0f);		
		}


		colorAccumulator += ( out_LUT.colors.values[i] / (float) spectral_lookup_point_num );
	}
	for( int i = 0; i < spectral_lookup_point_num; i++ )
	{
		// average color in distro will be 1.0, 1.0, 1.0 white	
		out_LUT.colors.values[i] /= colorAccumulator;
	}

	out_LUT.num_indices = spectral_lookup_point_num;
	return out_LUT;
}

void get_interpolated_LUT_value( spectral_LUT * spectral_LUT, float sample, float * wavelength_out, AtColor * color_out )
{
	const float f_index = spectral_LUT_findex_from_sample( sample );
	const float fraction = f_index - floor( f_index );

	const int index1 = (int) floor( f_index );
	const int index2 = index1 + 1;

	const float wavelength1 = spectral_LUT->wavelengths.values[ index1 ];
	const float wavelength2 = spectral_LUT->wavelengths.values[ index2 ];

	const AtColor color1 = spectral_LUT->colors.values[ index1 ];
	const AtColor color2 = spectral_LUT->colors.values[ index2 ];
	
	*wavelength_out = LERP(fraction, wavelength1, wavelength2);
	*color_out = AiColorLerp(fraction, color1, color2);
}

float dispersedIOR(float IOR, float dispersion, float wavelength )
{
	const float start_wavelength = CIE_data.start_wavelength;
	const float end_wavelength = CIE_data.end_wavelength;
	const float indexFraction = (wavelength - start_wavelength) / (end_wavelength - start_wavelength) ;
	
	const float curved_dispersion_value = pow((1.0f - indexFraction), 2.0f);

	return IOR + (curved_dispersion_value - 0.5f) * dispersion;
}
