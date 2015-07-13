/* 
 * JF Nested Dielectric by Jonah Friedman
 * 1.0.5
 * Copyright (c) 2014, Psyop Media Company, LLC and Jonah Friedman 
 * Open sourced under the 3 clause BSD license, see license.txt
 */



#if ( AI_VERSION_ARCH_NUM * 100 + AI_VERSION_MAJOR_NUM ) < 401
	#define SAMPLETYPE double
#else
	#define SAMPLETYPE float
#endif



// ---------------------------------------------------//
// - Enumerations 
// ---------------------------------------------------//

// - note: there are more enumerations related to spectral stuff in spectral_distribution.h

const char *  enum_brdfs[] =
{
	"stretched_phong",
	"cook_torrance",
	"ward_rayTangent",
	"ward_userTangent",
	NULL
};

enum brdfs
{
	b_stretched_phong,
	b_cook_torrance,
	b_ward_rayTangent,
	b_ward_userTangent,
};

const char *  enum_shadow_modes[] =
{
	"black",
	"transmit_only",
	"transmit_and_outer_fresnels",
	"transmit_and_all_fresnels",
	NULL
};

enum shadow_modes
{
	sh_black,
	sh_transmit_only,
	sh_transmit_and_outer_fresnels,
	sh_transmit_and_all_fresnels,
};

const char *  enum_specular_raytypes[] =
{
	"Glossy",
	"Reflected",
	NULL
};

enum specular_raytypes
{
	rt_glossy,
	rt_reflected,
};


// ---------------------------------------------------//
// - Definitions and data types  
// ---------------------------------------------------//

#define ZERO_EPSILON 0.000000005f

const int max_media_count = 32 + 1;  // media 0 is reserved for the air everything exists in, so to allow 32 media this has to be set to 33. 

typedef struct mediaIntStruct { int v[max_media_count]; } mediaIntStruct;
typedef struct mediaFloatStruct { float v[max_media_count]; } mediaFloatStruct;
typedef struct mediaBoolStruct { bool v[max_media_count]; } mediaBoolStruct;
typedef struct mediaAtColorStruct { AtColor v[max_media_count]; } mediaAtColorStruct;

// ray_ is tracking something to do wtih the ray tree
// caustic_ is tracking caustic behavior
// media_ are arrays about media
// 

typedef struct Ray_State_Datatype {
	bool    ray_monochromatic; // branching ray tree data, should be cached
	float   ray_wavelength;
	int     ray_TIRDepth; // branching ray tree data, should be cached
	int     ray_invalidDepth; // branching ray tree data, should be cached
	AtColor ray_energy; // branching ray tree data, should be cached

	bool caustic_behaviorSet; // branching ray tree data, should be cached 
	int  caustic_mode;
	bool caustic_refractDirect;
	bool caustic_refractIndirect;
	bool caustic_TIR;
	bool caustic_dispersion;
	bool caustic_specDirect;
	bool caustic_specIndirect;
	bool caustic_specInternal;

	mediaIntStruct   media_inside; // branching ray tree data, should be cached
	mediaIntStruct   shadow_media_inside;
	mediaFloatStruct media_iOR;
	mediaBoolStruct  media_disperse;
	mediaFloatStruct media_dispersion;
	mediaIntStruct   media_BTDF;
	mediaIntStruct   media_BRDF;
	
	mediaFloatStruct media_refractDirect;
	mediaFloatStruct media_refractIndirect;
	mediaFloatStruct media_specDirect;
	mediaFloatStruct media_specIndirect;

	mediaFloatStruct media_refractRoughnessU;
	mediaFloatStruct media_refractRoughnessV;
	mediaFloatStruct media_specRoughnessU;
	mediaFloatStruct media_specRoughnessV;

	mediaAtColorStruct media_transmission;

	mediaFloatStruct media_blurAnisotropicPoles;

	spectral_LUT * spectral_LUT_;
	float energy_cutoff;
	bool polarized;
	AtVector polarizationVector;

} Ray_State_Datatype;



typedef struct Ray_State_Cache_Datatype {
	bool    ray_monochromatic;
	int     ray_TIRDepth;
	int     ray_invalidDepth;
	AtColor ray_energy;
	bool caustic_behaviorSet; 
	mediaIntStruct   media_inside;
} Ray_State_Cache_Datatype;



struct ShaderData{
	AtSampler * dispersion_sampler;
	AtSampler * specular_sampler;
	AtSampler * refraction_sampler;
	AtSampler * russian_roullete_single_sampler;
	std::string aov_direct_refraction;
	std::string aov_indirect_refraction;
	std::string aov_direct_specular;
	std::string aov_indirect_specular;
	spectral_LUT spectral_LUT_;
	AtVector polarizationVector; 
};

const int base_sampler_seed = 5;




// ---------------------------------------------------//
// - utilities  
// ---------------------------------------------------//


void uncacheRayState( Ray_State_Datatype * RayState, Ray_State_Cache_Datatype * RayStateCache )
{
	memcpy(&RayState->media_inside, &RayStateCache->media_inside, sizeof(mediaIntStruct) );
	RayState->ray_monochromatic 	= RayStateCache->ray_monochromatic;
	RayState->caustic_behaviorSet 	= RayStateCache->caustic_behaviorSet;
	RayState->ray_TIRDepth 			= RayStateCache->ray_TIRDepth;
	RayState->ray_invalidDepth 		= RayStateCache->ray_invalidDepth;
	RayState->ray_energy 			= RayStateCache->ray_energy;
}

void cacheRayState( Ray_State_Datatype * RayState, Ray_State_Cache_Datatype * RayStateCache )
{
	memcpy(&RayStateCache->media_inside, &RayState->media_inside, sizeof(mediaIntStruct) );
	RayStateCache->ray_monochromatic 	= RayState->ray_monochromatic;
	RayStateCache->caustic_behaviorSet 	= RayState->caustic_behaviorSet;
	RayStateCache->ray_TIRDepth 		= RayState->ray_TIRDepth;
	RayStateCache->ray_invalidDepth 	= RayState->ray_invalidDepth;
	RayStateCache->ray_energy 			= RayState->ray_energy;
}






float snell (float n1, float n2, float cosThetaI)
{
	// computes snell's law, and returns cos(thetaT), the cosine of the transmitted ray. Cosine is generally more useful than the angle itself. 
	// returns -1.0f if the result is undefined, which means TIR. 
	float aSinThetaT = (sinf( acosf( cosThetaI ) ) * n1) / n2;
	if (aSinThetaT > 0.9999999f)
	{
		return (-1.0f);
	}
	else
	{
		float cosThetaT = cosf(asinf(aSinThetaT));
		return (cosThetaT);
	}
}

float fresnelEquations (float n1, float n2, float cosThetaI, float polarization = 0.5f, bool TIRisRefraction = false )
{
	// computes the fresnel equations, all both of em. Returns reflectance. 
	// n1 is medium IOR, n2 is material IOR, cosThetaI is the dot of the incident ray with the normal, and polarization is a blend between S or P polarization.

	// this fixes some NANs. They happen when the dot product comes out in a rare case that cosThetaI is outside the 1 to -1 range
	cosThetaI = (cosThetaI > 1.0f) ? 1.0f : cosThetaI;
	cosThetaI = (cosThetaI < -1.0f) ? -1.0f : cosThetaI;

	float cosThetaT = snell(n1, n2, cosThetaI);
	
	if (cosThetaT < -0.0f)
	{
		// Total Internal Reflection, reflectance is 1.0
		if (TIRisRefraction)
			return 0.0f;
		else
			return 1.0f;
	}
	else
	{
		// Not TIR, fresnel equations can commence
		const float fresnelS = (pow( (n1 * cosThetaI - n2 * cosThetaT ) / (n1 * cosThetaI + n2 * cosThetaT ), 2)) ;
		const float fresnelP = (pow( (n1 * cosThetaT - n2 * cosThetaI ) / (n1 * cosThetaT + n2 * cosThetaI ), 2)) ;

		const float returnvalue = (fresnelP * polarization) + (fresnelS * (1-polarization));

		if (returnvalue > 1.0f)
			return 1.0f;
		else
			return returnvalue;
	}
}
// char * const mediumMessage2 = "Inside 2 media!.";
// AiMsgWarning(mediumMessage2);

void parallelPark ( AtVector refractionVector, AtShaderGlobals* shaderGlobalsCopy )
{
	// Paralell Parking Maneuver
	AtVector reflectionSourceVector;
	refractionVector = -refractionVector ;
	AiReflect( &refractionVector, &shaderGlobalsCopy->N, &reflectionSourceVector ) ;
	reflectionSourceVector = -AiV3Normalize( reflectionSourceVector ) ;
	shaderGlobalsCopy->Nf = -shaderGlobalsCopy->Nf ;  // set the shading normal to be the regular normal, in order to do specs from behind things
	shaderGlobalsCopy->Ngf = -shaderGlobalsCopy->Ngf ;
	shaderGlobalsCopy->Rd = reflectionSourceVector ; // Set the ray direction to be the backwards refracted ray
}

float russianRoulette( AtSamplerIterator* rrSamplerIterator, float value, float probability)
{		   
	float outvalue = 0.0f;

	if (probability >= 1.0 || probability <= 0.0)
		return value;

	// Prepare a sampler to get a value to compare against.
	SAMPLETYPE sample[2];
	const bool foundSample = AiSamplerGetSample( rrSamplerIterator, sample );
	float fsample = (float) sample[0];
	
	while ( AiSamplerGetSample( rrSamplerIterator, sample ) ) 
	{
		// Exhaust the sampler. Good night, sampler. 
		// This seems necessary, having the unexhausted sampler caused some problems with the specular sampler. 
	}

	if (fsample < probability) 
	{
		// for example, if there's a 0.33 chance of tracing, we trace at 3x value (in the 0.33 chance that we trace at all).
		outvalue = value / probability;
	}
	else 
	{
		outvalue = 0.0f;
	}	

	return outvalue;
}



float refractiveRoughness( float roughness1, float roughness2, float n1, float n2)
{	
	float averageRoughness = (roughness1 + roughness2) / 2.0f ;

	float largerIOR = (n1 > n2) ? n1 : n2 ;
	float smallerIOR = (n1 < n2) ? n1 : n2 ; 

	float IORDifference = (largerIOR / smallerIOR) - 1.0f ;

	return ( averageRoughness * IORDifference );
}

float refractRoughnessConvert(float roughness)
{
	// squaring the roughness has been experimentally shown to do a pretty good job of matching 
	// the hijacked brdf roughnesses to the actual microfacet BTDF roughness.
	// So if you want roughness .5, you feed the hijacked BRDF roughness 0.5 ^ 2, which does a pretty good job. 
	
	return roughness * roughness;
	//return roughness;
}


void blurAnisotropicPoles( float *  roughnessU, float * roughnessV, float * blurAmount, AtVector * normal, AtVector * raw_tangent )
{
	const float averagedRoughness = ( *roughnessU + *roughnessV ) / 2.0f ;
	const float tDotN = AiV3Dot(AiV3Normalize( *raw_tangent ), *normal) ;
	const float blendValue = pow( std::abs( tDotN ),  ( 1.0f / *blurAmount ) - 1.0f ) ;  // power function: ^0 = 1 for fully blurry, ^ high for not blurry at all)

	*roughnessU = ( *roughnessU * ( 1.0f - blendValue) ) + ( averagedRoughness * ( blendValue) );
	*roughnessV = ( *roughnessV * ( 1.0f - blendValue) ) + ( averagedRoughness * ( blendValue) );

	return;
}



AtColor transmissionColor( AtColor * transmission, float depth)
{
	if (*transmission != AI_RGB_WHITE)
	{		
		return AiColorCreate( 
			pow( transmission->r, depth ), 
			pow( transmission->g, depth ), 
			pow( transmission->b, depth ) );
	}
	else
	{
		return AI_RGB_WHITE;
	}
}

AtColor transmissionOnSample( AtColor * transmission, AtScrSample * sample, bool traceHit )
{
	AtColor returnEnergy;

	const float depth = (float) sample->z;
	if (traceHit)
	{
		returnEnergy = transmissionColor( transmission, depth);
	}
	else
	{
		returnEnergy = AiColorCreate( 
			transmission->r >= 1.0f ? 1.0f: 0.0f, 
			transmission->g >= 1.0f ? 1.0f: 0.0f,
			transmission->b >= 1.0f ? 1.0f: 0.0f );	
	}

	return returnEnergy;
}


void updateMediaInsideLists(int m_cMatID, bool entering, mediaIntStruct * media_inside_list, bool reverse = false)
{
	// If we're entering an object, increment, if we're leaving, decrement. 
	// Reverse reverses the behavior. 

	if ( !reverse)
		media_inside_list->v[m_cMatID] = media_inside_list->v[m_cMatID] + (entering ? 1 : -1 );
	else
		media_inside_list->v[m_cMatID] = media_inside_list->v[m_cMatID] + (entering ? -1 : 1 );
}
