/* 
 * JF Nested Dielectric by Jonah Friedman
 * 1.0.4
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

const int base_sampler_seed = 5;

struct JFND_Shader_Data{
	AtSampler * dispersion_sampler;
	AtSampler * specular_sampler;
	AtSampler * refraction_sampler;
	AtSampler * russian_roullete_single_sampler;
	std::string aov_direct_refraction;
	std::string aov_indirect_refraction;
	std::string aov_direct_specular;
	std::string aov_indirect_specular;
	spectral_LUT spectral_LUT_ptr;
	AtVector polarizationVector;
	int gloss_samples;
	int refr_samples;


	~JFND_Shader_Data() {
		AiSamplerDestroy(this->dispersion_sampler);
		AiSamplerDestroy(this->specular_sampler);
		AiSamplerDestroy(this->refraction_sampler);
		AiSamplerDestroy(this->russian_roullete_single_sampler);
	}

	JFND_Shader_Data() {
		this->dispersion_sampler = NULL;
		this->specular_sampler = NULL;
		this->refraction_sampler = NULL;
		this->russian_roullete_single_sampler = NULL;
	}

	void update( AtNode* node ) {
		AiSamplerDestroy(this->dispersion_sampler);
		AiSamplerDestroy(this->specular_sampler);
		AiSamplerDestroy(this->refraction_sampler);
		AiSamplerDestroy(this->russian_roullete_single_sampler);

		AtNode* render_options = AiUniverseGetOptions();

		this->refr_samples = AiNodeGetInt(render_options, "GI_refraction_samples");
		this->gloss_samples = AiNodeGetInt(render_options, "GI_glossy_samples");

		this->dispersion_sampler = AiSamplerSeeded( base_sampler_seed + 1, this->refr_samples, 2 );	
		this->specular_sampler = AiSamplerSeeded( base_sampler_seed + 2, this->gloss_samples, 2);
		this->refraction_sampler = AiSamplerSeeded( base_sampler_seed + 3, this->refr_samples, 2 );
		this->russian_roullete_single_sampler = AiSamplerSeeded( base_sampler_seed + 7, 1, 2 );

		this->aov_direct_refraction = AiNodeGetStr(node, "aov_direct_refraction");
		this->aov_indirect_refraction = AiNodeGetStr(node, "aov_indirect_refraction");
		this->aov_direct_specular = AiNodeGetStr(node, "aov_direct_specular");
		this->aov_indirect_specular = AiNodeGetStr(node, "aov_indirect_specular");

		this->refr_samples *= this->refr_samples;
		this->gloss_samples *= this->gloss_samples;

		const bool do_disperse = AiNodeGetBool(node, "disperse");
		if (do_disperse)
		{
			const int spectrum_selection = AiNodeGetInt(node, "spectral_distribution");
			const int gamut_selection = AiNodeGetInt(node, "spectral_gamut");
			const float saturation = AiNodeGetFlt(node, "spectral_saturation");
			const bool clamp = AiNodeGetBool(node, "spectral_clamp_negative_colors");
			this->spectral_LUT_ptr = build_nonuniform_spectral_LUT(spectrum_selection, gamut_selection, saturation, clamp, AI_RGB_WHITE);
		}
		const bool polarize = AiNodeGetBool(node, "polarize");
		if (polarize)
		{
			const float polarizerRotation = AiNodeGetFlt(node, "polarization_angle") * (float) AI_PI;
			AtMatrix cameramatrix;
			AiNodeGetMatrix( AiUniverseGetCamera(), "matrix", cameramatrix );
			AtVector pfilterInCameraSpace;
			AiV3Create(pfilterInCameraSpace, sinf(polarizerRotation), cosf(polarizerRotation), 0.0f);
			AtVector pfilterInWorldSpace;
			AiM4VectorByMatrixMult(&pfilterInWorldSpace, cameramatrix, &pfilterInCameraSpace) ;

			this->polarizationVector = AiV3Normalize(pfilterInWorldSpace);
		}
		else
		{
			this->polarizationVector = AI_V3_ZERO;
		}

	}

};



const int max_media_count = 32 + 1;  // media 0 is reserved for the air everything exists in, so to allow 32 media this has to be set to 33. 

typedef struct mediaIntStruct { int v[max_media_count]; } mediaIntStruct;
typedef struct mediaFloatStruct { float v[max_media_count]; } mediaFloatStruct;
typedef struct mediaBoolStruct { bool v[max_media_count]; } mediaBoolStruct;
typedef struct mediaAtColorStruct { AtColor v[max_media_count]; } mediaAtColorStruct;

// ray_ is tracking something to do wtih the ray tree
// caustic_ is tracking caustic behavior
// media_ are arrays about media
// 



typedef struct Ray_State_Cache_Datatype {
	bool    ray_monochromatic;
	int     ray_TIRDepth;
	int     ray_invalidDepth;
	AtColor ray_energy;
	bool caustic_behaviorSet; 
	mediaIntStruct   media_inside;
} Ray_State_Cache_Datatype;



typedef struct Ray_State_Datatype {
	bool    ray_monochromatic; // branching ray tree data, should be cached
	float   ray_wavelength;
	int     ray_TIRDepth; // branching ray tree data, should be cached
	int     ray_invalidDepth; // branching ray tree data, should be cached
	AtColor ray_energy; // branching ray tree data, should be cached
	AtColor ray_energy_photon; // branching ray tree data, should be cached

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

	spectral_LUT * spectral_LUT_ptr;
	float energy_cutoff;
	bool polarized;
	AtVector polarizationVector;

	void init(JFND_Shader_Data *data, AtShaderGlobals * sg, AtNode * node) {
		this->ray_monochromatic = false;
		this->ray_wavelength = 0.0f;
		this->ray_TIRDepth = 0;
		this->ray_invalidDepth = 0;
		this->ray_energy = AI_RGB_WHITE;
		this->ray_energy_photon = AI_RGB_WHITE;

		// Array initialization - possibly unnecessary except mediaInside? TO DO: investigate
		for( int i = 0; i < max_media_count; i++ )
		{
			this->media_inside.v[i] = 0;
			this->media_iOR.v[i] = 1.0f;
			this->media_disperse.v[i] = false;
			this->media_dispersion.v[i] = 0.0f;
			this->media_BTDF.v[i] = 0;
			this->media_BRDF.v[i] = 0;

			this->media_refractDirect.v[i] = 1.0f;
			this->media_refractIndirect.v[i] = 1.0f;
			this->media_specDirect.v[i] = 1.0f;
			this->media_specIndirect.v[i] = 1.0f;

			this->media_refractRoughnessU.v[i] = 0.0f;
			this->media_refractRoughnessV.v[i] = 0.0f;
			this->media_specRoughnessU.v[i] = 0.0f;
			this->media_specRoughnessV.v[i] = 0.0f;
			this->media_transmission.v[i] = AI_RGB_WHITE;

			this->media_blurAnisotropicPoles.v[i] = 0.0f;
		}

		// caustics behavior values get set once a ray enters a medium from a diffuse ray
		this->caustic_behaviorSet = false;
		this->caustic_mode = 0;
		this->caustic_refractDirect = false;
		this->caustic_refractIndirect = false;
		this->caustic_TIR = false;
		this->caustic_specInternal = false;
		this->caustic_dispersion = false;
		this->caustic_specDirect = false;
		this->caustic_specIndirect = false;

		// Set once values - values set at initialization for all descendent rays
		this->energy_cutoff = pow(10.0f, 16);
		this->polarized = false;
		AiV3Create(this->polarizationVector, 0.0f, 1.0f, 0.0f);
	}
	
	void setEnergyCutoff( float energy_cutoff_exponent)
	{
		this->energy_cutoff = pow(10.0f, energy_cutoff_exponent );
	}

	void setPolarization(bool polarize, AtVector polarizationVector)
	{
		this->polarized = polarize;
		this->polarizationVector = polarizationVector;
	}

	void uncacheRayState(Ray_State_Cache_Datatype * RayStateCache )
	{
		memcpy(&this->media_inside, &RayStateCache->media_inside, sizeof(mediaIntStruct) );
		this->ray_monochromatic 	= RayStateCache->ray_monochromatic;
		this->caustic_behaviorSet 	= RayStateCache->caustic_behaviorSet;
		this->ray_TIRDepth 			= RayStateCache->ray_TIRDepth;
		this->ray_invalidDepth 		= RayStateCache->ray_invalidDepth;
		this->ray_energy 			= RayStateCache->ray_energy;
	}

	void cacheRayState( Ray_State_Cache_Datatype * RayStateCache )
	{
		memcpy(&RayStateCache->media_inside, &this->media_inside, sizeof(mediaIntStruct) );
		RayStateCache->ray_monochromatic 	= this->ray_monochromatic;
		RayStateCache->caustic_behaviorSet 	= this->caustic_behaviorSet;
		RayStateCache->ray_TIRDepth 		= this->ray_TIRDepth;
		RayStateCache->ray_invalidDepth 	= this->ray_invalidDepth;
		RayStateCache->ray_energy 			= this->ray_energy;
	}

	void setMediaIOR(int i, float ior)
	{
		this->media_iOR.v[i] = ior;
	}

	void setMediaDispersion(int i, bool disperse, float dispersion)
	{
		this->media_disperse.v[i] = disperse;
		this->media_dispersion.v[i] = dispersion;
	}
		
	void setMediaSpecular(int i, float direct, float indirect, int BRDF, float roughnessU, 
		float roughnessV)
	{
		this->media_BRDF.v[i] = BRDF;
		this->media_specRoughnessU.v[i] = roughnessU;
		this->media_specDirect.v[i] = direct;
		this->media_specIndirect.v[i] = indirect;
		if (BRDF >= 2)
			this->media_specRoughnessV.v[i] = roughnessV;
		else
			this->media_specRoughnessV.v[i] = 0;
	}

	void setRefractionSettings(int i, float direct, float indirect, int BTDF, float roughnessU, 
		float roughnessV, AtColor transmission, float transmissionScale)
	{
		this->media_refractDirect.v[i] = direct;
		this->media_refractIndirect.v[i] = indirect;
		this->media_BTDF.v[i] = BTDF;
		this->media_refractRoughnessU.v[i] = roughnessU;
		if (BTDF >= 2)
			this->media_refractRoughnessV.v[i] = roughnessV;
		else
			this->media_refractRoughnessV.v[i] = 0.0f;
		
		const float tScale = 1.0f / transmissionScale;
		const AtColor scaledTransmission = 
		{
			pow(transmission.r, tScale),
			pow(transmission.g, tScale),
			pow(transmission.b, tScale)
		};
		this->media_transmission.v[i] = scaledTransmission;
	}

	void setAnisotropySettings(int i, float blurAnisotropicPoles) 
	{
		this->media_blurAnisotropicPoles.v[i] = blurAnisotropicPoles;
	}

} Ray_State_Datatype;



// ---------------------------------------------------//
// - utilities  
// ---------------------------------------------------//


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





bool rayIsPhoton(AtShaderGlobals * sg) {
	bool result = false;
	AiStateGetMsgBool("photon_ray_type", &result);
	return result;
}









typedef struct InterfaceInfo {
	int currentID; // the ID of the material currently being evaluated
	int startingMedium;
	int startingMediumSecondary;
	int m1;
	int m2;
	int m_higherPriority;
	float n1;
	float n2;
	float polarizationTerm;
	AtColor t1;
	AtColor t2;
	bool entering;
	bool validInterface;
	bool mediaExit;
	bool mediaEntrance;

	Ray_State_Datatype * rs;

	InterfaceInfo(Ray_State_Datatype * RayState, int currentID, AtShaderGlobals * sg) 
	{
		this->currentID = currentID;
		this->startingMedium = 0;
		this->startingMediumSecondary = 0;
		this->m1 = 0;
		this->m2 = 0;
		this->m_higherPriority = 0;

		this->n1 = 0.0f;
		this->n2 = 0.0f;

		this->polarizationTerm = 0.5f;

		this->t1 = AI_RGB_WHITE;
		this->t2 = AI_RGB_WHITE;
		this->entering = false;
		this->validInterface = false;
		this->mediaExit = false;
		this->mediaEntrance = false;

		this->rs = RayState;

		bool shadowRay = (sg->Rt == AI_RAY_SHADOW);
		this->setCurrentMediaInfo(shadowRay);
		this->setIntersectionIDs(sg);
		this->setInterfaceProperties();
		this->setPolarizationTerm(sg);
	}


	private: 

	void setCurrentMediaInfo(bool shadowRay) 
	{
		// ---------------------------------------------------//
		// - media logic     
		// 		(determine what media we're in based on mediaInside arrays. )
		// ---------------------------------------------------//

		/*
		 * info.startingMedium is the medium we were already inside when we started tracing
		 * info.startingMediumSecondary is the next medium we will be inside, if it turns out we're leaving this one
		 *
		 * example: we're inside glass and water, glass is the priority. startingMedium is glass, startingMediumSecondary is water.
		 * This is all evaulated without any knowledge of what surface we're evaluating, this is just parsing the mediaInside arrays. 
		 */

		mediaIntStruct * media_inside_ptr = shadowRay ? &this->rs->shadow_media_inside : &this->rs->media_inside;

		for( int i = 0; i < max_media_count; i++ )
		{	
			if ( media_inside_ptr->v[i] > 0 ) // if we find a medium that we're inside
			{
				if ( this->startingMedium == 0 ) // ..and we havent found our current medium yet
				{
					this->startingMedium = i;  // then we've found our medium.
					if (media_inside_ptr->v[i] > 1) // ...and if we're in more than one of these
					{
						this->startingMediumSecondary = i; // then our second medium is the same medium.
						break;
					}
				}
				else if ( this->startingMediumSecondary == 0 )  // ..and we've already found our current medium
				{
					this->startingMediumSecondary = i; // ..then we've found our second medium
					break;
				}
			}
		}
	}

	void setIntersectionIDs(AtShaderGlobals* sg) 
	{
		/*
		 * m = medium ID, n = IOR, T = transmission
		 * 1 is the previous medium, 2 is the next, as seen in n1 and n2 of refraction diagrams
		 * such diagrams: http://en.wikipedia.org/wiki/File:Snells_law.svg
		 */

		this->entering = ( AiV3Dot(sg->Ng, sg->Rd) < 0.0f ) ;
		
		if (this->entering) // (entering a possible medium)
		{
			if (this->startingMedium == 0)
			{
				// entering the first medium, so no media information is available
				this->validInterface = true;
				this->mediaEntrance = true;
				this->m1 = 0;
				this->m2 = this->currentID;
			}
			else if ( this->currentID < this->startingMedium )
			{
				// entering a new highest priority medium
				this->validInterface = true;
				this->m1 = this->startingMedium;
				this->m2 = this->currentID;
			}
			else if ( this->currentID >= this->startingMedium )
			{
				// entering a medium that is not higher priority, interface is invalid
				// or entering the same medium we're already in, interface is also invalid
				this->m1 = this->m2 = this->startingMedium;
			}
		}
		else // (leaving a possible medium)
		{
			if (this->startingMedium == this->currentID)
			{
				// leaving the highest priority medium, which may be 0 (meaning no medium)
				this->validInterface = true;
				this->m1 = this->startingMedium;
				this->m2 = this->startingMediumSecondary;
				if (this->startingMediumSecondary == 0)
				{
					// and there is no second medium, meaning we're leaving all media
					this->mediaExit = true;
				}
			}
			else
			{
				// leaving a medium that is not the current medium, interface invalid
				this->m1 = this->m2;
				this->m1 = this->startingMedium;
			}
		}

		if (this->m1 == 0) 
			this->m_higherPriority = this->m2;
		else if (this->m2 == 0) 
			this->m_higherPriority = this->m1;
		else if (this->m1 < this->m2) 
			this->m_higherPriority = this->m1;
		else if (this->m1 > this->m2) 
			this->m_higherPriority = this->m2;
	}

	void setInterfaceProperties() 
	{
		// - Interface Logic
		// 		(Disperse and multisample)

		if ( this->validInterface )
		{	
			if (this->rs->ray_monochromatic)
			{
				this->n1 = dispersedIOR( this->rs->media_iOR.v[this->m1], this->rs->media_dispersion.v[this->m1], this->rs->ray_wavelength );
				this->n2 = dispersedIOR( this->rs->media_iOR.v[this->m2], this->rs->media_dispersion.v[this->m2], this->rs->ray_wavelength );
			}
			else
			{
				this->n1 = this->rs->media_iOR.v[this->m1];
				this->n2 = this->rs->media_iOR.v[this->m2];
			}


			if (this->n1 == this->n2)
			{			
				this->validInterface = false;
			}
		}

		this->t1 = this->rs->media_transmission.v[this->m1];
		this->t2 = this->rs->media_transmission.v[this->m2];
	}


	void setPolarizationTerm(AtShaderGlobals * sg) 
	{
		if (this->rs->polarized)
		{
			this->polarizationTerm = std::abs( AiV3Dot(sg->Nf, this->rs->polarizationVector ));
		}
		else 
		{
			this->polarizationTerm = 0.5f;
		}
	}

	public:
	bool doBlurryRefraction() 
	{
		bool do_blurryRefraction = false;
		if ( this->rs->media_refractRoughnessU.v[this->m1] > ZERO_EPSILON ||
			 this->rs->media_refractRoughnessV.v[this->m1] > ZERO_EPSILON ||
			 this->rs->media_refractRoughnessU.v[this->m2] > ZERO_EPSILON ||
			 this->rs->media_refractRoughnessV.v[this->m2] > ZERO_EPSILON )
		{
			do_blurryRefraction = true;
		}
		return do_blurryRefraction;
	}

	bool setupDispersion(JFND_Shader_Data * data ) 
	{
		bool do_disperse = false;

		if (!this->rs->ray_monochromatic)
		{
			// already monochromatic means don't disperse, just keep the existing wavelength and trace that way
			if ( this->rs->media_disperse.v[this->currentID] )
				this->rs->spectral_LUT_ptr = &data->spectral_LUT_ptr;

			do_disperse = this->rs->media_disperse.v[this->m2]; // only disperse if the medium we're entering is dispersey, and the ray is not already monochromatic
		}

		return do_disperse;
	}

	AtColor getTransmissionColor(AtShaderGlobals * sg)
	{
		float depth = (float) sg->Rl;
		if (this->t1 != AI_RGB_WHITE)
		{		
			return AiColorCreate( 
				pow(this->t1.r, depth), 
				pow(this->t1.g, depth), 
				pow(this->t1.b, depth));
		}
		else
		{
			return AI_RGB_WHITE;
		}
	}

	AtColor getShadowTransparency(AtShaderGlobals * sg, int shadowMode) 
	{
		float fresnelTerm = 0.0f;
		switch ( shadowMode )
		{
			case sh_black:
				return AI_RGB_BLACK;
			case sh_transmit_only:		
				return this->getTransmissionColor(sg);
			case sh_transmit_and_outer_fresnels:
				if (this->validInterface && this->entering)
				{					
					fresnelTerm = fresnelEquations (this->n1, this->n2, AiV3Dot(sg->Nf, -sg->Rd), this->polarizationTerm, false);
				}
				return this->getTransmissionColor(sg) * (AI_RGB_WHITE - fresnelTerm);
			case sh_transmit_and_all_fresnels:
				if (this->validInterface)
				{
					fresnelTerm = fresnelEquations (this->n1, this->n2, AiV3Dot(sg->Nf, -sg->Rd), this->polarizationTerm, false);
				}
				return this->getTransmissionColor(sg) * (AI_RGB_WHITE - fresnelTerm);
		}
		return AI_RGB_WHITE;
	}

	float getSpecRoughnessU() 
	{
		return (this->rs->media_specRoughnessU.v[this->m2] + this->rs->media_specRoughnessU.v[this->m1]) / 2.0f;
	}

	float getSpecRoughnessV()
	{
		return (this->rs->media_specRoughnessV.v[this->m2] + this->rs->media_specRoughnessV.v[this->m1]) / 2.0f;		
	}

	float getRefrRoughnessU()
	{
		return refractiveRoughness( this->rs->media_refractRoughnessU.v[this->m1], this->rs->media_refractRoughnessU.v[this->m2], this->n1, this->n2 );
	}
	
	float getRefrRoughnessV() 
	{
		return refractiveRoughness( this->rs->media_refractRoughnessV.v[this->m1], this->rs->media_refractRoughnessV.v[this->m2], this->n1, this->n2 );
	}

	int getBTDFType() 
	{
		return this->rs->media_BTDF.v[this->m_higherPriority];
	}

	bool refractionNeedsUserTangent() 
	{
		return this->getBTDFType() == b_ward_userTangent;
	}

	void* getBTDFData(AtShaderGlobals *ppsg, float refr_roughnessU, float refr_roughnessV, AtVector customTangentVector) 
	{		
		AtVector tangentSourceVector, uTangent, vTangent;
		switch ( this->getBTDFType() )
		{
			case b_stretched_phong:
				// Stretched Phong
				return AiStretchedPhongMISCreateData(ppsg, (0.5f / SQR(refr_roughnessU) - 0.5f));
				break;
			case b_cook_torrance:
				// Cook Torrance
				return AiCookTorranceMISCreateData(ppsg, &AI_V3_ZERO, &AI_V3_ZERO, refr_roughnessU, refr_roughnessU);
				break;
			case b_ward_rayTangent:
				// Ward with refraction-derivitive tangents
				tangentSourceVector = AiV3Normalize(ppsg->Rd);
				blurAnisotropicPoles(&refr_roughnessU, &refr_roughnessV, &this->rs->media_blurAnisotropicPoles.v[this->m_higherPriority], &ppsg->N, &tangentSourceVector);
				uTangent = AiV3Cross(ppsg->Nf, AiV3Normalize(ppsg->Rd)); 
				vTangent = AiV3Cross(ppsg->Nf, uTangent);
				return AiWardDuerMISCreateData(ppsg, &uTangent, &vTangent, refr_roughnessU, refr_roughnessV); 
				break;
			case b_ward_userTangent:
				// Ward with user tangents
				tangentSourceVector = AiV3Normalize( customTangentVector );
				blurAnisotropicPoles(&refr_roughnessU, &refr_roughnessV, &this->rs->media_blurAnisotropicPoles.v[this->m_higherPriority], &ppsg->N, &tangentSourceVector);
				uTangent = AiV3Cross(ppsg->Nf, tangentSourceVector); 
				vTangent = AiV3Cross(ppsg->Nf, uTangent);
				return AiWardDuerMISCreateData( ppsg, &uTangent, &vTangent, refr_roughnessU, refr_roughnessV ) ; 
				break;
		}
		return NULL;
	}

	AtVector getBTDFSample( void* btdf_data, SAMPLETYPE refraction_sample[2]) {
		switch ( this->getBTDFType() )
		{
			case b_stretched_phong:
				return AiStretchedPhongMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);
			case b_cook_torrance:
				return AiCookTorranceMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);
			case b_ward_rayTangent:
				return AiWardDuerMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);
			case b_ward_userTangent:
				return AiWardDuerMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);
		}
		return AI_V3_ZERO;
	}

	bool directRefractionNeedsUserTangent(int dr_btdf) 
	{
		return dr_btdf == b_ward_userTangent;
	}


	void getDirectRefractionBTDFs(int dr_btdf, AtShaderGlobals *ppsg, float dr_roughnessU, float dr_roughnessV, 
		float drs_roughnessU, float drs_roughnessV, AtVector customTangentVector, void **btdfA, void **btdfB) 
	{
		AtVector tangentSourceVector, uTangent, vTangent;
		switch ( dr_btdf )
		{
			case b_stretched_phong:
				// Stretched Phong
				*btdfA = AiStretchedPhongMISCreateData(ppsg, (0.5f / SQR(dr_roughnessU) - 0.5f));
				*btdfB = AiStretchedPhongMISCreateData(ppsg, (0.5f / SQR(drs_roughnessU) - 0.5f));
				break;
			case b_cook_torrance:
				// Cook Torrance
				*btdfA = AiCookTorranceMISCreateData(ppsg, &AI_V3_ZERO, &AI_V3_ZERO, dr_roughnessU, dr_roughnessU);
				*btdfB = AiCookTorranceMISCreateData(ppsg, &AI_V3_ZERO, &AI_V3_ZERO, drs_roughnessU, drs_roughnessU);
				break;
			case b_ward_rayTangent:
				// Ward with refraction-derivitive tangents
				tangentSourceVector = AiV3Normalize(ppsg->Rd);
				blurAnisotropicPoles(&dr_roughnessU, &dr_roughnessV, &this->rs->media_blurAnisotropicPoles.v[this->m_higherPriority], &ppsg->Nf, &tangentSourceVector);
				blurAnisotropicPoles(&drs_roughnessU, &drs_roughnessV, &this->rs->media_blurAnisotropicPoles.v[this->m_higherPriority], &ppsg->Nf, &tangentSourceVector);
				uTangent = AiV3Cross(ppsg->Nf, tangentSourceVector ); 
				vTangent = AiV3Cross(ppsg->Nf, uTangent);
				*btdfA = AiWardDuerMISCreateData(ppsg, &uTangent, &vTangent, dr_roughnessU, dr_roughnessV); 
				*btdfB = AiWardDuerMISCreateData(ppsg, &uTangent, &vTangent, drs_roughnessU, drs_roughnessV); 
				break;
			case b_ward_userTangent:
				// Ward with user tangents
				tangentSourceVector = AiV3Normalize( customTangentVector );
				blurAnisotropicPoles(&dr_roughnessU, &dr_roughnessV, &this->rs->media_blurAnisotropicPoles.v[this->m_higherPriority], &ppsg->Nf, &tangentSourceVector);
				blurAnisotropicPoles(&drs_roughnessU, &drs_roughnessV, &this->rs->media_blurAnisotropicPoles.v[this->m_higherPriority], &ppsg->Nf, &tangentSourceVector);
				uTangent = AiV3Cross(ppsg->Nf, tangentSourceVector ); 
				vTangent = AiV3Cross(ppsg->Nf, uTangent);					
				*btdfA = AiWardDuerMISCreateData( ppsg, &uTangent, &vTangent, dr_roughnessU, dr_roughnessV ) ; 
				*btdfB = AiWardDuerMISCreateData( ppsg, &uTangent, &vTangent, drs_roughnessU, drs_roughnessV ) ; 
				break;
		}
	}

	void directRefractionSampleLights(AtShaderGlobals * sg, int dr_btdf, bool refract_skydomes, bool two_lobes, 
		void* btdf_data_direct, void* btdf_data_direct2, AtColor * acc_refract_direct, 
		AtColor * acc_refract_direct_second) 
	{
		AiLightsPrepare(sg);
		AiStateSetMsgBool("opaqueShadowMode", true);
		while (AiLightsGetSample(sg))
		{
			if ( refract_skydomes || !AiNodeIs( sg->Lp,"skydome_light" ))
			{
				switch ( dr_btdf )
				{
					case b_stretched_phong:
						*acc_refract_direct += AiEvaluateLightSample(sg, btdf_data_direct, AiStretchedPhongMISSample, AiStretchedPhongMISBRDF, AiStretchedPhongMISPDF);
						if (two_lobes)
							*acc_refract_direct_second += AiEvaluateLightSample(sg, btdf_data_direct2, AiStretchedPhongMISSample, AiStretchedPhongMISBRDF, AiStretchedPhongMISPDF);
						break;
					case b_cook_torrance:
						*acc_refract_direct += AiEvaluateLightSample(sg, btdf_data_direct, AiCookTorranceMISSample, AiCookTorranceMISBRDF, AiCookTorranceMISPDF);
						if (two_lobes)
							*acc_refract_direct_second += AiEvaluateLightSample(sg, btdf_data_direct2, AiCookTorranceMISSample, AiCookTorranceMISBRDF, AiCookTorranceMISPDF);
						break;
					case b_ward_rayTangent:
						*acc_refract_direct += AiEvaluateLightSample(sg, btdf_data_direct, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
						if (two_lobes)
							*acc_refract_direct_second += AiEvaluateLightSample(sg, btdf_data_direct2, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
						break;
					case b_ward_userTangent:
						*acc_refract_direct += AiEvaluateLightSample(sg, btdf_data_direct, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
						if (two_lobes)
							*acc_refract_direct_second += AiEvaluateLightSample(sg, btdf_data_direct2, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
						break;
				}
			}
		}
		AiStateSetMsgBool("opaqueShadowMode", false);
	}





} InterfaceInfo;



typedef struct TraceSwitch
{
	bool refr_ind;
	bool refr_dir;
	bool spec_ind;
	bool spec_dir;
	bool TIR;

	TraceSwitch(InterfaceInfo * iinfo) 
	{
		Ray_State_Datatype * rs = iinfo->rs;
		this->refr_ind = rs->media_refractIndirect.v[iinfo->m1] > ZERO_EPSILON 
			&& rs->media_refractIndirect.v[iinfo->m2] > ZERO_EPSILON;
		this->refr_dir = rs->media_refractDirect.v[iinfo->currentID] > ZERO_EPSILON 
			&& iinfo->mediaExit; // media exit is AND'd in later. 
		this->spec_ind = rs->media_specIndirect.v[iinfo->m_higherPriority] > ZERO_EPSILON;
		this->spec_dir = iinfo->mediaEntrance && rs->media_specDirect.v[iinfo->m_higherPriority] > ZERO_EPSILON;
		this->TIR = this->refr_ind || this->spec_ind;
	}


	void setTraceNone() 
	{
		this->refr_ind = false;
		this->refr_dir = false;
		this->spec_ind = false;
		this->spec_dir = false;
		this->TIR = false;
	}

	void setInternalReflections(InterfaceInfo * iinfo, bool internal_reflections_enabled) 
	{
		bool internal = !iinfo->mediaEntrance;
		if (internal && !internal_reflections_enabled)
			this->spec_ind = false;
	}

	void setPathtracedCaustic(InterfaceInfo * iinfo) 
	{
		Ray_State_Datatype * rs = iinfo->rs;
		this->refr_dir = this->refr_dir && rs->caustic_refractDirect;
		this->refr_ind = this->refr_ind && true;
		this->spec_ind = this->spec_ind && rs->caustic_specIndirect;;
		this->spec_dir = this->spec_dir && rs->caustic_specDirect;
		this->TIR = this->TIR && rs->caustic_TIR && (rs->caustic_refractDirect || rs->caustic_refractIndirect);

		if (iinfo->mediaExit)
		{
			// Indirect refraction is only needed on mediaExit if we actually want caustics from indirect refractions.
			// Otherwise they are carriers for direct refraction. 
			this->refr_ind = this->refr_ind && rs->caustic_refractIndirect;
		}
		else
		{
			// Indrect refraction is needed if we want caustics from indirect refraction, or direct refraction, 
			this->refr_ind = this->refr_ind && (rs->caustic_refractDirect || rs->caustic_refractIndirect);
		}

		if (iinfo->mediaEntrance)
		{
			// if we have an inside out surface, TIR caustics will only be traced if indirect specular is on. 
			this->TIR = this->spec_ind;
		}
		else
		{
			this->spec_ind = this->spec_ind && rs->caustic_specInternal;
		}

	}

	void setPhotonCaustic(InterfaceInfo * iinfo) 
	{
		Ray_State_Datatype * rs = iinfo->rs;
		this->refr_dir = false;
		// in photon land, direct refractions means the light refracts straight through.
		// indirect refractions means.. very little. 
		// to the renderer this would always be indirect refraction though. 
		this->refr_ind = this->refr_ind && rs->caustic_refractDirect;

		// in photon land, direct specular caustics will mean reflections off the outside of the surface.
		// indirect specular will mean the caustics will beounce off interior surfaces too. 
		// to the renderer this would always be indirect specular though. 
		this->spec_ind = this->spec_ind && (
			(rs->caustic_specDirect && iinfo->entering)
			|| (rs->caustic_specIndirect)
		);

		this->spec_dir = false;
		this->TIR = this->TIR && rs->caustic_TIR;
	}

	bool traceAnyRefraction() 
	{
		return this->refr_ind || this->refr_dir;
	}

	bool traceAnySpecular()
	{
		return this->spec_ind || this->spec_dir;
	}

	bool traceAnything()
	{
		return this->traceAnyRefraction() || this->traceAnySpecular();
	}
};
