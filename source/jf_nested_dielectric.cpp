/* 
 * JF Nested Dielectric by Jonah Friedman
 * 1.0.5
 * Copyright (c) 2014, Psyop Media Company, LLC and Jonah Friedman 
 * Open sourced under the 3 clause BSD license, see license.txt
 */

#include <ai.h>
#include <cstring>
#include <string>
#include <math.h>
#include "spectral_distributions.h"
#include "jf_nested_dielectric.h"

AI_SHADER_NODE_EXPORT_METHODS(jf_nested_dielectric_methods);

enum jf_nested_dielectricParams
{
	p_mediumPriority,
	p_mediumIOR,
	p_mediumTransmittance,
	p_mediumTransmittance_scale,

	p_polarize,
	p_polarization_angle,

	p_refraction_scale,
	p_btdf,
	p_refraction_roughness_u,
	p_refraction_roughness_v,

	p_indirect_refraction,

	p_direct_refraction,
	p_dr_roughnessOffset,
	p_dr_roughnessDepthMultiplier,
	p_dr_roughnessDepthAdder,
	p_dr_use_refraction_btdf,
	p_dr_btdf,
	p_dr_roughness_u,
	p_dr_roughness_v,
	p_dr_second_scale,
	p_dr_second_roughnessMultiplier,

	p_disperse,
	p_dispersion,
	p_spectral_distribution,
	p_spectral_gamut,
	p_spectral_clamp_negative_colors,
	p_spectral_saturation,

	p_specular_scale,
	p_brdf,
	p_direct_specular,
	p_indirect_specular,
	p_specular_roughness_u,
	p_specular_roughness_v,

	p_specular_ray_type,

	p_blur_anisotropic_poles,

	p_enable_internal_reflections,
	p_russian_roulette_probability,
	p_energy_cutoff_exponent,

	p_ward_tangent,

	p_refract_skydomes,
	p_refract_skies,

	p_reflect_skydomes,
	p_reflect_skies,

	p_emission_at_interfaces,
	p_emission_at_exits,
	p_emission_at_entrances,

	p_shadow_mode,
	p_caustic_mode,

	p_caustic_refractions_direct,
	p_caustic_dr_roughness_offset,
	p_caustic_refractions_indirect,
	p_caustic_TIR,
	p_caustic_internal_speculars,
	p_caustic_dispersion,

	p_caustic_specular_direct,
	p_caustic_specular_indirect,

	p_caustic_max_distance,

	p_aov_direct_refraction,
	p_aov_indirect_refraction,
	p_aov_direct_specular,
	p_aov_indirect_specular
};

node_parameters
{
	AiParameterINT("mediumPriority", 0);
	AiParameterFLT("mediumIOR", 1.33f);
	AiParameterRGB("mediumTransmittance", 1.0f, 1.0f, 1.0f);
	AiParameterFLT("mediumTransmittance_scale", 1.0f);

	AiParameterBOOL("polarize", false);
	AiParameterFLT("polarization_angle", 0.0f);

	AiParameterFLT("refraction_scale", 1.0f);
	AiParameterENUM("btdf", b_cook_torrance, enum_brdfs);
	AiParameterFLT("refraction_roughness_u", 0.00f);
	AiParameterFLT("refraction_roughness_v", 0.00f);

	AiParameterFLT("indirect_refraction", 1.0f);

	AiParameterFLT("direct_refraction", 1.0f);
	AiParameterFLT("dr_roughnessOffset", 0.05f);
	AiParameterFLT("dr_roughnessDepthMultiplier", 1.0f);
	AiParameterFLT("dr_roughnessDepthAdder", 0.0f);
	AiParameterBOOL("dr_use_refraction_btdf", true);
	AiParameterENUM("dr_btdf", b_cook_torrance, enum_brdfs);
	AiParameterFLT("dr_roughness_u", 0.00f);
	AiParameterFLT("dr_roughness_v", 0.00f);
	AiParameterFLT("dr_second_scale", 0.0f);
	AiParameterFLT("dr_second_roughnessMultiplier", 2.0f);

	AiParameterBOOL("disperse", false);
	AiParameterFLT("dispersion", 0.03f);
	AiParameterENUM("spectral_distribution", s_daylight, enum_spectra);
	AiParameterENUM("spectral_gamut", g_sRGB, enum_gamuts);
	AiParameterBOOL("spectral_clamp_negative_colors", true);
	AiParameterFLT("spectral_saturation", 1.00f);

	AiParameterFLT("specular_scale", 1.0f);
	AiParameterENUM("brdf", b_cook_torrance, enum_brdfs);
	AiParameterFLT("direct_specular", 1.0f);
	AiParameterFLT("indirect_specular", 1.0f);
	AiParameterFLT("specular_roughness_u", 0.00f);
	AiParameterFLT("specular_roughness_v", 0.00f);

	AiParameterENUM("specular_ray_type", rt_glossy, enum_specular_raytypes);
	
	AiParameterFLT("blur_anisotropic_poles", 0.0f);

	AiParameterBOOL("enable_internal_reflections", true);
	AiParameterFLT("russian_roulette_probability", 0.5f);
	AiParameterINT("energy_cutoff_exponent", -5);

	AiParameterVEC("ward_tangent", 0.0f, 1.0f, 0.0f );

	AiParameterBOOL("refract_skydomes", false);
	AiParameterBOOL("refract_skies", true);	

	AiParameterBOOL("reflect_skydomes", true);
	AiParameterBOOL("reflect_skies", true);

	AiParameterRGB("emission_at_interfaces", 0.0f, 0.0f, 0.0f );
	AiParameterRGB("emission_at_exits", 0.0f, 0.0f, 0.0f );
	AiParameterRGB("emission_at_entrances", 0.0f, 0.0f, 0.0f );

	AiParameterENUM("shadow_mode", sh_transmit_and_outer_fresnels, enum_shadow_modes);

	AiParameterINT("caustic_mode", 0);
	AiParameterBOOL("caustic_refractions_direct", true);
	AiParameterFLT("caustic_dr_roughness_offset", 0.1f);
	AiParameterBOOL("caustic_refractions_indirect", false);
	AiParameterBOOL("caustic_TIR", false);
	AiParameterBOOL("caustic_internal_speculars", false);
	AiParameterBOOL("caustic_dispersion", false);

	AiParameterBOOL("caustic_specular_direct", false);
	AiParameterBOOL("caustic_specular_indirect", false);

	AiParameterFLT("caustic_max_distance", -1.0f);

	AiParameterSTR("aov_direct_refraction", "");
	AiParameterSTR("aov_indirect_refraction", "Refraction");
	AiParameterSTR("aov_direct_specular", "Direct_Specular");
	AiParameterSTR("aov_indirect_specular", "Indirect_Specular");
}



// ---------------------------------------------------//
// - Arnold Shader Code 
// ---------------------------------------------------//

node_initialize
{
	ShaderData* data = new ShaderData;
	AiNodeSetLocalData(node,data);
	data->dispersion_sampler = NULL;
	data->specular_sampler = NULL;
	data->refraction_sampler = NULL;
	data->russian_roullete_single_sampler = NULL;
};

node_finish
{
	if (AiNodeGetLocalData(node) != NULL)
	{
		ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

		AiSamplerDestroy(data->dispersion_sampler);
		AiSamplerDestroy(data->specular_sampler);
		AiSamplerDestroy(data->refraction_sampler);
		AiSamplerDestroy(data->russian_roullete_single_sampler);

		AiNodeSetLocalData(node, NULL);
		delete data;
	}
}


node_update
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	AiSamplerDestroy(data->dispersion_sampler);
	AiSamplerDestroy(data->specular_sampler);
	AiSamplerDestroy(data->refraction_sampler);
	AiSamplerDestroy(data->russian_roullete_single_sampler);

	data->dispersion_sampler = AiSamplerSeeded( base_sampler_seed + 1, AiNodeGetInt(AiUniverseGetOptions(), "GI_refraction_samples"), 2 );	
	data->specular_sampler = AiSamplerSeeded( base_sampler_seed + 2, AiNodeGetInt(AiUniverseGetOptions(), "GI_glossy_samples"), 2 );	
	data->refraction_sampler = AiSamplerSeeded( base_sampler_seed + 3, AiNodeGetInt(AiUniverseGetOptions(), "GI_refraction_samples"), 2 );
	data->russian_roullete_single_sampler = AiSamplerSeeded( base_sampler_seed + 7, 1, 2 );

	data->aov_direct_refraction = AiNodeGetStr(node, "aov_direct_refraction");
	data->aov_indirect_refraction = AiNodeGetStr(node, "aov_indirect_refraction");
	data->aov_direct_specular = AiNodeGetStr(node, "aov_direct_specular");
	data->aov_indirect_specular = AiNodeGetStr(node, "aov_indirect_specular");

	const bool do_disperse = AiNodeGetBool(node, "disperse");
	if (do_disperse)
	{
		const int spectrum_selection = AiNodeGetInt(node, "spectral_distribution");
		const int gamut_selection = AiNodeGetInt(node, "spectral_gamut");
		const float saturation = AiNodeGetFlt(node, "spectral_saturation");
		const bool clamp = AiNodeGetBool(node, "spectral_clamp_negative_colors");
		data->spectral_LUT_ = build_nonuniform_spectral_LUT(spectrum_selection, gamut_selection, saturation, clamp, AI_RGB_WHITE);
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

		data->polarizationVector = AiV3Normalize(pfilterInWorldSpace);
	}
	else
	{
		data->polarizationVector = AI_V3_ZERO;
	}
}




shader_evaluate
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	Ray_State_Datatype *RayState;
	Ray_State_Cache_Datatype RayStateCache;

	// ---------------------------------------------------//
	// - Array and messaging preprocess 
	// ---------------------------------------------------//

	bool msgs_are_valid = false;
	AiStateGetMsgBool("msgs_are_valid", &msgs_are_valid);

	if ( msgs_are_valid )
	{
		void * rayState_ptr;
		AiStateGetMsgPtr("rayState_ptr", &rayState_ptr);
		RayState = static_cast<Ray_State_Datatype*>( rayState_ptr );

		cacheRayState( RayState, &RayStateCache );
	}
	else
	{		
		RayState = static_cast<Ray_State_Datatype*>( AiShaderGlobalsQuickAlloc(sg, sizeof( Ray_State_Datatype ) ) );
		AiStateSetMsgPtr("rayState_ptr", RayState);

		RayState->ray_monochromatic = false;
		RayState->ray_wavelength = 0.0f;
		RayState->ray_TIRDepth = 0;
		RayState->ray_invalidDepth = 0;
		RayState->ray_energy = AI_RGB_WHITE;

		// Array initialization - possibly unnecessary except mediaInside? TO DO: investigate
		for( int i = 0; i < max_media_count; i++ )
		{
			RayState->media_inside.v[i] = 0;
			RayState->media_iOR.v[i] = 1.0f;
			RayState->media_disperse.v[i] = false;
			RayState->media_dispersion.v[i] = 0.0f;
			RayState->media_BTDF.v[i] = 0;
			RayState->media_BRDF.v[i] = 0;

			RayState->media_refractDirect.v[i] = 1.0f;
			RayState->media_refractIndirect.v[i] = 1.0f;
			RayState->media_specDirect.v[i] = 1.0f;
			RayState->media_specIndirect.v[i] = 1.0f;

			RayState->media_refractRoughnessU.v[i] = 0.0f;
			RayState->media_refractRoughnessV.v[i] = 0.0f;
			RayState->media_specRoughnessU.v[i] = 0.0f;
			RayState->media_specRoughnessV.v[i] = 0.0f;
			RayState->media_transmission.v[i] = AI_RGB_WHITE;

			RayState->media_blurAnisotropicPoles.v[i] = 0.0f;
		}

		// caustics behavior values get set once a ray enters a medium from a diffuse ray
		RayState->caustic_behaviorSet = false;
		RayState->caustic_mode = 0;
		RayState->caustic_refractDirect =
			RayState->caustic_refractIndirect =
			RayState->caustic_TIR =
			RayState->caustic_specInternal =
			RayState->caustic_dispersion =
			RayState->caustic_specDirect =
			RayState->caustic_specIndirect =
			false;

		// Set once values - values set at initialization for all descendent rays
		RayState->energy_cutoff = pow(10.0f, (float) AiShaderEvalParamInt(p_energy_cutoff_exponent) );		
		RayState->polarized = AiShaderEvalParamBool(p_polarize);
		RayState->polarizationVector = data->polarizationVector;
	}
	AiStateSetMsgBool("msgs_are_valid", true); // Any child rays from this will find valid messages. 


	mediaIntStruct * media_inside_ptr;
	
	if (sg->Rt == AI_RAY_SHADOW)
	{
		/*
		 * Shadow rays need their own mediaInside tracking.
		 * The shadow_media_inside struct is initialized by copying the media_inside struct
		 * Under the following conditions the struct is re-initialized:
		 *   The previous sample evaluated was not a shadow ray
		 *   The transparency index decreased or stayed the same, meaning the previous was a different shadow ray
		 *   The shadow list has never been initialzied
		 */

		media_inside_ptr = &RayState->shadow_media_inside;

		bool shadowlist_is_valid = false;
		AiStateGetMsgBool("shadowlist_is_valid", &shadowlist_is_valid);


		bool transp_index_reset = false;
		int prev_transp_index;
		if (AiStateGetMsgInt( "prev_transp_index", &prev_transp_index))
		{
			if ((int) sg->transp_index <= prev_transp_index)
			{
				transp_index_reset = true;	
			}
		}

		if ( !shadowlist_is_valid || transp_index_reset )  // if there has been another kind of ray, or the transp_index did not increase
		{
			// intialize the shadow media inside list
			memcpy(&RayState->shadow_media_inside, &RayState->media_inside, sizeof(mediaIntStruct) );
		}

		AiStateSetMsgInt("prev_transp_index", sg->transp_index);
		AiStateSetMsgBool("shadowlist_is_valid", true);
	}
	else
	{
		media_inside_ptr = &RayState->media_inside;

		AiStateSetMsgBool("shadowlist_is_valid", false);
	}


	/*
	 * --------------------------------------------------- *
	 * - shader parameters setup 
	 * --------------------------------------------------- *
	 *
	 * m_cMatID is the medium ID of the material we're currently evaluating, regardless of
	 * whether or not the interface is valid or what media we're inside.
	 *
	 * The medium IDs all get 1 added to them. This means 0 is reserved for "no medium". This simplifies other code, but it means
	 * that if the user gives a medium an ID of 0, internally in the shader it has an ID of 1. 
	 */

	const int m_cMatID = AiShaderEvalParamInt(p_mediumPriority) + 1;

	if (media_inside_ptr->v[m_cMatID] < -10 )
	{
		/* 
		 * Sometimes two pieces of geometry overlapping cause this to go crazy, with values like -60
		 * I think -10 exceeds any plausible correct value of this
		 * Really -1 shows that things are modeled improperly.
		 * In any case, this fixes it and the warning is helpful.
		 */

		char * const overlapping_surfaces_message = "JF Nested Dielectric: Crazy values in media lists, you may have some perfectly overlapping surfaces.";
		AiMsgWarning(overlapping_surfaces_message);
		sg->out.RGBA = AI_RGBA_BLACK;
		return; 
	}
	if ( media_inside_ptr->v[m_cMatID] > 30)
	{
		/* 
		 * Sometimes two pieces of geometry overlapping cause this to go crazy, with values like -60
		 * I think -10 exceeds any plausible correct value of this
		 * Really -1 shows that things are modeled improperly.
		 * In any case, this fixes it and the warning is helpful.
		 */

		char * const overlapping_surfaces_message = "JF Nested Dielectric: Crazy positive values in media lists, you may have some perfectly overlapping surfaces.";
		AiMsgWarning(overlapping_surfaces_message);
		sg->out.RGBA = AI_RGBA_BLACK;
		return; 
	}

	if (m_cMatID >= max_media_count || m_cMatID < 0)
	{
		char * const priority_error_message = "JF Nested Dielectric: Medium priority must be between 0 and 32!";
		AiMsgError(priority_error_message);
		return;
	}


	{			
		const AtColor _cMatMediumTransmission_color = AiShaderEvalParamRGB(p_mediumTransmittance);
		const float _cMatMediumTransmission_scale = 1.0f / AiShaderEvalParamFlt(p_mediumTransmittance_scale);
		const AtColor cMatMediumTransmission = 
		{
			pow(_cMatMediumTransmission_color.r, _cMatMediumTransmission_scale),
			pow(_cMatMediumTransmission_color.g, _cMatMediumTransmission_scale),
			pow(_cMatMediumTransmission_color.b, _cMatMediumTransmission_scale)
		};
		RayState->media_iOR.v[m_cMatID] = AiShaderEvalParamFlt(p_mediumIOR);
		RayState->media_disperse.v[m_cMatID] =  AiShaderEvalParamBool(p_disperse);
		RayState->media_dispersion.v[m_cMatID] = AiShaderEvalParamFlt(p_dispersion);
		RayState->media_transmission.v[m_cMatID] = cMatMediumTransmission;
		RayState->media_BTDF.v[m_cMatID] = AiShaderEvalParamEnum(p_btdf);
		RayState->media_BRDF.v[m_cMatID] = AiShaderEvalParamEnum(p_brdf);
		RayState->media_refractRoughnessU.v[m_cMatID] = refractRoughnessConvert( AiShaderEvalParamFlt(p_refraction_roughness_u) );
		RayState->media_specRoughnessU.v[m_cMatID] = AiShaderEvalParamFlt(p_specular_roughness_u);

		const float _refractionScale = AiShaderEvalParamFlt(p_refraction_scale);
		RayState->media_refractDirect.v[m_cMatID] = AiShaderEvalParamFlt(p_direct_refraction) * _refractionScale;
		RayState->media_refractIndirect.v[m_cMatID] = AiShaderEvalParamFlt(p_indirect_refraction) * _refractionScale;	

		const float _specularScale = AiShaderEvalParamFlt(p_specular_scale);
		RayState->media_specDirect.v[m_cMatID] = AiShaderEvalParamFlt(p_direct_specular) * _specularScale;
		RayState->media_specIndirect.v[m_cMatID] = AiShaderEvalParamFlt(p_indirect_specular) * _specularScale;

		if (RayState->media_BTDF.v[m_cMatID] >= 2)
			RayState->media_refractRoughnessV.v[m_cMatID] = refractRoughnessConvert( AiShaderEvalParamFlt(p_refraction_roughness_v) );
		else
			RayState->media_refractRoughnessV.v[m_cMatID] = 0.0f;

		if (RayState->media_BRDF.v[m_cMatID] >= 2)
			RayState->media_specRoughnessV.v[m_cMatID] = AiShaderEvalParamFlt(p_specular_roughness_v);
		else
			RayState->media_specRoughnessV.v[m_cMatID] = 0;

		RayState->media_blurAnisotropicPoles.v[m_cMatID] = AiShaderEvalParamFlt(p_blur_anisotropic_poles);
	}

	// ---------------------------------------------------//
	// - polarization logic     
	// 		(determine the polarization term for use in the fresnel equations. )
	// ---------------------------------------------------//

	float cPolarizationTerm = 0.5f;

	if (RayState->polarized)
	{
		cPolarizationTerm = std::abs( AiV3Dot(sg->Nf, RayState->polarizationVector ));
	}

	// ---------------------------------------------------//
	// - media logic     
	// 		(determine what media we're in based on mediaInside arrays. )
	// ---------------------------------------------------//

	/*
	 * startingMedium is the medium we were already inside when we started tracing
	 * startingMediumSecondary is the next medium we will be inside, if it turns out we're leaving this one
	 *
	 * example: we're inside glass and water, glass is the priority. startingMedium is glass, startingMediumSecondary is water.
	 * This is all evaulated without any knowledge of what surface we're evaluating, this is just parsing the mediaInside arrays. 
	 */

	int startingMedium = 0, startingMediumSecondary = 0; 

	for( int i = 0; i < max_media_count; i++ )
	{
	
	if ( media_inside_ptr->v[i] > 0 ) // if we find a medium that we're inside
		{
			if ( startingMedium == 0 ) // ..and we havent found our current medium yet
			{
				startingMedium = i;  // then we've found our medium.
				if (media_inside_ptr->v[i] > 1) // ...and if we're in more than one of these
				{
					startingMediumSecondary = i; // then our second medium is the same medium.
					break;
				}
			}
			else if ( startingMediumSecondary == 0 )  // ..and we've already found our current medium
			{
				startingMediumSecondary = i; // ..then we've found our second medium
				break;
			}
		}
	}

	// ---------------------------------------------------//
	// - interface logic 
	//		(Figure out the properties and validity of the interface we're tracing. )
	// ---------------------------------------------------//
		
	const bool entering = ( AiV3Dot(sg->Ng, sg->Rd) < 0.0f ) ;

	/*
	 * m = medium ID, n = IOR, T = transmission
	 * 1 is the previous medium, 2 is the next, as seen in n1 and n2 of refraction diagrams
	 * such diagrams: http://en.wikipedia.org/wiki/File:Snells_law.svg
	 */

	int m1 = 0, m2 = 0;
	float n1 = 1.0f, n2 = 1.0f;
	AtColor t1;
	AtColor t2;

	bool validInterface = false, mediaExit = false, mediaEntrance = false;
	
	if (entering) // (entering a possible medium)
	{
		if (startingMedium == 0)
		{
			// entering the first medium, so no media information is available
			validInterface = mediaEntrance = true;
			m1 = 0;
			m2 = m_cMatID;
		}
		else if ( m_cMatID < startingMedium )
		{
			// entering a new highest priority medium
			validInterface = true;
			m1 = startingMedium;
			m2 = m_cMatID;
		}
		else if ( m_cMatID >= startingMedium )
		{
			// entering a medium that is not higher priority, interface is invalid
			// or entering the same medium we're already in, interface is also invalid
			m1 = m2 = startingMedium;
		}
	}
	else // (leaving a possible medium)
	{
		if (startingMedium == m_cMatID)
		{
			// leaving the highest priority medium, which may be 0 (meaning no medium)
			validInterface = true;
			m1 = startingMedium;
			m2 = startingMediumSecondary;
			if (startingMediumSecondary == 0)
			{
				// and there is no second medium, meaning we're leaving all media
				mediaExit = true;
			}
		}
		else
		{
			// leaving a medium that is not the current medium, interface invalid
			m1 = m2 = startingMedium;
		}
	}

	int m_higherPriority;
	if (m1 == 0) 
		m_higherPriority = m2;
	else if (m2 == 0) 
		m_higherPriority = m1;
	else if (m1 < m2) 
		m_higherPriority = m1;
	else if (m1 > m2) 
		m_higherPriority = m2;

	// - Interface Logic
	// 		(Disperse and multisample)

	bool do_blurryRefraction = false;
	if ( RayState->media_refractRoughnessU.v[m1] > ZERO_EPSILON ||
		 RayState->media_refractRoughnessV.v[m1] > ZERO_EPSILON ||
		 RayState->media_refractRoughnessU.v[m2] > ZERO_EPSILON ||
		 RayState->media_refractRoughnessV.v[m2] > ZERO_EPSILON )
	{
		do_blurryRefraction = true;
	}
	bool do_disperse = false;
	if (!RayState->ray_monochromatic)
	{
		// already monochromatic means don't disperse, just keep the existing wavelength and trace that way
		if ( RayState->media_disperse.v[m_cMatID] )
			RayState->spectral_LUT_ = &data->spectral_LUT_; // copy over the pointer to the spectral LUT

		if (RayState->media_disperse.v[m2])
			do_disperse = RayState->media_disperse.v[m2];
	}
	const bool do_multiSampleRefraction = do_disperse || do_blurryRefraction;


	if ( validInterface )
	{	
		if (RayState->ray_monochromatic)
		{
			n1 = dispersedIOR( RayState->media_iOR.v[m1], RayState->media_dispersion.v[m1], RayState->ray_wavelength );
			n2 = dispersedIOR( RayState->media_iOR.v[m2], RayState->media_dispersion.v[m2], RayState->ray_wavelength );
		}
		else
		{
			n1 = RayState->media_iOR.v[m1];
			n2 = RayState->media_iOR.v[m2];
		}


		if (n1 == n2)
		{			
			validInterface = false;
		}
	}


	t1 = RayState->media_transmission.v[m1];
	t2 = RayState->media_transmission.v[m2];

	AtColor cTransmission = transmissionColor(&t1, (float) sg->Rl) ;

	// ---------------------------------------------------//
	// - Shadow rays
	// ---------------------------------------------------//
	
	if (sg->Rt == AI_RAY_SHADOW)
	{
		AtColor transparency;
		float fresnelTerm = 0.0f;

		int shadowMode;
		bool opaqueShadowMode = false;
		AiStateGetMsgBool("opaqueShadowMode", &opaqueShadowMode);
		if (opaqueShadowMode) 
			shadowMode = 0;
		else
			shadowMode = AiShaderEvalParamEnum(p_shadow_mode);
		
		switch ( shadowMode )
		{
			case sh_black:
				transparency = AI_RGB_BLACK;
				break;
			case sh_transmit_only:		
				transparency = cTransmission;
				break;
			case sh_transmit_and_outer_fresnels:
				if (validInterface && entering)
				{					
					fresnelTerm = fresnelEquations (n1, n2,  AiV3Dot(sg->Nf, -sg->Rd), cPolarizationTerm, false);
				}
				transparency = cTransmission * (AI_RGB_WHITE - fresnelTerm);
				break;
			case sh_transmit_and_all_fresnels:
				if (validInterface)
				{
					fresnelTerm = fresnelEquations (n1, n2,  AiV3Dot(sg->Nf, -sg->Rd), cPolarizationTerm, false);
				}
				transparency = cTransmission * (AI_RGB_WHITE - fresnelTerm);
				break;
		}

		AiShaderGlobalsApplyOpacity(sg, AI_RGB_WHITE - transparency);
		if (sg->out_opacity != AI_RGB_WHITE)
			updateMediaInsideLists(m_cMatID, entering, media_inside_ptr, false);

		return;
	} else {
		RayState->ray_energy *= cTransmission;
	}


	// ---------------------------------------------------//
	// Main Ray Tracing
	// ---------------------------------------------------//
	
	AtColor acc_refract_indirect = AI_RGB_BLACK;
	AtColor acc_refract_direct = AI_RGB_BLACK;
	AtColor acc_refract_direct_second = AI_RGB_BLACK;
	AtColor acc_spec_indirect = AI_RGB_BLACK;
	AtColor acc_spec_direct = AI_RGB_BLACK;

	if ( validInterface )
	{
		bool trace_refract_indirect = RayState->media_refractIndirect.v[m1] > ZERO_EPSILON && RayState->media_refractIndirect.v[m2] > ZERO_EPSILON;
		bool trace_refract_direct = RayState->media_refractDirect.v[m_cMatID] > ZERO_EPSILON ; // media exit is AND'd in later. 
		bool trace_spec_indirect = RayState->media_specIndirect.v[m_higherPriority] > ZERO_EPSILON;
		bool trace_spec_direct = mediaEntrance && RayState->media_specDirect.v[m_higherPriority] > ZERO_EPSILON;
		bool trace_TIR = trace_refract_indirect || trace_spec_indirect;

		float overallResultScale = 1.0f;

		bool photon_ray_type = false;
		AiStateGetMsgBool("photon_ray_type", &photon_ray_type);

		bool traceCaust_photon = false;
		bool traceCaust_pathtraced = false;

		if (!mediaEntrance && !AiShaderEvalParamBool(p_enable_internal_reflections) )
			trace_spec_indirect = false;

		bool causticPath = photon_ray_type || sg->Rr_diff > 0 ;
		if ( causticPath ) 
		{
			// caustic controls

			if (sg->Rt == AI_RAY_DIFFUSE)
			{
				const float caustic_max_distance = AiShaderEvalParamFlt(p_caustic_max_distance);
				if (caustic_max_distance > ZERO_EPSILON)
				{
					const float raylength = (float) sg->Rl;
					if ( raylength > caustic_max_distance ) 
					{
						trace_refract_indirect = trace_refract_direct = trace_spec_indirect = trace_spec_direct = trace_TIR = false;
					}
					else 
					{
						overallResultScale = 1.0f - ( raylength / caustic_max_distance );
					}
				}
			}

			if (mediaEntrance || !RayState->caustic_behaviorSet)
			{
				RayState->caustic_behaviorSet = true;
				RayState->caustic_mode = AiShaderEvalParamInt(p_caustic_mode);
				RayState->caustic_refractDirect = AiShaderEvalParamBool(p_caustic_refractions_direct);
				RayState->caustic_refractIndirect = AiShaderEvalParamBool(p_caustic_refractions_indirect);
				RayState->caustic_TIR = AiShaderEvalParamBool(p_caustic_TIR);
				RayState->caustic_specInternal = AiShaderEvalParamBool(p_caustic_internal_speculars);
				RayState->caustic_dispersion = AiShaderEvalParamBool(p_caustic_dispersion);
				RayState->caustic_specDirect = AiShaderEvalParamBool(p_caustic_specular_direct);
				RayState->caustic_specIndirect = AiShaderEvalParamBool(p_caustic_specular_indirect);
			}

			traceCaust_photon = (photon_ray_type && RayState->caustic_mode == 1);
			traceCaust_pathtraced = (sg->Rr_diff > 0 && RayState->caustic_mode == 0);
			causticPath = traceCaust_pathtraced || traceCaust_photon;

			if (traceCaust_pathtraced) //means we're in the right kind of caustics for the right kind of ray.
			{
				do_disperse = do_disperse && RayState->caustic_dispersion;

				trace_refract_direct = trace_refract_direct && RayState->caustic_refractDirect;
				trace_spec_indirect  = trace_spec_indirect  && RayState->caustic_specIndirect;
				trace_spec_direct    = trace_spec_direct    && RayState->caustic_specDirect;
				trace_TIR = trace_TIR && RayState->caustic_TIR && (RayState->caustic_refractDirect || RayState->caustic_refractIndirect); // TO DO: does that really work?

				if (mediaExit)
				{
					// Indirect refraction is only needed on mediaExit if we actually want caustics from indirect refractions.
					// Otherwise they are carriers for direct refraction. 
					trace_refract_indirect = trace_refract_indirect && RayState->caustic_refractIndirect;
				}
				else
				{
					// Indrect refraction is needed if we want caustics from indirect refraction, or direct refraction, 
					trace_refract_indirect = trace_refract_indirect && (RayState->caustic_refractDirect || RayState->caustic_refractIndirect);
				}

				if (mediaEntrance)
				{
					// if we have an inside out surface, TIR caustics will only be traced if indirect specular is on. 
					trace_TIR = trace_spec_indirect;
				}
				else
				{
					trace_spec_indirect = trace_spec_indirect && RayState->caustic_specInternal;
				}
			}
			else if (traceCaust_photon)
			{
				// in photon land, direct refractions means the light refracts straight through.
				// indirect refractions means.. very little. 
				// to the renderer this would always be indirect refraction though. 
				trace_refract_indirect = trace_refract_indirect && RayState->caustic_refractDirect;

				// in photon land, direct specular caustics will mean reflections off the outside of the surface.
				// indirect specular will mean the caustics will beounce off interior surfaces too. 
				// to the renderer this would always be indirect specular though. 
				trace_spec_indirect = trace_spec_indirect && (
					(RayState->caustic_specDirect && entering)
					|| (RayState->caustic_specIndirect)
					);

				trace_TIR = trace_TIR && RayState->caustic_TIR;
				trace_refract_direct = trace_spec_direct = false;
				do_disperse = do_disperse && RayState->caustic_dispersion;
			}
			else
			{
				trace_refract_indirect 
				= trace_refract_direct 
				= trace_spec_indirect 
				= trace_spec_direct 
				= trace_TIR
				= false;
			}
		}

		bool energySignificant;

		if (   std::abs( RayState->ray_energy.r ) < RayState->energy_cutoff
			&& std::abs( RayState->ray_energy.g ) < RayState->energy_cutoff
			&& std::abs( RayState->ray_energy.b ) < RayState->energy_cutoff
			)
		{
			energySignificant = false;
		}
		else
		{
			energySignificant = true;
		}

		trace_refract_direct = trace_refract_direct && mediaExit;

		const bool traceAnyRefraction = trace_refract_indirect || trace_refract_direct;
		const bool traceAnySpecular = trace_spec_indirect || trace_spec_direct;
		const bool traceAnything = (traceAnyRefraction || traceAnySpecular) && energySignificant;

		// ---------------------------------------------------//
		// Main Ray Tracing
		// Valid interfaces
		// ---------------------------------------------------//

		// decision point- trace anything
		if ( traceAnything ) 
		{
			AtRay ray;
			AtScrSample sample;
			
			float fresnelTerm = fresnelEquations (n1, n2,  AiV3Dot(sg->Nf, -sg->Rd), cPolarizationTerm, true);

			bool do_TIR = false;
			AtColor TIR_color = AI_RGB_BLACK;

			AtVector uTangent; 
			AtVector vTangent;
			AtVector tangentSourceVector;

			float spec_roughnessU = (RayState->media_specRoughnessU.v[m2] + RayState->media_specRoughnessU.v[m1]) / 2.0f;
			float spec_roughnessV = (RayState->media_specRoughnessV.v[m2] + RayState->media_specRoughnessV.v[m1]) / 2.0f;

			float refr_roughnessU = refractiveRoughness( RayState->media_refractRoughnessU.v[m1], RayState->media_refractRoughnessU.v[m2], n1, n2 );
			float refr_roughnessV = refractiveRoughness( RayState->media_refractRoughnessV.v[m1], RayState->media_refractRoughnessV.v[m2], n1, n2 );
				
			if ( traceAnyRefraction )
			{
				int refractSamplesTaken = 0;

				// ---------------------------------------------------//
				// Refraction
				// Samplers
				// ---------------------------------------------------//

				AtSamplerIterator* dispersionIterator = AiSamplerIterator( data->dispersion_sampler, sg) ;
				AtSamplerIterator* refractionIterator = AiSamplerIterator( data->refraction_sampler, sg);

				SAMPLETYPE dispersion_sample[2];
				SAMPLETYPE refraction_sample[2];

				// ---------------------------------------------------//
				// Refraction
				// BTDF preprocessing
				// ---------------------------------------------------//


				if (causticPath)
				{
					const float causticDRRoughnessOffset = refractRoughnessConvert( AiShaderEvalParamFlt(p_caustic_dr_roughness_offset) );
					refr_roughnessU += causticDRRoughnessOffset;
					refr_roughnessV += causticDRRoughnessOffset;
				}

				AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, &sg->Rd, AI_BIG, sg); 				
				const bool refracted = AiRefractRay(&ray, &sg->Nf, n1, n2, sg);

				AtShaderGlobals ppsg = *sg;
				if (!do_disperse)
					parallelPark(ray.dir, &ppsg);

				// decision point- indirect refraction
				if ( trace_refract_indirect )
				{
					const float cMediaIndirectRefractionProduct = RayState->media_refractIndirect.v[m1] * RayState->media_refractIndirect.v[m2];
					void * btdf_data = NULL;
					if (!do_disperse && do_blurryRefraction)
					{
						switch ( RayState->media_BTDF.v[m_higherPriority] )
						{
							case b_stretched_phong:
								// Stretched Phong
								btdf_data = AiStretchedPhongMISCreateData(&ppsg, (0.5f / SQR(refr_roughnessU) - 0.5f));
								break;
							case b_cook_torrance:
								// Cook Torrance
								btdf_data = AiCookTorranceMISCreateData(&ppsg, &AI_V3_ZERO, &AI_V3_ZERO, refr_roughnessU, refr_roughnessU);
								break;
							case b_ward_rayTangent:
								// Ward with refraction-derivitive tangents
								tangentSourceVector = AiV3Normalize(ppsg.Rd);
								blurAnisotropicPoles(&refr_roughnessU, &refr_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->N, &tangentSourceVector);
								uTangent = AiV3Cross( ppsg.Nf, tangentSourceVector ); 
								vTangent = AiV3Cross(ppsg.Nf, uTangent);
								btdf_data = AiWardDuerMISCreateData(&ppsg, &uTangent, &vTangent, refr_roughnessU, refr_roughnessV); 
								break;
							case b_ward_userTangent:
								// Ward with user tangents
								tangentSourceVector = AiV3Normalize( AiShaderEvalParamVec(p_ward_tangent) );
								blurAnisotropicPoles(&refr_roughnessU, &refr_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->N, &tangentSourceVector);
								uTangent = AiV3Cross( ppsg.Nf, tangentSourceVector ); 
								vTangent = AiV3Cross(ppsg.Nf, uTangent);
								btdf_data = AiWardDuerMISCreateData( &ppsg, &uTangent, &vTangent, refr_roughnessU, refr_roughnessV ) ; 
								break;
						}
					}
					// ---------------------------------------------------//
					// Refraction - Indirect
					// Indirect REfraction Sampling
					// ---------------------------------------------------//
					
					AtRay dispersalRay;
					float dispersal_seed = -1.0f;
					int dispersed_TIR_samples = 0;
					bool dispersed_TIR = false;

					const bool refract_skies = AiShaderEvalParamBool(p_refract_skies);

					while ( AiSamplerGetSample(refractionIterator, refraction_sample) )
					{
						AtColor monochromaticColor = AI_RGB_WHITE;
						bool refracted_dispersion = true;
						if (refracted || do_disperse)
						{
							if (do_disperse)
							{
								AiSamplerGetSample(dispersionIterator, dispersion_sample);

								if (dispersal_seed < 0.0f)
								{
									// The job of a dispersal seed is to fix any correlations, but still allow stratefied sampling to work in batches of samples.
									dispersal_seed =  ( std::abs( sg->sx + sg->sy ) * 113 + (float) dispersion_sample[1] ) * 3.456f  ;
								}
								const float LUT_value = fmod( (float) (dispersal_seed + dispersion_sample[0]), 1.0f );

								float cWavelength; 
								get_interpolated_LUT_value( RayState->spectral_LUT_, LUT_value, &cWavelength, &monochromaticColor);
								RayState->ray_monochromatic = true;
								RayState->ray_wavelength = cWavelength;

								const float n1_dispersed = dispersedIOR(n1, RayState->media_dispersion.v[m1], cWavelength);
								const float n2_dispersed = dispersedIOR(n2, RayState->media_dispersion.v[m2], cWavelength);

								AiMakeRay(&dispersalRay, AI_RAY_REFRACTED, &sg->P, &sg->Rd, AI_BIG, sg); 								
								refracted_dispersion = AiRefractRay(&dispersalRay, &sg->Nf, n1_dispersed, n2_dispersed, sg);
								ray.dir = dispersalRay.dir;

								if ( do_blurryRefraction )
								{									
									parallelPark(ray.dir, &ppsg);
									switch ( RayState->media_BTDF.v[m_higherPriority] )
									{
										case b_stretched_phong:
											// Stretched Phong
											btdf_data = AiStretchedPhongMISCreateData(&ppsg, (0.5f / SQR(refr_roughnessU) - 0.5f));
											break;
										case b_cook_torrance:
											// Cook Torrance
											btdf_data = AiCookTorranceMISCreateData(&ppsg, &AI_V3_ZERO, &AI_V3_ZERO, refr_roughnessU, refr_roughnessU);
											break;
										case b_ward_rayTangent:
											// Ward with refraction-derivitive tangents
											tangentSourceVector = AiV3Normalize(ppsg.Rd);
											blurAnisotropicPoles(&refr_roughnessU, &refr_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->N, &tangentSourceVector);
											uTangent = AiV3Cross( ppsg.Nf, tangentSourceVector ); 
											vTangent = AiV3Cross(ppsg.Nf, uTangent);
											btdf_data = AiWardDuerMISCreateData(&ppsg, &uTangent, &vTangent, refr_roughnessU, refr_roughnessV); 
											break;
										case b_ward_userTangent:
											// Ward with user tangents
											tangentSourceVector = AiV3Normalize( AiShaderEvalParamVec(p_ward_tangent) );
											blurAnisotropicPoles(&refr_roughnessU, &refr_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->N, &tangentSourceVector);
											uTangent = AiV3Cross( ppsg.Nf, tangentSourceVector ); 
											vTangent = AiV3Cross(ppsg.Nf, uTangent);
											btdf_data = AiWardDuerMISCreateData( &ppsg, &uTangent, &vTangent, refr_roughnessU, refr_roughnessV ) ; 
											break;
									}
								}
							}

							if ( do_blurryRefraction )
							{
								switch ( RayState->media_BTDF.v[m_higherPriority] )
								{
									case b_stretched_phong:
										ray.dir = AiStretchedPhongMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);
										break;
									case b_cook_torrance:
										ray.dir = AiCookTorranceMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);
										break;
									case b_ward_rayTangent:
										ray.dir = AiWardDuerMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);
										break;
									case b_ward_userTangent:
										ray.dir = AiWardDuerMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);
										break;
								}
							}

							if (!refracted_dispersion)
							{
								TIR_color += monochromaticColor;
								dispersed_TIR = true;
								dispersed_TIR_samples++;
								refractSamplesTaken++ ;
							}
							else if ((AiV3Dot(ray.dir,sg->Nf) < 0.0f) && (ray.dir != AI_V3_ZERO))
							{
								updateMediaInsideLists(m_cMatID, entering, media_inside_ptr, false);
								const AtColor weight = (1.0f - fresnelTerm) 
										* monochromaticColor
										* cMediaIndirectRefractionProduct
										* overallResultScale;

								const AtColor energyCache = RayState->ray_energy;
								RayState->ray_energy *= weight;

								refractSamplesTaken++ ;
								const bool tracehit = AiTrace(&ray, &sample);
								if (tracehit || refract_skies) 
								{
									acc_refract_indirect += sample.color * weight * transmissionOnSample(&t2, &sample, tracehit );
								}

								RayState->ray_energy = energyCache;
								updateMediaInsideLists(m_cMatID, entering, media_inside_ptr, true);
							}

							if (!do_multiSampleRefraction)
								break;
						}
					}
					
					/*
					 * Refraction - TIR
					 * Total internal reflection sampling
					 * TIR counts as refraction rays, and depth is counted seperately by RayState->ray_TIRDepth
					 * Actually occurs in the specular loop below though. 
					 */

					// decision point - tir
					if (do_disperse)
					{
						if (dispersed_TIR)
						{
							TIR_color *= (float) AiSamplerGetSampleInvCount(refractionIterator);
							do_TIR = true;
						}
					}
					else
					{						
						if ( trace_TIR && !refracted && !do_disperse )
						{
							TIR_color = AI_RGB_WHITE;
							do_TIR = true;
						}
					}

					if (refractSamplesTaken > 0) 
						acc_refract_indirect /= (float) refractSamplesTaken;
					else 
						acc_refract_indirect = AI_RGB_BLACK;  //to do: necessary line?
				}

				// Exhaust the sampler. Good night, sampler. This seems necessary, having the unexhausted sampler caused some problems with the RR sampler. 
				while ( AiSamplerGetSample( dispersionIterator, dispersion_sample ) ){}

				// ---------------------------------------------------//
				// Refraction - Direct
				// Direct Refraction Sampling
				// ---------------------------------------------------//

				// decision point- direct refraction
				if (trace_refract_direct && refracted )
				{
					const bool use_refraction_btdf = AiShaderEvalParamBool( p_dr_use_refraction_btdf );
					float dr_roughnessU;
					float dr_roughnessV;
					int dr_btdf;

					if (use_refraction_btdf)
					{
						dr_roughnessU = refr_roughnessU;
						dr_roughnessV = refr_roughnessV;
						dr_btdf = RayState->media_BTDF.v[m_cMatID];
					}
					else
					{
						// manual mode
						dr_roughnessU = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughness_u ) );
						dr_roughnessV = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughness_v ) );
						dr_btdf =  AiShaderEvalParamEnum( p_dr_btdf );
					}

					// offsets and depth modification
					const float dr_roughnessOffset = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughnessOffset ) );
					const float dr_roughnessDepthAdder = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughnessDepthAdder ) );
					const float dr_roughnessDepthMultiplier = AiShaderEvalParamFlt( p_dr_roughnessDepthMultiplier );

					dr_roughnessU =  (dr_roughnessU) * pow(dr_roughnessDepthMultiplier, (float) sg->Rr) + (dr_roughnessDepthAdder * (float) sg->Rr) + dr_roughnessOffset;
					dr_roughnessV =  (dr_roughnessV) * pow(dr_roughnessDepthMultiplier, (float) sg->Rr) + (dr_roughnessDepthAdder * (float) sg->Rr) + dr_roughnessOffset;

					float drs_roughnessU;
					float drs_roughnessV;
					float dr_first_scale = RayState->media_refractDirect.v[m_cMatID];
					float dr_second_scale = AiShaderEvalParamFlt( p_dr_second_scale );
					if (dr_second_scale > ZERO_EPSILON)
					{
						const float dr_second_roughnessMultiplier = AiShaderEvalParamFlt( p_dr_second_roughnessMultiplier );
						drs_roughnessU = dr_roughnessU * dr_second_roughnessMultiplier;
						drs_roughnessV = dr_roughnessV * dr_second_roughnessMultiplier;
						dr_first_scale = (1.0f - dr_second_scale) * RayState->media_refractDirect.v[m_cMatID];
						dr_second_scale = dr_second_scale * RayState->media_refractDirect.v[m_cMatID];
					}

					void * btdf_data_direct = NULL;
					void * btdf_data_direct2 = NULL;

					switch ( dr_btdf )
					{
						case b_stretched_phong:
							// Stretched Phong
							btdf_data_direct = AiStretchedPhongMISCreateData(&ppsg, (0.5f / SQR(dr_roughnessU) - 0.5f));
							btdf_data_direct2 = AiStretchedPhongMISCreateData(&ppsg, (0.5f / SQR(drs_roughnessU) - 0.5f));
							break;
						case b_cook_torrance:
							// Cook Torrance
							btdf_data_direct = AiCookTorranceMISCreateData(&ppsg, &AI_V3_ZERO, &AI_V3_ZERO, dr_roughnessU, dr_roughnessU);
							btdf_data_direct2 = AiCookTorranceMISCreateData(&ppsg, &AI_V3_ZERO, &AI_V3_ZERO, drs_roughnessU, drs_roughnessU);
							break;
						case b_ward_rayTangent:
							// Ward with refraction-derivitive tangents
							tangentSourceVector = AiV3Normalize(ppsg.Rd);
							blurAnisotropicPoles(&dr_roughnessU, &dr_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->Nf, &tangentSourceVector);
							blurAnisotropicPoles(&drs_roughnessU, &drs_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->Nf, &tangentSourceVector);
							uTangent = AiV3Cross( ppsg.Nf, tangentSourceVector ); 
							vTangent = AiV3Cross(ppsg.Nf, uTangent);
							btdf_data_direct = AiWardDuerMISCreateData(&ppsg, &uTangent, &vTangent, dr_roughnessU, dr_roughnessV); 
							btdf_data_direct2 = AiWardDuerMISCreateData(&ppsg, &uTangent, &vTangent, drs_roughnessU, drs_roughnessV); 
							break;
						case b_ward_userTangent:
							// Ward with user tangents
							tangentSourceVector = AiV3Normalize( AiShaderEvalParamVec(p_ward_tangent) );
							blurAnisotropicPoles(&dr_roughnessU, &dr_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->Nf, &tangentSourceVector);
							blurAnisotropicPoles(&drs_roughnessU, &drs_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->Nf, &tangentSourceVector);
							uTangent = AiV3Cross( ppsg.Nf, tangentSourceVector ); 
							vTangent = AiV3Cross(ppsg.Nf, uTangent);					
							btdf_data_direct = AiWardDuerMISCreateData( &ppsg, &uTangent, &vTangent, dr_roughnessU, dr_roughnessV ) ; 
							btdf_data_direct2 = AiWardDuerMISCreateData( &ppsg, &uTangent, &vTangent, drs_roughnessU, drs_roughnessV ) ; 
							break;
					}
					AiStateSetMsgBool("opaqueShadowMode", true);
					AiLightsPrepare(&ppsg);
					const bool refract_skydomes = AiShaderEvalParamBool(p_refract_skydomes);
					while (AiLightsGetSample(&ppsg))
					{
						if ( refract_skydomes || !AiNodeIs( ppsg.Lp,"skydome_light" ))
						{
							switch ( dr_btdf )
							{
								case b_stretched_phong:
									acc_refract_direct += AiEvaluateLightSample(&ppsg, btdf_data_direct, AiStretchedPhongMISSample, AiStretchedPhongMISBRDF, AiStretchedPhongMISPDF);
									if (dr_second_scale > ZERO_EPSILON)
									{
										acc_refract_direct_second += AiEvaluateLightSample(&ppsg, btdf_data_direct2, AiStretchedPhongMISSample, AiStretchedPhongMISBRDF, AiStretchedPhongMISPDF);
									}
									break;
								case b_cook_torrance:
									acc_refract_direct += AiEvaluateLightSample(&ppsg, btdf_data_direct, AiCookTorranceMISSample, AiCookTorranceMISBRDF, AiCookTorranceMISPDF);
									if (dr_second_scale > ZERO_EPSILON)
									{
										acc_refract_direct_second += AiEvaluateLightSample(&ppsg, btdf_data_direct2, AiCookTorranceMISSample, AiCookTorranceMISBRDF, AiCookTorranceMISPDF);
									}
									break;
								case b_ward_rayTangent:
									acc_refract_direct += AiEvaluateLightSample(&ppsg, btdf_data_direct, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
									if (dr_second_scale > ZERO_EPSILON)
									{
										acc_refract_direct_second += AiEvaluateLightSample(&ppsg, btdf_data_direct2, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
									}
									break;
								case b_ward_userTangent:
									acc_refract_direct += AiEvaluateLightSample(&ppsg, btdf_data_direct, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
									if (dr_second_scale > ZERO_EPSILON)
									{
										acc_refract_direct_second += AiEvaluateLightSample(&ppsg, btdf_data_direct2, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
									}
									break;
							}
						}
						AiStateSetMsgBool("opaqueShadowMode", false);
					}
					acc_refract_direct *= dr_first_scale * (1.0f - fresnelTerm) * overallResultScale;
					acc_refract_direct_second *= dr_second_scale * (1.0f - fresnelTerm) * overallResultScale;
				}				
			}


			// ---------------------------------------------------//
			// (Weak) Specular / Reflection
			// ---------------------------------------------------//

			const float rrProbability = AiShaderEvalParamFlt(p_russian_roulette_probability);
			const bool rrActive = rrProbability > ZERO_EPSILON && (sg->Rr_refr > 0) && !do_TIR ;

			if ( rrActive )
			{
				AtSamplerIterator* russian_roullete_iterator = AiSamplerIterator( data->russian_roullete_single_sampler, sg );
				fresnelTerm = russianRoulette( russian_roullete_iterator, fresnelTerm, rrProbability);
			}

			if (do_TIR)
				fresnelTerm = 1.0f;

			if (fresnelTerm > ZERO_EPSILON)
			{
				// reset data to do reflection
				if (do_disperse)
				{
					// reflected rays off a dispersive medium aren't monochromatic
					RayState->ray_monochromatic = false;
				}

			// decision point- any specular
			if ( traceAnySpecular || do_TIR ) 
			{
					SAMPLETYPE specular_sample[2];
					AtSamplerIterator* specularIterator = AiSamplerIterator( data->specular_sampler, sg);
					AtRay specularRay;
					AtUInt32 rayType;

					if (AiShaderEvalParamInt(p_specular_ray_type) == 0)  
						rayType = AI_RAY_GLOSSY;
					else if (AiShaderEvalParamInt(p_specular_ray_type) == 1)  
						rayType = AI_RAY_REFLECTED;
					if (do_TIR) 
						rayType = AI_RAY_REFRACTED;

					AiMakeRay(&specularRay, rayType, &sg->P, NULL, AI_BIG, sg);

					void * brdf_data; 

					switch ( RayState->media_BRDF.v[m_higherPriority] )
					{						
						case b_stretched_phong:
							// Stretched Phong
							brdf_data = AiStretchedPhongMISCreateData(sg, (0.5f / SQR(spec_roughnessU) - 0.5f));
							break;
						case b_cook_torrance:
							// Cook Torrance
							brdf_data = AiCookTorranceMISCreateData(sg, &AI_V3_ZERO, &AI_V3_ZERO, spec_roughnessU, spec_roughnessU);
							break;
						case b_ward_rayTangent:
							// Ward with refraction-derivitive tangents
							blurAnisotropicPoles(&spec_roughnessU, &spec_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->Nf, &tangentSourceVector);
							AiV3Cross(uTangent, sg->N, sg->Rd); 
							AiV3Cross(vTangent, sg->N, uTangent);
							brdf_data = AiWardDuerMISCreateData(sg, &uTangent, &vTangent, spec_roughnessU, spec_roughnessV); 
							break;
						case b_ward_userTangent:
							// Ward with user tangents
							tangentSourceVector = AiV3Normalize( AiShaderEvalParamVec(p_ward_tangent) );
							blurAnisotropicPoles(&spec_roughnessU, &spec_roughnessV, &RayState->media_blurAnisotropicPoles.v[m_higherPriority], &sg->Nf, &tangentSourceVector);
							AiV3Cross( uTangent, sg->N, tangentSourceVector ) ;
							AiV3Cross( vTangent, sg->N, uTangent ) ;
							brdf_data = AiWardDuerMISCreateData( sg, &vTangent, &uTangent, spec_roughnessU, spec_roughnessV ) ; 
							break;
					}

					// ---------------------------------------------------//
					// Indirect Specular / Reflection
					// ---------------------------------------------------//

					// decision point- indirect specular
					if ( trace_spec_indirect || do_TIR )
					{
						const float weight = fresnelTerm * RayState->media_specIndirect.v[m1] * overallResultScale;
						const AtColor energyCache = RayState->ray_energy;						
						const bool reflect_skies = AiShaderEvalParamBool(p_reflect_skies);

						if ( do_TIR )
						{
							RayState->ray_TIRDepth++ ;

							if (RayState->ray_TIRDepth < 50 && specularRay.refr_bounces > 1)
								specularRay.refr_bounces-- ;

							RayState->ray_energy *= TIR_color;
						} 
						else
						{
							RayState->ray_energy *= weight;		
						}
						if (spec_roughnessU > ZERO_EPSILON || spec_roughnessV > ZERO_EPSILON)
						{
							while ( AiSamplerGetSample(specularIterator, specular_sample) )
							{
								switch ( RayState->media_BRDF.v[m_higherPriority] )
								{
									case b_stretched_phong:
										specularRay.dir = AiStretchedPhongMISSample(brdf_data, (float) specular_sample[0], (float) specular_sample[1]);
										break;
									case b_cook_torrance:
										specularRay.dir = AiCookTorranceMISSample(brdf_data, (float) specular_sample[0], (float) specular_sample[1]);
										break;
									case b_ward_rayTangent:
										specularRay.dir = AiWardDuerMISSample(brdf_data, (float) specular_sample[0], (float) specular_sample[1]);
										break;
									case b_ward_userTangent:
										specularRay.dir = AiWardDuerMISSample(brdf_data, (float) specular_sample[0], (float) specular_sample[1]);
										break;
								}
								
								if (AiV3Dot(specularRay.dir,sg->Nf) > ZERO_EPSILON )
								{									

									const bool tracehit = AiTrace(&specularRay, &sample);
									if (tracehit || reflect_skies) 
									{
										if (do_TIR) 
											acc_spec_indirect += sample.color * TIR_color * transmissionOnSample(&t1, &sample, tracehit );
										else 
											acc_spec_indirect += sample.color * weight * transmissionOnSample(&t1, &sample, tracehit );
									}
								}
							}
							acc_spec_indirect *= (float) AiSamplerGetSampleInvCount(specularIterator);

						}
						else
						{
							AiReflectRay(&specularRay, &sg->Nf, sg);
							const bool tracehit = AiTrace(&specularRay, &sample);
							if (tracehit || reflect_skies) 
							{
								if (do_TIR) 
									acc_spec_indirect += sample.color * TIR_color * overallResultScale * transmissionOnSample(&t1, &sample, tracehit );
								else 
									acc_spec_indirect = sample.color * fresnelTerm * RayState->media_specIndirect.v[m1] * overallResultScale * transmissionOnSample(&t1, &sample, tracehit );
							}
						}
						RayState->ray_energy = energyCache;						
					}

					// ---------------------------------------------------//
					// Direct Specular / Reflection
					// ---------------------------------------------------//

					// decision point- direct specular
					if ( trace_spec_direct )
					{
						AtColor weight = AI_RGB_WHITE * fresnelTerm * RayState->media_specDirect.v[m_higherPriority] * overallResultScale;
						if (do_TIR)
							weight *= TIR_color;
						
						AiStateSetMsgBool("opaqueShadowMode", true);						
						AiLightsPrepare(sg); 
						const bool reflect_skydomes = AiShaderEvalParamBool(p_reflect_skydomes);
						while ( AiLightsGetSample(sg) ) // loop over the lights to compute direct effects
						{
							float l_weight = AiLightGetSpecular(sg->Lp);

							if ( 
								(reflect_skydomes || !AiNodeIs( sg->Lp,"skydome_light" )) && 
								l_weight > ZERO_EPSILON
								)
							{
								switch ( RayState->media_BRDF.v[m_higherPriority] )
								{
									case b_stretched_phong:
										acc_spec_direct += l_weight * AiEvaluateLightSample(sg, brdf_data, AiStretchedPhongMISSample, AiStretchedPhongMISBRDF, AiStretchedPhongMISPDF);
										break;
									case b_cook_torrance:
										acc_spec_direct += l_weight * AiEvaluateLightSample(sg, brdf_data, AiCookTorranceMISSample, AiCookTorranceMISBRDF, AiCookTorranceMISPDF);
										break;
									case b_ward_rayTangent:
										acc_spec_direct += l_weight * AiEvaluateLightSample(sg, brdf_data, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
										break;
									case b_ward_userTangent:
										acc_spec_direct += l_weight * AiEvaluateLightSample(sg, brdf_data, AiWardDuerMISSample, AiWardDuerMISBRDF, AiWardDuerMISPDF);
										break;
								}								
							}
						}
						AiStateSetMsgBool("opaqueShadowMode", false);
						acc_spec_direct *= weight;
					}
				}
			}
		}

	} // End of ordinary (non-shadow) ray section

	// ---------------------------------------------------//
	// - Invalid Interface Tracing
	// ---------------------------------------------------//

	AtRGBA invalidInterfaceResult = AI_RGBA_BLACK;

	if ( !validInterface )
	{
		if ( RayState->ray_invalidDepth < 70 )
		{
			AtRay ray;
			AtScrSample sample;

			RayState->ray_invalidDepth ++;
			updateMediaInsideLists(m_cMatID, entering, media_inside_ptr, false);

			AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, NULL, AI_BIG, sg);
			ray.dir = sg->Rd;
			ray.level --;
			ray.refr_bounces --;
			const bool tracehit = AiTrace(&ray, &sample);

			AiRGBtoRGBA( sample.color * transmissionOnSample(&t2, &sample, tracehit ), invalidInterfaceResult );
			invalidInterfaceResult.a = sample.alpha;			
		}
		else
		{
			char * const overlapping_surfaces_message = "JF Nested Dielectric: Crazy numbers of invalid interfaces. Some geo may be overlapping.";
			AiMsgWarning(overlapping_surfaces_message);
			sg->out.RGBA = AI_RGBA_BLACK;
			return;
		}
	}

	// ---------------------------------------------------//
	// - Emission
	// ---------------------------------------------------//

	AtColor emission = AI_RGB_BLACK;
	if (validInterface)
	{
		AiStateSetMsgBool("opaqueShadowMode", true);
		emission += AiShaderEvalParamRGB(p_emission_at_interfaces);
		if (mediaExit)
		{
			emission += AiShaderEvalParamRGB(p_emission_at_exits);
		}
		else if (mediaEntrance)
		{
			emission += AiShaderEvalParamRGB(p_emission_at_entrances);
		}
		AiStateSetMsgBool("opaqueShadowMode", false);
	}

	// ---------------------------------------------------//
	// - Output and AOVs
	// ---------------------------------------------------//	

	if ( validInterface )
	{
		AiAOVSetRGBA(sg, data->aov_direct_refraction.c_str(), AiRGBtoRGBA( acc_refract_direct + acc_refract_direct_second ));			
		AiAOVSetRGBA(sg, data->aov_direct_specular.c_str(), AiRGBtoRGBA( acc_spec_direct ));			
		AiAOVSetRGBA(sg, data->aov_indirect_refraction.c_str(), AiRGBtoRGBA( acc_refract_indirect ));			
		AiAOVSetRGBA(sg, data->aov_indirect_specular.c_str(), AiRGBtoRGBA( acc_spec_indirect ));			
		AtRGBA outcolor;
		AiRGBtoRGBA(
			( acc_refract_indirect 
			+ acc_refract_direct 
			+ acc_refract_direct_second 
			+ acc_spec_indirect
			+ acc_spec_direct
			+ emission ), 
			outcolor) ;
		sg->out.RGBA = outcolor;
	}
	else
	{
		sg->out.RGBA = invalidInterfaceResult;
	}

	// ---------------------------------------------------//
	// - Array and messaging post process 
	// Reset all messages to the previous state. If there was no previous state, Mark them invalid. 
	// ---------------------------------------------------//

	AiStateSetMsgRGB("photon_energy",RayState->ray_energy); // this doesn't seem right, heberlein caustics related only

	if (msgs_are_valid)
	{
		AiStateSetMsgBool("msgs_are_valid", true);
		uncacheRayState( RayState, &RayStateCache );
	}
	else
	{
		AiStateSetMsgBool("msgs_are_valid", false);
	}
}

node_loader
{
	if (i > 0)
		return false;

	node->methods		   = jf_nested_dielectric_methods;
	node->output_type  = AI_TYPE_RGBA;
	node->name					  = "jf_nested_dielectric";
	node->node_type  = AI_NODE_SHADER;
	strcpy(node->version, AI_VERSION);
	return true;
}
