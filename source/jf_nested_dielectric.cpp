/* 
 * JF Nested Dielectric by Jonah Friedman
 * 1.0.4
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
	JFND_Shader_Data* data = new JFND_Shader_Data();
	AiNodeSetLocalData(node, data);
};

node_finish
{
	if (AiNodeGetLocalData(node) != NULL)
	{
		JFND_Shader_Data *data = (JFND_Shader_Data*) AiNodeGetLocalData(node);
		AiNodeSetLocalData(node, NULL);
		delete data;
	}
}


node_update
{
	AtNode* render_options = AiUniverseGetOptions();

	JFND_Shader_Data *data = (JFND_Shader_Data*)AiNodeGetLocalData(node);
	data->update(node);
}




shader_evaluate
{
	JFND_Shader_Data *data = (JFND_Shader_Data*)AiNodeGetLocalData(node);
	Ray_State *RayState;
	Ray_State_Cache RayStateCache;

	// ---------------------------------------------------//
	// - Array and messaging preprocess 
	// ---------------------------------------------------//

	bool msgs_are_valid = false;
	AiStateGetMsgBool("msgs_are_valid", &msgs_are_valid);

	if ( msgs_are_valid )
	{
		void * rayState_ptr;
		AiStateGetMsgPtr("rayState_ptr", &rayState_ptr);
		RayState = static_cast<Ray_State*>( rayState_ptr );
		RayState->cacheRayState( &RayStateCache);
	}
	else
	{		
		RayState = static_cast<Ray_State*>( AiShaderGlobalsQuickAlloc(sg, sizeof( Ray_State ) ) );
		RayState->init(data, sg, node);

		RayState->setEnergyCutoff( (float) AiShaderEvalParamInt(p_energy_cutoff_exponent));
		RayState->setPolarization(AiShaderEvalParamBool(p_polarize), data->polarizationVector);

		AiStateSetMsgPtr("rayState_ptr", RayState);
	}
	AiStateSetMsgBool("msgs_are_valid", true); // Any child rays from this will find valid messages. 


	MediaIntStruct * media_inside_ptr;
	
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
			memcpy(&RayState->shadow_media_inside, &RayState->media_inside, sizeof(MediaIntStruct) );
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
		RayState->setMediaIOR(m_cMatID, AiShaderEvalParamFlt(p_mediumIOR));
		RayState->setMediaDispersion(m_cMatID, AiShaderEvalParamBool(p_disperse), AiShaderEvalParamFlt(p_dispersion));	
		RayState->setAnisotropySettings(m_cMatID, AiShaderEvalParamFlt(p_blur_anisotropic_poles));
			
		const float sScale = AiShaderEvalParamFlt(p_specular_scale);
		const float sDirect = AiShaderEvalParamFlt(p_direct_specular) * sScale;
		const float sIndirect = AiShaderEvalParamFlt(p_indirect_specular) * sScale;
		RayState->setMediaSpecular(m_cMatID, sDirect, sIndirect, AiShaderEvalParamEnum(p_brdf), 
			AiShaderEvalParamFlt(p_specular_roughness_u), AiShaderEvalParamFlt(p_specular_roughness_v));

		const float rScale = AiShaderEvalParamFlt(p_refraction_scale);
		const float rDirect = AiShaderEvalParamFlt(p_direct_refraction) * rScale;
		const float rIndirect = AiShaderEvalParamFlt(p_indirect_refraction) * rScale;
		const float rURough = refractRoughnessConvert( AiShaderEvalParamFlt(p_refraction_roughness_u) );
		const float rVRough = refractRoughnessConvert( AiShaderEvalParamFlt(p_refraction_roughness_v) );
		RayState->setRefractionSettings(m_cMatID, rDirect, rIndirect, AiShaderEvalParamEnum(p_btdf), rURough, rVRough, 
			AiShaderEvalParamRGB(p_mediumTransmittance), AiShaderEvalParamFlt(p_mediumTransmittance_scale));
	}


	// ---------------------------------------------------//
	// - get interface info     
	// ---------------------------------------------------//

	InterfaceInfo iinfo = InterfaceInfo( RayState, m_cMatID, sg);

	bool do_blurryRefraction = iinfo.doBlurryRefraction();
	bool do_disperse = iinfo.setupDispersion(data);
	const bool do_multiSampleRefraction = do_disperse || do_blurryRefraction;

	// ---------------------------------------------------//
	// - Shadow rays
	// ---------------------------------------------------//
	
	if (sg->Rt == AI_RAY_SHADOW)
	{
		AtColor transparency = iinfo.getShadowTransparency(sg, AiShaderEvalParamEnum(p_shadow_mode));
		AiShaderGlobalsApplyOpacity(sg, AI_RGB_WHITE - transparency);
		if (sg->out_opacity != AI_RGB_WHITE)
			updateMediaInsideLists(m_cMatID, iinfo.entering, media_inside_ptr, false);
		return;
	} else {
		AtColor cTransmission = iinfo.getTransmissionColor(sg);
		RayState->ray_energy *= cTransmission;
		RayState->ray_energy_photon *= cTransmission;
	}


	// ---------------------------------------------------//
	// param caching
	//     caching select parameters if it improves the code to do so
	// ---------------------------------------------------//
	const AtVector pval_custom_tangent = AiShaderEvalParamVec(p_ward_tangent);
	const int pval_specular_ray_type = AiShaderEvalParamInt(p_specular_ray_type);
	const bool pval_enable_internal_reflections = AiShaderEvalParamBool(p_enable_internal_reflections);

	// ---------------------------------------------------//
	// Main Ray Tracing
	// ---------------------------------------------------//
	
	AtColor acc_refract_indirect = AI_RGB_BLACK;
	AtColor acc_refract_direct = AI_RGB_BLACK;
	AtColor acc_refract_direct_second = AI_RGB_BLACK;
	AtColor acc_spec_indirect = AI_RGB_BLACK;
	AtColor acc_spec_direct = AI_RGB_BLACK;

	if ( iinfo.validInterface )
	{
		TraceSwitch traceSwitch = TraceSwitch(&iinfo);
		traceSwitch.setInternalReflections(&iinfo, pval_enable_internal_reflections);

		float overallResultScale = 1.0f;
		const bool photon_ray_type = rayIsPhoton(sg);
		const bool causticPath = photon_ray_type || sg->Rr_diff > 0 ;
		if ( causticPath ) 
		{
			if (iinfo.mediaEntrance || !RayState->caustic_behaviorSet)
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

			// caustic controls
			if (sg->Rt == AI_RAY_DIFFUSE)
			{
				const float caustic_max_distance = AiShaderEvalParamFlt(p_caustic_max_distance);
				if (caustic_max_distance > ZERO_EPSILON)
				{
					const float raylength = (float) sg->Rl;
					if ( raylength > caustic_max_distance ) 
					{
						traceSwitch.setTraceNone();
					}
					else 
					{
						overallResultScale = 1.0f - ( raylength / caustic_max_distance );
					}
				}
			}


			if ((sg->Rr_diff > 0 && RayState->caustic_mode == 0)) //means we're in the right kind of caustics for the right kind of ray.
			{
				do_disperse = do_disperse && RayState->caustic_dispersion;
				traceSwitch.setPathtracedCaustic(&iinfo);
			}
			else if ((photon_ray_type))
			{
				do_disperse = do_disperse && RayState->caustic_dispersion;
				traceSwitch.setPhotonCaustic(&iinfo);
			}
			else
			{
				traceSwitch.setTraceNone();
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

		// ---------------------------------------------------//
		// Main Ray Tracing
		// Valid interfaces
		// ---------------------------------------------------//

		// decision point- trace anything
		if ( energySignificant && traceSwitch.traceAnything() ) 
		{
			AtRay ray;
			AtScrSample sample;
			
			float fresnelTerm = fresnelEquations (iinfo.n1, iinfo.n2,  AiV3Dot(sg->Nf, -sg->Rd), iinfo.polarizationTerm, true);

			bool do_TIR = false;
			AtColor TIR_color = AI_RGB_BLACK;

			AtVector uTangent; 
			AtVector vTangent;
			AtVector tangentSourceVector;

			if ( traceSwitch.traceAnyRefraction() )
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

				float refr_roughnessU = iinfo.getRefrRoughnessU();
				float refr_roughnessV = iinfo.getRefrRoughnessV();

				if (causticPath)
				{
					const float causticDRRoughnessOffset = refractRoughnessConvert( AiShaderEvalParamFlt(p_caustic_dr_roughness_offset) );
					refr_roughnessU += causticDRRoughnessOffset;
					refr_roughnessV += causticDRRoughnessOffset;
				}

				AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, &sg->Rd, AI_BIG, sg); 				
				const bool refracted = AiRefractRay(&ray, &sg->Nf, iinfo.n1, iinfo.n2, sg);

				AtVector old_Rd = sg->Rd;
				AtShaderGlobals ppsg = *sg;
				if (!do_disperse)
					parallelPark(ray.dir, &ppsg);

				// decision point- indirect refraction
				// if ( trace_refract_indirect )
				if ( traceSwitch.refr_ind )
				{
					void * btdf_data = NULL;
					if (!do_disperse && do_blurryRefraction)
					{
						btdf_data = iinfo.getBTDFData(&ppsg, refr_roughnessU, refr_roughnessV, pval_custom_tangent);
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
								get_interpolated_LUT_value( RayState->spectral_LUT_ptr, LUT_value, &cWavelength, &monochromaticColor);
								RayState->ray_monochromatic = true;
								RayState->ray_wavelength = cWavelength;

								const float n1_dispersed = dispersedIOR(iinfo.n1, RayState->media_dispersion.v[iinfo.m1], cWavelength);
								const float n2_dispersed = dispersedIOR(iinfo.n2, RayState->media_dispersion.v[iinfo.m2], cWavelength);

								AiMakeRay(&dispersalRay, AI_RAY_REFRACTED, &sg->P, &sg->Rd, AI_BIG, sg); 								
								refracted_dispersion = AiRefractRay(&dispersalRay, &sg->Nf, n1_dispersed, n2_dispersed, sg);
								ray.dir = dispersalRay.dir;

								if ( do_blurryRefraction )
								{									
									parallelPark(ray.dir, &ppsg);
									btdf_data = iinfo.getBTDFData(&ppsg, refr_roughnessU, refr_roughnessV, pval_custom_tangent);
								}
							}

							if ( do_blurryRefraction )
							{
								ray.dir = iinfo.getBTDFSample(btdf_data, refraction_sample);
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
								updateMediaInsideLists(m_cMatID, iinfo.entering, media_inside_ptr, false);
								const AtColor weight = (1.0f - fresnelTerm) 
									* monochromaticColor
									* RayState->media_refractIndirect.v[iinfo.m1]
									* RayState->media_refractIndirect.v[iinfo.m2]
									* overallResultScale;

								const AtColor energyCache = RayState->ray_energy;
								RayState->ray_energy *= weight;
								const AtColor energyCache_photon = RayState->ray_energy_photon;
								RayState->ray_energy_photon *= weight;

								if (sg->Rt == AI_RAY_CAMERA && (do_disperse || do_blurryRefraction))								
									RayState->ray_energy_photon /= (float) data->refr_samples;

								AiStateSetMsgRGB("photon_energy",RayState->ray_energy_photon);

								refractSamplesTaken++ ;
								const bool tracehit = AiTrace(&ray, &sample);
								if (tracehit || refract_skies) 
								{
									acc_refract_indirect += sample.color * weight * transmissionOnSample(&iinfo.t2, &sample, tracehit );
								}

								RayState->ray_energy = energyCache;
								RayState->ray_energy_photon = energyCache_photon;

								updateMediaInsideLists(m_cMatID, iinfo.entering, media_inside_ptr, true);
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
						// if ( trace_TIR && !refracted && !do_disperse )
						if ( traceSwitch.TIR && !refracted && !do_disperse )
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
				// if (trace_refract_direct && refracted )
				if (traceSwitch.refr_dir && refracted )
				{
					// offsets and depth modification
					const float rOffset = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughnessOffset ) );
					const float rDepthAdder = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughnessDepthAdder ) );
					const float rDepthMultiplier = AiShaderEvalParamFlt( p_dr_roughnessDepthMultiplier );

					const bool use_refr_settings = AiShaderEvalParamBool( p_dr_use_refraction_btdf );
					int dr_btdf 		= use_refr_settings ? iinfo.getBTDFType() : AiShaderEvalParamEnum( p_dr_btdf );
					float dr_roughnessU = use_refr_settings ? refr_roughnessU : refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughness_u ) );
					float dr_roughnessV = use_refr_settings ? refr_roughnessV : refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughness_v ) );
					dr_roughnessU *= pow(rDepthMultiplier, (float) sg->Rr) + (rDepthAdder * (float) sg->Rr) + rOffset;
					dr_roughnessV *= pow(rDepthMultiplier, (float) sg->Rr) + (rDepthAdder * (float) sg->Rr) + rOffset;

					float dr_first_scale = RayState->media_refractDirect.v[m_cMatID];
					float dr_second_scale = AiShaderEvalParamFlt( p_dr_second_scale ) * dr_first_scale;
					dr_first_scale -= dr_second_scale;

					const bool two_lobes = dr_second_scale > ZERO_EPSILON;
					const float drs_rDepthMultiplier = two_lobes ? AiShaderEvalParamFlt( p_dr_second_roughnessMultiplier ) : 0.0f;
					const float drs_roughnessU = dr_roughnessU * drs_rDepthMultiplier;
					const float drs_roughnessV = dr_roughnessV * drs_rDepthMultiplier;

					void * btdf_data_direct = NULL;
					void * btdf_data_direct2 = NULL;
					iinfo.getDirectRefractionBTDFs(dr_btdf, &ppsg, dr_roughnessU, dr_roughnessV, drs_roughnessU, drs_roughnessV, pval_custom_tangent, &btdf_data_direct, &btdf_data_direct2);

					const bool refract_skydomes = AiShaderEvalParamBool(p_refract_skydomes);
					AiLightsPrepare(&ppsg);
					AiStateSetMsgBool("opaqueShadowMode", true);

					iinfo.directRefractionSampleLights(&ppsg, dr_btdf, refract_skydomes, two_lobes, btdf_data_direct, btdf_data_direct2, &acc_refract_direct, &acc_refract_direct_second);

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
				if ( traceSwitch.traceAnySpecular() || do_TIR ) 
				{
					float spec_roughnessU = iinfo.getSpecRoughnessU();
					float spec_roughnessV = iinfo.getSpecRoughnessV();

					SAMPLETYPE specular_sample[2];
					AtSamplerIterator* specularIterator = AiSamplerIterator( data->specular_sampler, sg);
					AtRay specularRay;
					AtUInt32 rayType;

					if (pval_specular_ray_type == 0)  
						rayType = AI_RAY_GLOSSY;
					else if (pval_specular_ray_type == 1)  
						rayType = AI_RAY_REFLECTED;
					if (do_TIR) 
						rayType = AI_RAY_REFRACTED;

					AiMakeRay(&specularRay, rayType, &sg->P, NULL, AI_BIG, sg);


					void * brdf_data; 
					int spec_brdf = RayState->media_BRDF.v[iinfo.m_higherPriority];
					bool sharp_reflection = (spec_roughnessU < ZERO_EPSILON && spec_roughnessV < ZERO_EPSILON);
					if (sharp_reflection)
						spec_brdf = b_stretched_phong;
					switch ( spec_brdf )
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
							blurAnisotropicPoles(&spec_roughnessU, &spec_roughnessV, &RayState->media_blurAnisotropicPoles.v[iinfo.m_higherPriority], &sg->Nf, &tangentSourceVector);
							AiV3Cross(uTangent, sg->N, sg->Rd); 
							AiV3Cross(vTangent, sg->N, uTangent);
							brdf_data = AiWardDuerMISCreateData(sg, &uTangent, &vTangent, spec_roughnessU, spec_roughnessV); 
							break;
						case b_ward_userTangent:
							// Ward with user tangents
							tangentSourceVector = AiV3Normalize( pval_custom_tangent );
							blurAnisotropicPoles(&spec_roughnessU, &spec_roughnessV, &RayState->media_blurAnisotropicPoles.v[iinfo.m_higherPriority], &sg->Nf, &tangentSourceVector);
							AiV3Cross( uTangent, sg->N, tangentSourceVector ) ;
							AiV3Cross( vTangent, sg->N, uTangent ) ;
							brdf_data = AiWardDuerMISCreateData( sg, &vTangent, &uTangent, spec_roughnessU, spec_roughnessV ) ; 
							break;
					}

					// ---------------------------------------------------//
					// Indirect Specular / Reflection
					// ---------------------------------------------------//

					// decision point- indirect specular
					// if ( trace_spec_indirect || do_TIR )
					if ( traceSwitch.spec_ind || do_TIR )
					{
						const float weight = fresnelTerm * RayState->media_specIndirect.v[iinfo.m1] * overallResultScale;
						const AtColor energyCache = RayState->ray_energy;						
						const AtColor energyCache_photon = RayState->ray_energy_photon;						
						const bool reflect_skies = AiShaderEvalParamBool(p_reflect_skies);

						if ( do_TIR )
						{
							RayState->ray_TIRDepth++ ;

							if (RayState->ray_TIRDepth < 50 && specularRay.refr_bounces > 1)
								specularRay.refr_bounces-- ;

							RayState->ray_energy *= TIR_color;
							RayState->ray_energy_photon *= TIR_color;
						} 
						else
						{
							RayState->ray_energy *= weight;		
							RayState->ray_energy_photon *= weight;		
						}


						if (sg->Rt == AI_RAY_CAMERA)
							RayState->ray_energy_photon /= (float) data->gloss_samples;

						AiStateSetMsgRGB("photon_energy",RayState->ray_energy_photon);

						while ( AiSamplerGetSample(specularIterator, specular_sample) )
						{
							if (sharp_reflection) 
							{
								specular_sample[0] = 0.5f;
								specular_sample[1] = 0.5f;
							}
							switch ( spec_brdf )
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
										acc_spec_indirect += sample.color * TIR_color * transmissionOnSample(&iinfo.t1, &sample, tracehit );
									else 
										acc_spec_indirect += sample.color * weight * transmissionOnSample(&iinfo.t1, &sample, tracehit );
								}
							}
							if (sharp_reflection)
								break;
						}
						if (!sharp_reflection)
							acc_spec_indirect *= (float) AiSamplerGetSampleInvCount(specularIterator);

						RayState->ray_energy = energyCache;						
						RayState->ray_energy_photon = energyCache_photon;						
					}

					// ---------------------------------------------------//
					// Direct Specular / Reflection
					// ---------------------------------------------------//

					// decision point- direct specular
					// if ( trace_spec_direct )
					if ( traceSwitch.spec_dir )
					{
						AtColor weight = AI_RGB_WHITE * fresnelTerm * RayState->media_specDirect.v[iinfo.m_higherPriority] * overallResultScale;
						if (do_TIR)
							weight *= TIR_color;
						
						AiStateSetMsgBool("opaqueShadowMode", true);						
						AiLightsPrepare(sg); 
						const bool reflect_skydomes = AiShaderEvalParamBool(p_reflect_skydomes);
						while ( AiLightsGetSample(sg) ) // loop over the lights to compute direct effects
						{
							float l_weight = AiLightGetSpecular(sg->Lp);

							if ( 
								(reflect_skydomes || !AiNodeIs( sg->Lp,"skydome_light" )) && l_weight > ZERO_EPSILON
								)
							{
								switch ( RayState->media_BRDF.v[iinfo.m_higherPriority] )
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

	AtRGBA validInterfaceResult = AI_RGBA_BLACK;

	if ( !iinfo.validInterface )
	{
		if ( RayState->ray_invalidDepth < 70 )
		{
			AtRay ray;
			AtScrSample sample;

			RayState->ray_invalidDepth ++;
			updateMediaInsideLists(m_cMatID, iinfo.entering, media_inside_ptr, false);

			AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, NULL, AI_BIG, sg);
			ray.dir = sg->Rd;
			ray.level --;
			ray.refr_bounces --;
			const bool tracehit = AiTrace(&ray, &sample);

			AiRGBtoRGBA( sample.color * transmissionOnSample(&iinfo.t2, &sample, tracehit ), validInterfaceResult );
			validInterfaceResult.a = sample.alpha;			
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
	if (iinfo.validInterface)
	{
		AiStateSetMsgBool("opaqueShadowMode", true);
		emission += AiShaderEvalParamRGB(p_emission_at_interfaces);
		if (iinfo.mediaExit)
		{
			emission += AiShaderEvalParamRGB(p_emission_at_exits);
		}
		else if (iinfo.mediaEntrance)
		{
			emission += AiShaderEvalParamRGB(p_emission_at_entrances);
		}
		AiStateSetMsgBool("opaqueShadowMode", false);
	}

	// ---------------------------------------------------//
	// - Output and AOVs
	// ---------------------------------------------------//	

	if ( iinfo.validInterface )
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
		sg->out.RGBA = validInterfaceResult;
	}

	// ---------------------------------------------------//
	// - Array and messaging post process 
	// Reset all messages to the previous state. If there was no previous state, Mark them invalid. 
	// ---------------------------------------------------//

	AiStateSetMsgRGB("photon_energy",AI_RGB_BLACK); // this doesn't seem right, heberlein caustics related only

	if (msgs_are_valid)
	{
		AiStateSetMsgBool("msgs_are_valid", true);
		RayState->uncacheRayState(&RayStateCache);
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
