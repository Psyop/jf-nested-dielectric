/* 
 * JF Nested Dielectric by Jonah Friedman
 * 1.0.4
 * Copyright (c) 2014, Psyop Media Company, LLC and Jonah Friedman 
 * Open sourced under the 3 clause BSD license, see license.txt
 */

#define NOMINMAX // lets you keep using std::min

#include <ai.h>
#include <cstring>
#include <string>
#include <math.h>
#include "spectral_distributions.h"
#include "jf_nested_dielectric.h"

AI_SHADER_NODE_EXPORT_METHODS(jf_nested_dielectric_methods);

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
    Ray_State *rayState;

    // ---------------------------------------------------//
    // - Array and messaging preprocess 
    // ---------------------------------------------------//

    bool msgs_are_valid = false;
    AiStateGetMsgBool(JFND_MSG_VALID_BOOL, &msgs_are_valid);

    if ( msgs_are_valid )
    {
        void * rayState_ptr;
        AiStateGetMsgPtr(JFND_MSG_RAYSTATE_PTR, &rayState_ptr);
        rayState = static_cast<Ray_State*>( rayState_ptr );
        rayState->cacheRayState();
    }
    else
    {       
        rayState = static_cast<Ray_State*>( AiShaderGlobalsQuickAlloc(sg, sizeof( Ray_State ) ) );
        rayState->init(data, sg, node, data->polarizationVector);
        AiStateSetMsgPtr(JFND_MSG_RAYSTATE_PTR, rayState);
    }
    AiStateSetMsgBool(JFND_MSG_VALID_BOOL, true); // Any child rays from this will find valid messages. 

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

        media_inside_ptr = &rayState->shadow_media_inside;

        bool shadowlist_is_valid = false;
        AiStateGetMsgBool(JFND_MSG_SHADOW_VALID_BOOL, &shadowlist_is_valid);


        bool transp_index_reset = false;
        int prev_transp_index;
        if (AiStateGetMsgInt(JFND_MSG_PREV_TRANSP_INT, &prev_transp_index))
        {
            if ((int) sg->transp_index <= prev_transp_index)
            {
                transp_index_reset = true;  
            }
        }

        if ( !shadowlist_is_valid || transp_index_reset )  // if there has been another kind of ray, or the transp_index did not increase
        {
            // intialize the shadow media inside list
            memcpy(&rayState->shadow_media_inside, &rayState->media_inside, sizeof(MediaIntStruct) );
        }

        AiStateSetMsgInt(JFND_MSG_PREV_TRANSP_INT, sg->transp_index);
        AiStateSetMsgBool(JFND_MSG_SHADOW_VALID_BOOL, true);
    }
    else
    {
        media_inside_ptr = &rayState->media_inside;
        AiStateSetMsgBool(JFND_MSG_SHADOW_VALID_BOOL, false);
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

        const char *overlapping_surfaces_message = "JF Nested Dielectric: Crazy values in media lists, you may have some perfectly overlapping surfaces.";
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

        const char *overlapping_surfaces_message = "JF Nested Dielectric: Crazy positive values in media lists, you may have some perfectly overlapping surfaces.";
        AiMsgWarning(overlapping_surfaces_message);
        sg->out.RGBA = AI_RGBA_BLACK;
        return; 
    }

    if (m_cMatID >= max_media_count || m_cMatID < 0)
    {
        const char *priority_error_message = "JF Nested Dielectric: Medium priority must be between 0 and 32!";
        AiMsgError(priority_error_message);
        return;
    }

    // ---------------------------------------------------//
    // - get interface info     
    // ---------------------------------------------------//

    rayState->readBasicMatParameters(sg, node, m_cMatID);
    InterfaceInfo iinfo = InterfaceInfo(rayState, m_cMatID, sg);

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
        rayState->ray_energy *= cTransmission;
        rayState->ray_energy_photon *= cTransmission;
    }

    // ---------------------------------------------------//
    // param caching
    //     caching select parameters if it improves the code to do so
    // ---------------------------------------------------//
    const AtVector pval_custom_tangent = AiShaderEvalParamVec(p_ward_tangent);
    const int pval_specular_ray_type = AiShaderEvalParamInt(p_specular_ray_type);

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
        TraceSwitch traceSwitch = TraceSwitch(&iinfo, AiShaderEvalParamBool(p_enable_internal_reflections));
        float overallResultScale = 1.0f;
        const bool isPhoton = rayIsPhoton(sg);
        const bool causticPath = isPhoton || sg->Rr_diff > 0 ;
        if ( causticPath ) 
        {
            if (iinfo.mediaEntrance || !rayState->caustic_behaviorSet)
            {
                rayState->readCausticMatParameters(sg, node);
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


            if ((sg->Rr_diff > 0 && rayState->caustic_mode == 0)) //means we're in the right kind of caustics for the right kind of ray.
            {
                do_disperse = do_disperse && rayState->caustic_dispersion;
                traceSwitch.setPathtracedCaustic(&iinfo);
            }
            else if (isPhoton)
            {
                do_disperse = do_disperse && rayState->caustic_dispersion;
                traceSwitch.setPhotonCaustic(&iinfo);
            }
            else
            {
                traceSwitch.setTraceNone();
            }
        }       

        bool energySignificant;

        if (   std::abs( rayState->ray_energy.r ) < rayState->energy_cutoff
            && std::abs( rayState->ray_energy.g ) < rayState->energy_cutoff
            && std::abs( rayState->ray_energy.b ) < rayState->energy_cutoff
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
        float dispersion_sample[2], refraction_sample[2], specular_sample[2];

        // decision point- trace anything
        if ( energySignificant && traceSwitch.traceAnything() ) 
        {
            AtRay ray;
            AtScrSample sample;
            
            float fresnelTerm = fresnelEquations (iinfo.n1, iinfo.n2,  AiV3Dot(sg->Nf, -sg->Rd), iinfo.polarizationTerm, true);

            bool do_TIR = false;
            AtColor TIR_color = AI_RGB_BLACK;

            if ( traceSwitch.traceAnyRefraction() )
            {
                int refractSamplesTaken = 0;

                // ---------------------------------------------------//
                // Refraction
                // Samplers
                // ---------------------------------------------------//

                AtSamplerIterator* dispersionIterator = AiSamplerIterator( data->dispersion_sampler, sg) ;
                AtSamplerIterator* refractionIterator = AiSamplerIterator( data->refraction_sampler, sg);

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

                AtVector sgrd_cache = sg->Rd;
                AtShaderGlobals ppsg = *sg;
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
                    const bool refract_skies = AiShaderEvalParamBool(p_refract_skies);

                    while ( AiSamplerGetSample(refractionIterator, refraction_sample) )
                    {
                        if (!refracted && !do_disperse)
                            continue;

                        AtColor monochromeColor = AI_RGB_WHITE;
                        bool dispersion_sample_TIR = false;
                        if (do_disperse)
                        {
                            AiSamplerGetSample(dispersionIterator, dispersion_sample);
                            if (dispersal_seed < 0.0f)
                            {   // The job of a dispersal seed is to fix any correlations, but still allow stratefied sampling to work in batches of samples.
                                dispersal_seed =  ( std::abs( sg->sx + sg->sy ) * 113 + (float) dispersion_sample[1] ) * 3.456f  ;
                            }
                            float n1_disp, n2_disp;
                            iinfo.disperse((float) (dispersal_seed + dispersion_sample[0]), &n1_disp, &n2_disp, &monochromeColor);

                            AiMakeRay(&dispersalRay, AI_RAY_REFRACTED, &sg->P, &sg->Rd, AI_BIG, sg);                                
                            if (!AiRefractRay(&dispersalRay, &sg->Nf, n1_disp, n2_disp, sg))
                            {   // TIR
                                dispersion_sample_TIR = true;
                                TIR_color += monochromeColor;
                                dispersed_TIR_samples++;
                                refractSamplesTaken++ ;
                            }
                            ray.dir = dispersalRay.dir;
                        }


                        if ( do_blurryRefraction ) {
                            if ( do_disperse )
                            {   // redo ppsg creation
                                ppsg.Rd = sgrd_cache;
                                parallelPark(ray.dir, &ppsg);
                                btdf_data = iinfo.getBTDFData(&ppsg, refr_roughnessU, refr_roughnessV, pval_custom_tangent);
                            }
                            ray.dir = iinfo.getBTDFSample(btdf_data, refraction_sample);
                        }

                        if ((AiV3Dot(ray.dir,sg->Nf) < 0.0f) && (ray.dir != AI_V3_ZERO) && !dispersion_sample_TIR)
                        {
                            updateMediaInsideLists(m_cMatID, iinfo.entering, media_inside_ptr, false);
                            const AtColor weight = (1.0f - fresnelTerm) * monochromeColor 
                                * overallResultScale * iinfo.getIndirectRefractionWeight();

                            const AtColor energyCache = rayState->updateEnergyReturnOrig(weight);
                            const AtColor energyCache_photon = rayState->updatePhotonEnergyReturnOrig(weight);

                            if (sg->Rt == AI_RAY_CAMERA && (do_disperse || do_blurryRefraction))                                
                                rayState->ray_energy_photon /= (float) data->refr_samples;

                            AiStateSetMsgRGB("photon_energy",rayState->ray_energy_photon);

                            refractSamplesTaken++ ;
                            const bool tracehit = AiTrace(&ray, &sample);
                            if (tracehit || refract_skies) 
                                acc_refract_indirect += sample.color * weight * transmissionOnSample(&iinfo.t2, &sample, tracehit );

                            rayState->resetEnergyCache(energyCache);
                            rayState->resetPhotonEnergyCache(energyCache_photon);

                            updateMediaInsideLists(m_cMatID, iinfo.entering, media_inside_ptr, true);
                        }

                        if (!do_multiSampleRefraction)
                            break;
                    }
                    
                    /*
                     * Refraction - TIR
                     * Total internal reflection sampling
                     * TIR counts as refraction rays, and depth is counted seperately by rayState->ray_TIRDepth
                     * Actually occurs in the specular loop below though. 
                     */

                    // decision point - tir
                    if (do_disperse)
                    {
                        if (dispersed_TIR_samples > 0)
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
                    int dr_btdf         = use_refr_settings ? iinfo.getBTDFType() : AiShaderEvalParamEnum( p_dr_btdf );
                    float dr_roughnessU = use_refr_settings ? refr_roughnessU : refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughness_u ) );
                    float dr_roughnessV = use_refr_settings ? refr_roughnessV : refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughness_v ) );
                    dr_roughnessU *= pow(rDepthMultiplier, (float) sg->Rr) + (rDepthAdder * (float) sg->Rr) + rOffset;
                    dr_roughnessV *= pow(rDepthMultiplier, (float) sg->Rr) + (rDepthAdder * (float) sg->Rr) + rOffset;

                    float dr_first_scale = rayState->media_refractDirect.v[m_cMatID];
                    float dr_second_scale = AiShaderEvalParamFlt( p_dr_second_scale ) * dr_first_scale;
                    dr_first_scale -= dr_second_scale;

                    const bool two_lobes = dr_second_scale > ZERO_EPSILON;
                    const float drs_rDepthMultiplier = two_lobes ? AiShaderEvalParamFlt( p_dr_second_roughnessMultiplier ) : 0.0f;
                    const float drs_roughnessU = dr_roughnessU * drs_rDepthMultiplier;
                    const float drs_roughnessV = dr_roughnessV * drs_rDepthMultiplier;

                    void * btdf_data_direct = NULL;
                    void * btdf_data_direct2 = NULL;
                    iinfo.getDirectRefractionBTDFs(dr_btdf, &ppsg, dr_roughnessU, dr_roughnessV, 
                        drs_roughnessU, drs_roughnessV, pval_custom_tangent, &btdf_data_direct, &btdf_data_direct2);

                    const bool refract_skydomes = AiShaderEvalParamBool(p_refract_skydomes);
                    AiLightsPrepare(&ppsg);
                    AiStateSetMsgBool("opaqueShadowMode", true);

                    iinfo.directRefractionSampleLights(&ppsg, dr_btdf, refract_skydomes, two_lobes, 
                        btdf_data_direct, btdf_data_direct2, &acc_refract_direct, &acc_refract_direct_second);

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
                    rayState->ray_monochromatic = false;
                }

                // decision point- any specular
                if ( traceSwitch.traceAnySpecular() || do_TIR ) 
                {
                    float spec_roughnessU = iinfo.getSpecRoughnessU();
                    float spec_roughnessV = iinfo.getSpecRoughnessV();
                    bool sharp_reflection = (spec_roughnessU < ZERO_EPSILON && spec_roughnessV < ZERO_EPSILON);

                    AtSamplerIterator* specularIterator = AiSamplerIterator( data->specular_sampler, sg);
                    AtRay specularRay;

                    AtUInt32 rayType;
                    if (do_TIR) 
                        rayType = AI_RAY_REFRACTED;
                    else if (pval_specular_ray_type == 0)  
                        rayType = AI_RAY_GLOSSY;
                    else if (pval_specular_ray_type == 1)  
                        rayType = AI_RAY_REFLECTED;

                    AiMakeRay(&specularRay, rayType, &sg->P, NULL, AI_BIG, sg);

                    int spec_brdf = rayState->media_BRDF.v[iinfo.m_higherPriority];
                    if (sharp_reflection)
                        spec_brdf = b_cook_torrance;

                    AtVector tangentSourceVector, uTangent, vTangent;
                    switch ( spec_brdf )
                    {
                        case b_cook_torrance:
                            // Cook Torrance
                            uTangent = AI_V3_ZERO;
                            vTangent = AI_V3_ZERO;
                            break;
                        case b_cook_torrance_ray_tangent:
                            // Ward with refraction-derivitive tangents
                             // to do: tangentSourceVector is uninitialized? Should this be the eye ray or something?
                            blurAnisotropicPoles(&spec_roughnessU, &spec_roughnessV, &rayState->media_blurAnisotropicPoles.v[iinfo.m_higherPriority], &sg->Nf, &tangentSourceVector);
                            AiV3Cross(uTangent, sg->N, sg->Rd); 
                            AiV3Cross(vTangent, sg->N, uTangent);
                            break;
                        case b_cook_torrance_user_tangent:
                            // Ward with user tangents
                            tangentSourceVector = AiV3Normalize( pval_custom_tangent );
                            blurAnisotropicPoles(&spec_roughnessU, &spec_roughnessV, &rayState->media_blurAnisotropicPoles.v[iinfo.m_higherPriority], &sg->Nf, &tangentSourceVector);
                            AiV3Cross( uTangent, sg->N, tangentSourceVector ) ;
                            AiV3Cross( vTangent, sg->N, uTangent ) ;
                            break;
                    }
                    void * brdf_data = AiCookTorranceMISCreateData(sg, &AI_V3_ZERO, &AI_V3_ZERO, spec_roughnessU, spec_roughnessU);

                    // ---------------------------------------------------//
                    // Indirect Specular / Reflection
                    // ---------------------------------------------------//

                    // decision point- indirect specular
                    // if ( trace_spec_indirect || do_TIR )
                    if ( traceSwitch.spec_ind || do_TIR )
                    {
                        const float weight = fresnelTerm * rayState->media_specIndirect.v[iinfo.m1] * overallResultScale;
                        const bool reflect_skies = AiShaderEvalParamBool(p_reflect_skies);
                        AtColor energyCache;
                        AtColor energyCache_photon; 

                        if ( do_TIR )
                        {
                            rayState->ray_TIRDepth++;
                            if (rayState->ray_TIRDepth < 50 && specularRay.refr_bounces > 1)
                                specularRay.refr_bounces-- ;
                            energyCache = rayState->updateEnergyReturnOrig(TIR_color);
                            energyCache_photon = rayState->updatePhotonEnergyReturnOrig(TIR_color);
                        } 
                        else
                        {
                            energyCache = rayState->updateEnergyReturnOrig(weight * AI_RGB_WHITE);
                            energyCache_photon = rayState->updatePhotonEnergyReturnOrig(weight * AI_RGB_WHITE); 
                        }


                        if (sg->Rt == AI_RAY_CAMERA)
                            rayState->ray_energy_photon /= (float) data->gloss_samples;

                        AiStateSetMsgRGB("photon_energy",rayState->ray_energy_photon);

                        while ( AiSamplerGetSample(specularIterator, specular_sample) )
                        {
                            if (sharp_reflection) 
                            {
                                specular_sample[0] = 0.5f;
                                specular_sample[1] = 0.5f;
                            }

                            specularRay.dir = AiCookTorranceMISSample(brdf_data, (float) specular_sample[0], (float) specular_sample[1]);
                                
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

                        rayState->resetEnergyCache(energyCache);
                        rayState->resetPhotonEnergyCache(energyCache_photon);               
                    }

                    // ---------------------------------------------------//
                    // Direct Specular / Reflection
                    // ---------------------------------------------------//

                    // decision point- direct specular
                    // if ( trace_spec_direct )
                    if ( traceSwitch.spec_dir )
                    {
                        AtColor weight = AI_RGB_WHITE * fresnelTerm * rayState->media_specDirect.v[iinfo.m_higherPriority] * overallResultScale;
                        if (do_TIR)
                            weight *= TIR_color;
                        
                        AiStateSetMsgBool("opaqueShadowMode", true);                        
                        AiLightsPrepare(sg); 
                        const bool reflect_skydomes = AiShaderEvalParamBool(p_reflect_skydomes);
                        while ( AiLightsGetSample(sg) ) // loop over the lights to compute direct effects
                        {
                            float l_weight = AiLightGetSpecular(sg->Lp);
                            if ( (reflect_skydomes || !AiNodeIs( sg->Lp,"skydome_light" )) && l_weight > ZERO_EPSILON)
                            {
                                acc_spec_direct += l_weight * AiEvaluateLightSample(sg, brdf_data, AiCookTorranceMISSample, AiCookTorranceMISBRDF, AiCookTorranceMISPDF);                         
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
        if ( rayState->ray_invalidDepth < 70 )
        {
            AtRay ray;
            AtScrSample sample;

            rayState->ray_invalidDepth ++;
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
            const char *overlapping_surfaces_message = "JF Nested Dielectric: Crazy numbers of invalid interfaces. Some geo may be overlapping.";
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
        AiStateSetMsgBool(JFND_MSG_VALID_BOOL, true);
        rayState->uncacheRayState();
    }
    else
    {
        AiStateSetMsgBool(JFND_MSG_VALID_BOOL, false);
    }
}

node_loader
{
    if (i > 0)
        return false;

    node->methods          = jf_nested_dielectric_methods;
    node->output_type  = AI_TYPE_RGBA;
    node->name                    = "jf_nested_dielectric";
    node->node_type  = AI_NODE_SHADER;
    strcpy(node->version, AI_VERSION);
    return true;
}
