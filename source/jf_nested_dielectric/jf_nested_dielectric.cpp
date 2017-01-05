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
    AiParameterFLT("dr_roughnessDepthMultiplier", 0.0f); // to do: rename to addOnReflect
    AiParameterFLT("dr_roughnessDepthAdder", 0.0f); // to do: rename to addOnRefract
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
    Ray_State_Cache rayStateCache;

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
        rayState->cacheRayState( &rayStateCache);
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
            rayState->lastShadowDepth = 0.0f;
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
     * media_id is the medium ID of the material we're currently evaluating, regardless of
     * whether or not the interface is valid or what media we're inside.
     *
     * The medium IDs all get 1 added to them. This means 0 is reserved for "no medium". This simplifies other code, but it means
     * that if the user gives a medium an ID of 0, internally in the shader it has an ID of 1. 
     */

    const int media_id = AiShaderEvalParamInt(p_mediumPriority) + 1;

    if (media_inside_ptr->v[media_id] < -10 )
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
    if ( media_inside_ptr->v[media_id] > 30)
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

    if (media_id >= max_media_count || media_id < 0)
    {
        const char *priority_error_message = "JF Nested Dielectric: Medium priority must be between 0 and 32!";
        AiMsgError(priority_error_message);
        return;
    }

    // ---------------------------------------------------//
    // - get interface info     
    // ---------------------------------------------------//

    rayState->readBasicMatParameters(sg, node, media_id);
    InterfaceInfo iinfo = InterfaceInfo(rayState, media_id, sg);

    const bool do_blurryRefraction = iinfo.doBlurryRefraction();
    bool do_disperse = iinfo.setupDispersion(data);
    // const bool do_multiSampleRefraction = do_disperse || do_blurryRefraction;

    // ---------------------------------------------------//
    // - Shadow rays
    // ---------------------------------------------------//
    
    if (sg->Rt == AI_RAY_SHADOW)
    {
        float depth_portion = (float) sg->Rl - rayState->lastShadowDepth;
        rayState->lastShadowDepth = (float) sg->Rl;
        AtColor transparency = iinfo.getShadowTransparency(sg, depth_portion, AiShaderEvalParamEnum(p_shadow_mode));
        sg->out_opacity = AI_RGB_WHITE - transparency;
        if (sg->out_opacity != AI_RGB_WHITE)
            updateMediaInsideLists(media_id, iinfo.entering, media_inside_ptr, false);
        return;
    } else {
        AtColor cTransmission = iinfo.getTransmissionColor(sg, (float) sg->Rl);
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
        const bool causticPath = sg->Rr_diff > 0 ;
        if ( causticPath ) 
        {
            // to do: does not currently apply to photon caustics, only applies to path traced
            // means no cuastics settings in JFND apply. 
            // need to know we're making photons to make this work properly. 
            if (iinfo.mediaEntrance || !rayState->caustic_behaviorSet)
                rayState->readCausticMatParameters(sg, node);

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
            else
            {
                traceSwitch.setTraceNone();
            }
        }

        // ---------------------------------------------------//
        // Main Ray Tracing
        // Valid interfaces
        // ---------------------------------------------------//


        float dispersion_sample[2], refraction_sample[2], specular_sample[2];
        // decision point- trace anything
        if (traceSwitch.traceAnything()) 
        {
            AtRay ray;
            AtScrSample sample;
            
            float fresnelTerm = fresnelEquations(iinfo.n1, iinfo.n2,  AiV3Dot(sg->Nf, -sg->Rd), 
                iinfo.polarizationTerm, true);

            const bool refrStronger = fresnelTerm < 0.5;
            const bool optoEnergySignificant = std::abs( rayState->ray_energy.r ) > rayState->energy_cutoff 
                                        || std::abs( rayState->ray_energy.g ) > rayState->energy_cutoff 
                                        || std::abs( rayState->ray_energy.b ) > rayState->energy_cutoff;
            bool optoAllowBranch = optoEnergySignificant;

            const float rrProbability = AiShaderEvalParamFlt(p_russian_roulette_probability);
            float rrRefrStrength = 1.0;
            float rrSpecStrength = 1.0;
            if ( optoAllowBranch && rrProbability > ZERO_EPSILON && sg->Rr_refr > 0 )
            {
                // const float effectiveProbability = (AiColorToGrey(rayState->ray_energy) + 0.01f) * rrProbability * 10.0f;
                const float energyProb = (AiColorToGrey(rayState->ray_energy) + 0.01f) * 4.0;
                const float averageProb = (energyProb + rrProbability) / 2.0f;

                AtSamplerIterator* rr_iterator = AiSamplerIterator( data->russian_roullete_single_sampler, sg );
                const float rrResult = russianRoulette( rr_iterator, 1.0, averageProb);
                float dummy[2];
                while ( AiSamplerGetSample(rr_iterator, dummy) ) {} // exhaust the sampler. 

                if (rrResult == 0.0)
                    optoAllowBranch = false;
                else if (refrStronger) 
                    rrSpecStrength *= rrResult; // refr is stronger, spec was on the chopping block. 
                else
                    rrRefrStrength *= rrResult; // spec is stronger, refr was on the chopping block. 
            }

            const bool optoAllowRefr = optoAllowBranch || refrStronger;
            bool optoAllowRefl = optoAllowBranch || !refrStronger;

            bool do_TIR = false;
            AtColor TIR_color = AI_RGB_BLACK;

            if ( traceSwitch.traceAnyRefraction() && optoAllowRefr )
            {
                // ---------------------------------------------------//
                // Refraction
                // Samplers
                // ---------------------------------------------------//

                AtSamplerIterator* dispersionIt = AiSamplerIterator(data->dispersion_sampler, sg) ;
                AtSamplerIterator* refractionIt = AiSamplerIterator(data->refraction_sampler, sg);

                // ---------------------------------------------------//
                // Refraction
                // BTDF preprocessing
                // ---------------------------------------------------//

                float refr_roughnessU, refr_roughnessV;
                iinfo.getRefrRoughness(refr_roughnessU, refr_roughnessV);

                if (causticPath)
                {
                    const float causticDRRoughnessOffset = refractRoughnessConvert(AiShaderEvalParamFlt(p_caustic_dr_roughness_offset));
                    refr_roughnessU += causticDRRoughnessOffset;
                    refr_roughnessV += causticDRRoughnessOffset;
                }

                AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, &sg->Rd, AI_BIG, sg);                 
                const bool tir = !AiRefractRay(&ray, &sg->Nf, iinfo.n1, iinfo.n2, sg);
                if (tir)
                    optoAllowRefl = true;

                const AtVector cache_N   = sg->N;
                const AtVector cache_Nf  = sg->Nf;
                const AtVector cache_Ngf = sg->Ngf;
                const AtVector cache_Rd  = sg->Rd;
                sg->Nf *= -1; // flip the forward facing normals
                sg->Ngf *= -1; // flip the forward facing normals
                sg->Rd = parallelPark(ray.dir, cache_N); 

                // decision point- indirect refraction
                if ( traceSwitch.refr_ind )
                {

                    void * btdf_data = NULL;
                    if (!do_disperse && do_blurryRefraction)
                    {
                        // we're not dispersing, and we are doing blurry refraction. 
                        // No per sample btdf needed, we just make one here.  
                        btdf_data = iinfo.getRefrBRDFData(sg, refr_roughnessU, refr_roughnessV, pval_custom_tangent);
                    }
                    // ---------------------------------------------------//
                    // Refraction - Indirect
                    // Indirect REfraction Sampling
                    // ---------------------------------------------------//
                    


                    AtRay dispersalRay;
                    float dispersal_seed = -1.0f;
                    int dispersed_TIR_samples = 0;
                    int refractSamplesTaken = 0;
                    const bool refract_skies = AiShaderEvalParamBool(p_refract_skies);
                    const bool skip_sampling = tir && !do_disperse;

                    while ( AiSamplerGetSample(refractionIt, refraction_sample) )
                    {
                        if (skip_sampling)
                            continue;
                        
                        AtColor monochromeColor = AI_RGB_WHITE;
                        bool dispersion_sample_TIR = false;

                        if (do_disperse)
                        {
                            AiSamplerGetSample(dispersionIt, dispersion_sample);
                            if (dispersal_seed < 0.0f)
                            {   // The job of a dispersal seed is to fix any correlations, 
                                // but still allow stratefied sampling to work in batches of samples.
                                dispersal_seed =  ( std::abs( sg->sx + sg->sy ) * 113 + (float) dispersion_sample[1] ) * 3.456f  ;
                            }
                            float n1_disp, n2_disp;
                            iinfo.disperse(dispersal_seed + dispersion_sample[0], &n1_disp, &n2_disp, &monochromeColor);

                            sg->Rd = cache_Rd;
                            AiMakeRay(&dispersalRay, AI_RAY_REFRACTED, &sg->P, &cache_Rd, AI_BIG, sg);
                            dispersion_sample_TIR = !AiRefractRay(&dispersalRay, &cache_N, n1_disp, n2_disp, sg);
                            ray.dir = dispersalRay.dir; // note: must happen after AiRefractRay
                            if (dispersion_sample_TIR)
                            {   // TIR
                                TIR_color += monochromeColor;
                                dispersed_TIR_samples++;
                                refractSamplesTaken++ ;
                            }

                            if ( do_blurryRefraction ) 
                            {
                                sg->Rd = parallelPark(ray.dir, cache_N); 
                                btdf_data = iinfo.getRefrBRDFData(sg, refr_roughnessU, refr_roughnessV, pval_custom_tangent);
                            }
                        } 

                        if ( do_blurryRefraction )
                            ray.dir = AiCookTorranceMISSample(btdf_data, (float) refraction_sample[0], (float) refraction_sample[1]);

                        if (AiV3Dot(ray.dir, cache_Nf) < 0 && !dispersion_sample_TIR)
                        {
                            updateMediaInsideLists(media_id, iinfo.entering, media_inside_ptr, false);
                            const AtColor weight = (1.0f - fresnelTerm) * monochromeColor 
                                * overallResultScale * iinfo.getIndirectRefractionWeight() * rrRefrStrength;

                            const AtColor cache_energy = rayState->updateEnergyReturnOrig(weight);
                            const AtColor cache_photonEnergy = rayState->updatePhotonEnergyReturnOrig(weight);
                            const float drDepthAdderRefr = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughnessDepthAdder ) );
                            const float drAccumRoughness = rayState->updateDRAccumRoughnessReturnOrig(drDepthAdderRefr);

                            if (sg->Rt == AI_RAY_CAMERA && (do_disperse || do_blurryRefraction)) // to do: only needs to happen for photon rays
                                rayState->ray_energy_photon /= (float) data->refr_samples;

                            AiStateSetMsgRGB(JFND_MSG_PHOTON_RGB,rayState->ray_energy_photon);

                            refractSamplesTaken++ ;
                            const bool tracehit = AiTrace(&ray, &sample);
                            if (tracehit || refract_skies) 
                                acc_refract_indirect += sample.color * weight * transmissionOnSample(&iinfo.t2, &sample, tracehit );

                            rayState->resetEnergyCache(cache_energy);
                            rayState->resetPhotonEnergyCache(cache_photonEnergy);
                            rayState->resetDRAccumRoughness(drAccumRoughness);

                            updateMediaInsideLists(media_id, iinfo.entering, media_inside_ptr, true);
                        }

                        if (!do_disperse && !do_blurryRefraction)
                            break; // one sample only
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
                            TIR_color *= (float) AiSamplerGetSampleInvCount(refractionIt);
                            do_TIR = true;
                        }
                    }
                    else
                    {                       
                        if (traceSwitch.TIR && tir)
                        {
                            TIR_color = AI_RGB_WHITE;
                            do_TIR = true;
                        }
                    }

                    if (refractSamplesTaken > 0) 
                        acc_refract_indirect /= (float) refractSamplesTaken;
                }

                // Exhaust the sampler. Good night, sampler. This seems necessary, having the unexhausted sampler caused some problems with the RR sampler. 
                while ( AiSamplerGetSample(dispersionIt, dispersion_sample) ){} // to do: not necessary when not dispersed?

                // ---------------------------------------------------//
                // Refraction - Direct
                // Direct Refraction Sampling
                // ---------------------------------------------------//

                // decision point- direct refraction
                if (traceSwitch.refr_dir && !tir)
                {
                    // offsets and depth modification
                    const float rFlatOffset = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughnessOffset ) );
                    const float rOffset = rFlatOffset + rayState->dr_accumRoughness; // to do: + the ray state depth value

                    const bool use_refr_settings = AiShaderEvalParamBool( p_dr_use_refraction_btdf );
                    int dr_btdf         = use_refr_settings ? iinfo.getRefrBRDFType() : AiShaderEvalParamEnum( p_dr_btdf );
                    float dr_roughnessU = use_refr_settings ? 
                        refr_roughnessU + rOffset : 
                        refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughness_u ) ) + rOffset;

                    float dr_roughnessV = dr_roughnessU;
                    if (dr_btdf != b_cook_torrance) 
                        dr_roughnessV = use_refr_settings ? 
                            refr_roughnessV  + rOffset : 
                            refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughness_v ) ) + rOffset;

                    // float dr_first_scale = rayState->media_refractDirect.v[media_id];
                    float dr_first_scale = rayState->media[media_id].refractDirect;
                    float dr_second_scale = AiShaderEvalParamFlt( p_dr_second_scale ) * dr_first_scale;
                    dr_first_scale -= dr_second_scale;

                    const bool two_lobes = dr_second_scale > ZERO_EPSILON;
                    const float drs_rDepthMultiplier = two_lobes ? AiShaderEvalParamFlt( p_dr_second_roughnessMultiplier ) : 0.0f;
                    const float drs_roughnessU = dr_roughnessU * drs_rDepthMultiplier;
                    const float drs_roughnessV = dr_roughnessV * drs_rDepthMultiplier;

                    void * btdf_data_direct = NULL;
                    void * btdf_data_direct2 = NULL;
                    iinfo.getDirectRefractionBTDFs(dr_btdf, sg, dr_roughnessU, dr_roughnessV, 
                        drs_roughnessU, drs_roughnessV, pval_custom_tangent, &btdf_data_direct, 
                        &btdf_data_direct2);

                    const bool refract_skydomes = AiShaderEvalParamBool(p_refract_skydomes);
                    AiLightsPrepare(sg);
                    AiStateSetMsgBool(JFND_MSG_OPQ_SHADOW_MODE_BOOL, true);

                    iinfo.directRefractionSampleLights(sg, dr_btdf, refract_skydomes, two_lobes, 
                        btdf_data_direct, btdf_data_direct2, &acc_refract_direct, &acc_refract_direct_second);

                    acc_refract_direct *= dr_first_scale * (1.0f - fresnelTerm) * overallResultScale * rrRefrStrength;
                    acc_refract_direct_second *= dr_second_scale * (1.0f - fresnelTerm) * overallResultScale * rrRefrStrength;
                }

                sg->N = cache_N;   
                sg->Nf = cache_Nf;   
                sg->Ngf = cache_Ngf;   
                sg->Rd = cache_Rd;
            }

            // ---------------------------------------------------//
            // Specular / Reflection / TIR
            // ---------------------------------------------------//

            // If we're doing TIR, either refraction was stronger (and so will == 1.0), or no 
            // russian roullete was done (and strength will also be 1.0. No RR strength needed on 
            // TIR. )
            if (do_TIR)
                fresnelTerm = 1.0f; 

            if (fresnelTerm > ZERO_EPSILON && optoAllowRefl)
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
                    float spec_roughnessU, spec_roughnessV;
                    iinfo.getSpecRoughness(spec_roughnessU, spec_roughnessV);
                    const bool sharp_reflection = (spec_roughnessU < ZERO_EPSILON && spec_roughnessV < ZERO_EPSILON);
                    void *brdf_data = iinfo.getSpecBRDFData(sg, spec_roughnessU, spec_roughnessV, pval_custom_tangent);

                    // ---------------------------------------------------//
                    // Indirect Specular / Reflection
                    // ---------------------------------------------------//

                    // decision point- indirect specular
                    if ( traceSwitch.spec_ind || do_TIR )
                    {
                        // const float weight = fresnelTerm * rayState->media_specIndirect.v[iinfo.m1] * overallResultScale * rrSpecStrength;
                        const float weight = fresnelTerm * rayState->media[iinfo.m1].specIndirect * overallResultScale * rrSpecStrength;
                        const bool reflect_skies = AiShaderEvalParamBool(p_reflect_skies);
                        AtSamplerIterator* specularIterator = AiSamplerIterator( data->specular_sampler, sg);

                        AtUInt32 rt = pval_specular_ray_type == 0 ? AI_RAY_GLOSSY : AI_RAY_REFLECTED;
                        if (do_TIR) 
                            rt = AI_RAY_REFRACTED;

                        AtRay specularRay;
                        AiMakeRay(&specularRay, rt, &sg->P, NULL, AI_BIG, sg);

                        AtColor cache_energy = rayState->updateEnergyReturnOrig(do_TIR ? TIR_color : weight * AI_RGB_WHITE); 
                        AtColor cache_photonEnergy = rayState->updatePhotonEnergyReturnOrig(do_TIR ? TIR_color : weight * AI_RGB_WHITE); 
                        const float drDepthAdderSpec = refractRoughnessConvert( AiShaderEvalParamFlt( p_dr_roughnessDepthMultiplier ) );
                        const float drAccumRoughness = rayState->updateDRAccumRoughnessReturnOrig(drDepthAdderSpec);

                        if ( do_TIR && rayState->ray_TIRDepth < JFND_MAX_TIR_DEPTH && specularRay.refr_bounces > 1)
                        {
                            rayState->ray_TIRDepth++;
                            specularRay.level--;
                            specularRay.refr_bounces--; 
                        }

                        if (sg->Rt == AI_RAY_CAMERA)
                            rayState->ray_energy_photon /= (float) data->gloss_samples;

                        AiStateSetMsgRGB(JFND_MSG_PHOTON_RGB, rayState->ray_energy_photon);

                        while ( AiSamplerGetSample(specularIterator, specular_sample) )
                        {
                            if (sharp_reflection) 
                            {
                                specular_sample[0] = 0.5f;
                                specular_sample[1] = 0.5f;
                            }

                            specularRay.dir = AiCookTorranceMISSample(brdf_data, (float) specular_sample[0], (float) specular_sample[1]);
                                
                            if (AiV3Dot(specularRay.dir, sg->Nf) > ZERO_EPSILON )
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

                        rayState->resetEnergyCache(cache_energy);
                        rayState->resetPhotonEnergyCache(cache_photonEnergy);               
                        rayState->resetDRAccumRoughness(drAccumRoughness);
                    }

                    // ---------------------------------------------------//
                    // Direct Specular / Reflection
                    // ---------------------------------------------------//

                    // decision point- direct specular
                    if ( traceSwitch.spec_dir )
                    {
                        // AtColor weight = AI_RGB_WHITE * fresnelTerm * rayState->media_specDirect.v[iinfo.m_higherPriority] * overallResultScale;
                        AtColor weight = AI_RGB_WHITE * fresnelTerm * rayState->media[iinfo.m_higherPriority].specDirect * overallResultScale;
                        if (do_TIR)
                            weight *= TIR_color;
                        
                        AiStateSetMsgBool(JFND_MSG_OPQ_SHADOW_MODE_BOOL, true);                        
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
                        AiStateSetMsgBool(JFND_MSG_OPQ_SHADOW_MODE_BOOL, false);
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
        if ( rayState->ray_invalidDepth < JFND_MAX_INVALID_DEPTH)
        {
            AtRay ray;
            AtScrSample sample;

            updateMediaInsideLists(media_id, iinfo.entering, media_inside_ptr, false);
            AiMakeRay(&ray, AI_RAY_REFRACTED, &sg->P, &sg->Rd, AI_BIG, sg);
            rayState->ray_invalidDepth++;
            ray.refr_bounces--;
            ray.level--;
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
        AiStateSetMsgBool(JFND_MSG_OPQ_SHADOW_MODE_BOOL, true);
        emission += AiShaderEvalParamRGB(p_emission_at_interfaces);
        if (iinfo.mediaExit)
        {
            emission += AiShaderEvalParamRGB(p_emission_at_exits);
        }
        else if (iinfo.mediaEntrance)
        {
            emission += AiShaderEvalParamRGB(p_emission_at_entrances);
        }
        AiStateSetMsgBool(JFND_MSG_OPQ_SHADOW_MODE_BOOL, false);
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

    AiStateSetMsgRGB(JFND_MSG_PHOTON_RGB,AI_RGB_BLACK); // this doesn't seem right, heberlein caustics related only

    if (msgs_are_valid)
    {
        AiStateSetMsgBool(JFND_MSG_VALID_BOOL, true);
        rayState->uncacheRayState(&rayStateCache);
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
