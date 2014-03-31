# ------------------------------------------------------------------
# JF Nested Dielectric By Jonah Friedman
# 1.0.2
# Copyright (c) 2014, Psyop Media Company, LLC and Jonah Friedman 
# Open sourced under the 3 clause BSD license, see license.txt

import pymel.core as pm
import mtoa.utils as utils
import mtoa.ui.ae.utils as aeUtils
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

class AEjf_nested_dielectricTemplate(ShaderAETemplate):

	convertToMayaStyle = True

	def setup(self):
		# Add the shader swatch to the AE
		self.addSwatch()

		self.beginScrollLayout()

		# Add a list that allows to replace the shader for other one
		self.addCustom('message', 'AEshaderTypeNew', 'AEshaderTypeReplace')


		self.beginLayout("Basic Dielectric Properties", collapse=True)

		self.beginLayout("Basics", collapse=False)
		self.addControl("mediumPriority", 					label = "Medium Priority" )
		self.addControl("mediumIOR", 						label = "Medium IOR" )		
		self.endLayout()

		self.beginLayout("Transmission", collapse=False)
		self.addControl("mediumTransmittance", 				label = "Medium Transmittance" )
		self.addControl("mediumTransmittance_scale", 		label = "Medium Transmittance Scale" )
		self.endLayout()

		self.beginLayout("Polarizing Filter", collapse=True)
		self.addControl("polarize", 						label = "Polarize" )
		self.addControl("polarization_angle", 				label = "Polarization Angle" )
		self.endLayout()

		self.endLayout() # end basics



		self.beginLayout("Refraction", collapse=True)


		self.addControl("refraction_scale", 				label = "Refraction Scale" )
		self.addControl("btdf", 							label = "BTDF" )
		self.addControl("refraction_roughness_u", 			label = "Roughness U" )
		self.addControl("refraction_roughness_v", 			label = "Roughness V" )

		self.addControl("indirect_refraction", 				label = "Indirect Refraction" )

		self.beginLayout("Direct Refraction", collapse=False)
		self.addControl("direct_refraction", 				label = "Direct Refraction" )


		self.addControl("dr_use_refraction_btdf", 			label = "Use Refraction BTDF" )
		self.addControl("dr_roughnessOffset", 				label = "Roughness Offset" )

		self.addControl("dr_roughnessDepthMultiplier", 		label = "Roughness Depth Multiplier" )
		self.addControl("dr_roughnessDepthAdder", 			label = "Roughness Depth Adder" )

		self.beginLayout("Direct Refraction Manual BTDF", collapse=True)
		self.addControl("dr_btdf", 							label = "BTDF" )
		self.addControl("dr_roughness_u", 					label = "Roughness U" )
		self.addControl("dr_roughness_v", 					label = "Roughness V" )
		self.endLayout()

		
		self.beginLayout("Direct Refraction Second Lobe", collapse=False)
		self.addControl("dr_second_scale", 					label = "Second Lobe Scale" )
		self.addControl("dr_second_roughnessMultiplier", 	label = "Roughness Multiplier" )
		self.endLayout()

		self.endLayout() # end Direct Refraction


		self.endLayout() # end refraction



		self.beginLayout("Dispersion", collapse=True)

		self.addControl("disperse", 						label = "Disperse" )
		self.addControl("dispersion", 						label = "Dispersion" )

		self.beginLayout("Spectral Conversion", collapse=False)
		self.addControl("spectral_distribution", 			label = "Spectral Distribution" )
		self.addControl("spectral_gamut", 					label = "Gamut" )
		self.addControl("spectral_clamp_negative_colors", 	label = "Clamp Negative Colors" )
		self.addControl("spectral_saturation", 				label = "Spectral Saturation" )
		self.endLayout()

		self.endLayout() # end dispersion



		self.beginLayout("Specular", collapse=True)

		self.beginLayout("BRDF", collapse=False)
		self.addControl("specular_scale", 					label = "Specular Scale" )
		self.addControl("brdf", 							label = "BRDF" )
		self.addControl("specular_roughness_u", 			label = "Roughness U" )
		self.addControl("specular_roughness_v", 			label = "Roughness V" )
		self.endLayout()

		self.beginLayout("Components", collapse=True)
		self.addControl("direct_specular", 					label = "Direct Specular" )
		self.addControl("indirect_specular", 				label = "Indirect Specular" )
		self.endLayout()

		self.beginLayout("Ray Type (for trace depth and ray switches)", collapse=True)
		self.addControl("specular_ray_type", 				label = "Specular Ray Type" )
		self.endLayout()

		self.endLayout()



		self.beginLayout("Other", collapse=True)

		self.beginLayout("Optimization", collapse=True)
		self.addControl("enable_internal_reflections", 		label = "Internal Reflections" )
		self.addControl("russian_roulette_probability", 	label = "Russian Roulette Probability" )
		self.addControl("energy_cutoff_exponent", 			label = "Energy Cutoff Exponent" )
		self.endLayout()

		self.beginLayout("Anisotropy", collapse=True)
		self.addControl("blur_anisotropic_poles", 			label = "Blur Anisotropic Poles" )
		self.addControl("ward_tangent", 					label = "Ward Tangent" )
		self.endLayout()

		self.beginLayout("Sky Behaviors", collapse=False)
		self.addControl("refract_skydomes", 				label = "Refract Skydome Lights" )
		self.addControl("refract_skies",					label = "Refract Skies (environment shaders)" )
		self.addControl("reflect_skydomes", 				label = "Reflect Skydome Lights" )
		self.addControl("reflect_skies", 					label = "Refract Skies (environment shaders)" )
		self.endLayout()

		self.beginLayout("Emission", collapse=False)
		self.addControl("emission_at_interfaces", 			label = "Emit at Interfaces" )
		self.addControl("emission_at_exits", 				label = "Emit at Exits" )
		self.addControl("emission_at_entrances", 			label = "Emit at Entrances" )
		self.endLayout()

		self.endLayout()


		self.beginLayout("Shadows and Caustics", collapse=True)

		self.addControl("shadow_mode", 						label = "Shadow Mode" )
		self.addControl("caustic_mode", 					label = "Caustic Mode" )

		self.addControl("caustic_refractions_direct", 		label = "Direct Refraction" )
		self.addControl("caustic_dr_roughness_offset", 		label = "Roughness Offset" )
		self.addControl("caustic_refractions_indirect", 	label = "Indirect Refraction" )
		self.addControl("caustic_TIR", 						label = "TIR" )
		self.addControl("caustic_internal_speculars", 		label = "Internal Reflections" )
		self.addControl("caustic_dispersion", 				label = "Dispersion" )

		self.addControl("caustic_specular_direct", 			label = "Direct Specular" )
		self.addControl("caustic_specular_indirect", 		label = "Indirect Specular" )

		self.addControl("caustic_max_distance", 			label = "Max Distance" )

		self.endLayout()
		

		self.beginLayout("AOVs", collapse=True)

		self.addControl("aov_direct_refraction", 			label = "aov_direct_refraction" )
		self.addControl("aov_indirect_refraction", 			label = "aov_indirect_refraction" )
		self.addControl("aov_direct_specular", 				label = "aov_direct_specular" )
		self.addControl("aov_indirect_specular", 			label = "aov_indirect_specular" )

		self.endLayout()


		pm.mel.AEdependNodeTemplate(self.nodeName)

		self.addExtraControls()
		self.endScrollLayout()