import pymel.core as pm
import mtoa.utils as utils
import mtoa.ui.ae.utils as aeUtils
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

class AEjf_photonTemplate(ShaderAETemplate):

	convertToMayaStyle = True

	def setup(self):

		self.beginScrollLayout()

		self.beginLayout("Caches", collapse=False)
		self.addControl("mode", 							label = "Mode" )
		self.addControl("frame", 							label = "Frame" )
		self.addControl("file_path", 						label = "File Path" )
		self.endLayout()

		self.beginLayout("Write", collapse=False)
		self.addControl("write_merge_photons", 				label = "Write Merge Photons" )
		self.addControl("write_remerge_photons", 			label = "Write Remerge Photons" )
		self.addControl("write_merge_radius", 				label = "Write Merge Radius" )		
		self.endLayout()

		self.beginLayout("Read", collapse=False)
		self.addControl("read_radius", 						label = "Read Radius" )
		self.addControl("diffuse_color", 					label = "Diffuse Color" )
		self.addControl("exposure", 						label = "Exposure" )
		self.addControl("refracted_intensity", 				label = "Refracted Intensity" )
		self.addControl("reflected_intensity", 				label = "Reflected Intensity" )
		self.endLayout()


		pm.mel.AEdependNodeTemplate(self.nodeName)

		self.addExtraControls()
		self.endScrollLayout()