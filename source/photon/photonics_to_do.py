
Ultimate way this works:

	*something* triggers a pre-render process
	this process puts jf_photon into write mode
	prerender process happens
	when finished, *something* triggers jf_photon to write
		perhaps an output driver

	perhaps:
		when a ray 
		create a camera from the light, with a spot light-emulating lens
		fire rays, check if the found object is a caustic emitter
			if so, trace and evaluate the shaders it finds





	jf_photon defaults (or is set into) read mode, same path





TO DO:


Caustics:

# Make jf_photon

# Writing to photon map
# 	Store these in some kind thread-safe container
# 		investigate thread-safe containers
# 		perhaps just a std::vector with a critical section





	# Demonstrate that writes photon maps
	# Demonstrate reading the photon maps
	# evaluate
	


Writing caustics:
	Investigate removing points outside of the bucket to avoid the grid lines



Reading Caustics:
	# Add a "read - Octree Visualization" mode - because it is neato spaghettio
	Differentiate between reflected and refracted caustics, give them seperate scales
	Diffuse color
	Exposure
	Normalize against radius
	Blackman-harris filtering




Data structure accelleration:
	# Control test: 
	# 	0:56
	# 	3.62 s build
	# Pointer (less copying vectors): 
	# 	0:52

	# size_t max_per_bucket = 100, size_t max_nesting = 16
	# 	0:51
	# 	2.63 s build

	# size_t max_per_bucket = 5, size_t max_nesting = 16
	# 	1:05
	# 	6.0 s build

	# size_t max_per_bucket = 20, size_t max_nesting = 32
	# 	0:51
	# 	3.8 s build

	# size_t max_per_bucket = 20, size_t max_nesting = 16
	# 	0:49
	# 	3.5 s build

	# size_t max_per_bucket = 20, size_t max_nesting = 8
	# 	0:49
	# 	1.0 s build

	# size_t max_per_bucket = 100, size_t max_nesting = 24
	# 	0:44
	# 	3.5 s build

	# Pass a pointer to not copy the lists so much
	# 	evaluate speedup
	Put photons into data structure container for contiguous memory reads
		evaluate speedup
