
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
	# Differentiate between reflected and refracted caustics, give them seperate scales
	# Diffuse color
	# Exposure
	# Normalize against radius
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

	was: 40 s
	changed filtering to avoid a sqrt: 37 s
	removed the minor optimization and got to 35s



	# Pass a pointer to not copy the lists so much
	# 	evaluate speedup
	Put photons into data structure container for contiguous memory reads
		evaluate speedup



reduce:
	309mb, 1:10 to build and write
	160mb, 1:27 to build and build octrees and write


'''
00:06:47  1968MB         | releasing resources
00:06:47  1965MB WARNING |  Writing Photon Cloud to: C:/Users/Jonah/development/photons/photon_map_prism_test2.bin 
00:06:47  1965MB         |    Building octree for thread 0:
00:06:53  2963MB         |    5692096 sub-structures. 2^11 division. 
00:06:54  2045MB         |    9574 -> 1097 kilophotons
00:06:54  2045MB         |    b mb of compiled photons.
00:06:54  1789MB         |    Building octree for thread 1:
00:06:59  2760MB         |    5574080 sub-structures. 2^11 division. 
00:07:00  1862MB         |    9211 -> 1115 kilophotons
00:07:00  1862MB         |    b mb of compiled photons.
00:07:00  1615MB         |    Building octree for thread 2:
00:07:06  2487MB         |    4926344 sub-structures. 2^11 division. 
00:07:07  1680MB         |    9453 -> 952 kilophotons
00:07:07  1680MB         |    b mb of compiled photons.
00:07:07  1427MB         |    Building octree for thread 3:
00:07:12  2235MB         |    4551760 sub-structures. 2^11 division. 
00:07:13  1488MB         |    9140 -> 958 kilophotons
00:07:13  1488MB         |    b mb of compiled photons.
00:07:13  1243MB         |    Building octree for thread 4:
00:07:18  2205MB         |    5526144 sub-structures. 2^11 division. 
00:07:19  1314MB         |    9129 -> 1064 kilophotons
00:07:19  1314MB         |    b mb of compiled photons.
00:07:19  1070MB         |    Building octree for thread 5:
00:07:24  1894MB         |    4674936 sub-structures. 2^11 division. 
00:07:25  1129MB         |    8768 -> 866 kilophotons
00:07:25  1129MB         |    b mb of compiled photons.
00:07:25   895MB         |    Building octree for thread 6:
00:07:30  1683MB         |    4471672 sub-structures. 2^11 division. 
00:07:31   952MB         |    8489 -> 811 kilophotons
00:07:31   952MB         |    b mb of compiled photons.
00:07:31   725MB         |    Building octree for thread 7:
00:07:36  1570MB         |    4833856 sub-structures. 2^11 division. 
00:07:37   784MB         |    8657 -> 847 kilophotons
00:07:37   784MB         |    b mb of compiled photons.
00:07:37   552MB WARNING |  Rereducing compiled photons:
00:07:43  2296MB         |    10523752 sub-structures. 2^11 division. 
00:07:44   734MB WARNING |  Rereduction: 7713 -> 3647 kilophotons
00:07:44   430MB         |  Photon cloud: 97.397102 mb, 3647 kilophotons.
00:07:44   430MB         |  Photon cloud: reduced from: 1933.971069 mb, 3647438 kilophotons.
00:07:44   426MB         | Arnold shutdown
