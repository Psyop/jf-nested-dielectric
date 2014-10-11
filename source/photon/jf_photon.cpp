#include <ai.h>

#include <fstream>
#include <vector>
#include <ctime>

//#include <algorithm>
//#include <iterator>
//#include <map>
//#include <cmath> 
//#include <Windows.h>


AI_SHADER_NODE_EXPORT_METHODS(jf_photon_methods);


const char * enum_modes[] =
{
	"Disabled",
	"Read",
	"Read and Visualize Octree",
	"Write",
	NULL
};

enum brdfs
{
	m_disabled,
	m_read,
	m_read_visualize,
	m_write,
};

typedef struct photon_type{
	AtColor energy;
	AtVector pos;
	AtUInt16 type;
	// AtVector dir;
	// float rl;
} photon_type;

typedef std::vector<photon_type> photon_cloud_type; // vector of photons with all the information a photon has

typedef std::vector<size_t> photon_list_type; // a list of photons in a photon cloud

typedef unsigned long long file_int;


float Log2( float n )  
{  
    // log(n)/log(2) is log2.  
    return (float) (log( n ) / log( 2.0f ));  
}

unsigned long long g_accelerator_count;

typedef class photon_accellerator_type{
	private:
		AtVector bounds_n; //bounds_n + {_len, _len, _len} is the positive bounds
		float _len;

		bool _has_sub_accells;
		unsigned char _recursion_level;
		unsigned char _max_recursion;
		unsigned short _photons_per_bucket_hint;
		// size_t _octant;

		photon_list_type * _photon_list;
		photon_accellerator_type * _sub_accells[8];

		photon_cloud_type * _target_photon_cloud;



	private:
		void init_photon_list() {
			_photon_list = new photon_list_type;
		}
		// Bounds creation

		void init_bounds(photon_cloud_type* photon_cloud) {
			_target_photon_cloud = photon_cloud;
			// _photons_per_bucket_hint = max_per_bucket;
			// _max_recursion = max_nesting;
			_recursion_level = 0;
			init_photon_list();

			AtVector measured_bounds_n = AI_V3_ZERO;
			AtVector measured_bounds_p = AI_V3_ZERO;

			for(size_t i = 0; i != _target_photon_cloud->size(); i++) {
				_photon_list->push_back(i);
				photon_type * photon = &_target_photon_cloud->at(i);
				if (i == 0) {
					measured_bounds_n = photon->pos;
					measured_bounds_p = photon->pos;
				} else {
					if (measured_bounds_n.x > photon->pos.x)
						measured_bounds_n.x = photon->pos.x;
					if (measured_bounds_n.y > photon->pos.y)
						measured_bounds_n.y = photon->pos.y;
					if (measured_bounds_n.z > photon->pos.z)
						measured_bounds_n.z = photon->pos.z;

					if (measured_bounds_p.x < photon->pos.x)
						measured_bounds_p.x = photon->pos.x;
					if (measured_bounds_p.y < photon->pos.y)
						measured_bounds_p.y = photon->pos.y;
					if (measured_bounds_p.z < photon->pos.z)
						measured_bounds_p.z = photon->pos.z;
				}
			}
			AtVector dim = measured_bounds_p - measured_bounds_n;

			_len = std::max(std::max(dim.x, dim.y), dim.z);
			bounds_n = measured_bounds_n;
		}

		void set_bounds(unsigned char octant, AtVector par_bounds_n, float par_len) {
			/*    0 means negative, 1 means positive
			 * 0: 000
			 * 1: 001
			 * 2: 010
			 * 3: 011
			 * 4: 100
			 * 5: 101
			 * 6: 110
			 * 7: 111
			*/
			float half_len = par_len/2.0f;
			AtVector offset_vector = {half_len, half_len, half_len};
			bool neg_x = octant % 2 == 0; //0, 2, 4
			bool neg_y = (octant/2) % 2 == 0; //0, 1, 4, 5
			bool neg_z = octant < 4;

			bounds_n = par_bounds_n;
			if (!neg_x) {
				bounds_n.x += half_len;
			}
			if (!neg_y) {
				bounds_n.y += half_len;
			}
			if (!neg_z) {
				bounds_n.z += half_len;
			}
			_len = half_len;
		}

		// Bounds checking
		bool within_bounds(AtVector photon_pos) {
			AtVector offset = {_len, _len, _len};
			AtVector bounds_p = bounds_n + offset;
			if (photon_pos.x >= bounds_n.x && 
				photon_pos.y >= bounds_n.y && 
				photon_pos.z >= bounds_n.z &&
				photon_pos.x < bounds_p.x &&
				photon_pos.y < bounds_p.y &&
				photon_pos.z < bounds_p.z
				) {
				return true;
			}
			return false;
		}


		bool within_range(const AtVector* photon_pos, float radius) {		
			AtVector offset = {_len, _len, _len};
			AtVector bounds_p = bounds_n + offset;
			if (photon_pos->x + radius > bounds_n.x && 
				photon_pos->y + radius > bounds_n.y && 
				photon_pos->z + radius > bounds_n.z &&
				photon_pos->x - radius < bounds_p.x &&
				photon_pos->y - radius < bounds_p.y &&
				photon_pos->z - radius < bounds_p.z
				) {
				return true;
			}
			return false;
		}


		void add_ID_to_bucket(size_t ID) {
			_photon_list->push_back(ID);
		}


		// Culling

		void cull_photons_in_bucket(photon_accellerator_type * accell, float radius, photon_cloud_type * cloud_out) {
			photon_list_type* photons = accell->_photon_list;
			if (accell->_len < radius) {
			// if (false) {
				photon_type refr_conglom;
				photon_type refl_conglom;
				refr_conglom.energy = AI_RGB_BLACK;
				refl_conglom.energy = AI_RGB_BLACK;
				refr_conglom.pos = AI_V3_ZERO;
				refl_conglom.pos = AI_V3_ZERO;

				size_t refr_conglom_count = 0;
				size_t refl_conglom_count = 0;

				for(size_t i = 0; i != photons->size(); i++) {
					photon_type * photon = &_target_photon_cloud->at( photons->at(i) );

					if (photon->type == AI_RAY_REFRACTED) {
						refr_conglom_count++;
						refr_conglom.energy += photon->energy;
						refr_conglom.pos += photon->pos;
					} else if (photon->type == AI_RAY_REFLECTED || photon->type == AI_RAY_GLOSSY) {
						refl_conglom_count++;
						refl_conglom.energy += photon->energy;
						refl_conglom.pos += photon->pos;
					}
				}

				if (refr_conglom_count > 0) {
					refr_conglom.pos /= (float) refr_conglom_count;
					refr_conglom.type = AI_RAY_REFRACTED;
					cloud_out->push_back(refr_conglom);
				}
				if (refl_conglom_count > 0) {
					refl_conglom.pos /= (float) refl_conglom_count;
					refl_conglom.type = AI_RAY_REFLECTED;
					cloud_out->push_back(refl_conglom);
				}
			} else {
				for(size_t i = 0; i != photons->size(); i++) {
					photon_type * photon = &_target_photon_cloud->at( photons->at(i) );
					cloud_out->push_back(*photon);
				}
			}
		}

		// void cull_photons_in_all_buckets(float radius, photon_cloud_type * cloud_out) {
		// 	for (size_t i = 0; i < all_nested_accells->size(); i++) {
		// 		photon_accellerator_type * accell = all_nested_accells->at(i);
		// 		if (accell->_has_sub_accells == false) {
		// 			// AiMsgWarning("Culling photons in subAccell %d", i);
		// 			cull_photons_in_bucket(all_nested_accells->at(i), radius, cloud_out);			
		// 		}
		// 	}
		// }


		// Build tree

		void initialize_child(photon_accellerator_type * child, unsigned char octant) {
			g_accelerator_count++;

			child->init_photon_list();
			child->set_bounds(octant, bounds_n, _len); //will set _pos, bounds_n and bounds_p

			// child->_octant = octant;
			child->_has_sub_accells = false;
			child->_recursion_level = _recursion_level + 1;
			child->_max_recursion = _max_recursion;
			child->_photons_per_bucket_hint = _photons_per_bucket_hint;			

			child->_target_photon_cloud = _target_photon_cloud;


			for (unsigned char i = 0; i < 8; i++) {
				child->_sub_accells[i] = NULL;
			}
		}

		// Recursive tree build
		void build_structure(bool cull = false, float cull_radius = 1.0f, photon_cloud_type * cull_cloud_out = NULL) {
			if (_photon_list->size() < _photons_per_bucket_hint || _recursion_level >= _max_recursion) {
				// Stop building. Either we're divided enough, or we're at our max recursion depth. 
				if (cull) {
					cull_photons_in_bucket(this, cull_radius, cull_cloud_out);
				}
				return;
			}
			// make 8 children and initialize them
			for (unsigned char i = 0; i < 8; i++) {
				_sub_accells[i] = new photon_accellerator_type;
				initialize_child(_sub_accells[i], i);
			}

			// distribute points
			for(size_t i = 0; i != _photon_list->size(); i++) {
				size_t photon_ID = _photon_list->at(i);
				for (unsigned char i = 0; i < 8; i++) {
					AtVector photon_position = _target_photon_cloud->at(photon_ID).pos;
					if (_sub_accells[i]->within_bounds(photon_position)) {
						_sub_accells[i]->add_ID_to_bucket(photon_ID);
						break;
					}
				}
			}
			_has_sub_accells = true;
			delete _photon_list;
			_photon_list = NULL;
			for (unsigned char i = 0; i < 8; i++) {
				_sub_accells[i]->build_structure(cull, cull_radius, cull_cloud_out);
				if (cull) {
					_sub_accells[i]->ripple_destroy();
					delete _sub_accells[i];
					_sub_accells[i] = NULL;
				}
			}
		}

	public:
		void get_photons_in_radius(photon_list_type* photon_list_out, const AtVector* pos, float radius) {
			if (!_has_sub_accells) {
				photon_list_out->insert(
					photon_list_out->end(), 
					_photon_list->begin(), 
					_photon_list->end()
					);
				return;
			}

			photon_accellerator_type * s_sub_accell = _sub_accells[0];
			std::vector<photon_accellerator_type *> accells_to_search;

			for (unsigned char i = 0; i < 8; i++) {
				if( _sub_accells[i]->within_range(pos, radius)) {
					accells_to_search.push_back(_sub_accells[i]);
				}
			}

			while (accells_to_search.size() > 0) {
				std::vector<photon_accellerator_type *> s_accells = accells_to_search;
				accells_to_search.clear();
				for (size_t i = 0; i < s_accells.size(); i++) {

					photon_accellerator_type* sub_accell = s_accells[i];

					if (sub_accell->_has_sub_accells) {
						// no points here, add sub accels to accells_to_search
						for (unsigned char g = 0; g < 8; g++) {
							if (sub_accell->_sub_accells[g]->within_range(pos, radius)) {
								accells_to_search.push_back(sub_accell->_sub_accells[g]);
							}
						}
					} else {
						// We have found points. 

						// boxes_count++;
						// photon_count += sub_accell->_photon_list->size();
						// max_depth = max( sub_accell->_recursion_level, max_depth);

						photon_list_out->insert(
							photon_list_out->end(), 
							s_accells[i]->_photon_list->begin(), 
							s_accells[i]->_photon_list->end()
							);
					}
				}
			}
			// AiMsgWarning("Found %d photons in %d buckets, max %d", photon_count, boxes_count, max_depth);
		}


		void build_while_culling(photon_cloud_type* photon_cloud, float radius, photon_cloud_type * cloud_out) {
			init_bounds(photon_cloud);
			g_accelerator_count = 1;

			int subdivisions = (int) (Log2(_len/radius) + 2.0f);
			subdivisions = std::min(subdivisions, 24);

			_photons_per_bucket_hint = 1;
			_max_recursion = (unsigned char) subdivisions;
			build_structure(true, radius, cloud_out);

			AiMsgInfo("  %d kilo-substructures. 8^%d max division. ", 
				g_accelerator_count/1000, subdivisions);

			// cull_photons_in_all_buckets(radius, cloud_out);		
		}


		void build(photon_cloud_type* photon_cloud, float radius_hint) {
			unsigned short photons_per_bucket_hint = 32;
			unsigned short max_nesting = 11;
			g_accelerator_count = 1;

			init_bounds(photon_cloud);

			unsigned short subdivisions_hint = (int) (Log2(_len/radius_hint) - 2.0f);
			unsigned short subdivisions = std::min(subdivisions_hint, max_nesting);

			_photons_per_bucket_hint = photons_per_bucket_hint;
			_max_recursion = (unsigned char) subdivisions;
			build_structure();

			AiMsgInfo("  %d  kilo-substructures. 8^%d max divisions (optimized for radius %f). ", 
				g_accelerator_count/1000, subdivisions, radius_hint);
		}


		void destroy_structure() {
			ripple_destroy();

			if (_photon_list != NULL) {
				delete _photon_list;
			}

			g_accelerator_count = 0;
		}


		void ripple_destroy() {
			if (_photon_list != NULL) {
				delete _photon_list;
			}

			for (unsigned char i = 0; i < 8; i++) {
				if (_sub_accells[i] != NULL) {
					_sub_accells[i]->ripple_destroy();
					delete _sub_accells[i];
				}
			}
		}
} photon_accellerator_type;


float blackman_harris(float distance, float radius) {
	float x = distance/radius;

	float a0 = 0.35875f;
	float a1 = 0.48829f;
	float a2 = 0.14128f;
	float a3 = 0.01168f;

	float weight  = a0 + a1*cos(1.0f * AI_PI * x) + a2*cos(2.0f * AI_PI * x) + a3*cos(4.0f * AI_PI * x);

	return weight;
}


struct ShaderData{
	std::string file_name;
	bool abort;
	AtArray * write_thread_clouds;
	photon_cloud_type * read_cloud;
	photon_accellerator_type * read_cloud_accelerator;
	float sampling_normalizer;
	std::string aov_refr_caustics;
	std::string aov_refl_caustics;
};


enum jf_photonParams {
	p_mode,
	p_file_path,
	p_frame,
	p_write_merge_photons,
	p_write_remerge_photons,
	p_write_merge_radius,
	p_read_radius,
	p_diffuse_color,
	p_exposure,
	p_refracted_intensity,
	p_reflected_intensity,
};

 
node_parameters {
	AiParameterEnum("mode", m_read, enum_modes);
	AiParameterStr("file_path", "c:/photonics.[Frame].caustics");
	AiParameterInt("frame", 0);
	AiParameterBool("write_merge_photons", true);
	AiParameterBool("write_remerge_photons", true);
	AiParameterFlt("write_merge_radius", 0.005f);
	AiParameterFlt("read_radius", 0.1f);
	AiParameterRGB("diffuse_color", 0.7f, 0.7f, 0.7f);
	AiParameterFlt("exposure", 0.0f);
	AiParameterFlt("refracted_intensity", 1.0f);
	AiParameterFlt("reflected_intensity", 1.0f);
}


size_t init_cloud_size = sizeof(photon_type) * 100000;


node_initialize {
	ShaderData* data = new ShaderData;
	AiNodeSetLocalData(node,data);
	data->abort = false;

	AtNode * render_options = AiUniverseGetOptions();

	int mode = AiNodeGetInt(node, "mode");

	std::string file_path = AiNodeGetStr(node, "file_path");
	int frame = AiNodeGetInt(node, "frame");
	char frame_string[4];
 	sprintf(frame_string, "%04d", frame);
	std::string token = "[Frame]";
	if (file_path.find(token) != std::string::npos) {	
		file_path.replace(file_path.find(token),token.length(),frame_string);
	}
	data->file_name = file_path;

	if (mode == m_write) {
		// AiMsgWarning("Setting up %d photon containers for the threads.", AiNodeGetInt(render_options, "threads"));
		data->write_thread_clouds = AiArrayAllocate(AiNodeGetInt(render_options, "threads"), 1, AI_TYPE_POINTER);
		for (AtUInt32 i = 0; i < data->write_thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = new photon_cloud_type;
			cloud->reserve(init_cloud_size);
			AiArraySetPtr(data->write_thread_clouds, i, cloud );
		}

		return;
	} 

	if (mode == m_read || mode == m_read_visualize) {
		//const char* char_path = AiNodeGetStr(node, "file_path");
		std::string string_path = data->file_name;
		std::ifstream infile (string_path.c_str(), std::ios::binary);


		AiMsgWarning("Reading Photon Cloud from: %s ", string_path.c_str());
		if (!infile.good()) {
			data->abort = true;
			AiMsgError("Unable to read file! Check for invalid paths or bad permissions or something.");
			return;
		}
		
		infile.seekg (0, infile.end);
		file_int length = infile.tellg();

		file_int photon_size = sizeof(photon_type);
		file_int num_photons = length/photon_size;
		file_int mb = 1024 * 1024;
		file_int chunk_size = (128 * mb) ;
		file_int num_photons_in_chunk = (chunk_size / photon_size) + 1;
		chunk_size = photon_size * num_photons_in_chunk; //No remainders here, please. 

		data->read_cloud = new photon_cloud_type;
		photon_cloud_type * read_cloud = data->read_cloud; //Alias of the cloud in data->read_cloud;

		for (file_int i = 0; i < length; i+= chunk_size) {
			file_int read_bytes = std::min(length - i, chunk_size);
			file_int read_photons = read_bytes / photon_size;

			infile.seekg (i, infile.beg);

			AiMsgInfo("  %d mb chunk, %d kilophotons.", (int) read_bytes/mb, read_photons/1000);
			photon_type* photon_array = new photon_type[read_photons];
			char * buffer = (char*)(photon_array);

			infile.read(buffer, read_bytes);
			read_cloud->insert(read_cloud->end(), &photon_array[0], &photon_array[read_photons]);

			delete photon_array;
		}


		AiMsgInfo("  Read %d mb, %d kilophotons.", length/mb, num_photons/1000);
		// photon_type* last_photon = &photon_array[num_photons - 1];
		// AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );

		const clock_t photon_process_time = clock();

		if (num_photons != read_cloud->size()) {
			AiMsgError("Error in photon read. %d in file, %d in memory.", num_photons, read_cloud->size());
		} else {
			AiMsgInfo("  All photons accounted for.");
		}
		// last_photon = &read_cloud->at(num_photons - 1);
		// AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );

		photon_accellerator_type * accel = new photon_accellerator_type;
		float radius_hint = AiNodeGetFlt(node, "read_radius");
		accel->build(read_cloud, radius_hint);
		data->read_cloud_accelerator = accel;

		AiMsgInfo("Octree completed: %f seconds", (float(clock() - photon_process_time) /  CLOCKS_PER_SEC));
	}
}

 
node_update {
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	if (data->abort) {
		return;
	}

	int mode = AiNodeGetInt(node, "mode");
	AtNode * render_options = AiUniverseGetOptions();

	if (mode == m_write) {
		for (AtUInt32 i = 0; i < data->write_thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, i));
			cloud->clear();
			cloud->reserve(init_cloud_size);
		}

		unsigned long long AA_samples = AiNodeGetInt(render_options, "AA_samples");
		unsigned long long total_samples = AA_samples * AA_samples * AiNodeGetInt(render_options, "xres") * AiNodeGetInt(render_options, "yres");
		unsigned long long expected_sampling_baseline = 16777216;
		data->sampling_normalizer = (float) ((double) expected_sampling_baseline / (double) total_samples);
		AiMsgWarning("Based on expected %d kilosamples, normalization factor is %f.", total_samples/1000, data->sampling_normalizer);

		return;
	}

	if (mode == m_read || mode == m_read_visualize) {
		// STUB
		// Still can't find anything that needs doing here
	}
}


node_finish {
	if (AiNodeGetLocalData(node) == NULL) {
		return;
	}

	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

	if (data->abort) {
		return;
	}

	int mode = AiNodeGetInt(node, "mode");

	if (mode == m_write) {
		// const char* char_path = AiNodeGetStr(node, "file_path");
		std::string string_path = data->file_name;
		std::ofstream outfile (string_path.c_str(), std::ios::binary);

		AiMsgWarning("Writing Photon Cloud to: %s ", string_path.c_str());
		if (!outfile.good()) {
			AiMsgError("Unable to write file! Check for invalid paths or bad permissions or something.");
			return;
		}

		float merge_radius = AiNodeGetFlt(node, "write_merge_radius");
		bool merge = AiNodeGetBool(node, "write_merge_photons");
		bool remerge = AiNodeGetBool(node, "write_remerge_photons");

		size_t photon_count = 0;
		size_t orig_photon_count = 0;

		if (merge) {
			AtUInt32 threads = data->write_thread_clouds->nelements;
			photon_cloud_type compiled_cloud;

			AiMsgWarning("Reducing (merging) photons from subclouds:");

			for (AtUInt32 i = 0; (i < threads) && (i < threads); i++) {
				AtUInt32 thread_ID = i;
				photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, thread_ID));
				
				if (cloud->size() > 0) {
					photon_accellerator_type octree;

					AiMsgInfo("  Building octree for thread %d:", thread_ID);
					size_t prev_size = compiled_cloud.size();

					octree.build_while_culling(cloud, merge_radius, &compiled_cloud);
					octree.destroy_structure();
					AiMsgInfo("  %d -> %d kilophotons", cloud->size()/1000, (compiled_cloud.size() - prev_size)/1000);
					AiMsgInfo("(%d mb of compiled photons)", (compiled_cloud.size() * sizeof(photon_type))/(1024*1024));
			
					orig_photon_count += cloud->size();

					cloud->clear();
					delete cloud;
				}
			}

			if (compiled_cloud.size() > 0) {
				if (remerge) {						
					AiMsgWarning("Rereducing compiled photons:");
					photon_cloud_type cloud_out;
					photon_accellerator_type octree;
					octree.build_while_culling(&compiled_cloud, merge_radius, &cloud_out);
					octree.destroy_structure();

					AiMsgInfo("  %d -> %d kilophotons (rereduction)", compiled_cloud.size()/1000, cloud_out.size()/1000);
					outfile.write((const char*)&cloud_out.at(0), (file_int) sizeof(photon_type) * (file_int) cloud_out.size());
					photon_count += cloud_out.size();
				} else {
					outfile.write((const char*)&compiled_cloud.at(0), (file_int) sizeof(photon_type) * (file_int) compiled_cloud.size());
					photon_count += compiled_cloud.size();
				}	

			}
		} else {
			// Naive Write
			AtUInt32 threads = data->write_thread_clouds->nelements;

			for (AtUInt32 i = 0; i < threads; i += 1) { 
				AtUInt32 thread_ID = i;
				photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, thread_ID));

				if (cloud->size() > 0) {
					outfile.write((const char*)&cloud->at(0), (file_int) sizeof(photon_type) * (file_int) cloud->size());
					photon_count += cloud->size();		
				}
				delete cloud;
			}

		}

		float cloud_mb = (float) (photon_count * sizeof(photon_type)) / (1024.0f*1024.0f);
		float orig_cloud_mb = (float) (orig_photon_count * sizeof(photon_type)) / (1024.0f*1024.0f);
		AiMsgWarning("Photon cloud: %f mb, %d kilophotons.", cloud_mb, photon_count/1000);
		if (merge) {
			AiMsgInfo("  reduced from: %f mb, %d kilophotons.", orig_cloud_mb, orig_photon_count/1000);
		}

		outfile.close();
		AiArrayDestroy(data->write_thread_clouds);
	}

	if ((mode == m_read || mode == m_read_visualize)) {
		AiMsgWarning("Photon Octree: Destroying...");
		data->read_cloud_accelerator->destroy_structure();
		delete data->read_cloud_accelerator;
		AiMsgWarning("Photon Cloud: Destroying...");
		delete data->read_cloud;
	}

	delete data;
}
 

shader_evaluate {
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
	int mode = AiShaderEvalParamEnum(p_mode);
	sg->out.RGBA = AI_RGBA_BLACK;

	if (data->abort) {
		return;
	}
	
	if (mode == m_write) {
		AtColor photon_energy = AI_RGB_BLACK;
		bool is_photon = AiStateGetMsgRGB( "photon_energy", &photon_energy );

		if (is_photon) {
			photon_type photon = {photon_energy * data->sampling_normalizer, sg->P, sg->Rt};
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, sg->tid));
			cloud->push_back(photon);
			sg->out.RGB = photon_energy;
			sg->out.RGBA.a = 1.0f;
		} else {
			sg->out.RGBA = AI_RGBA_RED / 5;
		}
	}	

	if (mode == m_read || mode == m_read_visualize) { 
		float radius = AiShaderEvalParamFlt(p_read_radius);

		if (mode == m_read) {
			AtColor kd = AiShaderEvalParamRGB(p_diffuse_color);
			float k_refracted = AiShaderEvalParamFlt(p_refracted_intensity);
			float k_reflected = AiShaderEvalParamFlt(p_reflected_intensity);

			photon_list_type photon_IDs;
			data->read_cloud_accelerator->get_photons_in_radius(&photon_IDs, &sg->P, radius);

			AtColor refr_energy = AI_RGB_BLACK;
			AtColor refl_energy = AI_RGB_BLACK;
			// float refr_weight = 0.0f;
			// float refl_weight = 0.0f;

			for(size_t i = 0; i != photon_IDs.size(); i++) {
				photon_type * photon = &data->read_cloud->at(photon_IDs[i]);
				float dist = AiV3Dist(sg->P, photon->pos);
				if (dist < radius) {
					float weight = blackman_harris(dist, radius);
					if (photon->type == AI_RAY_REFRACTED) {
						refr_energy += photon->energy * weight;
						// refr_weight += weight;
					} else if (photon->type == AI_RAY_REFLECTED || photon->type == AI_RAY_GLOSSY) {
						refl_energy += photon->energy * weight;
						// refl_weight += weight;
					}
				}
			}
			float exposure = (float) std::pow(2, AiShaderEvalParamFlt(p_exposure));
			float normalization_factor = (radius * radius * 10000000.0f);
			refr_energy /= normalization_factor;
			refl_energy /= normalization_factor;

			refr_energy *= kd * k_refracted * exposure;
			refl_energy *= kd * k_reflected * exposure;

			sg->out.RGB = refr_energy + refl_energy;

			sg->out.RGBA.a = 1.0f;
		} else if (mode == m_read_visualize) {
			photon_list_type photon_IDs;
			data->read_cloud_accelerator->get_photons_in_radius(&photon_IDs, &sg->P, radius);
			float out_value = ((float) photon_IDs.size() / 300.0f);
			sg->out.RGBA = AI_RGBA_WHITE * out_value;
			sg->out.RGBA.a = 1.0f;
		}		
	}
}

 
node_loader {
	if (i > 0)
		return false;

	node->methods = jf_photon_methods;
	node->output_type = AI_TYPE_RGBA;
	node->name = "jf_photon";
	node->node_type = AI_NODE_SHADER;
	strcpy_s(node->version, AI_VERSION);
	return true;
}