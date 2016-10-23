#include <algorithm>
#include <ai.h>
#include <fstream>
#include <vector>
#include <ctime>

AI_SHADER_NODE_EXPORT_METHODS(jf_photon_methods);

const char * enum_modes[] =
{
	"Disabled",
	"Read",
	"Read and Visualize Octree",
	"Write",
	NULL
};

enum modes
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
} photon_type;

typedef std::vector<photon_type> photon_cloud_type; // vector of photons with all the information a photon has

typedef std::vector<size_t> photon_list_type; // a list of photons in a photon cloud

typedef unsigned long long file_int;


float Log2( float n )  
{  
    // log(n)/log(2) is log2.  
    return (float) (log( n ) / log( 2.0f ));  
}

typedef class photon_octree_type{
	public:
		AtVector bounds_n; //bounds_n + {len, len, len} is the positive bounds
		float len;

	private:
		bool _has_sub_octrees;
		unsigned char _recursion_level;
		unsigned char _max_recursion;
		unsigned short _photons_per_bucket_hint;

		photon_list_type * _photon_list;
		photon_octree_type * _sub_octrees[8];
		unsigned long long * octree_count;

		photon_cloud_type * _target_photon_cloud;



	public:
		void init_bounds_from_vectors(AtVector* bounds_n_out, AtVector* bounds_p_out) {
			AtVector dim = *bounds_p_out - *bounds_n_out;
			len = std::max(std::max(dim.x, dim.y), dim.z);
			bounds_n = *bounds_n_out;

			// 1% bounds padding. Cheap way to defend againt certain artifacts that have cropped up due to floating point errors.. I think.
			len *= 1.02f;
			bounds_n -= dim * 0.01f;
		}

		void measure_bounds_from_cloud(photon_cloud_type* photon_cloud, AtVector* bounds_n_out, AtVector* bounds_p_out) {
			bool uninitialized = (*bounds_p_out == AI_V3_ZERO && *bounds_n_out == AI_V3_ZERO);

			for(size_t i = 0; i != photon_cloud->size(); i++) {
				photon_type * photon = &photon_cloud->at(i);
				if (uninitialized && i == 0) {
					*bounds_n_out = photon->pos;
					*bounds_p_out = photon->pos;
				} else {
					if (bounds_n_out->x > photon->pos.x)
						bounds_n_out->x = photon->pos.x;
					if (bounds_n_out->y > photon->pos.y)
						bounds_n_out->y = photon->pos.y;
					if (bounds_n_out->z > photon->pos.z)
						bounds_n_out->z = photon->pos.z;

					if (bounds_p_out->x < photon->pos.x)
						bounds_p_out->x = photon->pos.x;
					if (bounds_p_out->y < photon->pos.y)
						bounds_p_out->y = photon->pos.y;
					if (bounds_p_out->z < photon->pos.z)
						bounds_p_out->z = photon->pos.z;
				}
			}

			init_bounds_from_vectors(bounds_n_out, bounds_p_out);
		}


	private:
		void init_photon_list() {
			_photon_list = new photon_list_type;
		}
		// Bounds creation

		void init_top_tree(photon_cloud_type* photon_cloud, bool measure_bounds) {
			_target_photon_cloud = photon_cloud;
			octree_count = new unsigned long long;
			*octree_count = 1;

			_recursion_level = 0;
			init_photon_list();

			for(size_t i = 0; i != photon_cloud->size(); i++) {
				_photon_list->push_back(i);
				photon_type * photon = &photon_cloud->at(i); 
			}
			if (measure_bounds) {
				AtVector measured_bounds_n = AI_V3_ZERO;
				AtVector measured_bounds_p = AI_V3_ZERO;
				measure_bounds_from_cloud(photon_cloud, &measured_bounds_n, &measured_bounds_p);	
			}
		}

		void set_bounds(unsigned char octant, AtVector par_bounds_n, float parlen) {
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
			float halflen = parlen/2.0f;
			AtVector offset_vector = {halflen, halflen, halflen};
			bool neg_x = octant % 2 == 0; //0, 2, 4
			bool neg_y = (octant/2) % 2 == 0; //0, 1, 4, 5
			bool neg_z = octant < 4;

			bounds_n = par_bounds_n;
			if (!neg_x) {
				bounds_n.x += halflen;
			}
			if (!neg_y) {
				bounds_n.y += halflen;
			}
			if (!neg_z) {
				bounds_n.z += halflen;
			}
			len = halflen;
		}

		// Bounds checking
		bool within_bounds(AtVector photon_pos) {
			AtVector offset = {len, len, len};
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
			AtVector offset = {len, len, len};
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

		void cull_photons_in_bucket(photon_octree_type * tree, float radius, photon_cloud_type * cloud_out) {
			photon_list_type* photons = tree->_photon_list;
			if (tree->len < radius) {
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

		// Build tree

		void initialize_child(photon_octree_type * child, unsigned char octant) {
			child->init_photon_list();
			child->set_bounds(octant, bounds_n, len); //will set _pos, bounds_n and bounds_p

			child->octree_count = octree_count;
			++ *child->octree_count;

			child->_has_sub_octrees = false;
			child->_recursion_level = _recursion_level + 1;
			child->_max_recursion = _max_recursion;
			child->_photons_per_bucket_hint = _photons_per_bucket_hint;			

			child->_target_photon_cloud = _target_photon_cloud;

			for (unsigned char i = 0; i < 8; i++) {
				child->_sub_octrees[i] = NULL;
			}
		}

		// Recursive tree build
		void build_structure(bool cull = false, float cull_radius = 1.0f, photon_cloud_type * cull_cloud_out = NULL) {
			if (*octree_count > 1 && (_photon_list->size() < _photons_per_bucket_hint || _recursion_level >= _max_recursion)) {
				// Stop building. Either we're divided enough, or we're at our max recursion depth. 
				if (cull) {
					cull_photons_in_bucket(this, cull_radius, cull_cloud_out);
				}
				return;
			}

			unsigned long long octree_size = *octree_count * sizeof(photon_octree_type);
			unsigned long long accel_mb_estimate = (octree_size)/(1024 * 1024);
			if (!cull && accel_mb_estimate > (unsigned long) ( 16 * 1024 )) {
				// This is only a safety mechanism.
				AiMsgError("JF Photon: More than 16 gb of octree. Radius may be too small. Stopping construction.");
				return;
			}
			// make 8 children and initialize them
			for (unsigned char i = 0; i < 8; i++) {
				_sub_octrees[i] = new photon_octree_type;
				initialize_child(_sub_octrees[i], i);
			}

			// distribute points
			for(size_t i = 0; i != _photon_list->size(); i++) {
				size_t photon_ID = _photon_list->at(i);
				for (unsigned char i = 0; i < 8; i++) {
					AtVector photon_position = _target_photon_cloud->at(photon_ID).pos;
					if (_sub_octrees[i]->within_bounds(photon_position)) {
						_sub_octrees[i]->add_ID_to_bucket(photon_ID);
						break;
					}
				}
			}
			_has_sub_octrees = true;
			delete _photon_list;
			_photon_list = NULL;
			for (unsigned char i = 0; i < 8; i++) {
				_sub_octrees[i]->build_structure(cull, cull_radius, cull_cloud_out);
				if (cull) {
					_sub_octrees[i]->ripple_destroy();
					delete _sub_octrees[i];
					_sub_octrees[i] = NULL;
				}
			}
		}

	public:
		void get_photons_in_radius(photon_list_type* photon_list_out, const AtVector* pos, float radius) {
			if (!_has_sub_octrees) {
				photon_list_out->insert(
					photon_list_out->end(), 
					_photon_list->begin(), 
					_photon_list->end()
					);
				return;
			}

			photon_octree_type * s_sub_tree = _sub_octrees[0];
			std::vector<photon_octree_type *> octrees_to_search;

			for (unsigned char i = 0; i < 8; i++) {
				if( _sub_octrees[i]->within_range(pos, radius)) {
					octrees_to_search.push_back(_sub_octrees[i]);
				}
			}

			while (octrees_to_search.size() > 0) {
				std::vector<photon_octree_type *> s_octrees = octrees_to_search;
				octrees_to_search.clear();
				for (size_t i = 0; i < s_octrees.size(); i++) {

					photon_octree_type* sub_tree = s_octrees[i];

					if (sub_tree->_has_sub_octrees) {
						// no points here, add sub accels to octrees_to_search
						for (unsigned char g = 0; g < 8; g++) {
							if (sub_tree->_sub_octrees[g]->within_range(pos, radius)) {
								octrees_to_search.push_back(sub_tree->_sub_octrees[g]);
							}
						}
					} else {
						// We have found points. 
						photon_list_out->insert(
							photon_list_out->end(), 
							s_octrees[i]->_photon_list->begin(), 
							s_octrees[i]->_photon_list->end()
							);
					}
				}
			}
		}


		void build_while_culling(photon_cloud_type* photon_cloud, float radius, photon_cloud_type * cloud_out, AtVector* bounds_n_in = NULL, AtVector* bounds_p_in = NULL) {
			if (bounds_n_in != NULL && bounds_p_in != NULL) {
				init_top_tree(photon_cloud, false);
				init_bounds_from_vectors(bounds_n_in, bounds_p_in);
			} else {
				init_top_tree(photon_cloud, true);				
			}

			int subdivisions = (int) (Log2(len/radius) + 2.0f);
			subdivisions = std::min(subdivisions, 24);

			_photons_per_bucket_hint = 1;
			_max_recursion = (unsigned char) subdivisions;
			build_structure(true, radius, cloud_out);
		}

		void build(photon_cloud_type* photon_cloud, float radius_hint) {
			unsigned short photons_per_bucket_hint = 32;
			unsigned short max_nesting = 11;

			init_top_tree(photon_cloud, true);
			_has_sub_octrees = false;
			
			unsigned short subdivisions_hint = (int) (Log2(len/radius_hint) - 2.0f);
			unsigned short subdivisions = std::min(subdivisions_hint, max_nesting);

			_photons_per_bucket_hint = photons_per_bucket_hint;
			_max_recursion = (unsigned char) subdivisions;
			build_structure();

			AiMsgInfo("  Octree stats: %g kilo-substructures. 8^%d max divisions (optimized for radius %g). ", 
				*octree_count/1000.0f, subdivisions, radius_hint);
			unsigned long long octree_size = *octree_count * sizeof(photon_octree_type);
			unsigned long long lists_estimate = photon_cloud->size() * sizeof(size_t);
			AiMsgInfo("  Octree stats: %g mb. ", (octree_size + lists_estimate) / (1024.0f * 1024.0f));
		}


		void destroy_structure() {
			ripple_destroy();

			if (_photon_list != NULL) {
				delete _photon_list;
			}

			delete octree_count;
		}


		void ripple_destroy() {
			if (_photon_list != NULL) {
				delete _photon_list;
			}

			for (unsigned char i = 0; i < 8; i++) {
				if (_sub_octrees[i] != NULL) {
					_sub_octrees[i]->ripple_destroy();
					delete _sub_octrees[i];
				}
			}
		}
} photon_octree_type;


float blackman_harris(float distance, float radius) {
	float x = distance/radius;

	float a0 = 0.35875f;
	float a1 = 0.48829f;
	float a2 = 0.14128f;
	float a3 = 0.01168f;

	float weight  = a0 + a1*cos(1.0f * AI_PI * x) + a2*cos(2.0f * AI_PI * x) + a3*cos(4.0f * AI_PI * x);
	return weight;
}


typedef struct t_collapse_data_type {
	AtCritSec* crit_sec;
	AtUInt32 thread_ID;
	float merge_radius;
	photon_cloud_type* cloud_in;
	photon_cloud_type* cloud_out;
	AtVector bounds_n;
	AtVector bounds_p;
} t_collapse_data_type;


void photon_cloud_append( photon_cloud_type* cloud_out, photon_cloud_type* cloud_in ) {
	cloud_out->insert(cloud_out->end(), cloud_in->begin(), cloud_in->end());
}

unsigned int threaded_cull(void * data) {
	t_collapse_data_type* thread_data = static_cast<t_collapse_data_type*> (data);

	if (thread_data->cloud_in->size() == 0) {
		return 0;
	}

	photon_octree_type octree;
	photon_cloud_type reduced_cloud;
	octree.build_while_culling(thread_data->cloud_in, thread_data->merge_radius, &reduced_cloud, &thread_data->bounds_n, &thread_data->bounds_p);
	octree.destroy_structure();
	
	size_t prev_size = thread_data->cloud_in->size();
	size_t new_size = reduced_cloud.size();


	AiCritSecEnter(thread_data->crit_sec);
		AiMsgInfo("JF Photon: Finished collapse for thread-cloud %d: ", thread_data->thread_ID);
		AiMsgInfo("  T-%d: %g -> %g kilophotons", thread_data->thread_ID, prev_size/1000.0f, new_size/1000.0f);
		AiMsgInfo("  T-%d: (%d mb of photons)", thread_data->thread_ID, (new_size * sizeof(photon_type))/(1024*1024));
		photon_cloud_append( thread_data->cloud_out, &reduced_cloud);
		delete thread_data->cloud_in;
	AiCritSecLeave(thread_data->crit_sec);

	return 0;
}


struct ShaderData{
	std::string file_name;
	bool abort;
	AtArray * write_thread_clouds;
	float write_sampling_normalizer;
	bool read_cloud_owner;
	bool read_octree_owner;
	photon_cloud_type * read_cloud;
	photon_octree_type * read_cloud_octree;
	float read_radius;
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
	AiMetaDataSetInt(mds, NULL, "maya.id", 0xc2ce82);
	AiMetaDataSetStr(mds, NULL, "maya.classification", "texture/3d");
	AiMetaDataSetStr(mds, NULL, "maya.name", "jf_photon");
	AiMetaDataSetInt(mds, NULL, "maya.id", 0xc2ce82);


	AiParameterEnum("mode", m_read, enum_modes);
	AiParameterStr("file_path", "c:/photonics.[Frame].bin");
	AiParameterInt("frame", -1);
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
	if (frame == -1) {
		frame = int( AiNodeGetFlt(render_options, "frame") + 0.49f );
	}

	char frame_string[4];
 	sprintf(frame_string, "%04d", frame);
	std::string token = "[Frame]";
	if (file_path.find(token) != std::string::npos) {	
		file_path.replace(file_path.find(token),token.length(),frame_string);
	}
	data->file_name = file_path;

	if (mode == m_write) {
		data->write_thread_clouds = AiArrayAllocate(AiNodeGetInt(render_options, "threads"), 1, AI_TYPE_POINTER);
		for (AtUInt32 i = 0; i < data->write_thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = new photon_cloud_type;
			cloud->reserve(init_cloud_size);
			AiArraySetPtr(data->write_thread_clouds, i, cloud );
		}

		return;
	} 

	if (mode == m_read || mode == m_read_visualize) {
		std::string string_path = data->file_name;
		std::ifstream infile (string_path.c_str(), std::ios::binary);

		float read_radius = AiNodeGetFlt(node, "read_radius");

		data->read_cloud_owner = true;
		data->read_octree_owner = true;
		data->read_cloud = NULL;
		data->read_cloud_octree = NULL;
		data->read_radius = read_radius;

		AtNodeIterator * shader_iterator = AiUniverseGetNodeIterator(AI_NODE_SHADER);
		while (!AiNodeIteratorFinished(shader_iterator)) {
			AtNode *other_photon_shader = AiNodeIteratorGetNext(shader_iterator);
			if (AiNodeIs(other_photon_shader, "jf_photon")) {
				ShaderData *other_data = (ShaderData*)AiNodeGetLocalData(other_photon_shader);
				if (other_photon_shader == node) {
					continue;					
				}
				if (other_data == NULL) {
					continue;
				}
				if (other_data->file_name.compare(data->file_name) == 0 &&
					other_data->read_cloud_owner == true &&
					other_data->read_octree_owner == true
					) {

					AiMsgWarning("JF Photon: Cloud already in memory, reusing: %s", string_path.c_str());
					data->read_cloud_owner = false;
					data->read_cloud = other_data->read_cloud;
					
					float tolerance_factor = 4.0f;
					if (other_data->read_radius < (read_radius * tolerance_factor) &&
						other_data->read_radius > (read_radius / tolerance_factor)
						) {
						AiMsgWarning("JF Photon: Octree already exists, reusing!");
						data->read_cloud_octree = other_data->read_cloud_octree;
						data->read_octree_owner = false;
					}
					break;

				}

				
			}
		}

		if (data->read_cloud_owner == true) {

			AiMsgWarning("JF Photon: Reading from: %s ", string_path.c_str());
			if (!infile.good()) {
				data->abort = true;
				AiMsgError("JF Photon: Unable to read file! Check for invalid paths or bad permissions or something.");
				return;
			}
			
			infile.seekg (0, infile.end);
			file_int length = infile.tellg();

			if (length == 0) {
				data->abort = true;
				AiMsgWarning("JF Photon: File contains zero photons.");
				return;
			}


			file_int photon_size = sizeof(photon_type);
			file_int num_photons = length/photon_size;
			file_int mb = 1024 * 1024;
			file_int chunk_size = (256 * mb) ;
			file_int num_photons_in_chunk = (chunk_size / photon_size) + 1;
			chunk_size = photon_size * num_photons_in_chunk; //No remainders here, please. 

			data->read_cloud = new photon_cloud_type;
			data->read_cloud_owner = true;
			photon_cloud_type * read_cloud = data->read_cloud; //Alias of the cloud in data->read_cloud;

			for (file_int i = 0; i < length; i+= chunk_size) {
				file_int read_bytes = std::min(length - i, chunk_size);
				file_int read_photons = read_bytes / photon_size;

				if (read_bytes != 0) {
					infile.seekg (i, infile.beg);

					//AiMsgInfo("  %d mb chunk, %d kilophotons.", (int) read_bytes/mb, read_photons/1000);
					AiMsgInfo("  %d mb chunk, %d photons.", (int) read_bytes/mb, read_photons);
					photon_type* photon_array = new photon_type[read_photons];
					char * buffer = (char*)(photon_array);

					infile.read(buffer, read_bytes);
					read_cloud->insert(read_cloud->end(), &photon_array[0], &photon_array[read_photons]);

					delete photon_array;
				}
			}

			AiMsgInfo("  Read %d mb, %d kilophotons.", length/mb, num_photons/1000);
			// photon_type* last_photon = &photon_array[num_photons - 1];
			// AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );

			if (length % photon_size != 0) {
				AiMsgWarning("JF Photon: File size indicates a non-integer number of photons.");
				data->abort = true;
			}

			if (num_photons != read_cloud->size()) {
				AiMsgError("JF Photon: Error in photon read. %d in file, %d in memory.", num_photons, read_cloud->size());
				data->abort = true;
			} else {
				AiMsgInfo("  All photons accounted for.");
			}
			// last_photon = &read_cloud->at(num_photons - 1);
			// AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );
		}

		if (data->read_octree_owner == true) {
			const clock_t photon_process_time = clock();
			
			photon_octree_type * accel = new photon_octree_type;
			accel->build(data->read_cloud, data->read_radius);
			data->read_cloud_octree = accel;

			AiMsgInfo("JF Photon: Octree completed in %f seconds", (float(clock() - photon_process_time) /  CLOCKS_PER_SEC));
		}
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

		float merge_radius = AiNodeGetFlt(node, "write_merge_radius");
		bool merge = AiNodeGetBool(node, "write_merge_photons");
		if (merge_radius == 0.0f && merge) {
			AiMsgError("JF Photon: Zero radius, unable to merge photons. If you would like raw photon clouds please disable merging.");
			return;
		}

		for (AtUInt32 i = 0; i < data->write_thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, i));
			cloud->clear();
			cloud->reserve(init_cloud_size);
		}

		unsigned long long AA_samples = AiNodeGetInt(render_options, "AA_samples");
		unsigned long long total_samples = AA_samples * AA_samples * AiNodeGetInt(render_options, "xres") * AiNodeGetInt(render_options, "yres");
		unsigned long long expected_sampling_baseline = 16777216;
		data->write_sampling_normalizer = (float) ((double) expected_sampling_baseline / (double) total_samples);
		AiMsgWarning("JF Photon: Based on expected %d kilosamples, normalization factor is %f.", total_samples/1000, data->write_sampling_normalizer);

		return;
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
		std::string string_path = data->file_name;
		std::ofstream outfile (string_path.c_str(), std::ios::binary);

		AiMsgWarning("JF Photon: Writing to: %s ", string_path.c_str());
		if (!outfile.good()) {
			AiMsgError("JF Photon: Unable to write file! Check for invalid paths or bad permissions or something.");
			return;
		}

		float merge_radius = AiNodeGetFlt(node, "write_merge_radius");
		bool merge = AiNodeGetBool(node, "write_merge_photons");
		bool remerge = AiNodeGetBool(node, "write_remerge_photons");

		if (merge_radius == 0.0f && merge) {
			AiMsgError("JF Photon: Zero radius, unable to merge photons. If you would like raw photon clouds please disable merging.");
			return;
		}

		size_t photon_count = 0;
		size_t orig_photon_count = 0;


		if (merge) {
			AtUInt32 thread_count = data->write_thread_clouds->nelements;
			photon_cloud_type compiled_cloud;

			AiMsgWarning("JF Photon: Reducing (merging) photons from subclouds:");

			AtVector universal_bounds_n = AI_V3_ZERO;
			AtVector universal_bounds_p = AI_V3_ZERO;
			photon_octree_type octree;
			for (AtUInt32 i = 0; i < thread_count; i++) { // TO DO: remove this atrocity- (i < thread_count) && (i < thread_count);
				AtUInt32 thread_ID = i;
				photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, thread_ID));
				AtVector prev_bound = universal_bounds_n;
				octree.measure_bounds_from_cloud(cloud, &universal_bounds_n, &universal_bounds_p);
			}

			AtArray* threads = AiArrayAllocate(thread_count, 1, AI_TYPE_POINTER);
			AtArray* threads_data = AiArrayAllocate(thread_count, 1, AI_TYPE_POINTER);

			AtCritSec crit_sec;
			AiCritSecInit(&crit_sec);

			for (AtUInt32 thread_ID = 0; thread_ID < thread_count; thread_ID++) {
				t_collapse_data_type* t_collapse_data = new t_collapse_data_type;
				t_collapse_data->crit_sec = &crit_sec;
				t_collapse_data->thread_ID = thread_ID;
				t_collapse_data->merge_radius = merge_radius;
				t_collapse_data->cloud_in = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, thread_ID));
				t_collapse_data->cloud_out = &compiled_cloud;
				t_collapse_data->bounds_n = universal_bounds_n;
				t_collapse_data->bounds_p = universal_bounds_p;
				AiArraySetPtr(threads_data, thread_ID, t_collapse_data);
				
				orig_photon_count += t_collapse_data->cloud_in->size();

				void * thread_handle = AiThreadCreate(threaded_cull, t_collapse_data, AI_PRIORITY_NORMAL);
				AiArraySetPtr(threads, thread_ID, thread_handle);
			}

			for (AtUInt32 thread_ID = 0; thread_ID < thread_count; thread_ID++) {
				void * thread_handle = AiArrayGetPtr(threads, thread_ID);
				AiThreadWait(thread_handle);
				AiThreadClose(thread_handle);
				
				t_collapse_data_type * thread_data = static_cast<t_collapse_data_type*>(AiArrayGetPtr(threads_data, thread_ID));
				delete thread_data;
			}

			AiCritSecClose(&crit_sec);
			
			AiArrayDestroy(threads);
			AiArrayDestroy(threads_data);

			if (compiled_cloud.size() > 0) {
				if (remerge) {						
					AiMsgWarning("JF Photon: Rereducing:");
					photon_cloud_type cloud_out;
					photon_octree_type octree;
					octree.build_while_culling(&compiled_cloud, merge_radius, &cloud_out, &universal_bounds_n, &universal_bounds_p);
					octree.destroy_structure();

					AiMsgInfo("  %g -> %g kilophotons (rereduction)", compiled_cloud.size()/1000.0f, cloud_out.size()/1000.0f);
					// AiMsgInfo("  %d -> %d photons (rereduction)", compiled_cloud.size(), cloud_out.size());
					outfile.write((const char*)&cloud_out.at(0), (file_int) sizeof(photon_type) * (file_int) cloud_out.size());
					photon_count += cloud_out.size();
				} else {
					outfile.write((const char*)&compiled_cloud.at(0), (file_int) sizeof(photon_type) * (file_int) compiled_cloud.size());
					photon_count += compiled_cloud.size();
				}	

			}
		} else {
			// Naive Write
			AtUInt32 thread_count = data->write_thread_clouds->nelements;

			for (AtUInt32 i = 0; i < thread_count; i += 1) { 
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
		AiMsgWarning("JF Photon: Final tally: %g mb, %g kilophotons.", cloud_mb, photon_count/1000.0f);
		if (merge) {
			AiMsgInfo("  (reduced from: %g mb, %g kilophotons)", orig_cloud_mb, orig_photon_count/1000.0f);
		}

		outfile.close();
		AiArrayDestroy(data->write_thread_clouds);
	}

	if ((mode == m_read || mode == m_read_visualize)) {
		if (data->read_octree_owner) {
			AiMsgWarning("JF Photon: Destroying Octree...");
			data->read_cloud_octree->destroy_structure();
			delete data->read_cloud_octree;
		}		
		if (data->read_cloud_owner == true) {
			AiMsgWarning("JF Photon: Destroying Cloud...");
			delete data->read_cloud;
		}
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
			photon_type photon = {photon_energy * data->write_sampling_normalizer, sg->P, sg->Rt};
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
			data->read_cloud_octree->get_photons_in_radius(&photon_IDs, &sg->P, radius);

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
			data->read_cloud_octree->get_photons_in_radius(&photon_IDs, &sg->P, radius);
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
	strcpy(node->version, AI_VERSION);
	return true;
}
