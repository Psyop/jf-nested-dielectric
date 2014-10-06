#include <ai.h>


#include <algorithm>
#include <fstream>
#include <iterator>
#include <vector>
#include <map>

#include <ctime>
#include <Windows.h>


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
} photon_type;

typedef std::vector<photon_type> photon_cloud_type; // vector of photons with all the information a photon has

typedef std::vector<size_t> photon_list_type; // a list of photons in a photon cloud


typedef class photon_accellerator_type{
	AtVector bounds_n;
	AtVector bounds_p;

	// AtVector _pos;
	float _len;

	bool _has_sub_accells;

	size_t _recursion_level;
	size_t _max_recursion;
	size_t _max_per_bucket;
	size_t _octant;

	photon_list_type * _photon_list;
	photon_accellerator_type * _sub_accells[8];

	photon_cloud_type * _target_photon_cloud;

	void init() {
		_photon_list = new photon_list_type;
	}

	void set_bounds(size_t octant, AtVector par_bounds_n, float par_len) {
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
		bounds_p = bounds_n + offset_vector;
		_len = half_len;

	}

	bool within_bounds(AtVector photon_pos) {
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

	bool within_range(AtVector photon_pos, float radius) {
		if (photon_pos.x + radius >= bounds_n.x && 
			photon_pos.y + radius >= bounds_n.y && 
			photon_pos.z + radius >= bounds_n.z &&
			photon_pos.x - radius < bounds_p.x &&
			photon_pos.y - radius < bounds_p.y &&
			photon_pos.z - radius < bounds_p.z
			) {
			return true;
		}
		return false;
	}

	void add_ID_to_bucket(size_t ID) {
		_photon_list->push_back(ID);
	}

	void initialize_child(photon_accellerator_type * child, size_t octant) {
		child->init();
		child->set_bounds(octant, bounds_n, _len); //will set _pos, bounds_n and bounds_p

		child->_octant = octant;

		child->_has_sub_accells = false;
		child->_recursion_level = _recursion_level + 1;
		child->_max_recursion = _max_recursion;
		child->_max_per_bucket = _max_per_bucket;
		
		child->_target_photon_cloud = _target_photon_cloud;

		child->_sub_accells[0] = NULL;
		child->_sub_accells[1] = NULL;
		child->_sub_accells[2] = NULL;
		child->_sub_accells[3] = NULL;
		child->_sub_accells[4] = NULL;
		child->_sub_accells[5] = NULL;
		child->_sub_accells[6] = NULL;
		child->_sub_accells[7] = NULL;

		// if (octant != 0 && octant != 1)
		// 	AiMsgWarning("Found an octant in initialization");
	}


	void build_structure() {
		if (_photon_list->size() < _max_per_bucket || _recursion_level >= _max_recursion) {
			// Stop building. Either we're divided enough, or we're at our max recursion depth. 
			// if (_photon_list->size() > 0) {
			// 	AiMsgWarning("Recursion limit. %d points, %d depth, octant %d", _photon_list->size(), _recursion_level, _octant);
			// }
			return;
		}
		// make 8 children and initialize them
		for (size_t i = 0; i < 8; i++) {
			_sub_accells[i] = new photon_accellerator_type;
			initialize_child(_sub_accells[i], i);
		}

		// distribute points
		for(size_t i = 0; i != _photon_list->size(); i++) {
			size_t photon_ID = _photon_list->at(i);
			for (size_t i = 0; i < 8; i++) {
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
		for (size_t i = 0; i < 8; i++) {
			_sub_accells[i]->build_structure();
		}
	}
	

public:
	photon_list_type get_photons(AtVector pos, float radius) {
		photon_list_type found_photons;

		if (!_has_sub_accells) {
			return *_photon_list;
			// return _photon_list->size();
		}

		photon_accellerator_type * s_sub_accell = _sub_accells[0];
		std::vector<photon_accellerator_type *> accells_to_search;

		for (size_t i = 0; i < 8; i++) {
			if( _sub_accells[i]->within_range(pos, radius)) {
				accells_to_search.push_back(_sub_accells[i]);
			}
		}

		//bool exhausted = false;
		size_t photon_count = 0;
		size_t boxes_count = 0;
		size_t max_depth = 0;
		while (accells_to_search.size() > 0) {
			std::vector<photon_accellerator_type *> s_accells = accells_to_search;
			accells_to_search.clear();
			for (size_t i = 0; i < s_accells.size(); i++) {
				if (s_accells[i]->_has_sub_accells) {
					// no points here, add sub accels to accells_to_search
					for (size_t g = 0; g < 8; g++) {
						if (s_accells[i]->_sub_accells[g]->within_range(pos, radius)) {
							accells_to_search.push_back(s_accells[i]->_sub_accells[g]);
						}
					}
				} else {
					boxes_count++;
					photon_count += s_accells[i]->_photon_list->size();
					max_depth = max( s_accells[i]->_recursion_level, max_depth);

					found_photons.insert(found_photons.end(), s_accells[i]->_photon_list->begin(), s_accells[i]->_photon_list->end());
				}

			}
		}
		// AiMsgWarning("Found %d photons in %d buckets, max %d", photon_count, boxes_count, max_depth);
		return found_photons;		
	}

	photon_list_type get_photons_in_radius(AtVector pos, float radius) {
		photon_list_type found_photons = get_photons(pos, radius);
		// AiMsgWarning("Got back %d photons", found_photons.size());
		photon_list_type filtered_photons;
		for (size_t i = 0; i < found_photons.size(); i++) {
			photon_type* photon = &_target_photon_cloud->at(found_photons[i]);
			if (AiV3Dist(pos, photon->pos) < radius) { 
				filtered_photons.push_back(found_photons[i]);
			}
		}
		return filtered_photons;
		// return found_photons;
	}

	void build(photon_cloud_type* photon_cloud, size_t max_per_bucket = 20, size_t max_nesting = 16) {
		_target_photon_cloud = photon_cloud;
		_max_per_bucket = max_per_bucket;
		_max_recursion = max_nesting;
		_recursion_level = 0;
		_octant = -1;
		init();

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

		// _pos = measured_bounds_n;
		_len = max(max(dim.x, dim.y), dim.z);
		bounds_n = measured_bounds_n;
		bounds_p = measured_bounds_p;

		build_structure();
	}
	void destroy_structure() {
		for (size_t i = 0; i < 8; i++) {
			if (_sub_accells[i] != NULL) {
				_sub_accells[i]->destroy_structure();
			}
		}
		if (_photon_list != NULL) {
			delete _photon_list;		
		}
	}
} photon_accellerator_type;






struct ShaderData{
	AtArray * write_thread_clouds;
	photon_cloud_type * read_cloud;
	photon_accellerator_type * read_cloud_accelerator;
};



enum jf_photonParams {
	p_mode,
	p_file_path,
	p_read_radius,
};
 
node_parameters {
	AiParameterEnum("mode", m_read, enum_modes);
	AiParameterStr("file_path", "");
	AiParameterFlt("read_radius", 0.1f);
}

size_t init_cloud_size = sizeof(photon_type) * 100000;
 
node_initialize {
	ShaderData* data = new ShaderData;
	AiNodeSetLocalData(node,data);

	AtNode * renderOptions = AiUniverseGetOptions();

	int mode = AiNodeGetInt(node, "mode");
	
	if (mode == m_write) {
		// AiMsgWarning("Setting up %d photon containers for the threads.", AiNodeGetInt(renderOptions, "threads"));
		data->write_thread_clouds = AiArrayAllocate(AiNodeGetInt(renderOptions, "threads"), 1, AI_TYPE_POINTER);
		for (AtUInt32 i = 0; i < data->write_thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = new photon_cloud_type;
			cloud->reserve(init_cloud_size);
			AiArraySetPtr(data->write_thread_clouds, i, cloud );
		}
		return;
	} 

	if (mode == m_read || mode == m_read_visualize) {

		
		// Sleep(2000);

		std::string string_path = AiNodeGetStr(node, "file_path");
		std::ifstream infile (string_path, std::ios::binary);

		AiMsgWarning("Reading Photon Cloud from: %s ", string_path.c_str());
		if (!infile.good()) {
			AiMsgError("Unable to read file! Check for invalid paths or bad permissions or something.");
			return;
		}
		

		infile.seekg (0, infile.end);
		size_t length = infile.tellg();

		size_t max_len = 1024*1024*512;
		length = min(length, max_len);

		infile.seekg (0, infile.beg);

		size_t num_photons = length/sizeof(photon_type);

		data->read_cloud = new photon_cloud_type;
		data->read_cloud->reserve(length);

		// allocate memory:
		photon_type* photon_array = new photon_type[num_photons];
		

		char * buffer = (char*)(photon_array);
		infile.read(buffer, length);

		int mb = 1024 * 1024;
		AiMsgWarning("The file is %d long, length %d, photons %d", infile.gcount()/mb, length/mb, num_photons);
		// photon_type* last_photon = &photon_array[num_photons - 1];
		// AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );

		const clock_t photon_process_time = clock();

		photon_cloud_type * v_cloud = data->read_cloud; //Alias of the cloud in data->read_cloud;
		
		v_cloud->assign(photon_array, photon_array + num_photons);
		AiMsgWarning("Photon vector contains %d photons", v_cloud->size());
		//last_photon = &v_cloud->at(num_photons - 1);
		// AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );
		delete photon_array;

		photon_accellerator_type * accel = new photon_accellerator_type;
		accel->build(v_cloud);
		data->read_cloud_accelerator = accel;

		AiMsgWarning("Photons read and processed: %f seconds", (float( clock () - photon_process_time ) /  CLOCKS_PER_SEC));

	}
}
 
node_update {
	int mode = AiNodeGetInt(node, "mode");

	if (mode == m_write) {		
		ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

		for (AtUInt32 i = 0; i < data->write_thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, i));
			cloud->clear();
			cloud->reserve(init_cloud_size);
		}
		return;
	}

	if (mode == m_read || mode == m_read_visualize) {
		// STUB
	}
}

node_finish {
	if (AiNodeGetLocalData(node) == NULL) {
		return;
	}

	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
	int mode = AiNodeGetInt(node, "mode");

	if (mode == m_write) {
		std::string string_path = AiNodeGetStr(node, "file_path");
		std::ofstream outfile (string_path, std::ios::binary);

		AiMsgWarning("Writing Photon Cloud to: %s ", string_path.c_str());
		if (!outfile.good()) {
			AiMsgError("Unable to write file! Check for invalid paths or bad permissions or something.");
			return;
		}

		size_t photon_count = 0;
		for (AtUInt32 i = 0; i < data->write_thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, i));

			if (cloud->size() > 0) {				
				outfile.write((const char*)&cloud->at(0), sizeof(photon_type) * cloud->size());
				photon_count += cloud->size();
			}

			cloud->clear();
			delete cloud;			
		}

		float cloud_mb = (float) (photon_count * sizeof(photon_type)) / (1024.0f*1024.0f);
		AiMsgWarning("Photon cloud: %f mb, %d photons.", cloud_mb, photon_count);

		outfile.close();
		AiArrayDestroy(data->write_thread_clouds);
	}

	if ((mode == m_read || mode == m_read_visualize)) {
		AiMsgWarning("Deleting photons.");
		data->read_cloud_accelerator->destroy_structure();
		delete data->read_cloud_accelerator;
		delete data->read_cloud;
	}


	delete data;
}
 
shader_evaluate {
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
	int mode = AiShaderEvalParamEnum(p_mode);
	sg->out.RGBA = AI_RGBA_BLACK;
	if (mode == m_write) {
		AtColor photon_energy = AI_RGB_BLACK;
		bool is_photon = AiStateGetMsgRGB( "photon_energy", &photon_energy );

		if (is_photon) {
			photon_type photon = {photon_energy, sg->P, sg->Rt};
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, sg->tid));
			cloud->push_back(photon);
		}
		sg->out.RGBA = AI_RGBA_RED;
	}
	

	if (mode == m_read || mode == m_read_visualize) { 
		float radius = AiShaderEvalParamFlt(p_read_radius);



		if (mode == m_read) {
			photon_list_type photon_IDs = data->read_cloud_accelerator->get_photons_in_radius(sg->P, radius);

			AtColor energy = AI_RGB_BLACK;
			for(size_t i = 0; i != photon_IDs.size(); i++) {
				energy += data->read_cloud->at(photon_IDs[i]).energy;
			}

			sg->out.RGB = energy;
			sg->out.RGBA.a = 1.0f;
		} else if (mode == m_read_visualize) {
			photon_list_type photon_IDs = data->read_cloud_accelerator->get_photons(sg->P, radius);
			float out_value = ((float) photon_IDs.size() / 300.0f);
			sg->out.RGBA = AI_RGBA_WHITE * out_value;
			sg->out.RGBA.a = 1.0f;

		}
		
	}

}
 
node_loader {
	if (i > 0)
		return false;

	node->methods		   = jf_photon_methods;
	node->output_type  = AI_TYPE_RGBA;
	node->name					  = "jf_photon";
	node->node_type  = AI_NODE_SHADER;
	strcpy_s(node->version, AI_VERSION);
	return true;
}