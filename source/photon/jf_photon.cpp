#include <ai.h>
#include <vector>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <vector>
#include <map>

#include <ctime>


AI_SHADER_NODE_EXPORT_METHODS(jf_photon_methods);


const char * enum_modes[] =
{
	"disabled",
	"read",
	"write",
	NULL
};

enum brdfs
{
	m_disabled,
	m_read,
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
		
		child->_photon_list = new photon_list_type;
		child->_target_photon_cloud = _target_photon_cloud;

		// if (octant != 0 && octant != 1)
		// 	AiMsgWarning("Found an octant in initialization");
	}


	void build_structure() {
		if (_photon_list->size() < _max_per_bucket || _recursion_level >= _max_recursion) {
			// Stop building. Eitehr we're divided enough, or we're at our max recursion depth. 
			if (_photon_list->size() > 0) {
				AiMsgWarning("Recursion limit. %d points, %d depth, octant %d", _photon_list->size(), _recursion_level, _octant);
			}
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

		for (size_t i = 0; i < 8; i++) {
			_sub_accells[i]->build_structure();
		}
	}
	
public:
	void build(photon_cloud_type* photon_cloud, size_t max_per_bucket = 100, size_t max_nesting = 16) {
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
		_len = std::max(std::max(dim.x, dim.y), dim.z);
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

	if (mode == m_read) {


		std::string string_path = AiNodeGetStr(node, "file_path");
		std::ifstream infile (string_path, std::ios::binary);

		AiMsgWarning("Reading Photon Cloud from: %s ", string_path.c_str());
		if (!infile.good()) {
			AiMsgError("Unable to read file! Check for invalid paths or bad permissions or something.");
			return;
		}

		infile.seekg (0, infile.end);
		size_t length = infile.tellg();
		infile.seekg (0, infile.beg);

		size_t num_photons = length/sizeof(photon_type);

		data->read_cloud = new photon_cloud_type;
		data->read_cloud->reserve(length);

		// allocate memory:
		photon_type* photon_array = new photon_type[num_photons];

		char * buffer = (char*)(photon_array);
		infile.read(buffer, length);

		AiMsgWarning("The file is %d long, length %d, photons %d", infile.gcount(), length, num_photons);
		photon_type* last_photon = &photon_array[num_photons - 1];
		AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );


		const clock_t photon_process_time = clock();

		photon_cloud_type * v_cloud = data->read_cloud; //Alias of the cloud in data->read_cloud;
		
		v_cloud->assign(photon_array, photon_array + num_photons);
		AiMsgWarning("Photon vector contains %d photons", v_cloud->size());
		last_photon = &v_cloud->at(num_photons - 1);
		AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );
		delete photon_array;


		// for(std::vector<photon_type>::iterator it = v_cloud->begin(); it != v_cloud->end(); ++it) {
		// 	// AtVector myVec = it->pos;
		//     if (abs(it->pos.x) < 0.0001) {
		//     	AiMsgWarning("Photon with attitude!");
		//     }
		// }

		photon_accellerator_type * accel = new photon_accellerator_type;
		accel->build(v_cloud);

		AiMsgWarning("Photons read and processed: %f seconds", (float( clock () - photon_process_time ) /  CLOCKS_PER_SEC));

		
		//AiMsgWarning("Lower bounds: %f, %f, %f", bounds_n.x, bounds_n.y, bounds_n.z);		
		//AiMsgWarning("Upper bounds: %f, %f, %f", bounds_p.x, bounds_p.y, bounds_p.z);
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

	if (mode == m_read) {
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

	if (mode == m_read) {
		// STUB
	}


	delete data;
}
 
shader_evaluate {
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);
	int mode = AiShaderEvalParamEnum(p_mode);
	if (mode == m_write) {
		AtColor photon_energy = AI_RGB_BLACK;
		bool is_photon = AiStateGetMsgRGB( "photon_energy", &photon_energy );

		if (is_photon) {
			photon_type photon = {photon_energy, sg->P, sg->Rt};
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->write_thread_clouds, sg->tid));
			cloud->push_back(photon);
		}
	}

	sg->out.RGBA = AI_RGBA_RED;
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