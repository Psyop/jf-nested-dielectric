#include <ai.h>
#include <vector>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <vector>

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
	AtVector positon;
	// AtVector direction;
	AtUInt16 type;
} photon_type;


typedef std::vector<photon_type> photon_cloud_type;

struct ShaderData{
	AtArray * write_thread_clouds;
	photon_cloud_type * read_cloud;
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

		const clock_t photon_process_time = clock();

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


		// Data gets the photon array, and the number of photons
		data->read_cloud->assign(photon_array, photon_array + num_photons);
		AiMsgWarning("Photon vector contains %d photons", data->read_cloud->size());
		last_photon = &data->read_cloud->at(num_photons - 1);
		AiMsgWarning("Color of last photon %f %f %f", last_photon->energy.r, last_photon->energy.g, last_photon->energy.b );
		delete photon_array;

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