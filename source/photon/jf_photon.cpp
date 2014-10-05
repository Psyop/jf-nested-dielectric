#include <ai.h>
// #include <cstring>
// #include <math.h>
// #include <stdlib.h>
// #include <iostream>
// #include <fstream>
#include <vector>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <vector>

#include <windows.h>



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
	AtArray * thread_clouds;
	photon_cloud_type * photon_cloud;
	volatile int count;
	volatile bool locked;
	volatile bool initialized;
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

 
node_initialize {
	ShaderData* data = new ShaderData;
	AiNodeSetLocalData(node,data);
	int mode = AiNodeGetInt(node, "mode");
	if (mode == m_write) {
		AtNode * renderOptions = AiUniverseGetOptions();

		AiMsgWarning("Setting up %d photon containers for the threads.", AiNodeGetInt(renderOptions, "threads"));

		data->thread_clouds = AiArrayAllocate(AiNodeGetInt(renderOptions, "threads"), 1, AI_TYPE_POINTER);
		for (AtUInt32 i = 0; i < data->thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = new photon_cloud_type;
			cloud->reserve(sizeof(photon_type) * 100000);
			AiArraySetPtr(data->thread_clouds, i, cloud );
		}
	}
}
 
node_update {
	int mode = AiNodeGetInt(node, "mode");
	if (mode == m_write) {		
		ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

		for (AtUInt32 i = 0; i < data->thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->thread_clouds, i));
			cloud->clear();
			cloud->reserve(sizeof(photon_type) * 10000);
		}

		data->locked = false;
		data->count = 0;
		data->initialized = false;

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

		AiMsgWarning("Writing Photon Cloud to: %s ", string_path.c_str());

		std::ofstream outfile (string_path, std::ios::binary);

		if (!outfile.good()) {
			AiMsgError("Unable to write file! Check for invalid paths of bad permissions or something.");
			return;
		}

		size_t real_photon_count = 0;
		for (AtUInt32 i = 0; i < data->thread_clouds->nelements; i++) {
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->thread_clouds, i));
			size_t cloud_size = cloud->size();
			if (cloud_size > 0) {				
				outfile.write((const char*)&cloud->at(0), sizeof(photon_type) * cloud->size());
				real_photon_count += cloud_size;
			}
			cloud->clear();
			delete cloud;			
		}

		float megabytes = (float) (real_photon_count * sizeof(photon_type)) / (1024.0f*1024.0f);
		AiMsgWarning("Photon cloud: %f mb, %d photons.", megabytes, real_photon_count);

		outfile.close();
		AiArrayDestroy(data->thread_clouds);

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
			photon_cloud_type * cloud = static_cast<photon_cloud_type*>(AiArrayGetPtr(data->thread_clouds, sg->tid));
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