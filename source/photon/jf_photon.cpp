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
	AiMsgWarning("node init called");
	data->photon_cloud = new photon_cloud_type;
}
 
node_update {
	AtNode * renderOptions = AiUniverseGetOptions();

	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

	data->photon_cloud->clear();
	data->photon_cloud->reserve(sizeof(photon_type) * 100000);

	data->locked = false;
	data->count = 0;
	data->initialized = false;

}

node_finish {
	AiMsgWarning("node finish called");
	if (AiNodeGetLocalData(node) != NULL) {
		ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

		std::string string_path = AiNodeGetStr(node, "file_path");

		AiMsgWarning("Size is %d, lenght is %d, lenght sanity check is %d", sizeof(photon_type), data->photon_cloud->size(), data->count);
		AiMsgWarning("Writing file on node finish");
		AiMsgWarning(string_path.c_str());

		std::ofstream outfile (string_path, std::ios::binary);

		if (!outfile.good()) {
			AiMsgError("Unable to write file!");
		}

		outfile.write((const char*)&data->photon_cloud->at(0), sizeof(photon_type) * data->photon_cloud->size());
		outfile.close();
		
	}
}
 
shader_evaluate {
	ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

	AtColor photon_energy = AI_RGB_BLACK;
	bool is_photon = AiStateGetMsgRGB( "photon_energy", &photon_energy );

	if (is_photon) {
		photon_type photon = {photon_energy, sg->P, sg->Rt};

		int i = 0;
		while (data->locked){
			i++;
			Sleep(1); // TO DO: Cross platform way to do sleep
		}
		if (i > 10) {
			AiMsgWarning("Slept %d times", i); // TO DO: Better message about this
		}

		data->locked = true;
		data->initialized = true;
		data->photon_cloud->push_back(photon);
		data->count++;
		data->locked = false;
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