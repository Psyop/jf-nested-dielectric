#include <ai.h>
#include <cstring>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

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
}
 
node_update {
}
 
node_finish {
}
 
shader_evaluate {
	sg->out.RGBA = AI_RGBA_BLACK;
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