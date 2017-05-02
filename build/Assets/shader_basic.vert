#version 450

// VBO-ból érkezõ változók
in vec3 vs_in_pos;

// a pipeline-ban tovább adandó értékek
out vec3 vs_out_pos;

uniform mat4 VP;

void main()
{
	gl_Position = VP*vec4( vs_in_pos, 1 );

	vs_out_pos = vs_in_pos;
}