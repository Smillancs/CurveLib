#version 450

// VBO-b�l �rkez� v�ltoz�k
in vec3 vs_in_pos;

// a pipeline-ban tov�bb adand� �rt�kek
out vec3 vs_out_pos;

uniform mat4 VP;

void main()
{
	gl_Position = VP*vec4( vs_in_pos, 1 );

	vs_out_pos = vs_in_pos;
}