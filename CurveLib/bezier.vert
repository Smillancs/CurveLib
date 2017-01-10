#version 420

// VBO-ból érkezõ változók
in vec3 vs_in_pos;

// a pipeline-ban tovább adandó értékek
out block
{
	vec3 pos;
} Out;

uniform mat4 world;

void main()
{
	gl_Position = world*vec4( vs_in_pos, 1 );

	Out.pos = gl_Position.xyz;
}