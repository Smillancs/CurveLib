#version 130

// VBO-ból érkezõ változók
in vec3 vs_in_pos;
in vec3 vs_in_e;
in vec3 vs_in_n;
in vec3 vs_in_b;
in float vs_in_k;
in float vs_in_t;

// a pipeline-ban tovább adandó értékek
out vec3 vs_out_pos;
out vec3 vs_out_e;
out vec3 vs_out_n;
out vec3 vs_out_b;
out float vs_out_k;
out float vs_out_t;


// shader külsõ paraméterei - most a három transzformációs mátrixot külön-külön vesszük át
uniform mat4 world;
uniform mat4 worldIT;
uniform mat4 MVP;

void main()
{
	gl_Position = MVP * vec4( vs_in_pos, 1 );

	vs_out_pos = (world * vec4( vs_in_pos, 1 )).xyz;
	vs_out_e = (world * vec4( vs_in_e, 0 )).xyz;
	vs_out_n = (world * vec4( vs_in_n, 0 )).xyz;
	vs_out_b = (world * vec4( vs_in_b, 0 )).xyz;
	vs_out_k = vs_in_k;
	vs_out_t = vs_in_t;
}