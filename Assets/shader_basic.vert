#version 130

// VBO-ból érkezõ változók
in vec3 vs_in_pos;
in vec3 vs_in_col;

// a pipeline-ban tovább adandó értékek
out vec3 vs_out_pos;
out vec3 vs_out_col;


// shader külsõ paraméterei - most a három transzformációs mátrixot külön-külön vesszük át
uniform mat4 VP;

void main()
{
	gl_Position = VP * vec4( vs_in_pos, 1 );

	vs_out_pos = vs_in_pos;
	vs_out_col = vs_in_col;
}
