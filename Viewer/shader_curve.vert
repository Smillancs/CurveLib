#version 130

// VBO-b�l �rkez� v�ltoz�k
in vec3 vs_in_pos;
in vec3 vs_in_e;
in vec3 vs_in_n;
in vec3 vs_in_b;
in float vs_in_k;
in float vs_in_t;

// a pipeline-ban tov�bb adand� �rt�kek
out vec3 vs_out_pos;
out vec3 vs_out_e;
out vec3 vs_out_n;
out vec3 vs_out_b;
out float vs_out_k;
out float vs_out_t;


// shader k�ls� param�terei - most a h�rom transzform�ci�s m�trixot k�l�n-k�l�n vessz�k �t
uniform mat4 MVP;

void main()
{
	gl_Position = MVP * vec4( vs_in_pos, 1 );

	vs_out_pos = vs_in_pos;
	vs_out_e = vs_in_e;
	vs_out_n = vs_in_n;
	vs_out_b = vs_in_b;
	vs_out_k = vs_in_k;
	vs_out_t = vs_in_t;
}
