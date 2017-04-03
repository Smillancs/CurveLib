#version 420

in block
{
	vec3	pos;
	vec3	col;
	vec3	patch_coords;
	float	k;
	float	t;
} In;

out vec4 fs_out_col;

uniform vec3 pointData[max_length];

void main()
{
	int MAX_STEP = 32;
	bool found = false;
	
	for(int i=0; i<MAX_STEP && !found;++i)
	{
	}

	vec2 pos = gl_FragCoord;
	fs_out_col = vec4(In.k, 5*In.t, -5*In.t, 1);
}
