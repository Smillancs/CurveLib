#version 420
#define max_length 128

in block
{
	vec3	pos;
	vec3	col;
	vec3	patch_coords;
	float	k;
  float dK;
	float	t;
  float v;
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

	vec2 pos = gl_FragCoord.xy;
	//fs_out_col = vec4(In.k, 5*In.t, -5*In.t, 1);
  fs_out_col = vec4(In.dK,0,0,1);
}
