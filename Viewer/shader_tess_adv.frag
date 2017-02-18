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

void main()
{
	fs_out_col = vec4(In.k, 5*In.t, -5*In.t, 1);
}
