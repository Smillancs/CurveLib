#version 420

in vec3	vs_in_pos;
in vec3	vs_in_e;
in vec3	vs_in_n;
in vec3	vs_in_b;
in float	vs_in_k;
in float	vs_in_t;

out block
{
	vec3	pos;
	vec3	e;
	vec3	n;
	vec3	b;
	float	k;
	float	t;
} Out;

void main()
{
	Out.pos		= vs_in_pos;
	Out.e		= vs_in_e;
	Out.n		= vs_in_n;
	Out.b		= vs_in_b;
	Out.k		= vs_in_k;
	Out.t		= vs_in_t;
}
