#version 420

in vec3	vs_in_pos;

out block
{
	vec3	pos;
} Out;

void main()
{
	Out.pos		= vs_in_pos;
}
