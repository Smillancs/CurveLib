#version 130

// bemeneti változó - most a vertex shader-bõl (vagyis ottani out)
in vec4 vs_out_col;

// kimeneti változó - a fragment színe
out vec4 fs_out_col;

void main()
{
	fs_out_col = vs_out_col;
}
