#version 130

// pipeline-ból bejövõ per-fragment attribútumok
in vec3 vs_out_pos;
in vec3 vs_out_col;

// kimenõ érték - a fragment színe
out vec4 fs_out_col;

//
// uniform változók
//

void main()
{
	fs_out_col = vec4(vs_out_col, 1);
}
