#version 130

// pipeline-ból bejövõ per-fragment attribútumok
in vec3 vs_out_pos;
in vec3 vs_out_e;
in vec3 vs_out_n;
in vec3 vs_out_b;
in float vs_out_k;
in float vs_out_t;

// kimenõ érték - a fragment színe
out vec4 fs_out_col;

//
// uniform változók
//

void main()
{
	fs_out_col = vec4(vs_out_k, 5*vs_out_t, -5*vs_out_t, 1);
}
