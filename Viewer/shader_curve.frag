#version 130

// pipeline-b�l bej�v� per-fragment attrib�tumok
in vec3 vs_out_pos;
in vec3 vs_out_e;
in vec3 vs_out_n;
in vec3 vs_out_b;
in float vs_out_k;
in float vs_out_t;

// kimen� �rt�k - a fragment sz�ne
out vec4 fs_out_col;

//
// uniform v�ltoz�k
//

void main()
{
	fs_out_col = vec4(vs_out_k, 5*vs_out_t, -5*vs_out_t, 1);
}
