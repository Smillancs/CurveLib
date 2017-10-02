#version 130

// lokális változók: két tömb
uniform vec3 positions[4] = vec3[4](
	vec3(-1, -1, 0),
	vec3( 1, -1, 0),
	vec3(-1,  1, 0),
	vec3( 1,  1, 0)
);

uniform vec4 color = vec4(1,0,0,1);

uniform mat4 VP;

// kimeneti változó: a csúcspont színe
out vec4 vs_out_col;

void main()
{
	// gl_VertexID: https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/gl_VertexID.xhtml
	gl_Position = VP*vec4(positions[gl_VertexID],1);
	vs_out_col	= color;
}
