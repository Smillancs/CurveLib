#version 420

// pipeline-ból bejövõ per-fragment attribútumok
in block
{
	vec3 pos;
	vec3 n;
	vec2 tex;
} In;

// kimenõ érték - a fragment színe
out vec4 fs_out_col;

//
// uniform változók
//

// színtér tulajdonságok
uniform vec3 eye_pos = vec3(0, 15, 15);

// fénytulajdonságok
uniform vec3 light_pos = vec3( 0, 5, 0 );
uniform vec4 La = vec4(0.1f, 0.1f, 0.1f, 1);
uniform vec4 Ld = vec4(0.5f, 0.5f, 0.5f, 1);
uniform vec4 Ls = vec4(1, 1, 1, 1);

// anyagtulajdonságok
uniform vec4 Ka = vec4(1, 1, 1, 1);
uniform vec4 Kd = vec4(0.75f, 0.25f, 0.125f, 1);
uniform vec4 Ks = vec4(0, 1, 0, 1);
uniform float specular_power = 16;
uniform sampler2D texImage;

void main()
{
	//
	// ambiens szín számítása
	//
	vec4 ambient = La * Ka;

	//
	// diffúz szín számítása
	//

	/* segítség:
		- normalizálás: http://www.opengl.org/sdk/docs/manglsl/xhtml/normalize.xml
	    - skaláris szorzat: http://www.opengl.org/sdk/docs/manglsl/xhtml/dot.xml
	    - clamp: http://www.opengl.org/sdk/docs/manglsl/xhtml/clamp.xml
	*/
	vec3 normal = normalize( In.n );
	vec3 toLight = normalize(light_pos - In.pos);
	float di = clamp( dot( toLight, normal), 0.0f, 1.0f );
	vec4 diffuse = Ld*Kd*di;

	//
	// fényfoltképzõ szín
	//

	/* segítség:
		- reflect: http://www.opengl.org/sdk/docs/manglsl/xhtml/reflect.xml
		- power: http://www.opengl.org/sdk/docs/manglsl/xhtml/pow.xml
	*/
	float si = 0;
	if ( di > 0 )
	{
		vec3 e = normalize( eye_pos - In.pos );
		vec3 r = reflect( -toLight, normal );
		si = pow( clamp( dot(e, r), 0.0f, 1.0f ), specular_power );
	}
	vec4 specular = Ls*Ks*si;

	fs_out_col = (ambient + diffuse + specular ) * texture(texImage, In.tex.st);
}