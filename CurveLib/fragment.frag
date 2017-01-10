#version 400

// pipeline-ból bejövő per-fragment attribútumok
in vec3 vs_out_pos;
in vec3 vs_out_normal;
in vec3 vs_out_col;

// kimenő érték - a fragment színe
out vec4 fs_out_col;

//
// uniform változók
//
uniform vec3 eye_pos = vec3(15, 15, 15);
uniform int  siPower = 32;
uniform bool is_ambient = true;
uniform bool is_diffuse = true;
uniform bool is_specular = true;

// színtér tulajdonságok
uniform vec3 light_dir = normalize(vec3(1, -1, 1));

// fénytulajdonságok
uniform vec4 La;
uniform vec4 Ld;
uniform vec4 Ls;

// anyagtulajdonságok
uniform vec4 Ka;
uniform vec4 Kd;
uniform vec4 Ks;


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

	vec3 n = normalize(vs_out_normal);
	vec3 l = -light_dir;
	float di = clamp(dot(n, l), 0.0f, 1.0f);
	vec4 diffuse = Ld * Kd * di;
	/*
	
	vec3 n = normalize(light_dir-vs_out_pos);
	float di = clamp(dot(n,  normalize( vs_out_normal )), 0.0f, 1.0f);
	vec4 diffuse = Ld * Kd * di;
	
	*/
	//
	// fényfoltképző szín
	//

	/* segítség:
		- reflect: http://www.opengl.org/sdk/docs/manglsl/xhtml/reflect.xml
		- power: http://www.opengl.org/sdk/docs/manglsl/xhtml/pow.xml
	*/
	
	vec4 specular = vec4(0, 0, 0, 1);

	if(di > 0){
		vec3 r = reflect(light_dir, n);
		vec3 c = normalize(eye_pos - vs_out_pos);
		float si = pow(clamp(dot(r, c), 0.0f, 1.0f), siPower);
		specular = Ls * Ks * si;
	}

	//
	// a fragment végső színének meghatározása
	//
	if(!is_ambient) ambient   = vec4(0, 0, 0, 1);
	if(!is_diffuse) diffuse	  = vec4(0, 0, 0, 1);
	if(!is_specular) specular = vec4(0, 0, 0, 1);

	fs_out_col =  vec4(vs_out_col,1) * (ambient + diffuse + specular);
}