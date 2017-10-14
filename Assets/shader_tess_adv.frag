#version 420
#define max_length 128

in block
{
	vec3	pos;
	vec3	col;
	vec3	patch_coords;
	float	k;
  float dK;
	float	t;
  float v;
} In;

out vec4 fs_out_col;

uniform vec3 pointData[max_length];

uniform int redRole = 0;
uniform int greenRole = 0;
uniform int blueRole = 0;

uniform float scaleRed = 1;
uniform float scaleGreen = 1;
uniform float scaleBlue = 1;

void main()
{
	int MAX_STEP = 32;
	bool found = false;

	for(int i=0; i<MAX_STEP && !found;++i)
	{
	}

	vec2 pos = gl_FragCoord.xy;
  vec3 color = vec3(0);
  if(redRole == 1)
    color.x = In.v * scaleRed;
  if(redRole == 2)
    color.x = In.k * scaleRed;
  if(redRole == 3)
    color.x = In.dK * scaleRed;
  if(redRole == 4)
    color.x = In.t * scaleRed;
  if(redRole == 5)
    color.x = -In.t * scaleRed;

  if(greenRole == 1)
    color.y = In.v * scaleGreen;
  if(greenRole == 2)
    color.y = In.k * scaleGreen;
  if(greenRole == 3)
    color.y = In.dK * scaleGreen;
  if(greenRole == 4)
    color.y = In.t * scaleGreen;
  if(greenRole == 5)
    color.y = -In.t * scaleGreen;

  if(blueRole == 1)
    color.z = In.v * scaleBlue;
  if(blueRole == 2)
    color.z = In.k * scaleBlue;
  if(blueRole == 3)
    color.z = In.dK * scaleBlue;
  if(blueRole == 4)
    color.z = In.t * scaleBlue;
  if(blueRole == 5)
    color.z = -In.t * scaleBlue;

  //color.x = colorRoles[0];
  fs_out_col = vec4(color.x, color.y, color.z, 1);
	//fs_out_col = vec4(In.k, 5*In.t, -5*In.t, 1);
}
