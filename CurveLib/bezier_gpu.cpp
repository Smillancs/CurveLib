#include "Bezier.h"

#ifdef GPGPU


void BezierCurve::initOGL()
{
	// clear color set to blue-ish
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);

	glEnable(GL_CULL_FACE);		// turn on back-face culling
	glEnable(GL_DEPTH_TEST);	// enable depth-test
	
	// create a VAO
	glGenVertexArrays(1, &m_vaoID);
	// activate the new VAO m_vaoID
	glBindVertexArray(m_vaoID);

	// create a VBO
	glGenBuffers(1, &m_vboID);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboID); // activate the VBO m_vboID
											// load the data stored in array vert into the VBO (essentially: upload the data to the GPU)
	/*glBufferData(GL_ARRAY_BUFFER,	// allocate memory for the active VBO and set its data
		sizeof(vert),		// size of the VBO allocation, in bytes
		vert,				// load data into the VBO from this location of the system memory
		GL_STATIC_DRAW);	// we only want to store data into the VBO once (STATIC), and we want to use the VBO as a source for drawing our scene at each frame (DRAW)
							// for other usage flags see http://www.opengl.org/sdk/docs/man/xhtml/glBufferData.xml
							// {STREAM, STATIC, DYNAMIC} and {DRAW, READ, COPY}


							// activate the first general attribute 'channel' in the VAO
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(
		0,				// set the attributes of VAO channel 0
		3,				// this channel has 3 componenets
		GL_FLOAT,		// each of those componenets are floats
		GL_FALSE,		// do not normalize
		sizeof(Vertex),	// stride
		0				// channel 0`s data begins at the beginning of the VBO, no offset
	);

	// activate 'channel' idx 1
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(
		1,
		3,
		GL_FLOAT,
		GL_FALSE,
		sizeof(Vertex),
		(void*)(sizeof(glm::vec3)));

	// texture coordinates
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(
		2,
		2,
		GL_FLOAT,
		GL_FALSE,
		sizeof(Vertex),
		(void*)(2*sizeof(glm::vec3)));
		*/
	// create index buffer
/*	glGenBuffers(1, &m_ibID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
	*/
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	
	//
	// shader initialization
	//
	GLuint vs_ID = loadShader(GL_VERTEX_SHADER, "vertex.vert");
	GLuint fs_ID = loadShader(GL_FRAGMENT_SHADER, "fragment.frag");

	// create the shader container (program)
	shaders.push_back(glCreateProgram());

	// attach the vertex and fragment (pixel) shaders to the program
	glAttachShader(shaders[0], vs_ID);
	glAttachShader(shaders[0], fs_ID);

	// make correspondances between the VAO channels and the shader 'in' variables
	// IMPORTANT: do this prior to linking the programs!
	glBindAttribLocation(shaders[0],	// ID of the shader program from which we want to map a variable to a channel
		0,				// the VAO channel number we want to bind the variable to
		"vs_in_pos");	// the name of the variable in the shader
	glBindAttribLocation(shaders[0], 1, "vs_in_normal");
	glBindAttribLocation(shaders[0], 2, "vs_in_tex0");

	// link the shaders
	glLinkProgram(shaders[0]);

	// check the linking
	GLint infoLogLength = 0, result = 0;

	glGetProgramiv(shaders[0], GL_LINK_STATUS, &result);
	glGetProgramiv(shaders[0], GL_INFO_LOG_LENGTH, &infoLogLength);
	if (GL_FALSE == result)
	{
		std::vector<char> ProgramErrorMessage(infoLogLength);
		glGetProgramInfoLog(shaders[0], infoLogLength, NULL, &ProgramErrorMessage[0]);
		fprintf(stdout, "%s\n", &ProgramErrorMessage[0]);

		char* aSzoveg = new char[ProgramErrorMessage.size()];
		memcpy(aSzoveg, &ProgramErrorMessage[0], ProgramErrorMessage.size());

		std::cout << "[app.Init()] Sáder Huba panasza: " << aSzoveg << std::endl;

		delete aSzoveg;
	}

	// we can dispose of the vertex and fragment shaders
	glDeleteShader(vs_ID);
	glDeleteShader(fs_ID);

}

BezierCurve::BezierCurve() 
	: mainCoeff(1.0), cp({}), coeff({})
{
}

BezierCurve::BezierCurve(const std::vector<Eigen::Vector2f> &_cp)
	: cp(_cp), mainCoeff(1.0)
{
	RegenerateCoeff();
}


BezierCurve BezierCurve::operator=(const BezierCurve & _other)
{
	cp = _other.GetControlPoints();
	RegenerateCoeff();
	return *this;
}

void BezierCurve::RegenerateCoeff()
{
	coeff.clear();
	size_t n = cp.size();
	for (size_t i = 0; i < n; ++i)
		coeff.push_back(binomial(n, i));
	dirtyCoeff = false;
}

void BezierCurve::addControlPoint(const Eigen::Vector2f _cp)
{
	cp.push_back(_cp);
	dirtyCoeff = true;
}

std::vector<Eigen::Vector2f> BezierCurve::deCasteljauEval(const float t, const size_t deg)
{
	std::vector<Eigen::Vector2f> tempCp = cp;

	for (size_t j = 0; j < deg; ++j)
	{
		for (size_t i = 0; i < tempCp.size()-1; ++i)
		{
			tempCp[i] = (1.0f-t)*tempCp[i]+t*tempCp[i+1];
		}
		tempCp.pop_back();
	}

	return tempCp;
}

Eigen::Vector2f BezierCurve::BernsteinEval(const float t)
{
	Eigen::Vector2f temp = Eigen::Vector2f(0,0);
	size_t n = cp.size();
	for(size_t i = 0; i < n; ++i)
	{
		temp += coeff[i]*pow(t,i)*pow(1-t,n-i)*cp[i];
	}
	return mainCoeff*temp;
}

//todo cachelni az evalhoz
BezierCurve BezierCurve::Diff(const size_t order)
{
	BezierCurve hodo;
	size_t n = cp.size()-order;
	for(size_t i = 0; i < n; ++i)
	{
		Eigen::Vector2f lcp = {0,0};
		for (size_t j = 0; j < order; ++j)
		{
			lcp += binomial(order,j)*pow(-1,order-j)*cp[i+j];
		}
		hodo.addControlPoint(lcp);
	}
	hodo.RegenerateCoeff();
	hodo.mainCoeff = factorial(cp.size())/factorial(n);
	return hodo;
}
/* does it even make sense?
double BezierCurve::DiffEval(const double t)
{
	BezierCurve temp = Diff(1);
	return BernsteinEval(temp(t));
}
*/

BezierCurve BezierCurve::Subdivision()
{
	return BezierCurve();
}

BezierCurve BezierCurve::Elevation()
{
	return BezierCurve();
}

BezierCurve BezierCurve::Reduction()
{
	return BezierCurve();
}


void BezierCurve::ToExplicit()
{

}


void BezierCurve::ToParametric()
{

}

#endif