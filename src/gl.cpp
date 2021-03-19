#include "../include/gl.hpp"

// -------------- FRUSTUM -------------- \\

GL_Frustum GL::Frustum::extractPlanes(const GL_Matrix4& _m){  

    GL_Frustum out;

    const __m128 _1x = _mm_load_ps(&_m(0, 0));
    const __m128 _2x = _mm_load_ps(&_m(1, 0));
    const __m128 _3x = _mm_load_ps(&_m(2, 0));
    const __m128 _4x = _mm_load_ps(&_m(3, 0));
    
	util::array4f16a v;
    //plane 0
    {
        const __m128 p = _mm_add_ps(_4x, _1x);		
		_mm_store_ps(v.v, p);
        std::memcpy(out.data(), v.v, sizeof(float) * 4);
    }

    //plane 1
    {
        const __m128 p = _mm_sub_ps(_4x, _1x);
		_mm_store_ps(v.v, p);
        std::memcpy(out.data() + 4, v.v, sizeof(float) * 4);
    }

    //plane 2
    {
        const __m128 p = _mm_sub_ps(_4x, _2x);
		_mm_store_ps(v.v, p);
        std::memcpy(out.data() + 2*4, v.v, sizeof(float) * 4);
    }

    //plane 3
    {
        const __m128 p = _mm_add_ps(_4x, _2x);
		_mm_store_ps(v.v, p);
        std::memcpy(out.data() + 3*4, v.v, sizeof(float) * 4);
    }

    //plane 4
    {
        const __m128 p = _mm_add_ps(_4x, _3x);
		_mm_store_ps(v.v, p);
        std::memcpy(out.data() + 4*4, v.v, sizeof(float) * 4);
    }

    //plane 5
    {
        const __m128 p = _mm_sub_ps(_4x, _3x);
		_mm_store_ps(v.v, p);
        std::memcpy(out.data() + 5*4, v.v, sizeof(float) * 4);
    }

    out.rowwise().normalize(); //TODO necessary??
    return out;
}

bool GL::Frustum::isPointInside(const GL_Frustum& _planes, const Eigen::Vector3f& _point){
    for(size_t i = 0; i < 6; ++i){
        const float d = _planes.row(i).head(3).dot(_point) + _planes(i, 3);
        if(d <= 0.f)
            return false;
    }
    return true;
}

// -------------- SHADERPROGRAM -------------- \\

void GL::ShaderProgram::print(std::string _id, ShaderProgram::Status _compComp, ShaderProgram::Status _compVert,
	ShaderProgram::Status _compGeom, ShaderProgram::Status _compFrag, ShaderProgram::Status _link, std::string _errorLog) const{
	if (!printDebug) return;
	std::stringstream s;
	s << "\n   Shader: " << _id << std::endl
	<< "Compiling: "
	<< (_compComp == Status::failed ? " X |" : _compComp == Status::success ? " S |" : " - |")
	<< (_compVert == Status::failed ? " X |" : _compVert == Status::success ? " S |" : " - |")
	<< (_compGeom == Status::failed ? " X |" : _compGeom == Status::success ? " S |" : " - |")
	<< (_compFrag == Status::failed ? " X |" : _compFrag == Status::success ? " S |" : " - |") << std::endl
	<< "  Linking: " + std::string(_link == Status::failed ? "Failed!" : _link == Status::success ? "Success!" : " - ") << std::endl
	<< "Error Log: " << (_errorLog.empty() ? "empty" : _errorLog) << std::endl;

	spdlog::debug(s.str());

}

bool GL::ShaderProgram::compileFromFile(const std::string& _path) {
	bool cExists = true;
	bool vExists = true;
	bool gExists = true;
	bool fExists = true;

	std::stringstream p;
	p << std::filesystem::current_path().string() << "/../../" << _path;
	const std::string path = p.str();

	std::ifstream compB(path + ".comp");
	cExists = compB.good();

	std::ifstream vertB(path + ".vert");
	vExists = vertB.good();

	std::ifstream geomB(path + ".geom");
	gExists = geomB.good();

	std::ifstream fragB(path + ".frag");
	fExists = fragB.good();

	if(!cExists && !vExists && !gExists && !fExists){
		spdlog::error("no valid path for shader: {}", path);
		return false;
	}

	bool success = compile(
		(cExists ? std::string{ std::istreambuf_iterator<char>(compB), std::istreambuf_iterator<char>() } : "").c_str(),
		(vExists ? std::string{ std::istreambuf_iterator<char>(vertB), std::istreambuf_iterator<char>() } : "").c_str(),
		(gExists ? std::string{ std::istreambuf_iterator<char>(geomB), std::istreambuf_iterator<char>() } : "").c_str(),
		(fExists ? std::string{ std::istreambuf_iterator<char>(fragB), std::istreambuf_iterator<char>() } : "").c_str());

	compB.close();
	vertB.close();
	geomB.close();
	fragB.close();

	return success;
}

bool GL::ShaderProgram::compile(const char* _compute, const char* _vertex, const char* _geom, const char* _frag) {
	Status compStatus = Status::missing;
	Status vertStatus = Status::missing;
	Status geomStatus = Status::missing;
	Status fragStatus = Status::missing;
	Status linkStatus = Status::missing;

	//std::cout << _compute << std::endl;

	if (compute != -1) {
		glDeleteShader(compute);
		compute = -1;
	}
	if (vertex != -1) {
		glDeleteShader(vertex);
		vertex = -1;
	}
	if (geom != -1) {
		glDeleteShader(geom);
		geom = -1;
	}
	if (frag != -1) {
		glDeleteShader(frag);
		frag = -1;
	}
	if (program != -1) {
		glDeleteShader(program);
		program = -1;
	}

	//Compile Compute
	if (_compute != NULL && _compute[0] != '\0') {
		compute = glCreateShader(GL_COMPUTE_SHADER);
		glShaderSource(compute, 1, &_compute, nullptr);
		glCompileShader(compute);
		GLint isCompiled = 0;
		glGetShaderiv(compute, GL_COMPILE_STATUS, &isCompiled);
		if (isCompiled == GL_FALSE) {
			GLint maxLength = 0;
			glGetShaderiv(compute, GL_INFO_LOG_LENGTH, &maxLength);
			std::vector<GLchar> errorLog(maxLength);
			glGetShaderInfoLog(compute, maxLength, &maxLength, &errorLog[0]);
			glDeleteShader(compute);
			compStatus = Status::failed;
			print(id, compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
			return false;
		} else compStatus = Status::success;
	}

	//Compile Vertex
	if (_vertex != NULL && _vertex[0] != '\0') {
		vertex = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertex, 1, &_vertex, nullptr);
		glCompileShader(vertex);
		GLint isCompiled = 0;
		glGetShaderiv(vertex, GL_COMPILE_STATUS, &isCompiled);
		if (isCompiled == GL_FALSE) {
			GLint maxLength = 0;
			glGetShaderiv(vertex, GL_INFO_LOG_LENGTH, &maxLength);
			std::vector<GLchar> errorLog(maxLength);
			glGetShaderInfoLog(vertex, maxLength, &maxLength, &errorLog[0]);
			glDeleteShader(vertex);
			vertStatus = Status::failed;
			print(id, compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
			return false;
		} else vertStatus = Status::success;
	}

	//Compile Geom
	if (_geom != NULL && _geom[0] != '\0') {
		geom = glCreateShader(GL_GEOMETRY_SHADER);
		glShaderSource(geom, 1, &_geom, nullptr);
		glCompileShader(geom);
		GLint isCompiled = 0;
		glGetShaderiv(geom, GL_COMPILE_STATUS, &isCompiled);
		if (isCompiled == GL_FALSE) {
			GLint maxLength = 0;
			glGetShaderiv(geom, GL_INFO_LOG_LENGTH, &maxLength);
			std::vector<GLchar> errorLog(maxLength);
			glGetShaderInfoLog(geom, maxLength, &maxLength, &errorLog[0]);
			glDeleteShader(geom);
			geomStatus = Status::failed;
			print(id, compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
			return false;
		} else geomStatus = Status::success;
	}

	//Compile Frag
	if (_frag != NULL && _frag[0] != '\0') {
		frag = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(frag, 1, &_frag, nullptr);
		glCompileShader(frag);
		GLint isCompiled = 0;
		glGetShaderiv(frag, GL_COMPILE_STATUS, &isCompiled);
		if (isCompiled == GL_FALSE) {
			GLint maxLength = 0;
			glGetShaderiv(frag, GL_INFO_LOG_LENGTH, &maxLength);
			std::vector<GLchar> errorLog(maxLength);
			glGetShaderInfoLog(frag, maxLength, &maxLength, &errorLog[0]);
			glDeleteShader(frag);
			fragStatus = Status::failed;
			print(id, compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
			return false;
		} else fragStatus = Status::success;
	}

	//Link
	program = glCreateProgram();
	if (_compute != NULL && _compute[0] != '\0') glAttachShader(program, compute);
	if (_vertex != NULL && _vertex[0] != '\0') glAttachShader(program, vertex);
	if (_geom != NULL && _geom[0] != '\0') glAttachShader(program, geom);
	if (_frag != NULL && _frag[0] != '\0') glAttachShader(program, frag);

	glLinkProgram(program);

	GLint isLinked = 0;
	glGetProgramiv(program, GL_LINK_STATUS, (int*)&isLinked);
	if (isLinked == GL_FALSE) {
		GLint maxLength = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);
		std::vector<GLchar> errorLog(maxLength);
		glGetProgramInfoLog(program, maxLength, &maxLength, &errorLog[0]);
		if (compute != -1)glDeleteShader(compute);
		if (vertex != -1)glDeleteShader(vertex);
		if (geom != -1)glDeleteShader(geom);
		if (frag != -1)glDeleteShader(frag);
		if (program != -1) glDeleteProgram(program);
		linkStatus = Status::failed;

		print(id, compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
		return false;
	} else linkStatus = Status::success;

	if (_compute != NULL && _compute[0] != '\0')glDetachShader(program, compute);
	if (_vertex != NULL && _vertex[0] != '\0')glDetachShader(program, vertex);
	if (_geom != NULL && _geom[0] != '\0')glDetachShader(program, geom);
	if (_frag != NULL && _frag[0] != '\0')glDetachShader(program, frag);

	print(id, compStatus, vertStatus, geomStatus, fragStatus, linkStatus, "");

	unbind();
	return true;
}

GLuint GL::ShaderProgram::getHandle() const{
	return program;
}

GL::ShaderProgram::ShaderProgram(std::string _id) : id(_id) {}

GL::ShaderProgram::ShaderProgram() : ShaderProgram(""){}

GL::ShaderProgram::~ShaderProgram() {
	glDeleteProgram(program);
}

void GL::ShaderProgram::bind() const{
	glUseProgram(getHandle());
}

void GL::ShaderProgram::unbind() const{
	glUseProgram(0);
}

// -------------- Camera -------------- \\

GL::Camera::Camera(){

}

GL_Matrix4 GL::Camera::projection(float _n, float _f, float _l, float _r, float _b, float _t){
	//step 1
	const __m128 v11 = _mm_set_ps(_r, _t, _f, 0.f);
	const __m128 v12 = _mm_set_ps(_l, _b, _n, 0.f);
	const __m128 rp = _mm_add_ps(v11, v12);
	const __m128 rs = _mm_sub_ps(v11, v12);
	const float n2 = 2.f*_n;
	const float fn = -n2 * _f;

	//step 2
	const __m128 v21 = _mm_set_ps(n2, n2, fn, 0.f);
	__m128 rd1 = _mm_div_ps(v21, rs);
	__m128 rd2 = _mm_div_ps(rp, rs);

	util::array4f16a rf1;
	util::array4f16a rf2;
	_mm_store_ps(rf1.v, rd1);
	_mm_store_ps(rf2.v, rd2);

	GL_Matrix4 out;
	//rd1
	out(0, 0) = rf1[0];
	out(1, 1) = rf1[1];
	out(2, 3) = rf1[2];
	//rd2
	out(0, 2) = rf2[0];
	out(1, 2) = rf2[1];
	out(2, 2) = -rf2[2];

	out(3, 2) = -1.f;

	return out;
}

GL_Matrix4 GL::Camera::projection(float _n, float _f, float _r, float _t){
	const float fn = _f - _n;
	const __m128 v1 = _mm_set_ps(_n, _n, -(_f+_n), -2.f*_f*_n);
	const __m128 v2 = _mm_set_ps(_r, _t, fn, fn);
	const __m128 r = _mm_div_ps(v1, v2);

	util::array4f16a v;
	_mm_store_ps(v.v, r);

	GL_Matrix4 out;
	out(0, 0) = v[0];
	out(1, 1) = v[1];
	out(2, 2) = v[2];
	out(2, 3) = v[3];
	out(3, 2) = -1.f;

	return out;
}

GL_Matrix4 GL::Camera::orthographic(float _n, float _f, float _r, float _t){
	const float fn = _f - _n;
	const __m128 v1 = _mm_set_ps(1.f, 1.f, -2.f, _f+_n);
	const __m128 v2 = _mm_set_ps(_r, _t, fn, fn);
	const __m128 r = _mm_div_ps(v1, v2);

	util::array4f16a v;
	_mm_store_ps(v.v, r);

	GL_Matrix4 out;
	out(0, 0) = v[0];
	out(1, 1) = v[1];
	out(2, 2) = v[2];
	out(2, 3) = v[3];
	out(3, 3) = 1.f;

	return out;
}

GL_Matrix4 GL::Camera::lookAt(const Eigen::Vector3f& _eye, const Eigen::Vector3f& _target, const Eigen::Vector3f& _up){

	const Eigen::Vector3f f = (_eye - _target).normalized(); //forward
	const Eigen::Vector3f l = (_up.cross(f)).normalized(); //left
	const Eigen::Vector3f u = f.cross(l); //up
	

	const __m128 v1 = _mm_set_ps(-l[0], -u[0], -f[0], 0.f);
	const __m128 v2 = _mm_set_ps(-l[1], -u[1], -f[1], 0.f);
	const __m128 v3 = _mm_set_ps(-l[2], -u[2], -f[2], 0.f);

	const __m128 e1 = _mm_set_ps(_eye[0], _eye[0], _eye[0], 0.f);
	const __m128 e2 = _mm_set_ps(_eye[1], _eye[1], _eye[1], 0.f);
	const __m128 e3 = _mm_set_ps(_eye[2], _eye[2], _eye[2], 0.f);

	const __m128 rm1 = _mm_mul_ps(v1, e1);
	const __m128 rm2 = _mm_mul_ps(v2, e2);
	const __m128 rm3 = _mm_mul_ps(v3, e3);

	const __m128 ra1 = _mm_add_ps(rm1, rm2);
	const __m128 ra2 = _mm_add_ps(ra1, rm3);

	util::array4f16a v;
	_mm_store_ps(v.v, ra2);

	GL_Matrix4 out;
	out.setIdentity();
	out(0,0) = l[0];
	out(0,1) = l[1];
	out(0,2) = l[2];
	//out(0,3) = 0.f;

	out(1,0) = u[0];
	out(1,1) = u[1];
	out(1,2) = u[2];
	//out(1,3) = 0.f;

	out(2,0) = f[0];
	out(2,1) = f[1];
	out(2,2) = f[2];
	//out(2,3) = 0.f;

	out(3,0) = v[0];
	out(3,1) = v[1];
	out(3,2) = v[2];
	//out(3,3) = 1.f;

	return out;
}

GL_Matrix4 GL::Camera::rotate(float _angle_rad, const Eigen::Vector3f& _axis){
	const float c = std::cos(_angle_rad);
	const float s = std::sin(_angle_rad);
	const float cm1 = (1.f - c);

	const __m128 v = _mm_set_ps(_axis[0], _axis[1], _axis[2], 0.f);
	const __m128 sv = _mm_set_ps(s, s, s, 0.f);
	const __m128 cm1v = _mm_set_ps(cm1, cm1, cm1, 0.f);

	const __m128 vx = _mm_set_ps(_axis[0], _axis[0], _axis[0], 0.f);
	const __m128 vy = _mm_set_ps(_axis[1], _axis[1], _axis[1], 0.f);
	const __m128 vz = _mm_set_ps(_axis[2], _axis[2], _axis[2], 0.f);

	const __m128 r1 = _mm_mul_ps(cm1v, _mm_mul_ps(vx, v));
	const __m128 r2 = _mm_mul_ps(cm1v, _mm_mul_ps(vy, v));
	const __m128 r3 = _mm_mul_ps(cm1v, _mm_mul_ps(vz, v));
	const __m128 r4 = _mm_mul_ps(sv, v);

	util::array4f16a a1;
	_mm_store_ps(a1.v, r1);

	util::array4f16a a2;
	_mm_store_ps(a2.v, r2);

	util::array4f16a a3;
	_mm_store_ps(a3.v, r3);

	util::array4f16a a4;
	_mm_store_ps(a4.v, r4);

	GL_Matrix4 out;
	out(0, 0) = a1[0] + c;
	out(1, 0) = a1[1] + a4[2];
	out(2, 0) = a1[2] - a4[1];
	out(3, 0) = 0.f;

	out(0, 1) = a2[0] - a4[2];
	out(1, 1) = a2[1] + c;
	out(2, 1) = a2[2] + a4[0];
	out(3, 1) = 0.f;

	out(0, 2) = a3[0] + a4[1];
	out(1, 2) = a3[1] - a4[0];
	out(2, 2) = a3[2] - c;
	out(3, 2) = 0.f;

	out(0, 3) = 0.f;
	out(1, 3) = 0.f;
	out(2, 3) = 0.f;
	out(3, 3) = 1.f;

	return out;
		
}

/*
	code taken and adapted from https://github.com/Pascal-So/turbotrack
	with permission by Pascal Sommer, 2021
*/
Eigen::Vector3f GL::Camera::shoemake_projection(const Eigen::Vector2f& _mousePos, float _radius) {
	const float r2 = _radius * _radius;
	const float d2 = _mousePos.squaredNorm();

	if (d2 <= r2) {
		// sphere
		return {_mousePos[0], _mousePos[1], std::sqrt(r2 - d2)};
	} else {
		// scaled sphere
		const float factor = _radius / std::sqrt(d2);
		return {factor * _mousePos[0], factor *  _mousePos[1], 0};
	}
}

/*
	code taken and adapted from https://github.com/Pascal-So/turbotrack
	with permission by Pascal Sommer, 2021
*/
Eigen::Vector3f GL::Camera::holroyd_projection(const Eigen::Vector2f& _mousePos, float _radius) {
	const float r2 = _radius * _radius;
	const float d2 = _mousePos.squaredNorm();

	if (d2 <= r2 / 2) {
		// sphere
		return {_mousePos[0], _mousePos[1], std::sqrt(r2 - d2)};
	} else {
		// hyperbola
		return {_mousePos[0], _mousePos[1], r2 / 2 / std::sqrt(d2)};
	}
}

/*
	code taken and adapted from https://github.com/Pascal-So/turbotrack
	with permission by Pascal Sommer, 2021
*/
Eigen::Quaternionf GL::Camera::trackball_shoemake(const Eigen::Vector2f& _oldPos, const Eigen::Vector2f& _newPos, float _radius){
	const Eigen::Vector3f p1 = shoemake_projection(_oldPos, _radius);
	const Eigen::Vector3f p2 = shoemake_projection(_newPos, _radius);
	return Eigen::Quaternionf::FromTwoVectors(p1, p2);
}

/*
	code taken and adapted from https://github.com/Pascal-So/turbotrack
	with permission by Pascal Sommer, 2021
*/
Eigen::Quaternionf GL::Camera::trackball_holroyd(const Eigen::Vector2f& _oldPos, const Eigen::Vector2f& _newPos, float _radius){
	const Eigen::Vector3f p1 = holroyd_projection(_oldPos, _radius);
	const Eigen::Vector3f p2 = holroyd_projection(_newPos, _radius);
	return Eigen::Quaternionf::FromTwoVectors(p1, p2);
}