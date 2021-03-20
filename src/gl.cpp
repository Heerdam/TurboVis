#include "../include/gl.hpp"

// -------------- SHADERPROGRAM -------------- 

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

	id = _path;

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

/*
	code taken and adapted from https://github.com/Pascal-So/turbotrack
	with permission by Pascal Sommer, 2021
*/
Vec3 GL::Camera::shoemake_projection(const Vec2& _mousePos, float _radius) {
	const float r2 = _radius * _radius;
	const float d2 = glm::dot(_mousePos, _mousePos);

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
Vec3 GL::Camera::holroyd_projection(const Vec2& _mousePos, float _radius) {
	const float r2 = _radius * _radius;
	const float d2 = glm::dot(_mousePos, _mousePos);

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
Quat GL::Camera::trackball_shoemake(const Vec2& _oldPos, const Vec2& _newPos, float _radius){
	const Vec3 p1 = glm::normalize(shoemake_projection(_oldPos, _radius));
	const Vec3 p2 = glm::normalize(shoemake_projection(_newPos, _radius));
	return glm::rotation(p1, p2);
}

/*
	code taken and adapted from https://github.com/Pascal-So/turbotrack
	with permission by Pascal Sommer, 2021
*/
Quat GL::Camera::trackball_holroyd(const Vec2& _oldPos, const Vec2& _newPos, float _radius){
	const Vec3 p1 = glm::normalize(holroyd_projection(_oldPos, _radius));
	const Vec3 p2 = glm::normalize(holroyd_projection(_newPos, _radius));
	return glm::rotation(p1, p2);
}

//max 16k vertices
GL::ShapeRenderer::ShapeRenderer(uint32_t _maxVertices) { //shape == triangle
    //create buffers
    {
        glGenVertexArrays(2, VAO);

        glGenBuffers(2, VBO_POS);
		glGenBuffers(2, VBO_COL);
		glGenBuffers(2, VBO_NRM);
        glGenBuffers(2, EBO);

        for (uint32_t i = 0; i < 2; ++i) {
            glBindVertexArray(VAO[i]);

            //vbo_pos
            glBindBuffer(GL_ARRAY_BUFFER, VBO_POS[i]);
            glBufferStorage(GL_ARRAY_BUFFER, _maxVertices * 3 * sizeof(float), nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
            VBO_ptr_pos[i] = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, _maxVertices * 3 * sizeof(float), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (GLvoid*)(0));
            glEnableVertexAttribArray(0);

            //col
			glBindBuffer(GL_ARRAY_BUFFER, VBO_POS[i]);
            glBufferStorage(GL_ARRAY_BUFFER, _maxVertices * sizeof(float), nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
            VBO_ptr_col[i] = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, _maxVertices * 3 * sizeof(float), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));
            glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, 1 * sizeof(float), (GLvoid*)(0));
            glEnableVertexAttribArray(1);

			//normals
			glBindBuffer(GL_ARRAY_BUFFER, VBO_POS[i]);
            glBufferStorage(GL_ARRAY_BUFFER, _maxVertices * 3 * sizeof(float), nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
            VBO_ptr_nrm[i] = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, _maxVertices * 3 * sizeof(float), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));
			glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (GLvoid*)(0));
            glEnableVertexAttribArray(2);

            //ebo
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO[i]);
            glBufferStorage(GL_ELEMENT_ARRAY_BUFFER, _maxVertices * 4 * sizeof(uint32_t), nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
            EBO_ptr[i] = reinterpret_cast<uint16_t*>(glMapBufferRange(GL_ELEMENT_ARRAY_BUFFER, 0, _maxVertices * 4 * sizeof(uint32_t), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));

            glBindVertexArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
    }
    //compile shaders
    {
        const GLchar* vertex_shader =
            "#version 430 core\n"
            "layout (location = 0) in vec3 Position;\n"
            "layout (location = 1) in vec4 Color;\n"
			"layout (location = 2) in vec3 Normals;\n"
            "layout (location = 3) uniform mat4 ProjMtx;\n"
            "out vec4 Frag_Color;\n"
			"out vec3 Nrm;\n"
            "void main(){\n"
            "    Frag_Color = Color;\n"
            "    gl_Position = ProjMtx * vec4(Position.xy,0.f,1.f);\n"
			"    Nrm = Normals;\n"
            "}\n";

        const GLchar* fragment_shader =
            "#version 430 core\n"
            "in vec4 Frag_Color;\n"
			"in vec3 Nrm;\n"
            "out vec4 Out_Color;\n"
			"const vec3 ldir = normalize(vec3(1., 1., 1.));\n"
            "void main(){\n"
			"   const float diff = max(dot(-ldir, Nrm, 0.);"
            "   Out_Color = diff * Frag_Color;\n"
            "}\n";

		shader.id = "ShapeRenderer";
        shader.compile("", vertex_shader, "", fragment_shader);
    }
}

void GL::ShapeRenderer::render(const float* _camera) {
    shader.bind();
    glBindVertexArray(VAO[currIndex]);

    glUniformMatrix4fv(2, 1, GL_FALSE, _camera);

    glDrawElements(GL_TRIANGLES, currIdx, GL_UNSIGNED_SHORT, (void*)0);
    //glDrawArrays(GL_TRIANGLES, 0, 1);

    glBindVertexArray(0);
    shader.unbind();

    currIdx = currVert = 0;
    currIndex = (currIndex + 1) % 2;
}

void GL::ShapeRenderer::drawLine(const Vec3& _p1, const Vec3& _p2, uint32_t _segments, float _thickness, const Vec4& _col) {

	_thickness = std::max(0.1f, _thickness);

    const float l = glm::length(_p2 - _p1);
    const Vec3 dir = glm::normalize(_p2 - _p1);

    const float rad = glm::radians(360.f / _segments);
	std::vector<Vec3> cpos(_segments);
	auto rot = glm::rotation({ 0.f, 0.f, 0.f }, dir);

	//vertices
    for (uint32_t i = 0; i < _segments; ++i){		
       cpos.push_back(_p1 + (rot * (Vec3(std::cos(i * rad), std::sin(i * rad), 0.f) * _thickness))); 
	}

	for (uint32_t i = _segments; i < 2*_segments; ++i)
       cpos.push_back(cpos[i - _segments] * (dir * l));     

	//normals

	//color
	const int col = (static_cast<int>(_col[3] * 255.f) << 24) |
					(static_cast<int>(_col[2] * 255.f) << 16) |
					(static_cast<int>(_col[1] * 255.f) << 8) |
					static_cast<int>(_col[0] * 255.f);
	std::memset(VBO_ptr_col[currIndex] + currVert, col, cpos.size());

	//indices


}

void GL::ShapeRenderer::drawAABB(const Vec3& _low, const Vec3& _high, float _thickness, const Vec4& _col) {
	/*
    const Vec3 p1 = _low;
    const Vec3 p2 = Vec3(_high[0], _low[1]);
    const Vec3 p3 = _high;
    const Vec3 p4 = Vec3(_low[0], _high[1]);
	const Vec3 p5 = _high;
	const Vec3 p6 = _high;
	const Vec3 p7 = _high;
	const Vec3 p8 = _high;

	//bottom
    drawLine(p1, p2, _thickness, _col);
    drawLine(p2, p3, _thickness, _col);
    drawLine(p3, p4, _thickness, _col);
    drawLine(p4, p1, _thickness, _col);

	//top
	drawLine(p1, p2, _thickness, _col);
    drawLine(p2, p3, _thickness, _col);
    drawLine(p3, p4, _thickness, _col);
    drawLine(p4, p1, _thickness, _col);

	//sides
	drawLine(p1, p2, _thickness, _col);
    drawLine(p2, p3, _thickness, _col);
    drawLine(p3, p4, _thickness, _col);
    drawLine(p4, p1, _thickness, _col);

	//circles
	drawSphere();
	drawSphere();
	drawSphere();
	drawSphere();

	drawSphere();
	drawSphere();
	drawSphere();
	drawSphere();
	*/
}

void GL::ShapeRenderer::drawSphere(const Vec3& _centre, float _radius, uint32_t _subdivisions, const Vec4& _col) {
    if (!spheres[_subdivisions].has_value())
        spheres[_subdivisions] = {GL::util::Geometry::Icosahedron(_subdivisions)};

	// -------------- Vertices -------------- 
	{
		const Vector_af32& vrt = spheres[_subdivisions].value().first;

		const __m256 offset = _mm256_set_ps(_centre.x, _centre.y, _centre.z, _centre.x, _centre.y, _centre.z, 0.f, 0.f);
		const __m256 radius = _mm256_set_ps(_radius, _radius, _radius, _radius, _radius, _radius, 0.f, 0.f);

		util::array4f32a vals;

		for (size_t i = 0; i < vrt.size(); ++i) {
			const uint32_t vo = currVert * 3;
			//normals
			_mm256_store_ps(vals.v, vrt[i]);
			std::memcpy(VBO_ptr_nrm[currIndex] + vo, vals.v, sizeof(float) * 6);
			//pos
			_mm256_store_ps(vals.v, _mm256_fmadd_ps(vrt[i], radius, offset));
			std::memcpy(VBO_ptr_pos[currIndex] + vo, vals.v, sizeof(float) * 6);
			currVert += 2;
		}

		//color
		const int col = (static_cast<int>(_col[3] * 255.f) << 24) |
						(static_cast<int>(_col[2] * 255.f) << 16) |
						(static_cast<int>(_col[1] * 255.f) << 8) |
						static_cast<int>(_col[0] * 255.f);
		std::memset(VBO_ptr_col[currIndex] + currVert, col, vrt.size()*2);
	}

	// -------------- INDICES -------------- 
	{
		const Vector_aui16& idx = spheres[_subdivisions].value().second;
		util::array4ui16a vals;
		std::memset(vals.v, currIdx, sizeof(vals.v));

		const __m256i offset = _mm256_load_si256((__m256i const *)vals.v);

		for (size_t i = 0; i < idx.size()-1; ++i){
			_mm256_storeu_epi16(vals.v, _mm256_add_epi16(idx[i], offset));
			std::memcpy(EBO_ptr[currIndex] + currIdx + i*16, vals.v, sizeof(vals.v));
		}
		//last element
		_mm256_storeu_epi16(vals.v, _mm256_add_epi16(idx[idx.size()-1], offset));
		size_t lastEl = 15;
		for(; lastEl >= 0; lastEl--)
			if(vals.v[lastEl] != currIdx) break; //TODO: this is very fishy
		std::memcpy(EBO_ptr[currIndex] + currIdx + (idx.size() - 1)*16, vals.v, sizeof(vals.v));

		currIdx += static_cast<uint32_t>(spheres[_subdivisions].value().first.size())*2;
	}

}

std::pair<GL::Vector_af32, GL::Vector_aui16> GL::util::Geometry::Icosahedron(uint16_t _subdivisions) {
    const float X = .525731112119133606f;
    const float Z = .850650808352039932f;
    const float N = 0.f;

    struct Triangle {
        uint16_t vertex[3];
    };

    using TriangleList = std::vector<Triangle>;
    using VertexList = std::vector<Vec3>;

    const VertexList vertices =
    {
        {-X,N,Z}, {X,N,Z}, {-X,N,-Z}, {X,N,-Z},
        {N,Z,X}, {N,Z,-X}, {N,-Z,X}, {N,-Z,-X},
        {Z,X,N}, {-Z,X, N}, {Z,-X,N}, {-Z,-X, N}
    };

    const TriangleList triangles =
    {
        {0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
        {8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
        {7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
        {6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
    };

    using Lookup = std::map<std::pair<uint16_t, uint16_t>, uint16_t>;

    auto vertex_for_edge = [&](Lookup& lookup, VertexList& vertices, uint16_t first, uint16_t second) -> uint16_t {
        Lookup::key_type key(first, second);
        if (key.first > key.second)
            std::swap(key.first, key.second);

        auto inserted = lookup.insert({key, static_cast<uint16_t>(vertices.size())});
        if (inserted.second) {
            auto& edge0 = vertices[first];
            auto& edge1 = vertices[second];
            auto point = normalize(edge0 + edge1);
            vertices.push_back(point);
        }

        return inserted.first->second;
    };

    const auto subdivide = [&](VertexList& vertices, TriangleList triangles) -> TriangleList {
        Lookup lookup;
        TriangleList result;

        for (auto&& each : triangles) {
            std::array<uint16_t, 3> mid;
            for (int edge = 0; edge < 3; ++edge) {
                mid[edge] = vertex_for_edge(lookup, vertices,
                                            each.vertex[edge], each.vertex[(edge + 1) % 3]);
            }

            result.push_back({each.vertex[0], mid[0], mid[2]});
            result.push_back({each.vertex[1], mid[1], mid[0]});
            result.push_back({each.vertex[2], mid[2], mid[1]});
            result.push_back({mid[0], mid[1], mid[2]});
        }

        return result;
    };

    VertexList vert = vertices;
    TriangleList tria = triangles;

    for (uint16_t i = 0; i < _subdivisions; ++i) {
        tria = subdivide(vert, tria);
    }

	const bool isEven = vert.size()%2;
	if(!isEven)
		vert.emplace_back(0.f, 0.f, 0.f);

    Vector_af32 v_out;
    v_out.reserve(vert.size() / 2);
    for (size_t i = 0, j = 0; i < vert.size(); i += 2, ++j) {
        auto v1 = glm::normalize(vert[i]);
        auto v2 = glm::normalize(vert[i + 1]);
        v_out.push_back(_mm256_set_ps(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, 0.f, 0.f));
    }

    Vector_aui16 i_out;
	i_out.reserve((tria.size() * 3) / 16);

	std::vector<uint16_t, aligned_allocator<uint16_t, sizeof(uint16_t)>> buffer;
	buffer.reserve(16);
    for (size_t i = 0, j = 0; i < tria.size(); ++i) {       
		auto& t = tria[i];
		for(size_t j = 0; j < 3; ++j){

			if(buffer.size() == 16){
				i_out.push_back(_mm256_load_si256((__m256i const *)buffer.data()));
				buffer.clear();
			}

			buffer.push_back(t.vertex[j]);
		}
    }

	if(!buffer.empty())
		i_out.push_back(_mm256_load_si256((__m256i const *)buffer.data()));


    return { v_out, i_out };
}

std::pair<GL::Vector_af32, GL::Vector_aui16>Zylinder(uint16_t _subdivisions){

}