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

GL::ShapeRenderer::ShapeRenderer(uint32_t _maxShapesPerType) : maxShapes(_maxShapesPerType) { 
    {
        spheres = zyls = planes = 0;

        // ------------------- FULL SCREEN QUAD -------------------

        const float quad[] = {
            1.0f,  1.0f,  // top right
            1.0f, -1.0f,  // bottom right
            -1.0f, -1.0f,  // bottom left
            -1.0f,  1.0f,   // top left 
        };

        const int16_t quadI[] = {
            0, 1, 2,
            2, 3, 0
        };

        GLuint VBO, EBO;

        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);

        //vbo
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);  

        //ebo      
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadI), quadI, GL_STATIC_DRAW);

        glBindVertexArray(0);

        // ------------------- SSBO -------------------
        
        glGenBuffers(2, SSBO_spheres);
        glGenBuffers(2, SSBO_zyls);
        glGenBuffers(2, SSBO_planes);

        for (size_t i = 0; i < 2; ++i) {
            //spheres
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, SSBO_spheres[i]);
            glBufferStorage(GL_SHADER_STORAGE_BUFFER, _maxShapesPerType * sizeof(Sphere), nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
            SSBO_ptr_spheres[i] = reinterpret_cast<Sphere*>(glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, _maxShapesPerType * sizeof(Sphere), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

            //zylinders
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, SSBO_zyls[i]);
            glBufferStorage(GL_SHADER_STORAGE_BUFFER, _maxShapesPerType * sizeof(Zylinder), nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
            SSBO_ptr_zyls[i] = reinterpret_cast<Zylinder*>(glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, _maxShapesPerType * sizeof(Zylinder), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

            //planes
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, SSBO_planes[i]);
            glBufferStorage(GL_SHADER_STORAGE_BUFFER, _maxShapesPerType * sizeof(Plane), nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
            SSBO_ptr_planes[i] = reinterpret_cast<Plane*>(glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, _maxShapesPerType * sizeof(Plane), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);       
        }
        
    }

    //compile shaders
    shader.compileFromFile("shader/shape");
}

void GL::ShapeRenderer::render(const Camera& _cam) {

    shader.bind();

    
    
    //buffers
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, SSBO_spheres[currIndex]);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, SSBO_zyls[currIndex]);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, SSBO_planes[currIndex]);

    //uniforms
    glUniform3fv(10, 1, glm::value_ptr(_cam.position));
    glUniformMatrix4fv(11, 1, GL_FALSE, glm::value_ptr(_cam.combined));

    glUniform3fv(12, 1, glm::value_ptr(_cam.dir));
    glUniform3fv(13, 1, glm::value_ptr(_cam.right));
    glUniform3fv(14, 1, glm::value_ptr(_cam.up));

    glUniform2f(15, float(_cam.width), float(_cam.height));

    glUniform1i(16, spheres);
    glUniform1i(17, zyls);
    glUniform1i(18, planes);

    //glUniform1f(31, _cam.fov);

    //draw
    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, (void*)0);
    glBindVertexArray(0);

    shader.unbind();

    currIndex = (currIndex + 1) % 2;
    spheres = zyls = planes = 0;
}

void GL::ShapeRenderer::drawLine(const Vec3& _p1, const Vec3& _p2, float _radius, const Vec3& _col) {
    assert(zyls < maxShapes);
    Zylinder zyl { _radius, _p1, _p2, _col };
    std::memcpy(SSBO_ptr_zyls[currIndex] + zyls++, &zyl, sizeof(Zylinder));
}

void GL::ShapeRenderer::drawAABB(const Vec3& _low, const Vec3& _high, float _thickness, const Vec3& _col) {
  
}

void GL::ShapeRenderer::drawSphere(const Vec3& _centre, float _radius, const Vec3& _col) {
    assert(spheres < maxShapes);
    Sphere s { _radius, _centre, _col };
    std::memcpy(SSBO_ptr_spheres[currIndex] + spheres++, &s, sizeof(Sphere));
}

void GL::ShapeRenderer::drawAxisWidget(){
    drawSphere(Vec3(0.f), 20.f, Vec4(0.f, 1.f, 0.f, 1.f));
    drawLine(Vec3(0.f), Vec3(300.f, 0.f, 0.f), 20.f, Vec4(1.f, 0.f, 0.f, 1.f)); //x
    drawLine(Vec3(0.f), Vec3(0.f, 300.f, 0.f), 20.f, Vec4(0.f, 1.f, 0.f, 1.f)); //y
    drawLine(Vec3(0.f), Vec3(0.f, 0.f, 300.f), 20.f, Vec4(0.f, 0.f, 1.f, 1.f)); //z
}