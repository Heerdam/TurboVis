#ifndef SHADERPROGRAM_HPP
#define SHADERPROGRAM_HPP

#include "defines.hpp"

namespace GL {

    class ShaderProgram {
        enum class Status { success, failed, missing };

        GLint program = -1, compute = -1, vertex = -1, geom = -1, frag = -1;
        void print(Status /*_compute*/, Status /*_vertex*/, Status /*_geometry*/, Status /*_fragment*/, Status /*_link*/, const std::string& /*_log*/) const noexcept;

    public:
        ShaderProgram(const std::string& /*_id*/) noexcept;
        ShaderProgram() noexcept;
        ~ShaderProgram() noexcept;

        inline static std::string SHADERFOLDERPATH = "";

        std::string id;
        
        inline static std::function<void(const std::string& /*_msg*/)> cb_debug = 0;

        /*
            assumes the following extensions:
            compute shader: [FILENAME].comp
            vertex shader: [FILENAME].ver
            geometry shader: [FILENAME].geo
            fragment shader: [FILENAME].frag
        */
        bool compileFromFile(const std::string& /*_filename*/);
        bool compile(const char* /*_compute*/, const char* /*_vertex*/, const char* /*_geometry*/, const char* /*_fragment*/);

        [[nodiscard]] GLuint getHandle() const noexcept;

        void bind() const noexcept;
        void unbind() const noexcept;
    };

}

inline void GL::ShaderProgram::print(
        ShaderProgram::Status _compComp, 
        ShaderProgram::Status _compVert,
        ShaderProgram::Status _compGeom, 
        ShaderProgram::Status _compFrag, 
        ShaderProgram::Status _link, 
        const std::string& _errorLog
    ) const noexcept{

    if(cb_debug == 0) return;

    std::stringstream s;
    s << "\n   Shader: " << id << std::endl
    << "Compiling: "
    << (_compComp == Status::failed ? " X |" : _compComp == Status::success ? " S |" : " - |")
    << (_compVert == Status::failed ? " X |" : _compVert == Status::success ? " S |" : " - |")
    << (_compGeom == Status::failed ? " X |" : _compGeom == Status::success ? " S |" : " - |")
    << (_compFrag == Status::failed ? " X |" : _compFrag == Status::success ? " S |" : " - |") << std::endl
    << "  Linking: " + std::string(_link == Status::failed ? "Failed!" : _link == Status::success ? "Success!" : " - ") << std::endl
    << "Error Log: " << (_errorLog.empty() ? "empty" : _errorLog) << std::endl;

    cb_debug(s.str());
};

inline bool GL::ShaderProgram::compileFromFile(const std::string& _filename) {
    bool cExists = true;
    bool vExists = true;
    bool gExists = true;
    bool fExists = true;

    std::stringstream p;
    p << SHADERFOLDERPATH << _filename;
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
        spdlog::error("no valid path for shader: {}", _filename);
        return false;
    }

    id = _filename;

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

inline bool GL::ShaderProgram::compile(const char* _compute, const char* _vertex, const char* _geom, const char* _frag) {
    Status compStatus = Status::missing;
    Status vertStatus = Status::missing;
    Status geomStatus = Status::missing;
    Status fragStatus = Status::missing;
    Status linkStatus = Status::missing;

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
            print(compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
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
            print(compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
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
            print(compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
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
            print(compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
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

        print(compStatus, vertStatus, geomStatus, fragStatus, linkStatus, std::string(errorLog.begin(), errorLog.end()));
        return false;
    } else linkStatus = Status::success;

    if (_compute != NULL && _compute[0] != '\0')glDetachShader(program, compute);
    if (_vertex != NULL && _vertex[0] != '\0')glDetachShader(program, vertex);
    if (_geom != NULL && _geom[0] != '\0')glDetachShader(program, geom);
    if (_frag != NULL && _frag[0] != '\0')glDetachShader(program, frag);

    print(compStatus, vertStatus, geomStatus, fragStatus, linkStatus, "");

    unbind();
    return true;
}

inline GLuint GL::ShaderProgram::getHandle() const noexcept{
    return program;
}

inline GL::ShaderProgram::ShaderProgram(const std::string& _id) noexcept : id(_id) {}

inline GL::ShaderProgram::ShaderProgram() noexcept : ShaderProgram(""){}

inline GL::ShaderProgram::~ShaderProgram() noexcept{
    glDeleteProgram(program);
}

inline void GL::ShaderProgram::bind() const noexcept{
    glUseProgram(getHandle());
}

inline void GL::ShaderProgram::unbind() const noexcept{
    glUseProgram(0);
}

#endif //SHADERPROGRAM_HPP
