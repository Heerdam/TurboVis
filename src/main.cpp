
#include "../include/defines.hpp"
#include "../include/gui.hpp"
#include "../include/data.hpp"
#include "../include/hagedornrenderer.hpp"
#include "../include/bench.hpp"
#include "../include/util.hpp"
#include "../include/camera.hpp"
#include "../include/hdf5.hpp"
#include "../include/shaderprogram.hpp"

void GLAPIENTRY MessageCallback(GLenum /*source*/, GLenum type, GLuint /*id*/, GLenum severity, GLsizei /*length*/, const GLchar* message, const void* /*userParam*/) {
	if (type != GL_DEBUG_TYPE_ERROR) return;
	fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
		(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""),
		type, severity, message);
}

int main() {

    FilePathResolver::resolve();

    spdlog::set_level(spdlog::level::trace);
    if (!glfwInit()) {
        spdlog::error("Failed to init glfw. Shutting down...", true);
        return EXIT_FAILURE;
    }

    GL::ShaderProgram::SHADERFOLDERPATH = FilePathResolver::SHADERDIR();
    GL::ShaderProgram::cb_debug = [&](const std::string& _msg){
        spdlog::debug(_msg);
    };

    glfwWindowHint(GLFW_REFRESH_RATE, 144);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
    glfwWindowHint(GLFW_DEPTH_BITS, 24);
    glfwWindowHint(GLFW_STENCIL_BITS, 8);

    const int WIDTH = 1920;
    const int HEIGHT = 1080;
    const float HALFWIDTH = WIDTH * 0.5f;
    const float HALFHEIGHT = HEIGHT * 0.5f;

    auto window = glfwCreateWindow(WIDTH, HEIGHT, "TurboVis", NULL, NULL);
    if (!window) {
        spdlog::error("Failed to create window. Shutting down...", true);
        return EXIT_FAILURE;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        spdlog::error("ERROR:\tFailed to init GLAD. Shutting down...", true);
        return EXIT_FAILURE;
    }

    if (!glfwExtensionSupported("GL_ARB_buffer_storage")) {
        spdlog::error("extension GL_ARB_buffer_storage not supported");
        return EXIT_FAILURE;
    }

	glDebugMessageCallback(MessageCallback, 0);

    glEnable(GL_STENCIL_TEST);
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEPTH_TEST);
    // glEnable(GL_SRGB); //why error?

    // -------------- GUI --------------
    Gui::InputMultiplexer::init(window);
    Gui::FrontendGui gui(window, WIDTH, HEIGHT);
    double time = glfwGetTime();
    double deltaTime = 0.;
    unsigned long long frame = 0;
    size_t frames = 0, fps = 0, maxfps = 0;

    // -------------- CAMERA --------------
    GL::Camera camera = GL::Camera(int64_t(WIDTH), int64_t(HEIGHT), glm::radians(55.f), 0.01f, 5.f);
    camera.position = { 0.f, 0.f, 3.5f };
    camera.target = { 0.f, 0.f, 0.f };
    //camera.combined = glm::perspectiveFov(camera.fov, float(camera.width), float(camera.height), camera.near, camera.far);
    camera.combined = glm::lookAt(camera.position, camera.target, camera.upAxis);
    camera.comb =  glm::perspective(glm::radians(55.f), float(WIDTH)/float(HEIGHT), 0.1f, 100.f) * camera.combined;
    camera.update();

    // -------------- GIMBEL CAMERA --------------
    GL::Camera gCamera = GL::Camera(int64_t(250), int64_t(250), glm::radians(55.f), 0.1f, 1500.f);
    gCamera.position = { 0.f, 0.f, 100.f };
    gCamera.target = { 0.f, 0.f, 0.f };
    gCamera.update();

    bool RMB_down = false;

    double mX, mY;
    glfwGetCursorPos(window, &mX, &mY);
    Vec2 oldPosition = Vec2((float)mY - HALFHEIGHT, (float)mX - HALFWIDTH);
    Vec2 newPosition = oldPosition;

    Gui::InputMultiplexer::mouseButtonCallback([&](GLFWwindow*, int _button, int _action, int _mods)-> void{
        if(_button == GLFW_MOUSE_BUTTON_RIGHT && _mods == 0x0){
            RMB_down = _action == GLFW_PRESS ? true : _action == GLFW_RELEASE ? false : RMB_down;
        }
    });

    bool camKeys[] = {false, false, false};

    // -------------- RENDERER --------------
    //GL::ShapeRenderer shape (100);
    //GL::DepthBufferVisualizer depthr(camera);
    //GL::RaymarchTester rayM;

    //load file
    const auto file = IO::getExample();
    GL::HagedornRenderer<double, GL::Camera> hager(camera);
    hager.set(file);
    
    hager.steps = 500;

    Gui::InputMultiplexer::keyCallback([&](GLFWwindow* _window, int _key, int _scancode, int _action, int _mods)-> void {
        if(_key == GLFW_KEY_X && _action == GLFW_PRESS){
            hager.steps = 250;
            //hager.k += 0.001;
            std::cout << hager.k << std::endl;
            camera.hasMoved = true; 
        } 
        if(_key == GLFW_KEY_Z && _action == GLFW_PRESS){
            hager.steps += 250;
            std::cout << hager.steps << std::endl;
            camera.hasMoved = true;
        } 

        if(_key == GLFW_KEY_R)
            camera.hasMoved = true;

    });

    Gui::InputMultiplexer::cursorPosCallback([&](GLFWwindow*, double _xpos, double _ypos)-> void{
        oldPosition = newPosition;
        newPosition = Vec2(float(_xpos) - HALFWIDTH, (HEIGHT - float(_ypos) - HALFHEIGHT));

        
        if(RMB_down){
            camera.hasMoved = !(oldPosition == newPosition);
            const auto rot = GL::Camera::trackball_holroyd(oldPosition, newPosition, 250.f);
            const glm::mat4 RotationMatrix = glm::toMat4(glm::normalize(rot));

            camera.combined *= RotationMatrix;
            camera.comb *= RotationMatrix;
            //std::cout << glm::to_string(camera.combined) << std::endl;
            
            camera.position = glm::transpose(glm::mat3(camera.combined)) * Vec3(glm::column(camera.combined, 3)); //3
            camera.dir = glm::normalize(Vec3(-glm::row(camera.combined, 2))); //-2
            camera.right = glm::normalize(Vec3(glm::row(camera.combined, 0))); //0
            camera.up = glm::normalize(Vec3(glm::row(camera.combined, 1))); //1

            //std::cout << glm::to_string(camera.position) << std::endl;

            //camera.update();
            //camera.combined = glm::perspectiveFov(camera.fov, float(camera.width), float(camera.height), camera.near, camera.far);
            //camera.combined *= glm::lookAt(camera.position, camera.target, camera.upAxis);
        }
    });

    hager.start(camera);

    while (!glfwWindowShouldClose(window)) {
        const double ctime = glfwGetTime();

        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        glfwPollEvents();

        glEnable(GL_DEPTH_TEST);

        // -------------- FUNCTION --------------  
        //rayM.render(camera, uniforms);
        //hager.k += 1./60.;
        hager.render(camera);
        camera.hasMoved = false;
         //depthr.render();

/*
        // -------------- GIMBEL --------------  
        glViewport(0, 0, gCamera.width, gCamera.height);
        const Mat3 rot = Mat3(camera.dir, camera.up, camera.right);
        //shape.drawAxisWidget(rot);
        //shape.render(gCamera);
        glViewport(0, 0, WIDTH, HEIGHT);
       
*/
        // -------------- GUI --------------     
        Gui::RenderInfo info;
        info.width = WIDTH;
        info.height = HEIGHT;
        info.fps = fps;
        info.maxfps = maxfps;
        info.progress = hager.getProgress();
        info.steps = hager.steps;
        gui.draw(window, info);   

        if(hager.steps != info.steps){
            hager.steps = info.steps;
            camera.hasMoved = true;
        }    
          

        glfwSwapBuffers(window);

        frame++;
        if (ctime - time >= 1.) {
            fps = frame;
            maxfps = std::max(fps, maxfps);
            frame = 0;
            time = ctime;
        }
        deltaTime = glfwGetTime() - ctime;
        //if(uniforms.tt) 
            //uniforms.t += deltaTime;
    }
    hager.stop();
    glfwTerminate();

    return EXIT_SUCCESS;
}
