
#include "../include/defines.hpp"
#include "../include/gui.hpp"
#include "../include/data.hpp"
#include "../include/gl.hpp"
#include "../include/bench.hpp"
#include "../include/util.hpp"
#include "../include/camera.hpp"

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
    glfwSwapInterval(0);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        spdlog::error("ERROR:\tFailed to init GLAD. Shutting down...", true);
        return EXIT_FAILURE;
    }

    if (!glfwExtensionSupported("GL_ARB_buffer_storage")) {
        spdlog::error("extension GL_ARB_buffer_storage not supported");
        return EXIT_FAILURE;
    }

	//glDebugMessageCallback(MessageCallback, 0);

    glEnable(GL_STENCIL_TEST);
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEPTH_TEST);
    // glEnable(GL_SRGB); //why error?

    // -------------- GUI --------------
    Gui::InputMultiplexer::init(window);
    Gui::GLGui gui;
    gui.initGui(window, WIDTH, HEIGHT);
    try {
        gui.initGL((unsigned int)(100000), (unsigned int)(200000));
    } catch (const TurboGUI::TurboGuiException& e) {
        spdlog::error("{}", e.what());
    }

    bool show_demo_window = false;

    double time = glfwGetTime();
    double deltaTime = 0.;
    unsigned long long frame = 0;
    size_t frames = 0, fps = 0;

    // -------------- CAMERA --------------
    GL::Camera camera = GL::Camera(int64_t(WIDTH), int64_t(HEIGHT), glm::radians(55.f), 0.1f, 1500.f);
    camera.position = { 0.f, 0.f, 100.f };
    camera.target = { 0.f, 0.f, 0.f };
    camera.combined = glm::perspectiveFov(camera.fov, float(camera.width), float(camera.height), camera.near, camera.far);
    camera.combined *= glm::lookAt(camera.position, camera.target, camera.upAxis);
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

    Gui::InputMultiplexer::keyCallback([&](GLFWwindow* _window, int _key, int _scancode, int _action, int _mods)-> void {
        if(_key == GLFW_KEY_X)
            camKeys[0] = _action == GLFW_PRESS ? true :  _action == GLFW_RELEASE ? false : camKeys[0];
        if(_key == GLFW_KEY_Y)
            camKeys[1] = _action == GLFW_PRESS ? true :  _action == GLFW_RELEASE ? false : camKeys[1];
        if(_key == GLFW_KEY_Z)
            camKeys[2] = _action == GLFW_PRESS ? true :  _action == GLFW_RELEASE ? false : camKeys[2];
         
    });

    Gui::InputMultiplexer::cursorPosCallback([&](GLFWwindow*, double _xpos, double _ypos)-> void{
        oldPosition = newPosition;
        newPosition = Vec2((float)_xpos - HALFWIDTH, (HEIGHT - (float)_ypos - HALFHEIGHT));

        if(RMB_down){
            const auto rot = GL::Camera::trackball_holroyd(oldPosition, newPosition, 200.f);
            const glm::mat3 RotationMatrix = glm::toMat3(rot);

            //rotation around x axis
            if(camKeys[0]){
                camera.right = normalize(rot * camera.right);
                camera.up = normalize(glm::cross(camera.dir, camera.right));
            }

            //rotation around y axis
            if(camKeys[1]){
                camera.dir = normalize(rot * camera.dir);
                camera.right = normalize(glm::cross(camera.dir, camera.up));
            }

            //rotation around z axis
            if(camKeys[2]){
                camera.up = normalize(rot * camera.up);
                camera.dir = normalize(glm::cross(camera.up, camera.right));
            }

            //const Vec3 rpos = camera.position - camera.target;
            //camera.position = RotationMatrix * rpos  + camera.target;
            //camera.update();
            camera.combined = glm::perspectiveFov(camera.fov, float(camera.width), float(camera.height), camera.near, camera.far);
            camera.combined *= glm::lookAt(camera.position, camera.target, camera.upAxis);
        }
    });

    GL::ShapeRenderer shape (100);
    GL::DepthBufferVisualizer depthr(camera);

    while (!glfwWindowShouldClose(window)) {
        const double ctime = glfwGetTime();

        glClearColor(0.1f, 0.f, 0.25f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        glfwPollEvents();
        gui.ImGui_ImplGlfw_UpdateMousePosAndButtons(window);

        glEnable(GL_DEPTH_TEST);

        // -------------- FUNCTION --------------  


        // -------------- GIMBEL --------------  
        glViewport(0, 0, gCamera.width, gCamera.height);
        const Mat3 rot = Mat3(camera.dir, camera.up, camera.right);
        shape.drawAxisWidget(rot);
        shape.render(gCamera);
        glViewport(0, 0, WIDTH, HEIGHT);
        //depthr.render();

        // -------------- GUI --------------       
        gui.begin();
        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);
        gui.drawStatsWindow((uint32_t)fps);
        gui.draw();
        gui.sync();
        

        glfwSwapBuffers(window);

        frame++;
        if (ctime - time >= 1.0) {
            fps = frame;
            frame = 0;
            time = ctime;
        }
        deltaTime = glfwGetTime() - ctime;
    }
    glfwTerminate();

    return EXIT_SUCCESS;
}
