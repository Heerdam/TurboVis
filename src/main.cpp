
#include "../include/defines.hpp"
#include "../include/gui.hpp"
#include "../include/data.hpp"
#include "../include/gl.hpp"
#include "../include/bench.hpp"

void GLAPIENTRY MessageCallback(GLenum /*source*/, GLenum type, GLuint /*id*/, GLenum severity, GLsizei /*length*/, const GLchar* message, const void* /*userParam*/) {
	if (type != GL_DEBUG_TYPE_ERROR) return;
	fprintf(stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
		(type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""),
		type, severity, message);
}
int main() {

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
    const float HALFHEIGHTH = HEIGHT * 0.5f;

    auto window = glfwCreateWindow(WIDTH, HEIGHT, "TurboVis 0.1a", NULL, NULL);
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

	glDebugMessageCallback(MessageCallback, 0);

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
    GL::Camera camera;
    camera.position = { 250.f, 0.f, 0.f };
    camera.combined = glm::perspective(65.f, 1920.f / 1080.f, 0.1f, 1000.f);
    camera.combined *= glm::lookAt(camera.position, {0.f, 0.f, 0.f}, {0.f, 0.f, 1.f});

    bool RMB_down = false;

    double mX, mY;
    glfwGetCursorPos(window, &mX, &mY);
    Vec2 oldPosition = Vec2((float)mY - HALFHEIGHTH, (float)mX - HALFWIDTH);
    Vec2 newPosition = oldPosition;

    Gui::InputMultiplexer::mouseButtonCallback([&](GLFWwindow*, int _button, int _action, int _mods)-> void{
        if(_button == GLFW_MOUSE_BUTTON_RIGHT && _mods == 0x0){
            RMB_down = _action == GLFW_PRESS ? true : _action == GLFW_RELEASE ? false : RMB_down;
        }
    });

    Gui::InputMultiplexer::cursorPosCallback([&](GLFWwindow*, double _xpos, double _ypos)-> void{
        oldPosition = newPosition;
        newPosition = Vec2((HEIGHT - (float)_ypos - HALFHEIGHTH), (float)_xpos - HALFWIDTH);

        if(RMB_down){
            const auto rot = GL::Camera::trackball_shoemake(oldPosition, newPosition, 250.f);
            const glm::mat4 RotationMatrix = glm::toMat4(rot);
            camera.combined *= RotationMatrix;
        }
    });

    // -------------- SHADER --------------
    GL::ShaderProgram shader;
    shader.compileFromFile("shader/shader");

    // -------------- DATA --------------
    Data::FunctionData data;
	data.createDataFromFunction([](int16_t _x, int16_t _y, int16_t _z)-> float {
		return 1.f;
	}, Data::Gridsize::x4);

	auto vals = data.toBuffer();
	auto verts = data.createBuffers();

	std::vector<uint32_t> indices (vals.size());
	for(size_t i = 0; i < vals.size(); ++i)
		indices[i] = (uint32_t)i;

	GLuint vao, vbo, ebo;

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(), GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

	glGenBuffers(1, &ebo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(uint32_t), indices.data(), GL_STATIC_DRAW);
	

    glBindVertexArray(0);




    while (!glfwWindowShouldClose(window)) {
        const double ctime = glfwGetTime();

        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        glfwPollEvents();
        gui.ImGui_ImplGlfw_UpdateMousePosAndButtons(window);

        // -------------- GEOMETRY PASS --------------
        shader.bind();
        glBindVertexArray(vao);
        glUniformMatrix4fv(1, 1, false, glm::value_ptr(camera.combined));
        glDrawElements(GL_POINTS, (GLsizei)indices.size(), GL_UNSIGNED_INT, nullptr);   
        shader.unbind();

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