
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

    auto window = glfwCreateWindow(1920, 1080, "TurboVis 0.1a", NULL, NULL);
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
    Gui::Gui gui;
    gui.initGui(window, 1920, 1080);
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

	Data::FunctionData data;
	data.createDataFromFunction([](int16_t _x, int16_t _y, int16_t _z)-> float {
		return 1.f;
	}, Data::Gridsize::x4);

    // -------------- CAMERA --------------
    GL::Camera camera;
    camera.position = { 500.f, 0.f, 0.f };
    camera.size = { 1920.f, 1080.f };
    camera.combined.setIdentity();



    using vec3 = glm::vec3;

    //camera.combined *= GL::Camera::lookAt(camera.position, {0.f, 0.f, 0.f}, {0.f, 0.f, 1.f});
    vec3 v1 = vec3( -500.f, 0.f, 0.f );
    vec3 v2 = vec3(0.f, 0.f, 0.f);
    vec3 v3 = vec3(0.f, 0.f, 1.f);
    auto t1 = bench_function( glm::lookAt, v1, v2, v3);

    Eigen::Vector3f w1 = {-500.f, 0.f, 0.f};
    Eigen::Vector3f w2 = {0.f, 0.f, 0.f};
    Eigen::Vector3f w3 = {0.f, 0.f, 1.f};

    auto t2 = bench_function( GL::Camera::lookAt, w1, w2, w3);
    //std::cout << t1 << std::endl;
    std::cout << t2 << std::endl;
    //std::cout << la1 << std::endl;
    //std::cout << glm::to_string(la2) << std::endl;

    return 0;

    // -------------- SHADER --------------
    GL::ShaderProgram shader;
    shader.compileFromFile("shader/shader");

    // -------------- DATA --------------
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
        glUniformMatrix4fv(1, 1, false, camera.combined.data());
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