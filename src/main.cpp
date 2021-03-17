
#include "../include/defines.hpp"

#if NDEBUG
void GLAPIENTRY MessageCallback(GLenum /*source*/, GLenum type, GLuint /*id*/, GLenum severity, GLsizei /*length*/, const GLchar* message, const void* /*userParam*/) {
	if (type != GL_DEBUG_TYPE_ERROR) return;
	spdlog::error("GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n", (type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""), type, severity, message);
}
#endif

int main() {

    if (!glfwInit()) {
		spdlog::error("Failed to init glfw. Shutting down...", true);
		return EXIT_FAILURE;
	}

	//glfwWindowHint(GLFW_REFRESH_RATE, 144);
	glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#if NDEBUG
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
#endif

	glfwWindowHint(GLFW_DEPTH_BITS, 24);
	glfwWindowHint(GLFW_STENCIL_BITS, 8);

    auto window = glfwCreateWindow(1920, 1080, "TurboVis 0.1a", NULL, NULL);
    if(!window){
        spdlog::error("Failed to create window. Shutting down...", true);
		return EXIT_FAILURE;
	}

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		spdlog::error("ERROR:\tFailed to init GLAD. Shutting down...", true);
		return EXIT_FAILURE;
	}

    glEnable(GL_STENCIL_TEST);
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_SRGB);

    double time = glfwGetTime();
    double deltaTime = 0.;
	unsigned long long frame = 0;
    size_t frames = 0;

    while (!glfwWindowShouldClose(window)){

        const double ctime = glfwGetTime();

        glfwPollEvents();

        frame++;

			if (ctime - time >= 1.0) {
				frame = 0;
				time = ctime;
			}

        glfwSwapBuffers(window);
		deltaTime = glfwGetTime() - ctime;
    }
    
    
    return EXIT_SUCCESS;
}