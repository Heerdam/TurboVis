#ifndef GUI_HPP
#define GUI_HPP

#include "tb_gui.h"

namespace GL {
    class Uniforms;
}

namespace Gui {

    class GLGui : public TurboGUI::GUI {

        static inline bool g_MouseJustPressed[ImGuiMouseButton_COUNT] = {};

        const char* ImGui_ImplGlfw_GetClipboardText(void*);
        void ImGui_ImplGlfw_SetClipboardText(void*, const char*);
        static void ImGui_ImplGlfw_MouseButtonCallback(GLFWwindow*, int, int, int);
        static void ImGui_ImplGlfw_ScrollCallback(GLFWwindow*, double, double);
        static void ImGui_ImplGlfw_KeyCallback(GLFWwindow*, int, int, int, int);     
        static void ImGui_ImplGlfw_CharCallback(GLFWwindow*, unsigned int);

       public:
        void initGui(GLFWwindow*, size_t, size_t);
        void ImGui_ImplGlfw_UpdateMousePosAndButtons(GLFWwindow*);
    };

    class InputMultiplexer {

        size_t id = 1;

        std::vector<std::pair<size_t, std::function<void(GLFWwindow*, double /*xpos*/, double /*ypos*/)>>> cursorPosCallbacks;
        std::vector<std::pair<size_t, std::function<void(GLFWwindow*, int /*key*/, int /*scancode*/, int /*action*/, int /*mods*/)>>> keyCallbacks;
        std::vector<std::pair<size_t, std::function<void(GLFWwindow*, uint32_t /*codepoint*/)>>> charCallbacks;
        std::vector<std::pair<size_t, std::function<void(GLFWwindow*, int /*button*/, int /*action*/, int /*mods*/)>>> mouseButtonCallbacks;
        std::vector<std::pair<size_t, std::function<void(GLFWwindow*, double /*xoffset*/, double /*yoffset*/)>>> scrollCallbacks;

        InputMultiplexer(){}
        static InputMultiplexer* instance;

        //static inline std::function<void(GLFWwindow*, double /*xpos*/, double /*ypos*/)> cursorPos;
    public:
        static void init(GLFWwindow*);
        static size_t cursorPosCallback(std::function<void(GLFWwindow*, double /*xpos*/, double /*ypos*/)>);
        static size_t keyCallback(std::function<void(GLFWwindow*, int /*key*/, int /*scancode*/, int /*action*/, int /*mods*/)>);
        static size_t charCallback(std::function<void(GLFWwindow*, uint32_t /*codepoint*/)>);
        static size_t mouseButtonCallback(std::function<void(GLFWwindow*, int /*button*/, int /*action*/, int /*mods*/)>);
        static size_t scrollCallback(std::function<void(GLFWwindow*, double /*xoffset*/, double /*yoffset*/)>);
        static void remove(size_t /*id*/);
    };


    struct RenderInfo {
        size_t fps, maxfps, width, height;
        GL::Uniforms* uniforms;
    };

    //template<class Callback>
    class FrontendGui {

        // ----------------- TRANSIENT STATES -----------------
        bool window_settings = false;
        bool window_render = false;
        bool window_debug = false;

        // ----------------------------------------------------

        Gui::GLGui gui;

        void drawTopMenu(const RenderInfo&);
        void drawFPS(const RenderInfo&);
        void drawRenderWindow(RenderInfo&);

        void setStyle();

    public:
        FrontendGui(GLFWwindow*, size_t, size_t);
        void draw(GLFWwindow*, RenderInfo&);

        // ----------------- CALLBACKS -----------------


    };

}

#endif //GUI_HPP