#ifndef GUI_HPP
#define GUI_HPP

#include "tb_gui.h"

namespace Gui {

    class Gui : public TurboGUI::GUI {

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

}

#endif //GUI_HPP