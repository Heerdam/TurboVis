
#include "../include/gui.hpp"

const char* Gui::GLGui::ImGui_ImplGlfw_GetClipboardText(void* user_data) {
    return glfwGetClipboardString((GLFWwindow*)user_data);
}

void Gui::GLGui::ImGui_ImplGlfw_SetClipboardText(void* user_data, const char* text) {
    glfwSetClipboardString((GLFWwindow*)user_data, text);
}

void Gui::GLGui::ImGui_ImplGlfw_MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS && button >= 0 && button < IM_ARRAYSIZE(g_MouseJustPressed))
        g_MouseJustPressed[button] = true;
}

void Gui::GLGui::ImGui_ImplGlfw_ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    ImGuiIO& io = ImGui::GetIO();
    io.MouseWheelH += (float)xoffset;
    io.MouseWheel += (float)yoffset;
}

void Gui::GLGui::ImGui_ImplGlfw_KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    ImGuiIO& io = ImGui::GetIO();
    if (action == GLFW_PRESS)
        io.KeysDown[key] = true;
    if (action == GLFW_RELEASE)
        io.KeysDown[key] = false;

    // Modifiers are not reliable across systems
    io.KeyCtrl = io.KeysDown[GLFW_KEY_LEFT_CONTROL] || io.KeysDown[GLFW_KEY_RIGHT_CONTROL];
    io.KeyShift = io.KeysDown[GLFW_KEY_LEFT_SHIFT] || io.KeysDown[GLFW_KEY_RIGHT_SHIFT];
    io.KeyAlt = io.KeysDown[GLFW_KEY_LEFT_ALT] || io.KeysDown[GLFW_KEY_RIGHT_ALT];
#ifdef _WIN32
    io.KeySuper = false;
#else
    io.KeySuper = io.KeysDown[GLFW_KEY_LEFT_SUPER] || io.KeysDown[GLFW_KEY_RIGHT_SUPER];
#endif
}

void Gui::GLGui::ImGui_ImplGlfw_CharCallback(GLFWwindow* window, unsigned int c) {
    ImGuiIO& io = ImGui::GetIO();
    io.AddInputCharacter(c);
}

void Gui::GLGui::ImGui_ImplGlfw_UpdateMousePosAndButtons(GLFWwindow* _window) {
    // Update buttons
    ImGuiIO& io = ImGui::GetIO();
    for (int i = 0; i < IM_ARRAYSIZE(io.MouseDown); i++) {
        // If a mouse press event came, always pass it as "mouse held this frame", so we don't miss click-release events that are shorter than 1 frame.
        io.MouseDown[i] = g_MouseJustPressed[i] || glfwGetMouseButton(_window, i) != 0;
        g_MouseJustPressed[i] = false;
    }

    // Update mouse position
    const ImVec2 mouse_pos_backup = io.MousePos;
    io.MousePos = ImVec2(-FLT_MAX, -FLT_MAX);
#ifdef __EMSCRIPTEN__
    const bool focused = true;  // Emscripten
#else
    const bool focused = glfwGetWindowAttrib(_window, GLFW_FOCUSED) != 0;
#endif
    if (focused) {
        if (io.WantSetMousePos) {
            glfwSetCursorPos(_window, (double)mouse_pos_backup.x, (double)mouse_pos_backup.y);
        } else {
            double mouse_x, mouse_y;
            glfwGetCursorPos(_window, &mouse_x, &mouse_y);
            io.MousePos = ImVec2((float)mouse_x, (float)mouse_y);
        }
    }
}

void Gui::GLGui::initGui(GLFWwindow* _window, size_t _width, size_t _height) {
    //Fonts and size
	{
		ImGuiIO& io = ImGui::GetIO();
		io.Fonts->AddFontDefault();
		io.DisplaySize = ImVec2((float)_width, (float)_height);
	}
	//Inputs
	{
		ImGuiIO& io = ImGui::GetIO();
		io.BackendFlags |= ImGuiBackendFlags_HasMouseCursors;         // We can honor GetMouseCursor() values (optional)
		io.BackendFlags |= ImGuiBackendFlags_HasSetMousePos;          // We can honor io.WantSetMousePos requests (optional, rarely used)
		io.BackendPlatformName = "TurboGui";

		// Keyboard mapping. ImGui will use those indices to peek into the io.KeysDown[] array.
		io.KeyMap[ImGuiKey_Tab] = GLFW_KEY_TAB;
		io.KeyMap[ImGuiKey_LeftArrow] = GLFW_KEY_LEFT;
		io.KeyMap[ImGuiKey_RightArrow] = GLFW_KEY_RIGHT;
		io.KeyMap[ImGuiKey_UpArrow] = GLFW_KEY_UP;
		io.KeyMap[ImGuiKey_DownArrow] = GLFW_KEY_DOWN;
		io.KeyMap[ImGuiKey_PageUp] = GLFW_KEY_PAGE_UP;
		io.KeyMap[ImGuiKey_PageDown] = GLFW_KEY_PAGE_DOWN;
		io.KeyMap[ImGuiKey_Home] = GLFW_KEY_HOME;
		io.KeyMap[ImGuiKey_End] = GLFW_KEY_END;
		io.KeyMap[ImGuiKey_Insert] = GLFW_KEY_INSERT;
		io.KeyMap[ImGuiKey_Delete] = GLFW_KEY_DELETE;
		io.KeyMap[ImGuiKey_Backspace] = GLFW_KEY_BACKSPACE;
		io.KeyMap[ImGuiKey_Space] = GLFW_KEY_SPACE;
		io.KeyMap[ImGuiKey_Enter] = GLFW_KEY_ENTER;
		io.KeyMap[ImGuiKey_Escape] = GLFW_KEY_ESCAPE;
		io.KeyMap[ImGuiKey_KeyPadEnter] = GLFW_KEY_KP_ENTER;
		io.KeyMap[ImGuiKey_A] = GLFW_KEY_A;
		io.KeyMap[ImGuiKey_C] = GLFW_KEY_C;
		io.KeyMap[ImGuiKey_V] = GLFW_KEY_V;
		io.KeyMap[ImGuiKey_X] = GLFW_KEY_X;
		io.KeyMap[ImGuiKey_Y] = GLFW_KEY_Y;
		io.KeyMap[ImGuiKey_Z] = GLFW_KEY_Z;

        Gui::InputMultiplexer::mouseButtonCallback(ImGui_ImplGlfw_MouseButtonCallback);
        Gui::InputMultiplexer::scrollCallback(ImGui_ImplGlfw_ScrollCallback);
        Gui::InputMultiplexer::keyCallback(ImGui_ImplGlfw_KeyCallback);
        Gui::InputMultiplexer::charCallback(ImGui_ImplGlfw_CharCallback);
	}
}

Gui::InputMultiplexer* Gui::InputMultiplexer::instance = new Gui::InputMultiplexer();

void Gui::InputMultiplexer::init(GLFWwindow* _window){

    //glfwSetWindowUserPointer(_window, instance);
    auto cursorPos = [](GLFWwindow* _window, double _xpos, double _ypos)->void {
        for(auto& f : instance->cursorPosCallbacks)
            f.second(_window, _xpos, _ypos);
    };

    glfwSetCursorPosCallback(_window, cursorPos);

    auto key = [](GLFWwindow* _window, int _key, int _scancode, int _action, int _mods)->void {
        for(auto& f : instance->keyCallbacks)
            f.second(_window, _key, _scancode, _action, _mods);
    };

    glfwSetKeyCallback(_window, key);

    auto charc = [](GLFWwindow* _window, uint32_t _codepoint)->void {
        for(auto& f : instance->charCallbacks)
            f.second(_window, _codepoint);
    };

    glfwSetCharCallback(_window, charc);

    auto mb = [](GLFWwindow* _window, int _button, int _action, int _mods)->void {
        for(auto& f : instance->mouseButtonCallbacks)
            f.second(_window, _button, _action, _mods);
    };

    glfwSetMouseButtonCallback(_window, mb);

    auto scroll = [](GLFWwindow* _window, double _xoffset, double _yoffset)->void {
        for(auto& f : instance->scrollCallbacks)
            f.second(_window, _xoffset, _yoffset);
    };

    glfwSetScrollCallback(_window, scroll);

}

size_t Gui::InputMultiplexer::cursorPosCallback(std::function<void(GLFWwindow* _window, double _xpos, double _ypos)> _func){
    size_t tid = ++instance->id;
    tid = (tid | 0x80000000);
    instance->cursorPosCallbacks.emplace_back(tid, _func);
    return tid;
}

size_t Gui::InputMultiplexer::keyCallback(std::function<void(GLFWwindow* _window, int _key, int _scancode, int _action, int _mods)> _func){
    size_t tid = ++instance->id;
    tid = (tid | 0x40000000);
    instance->keyCallbacks.emplace_back(tid, _func);
    return tid;
}

size_t Gui::InputMultiplexer::charCallback(std::function<void(GLFWwindow* _window, uint32_t _codepoint)> _func){
    size_t tid = ++instance->id;
    tid = (tid | 0x20000000);
    instance->charCallbacks.emplace_back(tid, _func);
    return tid;
}

size_t Gui::InputMultiplexer::mouseButtonCallback(std::function<void(GLFWwindow* _window, int _button, int _action, int _mods)> _func){
    size_t tid = ++instance->id;
    tid = (tid | 0x10000000);
    instance->mouseButtonCallbacks.emplace_back(tid, _func);
    return tid;
}   

size_t Gui::InputMultiplexer::scrollCallback(std::function<void(GLFWwindow* _window, double _xoffset, double _yoffset)> _func){
    size_t tid = ++instance->id;
    tid = (tid | 0x8000000);
    instance->scrollCallbacks.emplace_back(tid, _func);
    return tid;
}

void Gui::InputMultiplexer::remove(size_t _id){
    size_t tid = (_id >> 27) & 0x1F;
    switch(tid){
        case 0x1:
        {
            for(auto it = instance->scrollCallbacks.begin(); it != instance->scrollCallbacks.end(); ++it){
                if((*it).first == _id){
                    instance->scrollCallbacks.erase(it);
                    return;
                }
            }
        }
        return;
        case 0x2:
        {
            for(auto it = instance->mouseButtonCallbacks.begin(); it != instance->mouseButtonCallbacks.end(); ++it){
                if((*it).first == _id){
                    instance->mouseButtonCallbacks.erase(it);
                    return;
                }
            }
        }
        return;
        case 0x4:
        {
            for(auto it = instance->charCallbacks.begin(); it != instance->charCallbacks.end(); ++it){
                if((*it).first == _id){
                    instance->charCallbacks.erase(it);
                    return;
                }
            }
        }
        return;
        case 0x8:
        {
            for(auto it = instance->keyCallbacks.begin(); it != instance->keyCallbacks.end(); ++it){
                if((*it).first == _id){
                    instance->keyCallbacks.erase(it);
                    return;
                }
            }
        }
        return;
        case 0x10:
        {
            for(auto it = instance->cursorPosCallbacks.begin(); it != instance->cursorPosCallbacks.end(); ++it){
                if((*it).first == _id){
                    instance->cursorPosCallbacks.erase(it);
                    return;
                }
            }
        }
        return;
    }
}
