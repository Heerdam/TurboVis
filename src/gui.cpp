
#include "../include/gui.hpp"
#include "../include/util.hpp"
#include "../include/gl.hpp"

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

Gui::FrontendGui::FrontendGui(GLFWwindow* _window, size_t _width, size_t _height) { 

    ImGuiIO& io = ImGui::GetIO();
    std::stringstream ss;
    ss << FilePathResolver::ASSETDIR() << "ARIAL.TTF";
    io.Fonts->AddFontFromFileTTF(ss.str().c_str(), 30);
	io.DisplaySize = ImVec2(float(_width), float(_height));  

    gui.initGui(_window, _width, _height);
    //setStyle();
    try {
        gui.initGL(uint32_t(100000), uint32_t(200000));
    } catch (const TurboGUI::TurboGuiException& e) {
        spdlog::error("{}", e.what());
    }
}

void Gui::FrontendGui::draw(GLFWwindow* _window, RenderInfo& _info) {
    gui.ImGui_ImplGlfw_UpdateMousePosAndButtons(_window);
    gui.begin();
    drawTopMenu(_info);
    drawFPS(_info);
    if(window_render)
        drawRenderWindow(_info);
    if(window_debug)
        ImGui::ShowMetricsWindow(&window_debug);
    gui.draw();
    gui.sync();
}

void Gui::FrontendGui::drawTopMenu(const RenderInfo& _info){
    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(_info.width - 200, 50), ImGuiCond_Always);
    ImGui::Begin("menu", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_MenuBar |
        ImGuiWindowFlags_AlwaysAutoResize);

    if (ImGui::BeginMenuBar()) {
        if(ImGui::BeginMenu("Menu")) {

            ImGui::MenuItem("New");
            ImGui::MenuItem("Open", "Ctrl+O");

            ImGui::EndMenu();
        }
        
        if(ImGui::Button("Settings"))
            window_settings = true;
        
        if(ImGui::Button("Render"))
            window_render = true;

        if(ImGui::Button("Debug"))
            window_debug = true;
        
    }
    
    ImGui::EndMenuBar();

    ImGui::End();
}

void Gui::FrontendGui::drawRenderWindow(RenderInfo& _info){
    ImGui::Begin("renderer", &window_render);
    ImGui::SliderFloat("tr_fac", &_info.uniforms->tr_fac, 0.0f, 100000.f);
    ImGui::SliderInt("steps", &_info.uniforms->steps, 0, 10000);
    ImGui::Checkbox("grayscale", &_info.uniforms->grayscale);
    ImGui::End();
}

void Gui::FrontendGui::drawFPS(const RenderInfo& _info){
    ImGui::SetNextWindowPos(ImVec2(_info.width - 200, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(200, 50), ImGuiCond_Always);
    ImGui::Begin("fps", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoNav | 
        ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoInputs);

        const float frac = std::clamp(0.65f / 144.f * float(_info.fps), 0.f, 0.65f);
        ImVec4 col = ImVec4(1.f - frac, frac, 0.f, 1.f);
        ImGui::TextColored(col, "FPS: %i [%i]", _info.fps, _info.maxfps);

    ImGui::End();
}

void Gui::FrontendGui::setStyle(){
    ImGuiStyle& style = ImGui::GetStyle();

    style.WindowRounding    = 2.0f;
    style.ScrollbarRounding = 3.0f;
    style.GrabRounding      = 2.0f;
    style.AntiAliasedLines  = true;
    style.AntiAliasedFill   = true;
    style.WindowRounding    = 2;
    style.ChildRounding     = 2;
    style.ScrollbarSize     = 16;
    style.ScrollbarRounding = 3;
    style.GrabRounding      = 2;
    style.ItemSpacing.x     = 10;
    style.ItemSpacing.y     = 4;
    style.IndentSpacing     = 22;
    style.FramePadding.x    = 6;
    style.FramePadding.y    = 4;
    style.Alpha             = 1.0f;
    style.FrameRounding     = 3.0f;

    style.Colors[ImGuiCol_Text]                     = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
    style.Colors[ImGuiCol_TextDisabled]             = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
    style.Colors[ImGuiCol_WindowBg]                 = ImVec4(0.86f, 0.86f, 0.86f, 1.00f);
    //style.Colors[ImGuiCol_ChildWindowBg]          = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    style.Colors[ImGuiCol_ChildBg]                  = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    style.Colors[ImGuiCol_PopupBg]                  = ImVec4(0.93f, 0.93f, 0.93f, 0.98f);
    style.Colors[ImGuiCol_Border]                   = ImVec4(0.71f, 0.71f, 0.71f, 0.08f);
    style.Colors[ImGuiCol_BorderShadow]             = ImVec4(0.00f, 0.00f, 0.00f, 0.04f);
    style.Colors[ImGuiCol_FrameBg]                  = ImVec4(0.71f, 0.71f, 0.71f, 0.55f);
    style.Colors[ImGuiCol_FrameBgHovered]           = ImVec4(0.94f, 0.94f, 0.94f, 0.55f);
    style.Colors[ImGuiCol_FrameBgActive]            = ImVec4(0.71f, 0.78f, 0.69f, 0.98f);
    style.Colors[ImGuiCol_TitleBg]                  = ImVec4(0.85f, 0.85f, 0.85f, 1.00f);
    style.Colors[ImGuiCol_TitleBgCollapsed]         = ImVec4(0.82f, 0.78f, 0.78f, 0.51f);
    style.Colors[ImGuiCol_TitleBgActive]            = ImVec4(0.78f, 0.78f, 0.78f, 1.00f);
    style.Colors[ImGuiCol_MenuBarBg]                = ImVec4(0.86f, 0.86f, 0.86f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarBg]              = ImVec4(0.20f, 0.25f, 0.30f, 0.61f);
    style.Colors[ImGuiCol_ScrollbarGrab]            = ImVec4(0.90f, 0.90f, 0.90f, 0.30f);
    style.Colors[ImGuiCol_ScrollbarGrabHovered]     = ImVec4(0.92f, 0.92f, 0.92f, 0.78f);
    style.Colors[ImGuiCol_ScrollbarGrabActive]      = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
    style.Colors[ImGuiCol_CheckMark]                = ImVec4(0.184f, 0.407f, 0.193f, 1.00f);
    style.Colors[ImGuiCol_SliderGrab]               = ImVec4(0.26f, 0.59f, 0.98f, 0.78f);
    style.Colors[ImGuiCol_SliderGrabActive]         = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
    style.Colors[ImGuiCol_Button]                   = ImVec4(0.71f, 0.78f, 0.69f, 0.40f);
    style.Colors[ImGuiCol_ButtonHovered]            = ImVec4(0.725f, 0.805f, 0.702f, 1.00f);
    style.Colors[ImGuiCol_ButtonActive]             = ImVec4(0.793f, 0.900f, 0.836f, 1.00f);
    style.Colors[ImGuiCol_Header]                   = ImVec4(0.71f, 0.78f, 0.69f, 0.31f);
    style.Colors[ImGuiCol_HeaderHovered]            = ImVec4(0.71f, 0.78f, 0.69f, 0.80f);
    style.Colors[ImGuiCol_HeaderActive]             = ImVec4(0.71f, 0.78f, 0.69f, 1.00f);
    //style.Colors[ImGuiCol_Column]                 = ImVec4(0.39f, 0.39f, 0.39f, 1.00f);
    //style.Colors[ImGuiCol_ColumnHovered]          = ImVec4(0.26f, 0.59f, 0.98f, 0.78f);
    //style.Colors[ImGuiCol_ColumnActive]           = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
    style.Colors[ImGuiCol_Separator]                = ImVec4(0.39f, 0.39f, 0.39f, 1.00f);
    style.Colors[ImGuiCol_SeparatorHovered]         = ImVec4(0.14f, 0.44f, 0.80f, 0.78f);
    style.Colors[ImGuiCol_SeparatorActive]          = ImVec4(0.14f, 0.44f, 0.80f, 1.00f);
    style.Colors[ImGuiCol_ResizeGrip]               = ImVec4(1.00f, 1.00f, 1.00f, 0.00f);
    style.Colors[ImGuiCol_ResizeGripHovered]        = ImVec4(0.26f, 0.59f, 0.98f, 0.45f);
    style.Colors[ImGuiCol_ResizeGripActive]         = ImVec4(0.26f, 0.59f, 0.98f, 0.78f);
    style.Colors[ImGuiCol_PlotLines]                = ImVec4(0.39f, 0.39f, 0.39f, 1.00f);
    style.Colors[ImGuiCol_PlotLinesHovered]         = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
    style.Colors[ImGuiCol_PlotHistogram]            = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
    style.Colors[ImGuiCol_PlotHistogramHovered]     = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
    style.Colors[ImGuiCol_TextSelectedBg]           = ImVec4(0.26f, 0.59f, 0.98f, 0.35f);
    //style.Colors[ImGuiCol_ModalWindowDarkening]   = ImVec4(0.20f, 0.20f, 0.20f, 0.35f);
    style.Colors[ImGuiCol_DragDropTarget]           = ImVec4(0.26f, 0.59f, 0.98f, 0.95f);
    style.Colors[ImGuiCol_NavHighlight]             = style.Colors[ImGuiCol_HeaderHovered];
    style.Colors[ImGuiCol_NavWindowingHighlight]    = ImVec4(0.70f, 0.70f, 0.70f, 0.70f);
}
