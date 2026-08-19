#ifndef TW_GLFW2_H
#define TW_GLFW2_H

// Adapter for the single-window GLFW2 API used by these examples.  It keeps
// callback plumbing and framebuffer-vs-window scaling in one place while all
// windowing operations are provided by the vendored GLFW 2.7.9 library.
#define GLFW_NO_GLU
#include <GL/glfw.h>
#include <AntTweakBar.h>

// GLFW2 has no window user pointer, so the adapter stores callback context here.
typedef struct TwGLFW2Window { void *user_pointer; } GLFWwindow;
typedef void (*TwGLFW2ErrorFun)(int, const char *);
typedef void (*TwGLFW2KeyFun)(GLFWwindow *, int, int, int, int);
typedef void (*TwGLFW2CharFun)(GLFWwindow *, unsigned int);
typedef void (*TwGLFW2MouseButtonFun)(GLFWwindow *, int, int, int);
typedef void (*TwGLFW2MousePosFun)(GLFWwindow *, double, double);
typedef void (*TwGLFW2ScrollFun)(GLFWwindow *, double, double);
typedef void (*TwGLFW2WindowSizeFun)(GLFWwindow *, int, int);

static GLFWwindow tw_glfw2_window;
static double tw_glfw2_scale_x = 1.0, tw_glfw2_scale_y = 1.0;
static int tw_glfw2_scale_ready;
static int tw_glfw2_wheel;
static TwGLFW2KeyFun tw_glfw2_key_fun;
static TwGLFW2CharFun tw_glfw2_char_fun;
static TwGLFW2MouseButtonFun tw_glfw2_button_fun;
static TwGLFW2MousePosFun tw_glfw2_pos_fun;
static TwGLFW2ScrollFun tw_glfw2_scroll_fun;
static TwGLFW2WindowSizeFun tw_glfw2_size_fun;

static int tw_glfw2_mods(void)
{
    int mods = 0;
    if (glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS || glfwGetKey(GLFW_KEY_RSHIFT) == GLFW_PRESS) mods |= 1;
    if (glfwGetKey(GLFW_KEY_LCTRL)  == GLFW_PRESS || glfwGetKey(GLFW_KEY_RCTRL)  == GLFW_PRESS) mods |= 2;
    if (glfwGetKey(GLFW_KEY_LALT)   == GLFW_PRESS || glfwGetKey(GLFW_KEY_RALT)   == GLFW_PRESS) mods |= 4;
    return mods;
}

static void tw_glfw2_measure_scale(void)
{
    if (!tw_glfw2_scale_ready) {
        int width = 1, height = 1;
        GLint viewport[4];
        glfwGetWindowSize(&width, &height);
        glGetIntegerv(GL_VIEWPORT, viewport);
        if (width > 0 && height > 0 && viewport[2] > 0 && viewport[3] > 0) {
            tw_glfw2_scale_x = (double)viewport[2] / width;
            tw_glfw2_scale_y = (double)viewport[3] / height;
        }
        tw_glfw2_scale_ready = 1;
    }
}

static void GLFWCALL tw_glfw2_key_cb(int key, int action)
{ if (tw_glfw2_key_fun) tw_glfw2_key_fun(&tw_glfw2_window, key, 0, action, tw_glfw2_mods()); }
static void GLFWCALL tw_glfw2_char_cb(int codepoint, int action)
{ if (tw_glfw2_char_fun && action == GLFW_PRESS) tw_glfw2_char_fun(&tw_glfw2_window, (unsigned int)codepoint); }
static void GLFWCALL tw_glfw2_button_cb(int button, int action)
{ if (tw_glfw2_button_fun) tw_glfw2_button_fun(&tw_glfw2_window, button, action, tw_glfw2_mods()); }
static void GLFWCALL tw_glfw2_pos_cb(int x, int y)
{ if (tw_glfw2_pos_fun) tw_glfw2_pos_fun(&tw_glfw2_window, (double)x, (double)y); }
static void GLFWCALL tw_glfw2_wheel_cb(int position)
{
    if (tw_glfw2_scroll_fun) tw_glfw2_scroll_fun(&tw_glfw2_window, 0.0, (double)(position - tw_glfw2_wheel));
    tw_glfw2_wheel = position;
}
static void GLFWCALL tw_glfw2_size_cb(int width, int height)
{
    tw_glfw2_measure_scale();
    if (tw_glfw2_size_fun)
        tw_glfw2_size_fun(&tw_glfw2_window,
                          (int)(width * tw_glfw2_scale_x + 0.5),
                          (int)(height * tw_glfw2_scale_y + 0.5));
}

static GLFWwindow *tw_glfw2_create_window(int width, int height, const char *title, void *monitor, void *share)
{
    (void)monitor; (void)share;
    if (!glfwOpenWindow(width, height, 8, 8, 8, 8, 24, 0, GLFW_WINDOW)) return NULL;
    glfwSetWindowTitle(title);
    glfwDisable(GLFW_AUTO_POLL_EVENTS);
    glfwEnable(GLFW_KEY_REPEAT);
    tw_glfw2_window.user_pointer = NULL;
    return &tw_glfw2_window;
}
static void tw_glfw2_get_window_size(GLFWwindow *window, int *width, int *height)
{ (void)window; glfwGetWindowSize(width, height); }
static void tw_glfw2_get_framebuffer_size(GLFWwindow *window, int *width, int *height)
{
    int w, h; (void)window; tw_glfw2_measure_scale(); glfwGetWindowSize(&w, &h);
    *width = (int)(w * tw_glfw2_scale_x + 0.5); *height = (int)(h * tw_glfw2_scale_y + 0.5);
}
static void tw_glfw2_get_content_scale(GLFWwindow *window, float *xscale, float *yscale)
{ (void)window; tw_glfw2_measure_scale(); *xscale = (float)tw_glfw2_scale_x; *yscale = (float)tw_glfw2_scale_y; }
static void tw_glfw2_set_bar_size(TwBar *bar, int width, int height)
{
    int size[2];
    tw_glfw2_measure_scale();
    size[0] = (int)(width * tw_glfw2_scale_x + 0.5);
    size[1] = (int)(height * tw_glfw2_scale_y + 0.5);
    TwSetParam(bar, NULL, "size", TW_PARAM_INT32, 2, size);
}
static void tw_glfw2_set_bar_position(TwBar *bar, int x, int y)
{
    int position[2];
    tw_glfw2_measure_scale();
    position[0] = (int)(x * tw_glfw2_scale_x + 0.5);
    position[1] = (int)(y * tw_glfw2_scale_y + 0.5);
    TwSetParam(bar, NULL, "position", TW_PARAM_INT32, 2, position);
}
static void tw_glfw2_set_bar_values_width(TwBar *bar, int width)
{
    int scaled_width;
    tw_glfw2_measure_scale();
    scaled_width = (int)(width * tw_glfw2_scale_x + 0.5);
    TwSetParam(bar, NULL, "valueswidth", TW_PARAM_INT32, 1, &scaled_width);
}
static void tw_glfw2_get_cursor_pos(GLFWwindow *window, double *x, double *y)
{ int ix, iy; (void)window; glfwGetMousePos(&ix, &iy); *x = ix; *y = iy; }
static int tw_glfw2_window_should_close(GLFWwindow *window)
{ (void)window; return !glfwGetWindowParam(GLFW_OPENED); }
static void tw_glfw2_set_window_should_close(GLFWwindow *window, int close)
{ (void)window; if (close) glfwCloseWindow(); }
static void tw_glfw2_set_window_user_pointer(GLFWwindow *window, void *pointer)
{ window->user_pointer = pointer; }
static void *tw_glfw2_get_window_user_pointer(GLFWwindow *window)
{ return window->user_pointer; }
static void tw_glfw2_set_key(GLFWwindow *window, TwGLFW2KeyFun fun)
{ (void)window; tw_glfw2_key_fun = fun; glfwSetKeyCallback(tw_glfw2_key_cb); }
static void tw_glfw2_set_char(GLFWwindow *window, TwGLFW2CharFun fun)
{ (void)window; tw_glfw2_char_fun = fun; glfwSetCharCallback(tw_glfw2_char_cb); }
static void tw_glfw2_set_button(GLFWwindow *window, TwGLFW2MouseButtonFun fun)
{ (void)window; tw_glfw2_button_fun = fun; glfwSetMouseButtonCallback(tw_glfw2_button_cb); }
static void tw_glfw2_set_pos(GLFWwindow *window, TwGLFW2MousePosFun fun)
{ (void)window; tw_glfw2_pos_fun = fun; glfwSetMousePosCallback(tw_glfw2_pos_cb); }
static void tw_glfw2_set_scroll(GLFWwindow *window, TwGLFW2ScrollFun fun)
{ (void)window; tw_glfw2_scroll_fun = fun; glfwSetMouseWheelCallback(tw_glfw2_wheel_cb); }
static void tw_glfw2_set_size(GLFWwindow *window, TwGLFW2WindowSizeFun fun)
{ (void)window; tw_glfw2_size_fun = fun; glfwSetWindowSizeCallback(tw_glfw2_size_cb); }

#define GLFW_KEY_ESCAPE GLFW_KEY_ESC
#define GLFW_KEY_DELETE GLFW_KEY_DEL
#define GLFW_KEY_PAGE_UP GLFW_KEY_PAGEUP
#define GLFW_KEY_PAGE_DOWN GLFW_KEY_PAGEDOWN
#define GLFW_REPEAT GLFW_PRESS
#define GLFW_MOD_SHIFT 1
#define GLFW_MOD_CONTROL 2
#define GLFW_MOD_ALT 4
#define GLFW_TRUE 1
#define glfwCreateWindow tw_glfw2_create_window
#define glfwMakeContextCurrent(window) ((void)(window))
#define glfwGetWindowSize tw_glfw2_get_window_size
#define glfwGetFramebufferSize tw_glfw2_get_framebuffer_size
#define glfwGetWindowContentScale tw_glfw2_get_content_scale
#define glfwGetCursorPos tw_glfw2_get_cursor_pos
#define glfwWindowShouldClose tw_glfw2_window_should_close
#define glfwSetWindowShouldClose tw_glfw2_set_window_should_close
#define glfwSetWindowUserPointer tw_glfw2_set_window_user_pointer
#define glfwGetWindowUserPointer tw_glfw2_get_window_user_pointer
#define glfwSetKeyCallback tw_glfw2_set_key
#define glfwSetCharCallback tw_glfw2_set_char
#define glfwSetMouseButtonCallback tw_glfw2_set_button
#define glfwSetCursorPosCallback tw_glfw2_set_pos
#define glfwSetScrollCallback tw_glfw2_set_scroll
#define glfwSetFramebufferSizeCallback tw_glfw2_set_size
#define glfwSetErrorCallback(callback) ((void)(callback))
#define glfwWindowHint(hint, value) glfwOpenWindowHint((hint), (value))
#define glfwSwapBuffers(window) (glfwSwapBuffers)()

#endif
